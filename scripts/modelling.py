import os
import flopy
from flopy.utils.cvfdutil import area_of_polygon
from scipy.spatial.distance import jaccard

import discretization
from discretization import leapfrog_name
import numpy as np
import parameterization
import boundaries
from flopy.utils.gridgen import Gridgen
import utils
from project_base import unbacked_dir, project_dir


def create_run_model(name, case='base', time='current', scenario='summer_average', tscenario='steady', init='fresh', hypothesis='null',
                     patch_multiplier=None, conduit_k=None):
    mf6exe = '/bin/modflow/mf6.exe'
    length_units = "m"
    time_units = "days"

    t_params = parameterization.get_tdis_params(tscenario)


    nouter, ninner = 500, 500
    hclose, rclose, relax = 1e-5, 1e-5, 0.97

    params = parameterization.get_hydraulic_params(case)

    ws = unbacked_dir.joinpath('models', name)
    if not os.path.exists(ws):
        os.makedirs(ws)

    sim = flopy.mf6.MFSimulation(
        sim_name=name, sim_ws=ws, exe_name=mf6exe
    )

    time_spd = [(t_params['perlen'], t_params['nstp'], t_params['tsmult']) for i in range(t_params['nper'])]

    flopy.mf6.ModflowTdis(
        sim, nper=t_params['nper'], perioddata=time_spd, time_units=time_units
    )

    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)

    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="ALL",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        under_relaxation="NONE",
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="BICGSTAB",
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=relax,
        filename="{}.ims".format(gwf.name),
    )

    sim.register_ims_package(ims, [gwf.name])

    if leapfrog_name is None:
        grid = discretization.get_grid_w_ibound()
        flopy.mf6.ModflowGwfdis(
            gwf,
            length_units=length_units,
            nlay=grid.nlay,
            nrow=grid.nrow,
            ncol=grid.ncol,
            delr=grid.delr,
            delc=grid.delc,
            top=grid.top,
            botm=grid.botm,
            idomain=grid.idomain,
            xorigin=grid.xoffset,
            yorigin=grid.yoffset,
            angrot=grid.angrot
        )
    else:
        grid = discretization.get_grid_w_ibound()

        temp_sim = flopy.mf6.MFSimulation.load(sim_name='mfsim',sim_ws=project_dir.joinpath("data", "leapfrog", leapfrog_name))
        temp_model = temp_sim.get_model(leapfrog_name)
        temp_disu = temp_model.get_package("disu")

        disu = flopy.mf6.ModflowGwfdisu(
            gwf,
            nodes=temp_disu.nodes.data,
            nja=temp_disu.nja.data,
            nvert=temp_disu.nvert.data,
            top=temp_disu.top.data,
            bot=temp_disu.bot.data,
            area=temp_disu.area.data,
            iac=temp_disu.iac.data,
            ja=temp_disu.ja.data,
            ihc=temp_disu.ihc.data,
            cl12=temp_disu.cl12.data,
            hwva=temp_disu.hwva.data,
            angldegx=temp_disu.angldegx.data,
            vertices=temp_disu.vertices.array,
            cell2d=temp_disu.cell2d.array,
        )

    k = parameterization.get_params_by_zone(case, 'hk', hypothesis, patch_multiplier, conduit_k)
    flopy.mf6.ModflowGwfnpf(
        gwf,
        icelltype=[0] * grid.nnodes,
        save_flows=True,
        save_specific_discharge=True,
        k=k,
        k33overk=False,
        k33=parameterization.get_params_by_zone(case, 'vk', hypothesis, patch_multiplier, conduit_k),
    )

    if init == 'fresh':
        strt = 0
    else:
        data = utils.get_final_results(init)
        strt = data['head']

    flopy.mf6.ModflowGwfic(gwf, strt=strt)

    pd = [(0, 0.7, 0.0, "trans", "concentration")]

    flopy.mf6.ModflowGwfbuy(gwf, packagedata=pd)

    chd_data = boundaries.get_chd_data('seafloor', time, scenario, nper=t_params['nper'])
    #fsave
    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=chd_data,
        pname="CHD-1",
        auxiliary="CONCENTRATION",
    )

    onshore_ghb = boundaries.get_ghb_data('onshore', case, time, scenario, nper=t_params['nper'])
    offshore_ghb = boundaries.get_ghb_data('offshore', case, time, scenario, nper=t_params['nper'])
    ghb_data = {i: onshore_ghb[i] + offshore_ghb[i] for i in range(t_params['nper'])}

    flopy.mf6.ModflowGwfghb(
        gwf,
        stress_period_data=ghb_data,
        pname="GHB-1",
        auxiliary=["CONCENTRATION"]
    )

    head_filerecord = f"{name}.hds"
    budget_filerecord = f"{name}.bud"
    saverecord = {i: [("HEAD", "LAST"), ("BUDGET", "FREQUENCY", t_params['frequency']), ("BUDGET", "LAST")] for i in range(t_params['nper'])}

    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=head_filerecord,
        budget_filerecord=budget_filerecord,
        saverecord=saverecord,
    )

    flopy.mf6.ModflowGwfsto(gwf,
                            ss=parameterization.get_params_by_zone(case, 'ss'),
                            sy=0.3,
                            iconvert=1,
                            transient=True)
    gwt = flopy.mf6.ModflowGwt(sim, modelname="trans", save_flows=True)

    imsgwt = flopy.mf6.ModflowIms(
        sim,
        print_option="ALL",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        under_relaxation="NONE",
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="BICGSTAB",
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=relax,
        filename="{}.ims".format(gwt.name),
    )

    sim.register_ims_package(imsgwt, [gwt.name])

    if leapfrog_name is None:
        flopy.mf6.ModflowGwtdis(
            gwt,
            length_units=length_units,
            nlay=grid.nlay,
            nrow=grid.nrow,
            ncol=grid.ncol,
            delr=grid.delr,
            delc=grid.delc,
            top=grid.top,
            botm=grid.botm,
            idomain=grid.idomain,
            xorigin=grid.xoffset,
            yorigin=grid.yoffset,
            angrot=grid.angrot
        )
    else:
        flopy.mf6.ModflowGwtdisu(
            gwt,
            nodes=temp_disu.nodes.data,
            nja=temp_disu.nja.data,
            nvert=temp_disu.nvert.data,
            top=temp_disu.top.data,
            bot=temp_disu.bot.data,
            area=temp_disu.area.data,
            iac=temp_disu.iac.data,
            ja=temp_disu.ja.data,
            ihc=temp_disu.ihc.data,
            cl12=temp_disu.cl12.data,
            hwva=temp_disu.hwva.data,
            angldegx=temp_disu.angldegx.data,
            vertices=temp_disu.vertices.array,
            cell2d=temp_disu.cell2d.array,
        )


    cnc_data = boundaries.get_cnc_data('seafloor', nper=t_params['nper'])
    flopy.mf6.ModflowGwtcnc(
        gwt,
        stress_period_data=cnc_data,
        pname="CNC-1",
    )

    # todo improve porosity
    flopy.mf6.ModflowGwtmst(gwt, porosity=0.5)

    if init == 'fresh':
        strt = 0
    else:
        data = utils.get_final_results(init)
        strt = data['conc']

    flopy.mf6.ModflowGwtic(gwt, strt=strt)

    flopy.mf6.ModflowGwtadv(gwt, scheme="UPSTREAM")

    flopy.mf6.ModflowGwtdsp(gwt, xt3d_off=True,
                            diffc=params['dm_coeff'],
                            alh=params['alpha'],
                            ath1=params['alpha'] * 0.1,
                            alv=params['alpha'] * 0.01)

    sourcerecarray = [
        ("GHB-1", "AUX", "CONCENTRATION"),
        ("CHD-1", "AUX", "CONCENTRATION")
    ]

    flopy.mf6.ModflowGwtssm(gwt, sources=sourcerecarray)

    saverecord = {i: [("CONCENTRATION", "FREQUENCY", t_params['frequency']), ("BUDGET", "LAST"), ("CONCENTRATION", "LAST")] for i in range(t_params['nper'])}

    flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord=f"{name}.cbc",
        concentration_filerecord=f"{name}.ucn",
        concentrationprintrecord=[
            ("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")
        ],
        saverecord=saverecord,
        printrecord=[("CONCENTRATION", "LAST"), ("BUDGET", "LAST")],
    )

    flopy.mf6.ModflowGwfgwt(
        sim, exgtype="GWF6-GWT6", exgmnamea=gwf.name, exgmnameb=gwt.name
    )

    sim.write_simulation(silent=False)
    success, buff = sim.run_simulation(silent=False)
    if not success:
        print(buff)

    return success

def load_model(name):
    ws = unbacked_dir.joinpath('models', name)
    sim = flopy.mf6.MFSimulation.load(sim_name=name, sim_ws=ws)
    return sim

if __name__ == "__main__":
    create_run_steady_state("steady_cbc", case='base', scenario='summer', tscenario='steady')
