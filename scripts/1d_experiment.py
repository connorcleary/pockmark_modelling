import numpy as np
import matplotlib.pyplot as plt
import flopy
from project_base import unbacked_dir, project_dir

def build_run_model(name, Kv, delz, thickness, flux):
    # discretize
    nlay = int(thickness / delz) + 1 # +1 for seafloor,
    nrow = 1
    ncol = 1
    delr = 1
    delc = 1
    top = 0.1
    botm = [-delz*i for i in range(nlay)]

    # build flow model
    sim = flopy.mf6.MFSimulation(sim_name=name, exe_name='mf6', version='mf6', sim_ws=unbacked_dir.joinpath('1d_models', name))

    tdis = flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1e6, 1e3, 1)], time_units="days")
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)
    dis = flopy.mf6.ModflowGwfdis(gwf, nlay=nlay, nrow=nrow, ncol=ncol,
                                  delr=delr, delc=delc, top=top, botm=botm, length_units="m")
    sto = flopy.mf6.ModflowGwfsto(gwf, save_flows=True, iconvert=1, ss=1e-5, sy=0.5, steady_state={0: True})

    ic = flopy.mf6.ModflowGwfic(gwf, strt=0.1)
    npf = flopy.mf6.ModflowGwfnpf(gwf, save_specific_discharge=True, icelltype=1, k=Kv)
    # constant head on seafloor
    chd_data = [(0, 0, 0, 0.1, 35.0)] # top layer
    chd = flopy.mf6.ModflowGwfchd(gwf, stress_period_data={0: chd_data}, auxiliary='CONCENTRATION', pname='CHD-1')
    # fixed flux on base (with concentration zero)
    wel_data = [(nlay-1, 0, 0, flux, 0.0)] # last layer
    wel = flopy.mf6.ModflowGwfwel(gwf, stress_period_data={0: wel_data}, auxiliary='CONCENTRATION', pname='WEL-1')

    # build transport model
    gwt = flopy.mf6.ModflowGwt(sim, modelname=name+'_gwt', save_flows=True)
    dis = flopy.mf6.ModflowGwtdis(gwt, nlay=nlay, nrow=nrow, ncol=ncol,
                                  delr=delr, delc=delc, top=top, botm=botm, length_units="m")
    ic = flopy.mf6.ModflowGwtic(gwt, strt=0.0)
    adv = flopy.mf6.ModflowGwtadv(gwt, scheme='TVD')
    dsp = flopy.mf6.ModflowGwtdsp(gwt, alh=0, ath1=0, atv=0, diffc=8.64e-5)
    cnc_data = [(0, 0, 0, 35.0)] # top layer
    # add constant concentration to base of model
    cnc_data.append((nlay-1, 0, 0, 0.0)) # bottom layer
    # constant concentration on seafloor
    cnc = flopy.mf6.ModflowGwtcnc(gwt, stress_period_data={0: cnc_data})
    mst = flopy.mf6.ModflowGwtmst(gwt, porosity=0.5)
    ssm_data = [("CHD-1", "AUX", "CONCENTRATION"), ("WEL-1", "AUX", "CONCENTRATION")]
    ssm = flopy.mf6.ModflowGwtssm(gwt, sources=ssm_data)

    # get concentration observations in all layers
    obs = []
    for ilay in range(nlay):
        obs.append([f'layer{ilay}', "concentration", (ilay, 0, 0)])
    flopy.mf6.ModflowUtlobs(gwt,
                            print_input=False,
                            continuous={f"{name+'_gwt'}.obs.csv": obs},
                            pname="OBS",
                            filename=f"{name+'_gwt'}.obs")

    # couple models

    # create solver
    ims_gwf = flopy.mf6.ModflowIms(sim, pname='ims_gwf', complexity='SIMPLE', filename=f'{name}.ims')
    ims_gwt = flopy.mf6.ModflowIms(sim, pname='ims_gwt', complexity='SIMPLE', filename=f'{name}_gwt.ims', linear_acceleration='BICGSTAB')

    # Register IMS packages with the correct model names
    sim.register_ims_package(ims_gwf, [name])  # GWF IMS first
    sim.register_ims_package(ims_gwt, [name + '_gwt'])

    # GWF-GWT exchange
    flopy.mf6.ModflowGwfgwt(sim, exgtype='GWF6-GWT6', exgmnamea=name, exgmnameb=name+'_gwt')

    # write and run
    sim.write_simulation()
    sim.run_simulation()

def run_scenarios():

    for thickness, lt in zip([5, 15], ['low', 'high']):
        for flux, lf in zip([1e-5, 2e-4], ['low', 'high']):
            for delz, ld in zip([0.01, 0.5], ['low', 'high']):
                name = f't{lt[0]}_f{lf[0]}_d{ld[0]}2'
                build_run_model(name, 1, delz, thickness, flux)

def plot_results():
    # plot eight subplots with the concentration profiles for delz=0.05 and delz=0.5, and the 0.05 averaged up to 0.5
    fig, axs = plt.subplots(1, 4, figsize=(8, 4), sharey='row')


    for i, (thickness, lt) in enumerate(zip([5, 15], ['low', 'high'])):
            for k, (flux, lf) in enumerate(zip([1e-5, 2e-4], ['low', 'high'])):
                ax = axs[i*2 + k]
                values = []
                for delz, ld in zip([0.01, 0.5], ['low', 'high']):
                    name = f't{lt[0]}_f{lf[0]}_d{ld[0]}2'
                    obs = np.loadtxt(unbacked_dir.joinpath('1d_models', name, f'{name}_gwt.obs.csv'), delimiter=',', skiprows=1)
                    # middle of each layer excluding top layer
                    nlay = obs.shape[1] - 1
                    center = np.array([-(delz*(i+0.5)) for i in range(nlay-1)])
                    # concat with [0.1] at top
                    center = np.concat([[0.5], center])
                    # plot 0.05 as line
                    # plot 0.5 as points
                    if delz == 0.5:
                        ax.scatter(obs[-1, 1:], center, label=f'dz = {delz} m', s=10, marker='x', color='k', linewidth=0.5, zorder=3)
                        print(obs[-1, 2])
                    else:
                        ax.plot(obs[-1, 1:], center, label=f'dz = {delz} m', lw=0.5, color='blue')

                    # add peru colored marker v at -0.25 m
                    if delz == 0.5:
                        ax.scatter(obs[-1, 2], -0.25, color='peru', marker='v', s=50, label='- 0.25 m', zorder=5, facecolor='none')

                    # add value at 0.25 to values for both models
                    if delz == 0.5:
                        values.append(obs[-1, 2])
                    else:
                        # get value at -0.25 m by linear interpolation
                        #V v = np.interp(-0.25, center, obs[-1, 2:])
                        values.append(0)

                ax.set_title(f'D={thickness} m, \n qz={flux} m/day', fontsize=7)
                ax.set_xlabel('Concentration [PSU]')
                ax.set_box_aspect(2)
                ax.set_xticks([0, 35])
                # put letter top left of each subplot
                ax.text(-0.1, 1.1, chr(97 + i*2 + k), transform=ax.transAxes, fontsize=10, fontweight='bold', va='top', ha='right')
                # add c error between 0.5 and 0.01 models at in bottom right of plot with white background
                # ax.text(0.95, 0.05, fr'$\epsilon_{-0.25}$ = {values[1]-values[0]:.2f} PSU', transform=ax.transAxes, fontsize=6, va='bottom', ha='right', backgroundcolor='white')


    axs[0].set_ylabel('Depth [m]')
    # add legend
    handles, labels = ax.get_legend_handles_labels()
    # put legend outside top right of plot
    ax.legend(handles, labels, loc='upper left', bbox_to_anchor=(1.05, 1), fontsize=7)

    plt.tight_layout()
    plt.savefig(unbacked_dir.joinpath('figures', f'1d_scenario_results.png'), dpi=600)



if __name__ == "__main__":
    # run_scenarios()
    plot_results()

