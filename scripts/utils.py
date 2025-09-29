import os
import flopy
import numpy as np
import pandas as pd
import discretization
from project_base import data_dir
import geopandas as gpd
import os
os.environ["SHAPE_RESTORE_SHX"] = "YES"
from shapely import LineString
from project_base import project_dir, unbacked_dir


def get_final_results(name, force=False):

    if force or not os.path.exists(unbacked_dir.joinpath('outputs', name, 'final_step.npz')):
        sim = flopy.mf6.MFSimulation.load(
            sim_ws=unbacked_dir.joinpath('models', name),
        )
        gwt = sim.get_model("trans")

        gwf = sim.get_model(name)
        bud = gwf.output.budget()
        qx, qy, qz = flopy.utils.postprocessing.get_specific_discharge(
            bud.get_data(text="DATA-SPDIS", totim=bud.get_times()[-1])[0],
            gwf,
        )
        head = gwt.output.head().get_data(totim=gwt.output.head().get_times()[-1])

        times = gwt.output.concentration().get_times()
        conc = gwt.output.concentration().get_data(totim=times[-1])
        cbud = gwt.output.concentration().get_data(totim=times[-1])

        os.makedirs(unbacked_dir.joinpath('outputs', name), exist_ok=True)
        np.savez_compressed(unbacked_dir.joinpath('outputs', name, 'final_step.npz'),
                            conc=conc, head=head, qx=qx, qy=qy, qz=qz)

        return {
                "conc": conc,
                "head": head,
                "qx": qx,
                "qy": qy,
                "qz": qz,
            }
    else:
        return np.load(unbacked_dir.joinpath('outputs', name, 'final_step.npz'))

def get_final_pockmark_results(name):

    if not os.path.exists(unbacked_dir.joinpath('outputs', name, 'pockmarks_final_step.h5')):
        data = get_final_results(name)
        pockmarks = discretization.get_pockmark_surfaces()
        df = pd.DataFrame(index=pockmarks.keys(), columns=["concentration", "discharge"])

        pockmark_connectivity = discretization.get_pockmark_connectivity()

        sim = flopy.mf6.MFSimulation.load(
            sim_ws=unbacked_dir.joinpath('models', name),
        )
        gwt = sim.get_model("trans")
        gwf = sim.get_model(name)
        bud = gwf.output.budget()

        cbud_available = True
        try:
            cbud = gwt.output.budget()
        except:
            cbud_available = False

        bud_data = bud.get_data(text=b'    FLOW-JA-FACE', totim=bud.get_times()[-1])[0][0][0]
        if cbud_available:
            cbud_data = cbud.get_data(text=b'    FLOW-JA-FACE', totim=cbud.get_times()[-1])[0][0][0]

        for pockmark in pockmarks:
            connectivity = np.array(pockmark_connectivity[pockmark])[:, 1]
            cells = np.array(pockmark_connectivity[pockmark])[:, 0]
            discharges = (bud_data[connectivity]) #  * (bud_data[connectivity] > 0))
            if cbud_available:
                concs = (data['conc'][0][0][cells]) # * (bud_data[connectivity] > 0))
                conc_flux = (cbud_data[connectivity] * (bud_data[connectivity] > 0))
            else:
                concs = np.ones_like(discharges)

            print("If new run need to recalculate with concentration budget")
            df.loc[pockmark, "concentration"] = np.sum([discharge*conc for discharge, conc in zip(discharges, concs)])/np.sum(discharges)
            df.loc[pockmark, "discharge"] = np.sum(discharges)

        df.to_hdf(unbacked_dir.joinpath('outputs', name, 'pockmarks_final_step.h5'), key='df', mode='w')
        return df
    else:
        return pd.read_hdf(unbacked_dir.joinpath('outputs', name, 'pockmarks_final_step.h5'))

def get_discharge_around_pockmarks(name, buffer=100):
    if not os.path.exists(unbacked_dir.joinpath('outputs', name, f'around_pockmarks_final_step_{buffer}.h5')):
        data = get_final_results(name)
        around_pockmarks = discretization.get_surface_around_pockmarks(buffer=buffer)
        df = pd.DataFrame(index=around_pockmarks.keys(), columns=["concentration", "discharge"])

        pockmark_connectivity = discretization.get_around_pockmark_connectivity(buffer=buffer)

        sim = flopy.mf6.MFSimulation.load(
            sim_ws=unbacked_dir.joinpath('models', name),
        )
        gwt = sim.get_model("trans")
        gwf = sim.get_model(name)
        bud = gwf.output.budget()

        cbud_available = True
        try:
            cbud = gwt.output.budget()
        except:
            cbud_available = False

        bud_data = bud.get_data(text=b'    FLOW-JA-FACE', totim=bud.get_times()[-1])[0][0][0]
        if cbud_available:
            cbud_data = cbud.get_data(text=b'    FLOW-JA-FACE', totim=cbud.get_times()[-1])[0][0][0]

        for pockmark in around_pockmarks:
            connectivity = np.array(pockmark_connectivity[pockmark])[:, 1]
            cells = np.array(pockmark_connectivity[pockmark])[:, 0]
            discharges = (bud_data[connectivity]) #  * (bud_data[connectivity] > 0))
            if cbud_available:
                concs = (data['conc'][0][0][cells]) # * (bud_data[connectivity] > 0))
                conc_flux = (cbud_data[connectivity] * (bud_data[connectivity] > 0))
            else:
                concs = np.ones_like(discharges)

            print("If new run need to recalculate with concentration budget")
            df.loc[pockmark, "concentration"] = np.sum([discharge*conc for discharge, conc in zip(discharges, concs)])/np.sum(discharges)
            df.loc[pockmark, "discharge"] = np.sum(discharges)

        df.to_hdf(unbacked_dir.joinpath('outputs', name, f'around_pockmarks_final_step_{buffer}.h5'), key='df', mode='w')
        return df
    else:
        return pd.read_hdf(unbacked_dir.joinpath('outputs', name, f'around_pockmarks_final_step_{buffer}.h5'))


def get_pockmark_discharge_and_recharge(name):

    if not os.path.exists(unbacked_dir.joinpath('outputs', name, 'pockmarks_discharge_recharge_over_time_step.h5')):
        pockmarks = discretization.get_pockmark_surfaces()
        pockmark_connectivity = discretization.get_pockmark_connectivity()

        sim = flopy.mf6.MFSimulation.load(
            sim_ws=unbacked_dir.joinpath('models', name),
        )
        gwt = sim.get_model("trans")
        gwf = sim.get_model(name)
        bud = gwf.output.budget()
        times = bud.get_times()

        results = {list(pockmarks.keys())[i]: pd.DataFrame(index=times, columns=["discharge", "recharge"]) for i in range(len(pockmarks.keys()))}
        if not os.path.exists(unbacked_dir.joinpath('outputs', name)):
            os.makedirs(unbacked_dir.joinpath('outputs', name))

        for t in times:
            bud_data = bud.get_data(text=b'    FLOW-JA-FACE', totim=t)[0][0][0]
            for pockmark in pockmarks:
                connectivity = np.array(pockmark_connectivity[pockmark])[:, 1]
                cells = np.array(pockmark_connectivity[pockmark])[:, 0]
                discharge = np.sum(bud_data[connectivity]*(bud_data[connectivity] > 0))
                recharge = np.sum(bud_data[connectivity]*(bud_data[connectivity] < 0))
                results[pockmark].loc[t, "discharge"] = discharge
                results[pockmark].loc[t, "recharge"] = recharge

        for pockmark in pockmarks:
            results[pockmark].to_hdf(unbacked_dir.joinpath('outputs', name, 'pockmarks_discharge_recharge_over_time_step.h5'), key=pockmark, mode='a')

        return results
    else:
        pockmarks = discretization.get_pockmark_surfaces()
        store = pd.HDFStore(unbacked_dir.joinpath('outputs', name, 'pockmarks_discharge_recharge_over_time_step.h5'), mode='r')
        return {pockmark: pd.read_hdf(unbacked_dir.joinpath('outputs', name, 'pockmarks_discharge_recharge_over_time_step.h5'), pockmark) for pockmark in pockmarks}

def get_total_seafloor_discharge(name):
    if not os.path.exists(unbacked_dir.joinpath('outputs', name, 'total_seafloor_discharge.npz')):
        data = get_final_results(name)
        non_pockmarks = discretization.get_non_pockmark_surfaces()
        qz = data['qz']
        discharges = qz[non_pockmarks['cells'].astype(int)]*non_pockmarks['areas']
        # assert np.all(discharges > 0)
        concs = data['conc'][0][0][non_pockmarks['cells'].astype(int)]
        total_discharge = np.sum(discharges*(discharges > 0))
        average_conc = np.sum(discharges*concs*(discharges > 0))/total_discharge
        np.savez_compressed(unbacked_dir.joinpath('outputs', name, 'total_seafloor_discharge.npz'),
                            discharge=total_discharge, conc=average_conc)
        return {
            "discharge": total_discharge,
            "conc": average_conc
        }

    else:
        return np.load(unbacked_dir.joinpath('outputs', name, 'total_seafloor_discharge.npz'))

def get_salinity_and_flux_at_measurement_cells(name):
    if not os.path.exists(unbacked_dir.joinpath('outputs', name, 'measurement_results.h5')):
        data = get_final_results(name)
        measurement_cells = discretization.get_measurement_cells()
        qz = data['qz']
        conc = data['conc']

        df = pd.DataFrame(index=measurement_cells.files[:int(len(measurement_cells)/2)],
                          columns=["discharge_center", "discharge_min", "discharge_average", "discharge_max",
                                   "conc_center", "conc_min", "conc_average", "conc_max"])

        for center, circle in zip(measurement_cells.files[:int(len(measurement_cells)/2)],
                                measurement_cells.files[int(len(measurement_cells)/2):]):

            conc_center = conc[0][0][measurement_cells[center][2].astype(int)]
            conc_min = np.min(conc[0][0][measurement_cells[circle][2].astype(int)])
            conc_max = np.max(conc[0][0][measurement_cells[circle][2].astype(int)])
            conc_average = np.mean(conc[0][0][measurement_cells[circle][2].astype(int)])

            discharge_center = qz[measurement_cells[center][0].astype(int)]
            discharge_min = np.min(qz[measurement_cells[circle][0].astype(int)])
            discharge_max = np.max(qz[measurement_cells[circle][0].astype(int)])
            discharge_average = np.mean(qz[measurement_cells[circle][0].astype(int)])

            df['conc_center'][center] = conc_center
            df['conc_min'][center] = conc_min
            df['conc_max'][center] = conc_max
            df['conc_average'][center] = conc_average
            df['discharge_center'][center] = discharge_center
            df['discharge_min'][center] = discharge_min
            df['discharge_max'][center] = discharge_max
            df['discharge_average'][center] = discharge_average

        df.to_hdf(unbacked_dir.joinpath('outputs', name, 'measurement_results.h5'), key='df', mode='w')
        return df
    else:
        return pd.read_hdf(unbacked_dir.joinpath('outputs', name, 'measurement_results.h5'))

def get_all_concentrations(name):
    if not os.path.exists(unbacked_dir.joinpath('outputs', name, 'all_results.npz')):
        sim = flopy.mf6.MFSimulation.load(
            sim_ws=unbacked_dir.joinpath('models', name),
        )
        gwt = sim.get_model("trans")
        times = gwt.output.concentration().get_times()
        concs = []
        for time in times:
            conc = gwt.output.concentration().get_data(totim=time)
            concs.append(conc)
        np.savez_compressed(unbacked_dir.joinpath('outputs', name, 'all_results.npz'), concs=concs)
        return concs
    else:
        return np.load(unbacked_dir.joinpath('outputs', name, 'all_results.npz'))['concs']

def get_max_waiwhetu_concentration_over_time(name, name_init, times = None):

    if times is None:
        times = np.arange(0, 31)
    if not os.path.exists(unbacked_dir.joinpath('outputs', name, 'max_waiwhetu_concentration.npz')):
        output = []
        data = get_all_concentrations(name)
        data0 = get_final_results(name_init)['conc']
        waiwhetu = discretization.get_zone_array() == 6

        output.append(np.max(data0[0][0]*waiwhetu))
        for conc in data:
            output.append(np.max(conc[0][0]*waiwhetu))

        np.savez_compressed(unbacked_dir.joinpath('outputs', name, 'max_waiwhetu_concentration.npz'), times=times, max_conc=output)
        return {
            "times": times,
            "max_conc": output
        }
    else:
        return np.load(unbacked_dir.joinpath('outputs', name, 'max_waiwhetu_concentration.npz'))

def calc_pockmark_and_model_areas(buffer=50):
    pockmarks_path = data_dir.joinpath('pockmarks6.shp')
    model_area_path = data_dir.joinpath('model_area0.shp')
    pockmarks = gpd.read_file(pockmarks_path)
    model_area = gpd.read_file(model_area_path)
    # lat lon
    pockmarks = pockmarks.set_crs(epsg=4326)
    pockmarks.to_crs(epsg=2193, inplace=True)
    model_area = model_area.set_crs(epsg=2193)
    pockmark_area = pockmarks.geometry.area.sum()
    model_area_summed = model_area.geometry.area.sum()
    print(f"pockmark area: {pockmark_area} m2")
    print(f"model area: {model_area_summed} m2")
    diff = model_area_summed - pockmark_area
    print(f"non pockmark area: {diff} m2")

    # print lengths of sides of model area
    model_area_exterior = model_area.geometry.exterior[0]
    lengths = []
    for i in range(len(model_area_exterior.xy[0])-1):
        lengths.append(LineString(np.array([model_area_exterior.xy[0][i:i+2], model_area_exterior.xy[1][i:i+2]]).T).length)
    print(f"lengths of sides of model area: {lengths}")

    df = pd.DataFrame(index=pockmarks.index, columns=['area', 'area_with_buffer_diff'])
    for i, row in pockmarks.iterrows():
        df.loc[i, 'area'] = row.geometry.area
        df.loc[i, 'area_with_buffer_diff'] = row.geometry.buffer(buffer).area - row.geometry.area

    df.to_csv(unbacked_dir.joinpath('outputs', f'pockmark_areas_with_{buffer}m_buffer.csv'))

def calc_inner_outer_comparison(name, buffer=50):
    areas = pd.read_csv(unbacked_dir.joinpath('outputs', f'pockmark_areas_with_{buffer}m_buffer.csv'), index_col=0)
    inner_discharge = get_final_pockmark_results(name)
    outer_discharge = get_discharge_around_pockmarks(name, buffer=buffer)
    df = pd.DataFrame(index=areas.index, columns=['inner_discharge', 'outer_discharge', 'inner_conc', 'outer_conc', 'area', 'area_with_buffer_diff',
                                                 'discharge_per_m2_inner', 'discharge_per_m2_outer'])
    for pockmark in areas.index:
        df.loc[pockmark, 'area'] = areas.loc[pockmark, 'area']
        df.loc[pockmark, 'area_with_buffer_diff'] = areas.loc[pockmark, 'area_with_buffer_diff']
        df.loc[pockmark, 'inner_discharge'] = inner_discharge.loc[f'{pockmark}', 'discharge']
        df.loc[pockmark, 'outer_discharge'] = outer_discharge.loc[f'{pockmark}', 'discharge']
        df.loc[pockmark, 'inner_conc'] = inner_discharge.loc[f'{pockmark}', 'concentration']
        df.loc[pockmark, 'outer_conc'] = outer_discharge.loc[f'{pockmark}', 'concentration']
        df.loc[pockmark, 'discharge_per_m2_inner'] = inner_discharge.loc[f'{pockmark}', 'discharge']/areas.loc[pockmark, 'area']
        df.loc[pockmark, 'discharge_per_m2_outer'] = outer_discharge.loc[f'{pockmark}', 'discharge']/areas.loc[pockmark, 'area_with_buffer_diff']
        # add percent difference between inner and outer discharge
        df.loc[pockmark, 'percent_difference_discharge'] = (df.loc[pockmark, 'discharge_per_m2_inner'] - df.loc[pockmark, 'discharge_per_m2_outer'])/df.loc[pockmark, 'discharge_per_m2_outer']*100
        # add percent difference between inner and outer concentration
        df.loc[pockmark, 'percent_difference_conc'] = (df.loc[pockmark, 'inner_conc'] - df.loc[pockmark, 'outer_conc'])/df.loc[pockmark, 'outer_conc']*100
        # also as factor
        df.loc[pockmark, 'factor_difference_discharge'] = df.loc[pockmark, 'discharge_per_m2_inner']/df.loc[pockmark, 'discharge_per_m2_outer']
        df.loc[pockmark, 'factor_difference_conc'] = df.loc[pockmark, 'inner_conc']/df.loc[pockmark, 'outer_conc']
        # also as a factor of reduction in concentration compared to outside
        df.loc[pockmark, 'factor_reduction_conc'] = df.loc[pockmark, 'outer_conc']/df.loc[pockmark, 'inner_conc']
        df.loc[pockmark, 'factor_reduction_discharge'] = df.loc[pockmark, 'outer_discharge']/df.loc[pockmark, 'inner_discharge']

    df.to_csv(unbacked_dir.joinpath('outputs', name, f'pockmark_inner_outer_comparison_{buffer}m_buffer.csv'))
    return df

if __name__=="__main__":
    # get_final_results("disu_test", True)
    # df = get_final_pockmark_results("stdy_thick")
    # get_total_seafloor_discharge("disu_test")
    # df = get_salinity_and_flux_at_measurement_cells("steady_cbc")
    # results = get_pockmark_discharge_and_recharge("dropping_t")
    calc_pockmark_and_model_areas()
    # results = get_discharge_around_pockmarks('1average', buffer=50)
    # calc_inner_outer_comparison('1average', buffer=50)