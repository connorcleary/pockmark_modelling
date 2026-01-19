import numpy as np
from shapely import LineString

import discretization
from discretization import leapfrog_name
import geopandas as gpd
import data_processing
from project_base import project_dir

def get_hydraulic_params(case='base'):
    if case=='base':
        return {
            'hkpu': 0.82,
            'vkpu': 0.006,
            'hkpd': 0.0005,
            'vkpd': 1.8e-5,
            'hkpl': 1.5,
            'vkpl': 0.0094,
            'hkuw': 1000,
            'vkuw': 0.1,
            'dm_coeff': 8.64e-5,
            'sspu': 2e-5,
            'sspd': 3.1e-5,
            'sspl': 5.5e-5,
            'ssuw': 3.1e-5,
            'alpha': 0,
            'hksea': 1000,
            'vksea': 1000,
            'sssea': 1,
        }

def get_boundary_values(time="current", scenario='summer', step=None):

    distance = 3200

    model_area = gpd.read_file(project_dir.joinpath('data', 'model_area0.shp'))
    length = np.max([LineString(np.array([model_area.exterior[0].xy[0][i:i+2], model_area.exterior[0].xy[1][i:i+2]]).T).length
             for i in range(len(model_area.exterior[0].xy[0])-1)])
    print(f"length: {length} m")
    bore_values = data_processing.find_summer_values()

    # dims = discretization.get_dimensions()
    if time=='current':
        seafloor = -0.14
    elif time=='SLR1':
        seafloor = 0.49
    elif time=='SLR2':
        seafloor = 1.45

    if scenario=='summer_average':
        gradient = (bore_values[time]['mcew_median']-bore_values[time]['somes_median'])/distance
        onshore = bore_values[time]['mcew_median']
    elif scenario=='summer_minimum':
        gradient = (bore_values[time]['mcew_min']-bore_values[time]['somes_min'])/distance
        onshore = bore_values[time]['mcew_min']
    elif scenario=='swi':
        gradient = -0.05/distance
        onshore = 1.5
    elif scenario=='swi_b':
        gradient = -0.05/distance
        onshore = -1
    elif scenario=="dropping_head":
        assert step is not None
        onshore = 2.75-0.25*(step+1)
        somes = 0.884*onshore+0.22
        gradient = (onshore-somes)/distance
        return {
            'seafloor': seafloor,
            'onshore': onshore,
            'offshore': onshore-gradient*length,
        }


    return {
        'seafloor': seafloor,
        'onshore': onshore,
        'offshore': onshore-gradient*length,
    }


def  get_tdis_params(scenario='steady'):
    if scenario=='steady':
        return {
            'perlen': 1e6,
            'nstp': 1e3,
            'tsmult': 1,
            'nper': 1,
            'frequency': 100,
        }
    if scenario=='transient':
        return {
            'perlen': 30,
            'nstp': 30,
            'tsmult': 1,
            'nper': 1,
            'frequency': 1,
        }
    if scenario=='slr2':
        return {
            'perlen': 365*105,
            'nstp': 105,
            'tsmult': 1,
            'nper': 1,
            'frequency': 5,
        }
    if scenario=='slr1':
        return {
            'perlen': 365*45,
            'nstp': 45,
            'tsmult': 1,
            'nper': 1,
            'frequency': 5,
        }
    if scenario=='short':
        return {
            'perlen': 365*105,
            'nstp': 105,
            'tsmult': 1,
            'nper': 1,
            'frequency': 5,
        }
    if scenario=='long':
        return {
            'perlen': 1e6,
            'nstp': 1e4,
            'tsmult': 1,
            'nper': 1,
            'frequency': 1e4,
        }
    if scenario=='dropping_head':
        return {
            'perlen': 365,
            'nstp': 73,
            'tsmult': 1,
            'nper': 15,
            'frequency': 10,
        }
    if scenario=="transient_l":
        return {
            'perlen': 3650,
            'nstp': 365,
            'tsmult': 1,
            'nper': 1,
            'frequency': 365,
            }

def get_params_by_zone(case, param, hypothesis='null', patch_multiplier=None, conduit_k=None):
    zones = discretization.get_zone_array()
    params = get_hydraulic_params(case)
    if param == 'hk':
        zone_values = [params['hkpu'], params['hkpd'], params['hkpl'], params['hkuw'], params['hksea']]
    elif param == 'vk':
        zone_values = [params['vkpu'], params['vkpd'], params['vkpl'], params['vkuw'], params['vksea']]
    elif param == 'ss':
        zone_values = [params['sspu'], params['sspd'], params['sspl'], params['ssuw'], params['sssea']]

    if leapfrog_name is None:
        zones[zones == -1] = zone_values[0]
        zones[zones == 1] = zone_values[0]
        zones[zones == 2] = zone_values[1]
        zones[zones == 3] = zone_values[2]
        zones[zones == 4] = zone_values[3]
    else:
        print("make sure that you check the leapfrog zone values")
        zones[zones == 4] = zone_values[4]
        zones[zones == 1] = zone_values[1]
        zones[zones == 2] = zone_values[2]
        zones[zones == 3] = zone_values[0]
        zones[zones == 6] = zone_values[3]

    if hypothesis == 'null':
        return zones

    elif hypothesis == 'high_pockmark_k':
        if patch_multiplier is None:
            patch_multiplier = 2
        distal_below_pockmarks = discretization.get_distal_below_pockmark_cells()
        if param == 'hk' or param == 'vk':
            zones[distal_below_pockmarks] = zones[distal_below_pockmarks]*patch_multiplier
        return zones

    elif hypothesis == "conduits":
        if conduit_k is None:
            conduit_k = 10
        distal_and_lower_below_conduits = discretization.get_distal_and_lower_below_conduit_cells()
        if param == 'hk' or param == 'vk':
            zones[distal_and_lower_below_conduits] = conduit_k
        return zones

if __name__ == '__main__':
    # print(get_hydraulic_params())dm    get_boundary_values()
    print(get_boundary_values(scenario='summer'))
    print(get_boundary_values(scenario='slr'))
    print(get_boundary_values(scenario='low_summer'))
    print(get_boundary_values(scenario='low_slr'))
    print(get_boundary_values(scenario='swi_summer'))
    print(get_boundary_values(scenario='swi_slr'))

