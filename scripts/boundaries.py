import discretization
import numpy as np
import parameterization
from discretization import leapfrog_name

def get_ghb_data(boundary, case, time, scenario, nper=1):

    ghb_data = {}
    hk = parameterization.get_params_by_zone(case, 'hk')
    boundary_arrays = discretization.get_boundary_arrays()

    if nper==1:
        boundary_values = parameterization.get_boundary_values(time, scenario)
        ghb_data_temp = []
        if leapfrog_name is None:
            cells = np.argwhere(boundary_arrays[boundary])
            for cell in cells:
                cell_hk = hk[cell[0], cell[1], cell[2]]
                ghb_data_temp.append((cell[0], cell[1], cell[2], boundary_values[boundary], cell_hk, 0.0))
        else:
            for cell in boundary_arrays[boundary]:
                cell_hk = hk[cell]
                ghb_data_temp.append((cell, boundary_values[boundary], cell_hk, 0.0))
        ghb_data[0] = ghb_data_temp
    else:
        for step in range(nper):
            boundary_values = parameterization.get_boundary_values(time, scenario, step)
            ghb_data_temp = []
            for cell in boundary_arrays[boundary]:
                cell_hk = hk[cell]
                ghb_data_temp.append((cell, boundary_values[boundary], cell_hk, 0.0))
            ghb_data[step] = ghb_data_temp

    return ghb_data

def get_chd_data(boundary, time, scenario, salinity=35.0, nper=1):

    boundary_arrays = discretization.get_boundary_arrays()
    chd_data = {}

    if nper==1:
        boundary_values = parameterization.get_boundary_values(time, scenario)
        chd_data_temp = []
        if leapfrog_name is None:
            for cell in cells:
                cells = np.argwhere(boundary_arrays[boundary])
                chd_data_temp.append((cell[0], cell[1], cell[2], boundary_values[boundary], salinity))
        else:
            for cell in boundary_arrays[boundary][0]:
                chd_data_temp.append((cell, boundary_values[boundary], salinity))
        chd_data[0] = chd_data_temp
    else:
        for step in range(nper):
            boundary_values = parameterization.get_boundary_values(time, scenario, step)
            chd_data_temp = []
            for cell in boundary_arrays[boundary][0]:
                chd_data_temp.append((cell, boundary_values[boundary], salinity))
            chd_data[step] = chd_data_temp

    return chd_data

def get_cnc_data(boundary, salinity=35.0, nper=1):
    boundary_arrays = discretization.get_boundary_arrays()
    cnc_data = {}

    for step in range(nper):
        cnc_data_temp = []
        if leapfrog_name is None:
            cells = np.argwhere(boundary_arrays[boundary])
            for cell in cells:
                cnc_data_temp.append((cell[0], cell[1], cell[2], salinity))
        else:
            for cell in boundary_arrays[boundary][0]:
                cnc_data_temp.append((cell, salinity))
        cnc_data[step] = cnc_data_temp

    # if leapfrog_name is None:
    #     cells = np.argwhere(boundary_arrays[boundary])
    #     for cell in cells:
    #         cnc_data.append((cell[0], cell[1], cell[2], salinity))
    # else:
    #     for cell in boundary_arrays[boundary][0]:
    #         cnc_data.append((cell, salinity))

    return cnc_data