import os.path
import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import shapely
# from networkx.classes import edges
from numpy.f2py.crackfortran import true_intent_list
from rasterio import ubyte
from shapely.geometry import (
    LineString,
    MultiLineString,
    MultiPoint,
    Point,
    Polygon,
)
import rasterio
from rasterio.plot import show
import flopy
import flopy.discretization as fgrid
import flopy.plot as fplot
from flopy.utils import GridIntersect
import pyproj
import geopandas as gpd
from flopy.utils import Raster
from pathlib import Path
import pickle
import pandas as pd
from pyproj import Transformer
from project_base import unbacked_dir, project_dir

# delc = 5
# delr = 5
# delv = 0.5
# top = 0.0
# bottom = -56.0
# nlay = -int((bottom-top)/delv)
# bottom = (np.arange(0, nlay)+1) * -delv
# width = 1795
# length = 2710
# angle = np.deg2rad(-22)
# nrows = int(width/delc)
# ncols = int(length/delr)

# delr_array = np.ones(ncols) * delr
# delc_array = np.ones(nrows) * delc
# top_array = np.ones((nrows, ncols)) * top
# botm = np.repeat(bottom[:, np.newaxis], nrows, axis=1)
# botm = np.repeat(botm[:, :, np.newaxis], ncols, axis=2)

# leapfrog_name =  "MODFLOW_Simulati" # TEST
# leapfrog_name = "thicker_seafloor" # NEW
leapfrog_name = "grid2"

# def get_dimensions():
#     return {
#         "width": width,
#         "length": length,
#     }

def get_grid(compact_top=0.1):
    if leapfrog_name is None:
        grid = fgrid.StructuredGrid(
            delc_array, delr_array, top=top_array, botm=botm, crs=2193, xoff=1758230, yoff=5431297, angrot=68
        )
        return grid

    grid = fgrid.UnstructuredGrid.from_gridspec(project_dir.joinpath("data", "leapfrog", leapfrog_name, f"{leapfrog_name}.gsf"))
    if compact_top is not None:
        zones = get_zone_array()
        seafloor = np.where(zones == 4)[0]
        grid.top[seafloor] = grid.botm[seafloor] + compact_top
        grid._top[seafloor] = grid.botm[seafloor] + compact_top
        grid.cell_thickness[seafloor] = compact_top
        grid.top_botm[0, seafloor] = grid.top[seafloor]
        grid.xyzcellcenters[2][0][seafloor] = grid.botm[seafloor] + compact_top/2
        grid.xyzvertices[2][0][seafloor] = grid.botm[seafloor] + compact_top
        grid.zcellcenters[0][seafloor] = grid.botm[seafloor] + compact_top/2
        grid.zvertices[0, seafloor] = grid.botm[seafloor] + compact_top

    return grid


def get_ibound():
    if leapfrog_name is None:
        dis_dir = project_dir.joinpath('data', 'discretization')
        dis_dir.mkdir(exist_ok=True)
        dis_name = f"nlay={nlay}_nrow={nrows}_ncol={ncols}"
        dis_dir_name = dis_dir.joinpath(dis_name)
        dis_dir_name.mkdir(exist_ok=True)
        if not dis_dir_name.joinpath('ibound.npy').exists():
            grid = get_grid()
            zcenters = grid.xyzcellcenters[2]

            ix = GridIntersect(grid, method="vertex")
            offshore = gpd.read_file("/data/offshore_area.shp")
            offshore_area = ix.intersect(offshore.values[0,1])

            bathymetry = Raster.load("/data/cropped_bath.tif")
            bath_data = bathymetry.resample_to_grid(
                grid,  band=bathymetry.bands[0], method="cubic"
            )
            base_waiwhetu = Raster.load("/data/Mid Waiwhetu - Upper Waiwhetu contacts.asc")
            base_waiwhetu_data = base_waiwhetu.resample_to_grid(
                grid, band=base_waiwhetu.bands[0], method="cubic",
            )
            ibound = np.ones((nlay, nrows, ncols))
            ibound[zcenters[:, :, :] < base_waiwhetu_data] = 0

            offshore = np.ones((nrows, ncols))*np.nan
            for cell in offshore_area.cellids:
                offshore[cell[0], cell[1]] = 1

            ibound[zcenters[:, :, :] > np.multiply(bath_data, offshore)] = 0
            np.save(dis_dir_name.joinpath('ibound.npy'), ibound)
        else:
            ibound = np.load(dis_dir_name.joinpath('ibound.npy'))
    else:
        ibound = 1

    return ibound

def get_cell_connectivity(name):
    if not os.path.exists(project_dir.joinpath('data', 'discretization', f'leapfrog_{leapfrog_name}', 'cell_connectivity.csv')):
        temp_sim = flopy.mf6.MFSimulation.load(sim_name='mfsim',
                                               sim_ws=unbacked_dir.joinpath('models', name))
        temp_model = temp_sim.get_model(name)
        temp_disu = temp_model.get_package("disu")
        ja = temp_disu.ja.array
        iac = temp_disu.iac.array
        cell_connectivity = []
        for cell in iac:
            cell_connectivity.append(ja[:cell])
            ja = ja[cell:]

        np.savetxt(project_dir.joinpath('data', 'discretization', f'leapfrog_{leapfrog_name}', 'cell_connectivity.csv'), cell_connectivity)
        return cell_connectivity
    return np.load(project_dir.joinpath('data', 'discretization', f'leapfrog_{leapfrog_name}', 'cell_connectivity.csv'))

def get_offshore():
    dis_dir = project_dir.joinpath('data', 'discretization')
    dis_dir.mkdir(exist_ok=True)
    dis_name = f"nlay={nlay}_nrow={nrows}_ncol={ncols}"
    dis_dir_name = dis_dir.joinpath(dis_name)
    dis_dir_name.mkdir(exist_ok=True)
    if not dis_dir_name.joinpath('offshore_area.npy').exists():
        grid = get_grid()
        ix = GridIntersect(grid, method="vertex")
        offshore = gpd.read_file("/data/gis/offshore_area.shp")
        offshore_area = ix.intersect(offshore.values[0,1])
        offshore = np.ones((nrows, ncols)) * np.nan
        for cell in offshore_area.cellids:
            offshore[cell[0], cell[1]] = 1
        np.save(dis_dir_name.joinpath('offshore_area.npy'), offshore)
    else:
        offshore = np.load(dis_dir_name.joinpath('offshore_area.npy'))

    return offshore

def get_grid_w_ibound() -> object:
    grid = get_grid()
    ibound = get_ibound()
    grid.idomain = ibound
    return grid


def get_zone_array():

    if leapfrog_name is None:
        dis_dir = project_dir.joinpath('data', 'discretization')
        dis_dir.mkdir(exist_ok=True)
        dis_name = f"nlay={nlay}_nrow={nrows}_ncol={ncols}"
        dis_dir_name = dis_dir.joinpath(dis_name)
        dis_dir_name.mkdir(exist_ok=True)
        if not dis_dir_name.joinpath('zone_array.npy').exists():
            grid = get_grid_w_ibound()
            base_petone_l = Raster.load("/data/gis/Upper Waiwhetu - Petone Marine lower contacts.asc")
            base_petone_m = Raster.load("/data/gis/Petone Marine lower - Petone Marine distal contacts.asc")
            base_petone_u = Raster.load("/data/gis/Petone Marine distal - Petone Marine upper contacts.asc")

            base_petone_l_data = base_petone_l.resample_to_grid(
                grid, band=base_petone_l.bands[0], method="cubic"
            )
            base_petone_m_data = base_petone_m.resample_to_grid(
                grid, band=base_petone_m.bands[0], method="cubic"
            )
            base_petone_u_data = base_petone_u.resample_to_grid(
                grid, band=base_petone_u.bands[0], method="cubic"
            )

            zcenters = grid.xyzcellcenters[2]
            zone_array = grid.idomain.copy()
            zone_array[zcenters > base_petone_u_data] = 1
            zone_array[zcenters < base_petone_u_data] = 2
            zone_array[zcenters < base_petone_m_data]  = 3
            zone_array[zcenters < base_petone_l_data]  = 4
            for row in range(grid.nrow):
                for col in range(grid.ncol):
                    zone_array[np.argmax(grid.idomain[:, row, col] > 0), row, col] = -1

            zone_array = np.multiply(zone_array, grid.idomain)
            np.save(dis_dir_name.joinpath('zone_array.npy'), zone_array)
        else:
            zone_array = np.load(dis_dir_name.joinpath('zone_array.npy'))

        return zone_array
    else:
        arr = np.genfromtxt(project_dir.joinpath("data", "leapfrog", leapfrog_name, f"{leapfrog_name}.mfs"), skip_header=3)
        return arr

def get_boundary_lines():

    model_area = gpd.read_file(project_dir.joinpath('data', 'model_area0.shp'))
    model_area = model_area.to_crs(epsg=2193)
    coord_list = model_area['geometry'][0].exterior.coords.xy
    y, order = np.unique(coord_list[1], return_index=True)
    x = [coord_list[0][i] for i in order]
    assert len(x) == len(y)
    assert len(x) == 4
    y_order = np.argsort(y)
    offshore_line = LineString([(x[y_order[0]], y[y_order[0]]+1), (x[y_order[1]], y[y_order[1]]+1)])
    onshore_line = LineString([(x[y_order[2]], y[y_order[2]]-1), (x[y_order[3]], y[y_order[3]]-1)])
    assert not offshore_line.intersects(onshore_line)

    return offshore_line, onshore_line

def get_intersected_cells(shapes, xvertices, yvertices, cell_indices, names = None, area = True):

    if names == None:
        names = range(len(shapes))
    data = {f'{names[i]}': [[], [], []] for i in range(len(shapes['geometry']))}

    for i, cell in enumerate(list(zip(xvertices, yvertices))):
        cell = shapely.convex_hull(Polygon(list(zip(cell[0], cell[1]))))
        for j, (shape, name) in enumerate(zip(shapes['geometry'], names)):
            if shapely.intersects(cell, shape):
                data[f'{name}'][0].append(cell_indices[i]) # seafloor cell
                data[f'{name}'][1].append(cell.area) # seafloor cell area
                data[f'{name}'][2].append(cell_indices[i]+len(cell_indices)) # cell below seafloor
            elif shape.within(cell):
                data[f'{name}'][0].append(cell_indices[i])  # seafloor cell
                data[f'{name}'][1].append(cell.area)  # seafloor cell area
                data[f'{name}'][2].append(cell_indices[i] + len(cell_indices))  # cell below seafloor


    return data

def get_boundary_intersection(grid, line):
    intersected_cells = []
    for i, cell in enumerate(list(zip(grid.xvertices, grid.yvertices))):
        cell = shapely.convex_hull(Polygon(list(zip(cell[0], cell[1]))))
        if line.intersects(cell):
            intersected_cells.append(i)

    assert len(intersected_cells) <= len(grid.xvertices)
    assert len(intersected_cells) == len(np.unique(intersected_cells))
    assert len(intersected_cells) > 0
    return intersected_cells


def get_boundary_arrays():

    if leapfrog_name is None:

        grid = get_grid_w_ibound()
        base = grid.idomain.copy()
        zone_array = get_zone_array()

        dis_dir = project_dir.joinpath('data', 'discretization')
        dis_dir.mkdir(exist_ok=True)
        dis_name = f"nlay={nlay}_nrow={nrows}_ncol={ncols}"
        dis_dir_name = dis_dir.joinpath(dis_name)
        dis_dir_name.mkdir(exist_ok=True)

        if not dis_dir_name.joinpath('onshore.npy').exists():
            onshore = np.zeros_like(base).astype(bool)
            onshore[zone_array == 4] = True
            onshore[:, :, :-1] = False
            np.save(dis_dir_name.joinpath('onshore.npy'), onshore)
        else:
            onshore = np.load(dis_dir_name.joinpath('onshore.npy'))

        if not dis_dir_name.joinpath('surface.npy').exists():
            surface = np.zeros_like(base).astype(bool)
            offshore = get_offshore()
            for row in range(grid.nrow):
                for col in range(grid.ncol):
                    if offshore[row, col] == 1.0:
                        surface[np.argmax(zone_array[:, row, col]==-1, axis=0), row, col] = True
            np.save(dis_dir_name.joinpath('surface.npy'), surface)
        else:
            surface = np.load(dis_dir_name.joinpath('surface.npy'))

        if not dis_dir_name.joinpath('offshore.npy').exists():
            offshore = np.zeros_like(base).astype(bool)
            offshore[zone_array == 4] = True
            offshore[:, :, 1:] = False
            np.save(dis_dir_name.joinpath('offshore.npy'), offshore)
        else:
            offshore = np.load(dis_dir_name.joinpath('offshore.npy'))
    else:
        print("make sure that you check the leapfrog zone values")
        dis_dir = project_dir.joinpath('data', 'discretization')
        dis_dir.mkdir(exist_ok=True)
        dis_name = f"leapfrog_{leapfrog_name}"
        dis_dir_name = dis_dir.joinpath(dis_name)
        dis_dir_name.mkdir(exist_ok=True)

        offshore_line, onshore_line = get_boundary_lines()

        grid = get_grid()

        if not dis_dir_name.joinpath('onshore.npy').exists():
            intersected_cells = get_boundary_intersection(grid, onshore_line)
            zone_array = get_zone_array()
            waiwhetu = np.where(zone_array == 6)
            onshore = np.intersect1d(intersected_cells , waiwhetu)
            np.save(dis_dir_name.joinpath('onshore.npy'), onshore)
        else:
            onshore = np.load(dis_dir_name.joinpath('onshore.npy'))

        if not dis_dir_name.joinpath('offshore.npy').exists():
            intersected_cells = get_boundary_intersection(grid, offshore_line)
            zone_array = get_zone_array()
            waiwhetu = np.where(zone_array == 6)
            offshore = np.intersect1d(intersected_cells, waiwhetu)
            np.save(dis_dir_name.joinpath('offshore.npy'), offshore)
        else:
            offshore = np.load(dis_dir_name.joinpath('offshore.npy'))

        if not dis_dir_name.joinpath('surface.npy').exists():
            zone_array = get_zone_array()
            surface = np.where(zone_array == 4)
            np.save(dis_dir_name.joinpath('surface.npy'), surface)
        else:
            surface = np.load(dis_dir_name.joinpath('surface.npy'))

    return {
            'onshore': onshore,
            'seafloor': surface,
            'offshore': offshore
            }

def get_start_conc(salinity=35.0):
    boundaries = get_boundary_arrays()
    seafloor = boundaries['seafloor']
    start = np.zeros_like(seafloor)
    start[seafloor] = salinity
    return start

def get_cells_below_seafloor(wrong_formulation=False):

    dis_dir = project_dir.joinpath('data', 'discretization')
    dis_dir.mkdir(exist_ok=True)
    dis_name = f"nlay={nlay}_nrow={nrows}_ncol={ncols}"
    dis_dir_name = dis_dir.joinpath(dis_name)
    dis_dir_name.mkdir(exist_ok=True)

    grid = get_grid_w_ibound()
    base = grid.idomain.copy()
    zone_array = get_zone_array()

    if not dis_dir_name.joinpath('below_seafloor.npy').exists():
        below_seafloor = np.zeros_like(base).astype(bool)
        offshore = get_offshore()
        for row in range(grid.nrow):
            for col in range(grid.ncol):
                if offshore[row, col] == 1.0:
                    if wrong_formulation:
                        below_seafloor[np.argmax(grid.idomain[:, row, col] > 0, axis=0)+1, row, col] = True
                    else:
                        below_seafloor[np.argmax(zone_array[:, row, col] > 0, axis=0), row, col] = True

        np.save(dis_dir_name.joinpath('below_seafloor.npy'), below_seafloor)
    else:
        below_seafloor = np.load(dis_dir_name.joinpath('below_seafloor.npy'))

    return below_seafloor

def get_pockmark_surfaces():
    """

    :return:
    dictionary with entries for each pockmark
    first list is the cells below seafloor
    second list is the area of the cells
    third list is the seafloor cells
    """

    if not os.path.exists(project_dir.joinpath('data', 'discretization', f'leapfrog_{leapfrog_name}', 'pockmarks_surfaces.npz')):
        pockmarks = gpd.read_file(project_dir.joinpath('data', 'pockmarks6_new.shp'))
        grid = get_grid()
        # zone_array = get_zone_array()
        assert len(np.unique(grid.ncpl)) == 1
        seafloor = np.arange(0, grid.ncpl[0])
        xvertices_temp = grid.xvertices
        yvertices_temp = grid.yvertices
        xvertices = [xvertices_temp[cell] for cell in seafloor]
        yvertices = [yvertices_temp[cell] for cell in seafloor]

        pockmark_surfaces = get_intersected_cells(pockmarks, xvertices, yvertices, seafloor)

        np.savez_compressed(project_dir.joinpath('data', 'discretization', f'leapfrog_{leapfrog_name}', 'pockmarks_surfaces.npz'), **pockmark_surfaces)
        return pockmark_surfaces

    return np.load(project_dir.joinpath('data', 'discretization', f'leapfrog_{leapfrog_name}', 'pockmarks_surfaces.npz'))

def get_surface_around_pockmarks(buffer=100):
    if not os.path.exists(project_dir.joinpath('data', 'discretization', f'leapfrog_{leapfrog_name}', f'around_pockmarks_surfaces_{buffer}.npz')):
        pockmarks = gpd.read_file(project_dir.joinpath('data', 'pockmarks6_new.shp'))
        around_pockmarks = pockmarks.copy()
        around_pockmarks['geometry'] = pockmarks.geometry.buffer(buffer)
        # remove original pockmark shapes
        around_pockmarks['geometry'] = around_pockmarks.difference(pockmarks.union_all())

        grid = get_grid()
        # zone_array = get_zone_array()
        assert len(np.unique(grid.ncpl)) == 1
        seafloor = np.arange(0, grid.ncpl[0])
        xvertices_temp = grid.xvertices
        yvertices_temp = grid.yvertices
        xvertices = [xvertices_temp[cell] for cell in seafloor]
        yvertices = [yvertices_temp[cell] for cell in seafloor]

        around_pockmark_surfaces = get_intersected_cells(around_pockmarks, xvertices, yvertices, seafloor)

        np.savez_compressed(project_dir.joinpath('data', 'discretization', f'leapfrog_{leapfrog_name}',  f'around_pockmarks_surfaces_{buffer}.npz'), **around_pockmark_surfaces)
        return around_pockmark_surfaces

    return np.load(project_dir.joinpath('data', 'discretization', f'leapfrog_{leapfrog_name}',  f'around_pockmarks_surfaces_{buffer}.npz'))

def get_pockmark_connectivity():

    if not os.path.exists(project_dir.joinpath('data', 'discretization', f'leapfrog_{leapfrog_name}', 'pockmarks_connectivity.npz')):

        pockmarks_cells = get_pockmark_surfaces()
        temp_sim = flopy.mf6.MFSimulation.load(sim_name='mfsim',
                                               sim_ws=project_dir.joinpath('data', 'leapfrog', leapfrog_name))
        temp_model = temp_sim.get_model(leapfrog_name)
        temp_disu = temp_model.get_package("disu")
        ja = temp_disu.ja.array
        iac = temp_disu.iac.array
        cell_connectivity = []
        cell_connection_indices = []

        zone_array = get_zone_array()
        # todo only get connection to aquitard cells
        cell_connection_indices.append([i for i in range(iac[0])])
        cell_connectivity.append(ja[:iac[0]])
        ja = ja[iac[0]:]
        for cell in iac[1:]:
            cell_connectivity.append(ja[:cell])
            cell_connection_indices.append([i for i in range(cell_connection_indices[-1][-1]+1, cell_connection_indices[-1][-1]+cell+1)])
            ja = ja[cell:]

        pockmarks_connectivity = {}
        for pockmark in pockmarks_cells.keys():
            pockmarks_connectivity[pockmark] = []
            for cell in np.array(pockmarks_cells[pockmark][0]).astype(int):
                for connection, index in zip(cell_connectivity[cell], cell_connection_indices[cell]):
                    if zone_array[connection] != 4:
                        pockmarks_connectivity[pockmark].append([connection, index])

        np.savez_compressed(project_dir.joinpath('data', 'discretization', f'leapfrog_{leapfrog_name}', 'pockmarks_connectivity.npz'),
            **pockmarks_connectivity)
        return pockmarks_connectivity

    return np.load(project_dir.joinpath('data', 'discretization', f'leapfrog_{leapfrog_name}', 'pockmarks_connectivity.npz'))

def get_around_pockmark_connectivity(buffer=100):

    if not os.path.exists(project_dir.joinpath('data', 'discretization', f'leapfrog_{leapfrog_name}', f'around_pockmarks_connectivity_{buffer}.npz')):

        pockmarks_cells = get_surface_around_pockmarks(buffer=buffer)
        temp_sim = flopy.mf6.MFSimulation.load(sim_name='mfsim',
                                               sim_ws=project_dir.joinpath('data', 'leapfrog', leapfrog_name))
        temp_model = temp_sim.get_model(leapfrog_name)
        temp_disu = temp_model.get_package("disu")
        ja = temp_disu.ja.array
        iac = temp_disu.iac.array
        cell_connectivity = []
        cell_connection_indices = []

        zone_array = get_zone_array()
        # todo only get connection to aquitard cells
        cell_connection_indices.append([i for i in range(iac[0])])
        cell_connectivity.append(ja[:iac[0]])
        ja = ja[iac[0]:]
        for cell in iac[1:]:
            cell_connectivity.append(ja[:cell])
            cell_connection_indices.append([i for i in range(cell_connection_indices[-1][-1]+1, cell_connection_indices[-1][-1]+cell+1)])
            ja = ja[cell:]

        pockmarks_connectivity = {}
        for pockmark in pockmarks_cells.keys():
            pockmarks_connectivity[pockmark] = []
            for cell in np.array(pockmarks_cells[pockmark][0]).astype(int):
                for connection, index in zip(cell_connectivity[cell], cell_connection_indices[cell]):
                    if zone_array[connection] != 4:
                        pockmarks_connectivity[pockmark].append([connection, index])

        np.savez_compressed(project_dir.joinpath('data', 'discretization', f'leapfrog_{leapfrog_name}', f'around_pockmarks_connectivity_{buffer}.npz'),
            **pockmarks_connectivity)
        return pockmarks_connectivity

    return np.load(project_dir.joinpath('data', 'discretization', f'leapfrog_{leapfrog_name}', f'around_pockmarks_connectivity_{buffer}.npz'))


def get_non_pockmark_surfaces():

    if not os.path.exists(project_dir.joinpath('data', 'discretization', f'leapfrog_{leapfrog_name}', 'non_pockmark_surfaces.npz')):
        pockmarks = get_pockmark_surfaces()
        grid = get_grid()
        zone_array = get_zone_array()
        seafloor = np.where(zone_array == 4)
        below_seafloor = seafloor + grid.ncpl[0]
        below_seafloor_petone_distal = below_seafloor[np.where(zone_array[below_seafloor.astype(int)] == 1)]
        pockmark_cells = []
        for pockmark in pockmarks.keys():
            pockmark_cells.extend(pockmarks[pockmark][2])
        non_pockmark_petone_distal_seafloor_cells = np.setdiff1d(below_seafloor_petone_distal, pockmark_cells)
        areas = []
        xvertices_temp = grid.xvertices
        yvertices_temp = grid.yvertices

        for cell in non_pockmark_petone_distal_seafloor_cells:
            areas.append(Polygon(list(zip(xvertices_temp[cell], yvertices_temp[cell]))).area)

        np.savez_compressed(project_dir.joinpath('data', 'discretization', f'leapfrog_{leapfrog_name}', 'non_pockmark_surfaces.npz'),
                            cells=non_pockmark_petone_distal_seafloor_cells, areas=areas)

        return {'cells': non_pockmark_petone_distal_seafloor_cells, 'areas': areas}

    return np.load(project_dir.joinpath('data', 'discretization', f'leapfrog_{leapfrog_name}', 'non_pockmark_surfaces.npz'))

def get_measurement_cells(plot=False):
    if not os.path.exists(project_dir.joinpath('data', 'discretization', f'leapfrog_{leapfrog_name}', 'measurement_cells.npz')):

        points = pd.read_csv(project_dir.joinpath('data', 'coordinates_from_harding_appendix.csv'))
        points.set_index('name', inplace=True)
        points.drop('sv4', inplace=True)
        grid = get_grid()
        # zone_array = get_zone_array()
        assert len(np.unique(grid.ncpl)) == 1
        seafloor = np.arange(0, grid.ncpl[0])
        xvertices_temp = grid.xvertices
        yvertices_temp = grid.yvertices
        xvertices = [xvertices_temp[cell] for cell in seafloor]
        yvertices = [yvertices_temp[cell] for cell in seafloor]

        names = [point for point in points.index]
        shapes = {'geometry': []}

        for i in range(len(names)):
            names.append(f"{names[i]}_circle")
        for point, data in points.iterrows():
            # todo add offset to harding points
            transformer = Transformer.from_crs("EPSG:4326", "EPSG:2193")
            if 'sv' in point:
                data.array[0] += 0.0016666
            nztm = transformer.transform(*data.array)
            shapes['geometry'].append(Point(nztm[1], nztm[0]))
        for i in range(len(points)):
            shapes['geometry'].append(shapes['geometry'][i].buffer(5))

        d = {'names': names, 'geometry': shapes['geometry']}
        gdf = gpd.GeoDataFrame(d, crs="EPSG:2197")

        cells = get_intersected_cells(gdf, xvertices, yvertices, seafloor, names)

        if plot:
            f, ax = plt.subplots()
            gdf.plot(ax=ax, color='none', edgecolors='red')
            thickness = rasterio.open(project_dir.joinpath('data', 'marine_thickness_masked_model.tif'))
            thck = show(thickness, ax=ax, cmap='grey', alpha=0.7)
            plt.show()

        np.savez_compressed(project_dir.joinpath('data', 'discretization', f'leapfrog_{leapfrog_name}', 'measurement_cells.npz'), **cells)
        return cells

    return np.load(project_dir.joinpath('data', 'discretization', f'leapfrog_{leapfrog_name}', 'measurement_cells.npz'))

def get_distal_below_pockmark_cells(pockmarks="refinement_areas"):
    if not os.path.exists(project_dir.joinpath('data', 'discretization', f'leapfrog_{leapfrog_name}', 'distal_below_pockmark_cells.npy')):
        if pockmarks is None:
            pockmarks = gpd.read_file(project_dir.joinpath('data', 'pockmarks.shp'))
        elif pockmarks == "refinement_areas":
            # include areas outside of pockmarks where I have refined. i.e. the smaller pockmark and area outside pockmark
            # with high salinity
            pockmarks = gpd.read_file(project_dir.joinpath('data', 'refinement_areas.shp'))

        names = range(len(pockmarks))
        data = []
        grid = get_grid()
        zone_array = get_zone_array()
        xvertices = grid.xvertices
        yvertices = grid.yvertices
        petone_distal = np.where(zone_array == 1)[0]

        assert len(np.unique(grid.ncpl)) == 1
        petone_distal_loc = petone_distal % grid.ncpl[0]
        petone_distal_unique, unique_index, unique_inverse = np.unique(petone_distal_loc, return_index=True, return_inverse=True)
        unique = petone_distal[unique_index]

        xvertices = [xvertices[cell] for cell in unique]
        yvertices = [yvertices[cell] for cell in unique]

        union = shapely.ops.unary_union([shape for shape in pockmarks['geometry']])

        for i, cell in enumerate(list(zip(xvertices, yvertices))):
            cell = shapely.convex_hull(Polygon(list(zip(cell[0], cell[1]))))
            if shapely.intersects(cell, union):
                data.extend(petone_distal[unique_inverse == i])

        np.save(project_dir.joinpath('data', 'discretization', f'leapfrog_{leapfrog_name}', 'distal_below_pockmark_cells.npy'), data)
        return data
    return np.load(project_dir.joinpath('data', 'discretization', f'leapfrog_{leapfrog_name}', 'distal_below_pockmark_cells.npy'))

def get_distal_and_lower_below_conduit_cells(non_measurement=True, all_measurements=True):

    if not os.path.exists(project_dir.joinpath('data', 'discretization', f'leapfrog_{leapfrog_name}', 'distal_and_lower_below_conduit_cells.npy')):
        measurement_cells = get_measurement_cells()
        known_cells = []

        known_cells.append(measurement_cells['sv1'][2])
        known_cells.append(measurement_cells['sv6'][2])
        known_cells.append(measurement_cells['mc22'][2])
        if all_measurements:
            # include places where there are not conduits
            known_cells.append(measurement_cells['mc21'][2])
            known_cells.append(measurement_cells['mc23'][2])


        if non_measurement:
            grid = get_grid()
            assert len(np.unique(grid.ncpl)) == 1
            seafloor = np.arange(0, grid.ncpl[0])
            xvertices = grid.xvertices
            yvertices = grid.yvertices
            xvertices = [xvertices[cell] for cell in seafloor]
            yvertices = [yvertices[cell] for cell in seafloor]
            shapes = gpd.read_file(project_dir.joinpath('data', 'abritrary_pockmark_locations.shp'))
            data = get_intersected_cells(shapes, xvertices, yvertices, seafloor, names=['edge_pockmark_2', 'sv4'])
            known_cells.append(data['edge_pockmark_2'][2])
            known_cells.append(data['sv4'][2])

        cells = []
        zone_array = get_zone_array()
        grid = get_grid()
        ncpl = grid.ncpl[0]
        petone_distal_lower = np.concatenate([np.where(zone_array == 1)[0], np.where(zone_array == 2)[0]])
        for cell in known_cells:
            cells.extend(petone_distal_lower[petone_distal_lower % ncpl == cell % ncpl])

        np.save(project_dir.joinpath('data', 'discretization', f'leapfrog_{leapfrog_name}', 'distal_and_lower_below_conduit_cells.npy'), cells)
        return cells

    cells = np.load(project_dir.joinpath('data', 'discretization', f'leapfrog_{leapfrog_name}', 'distal_and_lower_below_conduit_cells.npy'))
    return cells
# ix.plot_polygon(offshore_area, ax=ax, edgecolor="red")

# grid = get_grid_w_ibound()
# zone_array = get_zone_array()
# pmv = flopy.plot.PlotMapView(modelgrid=grid, layer=18)
# pmv.plot_array(zone_array[12])
# plt.show()

if __name__=="__main__":
    # grid = get_grid(compact_top=True)
    # get_pockmark_connectivity()
    # get_pockmark_surfaces()
    # cells = get_measurement_cells(True)
    # get_distal_below_pockmark_cells()
    get_distal_and_lower_below_conduit_cells()

