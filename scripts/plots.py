# from vtkmodules.numpy_interface.algorithms import square, cross

import discretization
from matplotlib import pyplot as plt
from matplotlib.colors import LightSource, to_rgb
import matplotlib as mpl
import numpy as np
import flopy
import matplotlib.cm as cm
import matplotlib.colors as colors
import parameterization
import os
#import pyvista as pv
from flopy.export.vtk import Vtk
import rasterio
from rasterio.plot import show
from pyproj import Transformer
from discretization import leapfrog_name
import geopandas as gpd
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
from shapely.geometry import Polygon, Point
import shapely
from matplotlib.colors import LinearSegmentedColormap
import utils
import pandas as pd
from matplotlib.lines import Line2D
from project_base import project_dir, unbacked_dir

def plot_unstructured_vectors(ax, cross_section, qx, qy, qz, section, hstride, vstride, zone_array, reverse=False):
    # order x_centers
    plot_elements = []
    assert len(np.unique(cross_section.mg.ncpl)) == 1

    pd = np.where(zone_array == 1)
    ids = np.array(list(cross_section.polygons.keys()))
    layer_ids = np.intersect1d(pd, ids)
    layer_ids = layer_ids[layer_ids//cross_section.mg.ncpl[0]%vstride==0]
    plot_ids = layer_ids[::hstride]

    qz = qz[plot_ids]
    qx = qx[plot_ids]
    qy = qy[plot_ids]

    # get list of distances along
    distances = []
    for i in range(1, len(section.coords.xy[0])):
        distances.append(section.line_locate_point(Point(section.coords.xy[0][i], section.coords.xy[1][i])))

    plot_is = [np.argwhere(ids == id)[0][0] for id in plot_ids]

    qh = []
    for count, (i, id) in enumerate(zip(plot_is, plot_ids)):
        # find segment
        iseg = np.argwhere(cross_section.xcenters[i] < distances)[-1]
        # project qh to the line
        center_coords = cross_section.xypts[id][0]
        seg_end = (section.coords.xy[0][int(iseg+1)], section.coords.xy[1][int(iseg+1)])
        seg_vector = (seg_end[0]-center_coords[0], seg_end[1]-center_coords[1])
        qh_temp = np.dot((qx[count], qy[count]), seg_vector)/np.linalg.norm(seg_vector)
        qh.append(qh_temp)

    xcenters = np.array(cross_section.xcenters)[plot_is]
    ycenters = np.mean(np.array(cross_section.elev).T[plot_ids], axis=1)

    if reverse:
        xcenters = ax.get_xlim()[1] - xcenters
        qh = -np.array(qh)

    ax.quiver(xcenters, ycenters, qh, qz, angles='xy', pivot='tail', color='white', alpha=1)
    print("finished pockmark")


def plot_geology(ls=None, savepath=None, unstructured=False):
    zone_array = discretization.get_zone_array()
    grid = discretization.get_grid_w_ibound()

    if not unstructured:
        f = plt.figure(figsize=(7, 4))
        gs = f.add_gridspec(1, 2, width_ratios=[1, 0.5])
        plt.rcParams.update({'font.size': 8})
        ax = f.add_subplot(gs[0], projection='3d')
        # ax.view_init(elev=30, azim=45, roll=15)
        if ls is None:
            ls = LightSource(285, -60)
        ax.set_box_aspect([4, 5, 2])

        x = np.repeat(grid.xvertices[np.newaxis, :, :], grid.nlay+1, axis=0)
        y = np.repeat(grid.yvertices[np.newaxis, :, :], grid.nlay+1, axis=0)
        z = grid.zverts_smooth
        filled = zone_array>0
        colors = np.empty_like(zone_array[:, :, :], dtype=object)
        colors[zone_array == 0] = "white"
        colors[zone_array == -1] = "white"
        colors[zone_array == 1] = "mediumseagreen"
        colors[zone_array == 2] = "slategray"
        colors[zone_array == 3] = "powderblue"
        colors[zone_array == 4] = "goldenrod"
        ax.voxels(x, y, z, filled, facecolors=colors, edgecolors=colors, alpha=1, lightsource=ls)

        ax_leg = f.add_subplot(gs[1])
        ax_leg.set_axis_off()

        pu = mpl.patches.Patch(color="mediumseagreen", label='Petone Marine upper, sand')
        pm = mpl.patches.Patch(color="slategray", label='Petone Marine distal, silt')
        pl = mpl.patches.Patch(color="powderblue", label='Petone Marine lower, sand')
        uw = mpl.patches.Patch(color="goldenrod", label='Upper Waiwhetu, gravel')

        ax_leg.legend(handles=[pu, pm, pl, uw], loc='center left', fontsize=7)

        ax.set_xlabel("nztm x [m]", labelpad=5.0)
        ax.set_ylabel("nztm y [m]", labelpad=10.0)
        ax.set_zlabel("Depth [m]")
        if savepath is None:
            savepath = "/home/connor/PycharmProjects/pockmarks/figures/geology.png"

        f.savefig(savepath, dpi=600)
    else:
        grid = discretization.get_grid(unstructured=True)
        zones = discretization.get_zone_array(unstructured=True)
        vtk = Vtk(modelgrid=grid)
        vtk.add_array(zones, "zones")
        grid = vtk.to_pyvista()

        p = pv.Plotter()
        p.add_mesh(grid)
        p.show()

        pass


def plot_geo_multi_ls():
    azimuth = np.linspace(240, 360, 8)
    elv = np.linspace(0, -90, 4)
    fig_dir = "/home/connor/PycharmProjects/pockmarks/figures/geo_ls_test"
    if not os.path.exists(fig_dir):
        os.makedirs(fig_dir)

    for az in azimuth:
        for el in elv:
            ls = LightSource(az, el)
            plot_geology(ls, savepath=f"{fig_dir}/geology_az{az}_el{el}.png")

def plot_boundary_conditions(test_boundaries=False):

    if leapfrog_name is None:
        grid = discretization.get_grid_w_ibound()
        f = plt.figure(figsize=(7, 4))
        gs = f.add_gridspec(1, 2, width_ratios=[1, 0.5])
        plt.rcParams.update({'font.size': 8})
        ax = f.add_subplot(gs[0], projection='3d')
        # ax.view_init(elev=30, azim=45, roll=15)
        ls = LightSource(315, 45)
        ax.set_box_aspect([4, 5, 2])

        x = np.repeat(grid.xvertices[np.newaxis, :, :], grid.nlay + 1, axis=0)
        y = np.repeat(grid.yvertices[np.newaxis, :, :], grid.nlay + 1, axis=0)
        z = grid.zverts_smooth
        filled_base = grid.idomain.astype(bool)
        zone_array = discretization.get_zone_array()

        boundary_arrays = discretization.get_boundary_arrays()
        filled_onshore = boundary_arrays['onshore']
        filled_surface = boundary_arrays['seafloor']
        filled_offshore = boundary_arrays['offshore']

        filled_base = np.multiply(filled_base, ~filled_onshore)
        ax.voxels(x, y, z, filled_onshore, facecolors=cm.viridis(0.0), edgecolors=colors.colorConverter.to_rgba(cm.viridis(0), alpha=0), alpha=0.5, lightsource=ls)

        # plot surface boundary cells
        filled_base = np.multiply(filled_base, ~filled_surface)
        ax.voxels(x, y, z, filled_surface, facecolors=cm.viridis(1.0), edgecolors=colors.colorConverter.to_rgba(cm.viridis(1.0), alpha=0), alpha=0.5, lightsource=ls)

        # plot offshore boundary cells
        filled_base = np.multiply(filled_base, ~filled_offshore)
        ax.voxels(x, y, z, filled_offshore, facecolors=cm.viridis(0.5), edgecolors=colors.colorConverter.to_rgba(cm.viridis(0.5), alpha=0), alpha=0.5, lightsource=ls)

        # plot non boundary cells
        ax.voxels(x, y, z, filled_base, facecolors="grey", edgecolors=colors.colorConverter.to_rgba("grey", alpha=0.), alpha=0.05, lightsource=ls)
        #
        ax_leg = f.add_subplot(gs[1])
        ax_leg.set_axis_off()

        on = mpl.patches.Patch(color=cm.viridis(0.0), alpha=0.5, label='Onshore CHD')
        off = mpl.patches.Patch(color=cm.viridis(0.5), alpha=0.5, label='Offshore CHD')
        surf = mpl.patches.Patch(color=cm.viridis(1.0), alpha=0.5, label='Seafloor CHD and CNC')


        ax_leg.legend(handles=[on, off, surf], loc='center left', fontsize=7)

        ax.set_xlabel("nztm x [m]", labelpad=5.0)
        ax.set_ylabel("nztm y [m]", labelpad=10.0)
        ax.set_zlabel("Depth [m]")
        f.savefig("/home/connor/PycharmProjects/pockmarks/figures/boundary_conditions.png", dpi=600)

    else:
        f = plt.figure(figsize=(8, 6))
        plt.rcParams.update({'font.size': 8})
        gs = f.add_gridspec(2, 1, height_ratios=[0.5, 0.4], hspace=0.3)
        gs_top = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[0], width_ratios = [0.5, 1], wspace=0.25)
        ax_plan = f.add_subplot(gs_top[0])
        ax_xsect = f.add_subplot(gs_top[1])
        gs_bot = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[1], width_ratios = [1, 0.8], wspace=0.25)
        ax_p_plan = f.add_subplot(gs_bot[0])
        ax_p_xsect = f.add_subplot(gs_bot[1])

        grid = discretization.get_grid(0.5)

        ax_plan.set_aspect('equal')

        mapview = flopy.plot.PlotMapView(ax=ax_plan, modelgrid=grid, layer=4)
        linecollection = mapview.plot_grid(lw=0.05)
        cross_section = gpd.read_file('/home/connor/PycharmProjects/pockmarks/data/cross_section.shp')
        coords = cross_section['geometry'][0].coords.xy
        ax_plan.plot(coords[0], coords[1], ls="--", color='teal')

        temp_line = [(i, j) for i, j in zip(coords[0], coords[1])]
        line = {"Line": temp_line}
        xsect = flopy.plot.PlotCrossSection(ax=ax_xsect, modelgrid=grid, line=line)
        xsect.plot_grid(lw=0.05)
        lims = ax_xsect.get_xlim()
        ax_xsect.set_xlim(lims[0]-50, lims[1]+50)
        boundaries = discretization.get_boundary_arrays()
        onshore = np.ones_like(grid.botm) * np.nan
        onshore[boundaries['onshore']] = 1
        offshore = np.ones_like(grid.botm) * np.nan
        offshore[boundaries['offshore']] = 1
        seafloor = np.ones_like(grid.botm) * np.nan
        seafloor[boundaries['seafloor']] = 1

        cmap_onshore = LinearSegmentedColormap.from_list("onshore", [cm.viridis(0.0), cm.viridis(0.0)])
        cmap_offshore = LinearSegmentedColormap.from_list("offshore", [cm.viridis(0.5), cm.viridis(0.5)])
        cmap_seafloor = LinearSegmentedColormap.from_list("seafloor", [cm.viridis(1.0), cm.viridis(1.0)])

        xsect.plot_array(onshore, cmap=cmap_onshore, alpha=1, zorder=5, label="Onshore CHD")
        xsect.plot_array(offshore, cmap=cmap_offshore, alpha=1, zorder=5, label="Offshore CHD")
        xsect.plot_array(seafloor, cmap=cmap_seafloor, alpha=1, zorder=5, label="Seafloor CHD and CNC")

        if test_boundaries:
            mapview.plot_array(onshore, cmap=cmap_onshore, alpha=1, zorder=5)
            mapview.plot_array(offshore, cmap=cmap_offshore, alpha=1, zorder=5)

        zoom_ll = (1758102, 5432107)
        zoom_dx = 250
        zoom_dy = 150

        patch = mpatches.Rectangle(xy=zoom_ll, width=zoom_dx, height=zoom_dy, edgecolor='olivedrab', facecolor='none', lw=1,
                                   zorder=12)
        ax_plan.add_patch(patch)

        ax_p_plan.set_aspect('equal')
        mapview_p = flopy.plot.PlotMapView(ax=ax_p_plan, modelgrid=grid, layer=0,
                                           extent=(zoom_ll[0], zoom_ll[0]+zoom_dx, zoom_ll[1], zoom_ll[1]+zoom_dy))
        linecollection_p = mapview_p.plot_grid(lw=0.1)
        ax_p_plan.plot(coords[0], coords[1], ls="--", color='teal')

        zoom_area = Polygon([(zoom_ll[0], zoom_ll[1]), (zoom_ll[0]+zoom_dx, zoom_ll[1]),
                                (zoom_ll[0]+zoom_dx, zoom_ll[1]+zoom_dy), (zoom_ll[0], zoom_ll[1]+zoom_dy)])
        zoom_line = shapely.intersection(zoom_area, cross_section['geometry'][0])

        coords = zoom_line.coords.xy
        temp_line1 = [(i, j) for i, j in zip(coords[0], coords[1])]
        line2 = {"Line": temp_line1}

        xsect_p = flopy.plot.PlotCrossSection(ax=ax_p_xsect, modelgrid=grid, line=line2)
        xsect_p.plot_grid(lw=0.1)

        anchor_b= (1758450, 5433570)
        anchor_c = (zoom_ll[0]+zoom_dx+50, zoom_ll[1]+zoom_dy+50)
        anchor_d = (1758195, 5432240)
        angle_b = np.arctan((temp_line[0][0]-temp_line[1][0])/(temp_line[0][1]-temp_line[1][1]))
        angle_d = np.arctan((temp_line1[0][0]-temp_line1[1][0])/(temp_line1[0][1]-temp_line1[1][1]))
        angle_b = -90-np.rad2deg(angle_b)
        angle_d = -90-np.rad2deg(angle_d)

        tb = ax_plan.text(anchor_b[0], anchor_b[1], "(b)", fontsize=8, ha='left', va='bottom',
                     rotation=angle_b, rotation_mode='anchor')
        tc = ax_plan.text(anchor_c[0], anchor_c[1], "(c)", fontsize=8, ha='left', va='bottom')
        td = ax_p_plan.text(anchor_d[0], anchor_d[1], "(d)", fontsize=8, ha='left', va='bottom',
                       rotation=angle_d, rotation_mode='anchor')
        for t in [tb, tc, td]:
            t.set_bbox(dict(facecolor='white', alpha=0.7, edgecolor='None'))

        ax_plan.set_ylabel("Nztm y [m]")
        ax_plan.set_xlabel("Nztm x [m]")
        ax_xsect.set_ylabel("Nzvd z [m]")
        ax_xsect.set_xlabel("Distance along cross section [m]")
        ax_p_plan.set_ylabel("Nztm y [m]")
        ax_p_plan.set_xlabel("Nztm x [m]")
        ax_p_xsect.set_ylabel("Nzvd z [m]")
        ax_p_xsect.set_xlabel("Distance along cross section [m]")

        ax_p_xsect.set_ylim(-35, -10)

        onshore_patch = mpatches.Patch(color=cm.viridis(0.0), alpha=0.5, label='Onshore CHD')
        offshore_patch = mpatches.Patch(color=cm.viridis(0.5), alpha=0.5, label='Offshore CHD')
        seafloor_patch = mpatches.Patch(color=cm.viridis(1.0), alpha=0.5, label='Seafloor CHD and CNC')

        ax_xsect.legend(handles=[onshore_patch,offshore_patch,seafloor_patch],loc='upper right', fontsize=6)
        ax_plan.text(0.9, 1.05, '(a)', fontsize=9, transform=ax_plan.transAxes)
        ax_xsect.text(0.95, 1.05, '(b)', fontsize=9, transform=ax_xsect.transAxes)
        ax_p_plan.text(0.95, 1.05, '(c)', fontsize=9, transform=ax_p_plan.transAxes)
        ax_p_xsect.text(0.9, 1.05, '(d)', fontsize=9, transform=ax_p_xsect.transAxes)

        # ax_p_xsect.set_box_aspect(ax_p_plan.get_box_aspect())
        # f.set_constrained_layout(True)
        plt.savefig(f"/home/connor/PycharmProjects/pockmarks/figures/{leapfrog_name}_boundary_conditions.png", dpi=600)

def plot_aquitard_thickness(modelgrid=False):

    if modelgrid:
        offshore = discretization.get_offshore()
        grid = discretization.get_grid_w_ibound()
        ibound = grid.idomain
        zone_data = discretization.get_zone_array()
        f, ax = plt.subplots()
        pmv = flopy.plot.PlotMapView(ax=ax, modelgrid=grid, layer=18)
        distal_thickness = (np.argmax(zone_data>2, axis=0) - np.argmax(zone_data==2, axis=0))*grid.delz[0][0][0]
        thickness = pmv.plot_array(np.multiply(distal_thickness, offshore), cmap="bone_r", vmin=0, vmax=20)
        plt.colorbar(thickness, label="Silt thickness [m]")
        ax.set_aspect('equal')
        ax.set_xlabel("nztm x [m]", labelpad=5.0)
        ax.set_ylabel("nztm y [m]", labelpad=10.0)
        f.savefig("/home/connor/PycharmProjects/pockmarks/figures/aquitard_thickness.png", dpi=600)
    else:
        mc21 = (-41.2416, 174.8891)
        mc22 = (-41.242, 174.8873)
        mc23 = (-41.2394, 174.8873)
        # mc14 = (-41.2418, 174.8896)
        sv1 = (-41.242, 174.886)
        sv6 = (-41.24787+0.0018, 174.8878)
        # points = np.genfromtxt("/home/connor/PycharmProjects/pockmarks/data/gis/coordinates.csv", delimiter=',', skip_header=1)[:, 1:]
        points = np.genfromtxt("/home/connor/PycharmProjects/pockmarks/data/gis/coordinates_from_harding_appendix.csv",
                               delimiter=',', skip_header=1)[:,1:]
        mc_names = ["MC2-1", "MC2-2", "MC2-3"]
        sv_names = ["SV1", "SV6"]
        mc_coords = [points[0], points[1], points[2]]
        sv_coords = [points[3], points[5]]

        transformer = Transformer.from_crs("EPSG:4326", "EPSG:2193")

        f, ax = plt.subplots(figsize=(5, 4))
        plt.rcParams.update({'font.size': 7})
        ax.set_aspect('equal')
        thickness = rasterio.open("/home/connor/PycharmProjects/pockmarks/data/gis/marine_thickness_masked_model.tif")
        thck = show(thickness, ax=ax, cmap='grey', alpha=0.7, zorder=5)
        first = True
        for name, coord in zip(mc_names, mc_coords):
            nztm = transformer.transform(*coord)
            ax.scatter(nztm[1], nztm[0], marker='v', color='peru', lw=0.75, s=15, label="Sediment Cores" if first else None, zorder=15, facecolors='none')
            first = False
        first = True
        for name, coord in zip(sv_names, sv_coords):
            coord[0] += +0.0016666
            nztm = transformer.transform(*coord)
            print(nztm)
            ax.scatter(nztm[1], nztm[0], marker='^', color='blue', lw=0.75,  s=15, label="Discharge Observations" if first else None, zorder=15, facecolors='none')
            first = False

        points[4][0] += +0.0016666
        print(transformer.transform(*points[4]))

        all_names = mc_names + sv_names
        all_coords = mc_coords + sv_coords
        ha = ['left', 'right', 'right', 'left', 'right', 'right']
        va = ['bottom', 'bottom', 'bottom', 'top', 'bottom', 'bottom']
        for name, coord, ha, va in zip(all_names, all_coords, ha, va):
            nztm = transformer.transform(*coord)
            if ha == 'left':
                x = nztm[1] + 30
            else:
                x = nztm[1] - 30
            if va == 'bottom':
                y = nztm[0] + 20
            else:
                y = nztm[0] - 20
            ax.text(x, y, name, fontsize=7, ha=ha, va=va, zorder=20)
        ax.set_ylabel("Nztm y [m]")
        ax.set_xlabel("Nztm x [m]")
        ax.legend(loc='upper right', fontsize=7)
        norm = mpl.colors.Normalize(vmin=thickness.statistics(1).min, vmax=thickness.statistics(1).max)
        mappable = plt.cm.ScalarMappable(cmap='grey', norm=norm)
        plt.colorbar(mappable, ax=ax, label="Petone Marine Distal (silt) thickness [m]", extend='both', alpha=0.7)
        plt.show()
        f.savefig("/home/connor/PycharmProjects/pockmarks/figures/aquitard_thickness_map.png", dpi=600)

def plot_aquitards_resistance():
    offshore = discretization.get_offshore()
    grid = discretization.get_grid_w_ibound()
    ibound = grid.idomain
    zone_data = discretization.get_zone_array()
    f, ax = plt.subplots()
    pmv = flopy.plot.PlotMapView(ax=ax, modelgrid=grid, layer=18)
    h_params = parameterization.get_hydraulic_params()
    resistance = np.zeros((grid.nrow, grid.ncol))

    for row in range(grid.nrow):
        for col in range(grid.ncol):
            resistance[row, col] += len(np.where(zone_data[:, row, col] ==1))*1/h_params['vkpu']
            resistance[row, col] += len(np.where(zone_data[:, row, col]==2))*1/h_params['vkpd']
            resistance[row, col] += len(np.where(zone_data[:, row, col]==3))*1/h_params['vkpl']
            resistance[row, col] = len(np.where(zone_data[:, row, col]>=1) and np.where(zone_data[:, row, col]<=3))*grid.delz[0][0][0]*resistance[row, col]

    r_plot = pmv.plot_array(resistance, cmap="bone_r")
    plt.colorbar(r_plot, label="Aquitard resistance [days]")
    ax.set_aspect('equal')
    ax.set_xlabel("nztm x [m]", labelpad=5.0)
    ax.set_ylabel("nztm y [m]", labelpad=10.0)
    f.savefig("/home/connor/PycharmProjects/pockmarks/figures/aquitard_resistance.png", dpi=600)


def plot_cross_section_boundaries(name):
    sim = flopy.mf6.MFSimulation.load(
        sim_ws=f"/home/connor/PycharmProjects/pockmarks/models/{name}/",
    )
    gwf = sim.get_model(name)
    grid = discretization.get_grid_w_ibound()
    xsect = flopy.plot.PlotCrossSection(model=gwf, line={"Row": grid.nrow//2})
    xsect.plot_ibound()
    xsect.plot_bc('CHD', color='yellow')
    xsect.plot_bc('GHB', color='purple')
    plt.show()

def plot_surface_fluxes(name, timestep=-1, wrong_formulation=False):
    grid = discretization.get_grid_w_ibound()

    sim = flopy.mf6.MFSimulation.load(
        sim_ws=f"/home/connor/PycharmProjects/pockmarks/models/{name}/",
    )
    gwt = sim.get_model("trans")
    conc = gwt.output.concentration().get_data(totim=gwt.output.concentration().get_times()[timestep])
    seafloor = discretization.get_cells_below_seafloor(wrong_formulation)
    conc = np.nansum(np.multiply(conc, seafloor), axis=0)
    gwf = sim.get_model(name)
    bud = gwf.output.budget()
    qx, qy, qz = flopy.utils.postprocessing.get_specific_discharge(
                    bud.get_data(text="DATA-SPDIS", totim=bud.get_times()[timestep])[0],
                    gwf,
                )
    qz = np.nansum(np.multiply(qz, seafloor), axis=0)
    f, axs = plt.subplots(1, 2, figsize=(8, 5))
    c = flopy.plot.PlotMapView(ax=axs[0], modelgrid=grid)
    q = flopy.plot.PlotMapView(ax=axs[1], modelgrid=grid)
    c_cm = c.plot_array(conc, cmap='viridis', vmin=0, vmax=35)
    q_cm = q.plot_array(qz, cmap='PuOr', norm=colors.SymLogNorm(linthresh=1e-6, linscale=1e-6, vmin=-1e-4, vmax=1e-4))
    plt.colorbar(ax=axs[0], mappable=c_cm)
    plt.colorbar(ax=axs[1], mappable=q_cm)
    plt.tight_layout()
    plt.savefig(f"/home/connor/PycharmProjects/pockmarks/figures/surface_fluxes_{name}.png", dpi=600)
    pass

def plot_cross_section_flux_salinity(name, simple=True):

    f, axs = plt.subplots(2, 1, figsize=(4, 6))

    data = utils.get_final_results(name)
    conc = data['conc']
    qx = data['qx']
    qy = data['qy']
    qz = data['qz']

    if leapfrog_name is None:
        grid = discretization.get_grid_w_ibound()
        line = {"Row": grid.nrow // 2}
    else:
        grid = discretization.get_grid_w_ibound()
        if not simple:
            cross_section = gpd.read_file('/home/connor/PycharmProjects/pockmarks/data/gis/cross_section.shp')
        else:
            cross_section = gpd.read_file('/home/connor/PycharmProjects/pockmarks/data/gis/cross_section_simple.shp')
        coords = cross_section['geometry'][0].coords.xy
        temp_line = [(i, j) for i, j in zip(coords[0], coords[1])]
        line = {"Line": temp_line}

    xsect = flopy.plot.PlotCrossSection(ax=axs[0], modelgrid=grid, line=line)
    xsect.plot_array(conc, cmap='viridis', vmax=35, zorder=2)
    # xsect.plot_ibound(zorder=1)
    xsect2 = flopy.plot.PlotCrossSection(ax=axs[1], modelgrid=grid, line=line)
    q_cm = xsect2.plot_array(qz, cmap='PuOr',
                                norm=colors.SymLogNorm(linthresh=1e-6, linscale=1e-6, vmin=-1e-4, vmax=1e-4))
    plt.savefig(f"/home/connor/PycharmProjects/pockmarks/figures/xsect_flux_salinity_{name}.png", dpi=600)





def plot_cross_sections_and_surface_fluxes_over_time(name):
    grid = discretization.get_grid_w_ibound()

    sim = flopy.mf6.MFSimulation.load(
        sim_ws=f"/home/connor/PycharmProjects/pockmarks/models/{name}/",
    )
    gwt = sim.get_model("trans")
    times = gwt.output.concentration().get_times()
    seafloor = discretization.get_cells_below_seafloor()
    f, axs = plt.subplots(2, 5, figsize=(10, 6))
    for i, ax in enumerate(axs.flat):
        conc = gwt.output.concentration().get_data(totim=times[i])
        conc = np.nansum(np.multiply(conc, seafloor), axis=0)
        c = flopy.plot.PlotMapView(ax=ax, modelgrid=grid)
        c.plot_array(conc, cmap='viridis', vmin=0, vmax=35)
        ax.set_title(f"t= {int(times[i])} days")
        ax.set_xticklabels([])
        ax.set_yticklabels([])

    f.savefig(f"/home/connor/PycharmProjects/pockmarks/figures/seafloor_conc_over_time_{name}.png", dpi=600)

    f, axs = plt.subplots(4, 2, figsize=(8, 10))
    for i, ax in enumerate(axs.flat):
        xsect = flopy.plot.PlotCrossSection(ax=ax, model=gwt, line={"Row": grid.nrow // 2})
        conc = gwt.output.concentration().get_data(totim=times[i])
        xsect.plot_array(conc, cmap='viridis', vmax=35, zorder=2)
        xsect.plot_ibound(zorder=1)
        ax.set_title(f"t= {int(times[i])} days")
        ax.set_xticklabels([])
        ax.set_yticklabels([])

    f.savefig(f"/home/connor/PycharmProjects/pockmarks/figures/xsect_conc_over_time_{name}.png", dpi=600)

def plot_cross_section_vectors(name, simple=True):

    f, ax = plt.subplots(figsize=(10, 6))
    if not os.path.exists(f"/home/connor/PycharmProjects/pockmarks/outputs/{name}/final_step.npz"):
        sim = flopy.mf6.MFSimulation.load(
            sim_ws=f"/home/connor/PycharmProjects/pockmarks/models/{name}/",
        )
        gwt = sim.get_model("trans")

        gwf = sim.get_model(name)
        bud = gwf.output.budget()
        qx, qy, qz = flopy.utils.postprocessing.get_specific_discharge(
            bud.get_data(text="DATA-SPDIS", totim=bud.get_times()[-1])[0],
            gwf,
        )
        times = gwt.output.concentration().get_times()
        conc = gwt.output.concentration().get_data(totim=times[-1])
        os.makedirs(f"/home/connor/PycharmProjects/pockmarks/outputs/{name}/", exist_ok=True)
        np.savez_compressed(f"/home/connor/PycharmProjects/pockmarks/outputs/{name}/final_step.npz",
                            conc=conc, qx=qx, qy=qy, qz=qz)
    else:
        data = np.load(f"/home/connor/PycharmProjects/pockmarks/outputs/{name}/final_step.npz")
        conc = data['conc']
        qx = data['qx']
        qy = data['qy']
        qz = data['qz']

    if leapfrog_name is None:
        grid = discretization.get_grid_w_ibound()
        line = {"Row": grid.nrow // 2}
    else:
        grid = discretization.get_grid_w_ibound()
        if not simple:
            cross_section = gpd.read_file('/home/connor/PycharmProjects/pockmarks/data/gis/cross_section.shp')
        else:
            cross_section = gpd.read_file('/home/connor/PycharmProjects/pockmarks/data/gis/cross_section_simple.shp')
        coords = cross_section['geometry'][0].coords.xy
        temp_line = [(i, j) for i, j in zip(coords[0], coords[1])]
        line = {"line": temp_line}

    xsect = flopy.plot.PlotCrossSection(ax=ax, modelgrid=grid, line=line)
    xsect.plot_array(conc, cmap='viridis', vmax=35, zorder=2)
    # xsect.plot_ibound(zorder=1)

    quiver = xsect.plot_vector(
        qx,
        qy,
        qz,
        hstep=10,
        kstep=4,
        normalize=True,
        color="white",
        zorder=15
    )

    plt.savefig(f"/home/connor/PycharmProjects/pockmarks/figures/xsect_vectors_{name}_simple={simple}.png", dpi=600)

    pass

def plot_surface_vectors(name):
    grid = discretization.get_grid_w_ibound()
    sim = flopy.mf6.MFSimulation.load(
        sim_ws=f"/home/connor/PycharmProjects/pockmarks/models/{name}/",
    )
    gwt = sim.get_model("trans")

    gwf = sim.get_model(name)
    bud = gwf.output.budget()
    qx, qy, qz = flopy.utils.postprocessing.get_specific_discharge(
        bud.get_data(text="DATA-SPDIS", totim=bud.get_times()[-1])[0],
        gwf,
    )
    seafloor = discretization.get_cells_below_seafloor()
    qx = np.nansum(np.multiply(qx, seafloor), axis=0)
    qy = np.nansum(np.multiply(qy, seafloor), axis=0)

    conc = gwt.output.concentration().get_data(totim=gwt.output.concentration().get_times()[-1])
    seafloor = discretization.get_cells_below_seafloor()
    conc = np.nansum(np.multiply(conc, seafloor), axis=0)

    f, ax = plt.subplots(figsize=(6, 8))
    c = flopy.plot.PlotMapView(ax=ax, modelgrid=grid)
    c.plot_array(conc, cmap='viridis', vmin=0, vmax=35)
    c.plot_vector(qx, qy, istep=20, jstep=20, normalize=True, color="white", zorder=15)

    plt.savefig(f"/home/connor/PycharmProjects/pockmarks/figures/surface_vectors_{name}.png", dpi=600)

def plot_pockmark(name, i, grid, ax_plan=None, ax_section=None, label_sections=False,
                  ylim=None, vectors=True, contours=False, section_line=True, discharge_on_plan=False):

    data = utils.get_final_results(name)
    pockmark = gpd.read_file('/home/connor/PycharmProjects/pockmarks/data/pockmarks6_new.shp', crs='EPSG:2193')['geometry'][i]
    section = gpd.read_file('/home/connor/PycharmProjects/pockmarks/data/pockmark_sections.shp')['geometry'][i]
    colors = ['tab:red', 'tab:orange', 'tab:brown', 'tab:purple', 'tab:blue', 'tab:pink']
    color = colors[i]

    bounds = shapely.union(pockmark, section).bounds
    if bounds[2] - bounds[0] > bounds[3] - bounds[1]:
        add = ((bounds[2] - bounds[0]) - (bounds[3] - bounds[1])) / 2
        square_extent = (bounds[0], bounds[2], bounds[1] - add, bounds[3] + add)
    else:
        add = ((bounds[3] - bounds[1]) - (bounds[2] - bounds[0])) / 2
        square_extent = (bounds[0] - add, bounds[2] + add, bounds[1], bounds[3])

    if ax_plan is not None:
        pockmark_plan = flopy.plot.PlotMapView(ax=ax_plan, modelgrid=grid, layer=1, extent=square_extent)
    # pockmark_plan.plot_grid()
        pockmark_plan.plot_array(data['conc'], vmin=0, vmax=35)

    if ax_section is not None:
        coords = section.coords.xy
        temp_line = [(i, j) for i, j in zip(coords[0], coords[1])]
        pockmark_section = flopy.plot.PlotCrossSection(ax=ax_section, modelgrid=grid, line={"Line": temp_line})
        ax_section.set_ylim(-30, -5)
        # hk = parameterization.get_params_by_zone('base', 'hk', 'high_pockmark_k')
        # pockmark_section.plot_array(np.log10(hk), cmap='coolwarm', zorder=10)
        # lims = ax_section.get_xlim()
        # ax_section.set_xlim(lims[0] - 50, lims[1] + 50)

    if ax_plan is not None and section_line:
        ax_plan.plot(coords[0], coords[1], ls="--", color=color)

        angle_start = np.arctan((coords[1][1] - coords[1][0]) / (coords[0][1] - coords[0][0]))
        angle_start = np.abs(np.degrees(angle_start))
        if coords[0][1] > coords[0][0]:
            if coords[1][1] < coords[1][0]:
                mult = 1
                angle_start = -angle_start
            else:
                mult = - 1
        else:
            if coords[1][1] > coords[1][0]:
                mult = - 1
                angle_start = -angle_start - 180
            else:
                mult = 1
                angle_start = 180 + angle_start

        anchor_start = (coords[0][0] + mult * 10, coords[1][0] + mult * 10)
        start_vec = [coords[0][1] - coords[0][0], coords[1][1] - coords[1][0]] / np.linalg.norm(
            [coords[0][1] - coords[0][0], coords[1][1] - coords[1][0]])
        start_point_along_line = (coords[0][0] + 20 * start_vec[0], coords[1][0] + 20 * start_vec[1])
        orthogonal = (1, -start_vec[0] / start_vec[1])
        orthogonal = orthogonal / np.linalg.norm(orthogonal)

        anchor_start = (
        start_point_along_line[0] + mult * 10 * orthogonal[0], start_point_along_line[1] + mult * 10 * orthogonal[1])

        tstart = ax_plan.text(anchor_start[0], anchor_start[1], f"{i + 1}", fontsize=6, ha='left', va='bottom',
                              rotation=angle_start, rotation_mode='anchor')

        angle_end = np.arctan((coords[1][-2] - coords[1][-1]) / (coords[0][-2] - coords[0][-1]))
        angle_end = np.abs(np.degrees(angle_end))
        if coords[0][-1] > coords[0][-2]:
            if coords[1][-1] < coords[1][-2]:
                mult = 1
                angle_end = -angle_end
            else:
                mult = - 1
        else:
            if coords[1][-1] > coords[1][-2]:
                mult = - 1
                angle_end = -angle_end - 180
            else:
                mult = 1
                angle_end = 180 + angle_end

        anchor_end = (coords[0][-1] - mult * 10, coords[1][-1] + mult * 10)
        end_vec = [coords[0][-2] - coords[0][-1], coords[1][-2] - coords[1][-1]] / np.linalg.norm(
            [coords[0][-2] - coords[0][-1], coords[1][-2] - coords[1][-1]])
        end_point_along_line = (coords[0][-1] + 20 * end_vec[0], coords[1][-1] + 20 * end_vec[1])
        orthogonal = (1, -end_vec[0] / end_vec[1])
        orthogonal = orthogonal / np.linalg.norm(orthogonal)
        anchor_end = (
        end_point_along_line[0] + mult * 10 * orthogonal[0], end_point_along_line[1] + mult * 10 * orthogonal[1])

        tstart = ax_plan.text(anchor_end[0], anchor_end[1], f"{i + 1}'", fontsize=6, ha='right', va='bottom',
                              rotation=angle_end, rotation_mode='anchor')

    zone_array = discretization.get_zone_array()

    if ax_section is not None:
        pockmark_section.plot_array(data['conc'], vmin=0, vmax=35)
        if vectors:
            plot_unstructured_vectors(ax_section, pockmark_section, data['qx'], data['qy'], data['qz'], section, 4, 4,
                                      zone_array)
        if contours:
            contours = pockmark_section.contour_array(data['conc'], levels=[35*0.05], linewidths=0.5, colors='red', linestyles='dashed')
            # ax_section.clabel(contours, inline=True, fontsize=5, fmt={35*0.05: r"5% salinity"})

    if ax_plan is not None:
        ax_plan.set_ylabel('Nztm y [m]')
        ax_plan.set_box_aspect(1)
        ax_plan.set_xlabel('Nztm x [m]', labelpad=15)

    if ax_section is not None:
        pockmarks_data = utils.get_final_pockmark_results(name)
        ax_section.set_ylabel('Nzvd z [m]')
        ax_section.set_box_aspect(1)
        ax_section.set_xlabel('Distance along section [m]')
        if ylim is None:
            ax_section.set_ylim(-30, -5)
        else:
            ax_section.set_ylim(ylim[0], ylim[1])

        if label_sections:
            ax_section.text(0.025, 0.975, f"{i + 1}", fontsize=6, transform=ax_section.transAxes, ha='left', va='top')
            ax_section.text(0.975, 0.975, f"{i + 1}'", fontsize=6, transform=ax_section.transAxes, ha='right', va='top')

        ax_section.text(0.025, 0.01, rf"{pockmarks_data.iloc[i]['discharge']:.2f} m$^3$/day",
                        fontsize=6, transform=ax_section.transAxes, color="white", ha='left', va='bottom')
        ax_section.text(0.975, 0.01, f"{pockmarks_data.iloc[i]['concentration']:.1f} PSU" if pockmarks_data.iloc[i]['discharge'] > 0 else "35 PSU",
                        fontsize=6, transform=ax_section.transAxes, color='white', ha='right', va='bottom')

    if discharge_on_plan:
        pockmarks_data = utils.get_final_pockmark_results(name)
        ax_plan.text(0.025, 0.01, rf"{pockmarks_data.iloc[i]['discharge']:.2f} m$^3$/day",
                        fontsize=6, transform=ax_plan.transAxes, color="black", ha='left', va='bottom')
        ax_plan.text(0.975, 0.01, f"{pockmarks_data.iloc[i]['concentration']:.1f} PSU" if pockmarks_data.iloc[i]['discharge'] > 0 else "35 PSU",
                        fontsize=6, transform=ax_plan.transAxes, color='black', ha='right', va='bottom')


def plot_salinity_at_pockmarks(name, ylim=None):
    # load pockmarks
    pockmarks = gpd.read_file('/home/connor/PycharmProjects/pockmarks/data/pockmarks6_new.shp', crs='EPSG:2193')
    sections = gpd.read_file('/home/connor/PycharmProjects/pockmarks/data/pockmark_sections.shp')
    # load grid
    grid = discretization.get_grid_w_ibound()
    # load concentration
    data = utils.get_final_results(name)
    conc = data['conc']
    # setup gridspec
    f = plt.figure(figsize=(7.5, 8.5))

    gs = f.add_gridspec(6, 4, hspace=0.5, wspace=0.5, height_ratios=[1, 1, 1, 1, 0, 0.1])
    plt.rcParams.update({'font.size': 7})
    ax_map = f.add_subplot(gs[:2, :2])

    axs_plan = []
    axs_section = []

    axs_plan.append(f.add_subplot(gs[0, 2]))
    axs_section.append(f.add_subplot(gs[0, 3]))
    axs_plan.append(f.add_subplot(gs[1, 2]))
    axs_section.append(f.add_subplot(gs[1, 3]))
    axs_plan.append(f.add_subplot(gs[2, 0]))
    axs_section.append(f.add_subplot(gs[2, 1]))
    axs_plan.append(f.add_subplot(gs[2, 2]))
    axs_section.append(f.add_subplot(gs[2, 3]))
    axs_plan.append(f.add_subplot(gs[3, 0]))
    axs_section.append(f.add_subplot(gs[3, 1]))
    axs_plan.append(f.add_subplot(gs[3, 2]))
    axs_section.append(f.add_subplot(gs[3, 3]))

    ax_cb = f.add_subplot(gs[-1, 1:3])

    colors = ['tab:red', 'tab:orange', 'tab:brown', 'tab:purple', 'tab:blue', 'tab:pink']

    # plot_map
    # mapview = flopy.plot.PlotMapView(ax=ax_map, modelgrid=grid, layer=4)
    base = rasterio.open("/home/connor/PycharmProjects/pockmarks/map/data/basemap.tif")
    model_area = gpd.read_file('/home/connor/PycharmProjects/pockmarks/data/model_area0.shp')
    extent = model_area.total_bounds
    show(base, ax=ax_map, alpha=0.7, zorder=0)
    ax_map.set_xlim(extent[0], extent[2])
    ax_map.set_ylim(extent[1], extent[3])
    ax_map.set_xlabel('Nztm x [m]')
    ax_map.set_ylabel('Nztm y [m]')
    # linecollection = mapview.plot_grid(lw=0.05)
    pockmarks.to_crs('EPSG:2193')
    pockmarks.plot(ax=ax_map, edgecolors='black', color=colors, zorder=12, alpha=1, lw=0.75)

    pockmarks_data = utils.get_final_pockmark_results(name)

    zone_array = discretization.get_zone_array()

    label_sides = ['left', 'left', 'left', 'right', 'right', 'right']

    for i, (color, pockmark, section, ax_plan, ax_section, label_side) in enumerate(zip(
        colors, pockmarks['geometry'], sections['geometry'], axs_plan, axs_section, label_sides
    )):
        bounds = shapely.union(pockmark, section).bounds
        if bounds[2]-bounds[0] > bounds[3]-bounds[1]:
            add = ((bounds[2] - bounds[0]) - (bounds[3] - bounds[1]))/2
            square_extent = (bounds[0], bounds[2], bounds[1]-add, bounds[3]+add)
        else:
            add = ((bounds[3] - bounds[1]) - (bounds[2] - bounds[0]))/2
            square_extent = (bounds[0]-add, bounds[2]+add, bounds[1], bounds[3])

        if label_side == 'left':
            ax_map.text(pockmark.bounds[0]-40, pockmark.centroid.y, f"{i+1}", fontsize=6, ha='right', va='center')
        else:
            ax_map.text(pockmark.bounds[2]+40, pockmark.centroid.y, f"{i+1}", fontsize=6, ha='left', va='center')


        pockmark_plan = flopy.plot.PlotMapView(ax=ax_plan, modelgrid=grid, layer=1, extent=square_extent)
        # pockmark_plan.plot_grid()
        pockmark_plan.plot_array(conc, vmin=0, vmax=35)
        coords = section.coords.xy

        temp_line = [(i, j) for i, j in zip(coords[0], coords[1])]

        pockmark_section = flopy.plot.PlotCrossSection(ax=ax_section, modelgrid=grid, line={"Line": temp_line})


        ax_plan.plot(coords[0], coords[1], ls="--", color=color)

        angle_start = np.arctan((coords[1][1] - coords[1][0]) / (coords[0][1] - coords[0][0]))
        angle_start = np.abs(np.degrees(angle_start))
        if coords[0][1] > coords[0][0]:
            if coords[1][1] < coords[1][0]:
                mult = 1
                angle_start = -angle_start
            else:
                mult = - 1
        else:
            if coords[1][1] > coords[1][0]:
                mult = - 1
                angle_start = -angle_start-180
            else:
                mult = 1
                angle_start = 180+angle_start


        anchor_start = (coords[0][0] + mult*10, coords[1][0] + mult*10)
        start_vec = [coords[0][1] - coords[0][0], coords[1][1] - coords[1][0]]/np.linalg.norm([coords[0][1] - coords[0][0], coords[1][1] - coords[1][0]])
        start_point_along_line = (coords[0][0] + 20*start_vec[0], coords[1][0] + 20*start_vec[1])
        orthogonal = (1, -start_vec[0]/start_vec[1])
        orthogonal = orthogonal/np.linalg.norm(orthogonal)

        anchor_start = (start_point_along_line[0] + mult * 10*orthogonal[0], start_point_along_line[1] + mult * 10*orthogonal[1])

        tstart = ax_plan.text(anchor_start[0], anchor_start[1], f"{i+1}", fontsize=6, ha='left', va='bottom',
                          rotation=angle_start, rotation_mode='anchor')

        angle_end = np.arctan((coords[1][-2] - coords[1][-1]) / (coords[0][-2] - coords[0][-1]))
        angle_end = np.abs(np.degrees(angle_end))
        if coords[0][-1] > coords[0][-2]:
            if coords[1][-1] < coords[1][-2]:
                mult = 1
                angle_end = -angle_end
            else:
                mult = - 1
        else:
            if coords[1][-1] > coords[1][-2]:
                mult = - 1
                angle_end = -angle_end- 180
            else:
                mult = 1
                angle_end = 180 + angle_end

        anchor_end = (coords[0][-1] - mult*10, coords[1][-1] + mult*10)
        end_vec = [coords[0][-2] - coords[0][-1], coords[1][-2] - coords[1][-1]] / np.linalg.norm(
            [coords[0][-2] - coords[0][-1], coords[1][-2] - coords[1][-1]])
        end_point_along_line = (coords[0][-1] + 20 * end_vec[0], coords[1][-1] + 20 * end_vec[1])
        orthogonal = (1, -end_vec[0] / end_vec[1])
        orthogonal = orthogonal / np.linalg.norm(orthogonal)
        anchor_end = (end_point_along_line[0] + mult * 10 * orthogonal[0], end_point_along_line[1] + mult * 10 * orthogonal[1])

        tstart = ax_plan.text(anchor_end[0], anchor_end[1], f"{i+1}'", fontsize=6, ha='right', va='bottom',
                              rotation=angle_end, rotation_mode='anchor')

        pockmark_section.plot_array(conc, vmin=0, vmax=35)
        if i in [3,4]:
            reverse=True
        else:
            reverse = False
        plot_unstructured_vectors(ax_section, pockmark_section, data['qx'], data['qy'], data['qz'], section,4, 4, zone_array, reverse=reverse)

        ax_plan.set_ylabel('Nztm y [m]')
        ax_section.set_ylabel('Nzvd z [m]')
        ax_section.set_ylim(-30, -5)

        ax_plan.set_box_aspect(1)
        ax_section.set_box_aspect(1)


        if i > 3:
            ax_plan.set_xlabel('Nztm x [m]', labelpad=15)
            ax_section.set_xlabel('Distance along section [m]')

        if reverse:
            # get old ticks and labels
            max = ax_section.get_xlim()[1]
            ax_section.invert_xaxis()
            # set new ticks as tick label - max
            ticks = ax_section.get_xticks()
            new_labels = [f"{max - tick:.0f}" for tick in ticks]
            ax_section.set_xticklabels(new_labels)

        ax_section.text(0.025, 0.975, f"{i+1}", fontsize=6, transform=ax_section.transAxes, ha='left', va='top')
        ax_section.text(0.975, 0.975, f"{i+1}'", fontsize=6, transform=ax_section.transAxes, ha='right', va='top')

        ax_section.text(0.025, 0.01, rf"{pockmarks_data.iloc[i]['discharge']:.2f} m$^3$/day",
                        fontsize=6, transform=ax_section.transAxes, color="white", ha='left', va='bottom')
        ax_section.text(0.975, 0.01, f"{pockmarks_data.iloc[i]['concentration']:.1f} PSU",
                        fontsize=6, transform=ax_section.transAxes, color='white', ha='right', va='bottom')


    norm = mpl.colors.Normalize(vmin=0, vmax=35)
    mappable = plt.cm.ScalarMappable(cmap='viridis', norm=norm)
    cb = plt.colorbar(mappable, cax=ax_cb, label='Salinity [PSU]', orientation='horizontal', fraction=1)

    non_pockmark_data = utils.get_total_seafloor_discharge(name)
    ax_map.text(0.025, 0.055, rf"{non_pockmark_data['discharge']:.2f} m$^3$/day",
                        fontsize=6, transform=ax_map.transAxes, color="white", ha='left', va='bottom')
    ax_map.text(0.025, 0.015, f"{non_pockmark_data['conc']:.1f} PSU",
                        fontsize=6, transform=ax_map.transAxes, color='white', ha='left', va='bottom')

    ax_map.text(0.975, 1.025, "(a)", fontsize=7, transform=ax_map.transAxes, ha='right', va='bottom')
    for ax, let in zip(axs_plan, ['(b)', '(d)', '(f)', '(h)', '(j)', '(l)']):
        ax.text(0.975, 1.025, let, fontsize=7, transform=ax.transAxes, ha='right', va='bottom')
    for ax, let in zip(axs_section, ['(c)', '(e)', '(g)', '(i)', '(k)', '(m)']):
        ax.text(0.975, 1.025, let, fontsize=7, transform=ax.transAxes, ha='right', va='bottom')


    for ax in np.array([ax_map, *axs_plan, *axs_section]):
        ax.fill_between(ax.get_xlim(), ax.get_ylim()[0], ax.get_ylim()[1], color='white', alpha=1, zorder=-1)

    figdir = unbacked_dir.joinpath('figures')
    figdir.mkdir(parents=True, exist_ok=True)
    plt.savefig(figdir.joinpath(f'{name}_pockmark_salinity.png'), dpi=600, transparent=True)

    # plot concentration at pockmarks
    # plan view
    # cross section

def plot_changes_and_examples(set, names, examples = ((0, 5), (1, 5), (4, 5), (5, 5)), color_bars_by_c=False, version=0):

    if version == 0:
        f = plt.figure(figsize=(8, 6))
        gs = f.add_gridspec(4,4, height_ratios=[1, 0.75, 0, 0.1], hspace=0.65, wspace=0.25)
        plt.rcParams.update({'font.size': 8})
        ax_discharges = f.add_subplot(gs[0, :2])
        ax_concentrations = f.add_subplot(gs[0, 2:])
        axs_pockmarks = [f.add_subplot(gs[1, i]) for i in range(4)]
        ax_cb = f.add_subplot(gs[-1, 1:3])

        pockmarks_discharge = []
        pockmarks_concentration = []
        total_discharge = []
        total_concentration = []
        for name in names:
            temp_pockmarks = utils.get_final_pockmark_results(name)
            pockmarks_discharge.append(temp_pockmarks['discharge'].sum())
            pockmarks_concentration.append((temp_pockmarks['concentration']*temp_pockmarks['discharge']).sum()/temp_pockmarks['discharge'].sum())
            temp_total = utils.get_total_seafloor_discharge(name)
            total_discharge.append(temp_total['discharge'])
            total_concentration.append(temp_total['conc'])

        norm = mpl.colors.Normalize(vmin=0, vmax=35)
        cmap = plt.cm.get_cmap('viridis')

        if color_bars_by_c:
            color_pockmark = cmap(norm(pockmarks_concentration))
            color_total = cmap(norm(total_concentration))
        else:
            # colors for pockmarks['tab:red', 'tab:orange', 'tab:brown', 'tab:purple', 'tab:blue', 'tab:pink']
            color_pockmark = ['tab:olive', 'tab:green', 'tab:cyan', 'tab:grey', 'orchid', 'chocolate']
            color_total = ['tab:olive', 'tab:green', 'tab:cyan', 'tab:grey', 'orchid', 'chocolate']

        X_axis = np.arange(len(names)/3)
        ax_discharges.bar(X_axis - 0.2, pockmarks_discharge, 0.4, color=color_pockmark , label='Pockmarks')
        ax_discharges.bar(X_axis + 0.2, total_discharge, 0.4, color=color_total,
                          alpha = 1 if color_bars_by_c else 0.3, label='Distributed')
        ax_discharges.set_xticklabels([None] + [f"{i+1}" for i in range(len(names))])
        ax_discharges.set_xlabel('Scenario')
        ax_discharges.set_ylabel(r'Discharge [m$^3$/day]')

        p_handle = mpatches.Patch(color='tab:grey', label='Pockmark')
        t_handle = mpatches.Patch(color='tab:grey', alpha=0.3, label='Distributed')

        ax_discharges.legend(handles=[p_handle, t_handle], loc='upper right', fontsize=6)

        conc_scens = [names[2], names[3], names[4], names[5]]
        conc_inits = [names[0], names[1], names[0], names[1]]
        conc_colors = ['tab:cyan', 'tab:grey', 'orchid', 'chocolate']

        for i, (scen, init, color) in enumerate(zip(conc_scens, conc_inits, conc_colors)):
            data = utils.get_max_waiwhetu_concentration_over_time(scen, init)
            ax_concentrations.plot(data['times'], data["max_conc"], label=f"Scenario {i+3}", color=color)

        ax_concentrations.legend(loc='upper left', fontsize=6)

        ax_concentrations.set_xlabel('Time [years]')
        ax_concentrations.yaxis.tick_right()
        ax_concentrations.yaxis.set_label_position("right")
        ax_concentrations.set_ylabel(r'C$_{\mathrm{max}}$ Upper Waiwhetu [PSU]')

        grid = discretization.get_grid_w_ibound()
        pockmark_colors = ['tab:red', 'tab:orange', 'tab:brown', 'tab:purple', 'tab:blue', 'tab:pink']
        for i, example in enumerate(examples):
            color = pockmark_colors[example[0]]
            plot_pockmark(names[example[0]], example[1], grid, color='black', ax_plan=None, ax_section=axs_pockmarks[i], ylim=-10)
            if i > 0:
                axs_pockmarks[i].set_yticklabels([])
                axs_pockmarks[i].set_ylabel('')
            axs_pockmarks[i].text(0.025, 0.975, f"Scenario {example[0]+1}", fontsize=6, transform=axs_pockmarks[i].transAxes, ha='left', va='top')
            axs_pockmarks[i].text(0.975, 0.975, f"Pockmark {example[1]+1}", fontsize=6, transform=axs_pockmarks[i].transAxes, ha='right', va='top')

        norm = mpl.colors.Normalize(vmin=0, vmax=35)
        mappable = plt.cm.ScalarMappable(cmap='viridis', norm=norm)
        cb = plt.colorbar(mappable, cax=ax_cb, label='C [PSU]', orientation='horizontal', fraction=1)

        plt.savefig(f"/home/connor/PycharmProjects/pockmarks/figures/changes_and_examples_{set}.png", dpi=600)

    elif version == 1:

        f = plt.figure(figsize=(6, 8))
        gs = f.add_gridspec(5, 3, height_ratios=[1, 1, 1, 0, 0.1], hspace=0.5, wspace=0.1)
        plt.rcParams.update({'font.size': 8})
        ax_discharges = f.add_subplot(gs[0, :])
        axs_pockmarks = [f.add_subplot(gs[1+i//3, i%3]) for i in range(6)]
        ax_cb = f.add_subplot(gs[-1, :])

        pockmarks_discharge = []
        for name in names:
            temp_pockmarks = utils.get_final_pockmark_results(name)
            pockmarks_discharge.append(temp_pockmarks['discharge'].sum())

        bar_colors = ['tab:green', 'tab:cyan', 'orchid']

        X_axis = 2*np.arange(int(len(names)/3))

        ax_discharges.bar(X_axis-0.6, pockmarks_discharge[:3], 0.6, color=bar_colors, alpha=1, label='Summer Average')
        ax_discharges.bar(X_axis, pockmarks_discharge[3:6], 0.6, color=bar_colors, alpha=0.75, label='Summer Minimum')
        ax_discharges.bar(X_axis+0.6, pockmarks_discharge[6:], 0.6, color=bar_colors, alpha=0.5, label='Seawater Intrusion')

        ax_discharges.set_xticks(X_axis)
        ax_discharges.set_xticklabels(['Present-day', '2070', '2130'])
        ax_discharges.set_ylabel(r'Pockmark discharge [m$^3$/day]')

        handle0 = mpatches.Patch(color='tab:grey', label='Summer Average')
        handle1 = mpatches.Patch(color='tab:grey', alpha=0.75, label='Summer Minimum')
        handle2 = mpatches.Patch(color='tab:grey', alpha=0.5, label='Seawater Intrusion')
        ax_discharges.legend(handles=[handle0, handle1, handle2], loc='upper right', fontsize=6)

        times = ['Present-day', '2070', '2130']
        grid = discretization.get_grid_w_ibound()
        for i, example in enumerate(examples[:3]):
            plot_pockmark(names[example[0]], example[1], grid, color='black', ax_plan=None, ax_section=axs_pockmarks[i],
                          ylim=-10, vectors=False, contours=True)
            if i > 0:
                axs_pockmarks[i].set_yticklabels([])
                axs_pockmarks[i].set_ylabel('')
            axs_pockmarks[i].text(0.025, 0.975, f"{times[i]}", fontsize=6, transform=axs_pockmarks[i].transAxes, ha='left', va='top')
            axs_pockmarks[i].text(0.975, 0.975, f"Average", fontsize=6, transform=axs_pockmarks[i].transAxes, ha='right', va='top')

        plt.savefig(f"/home/connor/PycharmProjects/pockmarks/figures/changes_and_examples_alternate_{set}.png", dpi=600)

def plot_examples_alternative_hypothesis(plan=[('1average', 1), ('hk1average', 1), ('c1average', 1)],
                                         section=[('c1average', 5), ('c9slr2_swi', 5)]):

    f = plt.figure(figsize=(6, 5))
    gs = f.add_gridspec(2, 6, height_ratios=[1, 1.5], hspace=0.5, wspace=0.1)
    plt.rcParams.update({'font.size': 8})
    axs_plan = [f.add_subplot(gs[0, int(2*i):int(2*i+2)]) for i in range(3)]
    axs_section = [f.add_subplot(gs[1, int(3*i):int(3*i+3)]) for i in range(2)]

    grid = discretization.get_grid_w_ibound()
    labels = ['Thin', 'Patchy', 'Conduits']
    for i, (name, pockmark) in enumerate(plan):
        plot_pockmark(name, pockmark, grid, ax_plan=axs_plan[i], ax_section=None, vectors=False, contours=False,
                      section_line=False, discharge_on_plan=True)
        if i > 0:
            axs_plan[i].set_yticklabels([])
            axs_plan[i].set_ylabel('')
        axs_plan[i].text(0.025, 0.975, f"{labels[i]}", fontsize=6, transform=axs_plan[i].transAxes, ha='left', va='top')
        axs_plan[i].text(0.975, 0.975, f"Average", fontsize=6, transform=axs_plan[i].transAxes, ha='right', va='top')

    labels = ['Average', 'Seawater Intrusion']
    for i, (name, pockmark) in enumerate(section):
        plot_pockmark(name, pockmark, grid, ax_plan=None, ax_section=axs_section[i], vectors=False, contours=False,
                      section_line=False, discharge_on_plan=False)
        if i > 0:
            axs_section[i].set_yticklabels([])
            axs_section[i].set_ylabel('')

        axs_section[i].set_ylim(-35, -10)

        axs_section[i].text(0.025, 0.975, f"Conduits", fontsize=6, transform=axs_section[i].transAxes, ha='left', va='top')
        axs_section[i].text(0.975, 0.975, f"{labels[i]}", fontsize=6, transform=axs_section[i].transAxes,
                         ha='right', va='top')

    for ax in axs_section:
        ax.fill_between(ax.get_xlim(), ax.get_ylim()[0], ax.get_ylim()[1], color='white', alpha=1, zorder=-1)

    plt.savefig(f"/home/connor/PycharmProjects/pockmarks/figures/examples_alternate_hypothesis_plan.png", dpi=600, transparent=True)

def plot_salinity_comparison(names, colors=['tab:green', 'orchid', 'chocolate'], ax=None, errors_on_modelled_conduits=False, labels=None, legend=True):
    if ax == None:
        f, ax = plt.subplots(1, 1, figsize=(2, 2))

    data = pd.read_csv("/home/connor/PycharmProjects/pockmarks/data/observations.csv", index_col=0,
                       dtype={'name': str, 'concentration_average': float, 'discharge_min': float,
                              'discharge_max': float})
    if labels == None:
        labels = ['Thin', 'Patchy', 'Conduits']

    for color, name, label in zip(colors, names, labels):
        modelled = utils.get_salinity_and_flux_at_measurement_cells(name)
        modelled.index = data.index
        for measurement in ['mc2-1', 'mc2-2', 'mc2-3']:
            cent = modelled.loc[measurement, 'conc_center']
            if not "Conduits" in label or errors_on_modelled_conduits :
                xerr = [cent-[modelled.loc[measurement, 'conc_min']], [modelled.loc[measurement, 'conc_max']-cent]]
            else:
                xerr = None
            ax.errorbar([cent], [data.loc[measurement, 'concentration_average']],
                        xerr=xerr,
                        fmt=' ', marker='o', fillstyle=None,
                        c=color,  markersize=4, zorder=4, capsize=2, linewidth=1, markerfacecolor='w', label=label if measurement == 'mc2-1' else None)

    ax.plot(np.linspace(0,35), np.linspace(0,35), color='black', ls='--', lw=0.5)
    ax.set_box_aspect(1)

    ax.text(0.85, 0.74, f"{'MC2-3'}", fontsize=6, transform=ax.transAxes, ha='left', va='center')
    ax.text(0.85, 0.6125, f"{'MC2-1'}", fontsize=6, transform=ax.transAxes, ha='left', va='center')
    ax.text(0.85, 0.145, f"{'MC2-2'}", fontsize=6, transform=ax.transAxes, ha='left', va='center')
    ax.set_ylabel('Observed salinity [PSU]')
    ax.set_xlabel('Modelled salinity [PSU]')

    if legend:
        ax.legend(fontsize=6)

def plot_flux_comparison(names, colors=['tab:green', 'orchid', 'chocolate'], ax=None, errors_on_modelled_conduits=False, labels = None, legend=False):
    if ax == None:
        f, ax = plt.subplots(1, 1, figsize=(2, 2))

    data = pd.read_csv("/home/connor/PycharmProjects/pockmarks/data/observations.csv", index_col=0,
                       dtype={'name': str, 'concentration_average': float, 'discharge_min': float,
                              'discharge_max': float})

    if labels == None:
        labels = ['Thin', 'Patchy', 'Conduits']

    for color, name, label in zip(colors, names, labels):
        modelled = utils.get_salinity_and_flux_at_measurement_cells(name)
        modelled.index = data.index
        for measurement in ['sv1', 'sv6']:
            xcent = modelled.loc[measurement, 'discharge_center']
            ycent = (data.loc[measurement, 'discharge_min']+data.loc[measurement, 'discharge_max'])/2
            if (not "Conduits" in label) or errors_on_modelled_conduits:
                xerr = [xcent - [modelled.loc[measurement, 'discharge_min']],
                              [modelled.loc[measurement, 'discharge_max'] - xcent]]
            else:
                xerr = None
            yerr = [[ycent - data.loc[measurement, 'discharge_min']], [data.loc[measurement, 'discharge_max'] - ycent]]
            ax.errorbar([xcent], [ycent],
                        xerr=xerr, yerr=yerr,
                        fmt=' ', marker='o', fillstyle=None,
                        c=color, markersize=4, zorder=4, capsize=2, linewidth=1, markerfacecolor='w',
                        label=label if measurement == 'sv1' else None)

    ax.plot(np.linspace(1e-6, 1e5), np.linspace(1e-6, 1e5), color='black', ls='--', lw=0.5)
    ax.set_box_aspect(1)
    ax.set_ylabel('Observed flux [m/day]')
    ax.set_xlabel('Modelled flux [m/day]')
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_ylim(1e-6, 1e5)
    ax.set_xlim(1e-6, 1e5)

    ax.text(0.65, 0.9, f"{'SV6'}", fontsize=6, transform=ax.transAxes, ha='left', va='center')
    ax.text(0.65, 0.85, f"{'SV1'}", fontsize=6, transform=ax.transAxes, ha='left', va='center')

    if legend:
        h, l = ax.get_legend_handles_labels()
        ax.legend(fontsize=6)

    return ax

def plot_dropping_head_fluxes(names, colors, labels, ax, pockmark=5):
    t_series = utils.get_pockmark_discharge_and_recharge(names[0])[f"{pockmark}"].index

    for name, color, label in zip(names, colors, labels):
        data = utils.get_pockmark_discharge_and_recharge(name)
        data = data[f"{pockmark}"].loc[t_series]
        ax.plot([i for i in range(15)], data['discharge'], color=color, ls=":", label=label)
        ax.plot([i for i in range(15)], data['recharge'], color=color, label=label)
    ax.set_yscale('symlog', linthresh=1e-2)
    ax.set_ylabel(f'Pockmark {pockmark} Flux [m$^3$/day]')
    ax.set_xlabel(r'$h_{\mathrm{on}}$ [m]')
    ax.set_xticks([i for i in range(15)][::2])
    ax.set_xticklabels([2.75-0.25*(i+1) for i in range(15)][::2])
    ax.set_box_aspect(1)
    discharge = Line2D([], [], color='grey', ls=':', label='Discharge')
    recharge = Line2D([], [], color='grey', label='Intrusion')
    ax.legend(handles=[discharge, recharge], loc='upper right', fontsize=6)


def plot_alternate_hypotheses_pockmarks_and_model_vs_observed(names, salinity_number=3, discharge_number=0):
    f = plt.figure(figsize=(8, 6.5))
    gs = f.add_gridspec(6, 3, width_ratios=[1,1,1.6], hspace=0.25, wspace=0.5)
    plt.rcParams.update({'font.size': 8})
    ax_salinity = f.add_subplot(gs[:3, -1])
    ax_discharges = f.add_subplot(gs[3:, -1])
    axs_plan = []
    axs_section = []
    for i, name in enumerate(names):
        axs_plan.append(f.add_subplot(gs[2*i:2*i+2, 0]))
        axs_section.append(f.add_subplot(gs[2*i:2*i+2, 1]))

    plot_salinity_comparison(names, ax=ax_salinity)
    ax_salinity.text(0.975, 1.025, f"(g)", fontsize=8, transform=ax_salinity.transAxes, ha='right', va='bottom')
    ax_discharges.text(0.975, 1.025, f"(h)", fontsize=8, transform=ax_discharges.transAxes, ha='right', va='bottom')
    plot_flux_comparison(names, ax=ax_discharges)

    grid = discretization.get_grid_w_ibound()
    for name, ax_plan, ax_section, letter in zip(names, axs_plan, axs_section, [['a', 'b'], ['c', 'd'], ['e', 'f']]):

        plot_pockmark(name, discharge_number, grid, ax_plan=None, ax_section=ax_section, vectors=False, contours=False,
                      section_line=False, discharge_on_plan=False)

        # add marker at 60.5, -11
        ax_section.scatter(60.5, -11+0.3, marker='^', color='blue', lw=0.75, s=15, facecolors='none')

        plot_pockmark(name, salinity_number, grid, ax_plan=ax_plan, ax_section=None, vectors=False, contours=False,
                      section_line=False, discharge_on_plan=True)

        ax_plan.scatter(1.758e6+298.5, 5.432e6+701.8-2, marker='v', color='darkgoldenrod', lw=0.75, s=15, facecolors='none')


        ax_plan.text(0.975, 1.025, f"({letter[0]})", fontsize=8, transform=ax_plan.transAxes, ha='right', va='bottom')
        ax_section.text(0.975, 1.025, f"({letter[1]})", fontsize=8, transform=ax_section.transAxes, ha='right', va='bottom')
        if not name == names[-1]:
            ax_plan.set_xticklabels([])
            ax_plan.set_xlabel('')
            ax_section.set_xticklabels([])
            ax_section.set_xlabel('')

    labels = ['Thin', 'Patchy', 'Conduits']
    for i in range(len(axs_plan)):
        axs_plan[i].text(0.975, 0.975, f"{labels[i]}", fontsize=6, transform=axs_plan[i].transAxes, ha='right', va='top')
        axs_section[i].text(0.975, 0.975, f"{labels[i]}", fontsize=6, transform=axs_section[i].transAxes, ha='right', va='top')


        axs_section[i].text(0.025, 0.975, f"Pockmark 1", fontsize=6, transform=axs_section[i].transAxes, ha='left', va='top')
        axs_plan[i].text(0.025, 0.975, f"Pockmark 4", fontsize=6, transform=axs_plan[i].transAxes, ha='left', va='top')

    legend_elements = []
    legend_elements.append(Line2D([0], [0], marker='^', color='w', label='SV1', markerfacecolor='none', markeredgecolor='blue', markersize=4, lw=0.5))
    legend_elements.append(Line2D([0], [0], marker='v', color='w', label='MC2-1', markerfacecolor='none', markeredgecolor='darkgoldenrod', markersize=4, lw=0.5))
    axs_plan[0].legend(handles=legend_elements, bbox_to_anchor=(0, 1.15, 1, 0.2), loc="lower left", ncol=2, fontsize=6)

    plt.tight_layout()
    plt.savefig(f"/home/connor/PycharmProjects/pockmarks/figures/alternate_hypotheses_pockmarks_and_model_vs_observed.png", dpi=600)
    # plt.show()

def plot_alternative_hypotheses_changes(names=None, pockmark_number=5, dropping_names=None):
    if names == None:
        names = [['1average', 'big1_swi', '3slr2'],
                 ['hk1ave_new', 'hkbig1_swi_n', "hkslr2_new"],
                 ['c1ave_new', 'cbig1_swi_n', "cslr2_new"]]
    if dropping_names == "test":
        dropping_names = ['dropping_t', 'hkdropping_t', 'cdropping']
    elif dropping_names == False:
        dropping_names = None
    else:
        dropping_names = ['dropping', 'hkdropping', 'cdropping']


    f = plt.figure(figsize=(8, 5))
    plt.rcParams.update({'font.size': 8})
    gs = gridspec.GridSpec(1, 2,  width_ratios=[1, 0.4], figure=f)
    gs_sections = gridspec.GridSpecFromSubplotSpec(3, 3, subplot_spec=gs[0], hspace=0.2, wspace=0.25)
    gs_results = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[1])

    axs_sections = [[f.add_subplot(gs_sections[i, j]) for j in range(3)] for i in range(3)]
    axs_results = [f.add_subplot(gs_results[i]) for i in range(2)]

    grid = discretization.get_grid_w_ibound()
    concepts = ["Thin", "Patchy", "Conduits"]
    scenarios = ["Present-day", "Depletion", " SLR 2130"]
    letters = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i']
    for i, (ax, name) in enumerate(zip(np.array(axs_sections).flatten(), np.array(names).flatten())):
        plot_pockmark(name, pockmark_number, grid, ax_plan=None, ax_section=ax, vectors=False, contours=True,
                      section_line=False, discharge_on_plan=False)

        ax.set_ylim(-35, -10)

        if i % 3 != 0:
            ax.set_yticklabels([])
            if not i == 4:
                ax.set_ylabel('')
        if i < 6:
            ax.set_xticklabels([])
            ax.set_xlabel('')
        if not i == 7:
            ax.set_xlabel('')

        ax.text(0.025, 0.975, f"{concepts[i//3]}", fontsize=6, transform=ax.transAxes, ha='left', va='top')
        ax.text(0.975, 0.975, f"{scenarios[i%3]}", fontsize=6, transform=ax.transAxes, ha='right', va='top')
        ax.text(0.975, 1.025, f"({letters[i]})", fontsize=8, transform=ax.transAxes, ha='right', va='bottom')
    X_axis = np.arange(2)
    for hypoth, offset, color in zip(np.array(names)[:, [0,-1]], [-0.2, 0, 0.2], ['tab:green', 'orchid', 'chocolate']):
        pockmarks_discharge = []
        for scena in hypoth:
            temp_pockmarks = utils.get_final_pockmark_results(scena)
            pockmarks_discharge.append(temp_pockmarks['discharge'].sum())
        axs_results[1].bar(X_axis + offset, pockmarks_discharge, 0.2, color=color)

    names = ["Thin", "Patchy", "Conduits"]
    handles = []
    for color, name in zip(['tab:green', 'orchid', 'chocolate'], names):
        handles.append(mpatches.Patch(color=color, label=name))

    axs_results[1].set_yscale('log')
    axs_results[1].set_xticks(X_axis)
    axs_results[1].set_xticklabels(["Present-day", "2130"])
    axs_results[1].set_xlabel("Year")
    axs_results[1].set_ylabel(r'Total Pockmark Discharge [m$^3$/day]')
    axs_results[1].set_box_aspect(1)
    axs_results[1].legend(handles=handles, loc='upper right', fontsize=6)
    axs_results[1].set_ylim(0, 100)

    if dropping_names != None:
        plot_dropping_head_fluxes(dropping_names, ["tab:green", "orchid", "chocolate"], ["Thin", "Patchy", "Conduit"], axs_results[0])
        axs_results[0].text(0.975, 1.025, f"(j)", fontsize=8, transform=axs_results[0].transAxes, ha='right',
                            va='bottom')
        axs_results[1].text(0.975, 1.025, f"(k)", fontsize=8, transform=axs_results[1].transAxes, ha='right',
                            va='bottom')
    else:
        axs_results[1].text(0.975, 1.025, f"(j)", fontsize=8, transform=axs_results[1].transAxes, ha='right',
                            va='bottom')
        axs_results[0].set_axis_off()
        lines = Line2D([], [], linewidth=0.5, color='red', linestyle='dashed', label="5% salinity")
        axs_results[0].legend(handles=[lines], loc='center left', fontsize=6)
        
    plt.tight_layout()
    plt.savefig(f"/home/connor/PycharmProjects/pockmarks/figures/alternate_hypotheses_changes.png", dpi=600)
    plt.show()

def plot_conduit_and_patches_locations():
    f, ax = plt.subplots(figsize=(5, 4))
    plt.rcParams.update({'font.size': 7})
    base = rasterio.open("/home/connor/PycharmProjects/pockmarks/map/data/basemap.tif")
    model_area = gpd.read_file('/home/connor/PycharmProjects/pockmarks/data/model_area0.shp')
    extent = model_area.total_bounds
    show(base, ax=ax, alpha=0.7, zorder=0)
    ax.set_xlim(extent[0], extent[2])
    ax.set_ylim(extent[1], extent[3])
    ax.set_xlabel('Nztm x [m]')
    ax.set_ylabel('Nztm y [m]')
    patches = gpd.read_file('/home/connor/PycharmProjects/pockmarks/data/refinement_areas.shp', crs='EPSG:2193')
    patches.to_crs('EPSG:2193')
    patches.plot(ax=ax, edgecolors='black', color="orchid", zorder=10, alpha=1, lw=0.75, label="Patches")
    points = np.genfromtxt("/home/connor/PycharmProjects/pockmarks/data/coordinates_from_harding_appendix.csv",
                           delimiter=',', skip_header=1)[:,1:]
    arbitrary = gpd.read_file('/home/connor/PycharmProjects/pockmarks/data/abritrary_pockmark_locations.shp', crs='EPSG:2193')
    arbitrary.to_crs('EPSG:2193')
    arbitrary.plot(ax=ax, edgecolors='black', color="chocolate", zorder=12, alpha=1, lw=0.75, markersize=15, label="Conduits")

    transformer = Transformer.from_crs("EPSG:4326", "EPSG:2193")

    for i, coord in enumerate(points):
        if i < 2:
            coord[0] += +0.0016666
        nztm = transformer.transform(*coord)
        ax.scatter(nztm[1], nztm[0], edgecolors='black', color="chocolate", zorder=12, alpha=1, lw=0.75, s=15)

    label_sides = ['left', 'left', 'left', 'right', 'right', 'right', 'right', 'right']
    numbers = [6, 5, 3, 4, 2, None, 1, None]
    for i, pockmark in enumerate(patches.geometry):
        if numbers[i] is None:
            continue
        elif label_sides[numbers[i]-1] == 'left':
            ax.text(pockmark.bounds[0] - 40, pockmark.centroid.y, numbers[i], fontsize=6, ha='right', va='center')
        else:
            ax.text(pockmark.bounds[2] + 40, pockmark.centroid.y, numbers[i], fontsize=6, ha='left', va='center')

    patch = mpatches.Patch(edgecolor='black', lw=0.75, facecolor='orchid', label='Patches')
    h, l = ax.get_legend_handles_labels()
    h.append(patch)
    ax.legend(handles=h, fontsize=6, loc='upper left')
    plt.savefig(f"/home/connor/PycharmProjects/pockmarks/figures/conduit_and_patches_locations.png", dpi=600)

def plot_multi_model_comparison():
    f, axs = plt.subplots(1, 2, figsize=(7, 3.5))
    plt.rcParams.update({'font.size': 8})
    plot_salinity_comparison(['1average', 'hk1ave_new', 'c1ave_new', 'hhhk_average'],
                             colors=['tab:green', 'orchid', 'chocolate', 'purple'],
                             labels=[r'Thin ($K_v = 1.8 \times 10^{-5}$ m/day)', r'Patchy ($K_v = 3.6 \times 10^{-5}$ m/day)', 'Conduits (K_v = 10 m/day)', r'Patchy ($K_v = 1.8 \times 10^{-3}$ m/day)'],
                            ax=axs[0], legend=False)
    plot_flux_comparison(['1average', 'hk1ave_new', 'c1ave_new', 'hhhk_average'],
                             colors=['tab:green', 'orchid', 'chocolate', 'purple'],
                             labels=[r'Thin ($K_v = 1.8 \times 10^{-5}$ m/day)', r'Patchy ($K_v = 3.6 \times 10^{-5}$ m/day)', 'Conduits ($K_v = 10$ m/day)', r'Patchy ($K_v = 1.8 \times 10^{-3}$ m/day)'],
                         ax=axs[1], legend=True)
    plt.tight_layout()
    plt.savefig(f"/home/connor/PycharmProjects/pockmarks/figures/multi_model_comparison_2.png", dpi=600)

def plot_longer_swi(models=None, names=None):
    if models == None:
        models = ['bl1_swi', 'hkbl1_swi_n', 'cbl1_swi_n',"hhhkbl1_swi_n"]
        names = [f'Thin {r'($K_v = 1.8 \times 10^{-5}$ m/day)'}', r'Patchy ($K_v = 3.6 \times 10^{-5}$ m/day)', 'Conduits ($K_v = 10$ m/day)', r'Patchy ($K_v = 1.8 \times 10^{-3}$ m/day)']
    f, axs = plt.subplots(1, 4, figsize=(7.5, 2.5))
    plt.rcParams.update({'font.size': 8})
    grid = discretization.get_grid_w_ibound()
    for ax, model, name in zip(axs, models, names):
        plot_pockmark(model, 5, grid, ax_plan=None, ax_section=ax, label_sections=False,
                      vectors=False, contours=True, section_line=True, discharge_on_plan=False, ylim=(-35, -5))
        ax.text(0.025, 0.975, f"{name}", fontsize=6, transform=ax.transAxes, ha='left', va='top')
        # ax.text(0.975, 0.975, f"SWI (10 years)", fontsize=6, transform=ax.transAxes, ha='right', va='top')

    axs[1].set_xlabel('')
    axs[2].set_xlabel('')
    axs[3].set_xlabel('')
    axs[1].set_ylabel('')
    axs[2].set_ylabel('')
    axs[3].set_ylabel('')
    axs[1].set_yticklabels([])
    axs[2].set_yticklabels([])
    axs[3].set_yticklabels([])

    plt.tight_layout()
    plt.savefig(f"/home/connor/PycharmProjects/pockmarks/figures/longer_swi.png", dpi=600)



if __name__=="__main__":
    # names = ['1average', '2slr1', '3slr2', '4minimum', '5slr1_min', '6slr2_min', '7swi', '8slr1_swi', '9slr2_swi']
    # examples  = ((0, 5), (1, 5), (2, 5))
    # plot_changes_and_examples("base", names=names, examples=examples, version=1)
    # plot_examples_alternative_hypothesis()
    # plot_aquitard_thickness()
    # plot_alternate_hypotheses_pockmarks_and_model_vs_observed(['1average', 'hk1ave_new', 'c1ave_new'])
    # plot_alternate_hypotheses_pockmarks_and_model_vs_observed(['1average'], 1)
    # plot_alternative_hypotheses_changes(pockmark_number=5, dropping_names=False)
    # plot_conduit_and_patches_locations()
   # plot_multi_model_comparison()
    # plot_longer_swi()
    plot_boundary_conditions()


