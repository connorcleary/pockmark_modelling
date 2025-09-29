import discretization
import pyvista as pv
from flopy.export.vtk import Vtk
from matplotlib.colors import ListedColormap
import numpy as np
from matplotlib import colors
import gemgis as gg
from scripts.discretization import leapfrog_name
import os
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

def get_vtk_from_leapfrog(name):
    if os.path.exists(f"/home/superuser/objective_3/data/leapfrog/{name}.vtk"):
        return pv.read(f"/home/superuser/objective_3/data/leapfrog/{name}.vtk")
    tempDict = gg.raster.read_msh(f"/home/superuser/objective_3/data/leapfrog/{name}.msh")
    tempMesh = gg.visualization.create_polydata_from_msh(tempDict)
    tempMesh.save(f"/home/superuser/objective_3/data/leapfrog/{name}.vtk")
    return tempMesh

def plot_geology(from_leapfrog_mesh=False):

    colors = ["goldenrod", "cornflowerblue", "slategray", "mediumseagreen"]
    names = ["Upper Waiwhetu", "Petone Marine Lower", "Petone Marine Distal", "Petone Marine Upper"]
    names_with_type = ["Upper Waiwhetu (gravel)", "Petone Marine Lower (sand)", "Petone Marine Distal (silt)", "Petone Marine Upper (sand)"]

    if not from_leapfrog_mesh:
        grid = discretization.get_grid(leapfrog_name='grid_w_seafloor')
        zones = discretization.get_zone_array(leapfrog_name='grid_w_seafloor')
        vtk = Vtk(modelgrid=grid, vertical_exageration=15)
        vtk.add_array(zones, "zones")
        grid = vtk.to_pyvista()
    # grid.plot(scalars="zones", show_edges=False, show_grid=False, notebook=False)

        mapping = np.linspace(grid["zones"].min(), grid["zones"].max(), 256)
        newcolors = np.empty((256, 4))
        newcolors[mapping >= 4] = colors.to_rgba(colors[0])
        newcolors[mapping < 4] = colors.to_rgba(colors[1])
        newcolors[mapping < 3] = colors.to_rgba(colors[2])
        newcolors[mapping < 2] = colors.to_rgba(colors[3])
    # Make the colormap from the listed colors
        my_colormap = ListedColormap(newcolors)

        grid.plot(scalars="zones", cmap=my_colormap, show_edges=False, show_grid=False, notebook=False, opacity=0.5)

    else:
        p = pv.Plotter(notebook=False, window_size=[800,600], off_screen=True)

        for name, color in zip(names, colors):
            mesh = get_vtk_from_leapfrog(name)
            mesh.points[:, 2] = mesh.points[:, 2] * 10
            p.add_mesh(mesh, color=color)
        #

        labels = dict(ztitle='Depth [m]', xtitle='Nztm x [m]', ytitle='Nztm y [m]',
                      bold=False, use_3d_text=False, font_size=16, n_xlabels=3, n_ylabels=3,
                      show_xlabels=True, show_ylabels=True, show_zlabels=True)
        # p.show_grid(color=None, **labels)
        p.set_background(color='white')
        # p.set_scale(zscale=15)
        #p.add_axes(**labels)
        p.enable_parallel_projection()
        p.camera_position = 'yz'
        p.camera.azimuth = -45
        p.camera.elevation = 35

        light = pv.Light(intensity=0.5)
        light.set_direction_angle(35, -45)
        p.add_light(light)
        # p.show()
        p.screenshot(f"/home/superuser/objective_3/figures/geology.png", transparent_background=True, scale=3)

        plt.rcParams.update({'font.size': 8})

        f,ax = plt.subplots(figsize=(3,1))
        ax.axis('off')
        legend_elements = [Patch(facecolor=colors[3-i], edgecolor=colors[3-i], label=names_with_type[3-i]) for i in range(4)]
        ax.legend(handles=legend_elements, loc='center', fontsize=8, markerfirst=False)
        plt.savefig(f"/home/superuser/objective_3/figures/geology_legend.png", dpi=600)

        f = plt.figure(figsize=(8,6))
        ax = f.add_subplot(111, projection='3d')
        ax.set_zlim(-50, 5)
        ax.set_xlim(1757300, 1759300)
        ax.set_ylim(5431350, 5433900)
        ax.view_init(azim=-45, elev=35)
        ax.set_proj_type('ortho')
        ax.set_box_aspect((np.ptp([1757300, 1759300]), np.ptp([5431350, 5433900]), 10*np.ptp([-50, 5])))
        ax.set_xlabel('Nztm x [m]', labelpad=10)
        ax.set_ylabel('Nztm y [m]', labelpad=10)
        ax.set_zlabel('Depth [m]')


        f.savefig(f"/home/superuser/objective_3/figures/geology_grid.png", dpi=600)








if __name__ == '__main__':
    plot_geology(from_leapfrog_mesh=False)