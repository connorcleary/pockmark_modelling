import numpy as np
from mpl_toolkits.mplot3d.proj3d import transform
from scipy.stats import norm, beta
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

def plot_geological_hypothese():
    fig, ax = plt.subplots(1, 1, figsize=(7, 2))


    x1 = np.linspace(0, 16, 100)
    x2 = np.linspace(16, 32, 100)
    x3 = np.linspace(32, 48, 100)

    dist = norm(0, 1.75)
    y1 = -dist.pdf(x1-8)
    y2 = -dist.pdf(x2-24)
    y3 = -dist.pdf(x3-40)

    ax.plot(x1, y1, c='k', lw=0.75)
    ax.plot(x2, y2, c='k', lw=0.75)
    ax.plot(x3, y3, c='k', lw=0.75)

    ax.axvline(16, c='k', ls="--", lw=0.75)
    ax.axvline(32, c='k', ls="--", lw=0.75)


    ax.set_ylim(-0.4, 0.2)
    ax.set_xlim(0, 48)

    ax.text(0.5, 0.185, "Thin Aquitard", ha='left', va='top', fontsize=6)
    ax.text(16.5, 0.185, "Patchy Aquitard", ha='left', va='top', fontsize=6)
    ax.text(32.5, 0.185, "Aquitard Conduits", ha='left', va='top', fontsize=6)

    ax.fill_between(np.concatenate([x1, x2, x3]), np.concatenate([y1, y2, y3]),
                    0.2*np.ones((300)), color='skyblue', alpha=0.3, ec='none')
    ax.fill_between(np.concatenate([x1, x2, x3]), np.concatenate([y1, y2, y3]),
                    -0.3*np.ones((300)), color='silver', alpha=1, ec='none')

    x_patch = np.linspace(20, 28, 100)
    ax.fill_between(x_patch,-dist.pdf(x_patch-24), -0.3*np.ones((100)), color='gainsboro', alpha=1, ec='none')

    x_center_conduit = np.linspace(39.9, 40.1, 10)
    y_center_conduit = -dist.pdf(x_center_conduit-40)
    ax.fill_between(x_center_conduit, y_center_conduit, -0.3*np.ones((10)), color='white', alpha=1, ec='white')

    x_edge_conduit = np.linspace(46, 46.2, 10)
    y_edge_conduit = -dist.pdf(x_edge_conduit-40)
    ax.fill_between(x_edge_conduit, y_edge_conduit, -0.3*np.ones((10)), color='white', alpha=1, ec='white')

    colors = [mcolors.to_rgba('silver'),
              mcolors.to_rgba('gainsboro'),
              mcolors.to_rgba('white')]

    colors_temp = [mcolors.to_hex(color, keep_alpha=True) for color in colors]

    cmap = mcolors.LinearSegmentedColormap.from_list("custom", colors_temp, N=3)
    mappable = plt.cm.ScalarMappable(cmap=cmap)
    cb = fig.colorbar(mappable, ax=ax, orientation='vertical', location='right', pad=0.02, shrink=1)
    cb.set_ticks([0, 1])
    cb.set_ticklabels(['Low K', "High K"], fontsize=6)

    ax.set_xticks([])
    ax.set_yticks([])

    ax.fill_between(ax.get_xlim(), ax.get_ylim()[0], ax.get_ylim()[1], color='white', alpha=1, zorder=-1)
    fig.tight_layout()
    fig.savefig("/home/superuser/objective_3/figures/geological_hypotheses.png", dpi=600, transparent=True)

def plot_discharge_typologies():
    fig, ax = plt.subplots(1, 1, figsize=(7, 2))

    ax.set_xlim(0, 54)
    ax.set_ylim(-0.45, 0.2)

    dist_deep = norm(0, 1.6)
    dist_shallow = norm(0, 1)

    x1 = np.linspace(0, 12, 100)
    y1 = -dist_deep.pdf(x1 - 6)

    x2 = np.linspace(12, 16, 100)
    y2 = -0.25*dist_shallow.pdf(x2 - 16)

    x3 = np.linspace(16, 20, 100)
    y3 = -0.25*dist_shallow.pdf(0)*np.ones(100)

    x4 = np.linspace(20, 24, 100)
    y4 = -0.25*dist_shallow.pdf(x4 - 20)

    x5 = np.linspace(24, 30, 100)
    y5 = np.zeros_like(x5)

    x6 = np.linspace(30, 42, 100)
    y6 = -dist_deep.pdf(x6 - 36)

    x7 = np.linspace(42, 46, 100)
    y7 = -0.25 * dist_shallow.pdf(x7 - 46)

    x8 = np.linspace(46, 50, 100)
    y8 = -0.25 * dist_shallow.pdf(0) * np.ones(100)

    x9 = np.linspace(50, 54, 100)
    y9 = -0.25 * dist_shallow.pdf(x9 - 50)

    x = np.concatenate([x1, x2, x3, x4, x5, x6, x7, x8, x9])
    y = np.concatenate([y1, y2, y3, y4, y5, y6, y7, y8, y9])
    ax.plot(x, y, c='k', lw=0.75)

    ax.fill_between(x, y, 0.2*np.ones_like(x), color='skyblue', alpha=0.3, ec='none')
    ax.fill_between(x, y, -0.35*np.ones_like(x), color='silver', alpha=1, ec='none')

    x_between_conduit = np.linspace(26.9, 27.1, 10)
    y_between_conduit = 0*np.ones((10))
    ax.fill_between(x_between_conduit, y_between_conduit, -0.35*np.ones((10)), color='white', alpha=1, ec='white')

    x_center_conduit = np.linspace(35.9, 36.1, 10)
    y_center_conduit = -dist_deep.pdf(x_center_conduit-36)
    ax.fill_between(x_center_conduit, y_center_conduit, -0.35*np.ones((10)), color='white', alpha=1, ec='white')

    x_edge_conduit = np.linspace(45.9, 46.1, 10)
    y_edge_conduit =  -0.25 * dist_shallow.pdf(x_edge_conduit - 46)
    ax.fill_between(x_edge_conduit, y_edge_conduit, -0.35*np.ones((10)), color='white', alpha=1, ec='white')

    gradient = np.gradient(100*y, x)

    start = 30
    end = 71
    stride = 10

    cmap = plt.cm.get_cmap('viridis')
    colornorm = mcolors.Normalize(vmin=0, vmax=35)
    color = cmap(colornorm(10))

    ax.quiver(x[start:end:stride], y[start:end:stride]-0.01,
              -gradient[start:end:stride]/np.sqrt(gradient[start:end:stride]**2+1), 1/np.sqrt(gradient[start:end:stride]**2+1),
              color=color, pivot='tip', angles='xy', scale=20)

    start = 160
    end = 210
    stride = 25

    color = cmap(colornorm(15))
    ax.quiver(x[start:end:stride], y[start:end:stride]-0.01,
              -gradient[start:end:stride]/np.sqrt(gradient[start:end:stride]**2+1), 1/np.sqrt(gradient[start:end:stride]**2+1),
              color=color, pivot='tip', angles='xy', scale_units='inches', scale=5)

    start = 219
    end = 280
    stride = 20

    color = cmap(colornorm(20))
    ax.quiver(x[start:end:stride], y[start:end:stride]-0.01,
              -gradient[start:end:stride]/np.sqrt(gradient[start:end:stride]**2+1), 1/np.sqrt(gradient[start:end:stride]**2+1),
              color=color, pivot='tip', angles='xy', scale_units='inches', scale=10)

    start = 310
    end = 360
    stride = 25

    color = cmap(colornorm(15))
    ax.quiver(x[start:end:stride], y[start:end:stride] - 0.01,
              -gradient[start:end:stride] / np.sqrt(gradient[start:end:stride] ** 2 + 1),
              1 / np.sqrt(gradient[start:end:stride] ** 2 + 1),
              color=color, pivot='tip', angles='xy', scale_units='inches', scale=5)

    color = cmap(colornorm(5))
    ax.quiver(x[450], y[450] + 0.01,
              -gradient[450] / np.sqrt(gradient[450] ** 2 + 1),
              1 / np.sqrt(gradient[450] ** 2 + 1),
              color=color, pivot='tail', angles='xy', scale_units='inches', scale=3)

    color = cmap(colornorm(5))
    ax.quiver(x[550], y[550] + 0.01,
              -gradient[550] / np.sqrt(gradient[550] ** 2 + 1),
              1 / np.sqrt(gradient[550] ** 2 + 1),
              color=color, pivot='tail', angles='xy', scale_units='inches', scale=3)

    color = cmap(colornorm(5))
    ax.quiver(x[-202], y[-202] + 0.01,
              -gradient[-202] / np.sqrt(gradient[-202] ** 2 + 1),
              1 / np.sqrt(gradient[-202] ** 2 + 1),
              color=color, pivot='tail', angles='xy', scale_units='inches', scale=3, zorder=10)

    ax.set_xticks([])
    ax.set_yticks([])

    fig.savefig("/home/superuser/objective_3/figures/discharge_typologies.png", dpi=600)


if __name__ == "__main__":
    plot_geological_hypothese()
    # plot_discharge_typologies()