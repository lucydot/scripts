import math

from six import next
from numpy import array
from itertools import cycle

from matplotlib import rc

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


default_colours = [[240, 163, 255], [0, 117, 220], [153, 63, 0], [76, 0, 92],
                   [66, 102, 0], [255, 0, 16], [157, 204, 0], [194, 0, 136],
                   [0, 51, 128], [255, 164, 5], [255, 255, 0], [255, 80, 5],
                   [94, 241, 242], [116, 10, 255], [153, 0, 0], [0, 153, 143],
                   [0, 92, 49], [43, 206, 72], [255, 204, 153], [148, 255, 181],
                   [143, 124, 0], [255, 168, 187], [128, 128, 128]]

default_fonts = ['Arial', 'Whitney Book', 'Helvetica',
                 'Liberation Sans', 'Andale Sans']


def get_tl_diagram(band_gap, defect_profiles, height=5, width=None, fonts=None,
                   chem_pot_label='', legend_on=True, ratio=(1, 1), emax=None,
                   x_min=None, x_max=None):
    """
    Plot transition level diagram
    Args:
        band_gap (float): The band gap of the material.
        defect_profiles (dict): Dictionary of {"symbol": profile}. The symbol
            would be things like "V_{Sn}" or "I_{Cl}". The profile is a numpy array
            of [fermi-levels, defect_formation_energies]. For example if your
            defect has a transition level 0.1 eV above the Fermi level with a
            formation energy of 2.1 eV and another transition level 1.5 eV
            above the Fermi level with a formation energy of 2.5 eV, the
            profile would be:

                profile = np.array([(-1, 0.1, 1.5, 3), (0.2, 1.5, 2.5, 2.5)])

            You will notice that I have added two more points to the profile,
            these are to ensure the line is drawn to below and above the valence
            band maximum and conduction band minimum, respecitvely. Please
            Slack me if this doesn't make sense.
        height (float): Printed figure height in inches
        width (float): Printed figure width in inches
        fonts (list of str): List of fonts
        legend_on (bool): Show legend
        ratio (2-tuple): Aspect ratio if width not specified
        emax (float): The maximum energy to display on the y-axis.
        x_min (float): The minimum energy to display on the x-axis. If this is
            less than zero then shading will indicate where the valence band
            maximum starts.
        x_max (float): The maximum energy to display on the x-axis. If this is
            greater than the band gap then shading will indicate where the
            conduction band minimum starts.
    """

    if width is None:
        width = height * ratio[0]/ratio[1]

    ticklabelsize = 22
    ticksize = 12
    linewidth = 1.5

    plt.figure(figsize=(width, height), facecolor="w", dpi=600)

    ax = plt.gca()
    ax.tick_params(width=linewidth, size=ticksize)
    ax.tick_params(which='major', size=ticksize, width=linewidth,
                   labelsize=ticklabelsize, pad=10)
    ax.tick_params(which='minor', size=ticksize/2, width=linewidth)

    ax.set_title(ax.get_title(), size=20)
    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(linewidth)

    labelsize = 22

    ax.set_xlabel(ax.get_xlabel(), size=labelsize)
    ax.set_ylabel(ax.get_ylabel(), size=labelsize)

    fonts = default_fonts if fonts is None else fonts + default_fonts

    rc('font', **{'family': 'sans-serif', 'sans-serif': fonts})
    rc('text', usetex=False)
    rc('pdf', fonttype=42)
    rc('mathtext', fontset='stixsans')

    rgb_colours = array(default_colours)/255.
    colours = cycle(rgb_colours)

    y_min = 0
    y_max = -float('inf')

    if x_min is None:
        x_min = 0.
    if x_max is None:
        x_max = band_gap

    for symbol, profile in defect_profiles:
        c = list(next(colours))
        ax.plot(profile[0], profile[1], marker='o',
                label=r"$\mathregular{{{}}}$".format(symbol),
                color=c, markeredgecolor=c, lw=2, markersize=12)
        y_max = max(profile[1]) if max(profile[1]) > y_max else y_max

    y_max = emax if emax else y_max * 1.1

    # Show colourful band edges
    if x_min != 0:
        ax.imshow([(0, 1), (0, 1)], cmap=plt.cm.Blues,
                  extent=(x_min, 0, y_min, y_max),
                  vmin=0, vmax=3,
                  interpolation="bicubic",
                  rasterized=True)

    if x_max != band_gap:
        ax.imshow([(1, 0), (1, 0)], cmap=plt.cm.Oranges,
                  extent=(band_gap, x_max, y_min, y_max),
                  vmin=0, vmax=3,
                  interpolation="bicubic",
                  rasterized=True)

    plt.xlim(x_min, x_max)
    plt.ylim((y_min, y_max))

    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    ax.set_aspect(((x1-x0)/(y1-y0)) * (height/width))

    # TODO: Scale these numbers based on axis ratio
    ax.xaxis.set_major_locator(ticker.MaxNLocator(4))
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    ax.yaxis.set_major_locator(ticker.MaxNLocator(4))
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))

    ax.set_ylabel('Formation Energy (eV)')
    ax.set_xlabel('Fermi Level (eV)')

    ax.tick_params(which='major', pad=10)
    ax.tick_params(which='minor', size=6)

    if legend_on:
        ncol = int(math.ceil(len(defect_profiles.keys())/10.))
        ax.legend(bbox_to_anchor=(1, 1), loc=2, frameon=False,
                  fontsize=20, numpoints=1, ncol=ncol)
    return plt