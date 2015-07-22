#! /usr/bin/env python

#
# Some default parameters for matplotlib
#

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.axes as ax

mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = ['Computer Modern Roman']
mpl.rcParams['text.usetex'] = True
mpl.rcParams['font.size'] = 28

mpl.rcParams['xtick.major.size'] = 10
mpl.rcParams['xtick.minor.size'] = 5
mpl.rcParams['ytick.major.size'] = 10
mpl.rcParams['ytick.minor.size'] = 5


mpl.rcParams['xtick.major.pad']='6'
#mpl.rcParams['legend.scatterpoints'] = 1

plt.axes().set_aspect(1)

# TODO: Set default edgecolor to none
