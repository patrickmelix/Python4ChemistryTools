import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, MultipleLocator, AutoLocator

matplotlib.rcParams['legend.numpoints'] = 1
matplotlib.rcParams['legend.fontsize'] = 18
matplotlib.rcParams['errorbar.capsize'] = 6
matplotlib.rcParams['axes.grid'] = True
matplotlib.rcParams['font.size'] = 18
matplotlib.rcParams['figure.figsize'] = (9.75, 5.85) #(10, 6)
matplotlib.rcParams['figure.dpi'] = 100
matplotlib.rcParams['savefig.dpi'] = 600
matplotlib.rcParams['legend.edgecolor'] = 'black'
#matplotlib.ticker.MaxNLocator.default_params['integer'] = True

cm = plt.get_cmap('gnuplot')
colors = ['blue', 'green', 'brown','red','black']
style = ['v','*','x','o','^','P','h','d','D','p','+']
msbig = 9
msmed = 8

eUnit = r"[$\mathrm{kJ}\,{\mathrm{mol}_{\mathrm{f.u.}}}^{-1}$]"
eUnitPresentation = r"[$\mathrm{kJ}\,\mathrm{mol}^{-1}$]"

kcal2kj = 4.184
