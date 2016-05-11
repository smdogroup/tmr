import sys
import numpy as np

# import matplotlib as mpl
# mpl.use("pgf")

# pgf_with_custom_preamble = {
# #    "font.family": "sansserif", # use serif/main font for text elements
# #    "text.usetex": True,    # use inline math for ticks
# #    "pgf.rcfonts": False,   # don't setup fonts from rc parameters
#     "pgf.preamble": [
#          "\\usepackage{sfmath}"]
# }
# mpl.rcParams.update(pgf_with_custom_preamble)

tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),    
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),    
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),    
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),    
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)] 


# Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.    
for i in range(len(tableau20)):    
    r, g, b = tableau20[i]    
    tableau20[i] = (r / 255., g / 255., b / 255.)    

import matplotlib.pylab as plt


# Get the name of the input file
fname = sys.argv[1]
fp = open(fname, 'r')
if fp:
    fp.readline()
    values = np.loadtxt(fp)

    # plt.figure()
    fig, ax1 = plt.subplots(figsize=(8,5.5))

    nnodes = values[:,2]
    fvals = values[:,3]
    fcorr = values[:,4]
    err = values[:,5]

    ax1.semilogx(1.0/np.sqrt(nnodes), fvals, '-s', color=tableau20[0], 
                 linewidth=2, markersize=7, label='function')
    ax1.semilogx(1.0/np.sqrt(nnodes), fcorr, '-o', color=tableau20[2], 
                 linewidth=2, markersize=7, label='corrected function')

    # Twin the axis
    ax2 = ax1.twinx()

    ax2.loglog([], [], '-s', color=tableau20[0], 
               linewidth=2, markersize=7, label='functional estimate')
    ax2.loglog([], [], '-o', color=tableau20[2], 
               linewidth=2, markersize=7, label='corrected functional')    
    ax2.loglog(1.0/np.sqrt(nnodes), err, '->', color=tableau20[4],  
               linewidth=2, markersize=7, label='error estimate')

    # ax1.set_ylim([1.0, 1.1])
    # ax1.set_ylim([0.95, 1.05])
    ax2.set_xlim([0.003, 0.05])
    ax1.set_xlim([0.003, 0.05])

    ax1.set_xlabel('Average $\\sf \Delta x$', 
                   fontweight='bold', fontsize=18, y=1.05)
    ax1.set_ylabel('Relative function value', 
                   fontweight='bold', fontsize=18)
    ax2.set_ylabel('Error estimate', 
                   fontweight='bold', fontsize=18)
    
    ax1.tick_params(axis='both', which='major', labelsize=15)
    ax2.tick_params(axis='both', which='major', labelsize=15)
    
    plt.gcf().subplots_adjust(bottom=0.13)
    ax2.legend(loc=4)
    
    split = fname.split('.')
    outname = split[0]+'.pdf'

    plt.savefig(outname)
