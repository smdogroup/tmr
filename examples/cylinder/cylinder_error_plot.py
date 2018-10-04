from __future__ import print_function
import tikzplots as tkz
import argparse
import numpy as np
import re

def parse_data_file(fname):
    with open(fname, 'r') as fp:
        lines = fp.readlines()

        # Read in the first line, and find the comma-separated values
        # in the header
        hline = lines[0]
        for index, h in enumerate(hline):
            if h == '=':
                hstr = hline[index+1:].split(',')

        # Strip away any white space
        header = []
        for h in hstr:
            header.append(h.strip())     

        data = []
        for line in lines[1:]:
            dline = []
            for entry in line.split():
                dline.append(float(entry))
            data.append(dline)

        return header, np.array(data)


# Create an argument parser to read in arguments from the commnad line
p = argparse.ArgumentParser()
p.add_argument('--files', nargs='+', type=str, help='List of files')
p.add_argument('--outfile', type=str, default='output.tex')
args = p.parse_args()

# Set the colors to use for each set of bars
colors = []
for i in range(10):
    colors.append('tableau%d'%(i))

tikzcolors = '''
\definecolor{tableau0}{RGB}{31,119,180}
\definecolor{tableau1}{RGB}{255,158,74}
\definecolor{tableau2}{RGB}{103,191,92}
\definecolor{tableau3}{RGB}{237,102,93}
\definecolor{tableau4}{RGB}{148,103,189}
\definecolor{tableau5}{RGB}{168,120,110}
\definecolor{tableau6}{RGB}{237,151,202}
\definecolor{tableau7}{RGB}{162,162,162}
\definecolor{tableau8}{RGB}{205,204,93}
\definecolor{tableau9}{RGB}{109,204,218}
'''

data = []
for fname in args.files:
    header, dat = parse_data_file(fname)
    data.append(dat)

# Plot the error on the y-axis
nnodes_index = header.index('nnodes')
fval_error_index = header.index('fval_error')

# Find the max value of y
xmin = 1e20
xmax = 0

ymin = 1e20
ymax = 0

# Look through all the data
for d in data:
    xmin = min(xmin, np.min(d[:, nnodes_index]))
    xmax = max(xmax, np.max(d[:, nnodes_index]))

    ymin = min(ymin, np.min(d[:, fval_error_index]))
    ymax = max(ymax, np.max(d[:, fval_error_index]))

# Round to the nearest multiple of 10
xmin = int(np.floor(np.log10(xmin)))
xmax = int(np.ceil(np.log10(xmax)))

ymin = int(np.floor(np.log10(ymin)))
ymax = int(np.ceil(np.log10(ymax)))

# Create a range
xticks = np.linspace(xmin, xmax, xmax - xmin + 1)
xtick_labels = []
for exp in range(xmin, xmax + 1, 1):
    xtick_labels.append('$10^{%d}$'%(exp))

yticks = np.linspace(ymin, ymax, ymax - ymin + 1)
ytick_labels = []
for exp in range(ymin, ymax + 1, 1):
    ytick_labels.append('$10^{%d}$'%(exp))

# The overall dimensions
xdim = 2.0
xscale = xdim/(xmax - xmin)

ydim = 1.75
yscale = ydim/(ymax - ymin)

# Get the header info
s = tkz.get_header()
s += tkz.get_begin_tikz(xdim=1.5, ydim=1.5, xunit='in', yunit='in')

s += tikzcolors

symbols = ['circle', 'square', 'triangle', 'delta', 'diamond']

for k, d in enumerate(data):
    xvals = np.log10(d[:, nnodes_index])
    yvals = np.log10(d[:, fval_error_index])
    
    s += tkz.get_2d_plot(xvals, yvals,
                         color=colors[k % 4],
                         symbol=symbols[k % 4 ],
                         symbol_size=0.03,
                         xscale=xscale, yscale=yscale, 
                         xmin=xmin, xmax=xmax,
                         ymin=ymin, ymax=ymax)

# Plot the axes
s += tkz.get_2d_axes(xmin, xmax, ymin, ymax,
                     xscale=xscale, yscale=yscale,
                     xticks=xticks, yticks=yticks,
                     xtick_labels=xtick_labels,
                     ytick_labels=ytick_labels,
                     tick_font='small',
                     tick_frac=0.01,
                     xlabel_offset=0.075,
                     label_font='large',
                     xlabel='Number of nodes',
                     ylabel_offset=0.165,
                     ylabel='Functional error')
                     

s += tkz.get_end_tikz()

fp = open(args.outfile, 'w')
fp.write(s)
fp.close()

