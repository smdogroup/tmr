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
p.add_argument('--labels', nargs='+', type=str, help='List of labels')
p.add_argument('--outfile', type=str, default='output.tex')
p.add_argument('--corrected', default=False, action='store_true')
p.add_argument('--error_offset', type=float, default=1e-2, help='Error order offset')
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
fval_corr_error_index = header.index('fval_corr_error')

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
    if args.corrected:
        ymin = min(ymin, np.min(d[:, fval_corr_error_index]))
    ymax = max(ymax, np.max(d[:, fval_error_index]))
    if args.corrected:
        ymax = max(ymax, np.max(d[:, fval_corr_error_index]))

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

s += r'\tikzstyle{dashed}= [dash pattern=on 6pt off 2pt]'
s += tikzcolors

symbols = ['circle', 'square', 'triangle', 'delta', 'diamond']

xerr = np.array([1.5e3, 4.5e4])
yerr = np.array([args.error_offset, 1.0])
order_list = [1, 2, 3, 4, 5, 6, 7]

for order in order_list:
    # Set the value of the order
    yerr[1] = yerr[0]*(xerr[0]/xerr[1])**(0.5*order)
    xvals = np.log10(xerr)
    yvals = np.log10(yerr)
    s += tkz.get_2d_plot(xvals, yvals,
                         line_dim='thin',
                         color='Gray!50', symbol=None,
                         xscale=xscale, yscale=yscale, 
                         xmin=xmin, xmax=xmax,
                         ymin=ymin, ymax=ymax)

    if ((xvals[1] >= xmin and xvals[1] <= xmax) and
        (yvals[1] >= ymin and yvals[1] <= ymax)):
        s += r'\draw[color=Gray!80, font=\scriptsize] (%e, %e) node[right] {%d};'%(
            xscale*xvals[1], yscale*yvals[1], order)

for k, d in enumerate(data):
    xvals = np.log10(d[:, nnodes_index])
    yvals = np.log10(d[:, fval_error_index])    
    s += tkz.get_2d_plot(xvals, yvals,
                         line_dim='very thick',
                         color=colors[k % len(colors)],
                         symbol=symbols[k % len(symbols)],
                         symbol_size=0.035,
                         xscale=xscale, yscale=yscale, 
                         xmin=xmin, xmax=xmax,
                         ymin=ymin, ymax=ymax)

    if args.corrected:
        yvals = np.log10(d[:, fval_corr_error_index])
        s += tkz.get_2d_plot(xvals, yvals,
                             line_dim='very thick, opacity=0.5, dashed',
                             color=colors[k % len(colors)],
                             symbol=symbols[k % len(symbols)],
                             symbol_size=0.035,
                             xscale=xscale, yscale=yscale, 
                             xmin=xmin, xmax=xmax,
                             ymin=ymin, ymax=ymax)

# Set the labels (lower-left corner)
if args.labels is not None:
    for k, label in enumerate(args.labels):
        x = xmin + 0.05*(xmax - xmin)
        y = ymin + 0.05*(ymax - ymin)*(len(args.labels)-k)
        length = 0.035*(xmax - xmin)
        s += tkz.get_legend_entry(x, y, length, label=label,
                                  font_size='small',
                                  line_dim='very thick',
                                  color=colors[k % len(colors)],
                                  symbol=symbols[k % len(symbols)],
                                  symbol_size=0.035,
                                  xscale=xscale, yscale=yscale)

# Plot the axes
s += tkz.get_2d_axes(xmin, xmax, ymin, ymax,
                     xscale=xscale, yscale=yscale,
                     xticks=xticks, yticks=yticks,
                     xtick_labels=xtick_labels,
                     ytick_labels=ytick_labels,
                     tick_font='normalsize',
                     tick_frac=0.01,
                     xlabel_offset=0.085,
                     label_font='Large',
                     xlabel='Number of nodes',
                     ylabel_offset=0.175,
                     ylabel='Functional error')
                     

s += tkz.get_end_tikz()

fp = open(args.outfile, 'w')
fp.write(s)
fp.close()

