from __future__ import print_function
import tikzplots as tkz
import argparse
import numpy as np

# Create an argument parser to read in arguments from the commnad line
p = argparse.ArgumentParser()
p.add_argument('--files', nargs='+', type=str, help='List of files')
p.add_argument('--labels', nargs='+', type=str, help='List of labels')
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
    data.append(np.loadtxt(fname))

# Find the max value of y
ymax = 0

# Look for all the data
bins_per_decade = 10
idx_min = data[0].shape[0]
idx_max = 0
for d in data:
    # Loop over all the bins
    for i in range(0, d.shape[0], bins_per_decade):
        flag = False
        for k in range(bins_per_decade):
            if d[i+k,3] > 0.01:
                flag = True
            if d[i+k,3] > ymax:
                ymax = d[i+k,3]

        # Set the new value for idx_min
        if flag and i < idx_min:
            idx_min = i
        if flag and i > idx_max:
            idx_max = i

idx_max = min(idx_max + bins_per_decade-1, data[0].shape[0]-1)

# Find the largest error value
x0 = int(np.ceil(np.log10(data[0][idx_min,0])))

# Find the smallest error value
x1 = int(np.floor(np.log10(data[0][idx_max,1])))

# Set the number of x values
xvalues = idx_max - idx_min + 1

# Create a range
xticks = xvalues*(np.linspace(x0, x1, x0 - x1 + 1) - x0)/(x1 - x0)

xtick_labels = []
for exp in range(x0, x1-1, -1):
    xtick_labels.append('$10^{%d}$'%(exp))

# Set the positions of the tick locations
if ymax < 10:
    ymax = int(np.ceil(ymax))
    yticks = np.linspace(0, ymax, ymax+1)
    ytick_labels = range(ymax+1)
elif ymax < 20:
    ymax = 2*int(np.ceil(ymax/2.0))
    yticks = np.linspace(0, ymax, ymax+1)
    ytick_labels = range(0, ymax+1, 2)
    yticks = np.linspace(0, ymax, ymax/2 + 1)
else:
    ymax = 5*int(np.ceil(ymax/5.0))
    yticks = np.linspace(0, ymax, ymax+1)
    ytick_labels = range(0, ymax+1, 5)
    yticks = np.linspace(0, ymax, ymax/5 + 1)

# Show the max/min value
xmin = 0
xmax = xvalues

ymin = 0

# The overall dimensions
xdim = 3.0
xscale = xdim/(2.5*xvalues)

ydim = 1.25
yscale = ydim/ymax

# Get the header info
s = tkz.get_header()
s += tkz.get_begin_tikz(xdim=3.5, ydim=2.0, xunit='in', yunit='in')

s += tikzcolors

# Create the plot background
for y in yticks:
    s += tkz.get_2d_plot([xmin, xmax], [y, y],
                          xscale=xscale, yscale=yscale,
                          color='gray', line_dim='thin',
                          xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)

for k, d in enumerate(data):
    bars = []
    ymin = 0
    for i in range(idx_min, idx_max+1):
        bars.append([d[i,3] + ymin])

    s += tkz.get_bar_chart(bars, color_list=[colors[k % len(colors)]], 
                            xscale=xscale, yscale=yscale, 
                            ymin=ymin, ymax=ymax)


# Set the labels (lower-left corner)
if args.labels is not None:
    x = xmin + 0.75*(xmax - xmin)
    y = ymin + (ymax - ymin)*(0.95 - 0.025*(len(args.labels)-1))

    xlen = 0.2*(xmax - xmin)
    ylen = 0.05*len(args.labels)*(ymax - ymin)
    
    x1 = x - 0.1*xlen
    x2 = x + 0.9*xlen
    y1 = y - 0.5*ylen
    y2 = y + 0.5*ylen

    s += r'\draw[color=white, fill=white] (%f,%f) rectangle (%f,%f);'%(
        xscale*x1, yscale*y1, xscale*x2, yscale*y2)


    for k, label in enumerate(args.labels):
        y = ymin + (ymax - ymin)*(0.95 - 0.05*k)

        xlen = 0.025*(xmax - xmin)
        ylen = 0.02*(ymax - ymin)
        x1 = x - 0.5*xlen
        x2 = x + 0.5*xlen
        y1 = y - 0.5*ylen
        y2 = y + 0.5*ylen
        line_dim='thick'
        color = colors[k % len(colors)]
        font_size = 'small'

        s += r'\draw[%s, color=%s, fill=%s, fill opacity=0.3]'%(
            line_dim, color, color)
        s += ' (%f, %f) rectangle (%f, %f);'%(
            xscale*x1, yscale*y1, xscale*x2, yscale*y2)
        s += r'\draw[font=\%s] (%f, %f) node[right] {%s};'%(
            font_size, xscale*(x + xlen), yscale*y, label)

# Plot the axes
s += tkz.get_2d_axes(xmin, xmax, ymin, ymax,
                     xscale=xscale, yscale=yscale,
                     xticks=xticks, yticks=yticks,
                     xtick_labels=xtick_labels,
                     ytick_labels=ytick_labels,
                     tick_font='small',
                     tick_frac=0.0125,
                     xlabel_offset=0.1,
                     xlabel='Element error',
                     ylabel_offset=0.075,
                     ylabel='Percentage')

s += tkz.get_end_tikz()

fp = open(args.outfile, 'w')
fp.write(s)
fp.close()
