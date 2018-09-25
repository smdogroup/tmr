from __future__ import print_function
import tikzplots as tkz
import argparse
import numpy as np

# Create an argument parser to read in arguments from the commnad line
p = argparse.ArgumentParser()
p.add_argument('--steps', type=int, default=5)
p.add_argument('--case', type=str, default='disk')
args = p.parse_args()

# Retrieve the number of steps
steps = args.steps
case = args.case

# Set the colors to use for each set of bars
colors = ['lightgray', 'ForestGreen']

data = []
for k in [0, steps-1]:
    data.append(np.loadtxt('results/%s_data%d.txt'%(case, k)))

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
xtick_labels = range(x0, x1-1, -1)

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
elif ymax < 30:
    ymax = 5*int(np.ceil(ymax/5.0))
    yticks = np.linspace(0, ymax, ymax+1)
    ytick_labels = range(0, ymax+1, 5)
    yticks = np.linspace(0, ymax, ymax/5 + 1)

print('xticks = ', xticks)
print('xtick_labels = ', xtick_labels)

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


# Plot the axes
s += tkz.get_2d_axes(xmin, xmax, ymin, ymax,
                     xscale=xscale, yscale=yscale,
                     xticks=xticks, yticks=yticks,
                     xtick_labels=xtick_labels,
                     ytick_labels=ytick_labels,
                     tick_font='small',
                     tick_frac=0.0125,
                     xlabel_offset=0.1,
                     xlabel='$\\log_{10}(\\text{Error})$', 
                     ylabel_offset=0.065,
                     ylabel='Percentage')

s += tkz.get_end_tikz()

fp = open('poisson_%s_plot.tex'%(case), 'w')
fp.write(s)
fp.close()

