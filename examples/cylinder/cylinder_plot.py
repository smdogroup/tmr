import tikzplots as tkz
import argparse
import numpy as np

# Create an argument parser to read in arguments from the commnad line
p = argparse.ArgumentParser()
p.add_argument('--steps', type=int, default=5)
args = p.parse_args()

# Retrieve the number of steps
steps = args.steps

# Set the colors to use for each set of bars
colors = ['BrickRed', 'NavyBlue', 'black', 'ForestGreen',
          'Violet', 'Magenta' ]

data = []
for k in range(steps):
    data.append(np.loadtxt('results/cylinder_data%d.txt'%(k))[12:56,:])

# Delta value
delta = 20.0

# Set the positions of the tick locations
yticks = [0, 5, 15, 20]

# Set the values for the ticks
yticks = np.linspace(0, delta*len(data), 4*len(data)+1)
ytick_labels = []
for i in range(len(data)):
    ytick_labels.extend([0, 5, 10, 15])
ytick_labels.append(20)

# Set the value of the x-ticks
xvalues = len(data[0])
xticks = range(0, xvalues+1, 4)
xtick_labels = np.linspace(np.log10(data[0][0,0]), 
                           np.log10(data[0][-1,1]), len(xticks))

# Show the max/min value
xmin = 0
xmax = xvalues
ymin = 0
ymax = delta*len(data)

# The overall dimensions
xdim = 3.0
xscale = xdim/(2*xvalues)

ydim = 1.5
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

for y in np.linspace(0, delta*(len(data)-1), len(data)):
    s += tkz.get_2d_plot([xmin, xmax], [y, y],
                          xscale=xscale, yscale=yscale,
                          color='gray', line_dim='thick',
                          xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)

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

for k, d in enumerate(data):
    bars = []
    ymin = delta*(len(data)-k-1)
    for i in range(d.shape[0]):
        bars.append([d[i,3] + ymin])

    s += tkz.get_bar_chart(bars, color_list=[colors[k % len(colors)]], 
                            xscale=xscale, yscale=yscale, 
                            ymin=ymin, ymax=ymax)

s += tkz.get_end_tikz()

fp = open('cylinder_plot.tex', 'w')
fp.write(s)
fp.close()

