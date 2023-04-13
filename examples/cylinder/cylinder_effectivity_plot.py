from __future__ import print_function
import tikzplots as tkz
import argparse
import numpy as np
import re


def parse_data_file(fname):
    with open(fname, "r") as fp:
        lines = fp.readlines()

        # Read in the first line, and find the comma-separated values
        # in the header
        hline = lines[0]
        for index, h in enumerate(hline):
            if h == "=":
                hstr = hline[index + 1 :].split(",")

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
p.add_argument("--files", nargs="+", type=str, help="List of files")
p.add_argument("--labels", nargs="+", type=str, help="List of labels")
p.add_argument("--outfile", type=str, default="output.tex")
p.add_argument("--plot", type=str, default="effectivity")
args = p.parse_args()

# Set the colors to use for each set of bars
colors = []
for i in range(10):
    colors.append("tableau%d" % (i))

tikzcolors = """
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
"""

data = []
for fname in args.files:
    try:
        header, dat = parse_data_file(fname)
    except:
        print("fname = ", fname)
    data.append(dat)

# Plot the error on the y-axis
nnodes_index = header.index("nnodes")
fval_eff_index = header.index("fval_effectivity")
indc_eff_index = header.index("indicator_effectivity")

# Find the max value of y
xmin = 1e20
xmax = 0

ymin = 0
ymax = 0

# Look through all the data
for d in data:
    xmin = min(xmin, np.min(d[:, nnodes_index]))
    xmax = max(xmax, np.max(d[:, nnodes_index]))

    if args.plot == "effectivity":
        ymax = max(ymax, np.max(d[:, fval_eff_index]))
        ymax = min(ymax, 100)
    else:
        ymax = max(ymax, np.max(d[:, indc_eff_index]))
        ymax = min(ymax, 500)

# Round to the nearest multiple of 10
xmin = int(np.floor(np.log10(xmin)))
xmax = int(np.ceil(np.log10(xmax)))

# Create a range
xticks = np.linspace(xmin, xmax, xmax - xmin + 1)
xtick_labels = []
for exp in range(xmin, xmax + 1, 1):
    xtick_labels.append("$10^{%d}$" % (exp))

# Set the positions of the tick locations
if ymax < 2.0:
    ymax_int = int(np.ceil(4.0 * ymax))
    ymax = ymax_int / 4.0
    yticks = np.linspace(0, ymax, ymax_int + 1)
    ytick_labels = yticks
elif ymax < 10:
    ymax = int(np.ceil(ymax))
    yticks = np.linspace(0, ymax, ymax + 1)
    ytick_labels = range(ymax + 1)
elif ymax < 20:
    ymax = 2 * int(np.ceil(ymax / 2.0))
    yticks = np.linspace(0, ymax, ymax + 1)
    ytick_labels = range(0, ymax + 1, 2)
    yticks = np.linspace(0, ymax, ymax / 2 + 1)
else:
    ymax = 5 * int(np.ceil(ymax / 5.0))
    yticks = np.linspace(0, ymax, ymax + 1)
    ytick_labels = range(0, ymax + 1, 5)
    yticks = np.linspace(0, ymax, ymax / 5 + 1)

# The overall dimensions
xdim = 2.0
xscale = xdim / (xmax - xmin)

ydim = 1.75
yscale = ydim / (ymax - ymin)

# Get the header info
s = tkz.get_header()
s += tkz.get_begin_tikz(xdim=1.5, ydim=1.5, xunit="in", yunit="in")

s += tikzcolors

symbols = ["circle", "square", "triangle", "delta", "diamond"]

for k, d in enumerate(data):
    xvals = np.log10(d[:, nnodes_index])
    if args.plot == "effectivity":
        yvals = d[:, fval_eff_index]
    else:
        yvals = d[:, indc_eff_index]

    s += tkz.get_2d_plot(
        xvals,
        yvals,
        line_dim="very thick",
        color=colors[k % 10],
        symbol=symbols[k % 4],
        symbol_size=0.035,
        xscale=xscale,
        yscale=yscale,
        xmin=xmin,
        xmax=xmax,
        ymin=ymin,
        ymax=ymax,
    )

# Set the labels (lower-right corner)
if args.labels is not None:
    for k, label in enumerate(args.labels):
        x = xmin + 0.75 * (xmax - xmin)
        y = ymin + 0.05 * (ymax - ymin) * (len(args.labels) - k)
        length = 0.035 * (xmax - xmin)
        s += tkz.get_legend_entry(
            x,
            y,
            length,
            label=label,
            font_size="small",
            line_dim="very thick",
            color=colors[k % 10],
            symbol=symbols[k % 4],
            symbol_size=0.035,
            xscale=xscale,
            yscale=yscale,
        )

if args.plot == "effectivity":
    title = "Effectivity"
else:
    title = "Indicator effectivity"

# Plot the axes
s += tkz.get_2d_axes(
    xmin,
    xmax,
    ymin,
    ymax,
    xscale=xscale,
    yscale=yscale,
    xticks=xticks,
    yticks=yticks,
    xtick_labels=xtick_labels,
    ytick_labels=ytick_labels,
    tick_font="normalsize",
    tick_frac=0.01,
    xlabel_offset=0.085,
    label_font="Large",
    xlabel="Number of nodes",
    ylabel_offset=0.175,
    ylabel=title,
)

s += tkz.get_end_tikz()

fp = open(args.outfile, "w")
fp.write(s)
fp.close()
