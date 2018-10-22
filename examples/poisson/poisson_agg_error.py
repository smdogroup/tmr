'''
Plot the aggregation error as a function of rho for the KS and p-norm
functionals
'''
from __future__ import print_function
import tikzplots as tkz
import numpy as np

def integrate(integrand):
    '''Integrate over equally spaced data'''
    sigma = [17.0/48.0, 59.0/48.0, 43.0/48, 49/48.9]
    r = len(sigma)

    integral = 0.0
    for i, s in enumerate(sigma):
        integral += s*integrand[i]
        integral += s*integrand[-1-i]

    for i in range(r, len(integrand)-r):
        integral += integrand[i]

    return integral

def get_disk_aggregate(functional, rho, R, n=1000):
    '''Evaluate the KS functional on a disk'''
    scale = 1.0/R**2
    r = np.linspace(0.0, R, n)

    a0 = 3920.0/363.0

    # This is the solution for constant order
    # phi = (0.1875*R**2 - 0.25*r**2 + (0.0625/R**2)*r**4)

    # Compute the solution scaled to r/R
    x = np.linspace(0.0, 1.0, n)
    phi = a0*R**2*(- (1.0/196)*x**14
                   + (1.0/24)*x**12
                   - (3.0/20)*x**10
                   + (5.0/16)*x**8
                   - (5.0/12)*x**6
                   + (3.0/8)*x**4
                   - (1.0/4)*x**2
                   + (363.0/3920))

    if functional == 'ks':
        ksmax = np.max(phi)
        integrand = r*np.exp(rho*scale*(phi - ksmax))
        kssum = 2*np.pi*(R/(n-1))*integrate(integrand)
        return scale*ksmax, scale*ksmax + np.log(kssum)/rho
    else:
        maxphi = scale*np.max(phi)
        integrand = r*np.power(np.fabs(scale*phi/maxphi), rho)
        psum = 2*np.pi*(R/(n-1))*integrate(integrand)
        return maxphi, maxphi*np.power(psum, 1.0/rho)

n = 100000
functional = 'ks'

xmin = -6
xmax = -3

m = 91
eps_vals = np.logspace(xmin, xmax, m)
ks_agg_err = np.zeros(m)
pnorm_agg_err = np.zeros(m)

with open('results/disk_aggregation_errors.dat', 'w') as fp:
    fp.write('Variables = rho, 1/rho, ks, pnorm\n')
    R = 100.0
    for k, eps in enumerate(eps_vals):
        rho = 1.0/eps
        max_value, ks = get_disk_aggregate('ks', rho, R, n)
        max_value, pnorm = get_disk_aggregate('pnorm', rho, R, n)
        ks_agg_err[k] = np.fabs(ks - max_value)
        pnorm_agg_err[k] = np.fabs(pnorm - max_value)
        fp.write('%25.16e %25.16e %25.16e %25.16e\n'%(
            rho, eps, ks_agg_err[k], pnorm_agg_err[k]))

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

ymax = max(np.max(ks_agg_err), np.max(pnorm_agg_err))
ymin = min(np.min(ks_agg_err), np.min(pnorm_agg_err))

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

xvals = np.log10(eps_vals)
yvals = np.log10(ks_agg_err)
s += tkz.get_2d_plot(xvals, yvals,
                     color=colors[0],
                     symbol=None,
                     xscale=xscale, yscale=yscale, 
                     xmin=xmin, xmax=xmax,
                     ymin=ymin, ymax=ymax)
s += tkz.get_2d_plot(xvals[::3], yvals[::3], line_dim=None,
                     color=colors[0],
                     symbol=symbols[0],
                     symbol_size=0.03,
                     xscale=xscale, yscale=yscale, 
                     xmin=xmin, xmax=xmax,
                     ymin=ymin, ymax=ymax)

xvals = np.log10(eps_vals)
yvals = np.log10(pnorm_agg_err)
s += tkz.get_2d_plot(xvals, yvals,
                     color=colors[1],
                     symbol=None,
                     xscale=xscale, yscale=yscale, 
                     xmin=xmin, xmax=xmax,
                     ymin=ymin, ymax=ymax)
 
s += tkz.get_2d_plot(xvals[::3], yvals[::3], line_dim=None,
                     color=colors[1],
                     symbol=symbols[1],
                     symbol_size=0.03,
                     xscale=xscale, yscale=yscale, 
                     xmin=xmin, xmax=xmax,
                     ymin=ymin, ymax=ymax)


# Set the labels (lower-left corner)
for k, label in enumerate(['KS', '$p$-norm']):
    x = xmin + 0.05*(xmax - xmin)
    y = ymin + (ymax - ymin)*(0.95 - 0.05*k)
    length = 0.035*(xmax - xmin)
    s += tkz.get_legend_entry(x, y, length, label=label,
                              font_size='small',
                              color=colors[k % len(colors)],
                              symbol=symbols[k % len(symbols)],
                              symbol_size=0.03,
                              xscale=xscale, yscale=yscale)

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
                     xlabel='1/$\\rho$',
                     ylabel_offset=0.165,
                     ylabel='Aggregation error')
                     

s += tkz.get_end_tikz()

fp = open('results/disk_exact_aggregation_error.tex', 'w')
fp.write(s)
fp.close()
