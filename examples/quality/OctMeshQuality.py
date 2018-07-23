import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import PercentFormatter
import matplotlib.colors as colors
import matplotlib.cm as cmx


# Configure the plotting environment
plt.rcParams['xtick.direction'] = 'out'
plt.rcParams['ytick.direction'] = 'out'

# Optionally set font to Computer Modern to avoid common missing font errors
params = {
    'axes.labelsize': 20,
    'legend.fontsize': 14,
    'xtick.labelsize': 20,
    'ytick.labelsize': 20,
    'text.usetex': True}
plt.rcParams.update(params)

# Latex math
plt.rcParams['text.latex.preamble'] = [r'\usepackage{sfmath}']
plt.rcParams['font.family'] = 'sans-serif'
# plt.rcParams['font.sans-serif'] = 'courier'
plt.rcParams['font.size'] = 18
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['lines.linewidth'] = 4
plt.rcParams['lines.color'] = 'r'

# Make sure everything is within the frame
plt.rcParams.update({'figure.autolayout': True})

# These are the "Tableau 20" colors as RGB.    
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),    
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),    
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),    
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),    
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]    
  
# Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.    
for i in range(len(tableau20)):    
    r, g, b = tableau20[i]    
    tableau20[i] = (r / 255., g / 255., b / 255.)

# Color scheme
hist_color = tableau20[0]

# Compute the relative size of the element in the forest
def computeShape(forest):
    octs = forest.getOctants()
    Xpts = forest.getPoints()
    conn = forest.getMeshConn()
    fshape = np.zeros(len(octs))
    jacobian_list={0:[1,2,4], 1:[3,0,5], 2:[0,3,6],
                   3:[2,1,7], 4:[6,5,0], 5:[4,7,1],
                   6:[7,4,2], 7:[5,6,3]}
    for i in range(len(octs)):
        npts = 8
        nodes = conn[npts*i:npts*(i+1)]
        pts = np.zeros((npts,3))
        pts[:, :] = Xpts[nodes[:], :]
        # Initialize the numerator
        num = 0.0
        # Compute the individual Ak matrices
        for j in range(npts):
            # Get the reference vertices of volume
            v = jacobian_list[j]
            Ak = np.array([pts[v[0],:]-pts[j,:], pts[v[1],:]-pts[j,:], \
                           pts[v[2],:]-pts[j,:] ])
            # Get determinant of local matrix alpha
            #alpha[j] = np.linalg.det(Ak)
            #sigma[j] = np.sum(np.dot(Ak.T, Ak).diagonal())
            num += np.sum(np.dot(Ak.T, Ak).diagonal())/np.linalg.det(Ak)**(2./3.)
        # Compute the shape metric 24/num
        fshape[i] = 24./num
        
    return fshape

# Compute the aspect ratio of each element in the forest
def computeAR(forest):
    octs = forest.getOctants()
    Xpts = forest.getPoints()
    conn = forest.getMeshConn()
    ar = np.zeros(len(octs))
    # List of vertices which create an edge on the octant
    edge_list = [[0, 1], [1, 3], [2, 3], [0, 2], [4, 5], [5, 7],
                 [6, 7], [4, 6], [1, 5], [3, 7], [2, 6], [0, 4]]
    for i in range(len(octs)):
        # Get the points
        npts = 8
        nodes = conn[npts*i:npts*(i+1)]
        pts = np.zeros((npts,3))
        pts[:, :] = Xpts[nodes[:], :]
        edge_lengths = np.zeros(12)
        # Compute the length of each edge on the octant
        for j, edge_pts in enumerate(edge_list):
            v1 = pts[edge_pts[0]]
            v2 = pts[edge_pts[1]]
            edge_lengths[j] = np.linalg.norm(v1-v2)
        # Compute the aspect ratio of the octant
        ar[i] = np.amax(edge_lengths)/np.amin(edge_lengths)

    return ar

# Compute the minimum angle of each element in the forest
def computeMinAngle(forest):
    octs = forest.getOctants()
    Xpts = forest.getPoints()
    conn = forest.getMeshConn()
    min_angs = np.zeros(len(octs))
    # List of verticies corresponding to interior angles
    angle_list = {0:[1, 2, 4], 1:[0, 5, 3],
                  2:[0, 6, 3], 3:[1, 2, 7],
                  4:[0, 5, 6], 5:[1, 4, 7],
                  6:[2, 4, 7], 7:[3, 5, 6]}
    for i in range(len(octs)):
        # Get the points
        npts = 8
        nodes = conn[npts*i:npts*(i+1)]
        pts = np.zeros((npts,3))
        pts[:, :] = Xpts[nodes[:], :]
        min_angle = 90.0
        for j in range(npts):
            node_neighbors = angle_list[j]
            for k in node_neighbors:
                for l in node_neighbors:
                    if k <> l:
                        # Compute 2 vectors from 3 unique points
                        vec1 = pts[k, :] - pts[j, :]
                        vec2 = pts[l, :] - pts[j, :]
                        # Make the vectors unit vectors
                        if (np.linalg.norm(vec1) <> np.float64(0.0)) & (np.linalg.norm(vec2) <> np.float64(0.0)):
                            vec1 /= np.linalg.norm(vec1)
                            vec2 /= np.linalg.norm(vec2)
                            # Compute the angle between the vectors as long as they
                            # are different
                            if np.array_equal(vec1, vec2) is False:
                                angle = (180.0/np.pi)*np.arccos(np.dot(vec1, vec2))
                                # Update min_angle if needed
                                if angle < min_angle:
                                    min_angle = angle
        min_angs[i] = min_angle

    return min_angs

# Write out mesh quality metrics to vtk file
def writeQualityToVtk(forest, ar, min_ang, fshape, fname='quality.vtk'):
    octs = forest.getOctants()
    Xpts = forest.getPoints()
    conn = forest.getMeshConn()
    npts = len(Xpts)
    nhex = len(octs)

    f = open(fname, 'w')
    f.write('# vtk DataFile Version 3.0\n')
    f.write('vtk output\nASCII\n')
    f.write('DATASET UNSTRUCTURED_GRID\n')

    # write out the points
    f.write('POINTS %d float\n'%(npts))
    for pt in Xpts:
        f.write('%e %e %e\n'%(pt[0], pt[1], pt[2]))
    
    # write out the mesh connectivity
    f.write('\nCELLS %d %d\n'%(nhex, nhex*9))
    for i in range(len(octs)):
        nodes = conn[8*i:8*(i+1)]
        f.write('8 %d %d %d %d %d %d %d %d\n'%(nodes[0], nodes[1], nodes[3], nodes[2], nodes[4], nodes[5], nodes[7], nodes[6]))

    # all hex
    f.write('\nCELL_TYPES %d\n'%(nhex))
    for i in range(nhex):
        f.write('%d\n'%(12))

    # write AR values
    f.write('\nCELL_DATA %d\n'%(nhex))
    f.write('SCALARS aspect_ratio float 1\n')
    f.write('LOOKUP_TABLE default\n')
    for elem_ar in ar:
        f.write('%e\n'%(elem_ar))

    # write minimum angle values
    f.write('SCALARS min_angle float 1\n')
    f.write('LOOKUP_TABLE default\n')
    for elem_ang in min_ang:
        f.write('%e\n'%(elem_ang))

    # Write the shape metric
    # write minimum angle values
    f.write('SCALARS shape_metric float 1\n')
    f.write('LOOKUP_TABLE default\n')
    for elem_shape in fshape:
        f.write('%e\n'%(elem_shape))
        
    f.close()

    return

def plotARHist(ar, xmax=None, fname='ar_hist.pdf'):
    hist_max = int(np.ceil(np.amax(ar)))
    nbins = 40*(hist_max-1)

    # Create the figure and set parameters
    fig, ax = plt.subplots(1, 1, figsize=(12, 9))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

    # Plot the histogram
    n, bins, patches = plt.hist(ar, bins=nbins, range=(1.0, hist_max),
                                density=True, stacked=False)
    widths = bins[:-1] - bins[1:]
    
    # Set the colors
    norm = colors.Normalize(ar.min(), ar.max())
    for b, p in zip(bins, patches):
        color = plt.cm.coolwarm(norm(b))
        p.set_facecolor(color)

    # Configure the axes
    if xmax:
        plt.xlim((1.0, xmax))
    else:
        plt.xlim((1.0, hist_max))
    plt.xlabel('Aspect Ratio')
    ax.yaxis.set_major_formatter(PercentFormatter(xmax=-1.0/widths[0], decimals=0))
    ax.tick_params(which='both', direction='out', top=False, right=False)
    plt.savefig(fname, bbox_inches='tight', pad_inches=0.05)

    return

def plotMinAngHist(min_ang, xmin=None, fname='angle_hist.pdf'):
    hist_max = 90.0
    nbins = 90

    # Create the figure and set parameters
    fig, ax = plt.subplots(1, 1, figsize=(12, 9))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

    # Plot the histogram
    n, bins, patches = plt.hist(min_ang, bins=nbins, range=(0.0, hist_max),
                                align='mid', density=True)
    # Set the colors
    norm = colors.Normalize(min_ang.min(), min_ang.max())
    for b, p in zip(bins, patches):
        color = plt.cm.coolwarm_r(norm(b))
        p.set_facecolor(color)

    # Configure the axes
    x_ticks = [0, 15, 30, 45, 60, 75, 90]
    plt.xticks(x_ticks, x_ticks)
    if xmin:
        plt.xlim((xmin, 90.0))
    else:
        plt.xlim((0.0, 90.0))
    plt.xlabel('Minimum Angle (deg.)')
    ax.yaxis.set_major_formatter(PercentFormatter(xmax=1, decimals=0))
    ax.tick_params(which='both', direction='out', top=False, right=False)
    plt.tight_layout()
    plt.savefig(fname, bbox_inches='tight', pad_inches=0.05)

    return

def plotShapeHist(fshape, xmin=None, fname='shape_hist.pdf'):
    hist_min = np.amin(fshape)
    nbins = 50
    #nbins = 40*(hist_min-1)

    # Create the figure and set parameters
    fig, ax = plt.subplots(1, 1, figsize=(12, 9))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

    # Plot the histogram
    n, bins, patches = plt.hist(fshape, bins=nbins, range=(xmin, 1.0), \
                                density=True)

    widths = bins[:-1] - bins[1:]
    
    # Set the colors
    norm = colors.Normalize(fshape.min(), fshape.max())
    for b, p in zip(bins, patches):
        color = plt.cm.coolwarm(norm(b))
        p.set_facecolor(color)
            
    plt.xlabel('Hexahedron Shape Metric')
    ax.yaxis.set_major_formatter(PercentFormatter(xmax=-1/widths[0], decimals=0))
    ax.tick_params(which='both', direction='out', top=False, right=False)
    plt.tight_layout()
    plt.savefig(fname, bbox_inches='tight', pad_inches=0.05,dpi=1000)

    return
