import sys, os
import numpy as np
from tacs import TACS

import matplotlib
matplotlib.use('Agg')

import matplotlib.pylab as plt
import matplotlib.tri as tri

for filename in sys.argv[1:]:
    outfile = None
    if not os.path.isfile(filename):
        print('File %s does not exist'%(filename))
        exit(0)
    else:
        outfile = os.path.splitext(filename)[0] + '.png'

    loader = TACS.FH5Loader()
    loader.loadData(filename)

    # Load the component types and data from the file
    comps, ltypes, ptr, conn = loader.getConnectivity()
    c_var_names, cdata = loader.getContinuousData()
    e_var_names, edata = loader.getElementData()

    # Get the number of elements
    num_elements = len(comps)

    # Get the element density
    rho = edata[:,e_var_names.split(',').index('dv1')]

    # Make a continuous element density value
    crho = np.zeros(cdata.shape[0])
    crho[conn] += rho

    # Normalize by the number of times the continuous density appears
    c = np.zeros(cdata.shape[0])
    c[conn] += np.ones(edata.shape[0])

    # Create a continuous density value
    crho /= c

    # Get the x/y coordinates
    x = cdata[:,0]
    y = cdata[:,1]

    # Assume that we have a mesh with all the same element type
    n = 2
    if ltypes[0] == TACS.QUAD_QUADRATIC_ELEMENT:
        n = 3
    elif ltypes[0] == TACS.QUAD_CUBIC_ELEMENT:
        n = 4
    elif ltypes[0] == TACS.QUAD_QUARTIC_ELEMENT:
        n = 5
    elif ltypes[0] == TACS.QUAD_QUINTIC_ELEMENT:
        n = 6

    # Specify as triangles
    triangles = []
    for elem in range(num_elements):
        for j in range(n-1):
            for i in range(n-1):
                triangles.append([conn[ptr[elem] + n*i + j],
                                  conn[ptr[elem] + n*i + j+1],
                                  conn[ptr[elem] + n*(i+1) + j+1]])
                triangles.append([conn[ptr[elem] + n*i + j],
                                  conn[ptr[elem] + n*(i+1) + j+1],
                                  conn[ptr[elem] + n*(i+1) + j]])

    # Create the triangles
    tri_obj = tri.Triangulation(x, y, triangles)

    # Plot the result as a figure
    fig, ax = plt.subplots()

    # Set the aspect ratio equal
    ax.set_aspect('equal')

    # Make sure that there are no ticks on either axis (these affect the bounding-box
    # and make extra white-space at the corners). Finally, turn off the axis so its
    # not visible.
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])
    ax.axis('off')

    # Set the number of levels to use.
    levels = np.linspace(0.0, 1.0, 26)

    # Create the contour plot
    ax.tricontourf(tri_obj, crho, levels, cmap='coolwarm', extend='max')

    # Save the figure. The bounding box is made tight to the figure, and the pading is
    # determined via the pad_inches argument. Transparent sets the background transparent.
    plt.savefig(outfile, dpi=500, transparent=True,
                bbox_inches='tight', pad_inches=0.01)

    # Close the figure
    plt.close()
