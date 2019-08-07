from __future__ import print_function, division
import numpy as np
from tmr import TMR
from egads4py import egads

# Read in the uCRM iges surfaces from final_surface.igs
ctx = egads.context()
filename = 'final_surface.igs'
oml = ctx.loadModel('final_surface.igs')

# Get the bodies and find the edge loops that are
# required to define
bodies = oml.getChildren()

edge_list = []

# Get the top and bottom surfaces
yloc = [0.0, 77.5, 276.5, 748.0]
surf_index = [0, 1, 0, 1, 2, 3, 5, 6]
edge_index = [1, 1, 3, 3, 3, 3, 3, 3]

# For each body in the shell
for s, e in zip(surf_index, edge_index):
    # Get the surface
    face = bodies[s].getChildren()[0]

    # Get the edge loop
    loop = face.getChildren()[0]

    # Get the edge
    edge = loop.getChildren()[e]

    # Get the edge
    edge_list.append(edge)

nctl = 25

top_curves = []
bottom_curves = []
for k in range(0, len(edge_list), 2):
    n = 250
    r, p = edge_list[k].getRange()
    xi = r[0] + (r[1] - r[0])*(0.5*(1.0 - np.cos(np.pi*np.linspace(0, 1.0, n))))
    x1 = np.zeros((n, 3))
    for i in range(n):
        x1[i], xt, xtt = edge_list[k].evaluate(xi[i])

    r, p = edge_list[k+1].getRange()
    xi = r[0] + (r[1] - r[0])*(0.5*(1.0 - np.cos(np.pi*np.linspace(0, 1.0, n))))
    # xi = xi[1:]
    x2 = np.zeros((n, 3))
    for i in range(n):
        x2[i], xt, xtt = edge_list[k+1].evaluate(xi[i])

    x = np.vstack((x1, x2))
    x[:,1] = yloc[k//2]

    interp = TMR.CurveInterpolation(x)
    interp.setNumControlPoints(nctl)
    curve = interp.createCurve(4)

    top, bottom = curve.split(0.5)
    top_curves.append(top)
    bottom_curves.append(bottom)

# Loft the curves
top_lofter = TMR.CurveLofter(top_curves)
top_surface = top_lofter.createSurface(2)

bottom_lofter = TMR.CurveLofter(bottom_curves)
bottom_surface = bottom_lofter.createSurface(2)

face_list = []

ku, kv, top_tu, ttv, wt, Xt = top_surface.getData()
ku, kv, bottom_tu, btv, wb, Xb = bottom_surface.getData()

for k in range(3):
    # Create the egads top surface
    top_tv = [ttv[k+1], ttv[k+1], ttv[k+2], ttv[k+2]]
    oclass = egads.SURFACE
    mtype = egads.BSPLINE
    bitflag = 2
    udegree = ku-1
    vdegree = kv-1
    nuknots = len(top_tu)
    nvknots = len(top_tv)
    ncpu = (len(top_tu) - ku)
    ncpv = (len(top_tv) - kv)
    idata = [bitflag,
             udegree, ncpu, nuknots,
             vdegree, ncpv, nvknots]
    w = wt[ncpu*k:ncpu*(k+2)]
    X = Xt[ncpu*k:ncpu*(k+2)]
    rdata = [top_tu, top_tv, X, w]
    top_surf = ctx.makeGeometry(oclass, mtype, rdata=rdata, idata=idata)

    # Create the egads bottom surface
    bottom_tv = [btv[k+1], btv[k+1], btv[k+2], btv[k+2]]
    oclass = egads.SURFACE
    mtype = egads.BSPLINE
    bitflag = 2
    udegree = ku-1
    vdegree = kv-1
    nuknots = len(bottom_tu)
    nvknots = len(bottom_tv)
    ncpu = (len(bottom_tu) - ku)
    ncpv = (len(bottom_tv) - kv)
    idata = [bitflag,
             udegree, ncpu, nuknots,
             vdegree, ncpv, nvknots]
    w = wb[ncpu*k:ncpu*(k+2)]
    X = Xb[ncpu*k:ncpu*(k+2)]
    rdata = [bottom_tu, bottom_tv, X, w]
    bottom_surf = ctx.makeGeometry(oclass, mtype, rdata=rdata, idata=idata)

    top_face = ctx.makeFace(top_surf, egads.SREVERSE,
                            rdata=[top_tu[0], top_tu[-1],
                                   top_tv[0], top_tv[-1]])
    bottom_face = ctx.makeFace(bottom_surf, egads.SFORWARD,
                               rdata=[bottom_tu[0], bottom_tu[-1],
                                      bottom_tv[0], bottom_tv[-1]])

    face_list.extend([top_face, bottom_face])

# Sew the faces together
model = ctx.sewFaces(face_list, toler=1e-3, manifold=False)

# Write out the model
model.saveModel('ucrm_9_oml.egads', overwrite=True)
model.saveModel('ucrm_9_oml.step', overwrite=True)
