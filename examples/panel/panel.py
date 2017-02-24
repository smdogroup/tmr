import numpy as np
from tmr import TMR


def create_panel(Lx, Ly, hz):
    nu = 4
    nv = 4
    x = np.linspace(-0.5*Lx, 0.5*Lx, nu)
    y = np.linspace(-0.5*Ly, 0.5*Ly, nv)

    pts = np.zeros((nu, nv, 3))
    for j in xrange(nv):
        for i in xrange(nu):
            pts[i,j,0] = x[i]
            pts[i,j,1] = y[j]

    panel_surf = TMR.BsplineSurface(pts)
    panel_surf.writeToVTK('panel_surface.vtk')


    # Create the blade stiffeners
    nz = 2
    z = np.linspace(0, -hz, nz)
    yloc = [-0.25*Ly, 0.25*Ly]

    blades = []
    for k in xrange(len(yloc)):
        pts = np.zeros((nu, nz, 3))
        for i in xrange(nu):
            for j in xrange(nz):
                pts[i,j,0] = x[i]
                pts[i,j,1] = yloc[k]
                pts[i,j,2] = z[j]

        blades.append(TMR.BsplineSurface(pts))
        blades[k].writeToVTK('blade%d_surface.vtk'%(k))

create_panel(1.0, 2.0, 0.1)
