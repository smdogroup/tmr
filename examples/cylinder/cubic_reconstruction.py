"""
This code tests the bi-cubic reconstruction performed within the
adaptive mesh refinement algorithm. This is used to smoothly
reconstrut the adjoint solution from the working mesh to the
refined mesh.
"""

import numpy as np
import matplotlib.pylab as plt


def func(x, y):
    """This is the true function"""

    return (x - 1) ** 3 + (y - 0.5) ** 5


def getShapeFuncs(xi, eta):
    """
    Compute the shape functions for the quadratic bi-Lagrange element.
    """

    # Compute the shape functions along each parametric direction
    nxi = np.array(
        [-0.5 * xi * (1.0 - xi), (1.0 - xi) * (1.0 + xi), 0.5 * xi * (1.0 + xi)]
    )
    neta = np.array(
        [-0.5 * eta * (1.0 - eta), (1.0 - eta) * (1.0 + eta), 0.5 * eta * (1.0 + eta)]
    )

    # Compute the derivative of the shape functions along each
    # parametric direction
    dxi = np.array([-0.5 + xi, -2.0 * xi, 0.5 + xi])
    deta = np.array([-0.5 + eta, -2.0 * eta, 0.5 + eta])

    # Evaluate the shape functions themselves
    N = np.outer(neta, nxi).flatten()
    Nxi = np.outer(neta, dxi).flatten()
    Neta = np.outer(deta, nxi).flatten()

    return N, Nxi, Neta


def getEnrichmentFuncs(xi, eta):
    """
    Return the enriched shape functions for the new element.
    """

    # Compute the enrichment function - note these are zero
    # at every parametric node location
    cxi = (1.0 + xi) * xi * (1.0 - xi)
    ceta = (1.0 + eta) * eta * (1.0 - eta)

    # Compute the derivative of the enrichment function
    dxi = 1.0 - 3.0 * xi * xi
    deta = 1.0 - 3.0 * eta * eta

    # Compute the enriched shape functions
    N = np.array(
        [cxi, cxi * eta, cxi * eta**2, ceta, ceta * xi, ceta * xi**2, cxi * ceta]
    )

    # Compute the derivative of the enriched shape functions
    Nxi = np.array(
        [dxi, dxi * eta, dxi * eta**2, 0.0, ceta, 2.0 * ceta * xi, dxi * ceta]
    )
    Neta = np.array(
        [0.0, cxi, 2.0 * cxi * eta, deta, deta * xi, deta * xi**2, cxi * deta]
    )

    return N, Nxi, Neta


def computeRecon(Xpts, uvals, uderiv):
    """
    Given the values at the nodes, and the derivative of the values at
    the nodes, compute the higher-order reconstruction term
    """

    # Form the additional 18 equations
    A = np.zeros((18, 7))
    b = np.zeros(18)

    c = 0
    for j in range(3):
        for i in range(3):
            # Set the parametric location
            xi = -1.0 + 1.0 * i
            eta = -1.0 + 1.0 * j

            # Compute the shape function
            N, Nxi, Neta = getShapeFuncs(xi, eta)

            # Compute the Jacobian transformation
            Xd = np.array(
                [
                    [np.dot(Nxi, Xpts[:, 0]), np.dot(Neta, Xpts[:, 0])],
                    [np.dot(Nxi, Xpts[:, 1]), np.dot(Neta, Xpts[:, 1])],
                ]
            )

            # Compute the inverse J = Xd^{-1}
            J = np.linalg.inv(Xd)

            # First, compute the contributions to the
            # righ-hand-side. The b vector contains the difference
            # between the prescribed derivative and the contribution
            # to the derivative from the quadratic shape function
            # terms

            # Set the derivative u,x
            b[c] = uderiv[3 * j + i, 0]
            b[c] -= J[0, 0] * np.dot(uvals, Nxi) + J[1, 0] * np.dot(uvals, Neta)

            # The derivative u,y
            b[c + 1] = uderiv[3 * j + i, 1]
            b[c + 1] -= J[0, 1] * np.dot(uvals, Nxi) + J[1, 1] * np.dot(uvals, Neta)

            # Now, set the terms from the enriched part of the
            # distribution
            N, Nxi, Neta = getEnrichmentFuncs(xi, eta)

            # For each row, cycle through and set the enrichment values
            for i in range(7):
                A[c, i] = J[0, 0] * Nxi[i] + J[1, 0] * Neta[i]
                A[c + 1, i] = J[0, 1] * Nxi[i] + J[1, 1] * Neta[i]

            c += 2

    # Scale the equations
    wvals = [0.5, 1.0, 0.5]
    weights = np.outer(wvals, wvals).flatten()

    for i in range(9):
        b[2 * i] *= weights[i]
        b[2 * i + 1] *= weights[i]
        A[2 * i, :] *= weights[i]
        A[2 * i + 1, :] *= weights[i]

    # Compute the solution to the least squares problem
    ubar, res, rank, s = np.linalg.lstsq(A, b)

    return ubar


def computeElemError(uvals, ubar, Xpts, nquad=5):
    """
    Compute the integral of the reconstruction error over the element
    """

    # Set the error variables to zero
    A = 0.0
    err1 = 0.0
    err2 = 0.0

    # Gauss quadrature scheme
    if nquad == 5:
        quadpts = [
            -0.906179845938664,
            -0.538469310105683,
            0.0,
            0.538469310105683,
            0.906179845938664,
        ]
        quadwts = [
            0.236926885056189,
            0.478628670499366,
            0.568888888888889,
            0.478628670499366,
            0.236926885056189,
        ]
    else:
        nquad = 6
        quadpts = [
            -0.9324695142031520278123016,
            -0.6612093864662645136613996,
            -0.2386191860831969086305017,
            0.2386191860831969086305017,
            0.6612093864662645136613996,
            0.9324695142031520278123016,
        ]
        quadwts = [
            0.1713244923791703450402961,
            0.3607615730481386075698335,
            0.4679139345726910473898703,
            0.4679139345726910473898703,
            0.3607615730481386075698335,
            0.1713244923791703450402961,
        ]

    for j in range(nquad):
        for i in range(nquad):
            # Set the parametric location
            xi = quadpts[i]
            eta = quadpts[j]

            # Compute the shape function
            N, Nxi, Neta = getShapeFuncs(xi, eta)

            # Compute the Jacobian transformation
            Xd = np.array(
                [
                    [np.dot(Nxi, Xpts[:, 0]), np.dot(Neta, Xpts[:, 0])],
                    [np.dot(Nxi, Xpts[:, 1]), np.dot(Neta, Xpts[:, 1])],
                ]
            )

            # Compute the inverse J = Xd^{-1}
            detJ = np.linalg.det(Xd)

            # Compute the x and y locations
            x = np.dot(N, Xpts[:, 0])
            y = np.dot(N, Xpts[:, 1])

            # Evaluate the value of the exact function
            uexact = func(x, y)

            # Now, set the terms from the enriched part of the
            # distribution
            Nr, nxi, neta = getEnrichmentFuncs(xi, eta)

            # Evaluate the approximate function and the corrected
            # function values
            u1 = np.dot(N, uvals)
            u2 = u1 + np.dot(Nr, ubar)

            A += detJ * quadwts[i] * quadwts[j]
            err1 += detJ * quadwts[i] * quadwts[j] * (uexact - u1) * (uexact - u1)
            err2 += detJ * quadwts[i] * quadwts[j] * (uexact - u2) * (uexact - u2)

    return A, err1, err2


def computeError(Xpts, nx, ny, dh=1e-30):

    # Parametric locations within the global mesh
    u = np.linspace(-1, 1, 2 * nx + 1)
    v = np.linspace(-1, 1, 2 * nx + 1)

    # Element-specific data
    elemXpts = np.zeros((9, 2))
    uvals = np.zeros(9)
    uderiv = np.zeros((9, 2))

    A = 0.0
    err1 = 0.0  # Interpolation error
    err2 = 0.0  # Reconstruction error

    # Loop over all the sub-elements in this mesh
    for j in range(ny):
        for i in range(nx):
            for jj in range(3):
                for ii in range(3):
                    # Set the parameter points
                    xi = u[2 * i + ii]
                    eta = v[2 * j + jj]

                    # Compute the nodal locations
                    N, Nxi, Neta = getShapeFuncs(xi, eta)
                    elemXpts[3 * jj + ii, 0] = np.dot(N, Xpts[:, 0])
                    elemXpts[3 * jj + ii, 1] = np.dot(N, Xpts[:, 1])

            # Set the values of the distribution at the nodes and the
            # derivative of u in the x/y directions at each of the
            # nodes in the element
            for ii in range(9):
                x = elemXpts[ii, 0]
                y = elemXpts[ii, 1]

                # Set the values of the function and its derivative
                uvals[ii] = func(x, y)
                uderiv[ii, 0] = func(x + 1j * dh, y).imag / dh
                uderiv[ii, 1] = func(x, y + 1j * dh).imag / dh

            # Compute the reconstruction
            ubar = computeRecon(elemXpts, uvals, uderiv)

            # Add the total error
            A0, e1, e2 = computeElemError(uvals, ubar, elemXpts)
            A += A0
            err1 += e1
            err2 += e2

    return np.sqrt(err1), np.sqrt(err2)


# Set the position of the initial element in the domain
globalXpts = 2.0 * np.array(
    [
        0.25,
        0.25,
        0.5,
        0.35,
        1.0,
        0.4,
        0.15,
        0.35,
        0.45,
        0.45,
        0.9,
        0.5,
        0.05,
        0.45,
        0.3,
        0.5,
        0.85,
        0.6,
    ]
).reshape((9, 2))

# Set the number of test cases to run
num = 20
interp_err = np.zeros(num)
recon_err = np.zeros(num)
nx = np.arange(1, num + 1)

# Set some reasonable colours
tableau20 = [
    (31, 119, 180),
    (174, 199, 232),
    (255, 127, 14),
    (255, 187, 120),
    (44, 160, 44),
    (152, 223, 138),
    (214, 39, 40),
    (255, 152, 150),
    (148, 103, 189),
    (197, 176, 213),
    (140, 86, 75),
    (196, 156, 148),
    (227, 119, 194),
    (247, 182, 210),
    (127, 127, 127),
    (199, 199, 199),
    (188, 189, 34),
    (219, 219, 141),
    (23, 190, 207),
    (158, 218, 229),
]


# Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.
for i in range(len(tableau20)):
    r, g, b = tableau20[i]
    tableau20[i] = (r / 255.0, g / 255.0, b / 255.0)


for i in range(num):
    print("Computing error ", i)
    interp_err[i], recon_err[i] = computeError(globalXpts, nx[i], nx[i])

alpha = 0.15
x = np.array([0.5, 0.025])

plt.figure()
plt.loglog(
    1.0 / (2 * nx + 1),
    interp_err,
    "-o",
    color=tableau20[0],
    linewidth=2,
    label="interpolation",
)
plt.loglog(
    1.0 / (2 * nx + 1),
    recon_err,
    "-s",
    color=tableau20[2],
    linewidth=2,
    label="reconstruction",
)
plt.loglog(x, (x**3) / alpha, "-", color=tableau20[6], linewidth=2, label="3rd order")
plt.loglog(x, (x**4) / alpha, "-", color=tableau20[8], linewidth=2, label="4th order")
plt.xlabel("Average $\\sf\\Delta x$", fontweight="bold", fontsize=18)
plt.ylabel("Error", fontweight="bold", fontsize=18)

plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.gcf().subplots_adjust(bottom=0.13)

plt.legend(loc=0)
plt.savefig("planar_reconstruction.pdf")
