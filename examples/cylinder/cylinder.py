from __future__ import print_function
from mpi4py import MPI
from tmr import TMR
from tacs import TACS, elements, constitutive, functions
import numpy as np
import argparse
import os
import ksFSDT

def integrate(integrand):
    sigma = [17.0/48.0, 59.0/48.0, 43.0/48, 49/48.9]
    r = len(sigma)

    integral = 0.0
    for i, s in enumerate(sigma):
        integral += s*integrand[i]
        integral += s*integrand[-1-i]

    for i in range(r, len(integrand)-r):
        integral += integrand[i]

    return integral

def cylinder_ks_functional(rho, t, E, nu, kcorr, ys,
                           L, R, alpha, beta, load, n=1000):
    A = np.zeros((5, 5))
    rhs = np.zeros(5)
    rhs[2] = -load
    ainv = 1.0/R

    # Compute the shear modulus
    G = 0.5*E/(1.0 + nu)

    Q11 = E/(1.0 - nu*nu)
    Q12 = nu*Q11
    Q22 = Q11
    Q33 = G

    # In-plane stiffness
    A11 = t*Q11
    A12 = t*Q12
    A22 = t*Q22
    A33 = t*G

    # Bending stiffness
    I = (t**3)/12.0
    D11 = Q11*I
    D12 = Q12*I
    D22 = Q22*I
    D33 = Q33*I

    # Shear stiffness
    bA11 = kcorr*t*G
    bA22 = kcorr*t*G

    # The first equation for u
    A[0,0] = -(A11*beta*beta + A33*alpha*alpha +
               D33*ainv*ainv*alpha*alpha)
    A[1,0] = -(A33 + A12)*alpha*beta
    A[2,0] = A12*beta*ainv
    A[3,0] = D33*ainv*alpha*alpha
    A[4,0] = D33*ainv*alpha*beta

    # The second equation for v
    A[0,1] = -(A12 + A33)*alpha*beta;
    A[1,1] = -(A33*beta*beta + A22*alpha*alpha + ainv*ainv*bA11 +
               D22*ainv*ainv*alpha*alpha)
    A[2,1] = (A22 + bA11)*ainv*alpha + D22*alpha*ainv*ainv*ainv
    A[3,1] = D12*ainv*alpha*beta
    A[4,1] = bA11*ainv + D22*ainv*alpha*alpha

    # The third equation for w
    A[0,2] = A12*beta*ainv
    A[1,2] = (bA11 + A22)*alpha*ainv + D22*alpha*ainv*ainv*ainv
    A[2,2] = -(bA11*alpha*alpha + bA22*beta*beta + A22*ainv*ainv +
               D22*ainv*ainv*ainv*ainv)
    A[3,2] = -bA22*beta - D12*beta*ainv*ainv
    A[4,2] = -bA11*alpha - D22*alpha*ainv*ainv

    # Fourth equation for theta
    A[0,3] = D33*ainv*alpha*alpha
    A[1,3] = D12*ainv*alpha*beta
    A[2,3] = -bA22*beta - D12*beta*ainv*ainv
    A[3,3] = -(D11*beta*beta + D33*alpha*alpha) - bA22
    A[4,3] = -(D12 + D33)*alpha*beta

    # Fifth equation for phi
    A[0,4] =  D33*ainv*alpha*beta
    A[1,4] =  bA11*ainv + D22*ainv*alpha*alpha
    A[2,4] = -bA11*alpha - D22*alpha*ainv*ainv
    A[3,4] = -(D33 + D12)*alpha*beta
    A[4,4] = -(D33*beta*beta + D22*alpha*alpha) - bA11

    # Solve for the coefficients
    ans = np.linalg.solve(A, rhs)

    U = ans[0]
    V = ans[1]
    W = ans[2]
    theta = ans[3]
    phi = ans[4]

    z1 = 0.5*t
    rz1 = 1.0 - z1*ainv
    a1 = -(beta*Q11*(U + z1*theta) + Q12*(alpha*(V*rz1 + z1*phi) - W*rz1*ainv))/ys
    b1 = -(beta*Q12*(U + z1*theta) + Q22*(alpha*(V*rz1 + z1*phi) - W*rz1*ainv))/ys
    c1 = Q33*(alpha*(U*rz1 + z1*theta) + beta*(V + z1*phi))/ys

    z2 = -0.5*t
    rz2 = 1.0 - z2*ainv
    a2 = -(beta*Q11*(U + z2*theta) + Q12*(alpha*(V*rz2 + z2*phi) - W*rz2*ainv))/ys
    b2 = -(beta*Q12*(U + z2*theta) + Q22*(alpha*(V*rz2 + z2*phi) - W*rz2*ainv))/ys
    c2 = Q33*(alpha*(U*rz2 + z2*theta) + beta*(V + z2*phi))/ys

    x = (a1**2 + b1**2 - a1*b1)
    y = 3.0*c1**2
    vm_max = max(np.sqrt(x*y/(x + y)), np.sqrt(x), np.sqrt(y))

    x = (a2**2 + b2**2 - a2*b2)
    y = 3.0*c2**2
    vm_max = max(vm_max, np.sqrt(x*y/(x + y)), np.sqrt(x), np.sqrt(y))

    x = np.linspace(0.0, L, n)
    y = np.linspace(0.0, 2.0*np.pi*R, n)
    S = np.zeros((n, n))

    tcoef1 = (a1**2 + b1**2 - a1*b1)
    tcoef2 = 3*c1**2
    bcoef1 = (a2**2 + b2**2 - a2*b2)
    bcoef2 = 3*c2**2

    for j in range(n):
        for i in range(n):
            sin2 = (np.sin(alpha*y[j])*np.sin(beta*x[i]))**2
            cos2 = (np.cos(alpha*y[j])*np.cos(beta*x[i]))**2
            vm1 = np.sqrt((tcoef1*sin2 + tcoef2*cos2))
            vm2 = np.sqrt((bcoef1*sin2 + bcoef2*cos2))
            S[i,j] = np.exp(rho*(max(vm1, vm2) - vm_max))

    # Integrate over the area
    h = (2.0*L*R*np.pi)/(n-1)**2
    ks_sum = 0.0
    for j in range(n-1):
        for i in range(n-1):
            ks_sum += (S[i,j] + S[i+1,j] + S[i,j+1] + S[i+1,j+1])
    ks_sum = 0.25*h*ks_sum

    return vm_max, vm_max + np.log(ks_sum)/rho

def disk_ks_functional(rho, t, E, nu, kcorr, ys, R, load, n=1000):
    '''
    Axisymmetric disk subject to a uniform pressure load
    '''
    D = ((t**3)/12.0)*E/(1.0 - nu**2)
    G = 0.5*E/(1.0 + nu)

    Q11 = E/(1.0 - nu**2)
    Q12 = nu*Q11
    Q22 = Q11

    # Compute the value of the radii
    r0 = np.linspace(0, R, n)
    S = np.zeros((n, 2))

    phi_over_r0 = load*(R**2/(16*D))
    dphidr0 = load*(R**2/(16*D))

    # Compute the value of the stress at the top surface
    z1 = 0.5*t
    err = z1*dphidr0
    ett = z1*phi_over_r0

    srr = Q11*err + Q12*ett
    stt = Q12*err + Q22*ett

    # Compute the von Mises at the top surface
    S[0,0] = np.sqrt(srr**2 + stt**2 - srr*stt)/ys

    # Compute the value of the stress at the bottom surface
    z2 = -0.5*t
    err = z2*dphidr0
    ett = z2*phi_over_r0

    srr = Q11*err + Q12*ett
    stt = Q12*err + Q22*ett

    # Compute the von Mises at the bottom surface
    S[0,1] = np.sqrt(srr**2 + stt**2 - srr*stt)/ys

    for i, r in enumerate(r0):
        if i == 0:
            continue

        # The value of the rotation
        phi = load*(R**3/(16*D))*(r/R)*(1 - (r/R)**2)
        dphidr = load*(R**2/(16*D))*(1 - 3*(r/R)**2)

        # Compute the value of the stress at the top surface
        z1 = 0.5*t
        err = z1*dphidr
        ett = z1*phi/r

        srr = Q11*err + Q12*ett
        stt = Q12*err + Q22*ett

        # Compute the von Mises at the top surface
        S[i,0] = np.sqrt(srr**2 + stt**2 - srr*stt)/ys

        # Compute the value of the stress at the bottom surface
        z2 = -0.5*t
        err = z2*dphidr
        ett = z2*phi/r

        srr = Q11*err + Q12*ett
        stt = Q12*err + Q22*ett

        # Compute the von Mises at the bottom surface
        S[i,1] = np.sqrt(srr**2 + stt**2 - srr*stt)/ys

    # Compute the maximum von Mises stress
    vm_mx = np.max(S)

    # Compute the contribution to the KS functional
    integrand = np.zeros(n)
    for i in range(n):
        integrand[i] = 2*np.pi*r0[i]*(np.exp(rho*(S[i,0] - vm_mx)) +
                                      np.exp(rho*(S[i,1] - vm_mx)))

    ks_sum = (R/(n-1))*integrate(integrand)

    return vm_mx, vm_mx + np.log(ks_sum)/rho

class CreateMe(TMR.QuadCreator):
    def __init__(self, bcs, case='cylinder'):
        TMR.QuadCreator.__init__(bcs)
        self.case = case
        return

    def createElement(self, order, quad):
        '''Create the element'''
        tmin = 0.0
        tmax = 1e6
        tnum = quad.tag
        stiff = ksFSDT.ksFSDT(ksweight, density, E, nu, kcorr, ys, t,
                              tnum, tmin, tmax)
        if self.case == 'cylinder':
            stiff.setRefAxis(np.array([0.0, 0.0, 1.0]))
        elif self.case == 'disk':
            stiff.setRefAxis(np.array([1.0, 0.0, 0.0]))
        elem = elements.MITCShell(order, stiff)
        return elem

def addFaceTraction(case, order, assembler, load):
    # Create the surface traction
    aux = TACS.AuxElements()

    # Get the element node locations
    nelems = assembler.getNumElements()

    if case == 'cylinder':
        for i in range(nelems):
            # Get the information about the given element
            elem, Xpt, vars0, dvars, ddvars = assembler.getElementData(i)

            # Loop over the nodes and create the traction forces in the
            # x/y/z directions
            nnodes = order*order
            tx = np.zeros(nnodes, dtype=TACS.dtype)
            ty = np.zeros(nnodes, dtype=TACS.dtype)
            tz = np.zeros(nnodes, dtype=TACS.dtype)

            # Set the components of the surface traction
            for j in range(nnodes):
                x = Xpt[3*j]
                y = Xpt[3*j+1]
                z = Xpt[3*j+2]
                theta = -R*np.arctan2(y, x)
                p = -load*np.sin(beta*z)*np.sin(alpha*theta)

                tx[j] = p*Xpt[3*j]/R
                ty[j] = p*Xpt[3*j+1]/R
                tz[j] = 0.0

            trac = elements.ShellTraction(order, tx, ty, tz)
            aux.addElement(i, trac)
    elif case == 'disk':
        nnodes = order*order
        tx = np.zeros(nnodes, dtype=TACS.dtype)
        ty = np.zeros(nnodes, dtype=TACS.dtype)
        tz = np.zeros(nnodes, dtype=TACS.dtype)
        tz[:] = load
        trac = elements.ShellTraction(order, tx, ty, tz)

        for i in range(nelems):
            aux.addElement(i, trac)

    return aux

def createRefined(case, forest, bcs, pttype=TMR.UNIFORM_POINTS):
    new_forest = forest.duplicate()
    new_forest.setMeshOrder(forest.getMeshOrder()+1, pttype)
    creator = CreateMe(bcs, case=case)
    return new_forest, creator.createTACS(new_forest)

def createProblem(case, forest, bcs, ordering, order=2, nlevels=2,
                  pttype=TMR.UNIFORM_POINTS):
    # Create the forest
    forests = []
    assemblers = []

    # Create the trees, rebalance the elements and repartition
    forest.balance(1)
    forest.setMeshOrder(order, pttype)
    forest.repartition()
    forests.append(forest)

    # Make the creator class
    creator = CreateMe(bcs, case=case)
    assemblers.append(creator.createTACS(forest, ordering))

    while order > 2:
        order = order-1
        forest = forests[-1].duplicate()
        forest.setMeshOrder(order, pttype)
        forest.balance(1)
        forests.append(forest)

        # Make the creator class
        creator = CreateMe(bcs, case=case)
        assemblers.append(creator.createTACS(forest, ordering))

    for i in range(nlevels-1):
        forest = forests[-1].coarsen()
        forest.setMeshOrder(2, pttype)
        forest.balance(1)
        forests.append(forest)

        # Make the creator class
        creator = CreateMe(bcs, case=case)
        assemblers.append(creator.createTACS(forest, ordering))

    # Create the multigrid object
    mg = TMR.createMg(assemblers, forests, omega=0.5)

    return assemblers[0], mg

# Set the communicator
comm = MPI.COMM_WORLD

# Create an argument parser to read in arguments from the commnad line
p = argparse.ArgumentParser()
p.add_argument('--case', type=str, default='cylinder')
p.add_argument('--steps', type=int, default=5)
p.add_argument('--htarget', type=float, default=10.0)
p.add_argument('--ordering', type=str, default='multicolor')
p.add_argument('--order', type=int, default=2)
p.add_argument('--ksweight', type=float, default=100.0)
p.add_argument('--uniform_refinement', action='store_true', default=False)
p.add_argument('--structured', action='store_true', default=False)
p.add_argument('--energy_error', action='store_true', default=False)
args = p.parse_args()

# Set the case type
case = args.case

# Parameters that define the geometry of the cylinder
t = 1.0 # 1 mm thickness
L = 100.0 # 100 mm length
R = 100.0/np.pi # Radius of the cylinder (in mm)

# Set the load to apply to the cylinder
load = 3.0 # 3 MPa
alpha = 4.0/R # Set the value of alpha
beta = 3.0*np.pi/L # Set the value of beta

# Set the material properties
density = 2700.0e-9 # kg/mm^3
E = 70e3 # 70 GPa
nu = 0.3 # Poisson's ratio
ys = 350.0 # 350 MPa
kcorr = 5.0/6.0 # The shear correction factor

# Set the type of ordering to use for this problem
ordering = args.ordering
ordering = ordering.lower()

# Set the value of the target length scale in the mesh
htarget = args.htarget

# Set the order
order = args.order

# Set the KS parameter
ksweight = args.ksweight

# Set the number of AMR steps to use
steps = args.steps

# Store whether to use a structured/unstructed mesh
structured = args.structured

# Set the load value to use depending on the case
exact_ks_functional = 0.0

# The boundary condition object
bcs = TMR.BoundaryConditions()

if case == 'cylinder':
    # Get the maximum stress and re-adjust the load
    vm_max, res = cylinder_ks_functional(1.0, t, E, nu, kcorr, ys,
        L, R, alpha, beta, load, n=10)
    load = load/vm_max

    # Compute the exact KS value
    one, exact_ks_functional = cylinder_ks_functional(ksweight,
        t, E, nu, kcorr, ys,
        L, R, alpha, beta, load, n=1000)

    # Load the geometry model
    geo = TMR.LoadModel('cylinder.stp')
    verts = geo.getVertices()
    edges = geo.getEdges()
    faces = geo.getFaces()

    # Create the simplified geometry with only faces
    geo = TMR.Model(verts, edges, [faces[0]])

    # Set the boundary conditions
    verts[0].setAttribute('Clamped')
    verts[1].setAttribute('Restrained')
    edges[0].setAttribute('Restrained')
    edges[2].setAttribute('Restrained')

    # Set the boundary conditions
    bcs.addBoundaryCondition('Clamped', [0, 1, 2, 5])
    bcs.addBoundaryCondition('Restrained', [0, 1, 5])
elif case == 'disk':
    # Get the maximum stress and re-adjust the load
    vm_max, res = disk_ks_functional(1.0, t, E, nu, kcorr, ys, R, load, n=10)
    load = load/vm_max

    # Compute the exact KS value
    # for n in [1000, 10000, 100000, 1000000, 10000000]:
    #     one, exact_ks_functional = disk_ks_functional(ksweight,
    #         t, E, nu, kcorr, ys, R, load, n=n)
    #     print('exact ks functional = ', exact_ks_functional)

    exact_ks_functional = 1.0372627248381336

    # Load the geometry model
    geo = TMR.LoadModel('cylinder.stp')
    verts = geo.getVertices()
    edges = geo.getEdges()
    faces = geo.getFaces()

    # Set the attributes
    verts[1].setAttribute('clamped')
    edges[2].setAttribute('clamped')

    # Set the boundary conditions
    bcs.addBoundaryCondition('clamped', [0, 1, 2, 3, 4, 5])

    # Create the new model
    geo = TMR.Model([verts[1]], [edges[2]], [faces[2]])

if comm.rank == 0:
    print('Exact KS functional = %25.15e'%(exact_ks_functional))

# Create the new mesh
mesh = TMR.Mesh(comm, geo)

# Set the meshing options
opts = TMR.MeshOptions()

# Set the mesh type
opts.mesh_type_default = TMR.UNSTRUCTURED
if structured:
    opts.mesh_type_default = TMR.STRUCTURED

opts.frontal_quality_factor = 1.25
opts.num_smoothing_steps = 50
opts.write_mesh_quality_histogram = 1
opts.triangularize_print_iter = 50000

# Create the surface mesh
mesh.mesh(htarget, opts)
mesh.writeToVTK('mesh.vtk')

# Create the corresponding mesh topology from the mesh-model
model = mesh.createModelFromMesh()
model.writeModelToTecplot('model.dat')
topo = TMR.Topology(comm, model)

# Create the quad forest and set the topology of the forest
depth = 0
if order == 2:
    depth = 1
forest = TMR.QuadForest(comm)
forest.setTopology(topo)
forest.setMeshOrder(order, TMR.UNIFORM_POINTS)
forest.createTrees(depth)

# Set the ordering to use
if ordering == 'rcm':
    ordering = TACS.PY_RCM_ORDER
elif ordering == 'multicolor':
    ordering = TACS.PY_MULTICOLOR_ORDER
else:
    ordering = TACS.PY_NATURAL_ORDER

# The target relative error
target_rel_err = 1e-4

# Open a log file to write
descript = 'unstructured'
if structured:
    descript = 'structured'
if args.uniform_refinement:
    descript += '_uniform'

# Create the log file and write out the header
log_fp = open('%s_order%d_%s.dat'%(case, order, descript), 'w')
s = 'Variables = iter, nelems, nnodes, fval, fcorr, abs_err, adjoint_corr, '
s += 'exact, fval_error, fval_corr_error, '
s += 'fval_effectivity, indicator_effectivity\n'
log_fp.write(s)

for k in range(steps):
    # Create the topology problem
    nlevs = min(3, depth+k+1)
    assembler, mg = createProblem(case, forest, bcs, ordering,
                                  order=order, nlevels=nlevs)
    aux = addFaceTraction(case, order, assembler, load)
    assembler.setAuxElements(aux)

    # Create the assembler object
    res = assembler.createVec()
    ans = assembler.createVec()

    use_direct_solve = False
    if use_direct_solve:
        mat = assembler.createFEMat()
        assembler.assembleJacobian(1.0, 0.0, 0.0, res, mat)

        # Create the direct solver
        pc = TACS.Pc(mat)
        pc.factor()
    else:
        # Solve the linear system
        mg.assembleJacobian(1.0, 0.0, 0.0, res)
        mg.factor()
        pc = mg
        mat = mg.getMat()

    gmres = TACS.KSM(mat, pc, 100, isFlexible=0)
    gmres.setMonitor(comm, freq=10)
    gmres.setTolerances(1e-14, 1e-30)
    gmres.solve(res, ans)
    ans.scale(-1.0)

    # Set the variables
    assembler.setVariables(ans)

    # Output for visualization
    flag = (TACS.ToFH5.NODES |
            TACS.ToFH5.DISPLACEMENTS |
            TACS.ToFH5.STRAINS |
            TACS.ToFH5.EXTRAS)
    f5 = TACS.ToFH5(assembler, TACS.PY_SHELL, flag)
    f5.writeToFile('results/solution%02d.f5'%(k))

    # Create the function
    func = functions.KSFailure(assembler, ksweight)
    func.setKSFailureType('continuous')
    fval = assembler.evalFunctions([func])[0]

    # Create the refined mesh
    forest_refined, assembler_refined = createRefined(case, forest, bcs)
    aux = addFaceTraction(case, order+1, assembler_refined, load)
    assembler_refined.setAuxElements(aux)

    if args.energy_error:
        # Compute the strain energy error estimate
        fval_corr = 0.0
        adjoint_corr = 0.0
        err_est, error = TMR.strainEnergyError(forest, assembler,
            forest_refined, assembler_refined)
    else:
        # Compute the adjoint
        res.zeroEntries()
        assembler.evalSVSens(func, res)

        # Compute the adjoint solution
        adjoint = assembler.createVec()
        gmres.solve(res, adjoint)
        adjoint.scale(-1.0)

        # Compute the adjoint and use adjoint-based refinement
        err_est, adjoint_corr, error = \
            TMR.adjointError(forest, assembler,
                             forest_refined, assembler_refined, adjoint)

        # Evaluate the function again on the refined mesh to obtain
        # the
        func_refined = functions.KSFailure(assembler_refined, ksweight)
        func_refined.setKSFailureType('continuous')
        fval_refined = assembler_refined.evalFunctions([func_refined])[0]

        # Compute the refined function value
        fval_corr = fval_refined + adjoint_corr

    f5_refine = TACS.ToFH5(assembler_refined, TACS.PY_SHELL, flag)
    f5_refine.writeToFile('results/solution_refined%02d.f5'%(k))

    # Compute the refinement from the error estimate
    low = -16
    high = 4
    bins_per_decade = 4
    nbins = bins_per_decade*(high - low)
    bounds = 10**np.linspace(high, low, nbins+1)
    bins = np.zeros(nbins+2, dtype=np.int)

    # Compute the mean and standard deviations of the log(error)
    ntotal = comm.allreduce(assembler.getNumElements(), op=MPI.SUM)
    mean = comm.allreduce(np.sum(np.log(error)), op=MPI.SUM)
    mean /= ntotal

    # Compute the standard deviation
    stddev = comm.allreduce(np.sum((np.log(error) - mean)**2), op=MPI.SUM)
    stddev = np.sqrt(stddev/(ntotal-1))

    # Get the total number of nodes
    nnodes = comm.allreduce(assembler.getNumOwnedNodes(), op=MPI.SUM)

    # Compute the error from the exact solution
    fval_error = np.fabs(exact_ks_functional - fval)
    fval_corr_error = np.fabs(exact_ks_functional - fval_corr)

    fval_effectivity = (fval_corr - fval)/(exact_ks_functional - fval)
    indicator_effectivity = err_est/np.fabs(exact_ks_functional - fval)

    # Write the log data to a file
    s = '%6d %6d %6d %20.15e %20.15e %20.15e %20.15e '
    s += '%20.15e %20.15e %20.15e %20.15e %20.15e\n'
    log_fp.write(
        s%(k, ntotal, nnodes, fval, fval_corr, err_est, adjoint_corr,
        exact_ks_functional, fval_error, fval_corr_error,
        fval_effectivity, indicator_effectivity))
    log_fp.flush()

    # Compute the bins
    for i in range(len(error)):
        if error[i] > bounds[0]:
            bins[0] += 1
        elif error[i] < bounds[-1]:
            bins[-1] += 1
        else:
            for j in range(len(bounds)-1):
                if (error[i] <= bounds[j] and
                    error[i] > bounds[j+1]):
                    bins[j+1] += 1

    # Compute the number of bins
    bins = comm.allreduce(bins, MPI.SUM)

    # Compute the sum of the bins
    total = np.sum(bins)

    # Print out the result
    if comm.rank == 0:
        print('fval      = ', fval)
        print('fval corr = ', fval_corr)
        print('estimate  = ', err_est)
        print('mean      = ', mean)
        print('stddev    = ', stddev)

        print('%10s   %10s   %12s   %12s'%(
            'high', 'low', 'bins', 'percentage'))
        print('%10.2e   %10s   %12d   %12.2f'%(
            bounds[0], ' ', bins[0], 100.0*bins[0]/total))
        for i in range(nbins):
            print('%10.2e   %10.2e   %12d   %12.2f'%(
                bounds[i], bounds[i+1], bins[i+1], 100.0*bins[i+1]/total))
        print('%10s   %10.2e   %12d   %12.2f'%(
            ' ', bounds[-1], bins[-1], 100.0*bins[-1]/total))

        # Set the data
        data = np.zeros((nbins, 4))
        for i in range(nbins-1, -1, -1):
            data[i,0] = bounds[i]
            data[i,1] = bounds[i+1]
            data[i,2] = bins[i+1]
            data[i,3] = 100.0*bins[i+1]/total
        np.savetxt('results/cylinder_data%d.txt'%(k), data)

    # Print out the error estimate
    assembler.setDesignVars(error)
    f5.writeToFile('results/error%02d.f5'%(k))

    # Perform the refinement
    if args.uniform_refinement:
        forest.refine()
    elif k < steps-1:
        # The refinement array
        refine = np.zeros(len(error), dtype=np.intc)

        # Compute the target relative error
        element_target_error = target_rel_err*fval/ntotal
        log_elem_target_error = np.log(element_target_error)

        # Determine the cutoff values
        cutoff = bins[-1]
        bin_sum = 0
        for i in range(len(bins)+1):
            bin_sum += bins[i]
            if bin_sum > 0.15*ntotal:
                cutoff = bounds[i]
                break

        log_cutoff = np.log(cutoff)

        if comm.rank == 0:
            print('elem_target_error =        %15.3e'%(
                element_target_error))
            print('log10(elem_target_error) = %15.3e'%(
                np.log10(element_target_error)))
            print('log10(mean) =              %15.3e'%(
                np.log10(np.exp(mean))))
            print('log10(stddev) =            %15.3e'%(
                np.log10(np.exp(stddev))))
            print('log_elem_target_error =    %15.3e'%(
                log_elem_target_error))
            print('mean =                     %15.3e'%(mean))
            print('stddev =                   %15.3e'%(stddev))
            print('cutoff =                   %15.3e'%(cutoff))

        # Element target error is still too high. Adapt based
        # solely on decreasing the overall error
        nrefine = 0
        for i, err in enumerate(error):
            # Compute the log of the error
            logerr = np.log(err)

            if logerr > log_cutoff:
                refine[i] = 1
                nrefine += 1

        # Refine the forest
        forest.refine(refine)
