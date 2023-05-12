'''
A simple example to test adaptive analysis with TACS/TMR. 

This example uses a square flat plate, fully fixed on all edges, with a uniform
pressure load. Static analysis is performed, followed by either uniform or 
adaptive refinement. Adaptive refinement is based on a selected output of 
interest. Various options for this example can be selected through command line 
arguments, see below. 

To run a simple adaptive analysis on the KS failure output, do the following:
1. run `python make_plate_geom.py` to generate the square plate geometry file
2. run `mpirun -np 4 python plate_adapt.py --niters 8 --adapt`
3. run `mpirun -np 4 python plate_adapt.py --niters 5`
4. compare the results of adaptive and uniform refinement between step 2 and 3
'''
import sys
import os
import argparse
from mpi4py import MPI
import numpy as np
from tmr import TMR
from tacs import TACS, elements, constitutive, functions, problems
from tmr.pytacsadapt import *

#===============================================================================
# Local functions
#===============================================================================
def elemCallBack(order, quad):
	# set constitutive properties
	rho = 2700.
	E = 70.e9
	nu = 0.33
	ys = 276.e6
	min_thickness = 0.001
	max_thickness = 0.01
	thickness = 0.01
	props = constitutive.MaterialProperties(rho=rho, E=E, nu=nu, ys=ys)
	stiff = constitutive.IsoShellConstitutive(props, t=thickness, tlb=min_thickness, tub=max_thickness)

	# create the element
	transform = elements.ShellNaturalTransform()
	element = None
	if order == 2:
		element = elements.Quad4Shell(transform, stiff)
	elif order == 3:
		element = elements.Quad9Shell(transform, stiff)
	elif order == 4:
		element = elements.Quad16Shell(transform, stiff)
	else:
		raise ValueError(f'orders {2}-{4} supported')
	return element

def probCallBack(name, comm, geomLoader, assembler, **kwargs):
	# assign keyword arguments
	args = kwargs.get('cmd_ln_args')

	# set the output flags
	flags = (TACS.OUTPUT_CONNECTIVITY |
			 TACS.OUTPUT_NODES |
			 TACS.OUTPUT_DISPLACEMENTS |
			 TACS.OUTPUT_STRAINS |
			 TACS.OUTPUT_STRESSES |
			 TACS.OUTPUT_EXTRAS)
	
	# create the output viewer
	f5Output = TACS.ToFH5(assembler, TACS.BEAM_OR_SHELL_ELEMENT, flags)
	for ind, comp_name in enumerate(geomLoader.getComponentDescripts()):
		f5Output.setComponentName(ind, comp_name)
	
	# set static problem options
	static_opts = {'subSpaceSize' : 10,
				   'nRestarts' : 1,
				   'L2Convergence' : 1.e-9,
				   'L2ConvergenceRel' : 1.e-12,
				   'useMonitor' : True,
				   'monitorFrequency' : 1}
	
	# create a static problem (meshLoader arg defaults to None)
	problem = problems.StaticProblem(name, assembler, comm, outputViewer=f5Output, options=static_opts)

	# add potential functions of interest to this problem
	problem.addFunction(funcName='ks_disp', funcHandle=functions.KSDisplacement, ftype='continuous', ksWeight=args.ksweight, direction=[0., 0., 1.])
	problem.addFunction(funcName='ks_fail', funcHandle=functions.KSFailure, ftype='continuous', ksWeight=args.ksweight)
	problem.addFunction(funcName='compliance', funcHandle=functions.Compliance)
	return problem

def initCallBack(coarse_space, **kwargs):
	# assign keyword arguments
	geom_file = kwargs.get('geom_file')
	geom_type = kwargs.get('geom_type')
	args = kwargs.get('cmd_ln_args')

	# load the geometry for the coarse model
	coarse_space.geomLoader = pyGeometryLoader(coarse_space.comm, geom_type)
	coarse_space.geomLoader.readGeometryFile(geom_file)

	# name the face and edges
	face_names = {0:'panel'}
	edge_names = {0:'south',1:'east',2:'north',3:'west'}
	coarse_space.geomLoader.nameGeometricEntities(face_names=face_names, edge_names=edge_names)

	# specify boundary conditions - use defaults for nums and vals to fully fix
	bc_names = ['south','east','north','west']
	bc_nums = {}
	bc_info = {'type' : 'edge',
			   'bc_names' : bc_names,
			   'bc_nums' : bc_nums}
	coarse_space.geomLoader.setBoundaryConditions(bc_info)

	# set the meshing options
	coarse_space.geomLoader.setMeshOptions(write_mesh_quality_histogram=0)
	if not args.adapt:
		coarse_space.geomLoader.setMeshOptions(mesh_type_default=TMR.STRUCTURED)

	# create the initial coarse mesh
	coarse_space.geomLoader.createMesh(h=args.hinit)

	# create the topology 
	coarse_space.geomLoader.createTopology(useMesh=True)

	# create the coarse forest
	if args.order < 4:
		coarse_space.geomLoader.createForest(order=args.order, interp=TMR.UNIFORM_POINTS)
	else: 
		coarse_space.geomLoader.createForest(order=args.order, interp=TMR.GAUSS_LOBATTO_POINTS)
	coarse_space.forest = coarse_space.geomLoader.getForest()
	return

def printroot(msg):
	if comm.rank == 0:
		print(f'{msg}')
	return

#===============================================================================
# Main
#===============================================================================
if __name__ == '__main__':
	# get the MPI communicator
	comm = MPI.COMM_WORLD

	# parse input arguments
	parser = argparse.ArgumentParser(description='Mesh refinement test case for a fully-clamped flat plate')
	parser.add_argument('--hinit', default=0.25, type=float, help='target mesh size used for initial mesh generation')
	parser.add_argument('--order', default=3, type=int, help='order of elements in the mesh')
	parser.add_argument('--p', default=-1.e2, type=float, help='uniform pressure load to apply')
	parser.add_argument('--output', default='ks_fail', type=str, help='name of output of interest')
	parser.add_argument('--ksweight', default=1.e6, type=float, help='weight factor (rho) used for KS-aggregated outputs')
	parser.add_argument('--niters', default=0, type=int, help='number of refinement iterations to perform')
	parser.add_argument('--adapt', default=False, action='store_true', help='boolean flag to apply adaptive/uniform mesh refinement')
	parser.add_argument('--err_tol', default=1.e-6, type=float, help='target output error criterion used for adaptation termination')
	parser.add_argument('--nmaxlevs', default=30, type=int, help='max number of refinement levels allowed during adaptation')
	parser.add_argument('--ndecrease', default=5, type=int, help='number of iterations needed to decrease the refinement threshold to 1')
	args = parser.parse_args()
	for key in vars(args):
		printroot(f'args.{key} = {str(getattr(args,key))}')

	# remove old files
	if comm.rank == 0:
		for item in os.listdir(os.getcwd()):
			fname, fext = os.path.splitext(item)
			if fext in ['.f5','.plt']:
				os.remove(os.path.join(os.getcwd(), item))

	# create the adaptive model
	model = pyTACSAdapt(comm, initCallBack, elemCallBack, probCallBack, args.output.lower(), error_tol=args.err_tol, num_max_ref_levels=args.nmaxlevs, num_decrease_iters=args.ndecrease)
	model.setOutputOfInterest(args.output.lower())

	# initialize the coarse-space model
	model.initializeCoarseSpace(geom_file='plate.step', geom_type='quad', cmd_ln_args=args)

	# loop over adaptive iterations
	for k in range(args.niters+1):
		printroot(f'starting iteration {k}')

		# set up the coarse-space model
		model.setupModel(model_type='coarse', cmd_ln_args=args)
		model.coarse.setAuxElements(geom_labels=['panel'], aux_type='pressure', faceIndex=0, p=args.p)
		printroot('coarse model set up')

		# solve the coarse primal problem
		model.solvePrimal(model_type='coarse', writeSolution=True)
		printroot('coarse primal problem solved')

		if args.adapt:
			# solve the coarse adjoint problem
			model.solveAdjoint(model_type='coarse', writeSolution=False)
			printroot('coarse adjoint problem solved')

			# set up the fine-space model
			model.createFineSpace(cmd_ln_args=args)
			model.fine.setAuxElements(geom_labels=['panel'], aux_type='pressure', faceIndex=0, p=args.p)
			printroot('fine model set up')

			# get the interpolated state in the fine space
			model.interpolateField(field_type='state', writeSolution=False)
			printroot('fine state interpolated')

			# do the high-order reconstruction of the adjoint (computes the 
			# difference between the reconstructed and interpolated fine-space adjoints)
			model.reconstructField(field_type='adjoint', compute_diff=True, writeSolution=False)
			printroot('fine adjoint reconstructed')

			# do the adjoint-based element-wise error estimation
			error_estimate, output_correction, element_errors = model.estimateOutputError()
			printroot('output error estimated')

			if error_estimate <= args.err_tol:
				# end early if error tolerance is satisfied
				printroot(f'target error tolerance met')
				break
			elif k < (args.niters):
				# do the adaptation if not yet satisfied
				printroot('applying adaptive refinement')
				model.adaptModelDecreasingThreshold(element_errors)

		# otherwise do uniform mesh refinement
		else:
			if k < (args.niters):
				printroot('applying uniform refinement')
				model.coarse.applyRefinement()

	# print out final information
	printroot(f"mesh_history: {model.mesh_history}")
	printroot(f"output_history: {model.output_history}")
	printroot(f"error_history: {model.error_history}")
	printroot(f"adaptation_history: {model.adaptation_history}")