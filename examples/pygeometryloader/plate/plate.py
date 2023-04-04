"""
A simple test of the functionality of the pyGeometryLoader class.

To run the defaults: mpirun -np 2 python plate.py
"""
# imports
import os
from mpi4py import MPI
from tmr import TMR
from tmr.pygeometryloader import pyGeometryLoader
from tacs import elements, constitutive
from make_plate_geom import makePlateGeom


# helper function for parallel printing
def pprint(msg):
    print(f"[{comm.rank}] {msg}")


# element callback function used to create TACS
def elemCallBack(order, quad):
    # set constitutive properties
    rho = 2700.0
    E = 70.0e9
    nu = 0.33
    ys = 276.0e6
    min_thickness = 0.001
    max_thickness = 0.01
    thickness = 0.0025
    props = constitutive.MaterialProperties(rho=rho, E=E, nu=nu, ys=ys)
    stiff = constitutive.IsoShellConstitutive(
        props, t=thickness, tlb=min_thickness, tub=max_thickness
    )

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
        raise ValueError(f"orders {2}-{4} supported")
    return element


# get the MPI communicator
comm = MPI.COMM_WORLD

# remove old files
if comm.rank == 0:
    for item in os.listdir(os.getcwd()):
        fname, fext = os.path.splitext(item)
        if fext in [".step", ".iges", ".dat", ".bdf"]:
            os.remove(os.path.join(os.getcwd(), item))

# create the geometry model
fname = "plate"
makeIGES = False
if makeIGES:
    fname += ".iges"
else:
    fname += ".step"
if comm.rank == 0:
    makePlateGeom(makeIGES=makeIGES)
comm.barrier()

# load the geometry model
geomLoader = pyGeometryLoader(comm, "quad")
geomLoader.readGeometryFile(fname)

# name the face and edges
face_names = {0: "panel/0", 1: "panel/1"}
edge_names = {3: "edges/root"}
geomLoader.nameGeometricEntities(
    face_names=face_names, edge_names=edge_names, writeToTecplot=True
)

# specify boundary conditions - use defaults for nums and vals to fully fix
bc_names = ["edges/root", "bogus_name"]
bc_info = {"type": "edge", "bc_names": bc_names}
geomLoader.setBoundaryConditions(bc_info)

# set the meshing options
geomLoader.setMeshOptions(write_mesh_quality_histogram=0)

# create a coarse mesh
geomLoader.createMesh(h=0.5, writeBDF=True)

# create the topology
geomLoader.createTopology(useMesh=True)

# create the quadforest
geomLoader.createForest(order=2, interp=TMR.UNIFORM_POINTS)

# create the assembler
geomLoader.createTACSAssembler(elemCallBack)

# test functions
globalIDs = [2]
compIDs = [0]
pprint(f"getModel(): {geomLoader.getModel()}")
pprint(f"getForest(): {geomLoader.getForest()}")
pprint(f"getAssembler(): {geomLoader.getAssembler()}")
pprint(f"getNumOwnedNodes(): {geomLoader.getNumOwnedNodes()}")
pprint(f"getNumTotalNodes(): {geomLoader.getNumTotalNodes()}")
pprint(f"getNumOwnedElements(): {geomLoader.getNumOwnedElements()}")
pprint(f"getNumTotalElements(): {geomLoader.getNumTotalElements()}")
pprint(f"getLocalNodes(): {geomLoader.getLocalNodes()}")
pprint(f"getGlobalNodes(): {geomLoader.getGlobalNodes()}")
pprint(f"getNodeRange(): {geomLoader.getNodeRange()}")
pprint(f"getGlobalToLocalNodeIDDict(): {geomLoader.getGlobalToLocalNodeIDDict()}")
pprint(
    f"getLocalNodeIDsFromGlobal(): {geomLoader.getLocalNodeIDsFromGlobal(globalIDs)}"
)
pprint(f"getGlobalToLocalElemIDDict(): {geomLoader.getGlobalToLocalElemIDDict()}")
pprint(
    f"getLocalElementIDsFromGlobal(): {geomLoader.getLocalElementIDsFromGlobal(globalIDs)}"
)
pprint(f"getLocalDepNodeConn(): {geomLoader.getLocalDepNodeConn()}")
pprint(f"getLocalElementConnectivity(): {geomLoader.getLocalElementConnectivity()}")
pprint(f"getGlobalElementConnectivity(): {geomLoader.getGlobalElementConnectivity()}")
pprint(f"getNumComponents(): {geomLoader.getNumComponents()}")
pprint(f"getComponentDescripts(): {geomLoader.getComponentDescripts()}")
pprint(f"getElementComponents(): {geomLoader.getElementComponents()}")
pprint(f"getConnectivityForComp(): {geomLoader.getConnectivityForComp(compIDs)}")
pprint(f"getLocalNodeIDsForComps(): {geomLoader.getLocalNodeIDsForComps(compIDs)}")
pprint(f"getGlobalNodeIDsForComps(): {geomLoader.getGlobalNodeIDsForComps(compIDs)}")
pprint(
    f"getLocalElementIDsForComps(): {geomLoader.getLocalElementIDsForComps(compIDs)}"
)
pprint(
    f"getGlobalElementIDsForComps(): {geomLoader.getGlobalElementIDsForComps(compIDs)}"
)
