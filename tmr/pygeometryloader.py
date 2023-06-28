#!/usr/bin/python
"""
pygeometryloader - A python-based geometry loader for setting up TACS and TMR

This python interface is designed to provide an easier interface to the
C-layer of TACS and TMR. This module uses the TMR library for its reading of 
geometry models, automated meshing, and translation to TACS and TMR objects
"""
# =============================================================================
# Imports
# =============================================================================
import os
import numpy as np
from mpi4py import MPI
from tmr import TMR
from tacs.utilities import BaseUI
from functools import wraps


# =============================================================================
# Class Method Decorators
# =============================================================================
def model_method(method):
    @wraps(method)
    def wrapped_method(self, *args, **kwargs):
        assert (
            self.geom_model is not None
        ), f"`{method.__name__}` is a TMR.Model-based method. It may only be used after `readGeometryFile()` has been called"
        return method(self, *args, **kwargs)

    return wrapped_method


def mesh_method(method):
    @wraps(method)
    def wrapped_method(self, *args, **kwargs):
        assert (
            self.mesh is not None
        ), f"`{method.__name__}` is a TMR.Mesh-based method. It may only be used after `createMesh()` has been called"
        return method(self, *args, **kwargs)

    return wrapped_method


def topology_method(method):
    @wraps(method)
    def wrapped_method(self, *args, **kwargs):
        assert (
            self.topo is not None
        ), f"`{method.__name__}` is a TMR.Topology-based method. It may only be used after `createTopology()` has been called"
        return method(self, *args, **kwargs)

    return wrapped_method


def forest_method(method):
    @wraps(method)
    def wrapped_method(self, *args, **kwargs):
        assert (
            self.forest is not None
        ), f"`{method.__name__}` is a TMR.Quad/OctForest-based method. It may only be used after `createForest()` has been called"
        return method(self, *args, **kwargs)

    return wrapped_method


def assembler_method(method):
    @wraps(method)
    def wrapped_method(self, *args, **kwargs):
        assert (
            self.assembler is not None
        ), f"`{method.__name__}` is a TACS.Assembler-based method. It may only be used after `createTACSAssembler()` has been called"
        return method(self, *args, **kwargs)

    return wrapped_method


# =============================================================================
# Auxiliary Classes
# =============================================================================
class QuadCreator(TMR.QuadCreator):
    """
    Creator class for quadtrees. The user should not use this class directly
    """

    def __init__(self, bcs, dvpn):
        TMR.QuadCreator.__init__(bcs, dvpn)
        self.elemCallBack = None

    def setElemCallBack(self, _elemCallBack):
        self.elemCallBack = _elemCallBack

    def createElement(self, order, quad):
        return self.elemCallBack(order, quad)


class OctCreator(TMR.OctCreator):
    """
    Creator class for octrees. The user should not use this class directly
    """

    def __init__(self, bcs, dvpn):
        TMR.OctCreator.__init__(bcs, dvpn)
        self.elemCallBack = None

    def setElemCallBack(self, _elemCallBack):
        self.elemCallBack = _elemCallBack

    def createElement(self, order, oct):
        return self.elemCallBack(order, oct)


# =============================================================================
# pyGeometryLoader
# =============================================================================
class pyGeometryLoader(BaseUI):
    def __init__(self, comm, geom_type, printDebug=False):
        """
        Parameters
        ----------
        comm : mpi4py.MPI.Intracomm
            The comm object on which to create the pyGeometryLoader object and
            the underlying TMR.Quad/OctForest and TACS.Assembler objects

        geom_type : str in ['quad','oct']
            The type of geometry/mesh that will be assumed by TMR

        printDebug : bool
            Boolean flag to turn on/off print statements for debugging
        """
        allowed_geom_types = ("quad", "oct")
        assert (
            geom_type in allowed_geom_types
        ), f"geom_type must be one of {allowed_geom_types}"

        # Set MPI communicator
        BaseUI.__init__(self, comm=comm)

        # Debug printing flag
        self.printDebug = printDebug

        # initialize class attributes
        self.geom_type = geom_type
        self.geom_model = None
        self.model_name = ""
        self.bcs = None
        self.mesh_options = None
        self.fs = None
        self.mesh = None
        self.mesh_model = None
        self.topo = None
        self.forest = None
        self.creator = None
        self.assembler = None

        # initialize data structures for components/labels
        self.nComp = 0
        self.compDescripts = []
        self.compIDtoLabel = {}
        return

    def readGeometryFile(self, geom_file, **kwargs):
        """
        Read in a STEP/IGES/EGADS file and make the TMR geometry model

        Parameters
        ----------
        geom_file : str
            Filename of the geometry file to load (STEP/IGES/EGADS)

        kwargs : keyword args
            arg=value pairs that will directly be passed to TMR.LoadModel

            See TMR.LoadModel for more details
        """
        # get the name of the geometry file thats being loaded in
        self.model_name = os.path.splitext(geom_file)[0]

        # load in the geometry model
        self.geom_model = TMR.LoadModel(fname=geom_file, **kwargs)
        if self.geom_type == "quad":
            # ensure no solids are in model for quad type
            verts = self.geom_model.getVertices()
            edges = self.geom_model.getEdges()
            faces = self.geom_model.getFaces()
            self.geom_model = TMR.Model(verts, edges, faces)
        return

    @model_method
    def nameGeometricEntities(
        self,
        vol_names=None,
        face_names=None,
        edge_names=None,
        vert_names=None,
        writeToTecplot=False,
        **kwargs,
    ):
        """
        Specify labels for geometric entities in the model (volumes, faces,
        edges, and vertices). These labels can be used to specify boundary
        conditions and apply loads.

        NOTE: Geometric labels are equivalent to component descriptions, and are
        used to group underlying mesh elements into components. Due to the quad/
        oct geometry discretization in TMR, only volume (oct) and face (quad)
        labels are valid for component classification. Edge and Vertex labels
        are not valid component identifiers.

        Parameters
        ----------
        vol_names : dict
            Dictionary of index:'name' to assign to specified volumes

        face_names : dict
            Dictionary of index:'name' to assign to specified faces

        edge_names : dict
            Dictionary of index:'name' to assign to specified edges

        vert_names : dict
            Dictionary of index:'name' to assign to specified vertices

        writeToTecplot : bool
            Boolean flag to write out the geometry and its labels to a Tecplot-
            formatted .dat ASCII file. Helpful when checking if the labels were
            set correctly.

        kwargs : keyword arguments
            arg=value pairs that will be passed to writeModelToTecplot()
        """
        # set the names for volumes
        if vol_names is not None:
            vols = self.geom_model.getVolumes()
            for ind, name in vol_names.items():
                vols[ind].setName(name)

        # set the names for faces
        if face_names is not None:
            faces = self.geom_model.getFaces()
            for ind, name in face_names.items():
                faces[ind].setName(name)

        # set the names for edges
        if edge_names is not None:
            edges = self.geom_model.getEdges()
            for ind, name in edge_names.items():
                edges[ind].setName(name)

        # set the names for vertices
        if vert_names is not None:
            verts = self.geom_model.getVertices()
            for ind, name in vert_names.items():
                verts[ind].setName(name)

        # update the component information
        self._updateComps()

        # write out to Tecplot if requested
        if writeToTecplot and self.comm.rank == 0:
            tec_name = self.model_name + ".dat"
            self.geom_model.writeModelToTecplot(fname=tec_name, **kwargs)
        return

    @model_method
    def setBoundaryConditions(self, bc_info):
        """
        Specify boundary conditions for the named/labeled geometric entities

        Parameters
        ----------
        bc_info : dict
            Dictionary of the boundary condition information:
            'type' : str in ['volume', 'face', 'edge', 'vertex']
            'bc_names' : list[str], names of the labels to apply the BC to
            'bc_nums' : dict, dictionary of 'label':list[int] which specify the
                        DOF ids that will be constrined for each label
            'bc_vals' : dict, dictionary of 'label':list[float] which specify
                        the values to specify for each DOF of the BC

            example:
            bc_info = {'type' : 'edge',
                       'bc_names' : ['root','tip'],
                       'bc_nums' : {'tip':[0,1,2]},
                       'bc_vals' : {}
                       }
            Applies boundary conditions to edges, with the labels 'root' and
            'tip'. Only DOFs 0,1,2 are constrained on 'tip', whereas all DOFs
            (default) are constrained on 'root'. The specified values assigned
            to each DOF can be set with bc_vals, but here the default (0.0) is
            used.
        """
        valid_bc_types = ("volume", "face", "edge", "vertex")

        # create the BoundaryConditions object if it does not exist yet
        if self.bcs is None:
            self.bcs = TMR.BoundaryConditions()

        # get the names of the geometric entities the BC will can be applied to
        entity_names = []
        bc_type = bc_info.get("type", None)
        assert (
            bc_type in valid_bc_types
        ), f"the type of boundary condition must be one of {valid_bc_types}"
        if bc_type == "volume":
            geom_entities = self.geom_model.getVolumes()
        elif bc_type == "face":
            geom_entities = self.geom_model.getFaces()
        elif bc_type == "edge":
            geom_entities = self.geom_model.getEdges()
        elif bc_type == "vertex":
            geom_entities = self.geom_model.getVertices()
        for entity in geom_entities:
            label = entity.getName()
            if label is not None:
                entity_names.append(label)

        # add the specified boundary conditions based on their name/label
        bc_names = bc_info.get("bc_names", [])
        bc_nums = bc_info.get("bc_nums", {})
        bc_vals = bc_info.get("bc_vals", {})
        for name in bc_names:
            if name in entity_names:
                nums = bc_nums.get(name, None)
                vals = bc_vals.get(name, None)
                self.bcs.addBoundaryCondition(name, nums, vals)
            else:
                # warn user their supplied label isn't in the model
                self._user_warning(
                    self.comm.rank,
                    f'skipping boundary condition definition for "{name}", label not found in model',
                )
        return

    def setMeshOptions(self, **kwargs):
        """
        Specify options that will be used by the TMR mesh generator

        Parameters
        ----------
        kwargs : keyword arguments
            arg=value pairs that will be set directly into the TMR.MeshOptions
            object attributes

            See TMR.MeshOptions for details
        """
        allowed_keys = (
            "mesh_type_default",
            "num_smoothing_steps",
            "frontal_quality_factor",
            "triangularize_print_level",
            "triangularize_print_iter",
            "reset_mesh_objects",
            "write_mesh_quality_histogram",
            "write_init_domain_triangle",
            "write_triangularize_intermediate",
            "write_pre_smooth_triangle",
            "write_post_smooth_triangle",
            "write_dual_recombine",
            "write_pre_smooth_quad",
            "write_post_smooth_quad",
            "write_quad_dual",
        )

        # create the mesh_options object if it does not yet exist
        if self.mesh_options is None:
            self.mesh_options = TMR.MeshOptions()

        # update the mesh_options attributes
        for key, val in kwargs.items():
            if key in allowed_keys:
                setattr(self.mesh_options, key, val)
            else:
                self._user_warning(
                    self.comm.rank, f'skipping "{key}", not a valid mesh option'
                )
        return

    def setElementFeatureSize(self, feature_type, **kwargs):
        """
        Set element sizing parameters used during meshing, assuming a specific
        kind of feature

        Parameters
        ----------
        feature_type : str
            String denoting the type of ElementFeatureSize to use

        kwargs : keyword args
            arg=value pairs that will directly be set into the constructor of
            the selected ElementFeatureSize

            See TMR.ElementFeatureSize for details on keyword arguments
        """
        allowed_feature_types = ("const", "linear", "box", "point")
        assert (
            feature_type in allowed_feature_types
        ), f"feature_type must be one of {allowed_feature_types}"

        # create the feature size object for the specified type
        if feature_type == "const":
            self.fs = TMR.ConstElementSize(**kwargs)
        elif feature_type == "linear":
            self.fs = TMR.LinearElementSize(**kwargs)
        elif feature_type == "box":
            self.fs = TMR.BoxElementSize(**kwargs)
        elif feature_type == "point":
            self.fs = TMR.PointElementSize(**kwargs)
        return

    @model_method
    def createMesh(self, h=1.0, writeBDF=False):
        """
        Mesh the geometry with specified mesh parameters

        Parameters
        ----------
        h : float
            Global element size target used during meshing
            (Only used if no ElementFeatureSize has been defined)

        writeBDF : bool
            Boolean flag to write out the mesh to a .bdf format
        """
        # create the mesh object if it does not yet exist
        if self.mesh is None:
            self.mesh = TMR.Mesh(self.comm, self.geom_model)

        # create the mesh
        self.mesh.mesh(h=h, opts=self.mesh_options, fs=self.fs)

        # create the model derived from the mesh
        self.mesh_model = self.mesh.createModelFromMesh()

        # write out the BDF if specified
        if writeBDF:
            bdf_name = self.model_name + ".bdf"
            self.mesh.writeToBDF(bdf_name, outtype=self.geom_type, bcs=self.bcs)
        return

    @model_method
    def createTopology(self, useMesh=False):
        """
        Set the topology of the model based on either the initial geometry or
        the mesh

        Parameters
        ----------
        useMesh : bool
            Boolean flag used to determine which model to use for the topology
        """
        if useMesh:
            assert (
                self.mesh_model is not None
            ), "Mesh must be created prior to setting the topology with the mesh model"
            self.topo = TMR.Topology(self.comm, self.mesh_model)
        else:
            self.topo = TMR.Topology(self.comm, self.geom_model)
        return

    @topology_method
    def createForest(self, nlevels=0, **kwargs):
        """
        Create the appropriate forest of quad/octree semi-structured
        partitioning of the full mesh

        Parameters
        ----------
        nlevels : int
            Number of levels of trees to make for the forest. Note that for each
            level the size of the mesh grows as 4x for quads and 8x for octs

        kwargs : keyword args
            arg=value pairs that will directly be set into the appropriate
            constructor for the forest

            See TMR.QuadForest/TMR.OctForest for details on keyword arguments
        """
        # create the forest for the specified type of geometry
        if self.geom_type == "quad":
            # create the forest of quadtrees
            self.forest = TMR.QuadForest(comm=self.comm, **kwargs)
        elif self.geom_type == "oct":
            # create the forest of octrees
            self.forest = TMR.OctForest(comm=self.comm, **kwargs)

        # set the topology for the forest
        self.forest.setTopology(self.topo)

        # create the trees for the forest
        self.forest.createTrees(nlevels)
        return

    @forest_method
    def createTACSAssembler(self, elemCallBack, dvpn=1, **kwargs):
        """
        Create and return the TACS.Assembler object used for analysis and
        optimization

        Parameters
        ----------
        elemCallBack : function handle
            Handle for function used to create the elements for the TACS Assembler
            See TMR.QuadCreator and TMR.OctCreator for details

        dvpn : int
            Number of design variables per node

        kwargs : keyword args
            arg=value pairs that will directly be set into
            Quad/OctCreator.createTACS()

        Returns
        -------
        assembler : TACS.Assembler object
        """
        # create the creator if it doesn't exist
        if self.creator is None:
            if self.geom_type == "quad":
                self.creator = QuadCreator(self.bcs, dvpn)
            elif self.geom_type == "oct":
                self.creator = OctCreator(self.bcs, dvpn)

            self.creator.setElemCallBack(elemCallBack)

        # create the TACSAssembler
        self.assembler = self.creator.createTACS(
            forest=self.forest, comp_names=self.compDescripts, **kwargs
        )
        return self.assembler

    def getModel(self, returnMeshModel=False):
        """
        Return the TMR.Model object

        Parameters
        ----------
        returnMeshModel : bool
            Boolean flag to select which model to return
        """
        if not returnMeshModel:
            return self.geom_model
        else:
            return self.mesh_model

    def getForest(self):
        """
        Return the TMR.Quad/OctForest object
        """
        return self.forest

    def getAssembler(self):
        """
        Return the TACS.Assembler object
        """
        return self.assembler

    @assembler_method
    def getNumOwnedNodes(self):
        """
        Return number of nodes owned by this processor.

        Returns
        -------
        numOwnedNodes : int
            The number of nodes owned locally by this proc
        """
        return self.assembler.getNumOwnedNodes()

    def getNumTotalNodes(self):
        """
        Return the total number of nodes in the mesh

        Returns
        -------
        numTotNodes : int
            The number of nodes in the entire mesh
        """
        # get the number of nodes owned locally on this proc
        numLocNodes = self.getNumOwnedNodes()

        # sum all proc contributions
        numTotNodes = self.comm.allreduce(numLocNodes, op=MPI.SUM)
        return numTotNodes

    @assembler_method
    def getNumOwnedElements(self):
        """
        Return number of elements owned by this processor.

        Returns
        -------
        numOwnedElems : int
            The number of elements owned locally by this proc
        """
        return self.assembler.getNumElements()

    def getNumTotalElements(self):
        """
        Return the total number of elements in the mesh

        Returns
        -------
        numTotElems : int
            The number of elements in the entire mesh
        """
        # get the number of nodes owned locally on this proc
        numLocElems = self.getNumOwnedElements()

        # sum all proc contributions
        numTotElems = self.comm.allreduce(numLocElems, op=MPI.SUM)
        return numTotElems

    @forest_method
    def getLocalNodes(self, nodeIDs=None):
        """
        Return x,y,z location for the nodes owned locally on this proc

        NOTE: The TMR.Quad/OctForest object contains the partitioned mesh

        NOTE: The mesh partition defined in a proc's forest and assembler may
        be slightly different. The forest mesh may include nodes on the mesh
        partition boundary that are not "owned" in the assembler.

        Parameters
        ----------
        nodeIDs : None or int or list[int]
            Local node IDs to return their node coordinates

        Returns
        -------
        coords : np.ndarray (2D)
            Coordinates of the local nodes
        """
        coords = self.forest.getPoints()
        if nodeIDs is not None:
            coords = coords[nodeIDs]
        return coords

    @mesh_method
    def getGlobalNodes(self, nodeIDs=None):
        """
        Return x,y,z location of for nodes in the global mesh

        NOTE: The TMR.Mesh object contains the entire global mesh, not the local
        partitioned mesh owned by a proc

        Parameters
        ----------
        nodeIDs : None or int or list[int]
            Global node IDs to return their node coordinates

        Returns
        -------
        coords : np.ndarray (2D)
            Coordinates of the global nodes
        """
        coords = self.mesh.getMeshPoints()
        if nodeIDs is not None:
            coords = coords[nodeIDs]
        return coords

    @forest_method
    def getNodeRange(self):
        """
        Get the range of global node numbers each processor owns

        Returns
        -------
        nodeRange : list[int]
            Node ranges defined for each proc
        """
        return self.forest.getNodeRange()

    @forest_method
    def getGlobalToLocalNodeIDDict(self):
        """
        Return a mapping of global to local node numbering for the nodes owned
        by this proc

        NOTE: This local numbering is for the nodes "owned" by this proc's
        forest and not its assembler. Therefore, the local IDs may not be
        directly compatible with the mesh stored in the assembler

        Returns
        -------
        globalToLocalNodeMap : dict
            Dictionary with globalNodeID:localNodeID key:value pairs for the nodes
            that are owned by this proc
        """
        # get the global node numbers owned by this proc
        nodeRange = self.getNodeRange()
        startID = nodeRange[self.comm.rank]
        endID = nodeRange[self.comm.rank + 1]
        global_ids = range(startID, endID)

        # create the global to local node map
        globalToLocalNodeMap = {
            gID: self.forest.getLocalNodeNumber(gID) for gID in global_ids
        }
        return globalToLocalNodeMap

    def getLocalNodeIDsFromGlobal(self, globalIDs):
        """
        Given a list of node IDs in global (non-partitioned) ordering returns
        the local (partitioned) node IDs for each processor. If a requested node
        is not included on this processor, an entry of -1 will be returned.

        Parameters
        ----------
        globalIDs : list[int]
            Global node IDs

        Returns
        -------
        localIDs : list[int]
            Local (partitioned) node ids
        """
        # get the global to local node id map for this proc
        globalToLocalNodeMap = self.getGlobalToLocalNodeIDDict()

        # get the local ids
        localIDs = [globalToLocalNodeMap.get(gID, -1) for gID in globalIDs]
        return localIDs

    def getGlobalToLocalElemIDDict(self):
        """
        Return a mapping of global to local element numbering for the elements
        owned by this proc

        Returns
        -------
        globalToLocalElemMap : dict
            Dictionary with globalElemID:localElemID key:value pairs for the
            elems that are owned by this proc
        """
        # get the number of elements owned on this proc
        numLocElems = self.getNumOwnedElements()

        # get the number of elems owned by each proc
        counts = self.comm.allgather(numLocElems)

        # compute this procs global elem offset
        offset = sum(counts[: self.comm.rank])

        # create the global to local elem map using the offset
        globalToLocalElemMap = {
            (localID + offset): localID for localID in range(numLocElems)
        }
        return globalToLocalElemMap

    def getLocalElementIDsFromGlobal(self, globalIDs):
        """
        Given a list of element IDs in global (non-partitioned) ordering returns
        the local (partitioned) elem IDs for each processor. If a requested elem
        is not included on this processor, an entry of -1 will be returned.

        Parameters
        ----------
        globalIDs : list[int]
            Global element IDs

        Returns
        -------
        localIDs : list[int]
            Local (partitioned) element IDs
        """
        # get the global to local elem id map for this proc
        globalToLocalElemMap = self.getGlobalToLocalElemIDDict()

        # get the local ids
        localIDs = [globalToLocalElemMap.get(gID, -1) for gID in globalIDs]
        return localIDs

    @forest_method
    def getLocalDepNodeConn(self):
        """
        Return the dependent node connectivity information for the locally owned
        partitioned mesh

        NOTE: Uses global node numbers

        Returns
        -------
        ptr : np.ndarray (1D)
            dependent node ids
        conn : np.ndarray (2D)
            connectivity of the dependent nodes
        weights : np.ndarray (1D)
            weights associated with the dependent nodes
        """
        return self.forest.getDepNodeConn()

    @forest_method
    def getLocalElementConnectivity(self):
        """
        Return the connectivity for the elements owned locally on this proc

        NOTE: Uses global node numbers

        NOTE: The TMR.Quad/OctForest object contains the partitioned mesh

        Returns
        -------
        conn : np.ndarray (2D)
            Connectivity of the elements local to this proc
        """
        return self.forest.getMeshConn()

    @mesh_method
    def getGlobalElementConnectivity(self):
        """
        Return the element connectivity for the entire global mesh

        NOTE: The TMR.Mesh object contains the entire global mesh, not the local
        partitioned mesh owned by a proc

        Returns
        -------
        conn : np.ndarray (2D)
            Element connectivity of the enitre global mesh
        """
        conn = None
        if self.geom_type == "quad":
            conn = self.mesh.getQuadConnectivity()
        elif self.geom_type == "oct":
            conn = self.mesh.getHexConnectivity()
        return conn

    def getNumComponents(self):
        """
        Return the number of user-defined components currently in the model

        Returns
        -------
        nComp : int
            the number of components defined in the model
        """
        return self.nComp

    def getComponentDescripts(self):
        """
        Return user-defined labels for each component

        Returns
        -------
        compDescripts : list[str]
            List of component descriptions
        """
        return self.compDescripts

    def getElementComponents(self):
        """
        Return the list of component IDs that each element in the mesh belongs to

        Returns
        -------
        compIDList : list[int]
            List of component IDs associated with each element in the global
            mesh. A component ID of -1 indicates an element is not associated
            with a component
        """
        nelems = self.getNumOwnedElements()
        compIDList = -1 * np.ones(nelems, dtype=int)
        for compID in self.compIDtoLabel.keys():
            # get the element IDs for each compID on this proc
            elemIDs = self.getLocalElementIDsForComps(compID)
            # set the comp IDs for these elements
            compIDList[elemIDs] = compID

        # assemble all compIDs from every proc
        compIDList = list(compIDList)
        compIDList = self.comm.allgather(compIDList)
        compIDList = [cID for sublist in compIDList for cID in sublist]

        # check that all elements belong to a component
        comps_not_found = [ind for ind, val in enumerate(compIDList) if val == -1]
        if comps_not_found:
            # warn user if some elements aren't in a any component
            self._user_warning(
                self.comm.rank,
                f"the following {len(comps_not_found)} elements were not associated with any component defined in the model: {comps_not_found}",
            )
        return compIDList

    def getConnectivityForComp(self, componentIDs):
        """
        Returns element connectivities associated with the specified components

        NOTE: Uses global node numbers

        Parameters
        ----------
        componentIDs : int or list[int]
            List of component ID numbers

        Returns
        -------
        conn : np.ndarray
            Element connectivity for the specified component IDs
        """
        # get the global element ids associated with the components
        elemIDs = self.getGlobalElementIDsForComps(componentIDs)

        # get the full global element connectivity
        global_conn = self.getGlobalElementConnectivity()

        # only retain the connectivity for associated elems
        conn = global_conn[elemIDs]
        return conn

    @forest_method
    def getLocalNodeIDsForComps(self, componentIDs):
        """
        Return the local (partitioned) node IDs on this proc that belong to the
        specified components

        Parameters
        ----------
        componentIDs : int or list[int]
            List of component ID numbers

        Returns
        -------
        nodeIDs : list[int]
            Local node IDs owned by this proc that belong to the specified
            components
        """
        if type(componentIDs) == int:
            componentIDs = [componentIDs]

        nodeIDs = set()
        for compID in componentIDs:
            # get the current component label
            comp_name = self.compIDtoLabel[compID]

            # get the local node ids associated with this label
            nodeIDs.update(self.forest.getNodesWithName(comp_name))

        return list(nodeIDs)

    def getGlobalNodeIDsForComps(self, componentIDs):
        """
        Returns a list of global node IDs associated with the specified
        components

        Parameters
        ----------
        componentIDs : int or list[int]
            List of component ID numbers

        Returns
        -------
        nodeIDs : list[int]
            Unique node IDs that belong to the given component IDs
        """
        # get the total element connectivity (global node ids) for the specified compIDs
        conn = self.getConnectivityForComp(componentIDs)

        # flatten and uniquify the nodeIDs from the connectivity info
        nodeIDs = list(set(np.ravel(conn)))
        return nodeIDs

    @forest_method
    def getLocalElementIDsForComps(self, componentIDs):
        """
        Return the local element IDs on this proc that belong to the specified
        components

        Parameters
        ----------
        componentIDs : int or list[int]
            List of component ID numbers

        Returns
        -------
        elemIDs : list[int]
            Local element IDs that belong to the given component IDs
        """
        if type(componentIDs) == int:
            componentIDs = [componentIDs]

        elemIDs = []
        for compID in componentIDs:
            # get the current component label
            comp_name = self.compIDtoLabel[compID]

            # get the geometric objects associated with this label
            geom_objs = None
            if self.geom_type == "quad":
                geom_objs = self.forest.getQuadsWithName(comp_name)
            elif self.geom_type == "oct":
                geom_objs = self.forest.getOctsWithName(comp_name)

            # get the element IDs associated with these objects
            for obj in geom_objs:
                # the tag attribute holds the LOCAL element id for this proc
                elemIDs.append(obj.tag)

        return elemIDs

    def getGlobalElementIDsForComps(self, componentIDs):
        """
        Returns a list of global element IDs associated with the specified
        components

        Parameters
        ----------
        componentIDs : int or list[int]
            List of component ID numbers

        Returns
        -------
        elemIDs : list[int]
            Global element IDs that belong to the given component IDs
        """
        # get the local elem ids for this proc
        localElemIDs = self.getLocalElementIDsForComps(componentIDs)

        # convert the local elem ids to global elem ids
        globalToLocalElemIDDict = self.getGlobalToLocalElemIDDict()
        localToGlobalElemIDDict = {
            val: key for key, val in globalToLocalElemIDDict.items()
        }
        globalElemIDs = [localToGlobalElemIDDict[locID] for locID in localElemIDs]

        # share global elem ids across all procs
        elemIDs = self.comm.allgather(globalElemIDs)
        elemIDs = [eid for sublist in elemIDs for eid in sublist]
        return elemIDs

    def _updateComps(self):
        """
        Add and/or update component IDs and descriptions based on the labels
        assigned to each geometric entity in the geometry model

        NOTE: The only valid component descriptions for Quad/Oct models are the
        labels associated with faces/volumes respectively. Edge and Vertex
        labels are never valid component descriptions. Since TMR assumes either
        Quad/Oct geometry exclusively, line and point elements are not supported.
        Likewise, TMR cannot mix Quad-based and Oct-based discretizations.
        """
        self.nComp = 0
        self.compDescripts = []
        comps_not_set = []

        # check volumes
        vols = self.geom_model.getVolumes()
        for vol in vols:
            label = vol.getName()
            if label is not None:
                if self.geom_type == "oct":
                    if (
                        self.compIDtoLabel.get(self.nComp)
                        and label != self.compIDtoLabel[self.nComp]
                    ):
                        self._user_warning(
                            self.comm.rank,
                            f'component ID {self.nComp}: description updated "{self.compIDtoLabel[self.nComp]}" ---> "{label}"',
                        )
                    self.compIDtoLabel[self.nComp] = label
                    self.compDescripts.append(label)
                    self.nComp += 1
                else:
                    comps_not_set.append(label)

        # check faces
        faces = self.geom_model.getFaces()
        for face in faces:
            label = face.getName()
            if label is not None:
                if self.geom_type == "quad":
                    if (
                        self.compIDtoLabel.get(self.nComp)
                        and label != self.compIDtoLabel[self.nComp]
                    ):
                        self._user_warning(
                            self.comm.rank,
                            f'component ID {self.nComp}: description updated "{self.compIDtoLabel[self.nComp]}" ---> "{label}"',
                        )
                    self.compIDtoLabel[self.nComp] = label
                    self.compDescripts.append(label)
                    self.nComp += 1
                else:
                    comps_not_set.append(label)

        # check edges
        edges = self.geom_model.getEdges()
        for edge in edges:
            label = edge.getName()
            if label is not None:
                comps_not_set.append(label)

        # check vertices
        verts = self.geom_model.getVertices()
        for vert in verts:
            label = vert.getName()
            if label is not None:
                comps_not_set.append(label)

        # warn user about labels that are not valid component descriptions
        if comps_not_set:
            self._user_warning(
                self.comm.rank,
                f"the following geometry labels are not valid component descriptions:\n\t{comps_not_set}",
            )
        return

    def _user_warning(self, rank, msg, rootOnly=True):
        """
        Writes a warning to the screen on the specified proc. The user should
        not use this function directly.

        Parameters
        ----------
        msg : str
            message to write to screen
        """
        if rootOnly and rank == 0:
            print(f"[{rank}] pyGeometryLoader Warning:\n\t{msg}")
        elif not rootOnly:
            print(f"[{rank}] pyGeometryLoader Warning:\n\t{msg}")
        return
