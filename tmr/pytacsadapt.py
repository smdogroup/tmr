#!/usr/bin/python
"""
pytacsadapt - A python-based interface for doing output-based error estimation 
and mesh adaptation with TACS and TMR.

This interface is designed to provide an easier-to-use interface to the C++ 
layer of TACS and TMR. 
"""

# =============================================================================
# Imports
# =============================================================================
import copy
from mpi4py import MPI
import numpy as np
from tmr import TMR
from tmr.pygeometryloader import *
from tacs.utilities import BaseUI


# =============================================================================
# pyApproximationSpace
# =============================================================================
class pyApproximationSpace(BaseUI):
    """
    This class defines an approximation space for the finite-element method that
    is being applied to the problem of interest. It consists of both TACS- and
    TMR-related objects which include: the model geometry, Quad/OctForest,
    finite-element assembler, finite-element problem, and copies of necessary
    information used in an adaptive solution.
    """

    def __init__(self, comm, _elemCallBack, _probCallBack):
        """
        Constructor

        Parameters
        ----------
        comm : mpi4py.MPI.Intracomm
            The comm object on which to create this approximation space

        _elemCallBack : function handle
            Handle to user-defined function for creating elements when making
            the TACS Assembler. See pyGeometryLoader.createTACSAssembler()
            for more details on the function prototype.

        _probCallBack : function handle
            Handle to user-defined function for creating the desired TACSProblem
            to be used in analysis. This function should create the TACSProblem
            and set an output to be used for error estimation and adaptation.
            It has the following prototype:
            _probCallBack(name, comm, geomLoader, assembler, **kwargs) -> TACSProblem
        """
        # Set MPI communicator
        BaseUI.__init__(self, comm=comm)

        # initialize the remaining class attributes
        self.elemCallBack = _elemCallBack
        self.probCallBack = _probCallBack
        self.geomLoader = None
        self.forest = None
        self.assembler = None
        self.problem = None
        self.state = None
        self.adjoint = None
        self.refine_iter = 0
        self.prob_name = ""
        self.output_name = ""
        return

    def createAssembler(self):
        """
        Create the TACS Assembler object and the state and adjoint vectors to be
        used during the adaptive analysis
        """
        self.assembler = self.geomLoader.createTACSAssembler(self.elemCallBack)
        self.state = self.assembler.createVec()
        self.adjoint = self.assembler.createVec()
        return

    def createProblem(self, prob_name, **kwargs):
        """
        Creates the TACS problem object to be used for analysis. To apply fixed
        loads to the problem, pyApproximationSpace.setAuxElements() can be called
        after pyApproximationSpace.createProblem() as many times as necessary.

        Parameters
        ----------
        prob_name : str
            Name to assign to the problem
        """
        # create the problem
        self.prob_name = prob_name
        self.problem = self.probCallBack(
            name=prob_name,
            comm=self.comm,
            geomLoader=self.geomLoader,
            assembler=self.assembler,
            **kwargs,
        )

        # check the output of interest is present in the problem
        output_names = self.problem.getFunctionKeys()
        assert (
            self.output_name in output_names
        ), f"Selected output of interest ({self.output_name}) is not in the list of outputs defined in the problem: {output_names}"
        return

    def setAuxElements(self, geom_labels, aux_type, **kwargs):
        """
        Creates and assigns the auxiliary elements to the assembler and problem.

        Applies identical aux elements to every element that belongs to the
        specified geometric entities.

        NOTE: Does not reset the aux elements before adding new ones. Previously
        added aux elements will persist unless manually reset.

        Parameters
        ----------
        geom_labels : list[str]
            Labels of geometric entities in the model with which the aux elements
            will be added to

        aux_type : str in ('traction', 'pressure', 'inertial', 'centrifugal')
            The type of aux element to add

        kwargs : keyword args
            arg=value pairs that will directly be passed to the aux element constructor
        """
        # check the aux_type is set properly
        valid_aux_types = ("traction", "pressure", "inertial", "centrifugal")
        aux_type = aux_type.lower()
        assert aux_type in valid_aux_types, f"aux_type must be one of {valid_aux_types}"

        for name in geom_labels:
            # get the elem ids associated with this geometric entity name
            entities = None
            if isinstance(self.forest, TMR.QuadForest):
                entities = self.forest.getQuadsWithName(name)
            elif isinstance(self.forest, TMR.OctForest):
                entities = self.forest.getOctsWithName(name)
            elem_ids = []
            for entity in entities:
                elem_ids.append(entity.tag)  # tag for elem id

            # create and add the aux elements
            for eid in elem_ids:
                # get the element object
                elem, _, _, _, _ = self.assembler.getElementData(eid)

                # create the specified aux element
                aux = None
                if aux_type == "traction":
                    aux = elem.createElementTraction(
                        kwargs.get("faceIndex"), kwargs.get("trac")
                    )
                elif aux_type == "pressure":
                    aux = elem.createElementPressure(
                        kwargs.get("faceIndex"), kwargs.get("p")
                    )
                elif aux_type == "inertial":
                    aux = elem.createElementInertialForce(kwargs.get("inertiaVec"))
                elif aux_type == "centrifugal":
                    aux = elem.createElementCentrifugalForce(
                        kwargs.get("omegaVec"),
                        kwargs.get("rotCenter"),
                        kwargs.get("firstOrder", False),
                    )

                # add to the aux elems in the problem
                self.problem.auxElems.addElement(eid, aux)

        # set the aux elems into the assembler too
        self.assembler.setAuxElements(self.problem.auxElems)
        return

    def applyRefinement(
        self, refine_indicator=None, num_min_levels=0, num_max_levels=TMR.MAX_LEVEL
    ):
        """
        Refine the mesh based on an element refinement indicator. Hanging-node
        h-refinement is carried out, followed by a 2:1 balancing such that there
        is only one level of difference between adjacent elements.

        Parameters
        ----------
        refine_indicator : None or np.ndarray (1D, int) [num_elements]
            Refinement indicator array, integer number of levels each element will
            be refined or coarsened. If None, uniform refinement of a single level
            is applied to the mesh.

        num_min_levels : int
            Minimum number of levels the forest is allowed to refine the elements

        num_max_levels : int
            Maximum number of levels the forest is allowed to refine the elements
        """
        self.forest.refine(
            refine_indicator, min_lev=num_min_levels, max_lev=num_max_levels
        )
        self.forest.balance(1)
        self.forest.repartition()
        self.forest.createNodes()
        self.refine_iter += 1
        return


# =============================================================================
# pyTACSAdapt
# =============================================================================
class pyTACSAdapt(BaseUI):
    """
    This class contains the coarse- and fine-space finite element approximations
    that can be used for output-based error estimation and adaptation. This is
    the high-level class that should be used for adaptive analysis.
    """

    def __init__(
        self, comm, _initCallBack, _elemCallBack, _probCallBack, _adapt_output, **kwargs
    ):
        """
        Constructor

        Parameters
        ----------
        comm : mpi4py.MPI.Intracomm
            The comm object on which to create the coarse- and fine-space models
            used in the adaptive analysis.

        _initCallBack : function handle
            Handle to user-defined function for initializing the adaptive model.
            This function should set up the coarse approximation space info by
            reading the geometry, making the initial mesh, and setting the
            initial Quad/OctForest.
            It has the following prototype:
            _initializeCallBack(coarse_space, **kwargs)

        _elemCallBack : function handle
            Handle to user-defined function for creating elements when making
            the TACS Assembler. See pyGeometryLoader.createTACSAssembler()
            for more details on the function prototype.

        _probCallBack : function handle
            Handle to user-defined function for creating the desired TACSProblem
            to be used in analysis. See pyApproximationSpace.__init__() for more
            details on the function prototype.

        _adapt_output : str
            Name of the output to be used for the adaptive analysis

        kwargs : keyword args
            arg=value pairs that can be used to specify adaptation parameters
        """
        # Set MPI communicator
        BaseUI.__init__(self, comm=comm)

        # initialize the class attributes
        self.initCallBack = _initCallBack
        self.coarse = pyApproximationSpace(comm, _elemCallBack, _probCallBack)
        self.fine = pyApproximationSpace(comm, _elemCallBack, _probCallBack)
        self.adapt_output = _adapt_output
        self.error_tol = kwargs.get("error_tol", None)
        self.num_min_ref_levels = kwargs.get("num_min_ref_levels", 0)
        self.num_max_ref_levels = kwargs.get("num_max_ref_levels", TMR.MAX_LEVEL)
        self.num_decrease_iters = kwargs.get("num_decrease_iters", None)
        self.mesh_history = {}
        self.output_history = {}
        self.error_history = {}
        self.adaptation_history = {}

        # set the output used for adaptation
        self.setOutputOfInterest(self.adapt_output)
        return

    def initializeCoarseSpace(self, **kwargs):
        """
        Initializes the coarse approximation space using the user-defined
        initCallBack() function.

        Parameters
        ----------
        kwargs : keyword args
            arg=value pairs that will directly be passed to self.initCallBack()
        """
        self.initCallBack(self.coarse, **kwargs)
        return

    def createFineSpace(self, **kwargs):
        """
        Creates the "fine" finite-element approximation space by copying the
        "coarse" approximation space and incrementing the approximation order
        by one.

        NOTE: Interpolation types (node-spacing) may need to be changed between
        the 'coarse' and 'fine' space for some high-order elements.

        Parameters
        ----------
        kwargs : keyword args
            arg=value pairs that will directly be passed to self.setupModel()
        """
        # copy the geometry info if necessary
        # geometry should be constant during the adaptive analysis
        if not self.fine.geomLoader:
            self.fine.geomLoader = copy.copy(self.coarse.geomLoader)

        # duplicate the forest, incrementing the order by 1
        # ensure interpolation type for high-order elements is properly set
        self.fine.forest = self.coarse.forest.duplicate()
        refined_order = self.coarse.forest.getMeshOrder() + 1
        if isinstance(self.fine.forest, TMR.QuadForest):
            if refined_order < 4:
                self.fine.forest.setMeshOrder(refined_order, TMR.UNIFORM_POINTS)
            else:
                self.fine.forest.setMeshOrder(refined_order, TMR.GAUSS_LOBATTO_POINTS)
        elif isinstance(self.fine.forest, TMR.OctForest):
            self.fine.forest.setMeshOrder(refined_order, TMR.UNIFORM_POINTS)

        # make sure the forest is balanced and shared with the geomLoader object
        self.fine.forest.balance(1)
        self.fine.geomLoader.forest = self.fine.forest

        # copy the refinement iteration counter
        self.fine.refine_iter = self.coarse.refine_iter

        # setup the fine-space model
        self.setupModel(model_type="fine", **kwargs)
        return

    def setOutputOfInterest(self, output_name):
        """
        Sets the output name used for adaptation in both the coarse-space and
        fine-space objects.

        Parameters
        ----------
        output_name : str
            Name of output to be used for output-based adaptation
        """
        setattr(self.coarse, "output_name", output_name)
        setattr(self.fine, "output_name", output_name)
        return

    def setupModel(self, model_type, **kwargs):
        """
        Sets up the model for the selected approximation space. After this
        function is called, the underlying finite-element assembler and problem
        in the model will be ready for use.

        Parameters
        ----------
        model_type : str in ('coarse', 'fine')
            Selector for which model to set up

        kwargs : keyword args
            arg=value pairs that will directly be passed to pyApproximationSpace.createProblem()
        """
        # select the appropriate model
        model = self._selectModel(model_type)

        # create the assembler
        model.createAssembler()

        # create and set up the problem used for analysis
        if "prob_name" not in kwargs:
            kwargs["prob_name"] = f"{model_type}_{model.refine_iter}"
        model.createProblem(**kwargs)

        # record the number of degrees of freedom for this model
        nnodes = model.assembler.getNumOwnedNodes()
        self.mesh_history[
            model.prob_name
        ] = model.assembler.getVarsPerNode() * model.comm.allreduce(nnodes, op=MPI.SUM)
        return

    def solvePrimal(self, model_type, writeSolution=False, **kwargs):
        """
        Solves the primal problem for the given model, copies the state to the
        model, and evaluates the output using the solved state. Optionally write
        out the solution data file.

        Parameters
        ----------
        model_type : str in ('coarse', 'fine')
            Selector for which model to solve the primal problem on

        writeSolution : bool
            Boolean flag for writing out the state field to a data file

        kwargs : keyword args
            arg=value pairs that will directly be passed to TACSProblem.solve()
        """
        # select the appropriate model
        model = self._selectModel(model_type)

        # solve the problem
        model.problem.solve(**kwargs)

        # copy the states
        model.problem.getVariables(model.state)

        # evaluate the output
        model.problem.evalFunctions(self.output_history, evalFuncs=[model.output_name])

        # write out the state field
        if writeSolution:
            model.problem.writeSolution(baseName=f"{model.prob_name}_state")
        return

    def solveAdjoint(self, model_type, writeSolution=False):
        """
        Solves the adjoint problem for the given model and copies the adjoint to
        the model. Optionally write out the solution data file.

        Parameters
        ----------
        model_type : str in ('coarse', 'fine')
            Selector for which model to solve the adjoint problem on

        writeSolution : bool
            Boolean flag for writing out the adjoint field to a data file
        """
        # select the appropriate model
        model = self._selectModel(model_type)

        # get the output sensitivity w.r.t. the states
        dJdu = model.assembler.createVec().getArray()
        model.problem.addSVSens([model.output_name], [dJdu])

        # solve for the adjoint variables
        model.adjoint.zeroEntries()
        model.problem.solveAdjoint(dJdu, model.adjoint)

        # write out the adjoint field
        if writeSolution:
            # overwrite the state field in the problem with the adjoint field
            model.problem.setVariables(model.adjoint)
            model.problem.writeSolution(baseName=f"{model.prob_name}_adjoint")
            # re-write the states back into the problem
            model.problem.setVariables(model.state)
        return

    def interpolateField(self, field_type, writeSolution=False):
        """
        Uses the information in the coarse-space model to interpolate the given
        field in the higher-order, fine-space model. The interpolated field
        is copied to the fine-space model. Optionally write out the solution
        data file.

        Parameters
        ----------
        field_type : str in ('state', 'adjoint')
            Selector for which field to interpolate in the fine-space model

        writeSolution : bool
            Boolean flag for writing out the selected field to a data file
        """
        field_type = self._selectField(field_type)

        # interpolate the selected field
        TMR.computeInterpSolution(
            self.coarse.forest,
            self.coarse.assembler,
            self.fine.forest,
            self.fine.assembler,
            getattr(self.coarse, field_type),
            getattr(self.fine, field_type),
        )

        # evaluate the fine-space output if the state is interpolated
        if field_type == "state":
            self.fine.problem.setVariables(getattr(self.fine, field_type))
            self.fine.problem.evalFunctions(
                self.output_history, evalFuncs=[self.fine.output_name]
            )

        # write out the selected field
        if writeSolution:
            # set the field in the fine problem
            self.fine.problem.setVariables(getattr(self.fine, field_type))
            self.fine.problem.writeSolution(
                baseName=f"{self.fine.prob_name}_{field_type}_interp"
            )
            # reset the state field in the fine problem
            self.fine.problem.setVariables(self.fine.state)
        return

    def reconstructField(self, field_type, compute_diff=False, writeSolution=False):
        """
        Uses the information in the coarse-space model to reconstruct the given
        field on the higher-order fine-space model. This enriches the field
        using the higher-order basis. The reconstructed field is copied to the
        fine model. Optionally store the difference between the reconstructed
        field and the interpolated field in the fine-space. Optionally write out
        the solution data file.

        Parameters
        ----------
        field_type : str in ('state', 'adjoint')
            Selector for which field to reconstruct using a higher-order basis
            in the fine-space model

        compute_diff : bool
            When True, computes and stores the difference between the reconstructed
            field and the interpolated field in the fine-space. Otherwise, just
            compute and store the reconstructed field in the fine-space.

        writeSolution : bool
            Boolean flag for writing out the selected field to a data file
        """
        field_type = self._selectField(field_type)

        # reconstruct the select field
        TMR.computeReconSolution(
            self.coarse.forest,
            self.coarse.assembler,
            self.fine.forest,
            self.fine.assembler,
            getattr(self.coarse, field_type),
            getattr(self.fine, field_type),
            compute_diff=compute_diff,
        )

        # write out the selected field
        if writeSolution:
            # set the field in the fine problem
            self.fine.problem.setVariables(getattr(self.fine, field_type))
            self.fine.problem.writeSolution(
                baseName=f"{self.fine.prob_name}_{field_type}_recon"
            )
            # reset the state field in the fine problem
            self.fine.problem.setVariables(self.fine.state)
        return

    def estimateOutputError(self):
        """
        Do the adjoint-based error estimation for the output of interest. Uses
        the current state and adjoint fields in the fine-space model for the
        estimation.

        Returns
        -------
        error_estimate : float
            The absolute error estimate in the output of interest

        output_correction : float
            The signed correction to the output of interest

        element_errors : np.ndarray (1D, float) [num_elements]
            The error in the output associated with each element in the mesh
        """
        # do the error estimation
        error_estimate, output_correction, element_errors = TMR.adjointError(
            self.coarse.forest,
            self.coarse.assembler,
            self.fine.forest,
            self.fine.assembler,
            self.fine.state,
            self.fine.adjoint,
        )

        # update histories
        self.output_history[
            f"{self.fine.prob_name}_{self.fine.output_name}"
        ] += output_correction
        self.error_history[f"adapt_iter_{self.fine.refine_iter}"] = error_estimate

        return error_estimate, output_correction, element_errors

    def adaptModelDecreasingThreshold(self, element_errors):
        """
        Refines elements in the coarse-space model based on a decreasing threshold
        adaptation strategy.

        Nemec, Marian, Michael Aftosmis, and Mathias Wintzer. "Adjoint-based
        adaptive mesh refinement for complex geometries." 46th AIAA Aerospace
        Sciences Meeting and Exhibit. 2008.

        Parameters
        ----------
        element_errors : np.ndarray (1D, float) [num_elements]
            The error in the output associated with each element in the mesh
        """
        # get the total number of elements
        nelems = len(element_errors)
        nelems_tot = self.comm.allreduce(nelems, op=MPI.SUM)

        # choose elements based on an equidistributed target error
        target_error = self.error_tol / nelems_tot
        error_ratio = element_errors / target_error
        refine_threshold = max(
            1.0, 2.0 ** (self.num_decrease_iters - self.coarse.refine_iter)
        )

        # record the refinement threshold:error_ratio pair for this iteration
        elem_counts = self.comm.allgather(nelems)
        total_error_ratios = np.zeros(nelems_tot)
        self.comm.Allgatherv(error_ratio, [total_error_ratios, elem_counts])
        self.adaptation_history[refine_threshold] = total_error_ratios

        # get the refinement indicator array
        refine = np.array(error_ratio >= refine_threshold, dtype=int)
        nref = np.count_nonzero(refine)
        nref = self.comm.allreduce(nref, op=MPI.SUM)

        # adapt the coarse-space model
        self.coarse.applyRefinement(
            refine.astype(np.intc),
            num_min_levels=self.num_min_ref_levels,
            num_max_levels=self.num_max_ref_levels,
        )
        return

    def _selectModel(self, model_type):
        """
        Helper function to check and select the appropriate model option
        """
        valid_model_types = ("coarse", "fine")
        model_type = model_type.lower()
        assert (
            model_type in valid_model_types
        ), f"model_type must be one of {valid_model_types}"
        return getattr(self, model_type)

    def _selectField(self, field_type):
        """
        Helper function to check and select the appropriate field option
        """
        valid_field_types = ("state", "adjoint")
        field_type = field_type.lower()
        assert (
            field_type in valid_field_types
        ), f"field_type must be one of {valid_field_types}"
        return field_type

    def printroot(self, msg):
        """
        Helper function to print a message on the root proc
        """
        if self.comm.rank == 0:
            print(f"{msg}")
        return

    def printproc(self, msg):
        """
        Helper function to print a message on all procs
        """
        print(f"[{self.comm.rank}] {msg}")
        return
