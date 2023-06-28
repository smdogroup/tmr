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
import h5py
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

    def writeField(self, field, field_name, reset=True):
        """
        Sets the problem variables with the supplied field and writes out the
        solution data file. Uses the field name to name the file.
        Optionally, will reset the problem variables with the current state
        field if specified.

        Parameters
        ----------
        field : TACS Vec or np.ndarray (1D, float)
            Vector to set as the problem variables when writing out the file

        field_name : str
            Name of the field to add to the end of the file name

        reset : bool
            Flag for resetting the problem variables with the current state
        """
        # set the field in the problem variables
        self.problem.setVariables(field)
        self.problem.writeSolution(baseName=f"{self.prob_name}_{field_name}")

        if reset:
            # re-write the states back into the problem
            self.problem.setVariables(self.state)
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
        self,
        comm,
        _initCallBack,
        _elemCallBack,
        _probCallBack,
        _adapt_output,
        adapt_params={},
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

        adapt_params : dict
            arg=value pairs that can be used to specify adaptation parameters
        """
        # Set MPI communicator
        BaseUI.__init__(self, comm=comm)

        # initialize the class attributes
        self.initCallBack = _initCallBack
        self.coarse = pyApproximationSpace(comm, _elemCallBack, _probCallBack)
        self.fine = pyApproximationSpace(comm, _elemCallBack, _probCallBack)
        self.adapt_output = _adapt_output

        # set the adaptation parameters
        valid_strategies = ("decreasing_threshold", "fixed_growth")
        self.adapt_strategy = adapt_params.get("adapt_strategy", "").lower()
        if self.adapt_strategy:
            assert (
                self.adapt_strategy in valid_strategies
            ), f"adaptation strategy must be one of {valid_strategies}"
            # get the general adaptation params
            self.num_min_ref_levels = adapt_params.get("num_min_ref_levels", 0)
            self.num_max_ref_levels = adapt_params.get(
                "num_max_ref_levels", TMR.MAX_LEVEL
            )
            self.error_tol = adapt_params.get("error_tol")
            assert (
                self.error_tol is not None
            ), f"error tolerance must be specified for adaptation"
            # get strategy-specific params
            if self.adapt_strategy == "decreasing_threshold":
                self.num_decrease_iters = adapt_params.get("num_decrease_iters")
                assert (
                    self.num_decrease_iters is not None
                ), f"number of decrease iterations must be specified for adaptation strategy `{self.adapt_strategy}`"
            if self.adapt_strategy == "fixed_growth":
                self.growth_refine_factor = np.clip(
                    adapt_params.get("growth_refine_factor", 0.0), 0.0, 1.0
                )
                self.growth_coarsen_factor = np.clip(
                    adapt_params.get("growth_coarsen_factor", 0.0), 0.0, 1.0
                )

        # initialize the storage for saving model info
        self.mesh_history = {}  # mesh size (ndof)
        self.output_history = {}  # outputs/corrected outputs
        self.error_history = {}  # error estimates
        self.adaptation_history = {
            "threshold": {},  # refine/coarsen thresholds
            "element_errors": {},
        }  # element-wise errors

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
        else:
            kwargs["prob_name"] += f"_{model_type}_{model.refine_iter}"
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
            model.writeField(model.state, "state")
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
            model.writeField(model.adjoint, "adjoint")
        return

    def evaluateResidual(self, model_type, vec=None, writeSolution=False):
        """
        Evaluates the residual of the specified model for a given perturbation.
        Optionally write out the residual to a data file.

        Parameters
        ----------
        model_type : str in ('coarse', 'fine')
            Selector for which model to evaluate the residual on

        vec : TACS Vec or np.ndarray (1D, float)
            Vector to add to rhs when evaluating the residual

        writeSolution : bool
            Boolean flag for writing out the residual field to a data file

        Returns
        -------
        res : TACS Vec
            Residual vector
        """
        # select the appropriate model
        model = self._selectModel(model_type)

        # create the residual vector
        res = model.assembler.createVec()
        res.zeroEntries()

        # evaluate the residual
        model.problem.getResidual(res, Fext=vec)

        # write out the residual field
        if writeSolution:
            model.writeField(res, "residual")

        return res

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
            self.fine.writeField(getattr(self.fine, field_type), f"{field_type}_interp")
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
            self.fine.writeField(getattr(self.fine, field_type), f"{field_type}_recon")
        return

    def estimateOutputError(self, writeSolution=False):
        """
        Do the adjoint-based error estimation for the output of interest. Uses
        the current state and adjoint fields in the fine-space model for the
        estimation. Optionally output the nodal error field to a data file.

        NOTE: Nodal error field is a scalar and is written to the "u" variable.

        Parameters
        ----------
        writeSolution : bool
            Boolean flag for writing out node error field to a data file

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
        (
            error_estimate,
            output_correction,
            element_errors,
            node_errors,
        ) = TMR.adjointError(
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

        # write out the nodal error field
        if writeSolution:
            vpn = self.fine.assembler.getVarsPerNode()
            node_error_vars = np.zeros(vpn * len(node_errors))
            node_error_vars[::vpn] = node_errors
            self.fine.writeField(field=node_error_vars, field_name="error")

        return error_estimate, output_correction, element_errors, node_errors

    def collectErrors(self, element_errors=None, node_errors=None):
        """
        Collects errors that are distributed across all the processors

        Parameters
        ----------
        element_errors : np.ndarray (1D, float)
            Distributed array of element errors

        node_errors : np.ndarray (1D, float)
            Distributed array of node errors

        Returns
        -------
        total_errors : np.ndarray (1D, float) or list[np.ndarray (1D, float)]
            Collected errors from all processors. If both element and nodal
            errors are collected, element errors are the first element in the
            returned list.
        """
        total_errors = []
        if element_errors is not None:
            # collect the error localized to the distributed elements
            elem_counts = self.comm.allgather(len(element_errors))
            element_errors_tot = np.empty(sum(elem_counts))
            self.comm.Allgatherv(element_errors, [element_errors_tot, elem_counts])
            total_errors.append(element_errors_tot)
        if node_errors is not None:
            # collect the error distributed to the nodes
            node_counts = self.comm.allgather(len(node_errors))
            node_errors_tot = np.empty(sum(node_counts))
            self.comm.Allgatherv(node_errors, [node_errors_tot, node_counts])
            total_errors.append(node_errors_tot)
        if len(total_errors) == 1:
            total_errors = total_errors[0]
        return total_errors

    def refineModel(self, element_errors=None):
        """
        Applies mesh refinement to the model using the selected adaptation
        strategy and a distribution of element-wise errors. If adaptation is not
        being done, uniform refinement is applied

        Parameters
        ----------
        element_errors : np.ndarray (1D, float) [num_elements]
            The error in the output associated with each element in the mesh
        """
        if self.adapt_strategy == "decreasing_threshold":
            self.adaptModelDecreasingThreshold(element_errors)
        elif self.adapt_strategy == "fixed_growth":
            self.adaptModelFixedGrowth(element_errors)
        else:
            # apply uniform refinement to the coarse model
            self.coarse.applyRefinement()
        return

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
        elem_counts = self.comm.allgather(nelems)
        nelems_tot = sum(elem_counts)

        # choose elements based on an equidistributed target error
        target_error = self.error_tol / nelems_tot
        error_ratio = element_errors / target_error
        refine_threshold = max(
            1.0, 2.0 ** (self.num_decrease_iters - self.coarse.refine_iter)
        )

        # record the refinement threshold:error_ratio pair for this iteration
        error_ratio_tot = np.zeros(nelems_tot)
        self.comm.Allgatherv(error_ratio, [error_ratio_tot, elem_counts])
        self.adaptation_history["threshold"][
            f"refine_{self.coarse.refine_iter}"
        ] = refine_threshold
        self.adaptation_history["element_errors"][
            f"adapt_iter_{self.coarse.refine_iter}"
        ] = error_ratio_tot

        # get the refinement indicator array
        adapt_indicator = np.array(error_ratio >= refine_threshold, dtype=int)
        nref = np.count_nonzero(adapt_indicator)
        nref = self.comm.allreduce(nref, op=MPI.SUM)

        # adapt the coarse-space model
        self.coarse.applyRefinement(
            adapt_indicator.astype(np.intc),
            num_min_levels=self.num_min_ref_levels,
            num_max_levels=self.num_max_ref_levels,
        )
        return

    def adaptModelFixedGrowth(self, element_errors):
        """
        Refines elements in the coarse-space model based on a fixed growth
        adaptation strategy.

        Parameters
        ----------
        element_errors : np.ndarray (1D, float) [num_elements]
            The error in the output associated with each element in the mesh
        """
        # get the elem counts
        nelems = len(element_errors)
        elem_counts = self.comm.allgather(nelems)
        nelems_tot = sum(elem_counts)
        nrefine = int(self.growth_refine_factor * nelems_tot)
        ncoarse = int(self.growth_coarsen_factor * nelems_tot)

        # gather all the element errors on each proc
        element_errors_tot = np.empty(nelems_tot)
        self.comm.Allgatherv(element_errors, [element_errors_tot, elem_counts])

        # compute the error ratio for all elements
        target_error = self.error_tol / nelems_tot
        error_ratio_tot = element_errors_tot / target_error

        # sort all the element error ratios in descending order
        sorted_inds = np.flip(np.argsort(error_ratio_tot))

        # create the adaptation indicator array
        adapt_indicator_tot = np.zeros(nelems_tot, dtype=int)
        refine_inds = None
        coarse_inds = None
        if nrefine > 0:
            refine_inds = sorted_inds[:nrefine]
            adapt_indicator_tot[refine_inds] = np.where(
                error_ratio_tot[refine_inds] >= 1.0, 1, adapt_indicator_tot[refine_inds]
            )
        if ncoarse > 0:
            coarse_inds = sorted_inds[-ncoarse:]
            adapt_indicator_tot[coarse_inds] = np.where(
                error_ratio_tot[coarse_inds] < 1.0, -1, adapt_indicator_tot[coarse_inds]
            )

        # update the adaptation history
        if refine_inds is not None:
            ind = 0
            while error_ratio_tot[refine_inds[ind]] >= 1.0:
                self.adaptation_history["threshold"][
                    f"refine_{self.coarse.refine_iter}"
                ] = error_ratio_tot[refine_inds[ind]]
                ind += 1
                if ind >= nrefine:
                    break
        if coarse_inds is not None:
            ind = ncoarse - 1
            while error_ratio_tot[coarse_inds[ind]] < 1.0:
                self.adaptation_history["threshold"][
                    f"coarsen_{self.coarse.refine_iter}"
                ] = error_ratio_tot[coarse_inds[ind]]
                ind -= 1
                if ind <= 0:
                    break
        self.adaptation_history["element_errors"][
            f"adapt_iter_{self.coarse.refine_iter}"
        ] = error_ratio_tot

        # scatter the adaptation indicator to each proc
        adapt_indicator = np.empty(nelems, dtype=int)
        self.comm.Scatterv([adapt_indicator_tot, elem_counts], adapt_indicator, root=0)

        # adapt the coarse-space model
        self.coarse.applyRefinement(
            adapt_indicator.astype(np.intc),
            num_min_levels=self.num_min_ref_levels,
            num_max_levels=self.num_max_ref_levels,
        )
        return

    def writeModelHistory(self, filename=""):
        """
        Writes out the model history information to an .hdf5 file

        Parameters
        ----------
        filename : str
            name of output file, should include the proper .hdf5 file extension
        """
        if self.comm.rank == 0:
            if not filename:
                filename = "model_history.hdf5"
            with h5py.File(filename, "w") as h5:
                # store the mesh histories
                h5["mesh_history/coarse"] = np.array(
                    [
                        self.mesh_history[key]
                        for key in self.mesh_history.keys()
                        if "coarse" in key
                    ]
                )
                h5["mesh_history/fine"] = np.array(
                    [
                        self.mesh_history[key]
                        for key in self.mesh_history.keys()
                        if "fine" in key
                    ]
                )

                # store the output histories
                h5["output_history/coarse"] = np.array(
                    [
                        self.output_history[key]
                        for key in self.output_history.keys()
                        if "coarse" in key
                    ]
                )
                h5["output_history/fine"] = np.array(
                    [
                        self.output_history[key]
                        for key in self.output_history.keys()
                        if "fine" in key
                    ]
                )
                h5["output_history"].attrs["output_name"] = self.adapt_output.replace(
                    "_", " "
                ).upper()

                # store the error estimate history
                h5["error_estimate_history"] = np.array(
                    list(self.error_history.values())
                )

                # store the adaptation refinement histories
                h5["adaptation_history/threshold/refine"] = np.array(
                    [
                        self.adaptation_history["threshold"][key]
                        for key in self.adaptation_history["threshold"].keys()
                        if "refine" in key
                    ]
                )
                h5["adaptation_history/threshold/coarsen"] = np.array(
                    [
                        self.adaptation_history["threshold"][key]
                        for key in self.adaptation_history["threshold"].keys()
                        if "coarsen" in key
                    ]
                )
                for key, val in self.adaptation_history["element_errors"].items():
                    h5[f"adaptation_history/element_errors/{key}"] = val
                if self.adapt_strategy:
                    h5["adaptation_history"].attrs[
                        "strategy"
                    ] = self.adapt_strategy.replace("_", " ").lower()
                    h5["adaptation_history"].attrs["error_tolerance"] = self.error_tol
                    h5["adaptation_history"].attrs[
                        "max_refine_levels"
                    ] = self.num_max_ref_levels
                    h5["adaptation_history"].attrs[
                        "min_refine_levels"
                    ] = self.num_min_ref_levels
                    if hasattr(self, "num_decrease_iters"):
                        h5["adaptation_history"].attrs[
                            "threshold_decrease_iters"
                        ] = self.num_decrease_iters
                    if hasattr(self, "growth_refine_factor"):
                        h5["adaptation_history"].attrs[
                            "growth_refine_factor"
                        ] = self.growth_refine_factor
                    if hasattr(self, "growth_coarsen_factor"):
                        h5["adaptation_history"].attrs[
                            "growth_coarsen_factor"
                        ] = self.growth_coarsen_factor
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
