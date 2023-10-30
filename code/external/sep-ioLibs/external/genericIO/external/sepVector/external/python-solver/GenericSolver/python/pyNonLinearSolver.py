from collections import deque
from math import isnan
import pySolver
import pyOperator as pyOp
from pyStepper import CvSrchStep, ParabolicStep
from pyStopper import BasicStopper
from pyProblem import ProblemLinearSymmetric
from pyLinearSolver import SymLCGsolver


# Beta functions
# grad=new gradient, grad0=old, dir=search direction
# From A SURVEY OF NONLINEAR CONJUGATE GRADIENT METHODS

def _betaFR(grad, grad0, dir, logger):
    """Fletcher and Reeves method"""
    # betaFR = sum(dprod(g,g))/sum(dprod(g0,g0))
    dot_grad = grad.dot(grad)
    dot_grad0 = grad0.dot(grad0)
    if dot_grad0 == 0.:  # Avoid division by zero
        beta = 0.
        if logger:
            logger.addToLog("Setting beta to zero since norm of previous gradient is zero!!!")
    else:
        beta = dot_grad / dot_grad0
    return beta


def _betaPRP(grad, grad0, dir, logger):
    """Polak, Ribiere, Polyak method"""
    # betaPRP = sum(dprod(g,g-g0))/sum(dprod(g0,g0))
    tmp1 = grad.clone()
    # g-g0
    tmp1.scaleAdd(grad0, 1.0, -1.0)
    dot_num = tmp1.dot(grad)
    dot_grad0 = grad0.dot(grad0)
    if dot_grad0 == 0.:  # Avoid division by zero
        beta = 0.
        if logger:
            logger.addToLog("Setting beta to zero since norm of previous gradient is zero!!!")
    else:
        beta = dot_num / dot_grad0
    return beta


def _betaHS(grad, grad0, dir, logger):
    """Hestenes and Stiefel"""
    # betaHS = sum(dprod(g,g-g0))/sum(dprod(d,g-g0))
    tmp1 = grad.clone()
    # g-g0
    tmp1.scaleAdd(grad0, 1.0, -1.0)
    dot_num = tmp1.dot(grad)
    dot_denom = tmp1.dot(dir)
    if dot_denom == 0.:  # Avoid division by zero
        beta = 0.
        if logger:
            logger.addToLog("Setting beta to zero since norm of denominator is zero!!!")
    else:
        beta = dot_num / dot_denom
    return beta


def _betaCD(grad, grad0, dir, logger):
    """Conjugate Descent"""
    # betaCD = -sum(dprod(g,g))/sum(dprod(d,g0))
    dot_num = grad.dot(grad)
    dot_denom = -grad0.dot(dir)
    if dot_denom == 0.:  # Avoid division by zero
        beta = 0.
        if logger:
            logger.addToLog("Setting beta to zero since norm of denominator is zero!!!")
    else:
        beta = dot_num / dot_denom
    return beta


def _betaLS(grad, grad0, dir, logger):
    """Liu and Storey"""
    # betaLS = -sum(dprod(g,g-g0))/sum(dprod(d,g0))
    tmp1 = grad.clone()
    # g-g0
    tmp1.scaleAdd(grad0, 1.0, -1.0)
    dot_num = tmp1.dot(grad)
    dot_denom = -grad0.dot(dir)
    if dot_denom == 0.:  # Avoid division by zero
        beta = 0.
        if logger:
            logger.addToLog("Setting beta to zero since norm of denominator is zero!!!")
    else:
        beta = dot_num / dot_denom
    return beta


def _betaDY(grad, grad0, dir, logger):
    """Dai and Yuan"""
    # betaDY = sum(dprod(g,g))/sum(dprod(d,g-g0))
    tmp1 = grad.clone()
    # g-g0
    tmp1.scaleAdd(grad0, 1.0, -1.0)
    dot_num = grad.dot(grad)
    dot_denom = tmp1.dot(dir)
    if dot_denom == 0.:  # Avoid division by zero
        beta = 0.
        if logger:
            logger.addToLog("Setting beta to zero since norm of denominator is zero!!!")
    else:
        beta = dot_num / dot_denom
    return beta


def _betaBAN(grad, grad0, dir, logger):
    """Bamigbola, Ali and Nwaeze"""
    # betaDY = sum(dprod(g,g-g0))/sum(dprod(g0,g-g0))
    tmp1 = grad.clone()
    # g-g0
    tmp1.scaleAdd(grad0, 1.0, -1.0)
    dot_num = tmp1.dot(grad)
    dot_denom = tmp1.dot(grad0)
    if dot_denom == 0.:  # Avoid division by zero
        beta = 0.
        if logger:
            logger.addToLog("Setting beta to zero since norm of denominator is zero!!!")
    else:
        beta = -dot_num / dot_denom
    return beta


def _betaHZ(grad, grad0, dir, logger):
    """Hager and Zhang"""
    # betaN = sum(dprod(g-g0-2*sum(dprod(g-g0,g-g0))*d/sum(dprod(d,g-g0)),g))/sum(dprod(d,g-g0))
    tmp1 = grad.clone()
    # g-g0
    tmp1.scaleAdd(grad0, 1.0, -1.0)
    # sum(dprod(g-g0,g-g0))
    dot_diff_g_g0 = tmp1.dot(tmp1)
    # sum(dprod(d,g-g0))
    dot_dir_diff_g_g0 = tmp1.dot(dir)
    if dot_dir_diff_g_g0 == 0.:  # Avoid division by zero
        beta = 0.
        if logger:
            logger.addToLog("Setting beta to zero since norm of denominator is zero!!!")
    else:
        # g-g0-2*sum(dprod(g-g0,g-g0))*d/sum(dprod(d,g-g0))
        tmp1.scaleAdd(dir, 1.0, -2.0 * dot_diff_g_g0 / dot_dir_diff_g_g0)
        # sum(dprod(g-g0-2*sum(dprod(g-g0,g-g0))*d/sum(dprod(d,g-g0)),g))
        dot_num = grad.dot(tmp1)
        # dot_num/sum(dprod(d,g-g0))
        beta = dot_num / dot_dir_diff_g_g0
    return beta


def _betaSD(grad, grad0, dir, logger):
    """Steepest descent"""
    beta = 0.
    return beta


class NLCGsolver(pySolver.Solver):
    """Non-Linear Conjugate Gradient and Steepest-Descent Solver object"""
    
    # Default class methods/functions
    def __init__(self, stoppr, stepper=None, beta_type="FR", logger=None):
        """Constructor for NLCG Solver"""
        # Calling parent construction
        super(NLCGsolver, self).__init__()
        # Defining stopper object
        self.stoppr = stoppr
        # Defining stepper object
        self.stepper = stepper if stepper is not None else ParabolicStep()
        # Beta function to use during the inversion
        self.beta_type = beta_type
        # Logger object to write on log file
        self.logger = logger
        # Overwriting logger of the Stopper object
        self.stoppr.logger = self.logger
        # print formatting
        self.iter_msg = "iter = %s, obj = %.5e, resnorm = %.2e, gradnorm = %.2e, feval = %d, geval = %d"
        return
    
    def __del__(self):
        """Default destructor"""
        return
    
    def beta_func(self, grad, grad0, dir):
        """Beta function interface"""
        beta_type = self.beta_type
        if beta_type == "FR":
            beta = _betaFR(grad, grad0, dir, self.logger)
        elif beta_type == "PRP":
            beta = _betaPRP(grad, grad0, dir, self.logger)
        elif beta_type == "HS":
            beta = _betaHS(grad, grad0, dir, self.logger)
        elif beta_type == "CD":
            beta = _betaCD(grad, grad0, dir, self.logger)
        elif beta_type == "LS":
            beta = _betaLS(grad, grad0, dir, self.logger)
        elif beta_type == "DY":
            beta = _betaDY(grad, grad0, dir, self.logger)
        elif beta_type == "BAN":
            beta = _betaBAN(grad, grad0, dir, self.logger)
        elif beta_type == "HZ":
            beta = _betaHZ(grad, grad0, dir, self.logger)
        elif beta_type == "SD":
            beta = _betaSD(grad, grad0, dir, self.logger)
        else:
            raise ValueError("ERROR! Requested Beta function type not existing")
        return beta
    
    def run(self, problem, verbose=False, restart=False):
        """Running NLCG or steppest-descent solver"""
        
        self.create_msg = verbose or self.logger
        
        # Resetting stopper before running the inversion
        self.stoppr.reset()
        
        if not restart:
            if self.create_msg:
                msg = 90 * "#" + "\n"
                msg += "\t\t\tNON-LINEAR %s SOLVER log file\n" % (
                    "STEEPEST-DESCENT" if self.beta_type == "SD" else "CONJUGATE GRADIENT")
                msg += "\tRestart folder: %s\n" % self.restart.restart_folder
                if self.beta_type != "SD":
                    msg += "\tConjugate method used: %s\n" % self.beta_type
                msg += 90 * "#" + "\n"
                if verbose:
                    print(msg.replace("log file", ""))
                if self.logger:
                    self.logger.addToLog(msg)
            
            # Setting internal vectors (model, search direction, and previous gradient vectors)
            prblm_mdl = problem.get_model()
            cg_mdl = prblm_mdl.clone()
            cg_dmodl = prblm_mdl.clone()
            cg_dmodl.zero()
            cg_grad0 = cg_dmodl.clone()
            
            # Other internal variables
            beta = 0.0
            iiter = 0
        else:
            # Retrieving parameters and vectors to restart the solver
            if self.create_msg:
                msg = "Restarting previous solver run from: %s" % self.restart.restart_folder
                if verbose:
                    print(msg)
                if self.logger:
                    self.logger.addToLog(msg)
            self.restart.read_restart()
            iiter = self.restart.retrieve_parameter("iter")
            self.stepper.alpha = self.restart.retrieve_parameter("alpha")
            initial_obj_value = self.restart.retrieve_parameter(
                "obj_initial")  # Retrieving initial objective function value
            cg_mdl = self.restart.retrieve_vector("cg_mdl")
            cg_dmodl = self.restart.retrieve_vector("cg_dmodl")
            cg_grad0 = self.restart.retrieve_vector("cg_grad0")
            # Setting the model and residuals to avoid residual twice computation
            problem.set_model(cg_mdl)
            prblm_mdl = problem.get_model()
            # Setting residual vector to avoid its unnecessary computation
            problem.set_residual(self.restart.retrieve_vector("prblm_res"))
        
        # Common variables unrelated to restart
        prev_mdl = prblm_mdl.clone().zero()
        
        while True:
            # Computing objective function
            obj0 = problem.get_obj(cg_mdl)  # Compute objective function value
            prblm_res = problem.get_res(cg_mdl)  # Compute residuals
            prblm_grad = problem.get_grad(cg_mdl)  # Compute the gradient
            if iiter == 0:
                initial_obj_value = obj0  # For relative objective function value
                # Saving objective function value
                self.restart.save_parameter("obj_initial", initial_obj_value)
                # iteration info
                if self.create_msg:
                    msg = self.iter_msg % (str(iiter).zfill(self.stoppr.zfill),
                                           obj0,
                                           problem.get_rnorm(cg_mdl),
                                           problem.get_gnorm(cg_mdl),
                                           problem.get_fevals(),
                                           problem.get_gevals())
                    # Writing on log file
                    if verbose:
                        print(msg)
                    if self.logger:
                        self.logger.addToLog(msg)
                # Check if either objective function value or gradient norm is NaN
                if isnan(obj0) or isnan(prblm_grad.norm()):
                    raise ValueError("ERROR! Either gradient norm or objective function value NaN!")
            if prblm_grad.norm() == 0.:
                print("Gradient vanishes identically")
                break
            
            # Saving results
            self.save_results(iiter, problem, force_save=False)
            # Keeping current inverted model
            prev_mdl.copy(prblm_mdl)
            
            if iiter >= 1:
                beta = self.beta_func(prblm_grad, cg_grad0, cg_dmodl)
                if beta < 0.:
                    if self.logger:
                        self.logger.addToLog("Beta negative setting to zero: beta value=%s!!!" % beta)
                    beta = 0.
            if self.beta_type != "SD":
                if self.logger:
                    self.logger.addToLog("beta coefficient: %s" % beta)
            
            # dmodl = beta*dmodl - grad
            cg_dmodl.scaleAdd(prblm_grad, beta, -1.0)
            # grad0 = grad
            cg_grad0.copy(prblm_grad)
            # Calling line search
            alpha, success = self.stepper.run(problem, cg_mdl, cg_dmodl, self.logger)
            if not success:
                if self.create_msg:
                    msg = "Stepper couldn't find a proper step size, will terminate solver"
                    if verbose:
                        print(msg)
                    # Writing on log file
                    if self.logger:
                        self.logger.addToLog(msg)
                problem.set_model(prev_mdl)
                break
            
            # Increasing iteration counter
            iiter = iiter + 1
            obj1 = problem.get_obj(cg_mdl)  # Compute objective function value
            # Redundant test on verifying convergence
            if obj0 <= obj1:
                if self.create_msg:
                    msg = "Objective function at new point greater or equal than previous one:\n\t" \
                          "obj_new = %.5e\tobj_cur = %.5e\n" \
                          "Potential issue in the stepper or in revaluation of objective function! Solver will stop!" \
                          % (obj1, obj0)
                    if verbose:
                        print(msg)
                    if self.logger:
                        self.logger.addToLog(msg)
                problem.set_model(prev_mdl)
                break
            
            # Saving current model and previous search direction in case of restart
            self.restart.save_parameter("iter", iiter)
            self.restart.save_parameter("alpha", alpha)
            self.restart.save_vector("cg_mdl", cg_mdl)
            self.restart.save_vector("cg_dmodl", cg_dmodl)
            self.restart.save_vector("cg_grad0", cg_grad0)
            # Saving data space vectors
            self.restart.save_vector("prblm_res", prblm_res)
            # iteration info
            if self.create_msg:
                msg = self.iter_msg % (str(iiter).zfill(self.stoppr.zfill),
                                       obj1,
                                       problem.get_rnorm(cg_mdl),
                                       problem.get_gnorm(cg_mdl),
                                       problem.get_fevals(),
                                       problem.get_gevals())
                if verbose:
                    print(msg)
                if self.logger:
                    self.logger.addToLog("\n" + msg)
            
            # Check if either objective function value or gradient norm is NaN
            if isnan(obj1) or isnan(prblm_grad.norm()):
                raise ValueError("ERROR! Either gradient norm or objective function value NaN!")
            if self.stoppr.run(problem, iiter, initial_obj_value, verbose):
                break
        
        # Writing last inverted model
        self.save_results(iiter, problem, force_save=True, force_write=True)
        if self.create_msg:
            msg = 90 * "#" + "\n"
            msg += "\t\t\tNON-LINEAR %s SOLVER log file end\n" % (
                "STEEPEST-DESCENT" if self.beta_type == "SD" else "CONJUGATE GRADIENT")
            msg += 90 * "#" + "\n"
            if verbose:
                print(msg.replace(" log file", ""))
            if self.logger:
                self.logger.addToLog(msg)
        
        # Clear restart object
        self.restart.clear_restart()
        
        return


class TNewtonsolver(pySolver.Solver):
    """Truncated Newton/Gauss-Newton solver object"""
    
    def __init__(self, stopper, niter_max, HessianOp, stepper=None, niter_min=None, warm_start=True, Newton_prefix=None,
                 logger=None):
        """
        Constructor for Truncated-Newton Solver.
        stoppr     = [no default] - stopper class; Stopper object necessary to terminate the solver
        niter_max  = [no default] - int; Maximum number of iterations for solving Newton system when starting with a zero initial model
        HessianOp  = [no default] - operator class; Operator to apply Hessian matrix onto model vector. The operator must contain a set_background function to set model vector on which the Hessian is computed. Note the symmetric solver will be used. Hence, this operator must be symmetric.
        stepper    = [CvSrch] - stepper class; Stepper object necessary to perform line-search step
        niter_min  = [niter_max] - int; Number of iterations for solving Newton system when linear inversion starts from previous search direction
        warm_start = [False] - boolean; If True, the linear Hessian inversion is started from the previous search direction if aligned with the current gradient
        logger 	  = [None] - logger class; Logger object to save inversion information at runtime
        """
        # Defining stopper
        self.stopper = stopper  # Stopper for non-linear problem
        # Setting maximum and minimum number of iterations
        self.niter_max = niter_max
        self.niter_min = niter_max
        # Setting linear inversion iterations
        if niter_min is not None:
            if niter_min <= niter_max:
                raise ValueError(
                    "niter_min of %d must be smaller or equal than niter_max of %d." % (niter_min, niter_max))
            self.niter_min = niter_min
        # Defining stepper object
        self.stepper = stepper if stepper is not None else CvSrchStep()
        # Warm starts requested?
        self.warm_start = warm_start
        # Hessian operator
        if HessianOp is not None:
            if "set_background" not in dir(HessianOp):
                raise AttributeError("Hessian operator must have a set_background function.")
        # Setting linear solver for solving Newton system and problem class
        StopLin = BasicStopper(niter=self.niter_max)
        self.lin_solver = SymLCGsolver(StopLin)
        self.NewtonPrblm = ProblemLinearSymmetric(HessianOp.domain.clone(), HessianOp.domain.clone(), HessianOp)
        return
    
    def run(self, problem, verbose=False, restart=False):
        """Running Truncated Newton solver"""
        return


class LBFGSsolver(pySolver.Solver):
    """L-BFGS (Limited-memory Broyden-Fletcher-Goldfarb-Shanno) Solver object"""
    
    def __init__(self, stopper, stepper=None, save_alpha=False, m_steps=None, H0=None, logger=None, save_est=False):
        """
		Constructor for LBFGS Solver:
		:param stopper    : Stopper, object to terminate the solver
		:param stepper    : Stepper, object to perform line-search step
		:param save_alpha : bool, Use previous step-length value as initial guess.
								Otherwise, the algorithm starts from an initial guess of 1.0 [False]
		:param m_steps    : int, Maximum number of steps to store to estimate the inverse Hessian (by default it runs BFGS method)
		:param H0         : Operator, initial estimated Hessian inverse (by default it assumes an identity operator)
		:param logger 	  : Logger, object to save inversion information at runtime
		:param save_est   : bool, save inverse Hessian estimate vectors (self.prefix must not be None) [False]
		"""
        # Calling parent construction
        super(LBFGSsolver, self).__init__()
        # Defining stopper object
        self.stopper = stopper
        # Defining stepper object
        self.stepper = stepper if stepper is not None else CvSrchStep()
        # Logger object to write on log file
        self.logger = logger
        # Overwriting logger of the Stopper object
        self.stopper.logger = self.logger
        # LBFGS-specific parameters
        self.save_alpha = save_alpha
        self.H0 = H0
        self.m_steps = m_steps
        self.save_est = save_est
        self.tmp_vector = None  # A copy of the model vector will be create when the function run is invoked
        # print formatting
        self.iter_msg = "iter = %s, obj = %.5e, resnorm = %.2e, gradnorm = %.2e, feval = %d, geval = %d"
    
    def save_hessian_estimate(self, index, iiter):
        """Function to save current vector of estimated Hessian inverse"""
        # index of the step to save
        if self.prefix is not None and self.save_est:
            step_filename = self.prefix + "step_vector_%s.H" % iiter
            grad_diff_filename = self.prefix + "grad_diff_vector_%s.H" % iiter
            self.step_vectors[index].writeVec(step_filename)
            self.grad_diff_vectors[index].writeVec(grad_diff_filename)
    
    def check_rho(self, denom_dot, step_index, iiter):
        """Function to check scaling factor of Hessian inverse estimate"""
        if denom_dot == 0.:
            if self.m_steps is not None:
                self.rho[step_index] = 0.
            else:
                self.rho.append(0.)
            msg = "Skipping update to estimated Hessian; y vector orthogonal to s vector at iteration %s" % iiter
            if self.logger:
                self.logger.addToLog(msg)
        elif denom_dot < 0.:
            if self.m_steps is not None:
                self.rho[step_index] = 0.
            else:
                self.rho.append(0.)
            msg = "Skipping update to estimated Hessian; not positive at iteration %s" % iiter
            if self.logger:
                self.logger.addToLog(msg)
        else:
            if self.m_steps is not None:
                self.rho[step_index] = 1.0 / denom_dot
            else:
                self.rho.append(1.0 / denom_dot)
            # Saving current update for inverse Hessian estimate (i.e., gradient-difference and model-step vectors)
            self.save_hessian_estimate(step_index, iiter)
        return
    
    # BFGSMultiply function
    def BFGSMultiply(self, dmodl, grad, iiter):
        """Function to apply approximated inverse Hessian"""
        # Array containing dot-products
        if self.m_steps is not None:
            alpha = [0.0] * self.m_steps
            # Handling of limited memory
            if iiter <= self.m_steps:
                initial_point = 0
            else:
                initial_point = iiter % self.m_steps
            # Right step list
            rloop = deque(range(0, min(iiter, self.m_steps)))
            rloop.reverse()
            # Rotate the list
            rloop.rotate(initial_point)
            # Left step list
            lloop = deque(range(0, min(iiter, self.m_steps)))
            # Rotate the list
            lloop.rotate(-initial_point)
        else:
            alpha = [0.0] * iiter
            rloop = deque(range(0, iiter))
            rloop.reverse()
            lloop = deque(range(0, iiter))
        # r = -grad
        dmodl.copy(grad)
        dmodl.scale(-1.0)
        # Apply right-hand series of operators
        for ii in rloop:
            # Check positivity, if not true skip the update
            if self.rho[ii] > 0.0:
                # alpha_i=rho_i*s_i'r
                alpha[ii] = self.rho[ii] * self.step_vectors[ii].dot(dmodl)
                # r=r-alpha_i*y_i
                dmodl.scaleAdd(self.grad_diff_vectors[ii], 1.0, -alpha[ii])
        # Comput center (If not provide Identity matrix is assumed)
        # r=H0r
        if self.H0 is not None:
            # Apply a forward of the initial Hessian estimate
            self.H0.forward(False, dmodl, self.tmp_vector)
            dmodl.copy(self.tmp_vector)
        # Apply left-hand series of operators
        for ii in lloop:
            # Check positivity, if not true skip the update
            if self.rho[ii] > 0.0:
                # beta=rhoiyi'r
                beta = self.rho[ii] * self.grad_diff_vectors[ii].dot(dmodl)
                dmodl.scaleAdd(self.step_vectors[ii], 1.0, alpha[ii] - beta)
        return
    
    def run(self, problem, verbose=False, restart=False):
        """Running LBFGS solver"""
        # Resetting stopper before running the inversion
        self.stopper.reset()
        # Preliminary variables for Hessian inverse estimation
        if self.m_steps is not None:
            self.step_vectors = [None] * self.m_steps  # s_i vectors
            self.grad_diff_vectors = [None] * self.m_steps  # y_i vectors
            self.rho = [None] * self.m_steps  # Scalar term necessary for Hessian inverse estimation
        else:
            self.step_vectors = []  # s_i vectors
            self.grad_diff_vectors = []  # y_i vectors
            self.rho = []  # Scalar term necessary for Hessian inverse estimation
        
        if not restart:
            msg = 90 * "#" + "\n"
            if self.m_steps is not None:
                msg += "Limited-memory Broyden-Fletcher-Goldfarb-Shanno (L-BFGS) algorithm log file\n"
                msg += "Maximum number of steps to be used for Hessian inverse estimation: %s \n" % self.m_steps
            else:
                msg = "Broyden-Fletcher-Goldfarb-Shanno (BFGS) algorithm log file\n"
            # Printing restart folder
            msg += "Restart folder: %s\n" % self.restart.restart_folder
            msg += 90 * "#" + "\n"
            if verbose:
                print(msg.replace("log file", ""))
            if self.logger:
                self.logger.addToLog(msg)
            
            # Setting internal vectors (model, search direction, and previous gradient vectors)
            prblm_mdl = problem.get_model()
            bfgs_mdl = prblm_mdl.clone()
            bfgs_dmodl = prblm_mdl.clone()
            bfgs_dmodl.zero()
            bfgs_grad0 = bfgs_dmodl.clone()
            
            # Other internal variables
            iiter = 0
        else:
            # Retrieving parameters and vectors to restart the solver
            msg = "Restarting previous solver run from: %s" % self.restart.restart_folder
            if verbose:
                print(msg)
            if self.logger:
                self.logger.addToLog(msg)
            self.restart.read_restart()
            iiter = self.restart.retrieve_parameter("iter")
            self.stepper.alpha = self.restart.retrieve_parameter("alpha")
            initial_obj_value = self.restart.retrieve_parameter("obj_initial")
            bfgs_mdl = self.restart.retrieve_vector("bfgs_mdl")
            bfgs_dmodl = self.restart.retrieve_vector("bfgs_dmodl")
            bfgs_grad0 = self.restart.retrieve_vector("bfgs_grad0")
            # Setting the model and residuals to avoid residual twice computation
            problem.set_model(bfgs_mdl)
            prblm_mdl = problem.get_model()
            # Setting residual vector to avoid its unnecessary computation
            problem.set_residual(self.restart.retrieve_vector("prblm_res"))
            # Retrieving Hessian inverse estimate
            self.rho = self.restart.retrieve_parameter("rho")
            for istep in range(iiter):
                if self.m_steps is not None:
                    if istep < self.m_steps:
                        self.grad_diff_files[istep] = self.restart.retrieve_vector("grad_diff_vectors%s" % istep)
                        self.step_files[istep] = self.restart.retrieve_vector("step_vectors%s.H" % istep)
                else:
                    self.grad_diff_files.append(self.restart.retrieve_vector("grad_diff_vectors%s" % istep))
                    self.step_files.append(self.restart.retrieve_vector("step_vectors%s" % istep))
        
        # Common variables unrelated to restart
        self.tmp_vector = bfgs_dmodl.clone()
        self.tmp_vector.zero()
        prev_mdl = prblm_mdl.clone().zero()
        
        # Inversion loop
        while True:
            # Computing objective function
            obj0 = problem.get_obj(bfgs_mdl)  # Compute objective function value
            prblm_res = problem.get_res(bfgs_mdl)  # Compute residuals
            prblm_grad = problem.get_grad(bfgs_mdl)  # Compute the gradient
            if iiter == 0:
                initial_obj_value = obj0  # For relative objective function value
                # Saving objective function value
                self.restart.save_parameter("obj_initial", initial_obj_value)
                msg = self.iter_msg % (str(iiter).zfill(self.stopper.zfill),
                                       obj0,
                                       problem.get_rnorm(bfgs_mdl),
                                       problem.get_gnorm(bfgs_mdl),
                                       problem.get_fevals(),
                                       problem.get_gevals())
                if verbose:
                    print(msg)
                # Writing on log file
                if self.logger:
                    self.logger.addToLog(msg)
                # Check if either objective function value or gradient norm is NaN
                if isnan(obj0) or isnan(prblm_grad.norm()):
                    raise ValueError("Either gradient norm or objective function value NaN!")
            if prblm_grad.norm() == 0.:
                print("Gradient vanishes identically")
                break
            
            # Saving results
            self.save_results(iiter, problem, force_save=False)
            # Saving current inverted model
            prev_mdl.copy(prblm_mdl)
            
            # Applying approximated Hessian inverse
            msg = "Appplying inverse Hessian estimate"
            if self.m_steps is not None:
                msg += "\nCurrent inverse dot-products of BFGS estimation vectors %s" \
                       % (self.rho[0:min(self.m_steps, iiter)])
            else:
                if len(self.rho) > 0:
                    msg += "\nCurrent inverse dot-products of BFGS estimation vectors %s" % self.rho
            if self.logger:
                self.logger.addToLog(msg)
            self.BFGSMultiply(bfgs_dmodl, prblm_grad, iiter)
            msg = "Done applying inverse Hessian estimate"
            if self.logger:
                self.logger.addToLog(msg)
            
            # grad0 = grad
            bfgs_grad0.copy(prblm_grad)
            # Calling line search
            alpha, success = self.stepper.run(problem, bfgs_mdl, bfgs_dmodl, self.logger)
            if not success:
                msg = "Stepper couldn't find a proper step size, will terminate solver"
                if verbose:
                    print(msg)
                # Writing on log file
                if self.logger:
                    self.logger.addToLog(msg)
                problem.set_model(prev_mdl)
                break
            
            obj1 = problem.get_obj(bfgs_mdl)  # Compute objective function value
            # Redundant test on verifying convergence
            if obj0 <= obj1:
                msg = "Objective function at new point greater or equal than previous one: obj_fun_old=%s obj_fun_new=%s\n" \
                      "Potential issue in the stepper or in revaluation of objective function!" % (obj0, obj1)
                if self.logger:
                    self.logger.addToLog(msg)
                problem.set_model(prev_mdl)
                raise ValueError(msg)
            
            # Compute new gradient
            prblm_grad = problem.get_grad(bfgs_mdl)
            # Compute updates for estimated Hessian inverse
            if self.m_steps is not None:
                # LBFGS
                step_index = iiter % self.m_steps  # Modulo to handle limited memory
                # yn+1=gn+1-gn
                self.grad_diff_vectors[step_index] = bfgs_grad0.clone()
                self.grad_diff_vectors[step_index].scaleAdd(prblm_grad, -1.0, 1.0)
                # sn+1=xn+1-xn = alpha * dmodl
                self.step_vectors[step_index] = bfgs_dmodl.clone()
                self.step_vectors[step_index].scale(alpha)
            else:
                # BFGS
                step_index = iiter
                # yn+1=gn+1-gn
                self.grad_diff_vectors.append(bfgs_grad0.clone())
                self.grad_diff_vectors[step_index].scaleAdd(prblm_grad, -1.0, 1.0)
                # sn+1=xn+1-xn = alpha * dmodl
                self.step_vectors.append(bfgs_dmodl.clone())
                self.step_vectors[step_index].scale(alpha)
            # rhon+1=1/yn+1'sn+1
            denom_dot = self.grad_diff_vectors[step_index].dot(self.step_vectors[step_index])
            # Checking rho
            self.check_rho(denom_dot, step_index, iiter)
            
            # Making first step-length value Hessian guess if not provided by user
            if iiter == 0 and self.H0 is None:
                self.restart.save_parameter("fist_alpha", alpha)
                self.H0 = pyOp.scalingOp(bfgs_dmodl, alpha)
                if self.logger:
                    self.logger.addToLog("First step-length value used as first Hessian inverse estimate!")
                self.stepper.alpha = 1.0
            
            # Increasing iteration counter
            iiter = iiter + 1
            
            # Using alpha = 1.0 after first iteration
            if iiter != 0 and not self.save_alpha:
                self.stepper.alpha = 1.0
            
            # Saving current model and previous search direction in case of restart
            self.restart.save_parameter("iter", iiter)
            self.restart.save_parameter("alpha", alpha)
            self.restart.save_vector("bfgs_mdl", bfgs_mdl)
            self.restart.save_vector("bfgs_dmodl", bfgs_dmodl)
            self.restart.save_vector("bfgs_grad0", bfgs_grad0)
            # Saving Inverse Hessian estimate for restart
            self.restart.save_parameter("rho", self.rho)
            self.restart.save_vector("grad_diff_vectors%s" % step_index, self.grad_diff_vectors[step_index])
            self.restart.save_vector("step_vectors%s" % step_index, self.step_vectors[step_index])
            # Saving data space vectors
            self.restart.save_vector("prblm_res", prblm_res)
            
            # iteration info
            msg = self.iter_msg % (str(iiter).zfill(self.stopper.zfill),
                                   obj1,
                                   problem.get_rnorm(bfgs_mdl),
                                   problem.get_gnorm(bfgs_mdl),
                                   problem.get_fevals(),
                                   problem.get_gevals())
            if verbose:
                print(msg)
            # Writing on log file
            if self.logger:
                self.logger.addToLog("\n" + msg)
            # Check if either objective function value or gradient norm is NaN
            if isnan(obj1) or isnan(prblm_grad.norm()):
                raise ValueError("Either gradient norm or objective function value NaN!")
            if self.stopper.run(problem, iiter, initial_obj_value, verbose):
                break
        
        # Writing last inverted model
        self.save_results(iiter, problem, force_save=True, force_write=True)
        msg = 90 * "#" + "\n"
        if self.m_steps is not None:
            msg += "Limited-memory Broyden-Fletcher-Goldfarb-Shanno (L-BFGS) algorithm log file end\n"
        else:
            msg += "Broyden-Fletcher-Goldfarb-Shanno (BFGS) algorithm log file end\n"
        msg += 90 * "#" + "\n"
        if verbose:
            print(msg.replace("log file ", ""))
        if self.logger:
            self.logger.addToLog(msg)
        self.restart.clear_restart()
        # Resetting inverse Hessian matrix
        self.H0 = None
        del self.tmp_vector
        self.tmp_vector = None
