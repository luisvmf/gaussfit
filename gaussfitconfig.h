//
//This is the configuration file for gaussfit parameters.
//
//

    /**
       Print all iterations to the console if 1;
    */
#define Debugiter 0
    /**
     * The line search algorithm.
     *  This parameter specifies a line search algorithm to be used by the
     *  L-BFGS routine.
     * Possible values are:
     * LBFGS_LINESEARCH_BACKTRACKING;
     * LBFGS_LINESEARCH_MORETHUENTE;
     * LBFGS_LINESEARCH_BACKTRACKING_ARMIJO;
     * LBFGS_LINESEARCH_BACKTRACKING_WOLFE;
     * LBFGS_LINESEARCH_BACKTRACKING_STRONG_WOLFE
	 */
#define gaussfit_linesearch LBFGS_LINESEARCH_BACKTRACKING_STRONG_WOLFE;
    /**
     * A parameter to control the accuracy of the line search routine.
     *  The default value is \c 0.9. If the function and gradient
     *  evaluations are inexpensive with respect to the cost of the
     *  iteration (which is sometimes the case when solving very large
     *  problems) it may be advantageous to set this parameter to a small
     *  value. A typical small value is \c 0.1. This parameter should be
     *  greater than the \ref ftol parameter (\c 1e-4) and smaller than
     *  \c 1.0.
     */
#define gaussfit_gtol 0.1;
    /**
     * A parameter to control the accuracy of the line search routine.
     *  The default value is \c 1e-4. This parameter should be greater
     *  than zero and smaller than \c 0.5.
     */
#define gaussfit_ftol 2.220446049250313e-9;
    /**
     * The number of corrections to approximate the inverse hessian matrix.
     *  The L-BFGS routine stores the computation results of previous \ref m
     *  iterations to approximate the inverse hessian matrix of the current
     *  iteration. This parameter controls the size of the limited memories
     *  (corrections). The default value is \c 6. Values less than \c 3 are
     *  not recommended. Large values will result in excessive computing time.
     */
#define gaussfit_m 20;
    /**
     * Distance for delta-based convergence test.
     *  This parameter determines the distance, in iterations, to compute
     *  the rate of decrease of the objective function. If the value of this
     *  parameter is zero, the library does not perform the delta-based
     *  convergence test. The default value is \c 0.
     */
	//XXX Changing gaussfit_past to 1 and increasing gaussfit_delta makes L-BFGS converge but the fit result is worse.
#define gaussfit_past 0;
   /**
     * Delta for convergence test.
     *  This parameter determines the minimum rate of decrease of the
     *  objective function. The library stops iterations when the
     *  following condition is met:
     *      (f' - f) / f < \ref delta,
     *  where f' is the objective value of \ref past iterations ago, and f is
     *  the objective value of the current iteration.
     *  The default value is \c 1e-5.
     */
#define gaussfit_delta 1e-5;
    /**
     * Epsilon for convergence test.
     *  This parameter determines the accuracy with which the solution is to
     *  be found. A minimization terminates when
     *      ||g|| < \ref epsilon * max(1, ||x||),
     *  where ||.|| denotes the Euclidean (L2) norm. The default value is
     *  \c 1e-5.
     */
#define gaussfit_epsilon 1e-8;
    /**
     * The maximum number of trials for the line search.
     *  This parameter controls the number of function and gradients evaluations
     *  per iteration for the line search routine. The default value is \c 40.
     */
#define gaussfit_max_linesearch 350;
    /**
     * The maximum number of iterations.
     *  The lbfgs() function terminates an optimization process with
     *  ::LBFGSERR_MAXIMUMITERATION status code when the iteration count
     *  exceedes this parameter. Setting this parameter to zero continues an
     *  optimization process until a convergence or error. The default value
     *  is \c 0.
     */
#define gaussfit_max_iterations 0;
    /**
     * The minimum step of the line search routine.
     *  The default value is \c 1e-20. This value need not be modified unless
     *  the exponents are too large for the machine being used, or unless the
     *  problem is extremely badly scaled (in which case the exponents should
     *  be increased).
     */
#define gaussfit_min_step 1e-20;
    /**
     * The maximum step of the line search.
     *  The default value is \c 1e+20. This value need not be modified unless
     *  the exponents are too large for the machine being used, or unless the
     *  problem is extremely badly scaled (in which case the exponents should
     *  be increased).
     */
#define gaussfit_max_step 1e+20;
    /**
     * A coefficient for the Wolfe condition.
     *  This parameter is valid only when the backtracking line-search
     *  algorithm is used with the Wolfe condition,
     *  ::LBFGS_LINESEARCH_BACKTRACKING_STRONG_WOLFE or
     *  ::LBFGS_LINESEARCH_BACKTRACKING_WOLFE .
     *  The default value is \c 0.9. This parameter should be greater
     *  the \ref ftol parameter and smaller than \c 1.0.
     */
#define gaussfit_wolfe 0.99;
    /**
     * The machine precision for floating-point values.
     *  This parameter must be a positive value set by a client program to
     *  estimate the machine precision. The line search routine will terminate
     *  with the status code (::LBFGSERR_ROUNDING_ERROR) if the relative width
     *  of the interval of uncertainty is less than this parameter.
     */
#define gaussfit_xtol 1e-20;
