# Penalized objective function
penalized_obj <- function(x, obj_fn, w, pen_fn, pen_par_id, diff_configs) {
    out <- obj_fn(x)
    
    if (!is.null(pen_par_id)) {
        out <- out + w * sum(pen_fn(x[pen_par_id]))
    }
    
    if (!is.null(diff_configs)) {
        pen_diff <- lapply(diff_configs, function(cfg) {
            x_mat <- matrix(x[cfg$mat], nrow = nrow(cfg$mat), ncol = ncol(cfg$mat))
            x_trans <- as.matrix(cfg$trans(x_mat))
            
            # Use pre-computed combn_idx and rescale_val
            diffs <- x_trans[cfg$combn_idx[1, ], , drop = FALSE] - 
                     x_trans[cfg$combn_idx[2, ], , drop = FALSE]
            
            sum(pen_fn(diffs), na.rm = TRUE) * cfg$rescale_val
        })
        out <- out + w * sum(unlist(pen_diff))
    }
    out
}

#' Penalized Parameter Estimation for Longitudinal CFA Models
#'
#' Performs penalized estimation on a lavaan model object by optimizing a
#' penalized objective function. The function extracts the objective function
#' from a lavaan model, applies a penalty function to specified parameters
#' or pairwise differences of parameters, and returns an updated model with
#' the optimized parameter estimates.
#'
#' @param x A fitted lavaan model object from which estimation components will
#'   be extracted.
#' @param w Numeric scalar. Penalty weight (multiplier) applied to the penalty
#'   terms.
#' @param pen_par_id Integer vector of parameter IDs to apply the penalty function
#'   directly to, in the same order as returned by `lavaan::coef()` and by
#'   [lavaan::partable()], with only the free elements.
#' @param pen_diff_id List of matrices containing parameter IDs. For each matrix,
#'   the penalty is applied to the pairwise differences of parameters in the same
#'   column indicated by the IDs. For matrices with names starting with "loading",
#'   the log transformation is applied before computing differences.
#' @param pen_fn A character string (`"l0a"` or `"alf"`) or a function that computes
#'   the penalty. Default is `"l0a"`.
#' @param pen_gr A function that computes the gradient of the penalty function.
#'   If `pen_fn` is `"l0a"` or `"alf"`, this is automatically set.
#' @param se Character string specifying the type of standard errors to compute.
#'   Options are `"none"` (default; no standard errors) or `"robust.huber.white"`
#'   (robust sandwich estimator using numerical Hessian and first-order information,
#'   which is the same as used in the `"mlr"` estimator).
#' @param opt_control A list of control parameters passed to [stats::nlminb()].
#'   Default includes `eval.max = 2e4`, `iter.max = 1e4`, and `abs.tol = 1e-20`.
#'
#' @section Warning:
#' The returned object is not fitted using standard ML. Standard errors reported
#' by `summary()` or `parameterEstimates()` will be missing unless
#' `se = "robust.huber.white"` was specified. Even then, they are based on an
#' experimental sandwich approximation and should be interpreted with caution.
#'
#' @return A lavaan model object updated with the penalized parameter estimates.
#'   The returned object includes an attribute `opt_info` containing the
#'   optimization information returned by `nlminb()`.
#'
#' @details
#' The function uses `nlminb()` to minimize a penalized objective function that
#' combines the standard lavaan objective function with a penalty term. Only the
#' parameter estimates and the log-likelihood should be interpreted. The
#' returned object was not "fitted" (`do.fit = FALSE`) to avoid users
#' interpreting the standard errors, which are generally not valid with
#' penalized estimation. The degrees of freedom may also be inaccurate. If the
#' optimization does not converge (convergence code != 0), a warning is issued.
#'
#' @seealso \code{\link[lavaan]{lavaan}}, \code{\link[stats]{nlminb}}
#'
#' @importFrom stats nlminb
#' @examples
#' library(lavaan)
#'
#' # Define a longitudinal factor model with PoliticalDemocracy data
#' model <- "
#'   dem60 =~ y1 + y2 + y3 + y4
#'   dem65 =~ y5 + y6 + y7 + y8
#'   dem60 ~~ dem65
#'   dem60 ~~ 1 * dem60
#'   dem65 ~~ NA * dem65
#'   dem60 ~ 0
#'   dem65 ~ NA * 1
#'   y1 ~~ y5
#'   y2 ~~ y6
#'   y3 ~~ y7
#'   y4 ~~ y8
#' "
#'
#' # Fit the model without constraints first to get parameter table
#' fit_un <- cfa(model, data = PoliticalDemocracy, std.lv = TRUE,
#'               meanstructure = TRUE, do.fit = FALSE)
#'
#' # Get parameter IDs
#' pt <- parTable(fit_un)
#' # Loadings
#' load_60 <- pt$free[pt$op == "=~" & pt$lhs == "dem60"]
#' load_65 <- pt$free[pt$op == "=~" & pt$lhs == "dem65"]
#' # Intercepts
#' int_60 <- pt$free[pt$op == "~1" & pt$lhs %in% c("y1", "y2", "y3", "y4")]
#' int_65 <- pt$free[pt$op == "~1" & pt$lhs %in% c("y5", "y6", "y7", "y8")]
#'
#' # Apply penalized estimation to penalize differences in loadings and intercepts
#' pen_fit <- penalized_est(
#'     x = fit_un,
#'     w = 0.03,
#'     pen_diff_id = list(
#'         loadings = rbind(load_60, load_65),
#'         intercepts = rbind(int_60, int_65)
#'     ),
#'     pen_fn = "l0a"
#' )
#'
#' # Compare parameter estimates
#' summary(pen_fit)
#'
#' @importFrom stats update
#' @importFrom utils modifyList
#' @export
penalized_est <- function(
    x,
    w,
    pen_par_id = NULL,
    pen_diff_id = NULL,
    pen_fn = "l0a",
    pen_gr = NULL,
    se = "none",
    opt_control = list()
) {
    # Define default control parameters
    control_defaults <- list(
        eval.max = 2e4,
        iter.max = 1e4,
        abs.tol = 1e-20
    )

    # Merge with user input
    control <- modifyList(control_defaults, opt_control)

    ff <- lavaan::lav_export_estimation(x)
    if (pen_fn %in% c("l0a", "alf")) {
        pen_gr <- switch(
            pen_fn,
            l0a = gr_l0a,
            alf = gr_alf
        )
        pen_fn <- get(pen_fn)
    }

    if (!is.function(pen_fn) && !pen_fn %in% c("l0a", "alf")) {
        stop("pen_fn must be 'l0a', 'alf', or a function.")
    }
    diff_configs <- NULL
    if (!is.null(pen_diff_id)) {
        diff_configs <- lapply(names(pen_diff_id), function(nm) {
            mat <- pen_diff_id[[nm]]
            is_loading <- grepl("^loading", nm)
            
            # Pre-assign transformations
            trans <- if (is_loading) log else identity
            gr_trans <- if (is_loading) function(x) 1/x else function(x) rep(1, length(x))
            
            # Pre-compute combinatorics and rescaling
            nrow_x <- nrow(mat)
            combn_idx <- combn(nrow_x, 2)
            rescale_val <- (nrow_x - 1) / ncol(combn_idx)
            
            # Pre-compute gradient row indices to avoid which() loops later
            grad_idx <- lapply(seq_len(nrow_x), function(i) {
                list(
                    idx1 = which(combn_idx[1, ] == i),
                    idx2 = which(combn_idx[2, ] == i)
                )
            })
            
            list(
                mat = mat, 
                trans = trans, 
                gr_trans = gr_trans,
                combn_idx = combn_idx, 
                rescale_val = rescale_val,
                grad_idx = grad_idx
            )
        })
    }
    f1 <- function(v) {
        penalized_obj(
            v,
            obj_fn = function(pars) {
                ff$objective_function(pars, lavaan_model = x)
            },
            w = w,
            pen_fn = pen_fn,
            pen_par_id = pen_par_id,
            diff_configs = diff_configs
        )
    }
    gr1 <- if (!is.null(pen_gr)) {
        function(v) {
            penalized_gr(
                v,
                gr_fn = function(pars) ff$gradient_function(pars, lavaan_model = x),
                w = w,
                pen_gr = pen_gr,
                pen_par_id = pen_par_id,
                diff_configs = diff_configs
            )
        }
    } else {
        NULL  # Let nlminb compute numerical gradient
    }
    opt <- nlminb(
        ff$starting_values,
        objective = f1,
        gradient = gr1,
        control = control
    )
    if (opt$convergence != 0) {
        warning(
            "Optimization did not converge. Try using better starting values, ",
            "or adjusting optimization control parameters."
        )
    }
    x_opt <- x@Options
    x_opt$start <- opt$par
    x_opt$do.fit <- FALSE
    x_opt$se <- "none"
    out <- lavaan::lavaan(
        lavaan::partable(x),
        slotOptions = x_opt,
        slotSampleStats = x@SampleStats,
        slotData = x@Data
        # do.fit = FALSE,
        # start = opt$par
    )
    out <- add_nlminb_info(out, opt)
    if (!se %in% c("none", "robust.huber.white")) {
        warning("se must be either 'none' or 'robust.huber.white'. ",
                "Defaulting to 'none'")
    }
    if (se == "robust.huber.white") {
        hess <- numDeriv::hessian(f1, opt$par)
        attr(out, "hessian") <- hess
        out <- add_vcov_pen(out, hess)
    }
    out
}

add_nlminb_info <- function(fit, opt) {
    fit@optim$x <- opt$par
    fit@optim$fx <- opt$objective
    fit@optim$iterations <- opt$iterations
    fit@optim$converged <- as.logical(1 - opt$convergence)
    fit@optim$control <- opt$control
    fit@optim$dx <- opt$gradient
    fit@optim$npar <- length(opt$par)
    fit
}

#' @importFrom lavaan lavInspect
add_vcov_pen <- function(fit, hess) {
    meat <- lavInspect(fit, "information.first.order")
    H_inv <- solve(hess)
    vc_out <- try(H_inv %*% meat %*% H_inv, silent = TRUE)
    if (inherits(vc_out, "try-error")) {
        vc_out <- NULL
    }
    fit@vcov$se <- "robust.huber.white"
    fit@vcov$vcov <- vc_out / lavInspect(fit, "nobs")
    fit@vcov$information <- "observed"
    fit@Options$se <- "robust.huber.white"
    fit@Options$information <- rep("observed", 2)
    fit@ParTable$se <- 0 * fit@ParTable$est
    fit@ParTable$se[which(fit@ParTable$free > 0)] <- sqrt(diag(
        fit@vcov$vcov
    ))
    fit
}

penalized_gr <- function(x, gr_fn, w, pen_gr, pen_par_id, diff_configs, ...) {
    out <- gr_fn(x)
    
    if (!is.null(pen_par_id)) {
        out <- out + w * hot_gr(x, pen_par_id, pen_gr, ...)
    }
    
    if (!is.null(diff_configs)) {
        pen_diff_gr <- lapply(diff_configs, function(cfg) {
            x_mat <- as.matrix(cfg$trans(matrix(x[cfg$mat], nrow = nrow(cfg$mat))))
            
            diffs <- x_mat[cfg$combn_idx[1, ], , drop = FALSE] - 
                     x_mat[cfg$combn_idx[2, ], , drop = FALSE]
                     
            grad_contribs <- pen_gr(diffs)
            grad <- matrix(0, nrow = nrow(x_mat), ncol = ncol(x_mat))
            
            # Loop is now incredibly fast because indices are pre-calculated
            for (i in seq_len(nrow(x_mat))) {
                g1 <- grad_contribs[cfg$grad_idx[[i]]$idx1, , drop = FALSE]
                g2 <- grad_contribs[cfg$grad_idx[[i]]$idx2, , drop = FALSE]
                
                grad[i, ] <- colSums(g1, na.rm = TRUE) - colSums(g2, na.rm = TRUE)
            }
            
            grad[which(is.na(x_mat))] <- NA
            grad_vec <- as.vector(grad) * cfg$rescale_val * cfg$gr_trans(x[cfg$mat])
            
            # Re-map back to the full parameter vector space (hot_gr equivalent inline)
            full_grad <- 0 * x
            full_grad[na.omit(as.numeric(cfg$mat))] <- na.omit(grad_vec)
            full_grad
        })
        out <- out + w * Reduce(`+`, pen_diff_gr)
    }
    out
}
