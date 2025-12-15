# Penalized objective function
penalized_obj <- function(
    x,
    obj_fn,
    w,
    pen_fn,
    pen_par_id = NULL,
    pen_diff_id = NULL
) {
    out <- obj_fn(x)
    if (!is.null(pen_par_id)) {
        out <- out + w * sum(pen_fn(x[pen_par_id]))
    }
    if (!is.null(pen_diff_id)) {
        trans_diff <- rep(list(identity), length(pen_diff_id))
        if (any(grepl("^loading", names(pen_diff_id)))) {
            trans_diff[[grep("^loading", names(pen_diff_id))]] <- log
        }
        pen_diff <- Map(
            function(mat, trans) {
                x_mat <- matrix(
                    x[mat],
                    nrow = nrow(mat),
                    ncol = ncol(mat)
                )
                composite_pair_loss(x_mat, fun = pen_fn, trans = trans)
            },
            mat = pen_diff_id,
            trans = trans_diff
        )
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
#' \dontrun{
#' library(lavaan)
#'
#' # Fit a longitudinal factor model using PoliticalDemocracy data
#' ind_mat <- cbind(c("y1", "y2", "y3", "y4"), c("y5", "y6", "y7", "y8"))
#' fit <- longcfa(ind_mat, lv_names = c("dem60", "dem65"), data = PoliticalDemocracy,
#'                long_equal = c("loadings", "intercepts"), lag_cov = TRUE)
#' # Obtain an unidentified model
#' mod_un <- longcfa_syntax(
#'     ind_mat, lv_names = c("dem60", "dem65"),
#'     lag_cov = TRUE,
#'     free_latvars = TRUE, free_latmeans = TRUE
#' )
#' fit_un <- cfa(mod_un, data = PoliticalDemocracy, do.fit = FALSE, std.lv = TRUE,
#'               start = fit)
#'
#' # Get parameter IDs for loadings
#' load_ids <- get_lav_par_id(fit_un, op = "=~", ind_matrix = ind_mat)
#' int_ids <- get_lav_par_id(fit_un, op = "~1", ind_matrix = ind_mat)
#'
#' # Apply penalized estimation with alignment loss
#' pen_fit <- penalized_est(
#'     x = fit_un,
#'     w = 0.1,
#'     pen_diff_id = list(cbind(t(load_ids), t(int_ids))),
#'     pen_fn = "alf"
#' )
#'
#' # Compare parameter estimates
#' cbind(coef(fit), coef(pen_fit))
#'
#' # Compare log-likelihoods
#' c("scalar invariance" = logLik(fit), "penalized" = logLik(pen_fit))
#' }
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
    f1 <- function(v) {
        penalized_obj(
            v,
            obj_fn = function(pars) {
                ff$objective_function(pars, lavaan_model = x)
            },
            w = w,
            pen_fn = pen_fn,
            pen_par_id = pen_par_id,
            pen_diff_id = pen_diff_id
        )
    }
    gr1 <- function(v) {
        penalized_gr(
            v,
            gr_fn = function(pars) ff$gradient_function(pars, lavaan_model = x),
            w = w,
            pen_gr = pen_gr,
            pen_par_id = pen_par_id,
            pen_diff_id = pen_diff_id
        )
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
    hess <- numDeriv::hessian(f1, opt$par)
    if (se == "robust.huber.white") {
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
    vc_out <- try(solve(hess) %*% meat %*% solve(hess), silent = TRUE)
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

penalized_gr <- function(
    x,
    gr_fn,
    w,
    pen_gr,
    pen_par_id = NULL,
    pen_diff_id = NULL,
    ...
) {
    out <- gr_fn(x)
    if (!is.null(pen_par_id)) {
        out <- out + w * hot_gr(x, pen_par_id, pen_gr, ...)
    }
    if (!is.null(pen_diff_id)) {
        trans_diff <- rep(list(identity), length(pen_diff_id))
        gr_trans_diff <- rep(list(function(x) 1), length(pen_diff_id))
        if (any(grepl("^loading", names(pen_diff_id)))) {
            trans_diff[[grep("^loading", names(pen_diff_id))]] <- log
            gr_trans_diff[[grep("^loading", names(pen_diff_id))]] <-
                function(x) 1 / x
        }
        pen_diff_gr <- Map(
            function(mat, trans, gr_trans) {
                hot_gr(
                    x,
                    mat,
                    gr_cpl,
                    gr_fun = pen_gr,
                    trans = trans,
                    gr_trans = gr_trans,
                    ...
                )
            },
            mat = pen_diff_id,
            trans = trans_diff,
            gr_trans = gr_trans_diff
        )
        out <- out + w * Reduce(`+`, pen_diff_gr)
    }
    out
}
