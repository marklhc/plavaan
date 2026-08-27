#' Effective Number of Parameters and Degrees of Freedom for Penalized Fits
#'
#' Computes the effective number of model parameters, and hence the
#' effective model degrees of freedom, of a penalized lavaan fit using the
#' soft-count property of the l0a penalty: a directly penalized parameter
#' with estimate `z` is counted as `l0a(z, eps) = z^2 / (z^2 + eps)`
#' effective parameters, which is close to 1 for large `|z|` and close to 0
#' for small `|z|`. For a difference penalty, each column contributes one
#' shared value (the invariance baseline) plus the l0a soft count of the
#' pairwise differences in that column.
#'
#' @param x A lavaan model object, typically the return value of
#'   [penalized_est()]. For such objects, the penalty specification is read
#'   from the object's \sQuote{penalized} attribute.
#' @param pen_par_id Integer vector of directly penalized free-parameter
#'   IDs, in the same order as returned by `lavaan::coef()`. For a fit from
#'   [penalized_est()], defaults to the value recorded in the object.
#' @param pen_diff_id Named list of integer matrices of free-parameter IDs
#'   for difference penalties (same structure as in [penalized_est()]). For
#'   a fit from [penalized_est()], defaults to the value recorded in the
#'   object.
#' @param pen_fn The penalty function: `"l0a"`, `"alf"`, or a function.
#'   For a fit from [penalized_est()], defaults to the value recorded in
#'   the object.
#' @param eps The epsilon used by the built-in penalties. For a fit from
#'   [penalized_est()], defaults to the value recorded in the object.
#' @param ... Additional arguments passed to a user-supplied `pen_fn`.
#'
#' @return A data frame of class `plavaan_efdf` with one row per penalty
#'   component (`"direct penalty"`, and one row per named block of
#'   `pen_diff_id`) and a final `"TOTAL"` row; the row names are the
#'   component labels. Columns:
#' \describe{
#'   \item{`npar`}{Number of (nominal) free parameters covered by the row.
#'     For the `TOTAL` row, `length(coef(x))`.}
#'   \item{`npar_effective`}{Effective number of parameters, with penalized
#'     parameters soft-counted by their penalty values.}
#'   \item{`df_saved`}{`npar - npar_effective`, the degrees of freedom saved
#'     by the penalty relative to the nominal model.}
#' }
#' The data frame also carries an `info` attribute with `n_stats` (number of
#' sample statistics), `df_model` (nominal model df, possibly negative),
#' `df_model_effective`, `w`, `eps`, and `pen_fn`.
#'
#' @details
#' For a fit returned by [penalized_est()], any argument left at its
#' default (`NULL` for `pen_par_id`/`pen_diff_id`, `"l0a"` for `pen_fn`,
#' `.01` for `eps`) is filled in from the object's \sQuote{penalized}
#' attribute; explicitly supplied values take precedence. Note that a value
#' equal to the default cannot be distinguished from an omitted one, so
#' passing, e.g., `eps = .01` explicitly will still be overridden by the
#' recorded value.
#'
#' An ID that appears in both `pen_par_id` and a `pen_diff_id` matrix is a
#' hard error, as the parameter would be double-counted.
#'
#' The effective-df interpretation is calibrated for the l0a penalty, which
#' takes the value 0 at the origin. With `pen_fn = "alf"`,
#' `alf(0) = eps^0.25 != 0`, so the soft counts do not reach 0 and the
#' effective df should be interpreted with care. A custom `pen_fn` is
#' evaluated as if it had l0a-like 0/1 behavior; the result is approximate.
#'
#' The nominal model df, `n_stats - length(coef(x))`, may be negative for
#' penalized (often under-identified) models; the effective model df,
#' `n_stats - npar_effective`, is the meaningful quantity. The number of
#' sample statistics `n_stats` is counted per group from the model's
#' parameter table.
#'
#' Fit indices based on the effective df (chi-square, CFI, RMSEA, AIC,
#' BIC) are obtained with `fitmeasures()` on a fit from
#' [penalized_est()], which freezes the parameters at the penalized
#' estimates and refits. This fit evaluation is experimental and requires
#' the fit to have been created with a non-\code{"none"} \code{test}
#' (e.g. `test = "Chisq"`); with the default `test = "none"`,
#' `fitmeasures()` is unavailable.
#'
#' @examples
#' library(lavaan)
#'
#' model <- "
#'   dem60 =~ y1 + y2 + y3 + y4
#'   dem65 =~ y5 + y6 + y7 + y8
#'   dem60 ~~ dem65
#' "
#' fit <- cfa(model, data = PoliticalDemocracy, std.lv = TRUE, do.fit = FALSE)
#' pt <- parTable(fit)
#' load_60 <- pt$free[pt$op == "=~" & pt$lhs == "dem60"]
#' load_65 <- pt$free[pt$op == "=~" & pt$lhs == "dem65"]
#'
#' pen <- penalized_est(
#'   fit, w = 0.03,
#'   pen_diff_id = list(loadings = rbind(load_60, load_65)),
#'   eps = .01,
#'   test = "Chisq"
#' )
#' effective_df(pen)
#' fitmeasures(pen, c("chisq", "df", "cfi"))
#'
#' @export
effective_df <- function(x, pen_par_id = NULL, pen_diff_id = NULL,
                         pen_fn = "l0a", eps = .01, ...) {
    spec <- plavaan_penalty_spec(x, pen_par_id, pen_diff_id, pen_fn, eps, ...)
    est <- plavaan_estimates(x)
    n_stats <- plavaan_n_stats(x)
    npar_total <- length(est)
    npar_eff <- plavaan_npar_eff(x, spec)
    penalty_fn <- plavaan_penalty_fn(spec)

    mk_row <- function(label, npar, npar_effective) {
        row <- data.frame(
            npar = npar,
            npar_effective = npar_effective,
            df_saved = npar - npar_effective,
            check.names = FALSE
        )
        rownames(row) <- label
        row
    }
    rows <- list()
    if (!is.null(spec$pen_par_id)) {
        direct <- spec$pen_par_id
        npar_effective_direct <- sum(penalty_fn(est[direct]))
        rows[["direct penalty"]] <-
            mk_row("direct penalty", length(direct), npar_effective_direct)
    }
    if (!is.null(spec$pen_diff_id)) {
        for (nm in names(spec$pen_diff_id)) {
            block <- spec$pen_diff_id[[nm]]
            k <- ncol(block)
            loss <- 0
            for (j in seq_len(k)) {
                loss <- loss +
                    composite_pair_loss(est[block[, j]], fun = penalty_fn)
            }
            npar_block <- length(unique(block[!is.na(block)]))
            rows[[nm]] <- mk_row(nm, npar_block, k + loss)
        }
    }
    rows[["TOTAL"]] <- mk_row("TOTAL", npar_total, npar_eff)
    out <- do.call(rbind, rows)
    class(out) <- c("plavaan_efdf", "data.frame")
    attr(out, "info") <- list(
        n_stats = n_stats,
        df_model = n_stats - npar_total,
        df_model_effective = n_stats - npar_eff,
        w = spec$w,
        eps = spec$eps,
        pen_fn = spec$pen_fn
    )
    out
}

#' Print Method for Effective Degrees of Freedom Tables
#'
#' Prints the component table from [effective_df()] followed by summary
#' lines with the number of sample statistics, the nominal and effective
#' model degrees of freedom, and the penalty settings.
#'
#' @param x An object of class `plavaan_efdf`, typically the return value
#'   of [effective_df()].
#' @param ... Passed to \code{print.data.frame()}.
#' @export
print.plavaan_efdf <- function(x, ...) {
    info <- attr(x, "info")
    print(as.data.frame(x), ...)
    pen_fn_label <- if (is.function(info$pen_fn)) "custom" else info$pen_fn
    w_label <- if (length(info$w) == 1 && is.na(info$w)) {
        "unknown (not a penalized_est() fit)"
    } else {
        info$w
    }
    df_note <- if (info$df_model < 0) {
        " (negative: the nominal model is under-identified; the effective df is the meaningful quantity)"
    } else {
        ""
    }
    cat(
        "\nn_stats (sample moments):  ", info$n_stats, "\n",
        "nominal model df:  ", info$df_model, df_note, "\n",
        "effective model df:  ", round(info$df_model_effective, 2), "\n",
        "penalty:  ", pen_fn_label, " (w = ", w_label,
        ", eps = ", info$eps, ")\n",
        sep = ""
    )
}

# Resolve the penalty specification from a fit's "penalized" attribute and
# (user-supplied) arguments, then validate it. Arguments left at their
# defaults (NULL / "l0a" / .01) are filled in from the attribute; explicit
# user values take precedence.
plavaan_penalty_spec <- function(x, pen_par_id = NULL, pen_diff_id = NULL,
                                 pen_fn = "l0a", eps = .01, ...) {
    dots <- list(...)
    attrs <- attr(x, "penalized")
    test <- "none"
    if (!is.null(attrs)) {
        if (is.null(pen_par_id)) pen_par_id <- attrs$pen_par_id
        if (is.null(pen_diff_id)) pen_diff_id <- attrs$pen_diff_id
        if (identical(pen_fn, "l0a")) pen_fn <- attrs$pen_fn
        if (identical(eps, 0.01)) eps <- attrs$eps
        w <- attrs$w
        # Default to "none" for objects created before `test` was recorded.
        if (!is.null(attrs$test)) test <- attrs$test
    } else {
        if (is.null(pen_par_id) && is.null(pen_diff_id)) {
            stop(
                "No penalty specification found: pass pen_par_id/pen_diff_id ",
                "or use a fit from penalized_est().",
                call. = FALSE
            )
        }
        w <- NA_real_
    }

    if (
        !is.function(pen_fn) &&
            !(
                is.character(pen_fn) &&
                    length(pen_fn) == 1 &&
                    !is.na(pen_fn) &&
                    pen_fn %in% c("l0a", "alf")
            )
    ) {
        stop("pen_fn must be 'l0a', 'alf', or a function.")
    }
    if (!is.numeric(eps) || length(eps) != 1 || !is.finite(eps) || eps <= 0) {
        stop("eps must be a positive numeric scalar.")
    }
    if (
        !is.null(pen_diff_id) &&
            (
                !is.list(pen_diff_id) ||
                    !all(vapply(pen_diff_id, is.matrix, logical(1)))
            )
    ) {
        stop("pen_diff_id must be a list of matrices.")
    }

    # Validate IDs the same way as random_start().
    n_par <- length(plavaan_estimates(x))
    pen_ids <- c(pen_par_id, unlist(pen_diff_id, use.names = FALSE))
    pen_ids <- pen_ids[!is.na(pen_ids)]
    if (length(pen_ids) > 0) {
        if (
            !is.numeric(pen_ids) || any(!is.finite(pen_ids)) ||
                any(pen_ids <= 0) || any(pen_ids != as.integer(pen_ids)) ||
                any(pen_ids > n_par)
        ) {
            stop("Penalty parameter IDs must be positive integer free-parameter indices.")
        }
    }
    if (!is.null(pen_par_id) && !is.null(pen_diff_id)) {
        overlap <- intersect(
            pen_par_id, unlist(pen_diff_id, use.names = FALSE)
        )
        overlap <- overlap[!is.na(overlap)]
        if (length(overlap) > 0) {
            stop(
                "Parameter ID(s) ",
                paste(overlap[seq_len(min(5L, length(overlap)))], collapse = ", "),
                if (length(overlap) > 5) " and more" else "",
                " appear in both pen_par_id and pen_diff_id; each parameter ",
                "may be penalized once.",
                call. = FALSE
            )
        }
    }

    list(
        w = w,
        pen_par_id = pen_par_id,
        pen_diff_id = pen_diff_id,
        pen_fn = pen_fn,
        eps = eps,
        test = test,
        dots = dots
    )
}

# Penalty function (as a closure) matching the resolved spec. Built-in
# penalties are bound to the spec's eps and ignore extra arguments, exactly
# as in penalized_est(); a custom function receives the spec's dots.
plavaan_penalty_fn <- function(spec) {
    if (is.function(spec$pen_fn)) {
        dots <- spec$dots
        function(z, ...) do.call(spec$pen_fn, c(list(z), dots))
    } else if (identical(spec$pen_fn, "alf")) {
        function(z, ...) alf(z, eps = spec$eps)
    } else {
        function(z, ...) l0a(z, eps = spec$eps)
    }
}

# Free-parameter estimates, in the same order as `lavaan::coef()` (i.e. the
# row order of parTable() filtered to free rows). parTable() is used rather
# than the coef() S3 generic because S3 dispatch from within this package's
# namespace cannot see lavaan's coef.lavaan method (it is registered in
# lavaan's namespace, not plavaan's).
plavaan_estimates <- function(x) {
    pt <- lavaan::parTable(x)
    as.numeric(pt$est[pt$free > 0])
}

# Effective number of parameters of the whole model: free parameters
# outside the penalty count in full; directly penalized parameters are
# soft-counted by their penalty values; each difference-penalty block
# contributes one shared value per column plus the soft count of the
# pairwise differences in each column. Summing the per-column
# composite_pair_loss() calls equals the full-matrix loss used in the
# optimization and is safe for structural NAs (NA differences are dropped
# inside composite_pair_loss()).
plavaan_npar_eff <- function(x, spec) {
    est <- plavaan_estimates(x)
    p <- length(est)
    direct <- spec$pen_par_id
    blocks <- spec$pen_diff_id
    P <- unique(c(direct, unlist(blocks, use.names = FALSE)))
    P <- P[!is.na(P)]
    penalty_fn <- plavaan_penalty_fn(spec)
    out <- p - length(P)
    if (length(direct)) {
        out <- out + sum(penalty_fn(est[direct]))
    }
    for (b in blocks) {
        k <- ncol(b)
        out <- out + k
        for (j in seq_len(k)) {
            out <- out + composite_pair_loss(est[b[, j]], fun = penalty_fn)
        }
    }
    out
}

# Number of sample statistics of the USER's fit. This is the quantity the
# effective degrees of freedom are measured against (df = n_stats -
# npar_effective), so it must equal the number of sample moments lavaan
# itself uses in the model test.
#
# Primary: lavaan::lav_pt_ndat(), lavaan's own count (exported in recent
# lavaan). It handles ordinal thresholds, composites, the correlation
# metric, multilevel, and group models, and by construction always matches
# lavaan's own model df.
#
# Fallback: for lavaan versions without lav_pt_ndat(), the count is derived
# structurally from the parameter table (see plavaan_n_stats_structural()).
plavaan_n_stats <- function(x) {
    pt <- lavaan::parTable(x)
    if (exists("lav_pt_ndat", envir = asNamespace("lavaan"))) {
        return(lavaan::lav_pt_ndat(pt))
    }
    plavaan_n_stats_structural(pt)
}

# Structural fallback for plavaan_n_stats(): count the sample moments per
# group from the parameter table. For a group with m observed variables of
# which c are continuous, this is the m * (m - 1) / 2 unique (co)variances,
# plus the c continuous variances (ordinal variances are fixed to 1 in the
# polychoric/biserial metric), plus the ordinal threshold parameters (one
# per category beyond the first, per ordinal variable), plus any
# continuous-variable intercepts.
#
# Two subtleties that break a naive k * (k + 1) / 2 count:
# - Threshold rows (op "|") carry the ordinal category counts; their rhs
#   ("t1", "t2", ...) are not variables and must be excluded, and each such
#   row is itself a sample statistic.
# - Ordinal "means" are absorbed into the thresholds, so only intercepts on
#   continuous observed variables are counted (WLSMV adds intercepts for
#   continuous variables automatically).
plavaan_n_stats_structural <- function(pt) {
    gr <- pt$group
    groups <- if (all(is.na(gr))) 1L else unique(gr[!is.na(gr)])
    n_stats <- 0
    for (g in groups) {
        sel <- if (length(groups) == 1 && is.numeric(groups)) {
            rep(TRUE, nrow(pt))
        } else {
            gr == g
        }
        # Threshold rows (op "|") hold category counts, not variables.
        nth <- pt$op != "|"
        allv <- unique(c(pt$lhs[sel & nth], pt$rhs[sel & nth]))
        # Intercept rows (op "~1") have an empty rhs; drop empty names.
        allv <- allv[nzchar(allv)]
        lat <- unique(pt$lhs[sel & pt$op == "=~"])
        obs <- setdiff(allv, lat)
        m <- length(obs)
        # Ordinal observed variables have threshold rows; the group's
        # threshold parameters equal the number of threshold rows.
        ord <- intersect(unique(pt$lhs[sel & pt$op == "|"]), obs)
        n_thresh <- sum(sel & pt$op == "|")
        # Intercept statistics on continuous observed variables.
        n_mean <- sum(sel & pt$op == "~1" & pt$lhs %in% setdiff(obs, ord))
        n_stats <- n_stats +
            m * (m - 1) / 2 + (m - length(ord)) + n_thresh + n_mean
    }
    n_stats
}
