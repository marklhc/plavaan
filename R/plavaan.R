#' Penalized lavaan Fit Objects (\sQuote{plavaan} Class)
#'
#' An S4 subclass of the \sQuote{lavaan} class, defined by
#' \code{setClass("plavaan", contains = "lavaan")}. Objects returned by
#' [penalized_est()] are \sQuote{plavaan} objects: they carry a
#' \sQuote{penalized} attribute recording the penalty specification, a
#' shared \sQuote{plavaan.cache} environment attribute (used to cache the
#' frozen refit), a \code{fitmeasures()} method that reports fit indices at
#' the effective degrees of freedom (see [effective_df()]), and a
#' \code{summary()} method that prints the frozen fit's summary together
#' with the effective parameter and degrees-of-freedom counts. Fit
#' evaluation (\code{fitmeasures()} and the chi-square test in
#' \code{summary()}) is experimental and only available when the fit was
#' created with a non-\code{"none"} \code{test} (see
#' [penalized_est()]).
#'
#' @name plavaan
#' @keywords internal
#' @import methods
#' @importFrom lavaan summary
NULL

setClass("plavaan", contains = "lavaan")

# Freeze all parameters at their current estimates and refit.
#
# The frozen refit is the basis for the fit measures reported for the
# plavaan class: with every parameter fixed, the model test statistic is
# the chi-square evaluated at the penalized estimates, and the (nominal)
# model degrees of freedom equal the number of sample statistics.
# plavaan_patch() then replaces them with the effective values.
#
# Verified mechanics (lavaan 0.7-2):
# - ptf$ustart must be backfilled with ptf$est. With ustart = NA,
#   lav_start() overwrites the start values with its default CFA/fabin
#   starts, and since the npar = 0 code path does
#   parTable$est <- parTable$start, the frozen fit would be evaluated at
#   the WRONG (fabin) parameters instead of the penalized estimates.
# - optim.force.converged = TRUE is required: the zero-free-parameter
#   branch of lav_step11_estoptim() sets converged = FALSE, and
#   lav_step14_test() skips the test computation unless
#   attr(x, "converged") is TRUE.
# - opts$test is set from the fit's recorded test specification (see the
#   `test` argument of penalized_est()): "none" skips the model test
#   (fit evaluation disabled), "Chisq" / "SatorraBentler" compute the
#   chi-square test that fitmeasures() and the summary read from the test
#   slot.
# - Pass both slotData and slotSampleStats from the original fit: for fits
#   created from sample statistics, @Data alone is insufficient (the refit
#   fails), and providing the original @SampleStats also avoids recomputing
#   the sample moments. Works for single-group, meanstructure, and group
#   models.
plavaan_freeze <- function(fit, test) {
    ptf <- lavaan::parTable(fit)
    ptf$free <- rep(0L, nrow(ptf))
    ptf$start <- ptf$est
    ptf$ustart <- ptf$est
    opts <- fit@Options
    opts$do.fit <- TRUE
    opts$test <- test
    opts$start <- NULL
    opts$optim.force.converged <- TRUE
    opts$se <- "none"
    lavaan::lavaan(
        ptf,
        slotOptions = opts,
        slotData = fit@Data,
        slotSampleStats = fit@SampleStats
    )
}

#' @importFrom stats pchisq
plavaan_patch <- function(frozen, fit, npar_eff) {
    n_stats <- plavaan_n_stats(fit)
    df_eff <- n_stats - npar_eff
    # A chi-square test is present only when fit evaluation was requested
    # (a non-"none" test); with test = "none" no p-value is computed.
    has_chisq <- any(vapply(
        frozen@test,
        function(t) identical(t$refdistr, "chisq"),
        logical(1)
    ))
    # A non-positive effective df means the effective model is still
    # under-identified, so the chi-square p-value is undefined (pchisq()
    # would return NaN with a warning). Report it as NA and warn once, but
    # only when a chi-square test was actually requested.
    if (has_chisq && df_eff <= 0) {
        warning(
            "The effective model degrees of freedom are non-positive (",
            round(df_eff, 2),
            "); the chi-square p-value is not available.",
            call. = FALSE
        )
    }
    for (nm in names(frozen@test)) {
        if (frozen@test[[nm]]$refdistr == "chisq") {
            frozen@test[[nm]]$df <- df_eff
            frozen@test[[nm]]$pvalue <- if (df_eff > 0) {
                pchisq(frozen@test[[nm]]$stat, df_eff, lower.tail = FALSE)
            } else {
                NA_real_
            }
        }
    }
    # n is the total number of observations. Use @loglik$ntotal rather
    # than lavInspect(fit, "nobs"), which is a per-group vector in group
    # models.
    ll <- frozen@loglik
    n <- ll$ntotal
    ll$npar <- npar_eff
    ll$AIC <- -2 * ll$loglik + 2 * npar_eff
    ll$BIC <- -2 * ll$loglik + npar_eff * log(n)
    ll$BIC2 <- -2 * ll$loglik + npar_eff * log(n / 2)
    frozen@loglik <- ll
    # The "npar" fitmeasure and the summary header read @Model@nx.free;
    # drive them with the rounded effective parameter count.
    frozen@Model@nx.free <- as.integer(round(npar_eff))
    # Report the original (penalized) run's optimizer status in the
    # summary: the frozen refit performs no optimization of its own, so
    # without this the summary would print "did not run (perhaps
    # do.fit = FALSE)?" (iterations <= 0), and -- since the frozen refit
    # is force-converged -- "ended normally" even when the penalized run
    # did not converge. (The 1L floor keeps the banner away in the
    # degenerate case of a zero-iteration penalized run.)
    frozen@optim$iterations <- max(1L, as.integer(fit@optim$iterations))
    frozen@optim$converged <- fit@optim$converged
    frozen
}

# Orchestrates the freeze-and-refit, with caching. The frozen object is
# cached in the shared "plavaan.cache" environment attribute (created by
# make_penalized_fit()), so repeated fitmeasures()/summary() calls do not
# refit. An environment is used rather than a plain attribute because R
# copies objects on write: `attr(fit, ...) <- value` inside a method would
# only modify a local copy, while all copies of the fit share the same
# environment.
plavaan_frozen <- function(fit) {
    cache <- attr(fit, "plavaan.cache")
    if (!is.null(cache) && !is.null(cache$frozen)) {
        return(cache$frozen)
    }
    spec <- plavaan_penalty_spec(fit)
    npar_eff <- plavaan_npar_eff(fit, spec)
    frozen <- tryCatch(plavaan_freeze(fit, spec$test), error = function(e) e)
    if (inherits(frozen, "error")) {
        stop(
            "Failed to compute fit measures for this penalized fit: ",
            conditionMessage(frozen),
            call. = FALSE
        )
    }
    frozen <- plavaan_patch(frozen, fit, npar_eff)
    if (!is.null(cache)) {
        cache$frozen <- frozen
    }
    frozen
}

#' Fit Measures for Penalized lavaan Fits
#'
#' \code{fitmeasures()} on a \sQuote{plavaan} object freezes all parameters
#' at the penalized estimates, refits the model, and reports the fit
#' measures with the effective degrees of freedom and effective parameter
#' count (see [effective_df()]). The frozen refit is cached in the shared
#' \sQuote{plavaan.cache} environment attribute, so repeated calls do not
#' refit.
#'
#' Fit evaluation is experimental and disabled by default: it requires the
#' fit to have been created with [penalized_est()] using a non-\code{"none"}
#' \code{test}. When \code{test = "none"} (the default), \code{fitmeasures()}
#' returns \code{NULL} with a message explaining how to enable it. When the
#' test is enabled, an experimental notice is shown on every call.
#'
#' @noRd
#' @importFrom lavaan fitmeasures
setMethod("fitmeasures", signature = "plavaan", function(
    object,
    fit_measures = "all",
    baseline_model = NULL,
    h1_model = NULL,
    fm_args = NULL,
    output = "vector",
    level = NULL,
    ...
) {
    spec <- plavaan_penalty_spec(object)
    if (identical(spec$test, "none")) {
        message(
            "Fit measures are not available for this penalized fit: fit ",
            "evaluation is experimental and disabled by default (test = ",
            "'none'). Re-run penalized_est() with test = 'Chisq' (or ",
            "'SatorraBentler') to enable fit measures and the chi-square test."
        )
        return(invisible(NULL))
    }
    message(
        "Fit evaluation for penalized fits is experimental; interpret fit ",
        "indices with caution."
    )
    # S4 does not forward a generic's named formals to a method declared as
    # (object, ...): arguments such as fit_measures, output, and level would
    # be matched to lavaan's fitmeasures() generic and then silently dropped
    # before reaching this method. Declare them explicitly (mirroring lavaan's
    # signature) and pass them through. fm_args uses a NULL sentinel so that
    # lavaan's own default is used when the caller omits it.
    call_args <- list(
        fit_measures = fit_measures,
        baseline_model = baseline_model,
        h1_model = h1_model,
        output = output,
        level = level
    )
    if (!is.null(fm_args)) {
        call_args$fm_args <- fm_args
    }
    do.call(lavaan::fitmeasures, c(list(plavaan_frozen(object)), call_args, list(...)))
})

#' Summary Method for Penalized lavaan Fits
#'
#' \code{summary()} on a \sQuote{plavaan} object. The original penalized
#' fit was not run with \code{do.fit = TRUE}, so the summary is computed
#' from a "frozen" refit with all parameters fixed at the penalized
#' estimates (see [effective_df()]), which makes fit indices available. A
#' message reports the effective number of parameters and degrees of
#' freedom.
#'
#' The chi-square test in the summary is part of the experimental fit
#' evaluation and is only shown when the fit was created with
#' [penalized_est()] using a non-\code{"none"} \code{test}; in that case
#' an experimental notice is also shown. With the default
#' \code{test = "none"}, the summary shows the model and parameter estimates
#' (at the effective parameter count) but no chi-square test.
#'
#' When the fit was created with [penalized_est()] using
#' \code{se = "robust.huber.white"}, the summary displays the standard
#' errors in the parameter estimates table (and the resulting z-values and
#' p-values), as for a plain lavaan fit. With the default \code{se = "none"},
#' no standard errors are shown.
#'
#' An S4 method is used (rather than an S3 method) because lavaan
#' registers S4 methods for \code{summary()} on the \sQuote{lavaan} class,
#' and S4 dispatch takes precedence for S4 objects: an S3
#' \code{summary.plavaan} would never be reached.
#'
#' @param object A \sQuote{plavaan} object, typically the return value of
#'   [penalized_est()].
#' @param ... Further arguments passed to the \code{summary} method for
#'   \sQuote{lavaan} objects.
#' @return The summary object of the frozen lavaan refit, as returned by
#'   the \code{summary} method for \sQuote{lavaan} objects.
#' @noRd
setMethod("summary", signature = "plavaan", function(object, ...) {
    frozen <- plavaan_frozen(object)
    # The frozen refit is run with se = "none" so that lavaan does not
    # (re)compute standard errors for the zero-free-parameter model. Its
    # ParTable nonetheless carries the original fit's standard errors: they
    # are inherited because plavaan_freeze() feeds parTable(fit) -- which
    # includes the se column -- into lavaan(). Restore the original se
    # setting so that summary() displays the standard errors, matching the
    # behaviour of a plain lavaan fit with se = "robust.huber.white". The
    # guard checks the frozen ParTable (what summary actually prints) so the
    # flag is only raised when SE values are genuinely present. Mutating the
    # local copy does not affect the cached frozen object (copy-on-write).
    frozen_se <- lavaan::parTable(frozen)$se
    if (
        identical(object@Options$se, "robust.huber.white") &&
            any(!is.na(frozen_se) & frozen_se > 0)
    ) {
        frozen@Options$se <- "robust.huber.white"
        frozen@vcov$se <- "robust.huber.white"
    }
    s <- summary(frozen, ...)
    spec <- plavaan_penalty_spec(object)
    npar_eff <- plavaan_npar_eff(object, spec)
    pen_fn_label <- if (is.function(spec$pen_fn)) "custom" else spec$pen_fn
    n_stats <- plavaan_n_stats(object)
    message(
        "Penalized fit (w = ", spec$w, ", eps = ", spec$eps,
        ", penalty = ", pen_fn_label, "): ",
        "effective npar = ", round(npar_eff, 2), ", effective df = ",
        round(n_stats - npar_eff, 2),
        " (nominal df = ", n_stats - length(plavaan_estimates(object)), ")."
    )
    if (!identical(spec$test, "none")) {
        message(
            "Fit evaluation for penalized fits is experimental; interpret ",
            "the chi-square test and fit indices with caution."
        )
    }
    s
})
