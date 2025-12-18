#' Composite Pairwise Loss Function
#'
#' Computes the total loss across all pairwise combinations of rows in a matrix.
#'
#' @param x A numeric vector, matrix, or data frame. If not a matrix, it will be
#'   coerced to one after applying the transformation function.
#' @param fun A function to compute the loss for each pairwise difference.
#'   The package supports the alignment loss (`alf`) and the approximate L0 penalty
#'   (`l0a`), but users can provide custom functions as well.
#' @param trans A transformation function to apply to `x` before computing
#'   pairwise differences. Default is `identity` (no transformation).
#' @param rescale Either `"df"` (default) to rescale the total loss by the degrees
#'   of freedom (number of rows - 1), or a numeric value (likely between 0 and 1)
#'   to multiply the total loss by.
#' @param ... Additional arguments passed to the loss function `fun`.
#'
#' @return A numeric scalar representing the sum of losses across all pairwise
#'   combinations of rows.
#'
#' @details
#' The function works by:
#' \enumerate{
#'   \item Applying the transformation function `trans` to the input `x`
#'   \item Converting the result to a matrix
#'   \item Generating all possible pairwise combinations of row indices
#'   \item Computing the difference between each pair of rows
#'   \item Applying the loss function `fun` to each difference
#'   \item Summing all the individual losses
#' }
#'
#' @examples
#' # Example with a simple matrix
#' x <- matrix(runif(12), nrow = 4)
#' composite_pair_loss(x, fun = alf)
#'
#' # Example with log transformation and L2 loss
#' composite_pair_loss(x, fun = function(x) x^2, trans = log)
#'
#' @importFrom utils combn
#' @export
composite_pair_loss <- function(x, fun, trans = identity, rescale = "df", ...) {
    x <- as.matrix(trans(x))
    nrow_x <- nrow(x)
    if (nrow_x < 2) {
        return(0)
    }
    combn_idx <- combn(nrow_x, 2)
    if (rescale == "df") {
        dof <- nrow_x - 1
        ncombn <- ncol(combn_idx)
        rescale <- dof / ncombn
    }
    if (!is.numeric(rescale)) {
        stop("rescale must be 'df' or a numeric value.")
    }
    out <- fun(x[combn_idx[1, ], ] - x[combn_idx[2, ], ], ...)
    sum(out, na.rm = TRUE) * rescale
}

#' @importFrom stats na.omit
hot_gr <- function(x, hot, fun, ...) {
    out <- 0 * x
    if (is.matrix(hot)) {
        x_hot <- matrix(x[hot], nrow = nrow(hot), ncol = ncol(hot))
    } else {
        x_hot <- x[hot]
    }
    out[na.omit(as.numeric(hot))] <- na.omit(fun(x_hot, ...))
    out
}

gr_cpl <- function(
    x,
    gr_fun,
    trans = identity,
    gr_trans = function(x) 1,
    rescale = "df"
) {
    x_mat <- as.matrix(trans(x))
    combn_idx <- combn(nrow(x_mat), 2)
    diffs <- x_mat[combn_idx[1, ], , drop = FALSE] -
        x_mat[combn_idx[2, ], , drop = FALSE]
    grad_contribs <- gr_fun(diffs)
    grad <- matrix(0, nrow = nrow(x_mat), ncol = ncol(x_mat))
    for (i in seq_len(nrow(x_mat))) {
        idx1 <- which(combn_idx[1, ] == i)
        idx2 <- which(combn_idx[2, ] == i)
        grad[i, ] <- colSums(
            grad_contribs[idx1, , drop = FALSE],
            na.rm = TRUE
        ) -
            colSums(grad_contribs[idx2, , drop = FALSE], na.rm = TRUE)
    }
    grad[which(is.na(x_mat))] <- NA
    if (rescale == "df") {
        dof <- nrow(x_mat) - 1
        ncombn <- ncol(combn_idx)
        rescale <- dof / ncombn
    }
    if (!is.numeric(rescale)) {
        stop("rescale must be 'df' or a numeric value.")
    }
    as.vector(grad) * rescale * gr_trans(as.vector(x))
}

#' Loss functions
#'
#' For small eps this provides a smooth,
#' numerically stable approximation of |x|^(1/2) (i.e. the square root of
#' the absolute value). The function is vectorized over x.
#'
#' @param x Numeric vector. Input values to transform.
#' @param eps Positive numeric scalar (default .001 for `alf()` and .01
#'   for `l0a()`). Small regularization constant to avoid
#'   non-differentiability and division-by-zero issues.
#' @return Numeric vector of the same length as x.
#' @name loss
NULL
#> NULL

#' @rdname loss
#'
#' @details The ALF, (x^2 + eps)^(1/4), is useful when a smooth
#'   surrogate for sqrt(|x|) is required (for optimization or
#'   regularization) while maintaining numerical stability near x = 0.
#'
#' @examples
#' alf(0)
#' alf(c(-4, -1, 0, 1, 4))
#' alf(0.5, eps = 1e-6)
#' @export
alf <- function(x, eps = .001) {
    (x^2 + eps)^.25
}

gr_alf <- function(v, eps = .001) {
    v / (2 * (v^2 + eps)^.75)
}

#' @rdname loss
#'
#' @details L0a, x^2/(x^2 + eps), is an approximation of the L0 penalty.
#'
#' @examples
#' l0a(0)
#' l0a(c(0, 1e-3, 0.1, 1))
#' l0a(c(-2, 0, 2), eps = 1e-4)
#' @export
l0a <- function(x, eps = .01) {
    x^2 / (x^2 + eps)
}

gr_l0a <- function(v, eps = .01) {
    2 * v * eps / (v^2 + eps)^2
}
