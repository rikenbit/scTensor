#' Build Signed Cell-Cell Interaction Graphs
#'
#' Constructs positive and negative cell-cell interaction matrices from
#' ligand/receptor expression and a ligand-receptor database with sign
#' annotations. The LR scores \eqn{X_{ijk} = L_{ik} R_{jk}} are aggregated
#' into \eqn{A^+} (activating) and \eqn{A^-} (inhibitory/death-inducing)
#' matrices based on the sign of each LR pair.
#'
#' @param lig_expr Matrix of ligand expression (cells x LR pairs), or
#'   (cell types x LR pairs) if pre-aggregated.
#' @param rec_expr Matrix of receptor expression (cells x LR pairs), same
#'   dimensions as \code{lig_expr}.
#' @param lr_sign Integer or numeric vector of length \code{ncol(lig_expr)}.
#'   +1 for activating LR pairs, -1 for inhibitory/death-inducing LR pairs.
#' @return A list with components:
#'   \describe{
#'     \item{A_pos}{Positive (activating) cell-cell interaction matrix (n x n)}
#'     \item{A_neg}{Negative (inhibitory) cell-cell interaction matrix (n x n)}
#'     \item{X_pos}{3D array of per-LR positive scores (sender x receiver x LR)}
#'     \item{X_neg}{3D array of per-LR negative scores (sender x receiver x LR)}
#'     \item{edge_table}{Data frame with columns: sender, receiver, lr_pair,
#'       score, sign}
#'   }
#' @export
#' @examples
#' # 3 cell types, 4 LR pairs (2 activating, 2 inhibitory)
#' lig <- matrix(c(1,0,0, 0,1,0, 1,1,0, 0,0,1), nrow = 3)
#' rec <- matrix(c(0,1,0, 1,0,1, 0,0,1, 1,0,0), nrow = 3)
#' signs <- c(1, 1, -1, -1)
#' result <- BuildSignedCCI(lig, rec, signs)
#' result$A_pos
#' result$A_neg
BuildSignedCCI <- function(lig_expr, rec_expr, lr_sign) {
  if (!is.matrix(lig_expr) || !is.matrix(rec_expr))
    stop("lig_expr and rec_expr must be matrices")
  if (nrow(lig_expr) != nrow(rec_expr) || ncol(lig_expr) != ncol(rec_expr))
    stop("lig_expr and rec_expr must have the same dimensions")

  n <- nrow(lig_expr)    # number of cell types
  K <- ncol(lig_expr)    # number of LR pairs

  if (length(lr_sign) != K)
    stop("lr_sign must have length equal to ncol(lig_expr)")
  if (!all(lr_sign %in% c(-1, 1)))
    stop("lr_sign must contain only +1 or -1")

  # Compute per-LR score tensor: X[i, j, k] = lig_expr[i, k] * rec_expr[j, k]
  X_pos <- array(0, dim = c(n, n, K))
  X_neg <- array(0, dim = c(n, n, K))

  for (k in seq_len(K)) {
    score_k <- outer(lig_expr[, k], rec_expr[, k])
    if (lr_sign[k] == 1) {
      X_pos[, , k] <- score_k
    } else {
      X_neg[, , k] <- score_k
    }
  }

  # Aggregate over LR dimension
  A_pos <- apply(X_pos, c(1, 2), sum)
  A_neg <- apply(X_neg, c(1, 2), sum)

  # Edge table
  edges <- list()
  for (k in seq_len(K)) {
    for (i in seq_len(n)) {
      for (j in seq_len(n)) {
        score <- lig_expr[i, k] * rec_expr[j, k]
        if (score > 0) {
          edges[[length(edges) + 1L]] <- data.frame(
            sender = i, receiver = j, lr_pair = k,
            score = score, sign = lr_sign[k],
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }

  if (length(edges) > 0) {
    edge_table <- do.call(rbind, edges)
  } else {
    edge_table <- data.frame(
      sender = integer(0), receiver = integer(0), lr_pair = integer(0),
      score = numeric(0), sign = integer(0),
      stringsAsFactors = FALSE
    )
  }

  list(
    A_pos = A_pos,
    A_neg = A_neg,
    X_pos = X_pos,
    X_neg = X_neg,
    edge_table = edge_table
  )
}
