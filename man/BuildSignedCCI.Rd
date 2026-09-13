\name{BuildSignedCCI}
\alias{BuildSignedCCI}
\title{Build Signed Cell-Cell Interaction Graphs}
\description{
Constructs positive and negative cell-cell interaction matrices from
ligand/receptor expression and a ligand-receptor database with sign
annotations. The LR scores \eqn{X_{ijk} = L_{ik} R_{jk}} are aggregated
into \eqn{A^+} (activating) and \eqn{A^-} (inhibitory/death-inducing)
matrices based on the sign of each LR pair.
}
\usage{
BuildSignedCCI(lig_expr, rec_expr, lr_sign)
}
\arguments{
\item{lig_expr}{
Matrix of ligand expression (cells x LR pairs), or
(cell types x LR pairs) if pre-aggregated.
}
\item{rec_expr}{
Matrix of receptor expression (cells x LR pairs), same
dimensions as \code{lig_expr}.
}
\item{lr_sign}{
Integer or numeric vector of length \code{ncol(lig_expr)}.
+1 for activating LR pairs, -1 for inhibitory/death-inducing LR pairs.
}
}
\value{
A list with components:
\describe{
  \item{A_pos}{Positive (activating) cell-cell interaction matrix (n x n)}
  \item{A_neg}{Negative (inhibitory) cell-cell interaction matrix (n x n)}
  \item{X_pos}{3D array of per-LR positive scores (sender x receiver x LR)}
  \item{X_neg}{3D array of per-LR negative scores (sender x receiver x LR)}
  \item{edge_table}{Data frame with columns: sender, receiver, lr_pair,
    score, sign}
}
}
\examples{
# 3 cell types, 4 LR pairs (2 activating, 2 inhibitory)
lig <- matrix(c(1,0,0, 0,1,0, 1,1,0, 0,0,1), nrow = 3)
rec <- matrix(c(0,1,0, 1,0,1, 0,0,1, 1,0,0), nrow = 3)
signs <- c(1, 1, -1, -1)
result <- BuildSignedCCI(lig, rec, signs)
result$A_pos
result$A_neg
}
