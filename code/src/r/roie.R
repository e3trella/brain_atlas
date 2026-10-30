ROIE <- function(crosstab) {
  ## Compute the Ro/e value from a contingency table.
  ## Args:
  ##   crosstab: contingency table of the given distribution.
  ## Returns:
  ##   The Ro/e matrix, with the same dimnames as crosstab.
  rowsum.matrix <- matrix(0, nrow = nrow(crosstab), ncol = ncol(crosstab))
  rowsum.matrix[, 1] <- rowSums(crosstab)
  colsum.matrix <- matrix(0, nrow = ncol(crosstab), ncol = ncol(crosstab))
  colsum.matrix[1, ] <- colSums(crosstab)
  allsum <- sum(crosstab)
  roie <- divMatrix(crosstab, rowsum.matrix %*% colsum.matrix / allsum)
  row.names(roie) <- row.names(crosstab)
  colnames(roie) <- colnames(crosstab)
  roie
}

divMatrix <- function(m1, m2) {
  ## Divide each element of two same-dimension matrices in turn.
  ## Args:
  ##   m1: the first matrix.
  ##   m2: the second matrix.
  ## Returns:
  ##   A matrix with the same dimension and dimnames as m1,
  ##   where result[i, j] = m1[i, j] / m2[i, j].
  dim_m1 <- dim(m1)
  dim_m2 <- dim(m2)
  if (sum(dim_m1 == dim_m2) == 2) {
    div.result <- matrix(rep(0, dim_m1[1] * dim_m1[2]), nrow = dim_m1[1])
    row.names(div.result) <- row.names(m1)
    colnames(div.result) <- colnames(m1)
    for (i in 1:dim_m1[1]) {
      for (j in 1:dim_m1[2]) {
        div.result[i, j] <- m1[i, j] / m2[i, j]
      }
    }
    return(div.result)
  } else {
    warning("The dimensions of m1 and m2 are different")
  }
}
