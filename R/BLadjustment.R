# Lawley's (1956) Bartlett correction term for the Likelihood Ratio test
# This is the procedure 'LAWLEY' as shown in Table 2 (page 13) of Jensen (1993)
# Jensen J (1993). A Historical Sketch and Some New Results on the Improved
# Log Likelihood Ratio Statistic. Scandinavian Journal of Statistics 20 (1):1-15
BLadjustment <- function(p, d2, d3, d4, d2d1, d2d2, d3d1) {
  din <- solve(d2)  # Fisher information inverse

  # Term a
  a <- 0
  for (i in 1:p) {
    for (j in 1:p) {
      term <- 0.25 * d4[i,j,,] - d3d1[i,j,,] + d2d2[i,j,,]
      a <- a + din[i,j] * sum(din * term)
    }
  }

  # Term b
  b <- 0
  for (i in 1:p) {
    for (j in 1:p) {
      for (k in 1:p) {
        for (f in 1:p) {
          for (g in 1:p) {
            for (h in 1:p) {
              term1 <-        d3[i,k,g] *   d3[j,f,h] / 6
              term2 <- 0.25 * d3[i,k,f] *   d3[j,g,h]
              term3 <-      - d3[i,k,g] * d2d1[j,h,f]
              term4 <-      - d3[i,k,f] * d2d1[j,h,g]
              term5 <-      d2d1[i,k,g] * d2d1[j,h,f]
              term6 <-      d2d1[i,k,f] * d2d1[j,h,g]

              total_term <- term1 + term2 + term3 + term4 + term5 + term6
              b <- b + din[i,j] * din[k,f] * din[g,h] * total_term
            }
          }
        }
      }
    }
  }

  return(a - b)
}

BLadjustment0 <- function(p, d2, d3, d4, d2d1, d2d2, d3d1) {
  din <- solve(d2)  # Fisher information inverse

  # Term a
  a <- 0
  for (i in 1:p) {
    for (j in 1:p) {
      for (k in 1:p) {
        for (f in 1:p) {
          term <- 0.25 * d4[i,j,k,f] - d3d1[i,j,k,f] + d2d2[i,j,k,f]
          a <- a + din[i,j] * din[k,f] * term
        }
      }
    }
  }

  # Term b
  b <- 0
  for (i in 1:p) {
    for (j in 1:p) {
      for (k in 1:p) {
        for (f in 1:p) {
          for (g in 1:p) {
            for (h in 1:p) {
              term1 <-        d3[i,k,g] *   d3[j,f,h] / 6
              term2 <- 0.25 * d3[i,k,f] *   d3[j,g,h]
              term3 <-      - d3[i,k,g] * d2d1[j,h,f]
              term4 <-      - d3[i,k,f] * d2d1[j,h,g]
              term5 <-      d2d1[i,k,g] * d2d1[j,h,f]
              term6 <-      d2d1[i,k,f] * d2d1[j,h,g]

              total_term <- term1 + term2 + term3 + term4 + term5 + term6
              b <- b + din[i,j] * din[k,f] * din[g,h] * total_term
            }
          }
        }
      }
    }
  }

  return(a - b)
}

