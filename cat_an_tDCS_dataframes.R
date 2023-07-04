library(data.table)

# here is data.table because of analysis very old, dplyr will be appeared later (or not)

#----cathodal dt--------------------------------------------------------

cathodal <-
  fread("tDCS_cathodal.csv",
        sep = ";",
        header = TRUE)

# View(cathodal)

#dt with rejected only
cathodal_rejected <- cathodal[Choice == "Rejected", ]

#dt with selected only
cathodal_selected <- cathodal[Choice == "Selected", ]


#dt with single difference insteas of repeated measured Evaluation

cathodal_dif <-
  cathodal[Rating == "R1", ]  #just for creating dataframe with certain number of rows


cathodal_dif[, 'Evaluation'] <-
  cathodal[Rating == "R2", Evaluation] -
  cathodal[Rating == "R1", Evaluation]
colnames(cathodal_dif)[colnames(cathodal_dif) == "Evaluation"] <-
  "Difference"
cathodal_dif[, 'Rating'] <- NULL

# View(cathodal_dif)

#----anodal dt--------------------------------------------------------

anodal <-
  fread("tDCS_anodal2.csv",
        sep = ";",
        header = TRUE)

# View(anodal)

#dt with rejected only
anodal_rejected <- anodal[Choice == "Rejected", ]

#dt with selected only
anodal_selected <- anodal[Choice == "Selected", ]


#dt with single difference insteas of repeated measured Evaluation

anodal_dif <-
  anodal[Rating == "R1", ]  #just for creating dataframe with certain number of rows


anodal_dif[, 'Evaluation'] <-  anodal[Rating == "R2", Evaluation] -
  anodal[Rating == "R1", Evaluation]
colnames(anodal_dif)[colnames(anodal_dif) == "Evaluation"] <-
  "Difference"
anodal_dif[, 'Rating'] <- NULL


#----   function for extracting t values and DFs (not mine)-----------------

pairwise.t.test.with.t.and.df <- function (x, g, p.adjust.method = p.adjust.methods, pool.sd = !paired, 
                                           paired = FALSE, alternative = c("two.sided", "less", "greater"), 
                                           ...) 
{
  if (paired & pool.sd) 
    stop("pooling of SD is incompatible with paired tests")
  DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(g)))
  g <- factor(g)
  p.adjust.method <- match.arg(p.adjust.method)
  alternative <- match.arg(alternative)
  if (pool.sd) {
    METHOD <- "t tests with pooled SD"
    xbar <- tapply(x, g, mean, na.rm = TRUE)
    s <- tapply(x, g, sd, na.rm = TRUE)
    n <- tapply(!is.na(x), g, sum)
    degf <- n - 1
    total.degf <- sum(degf)
    pooled.sd <- sqrt(sum(s^2 * degf)/total.degf)
    compare.levels <- function(i, j) {
      dif <- xbar[i] - xbar[j]
      se.dif <- pooled.sd * sqrt(1/n[i] + 1/n[j])
      t.val <- dif/se.dif
      if (alternative == "two.sided") 
        2 * pt(-abs(t.val), total.degf)
      else pt(t.val, total.degf, lower.tail = (alternative == 
                                                 "less"))
    }
    compare.levels.t <- function(i, j) {
      dif <- xbar[i] - xbar[j]
      se.dif <- pooled.sd * sqrt(1/n[i] + 1/n[j])
      t.val = dif/se.dif 
      t.val
    }       
  }
  else {
    METHOD <- if (paired) 
      "paired t tests"
    else "t tests with non-pooled SD"
    compare.levels <- function(i, j) {
      xi <- x[as.integer(g) == i]
      xj <- x[as.integer(g) == j]
      t.test(xi, xj, paired = paired, alternative = alternative, 
             ...)$p.value
    }
    compare.levels.t <- function(i, j) {
      xi <- x[as.integer(g) == i]
      xj <- x[as.integer(g) == j]
      t.test(xi, xj, paired = paired, alternative = alternative, 
             ...)$statistic
    }
    compare.levels.df <- function(i, j) {
      xi <- x[as.integer(g) == i]
      xj <- x[as.integer(g) == j]
      t.test(xi, xj, paired = paired, alternative = alternative, 
             ...)$parameter
    }
  }
  PVAL <- pairwise.table(compare.levels, levels(g), p.adjust.method)
  TVAL <- pairwise.table.t(compare.levels.t, levels(g), p.adjust.method)
  if (pool.sd) 
    DF <- total.degf
  else
    DF <- pairwise.table.t(compare.levels.df, levels(g), p.adjust.method)           
  ans <- list(method = METHOD, data.name = DNAME, p.value = PVAL, 
              p.adjust.method = p.adjust.method, t.value = TVAL, dfs = DF)
  class(ans) <- "pairwise.htest"
  ans
}


pairwise.table.t <- function (compare.levels.t, level.names, p.adjust.method) 
{
  ix <- setNames(seq_along(level.names), level.names)
  pp <- outer(ix[-1L], ix[-length(ix)], function(ivec, jvec) sapply(seq_along(ivec), 
                                                                    function(k) {
                                                                      i <- ivec[k]
                                                                      j <- jvec[k]
                                                                      if (i > j)
                                                                        compare.levels.t(i, j)               
                                                                      else NA
                                                                    }))
  pp[lower.tri(pp, TRUE)] <- pp[lower.tri(pp, TRUE)]
  pp
}





