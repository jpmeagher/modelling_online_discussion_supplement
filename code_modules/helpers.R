#' Helper functions for summarising posterior distributions and tables
#'
#' Used across figures_and_tables qmd documents

#' Posterior mean (standard deviation) summary
#'
#' @param x Numeric vector of posterior draws
#' @param digits Number of decimal places
#' @return Character string "mean (sd)"
mean_sd <- function(x, digits = 2) {
  n <- sum(!is.na(x))
  if (n == 0) {
    return("-")
  }
  m <- round(mean(x), digits = digits)
  s <- round(stats::sd(x), digits = digits)
  paste0(m, " (", s, ")")
}

#' Mean (standard error) summary
#'
#' @param x Numeric vector
#' @param digits Number of decimal places
#' @return Character string "mean (se)"
mean_se <- function(x, digits = 2) {
  n <- sum(!is.na(x))
  if (n == 0) {
    return("-")
  }
  m <- mean(x)
  s <- sd(x) / sqrt(n)
  paste0(
    sprintf(m, fmt = paste0("%#.", digits, "f")),
    " (",
    sprintf(s, fmt = paste0("%#.", digits, "f")),
    ")"
  )
}

#' Mean (credible interval) summary
#'
#' @param x Numeric vector
#' @param alpha Significance level for the interval
#' @param digits Number of decimal places
#' @return Character string "mean (lower, upper)"
mean_ci <- function(x, alpha = 0.05, digits = 2) {
  if (all(is.na(x))) {
    return("-")
  }
  m <- round(mean(x), digits = digits)
  quants <- round(
    quantile(x, prob = c(alpha / 2, 1 - alpha / 2)),
    digits = digits
  )
  paste0(m, " (", quants[1], ", ", quants[2], ")")
}

#' Sum (standard error) summary
#'
#' @param x Numeric vector
#' @param digits Number of decimal places
#' @return Character string "sum (se)"
sum_se <- function(x, digits = 2) {
  n <- sum(!is.na(x))
  if (n == 0) {
    return("-")
  }
  m <- sum(x)
  s <- stats::sd(x) * sqrt(n)
  paste0(
    sprintf(m, fmt = paste0("%#.", digits, "f")),
    " (",
    sprintf(s, fmt = paste0("%#.", digits, "f")),
    ")"
  )
}

#' Standard error
#'
#' @param x Numeric vector
#' @return Standard error of x
se <- function(x) {
  sd(x) / sqrt(length(x))
}
