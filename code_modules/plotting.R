#' Consistent ggplot2 theme for manuscript figures
#'
#' Based on theme_classic() with adjusted font sizes for publication
#'
#' @export
my_theme <-
  ggplot2::theme_classic() +
  ggplot2::theme(
    axis.text = ggplot2::element_text(size = 8),
    axis.title = ggplot2::element_text(size = 8),
    plot.title = ggplot2::element_text(size = 12),
    plot.subtitle = ggplot2::element_text(size = 10),
    legend.text = ggplot2::element_text(size = 10)
  )
