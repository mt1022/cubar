#' Bulmer model (1987, Nature)
#' @param uk_ratio u/k
#' @references
#' Bulmer M. Coevolution of codon usage and transfer RNA abundance. Nature. 1987
#' Feb 19-25;325(6106):728-30. doi: 10.1038/325728a0. PMID: 2434856.
model_b87 <- function(uk_ratio = 1){
    c_ratio <- seq(0.001, 0.999, 0.001)  # C_2 / C_1
    us_ratio <- uk_ratio * c_ratio / (1 - c_ratio)
    q <- 0.5 + us_ratio - 0.5 * sqrt(4 * us_ratio^2 + 1)
    p <- 1 - q
    plt <- data.table::data.table(c_ratio = c_ratio, qratio = q/p,
                                  type = 'Fix tRNA level')
    plt2 <- data.table::copy(plt)
    plt2[, qratio := c_ratio^2]
    plt2[, type := 'Fix codon usage']
    plt <- rbind(plt, plt2)
    ggplot2::ggplot(plt, ggplot2::aes(
        x = .data$c_ratio, y = .data$qratio, color = .data$type)) +
        ggplot2::geom_line() +
        ggplot2::scale_color_brewer(palette = 'Dark2') +
        ggplot2::labs(x = expression(italic(C[2]/C[1])),
                      y = expression(italic(q/p)), color = NULL)
}
