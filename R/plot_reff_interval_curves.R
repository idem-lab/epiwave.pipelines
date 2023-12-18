#' Plot interval curves for Reff
#'
#' @param filename output file name
#' @param df long form data frame of estimates from draws
#' @param dates infection dates sequence
#' @param start_date start date for plot
#' @param end_date end date for plot
#' @param jurisdictions vector of jurisdictions
#'
#' @importFrom dplyr filter group_by mutate
#' @importFrom ggdist curve_interval geom_lineribbon theme_ggdist
#' @importFrom ggplot2 aes facet_wrap geom_hline geom_line ggplot
#'  scale_fill_brewer scale_x_date scale_y_continuous xlab
#'
#' @return plot of interval curves for Reff
#' @export
plot_reff_interval_curves <- function (filename,
                                       df,
                                       dates,
                                       start_date,
                                       end_date,
                                       jurisdictions) {

    cat('plotting instanteneous reproduction number by infection date')
    ylab_name <- expression(R['eff']~from~'locally-acquired'~cases)

    dates_df <- data.frame(ymd = dates,
                           date_id = as.character(1:length(dates)))

    # rename value remove V and add date info.
    df2 <- do.call(rbind, lapply(
        unique(df$jur),
        function (x) {
            df2 <- df |>
                dplyr::filter(jur == x) |>
                dplyr::group_by(date_id) |>
                ggdist::curve_interval(
                    val,
                    .exclude = c('.draw', 'jur'),
                    .width = c(0.95)) |>
                dplyr::mutate(jur = as.factor(x))
        })) |>
        merge(dates_df)

    p <- ggplot2::ggplot(ggplot2::aes(x = ymd, y = val), data = df2) +
        ggdist::geom_lineribbon(ggplot2::aes(ymin = .lower,
                                             ymax = .upper, group = 1),
                                linewidth = 0.0000001,
                                col = 'light grey') +
        ggplot2::geom_line(ggplot2::aes(group = .draw),
                           alpha = 0.002, data = df) +
        ggplot2::scale_fill_brewer() +
        ggdist::theme_ggdist() +
        ggplot2::scale_y_continuous(name = ylab_name) +
        ggplot2::xlab('Date') +
        ggplot2::facet_wrap(~ jur, ncol = 2, scales = 'free') +
        ggplot2::scale_x_date(date_breaks = '1 week',
                              date_minor_breaks = '1 day',
                              date_labels = '%d %b',
                              limits = c(as.Date(start_date), end_date)) +
        ggplot2::geom_hline(yintercept = 1,
                            linetype = 'dotted',
                            colour = 'black')

    pdf(filename, width = 10, height = 12)
    print(p)
    dev.off()

    return(filename)
}
