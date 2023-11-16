

# colours for plotting
# blue <- 'steelblue3'
# purple <- '#C3A0E8'
#
# yellow <- RColorBrewer::brewer.pal(8, 'Set2')[6]
# blue_green <- colorRampPalette(c('blue', green))(10)[8]
# yellow_green <- colorRampPalette(c('yellow', green))(10)[8]
# orange <- RColorBrewer::brewer.pal(8, 'Set2')[2]
# pink <- RColorBrewer::brewer.pal(8, 'Set2')[4]
# fifo <- '#A8EB12'



#' takes in calculated values from timeseries samples, plot ribbon plot, adding
#' dates and states labels, and optionally add in observed data as overlaying
#' graphical elements
#'
#' @param simulations
#' @param dates
#' @param states
#' @param base_colour
#' @param start_date
#' @param case_forecast
#' @param type
#' @param case_validation_data
#'
#' @return
#' @export
#'
#' @examples
#'



p <- plot_interval_curves(reff_sims[[1]],
                     dates = days_infection,
                     start_date = as.Date('2023-10-01'),
                     end_date = as.Date('2023-11-01'),
                     jurisdictions = jurisdictions)


plot_interval_curves <- function(
        simulations,
        dates,
        start_date,
        end_date,
        jurisdictions = jurisdictions) {

    cat('plotting instanteneous reproduction number by infection date!')
    ylab_name <- expression(R['eff']~from~'locally-acquired'~cases)
    # ribbon_colour <- RColorBrewer::brewer.pal(8, 'Set2')[1]

    ### HERE
    #calculate mean and CI values for ribbon plot

    dates_df <- data.frame(ymd = dates,
                           date = as.character(1:length(dates)))

    df <- do.call(rbind, lapply(1:length(jurisdictions),
                                function (x) {
                                    dat <- simulations[ , -(1:5), x] |>
                                        as.data.frame() |>
                                        dplyr::mutate(.draw = row_number(),
                                                      jur = as.factor(jurisdictions[x]))
                                })) |>
        # df <- simulations[,-(1:5),1] |>
        # as.data.frame() |>
        # dplyr::mutate(.draw = row_number()) |>
        tidyr::pivot_longer(-c(.draw, jur),
                            names_to = 'date',
                            values_to = 'reff') |>
        dplyr::mutate(date = as.numeric(stringr::str_replace_all(date, 'V', ''))) |>
        merge(dates_df)



    # rename value remove V and add date info.
    df2 <- do.call(rbind, lapply(unique(df$jur),
                                 function (x) {
                                     df2 <- df |>
                                         dplyr::filter(jur == x) |>
                                         dplyr::group_by(date) |>
                                         curve_interval(reff,
                                                        .exclude = c('.draw', 'jur'),
                                                        .width = c(0.95)) |>
                                         dplyr::mutate(jur = as.factor(x))

                                 })) |>
        merge(dates_df)



    p <- ggplot(aes(x = ymd, y = reff), data = df2) +
        geom_lineribbon(aes(ymin = .lower, ymax = .upper, group = 1),
                        linewidth = .5) +
        # geom_line(aes(group = .draw), alpha = 0.005, data = df) +
        scale_fill_brewer() +
        theme_ggdist() +
        ggplot2::scale_y_continuous(name = ylab_name) +
        xlab('Date') +
        ggplot2::facet_wrap(~ jur, ncol = 2, scales = 'free') +
        ggplot2::scale_x_date(date_breaks = '1 week',
                              date_minor_breaks = '1 day',
                              date_labels = '%d%b',
                              limits = c(as.Date(start_date), end_date)) +
        ggplot2::geom_hline(yintercept = 1, linetype = 'dotted', colour = 'black')

    p
}



dates_df <- data.frame(ymd = days_infection,
                       date = as.character(1:length(days_infection)))


df <- reff_sims[[1]][,-(1:5),1] |>
    as.data.frame() |>
    dplyr::mutate(.draw = row_number(),
                  jur = 'SA') |>
    tidyr::pivot_longer(-c(.draw, jur),
                        names_to = 'date',
                        values_to = 'reff') |>
    dplyr::mutate(date = as.numeric(stringr::str_replace_all(date, 'V', ''))) |>
    merge(dates_df)

df2 <- df |>
    dplyr::filter(jur == 'SA') |>
    dplyr::group_by(date) |>
    curve_interval(reff, .exclude = c('.draw', 'jur'), .width = .95) |>
    dplyr::mutate(jur = as.factor('jur')) |>
    merge(dates_df)


ggplot() +
    geom_line(aes(ymd, reff, group = .draw), data = df, alpha = 0.005) +
    theme_ggdist() +
    ggplot2::scale_x_date(date_breaks = '1 week',
                          date_minor_breaks = '1 day',
                          date_labels = '%d%b',
                          limits = c(as.Date('2023-10-01'), as.Date('2023-11-01'))) +
    geom_line(aes(ymd, reff), data = df2)


ggplot(aes(x = ymd, y = reff), data = df2) +
    geom_lineribbon(aes(ymin = .lower, ymax = .upper, group = 1),
                    linewidth = .5) +
    # geom_line(aes(group = .draw), alpha = 0.005, data = df) +
    scale_fill_brewer() +
    theme_ggdist() +
    xlab('Date') +
    ggplot2::facet_wrap(~ jur, ncol = 2, scales = 'free') +
    ggplot2::scale_x_date(date_breaks = '1 week',
                          date_minor_breaks = '1 day',
                          date_labels = '%d%b',
                          limits = c(as.Date('2023-10-01'), as.Date('2023-11-01'))) +
    ggplot2::geom_hline(yintercept = 1, linetype = 'dotted', colour = 'black')
