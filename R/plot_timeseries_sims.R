#' Plot timeseries simulations
#'
#' @description Takes in calculated values from timeseries samples, plots ribbon
#'  plot, adds dates and states labels, and optionally adds in observed data
#'  as overlaying graphical elements.
#'
#' @param filename output file name
#' @param simulations simulations from greta model
#' @param type select between notification, infection, and reff
#' @param dates infection date sequence
#' @param states vector of jurisdictions
#' @param base_colour underlying color
#' @param start_date first date on x axis
#' @param end_date final date on x axis
#' @param case_validation_data case data to optionally overlay
#' @param infection_nowcast logical indicating whether we should annotate
#'  nowcast horizon
#' @param case_forecast logical whether we should forecast cases
#' @param valid_mat matrix indicating which case data are valid
#' @param nowcast_start start date for beginning of nowcast window
#' @param dim_sim whether simulations are in one or two dimensions
#' @param reff_ylim limits for y axis on reff plot
#'
#' @importFrom lubridate dmonths
#' @importFrom RColorBrewer brewer.pal
#' @importFrom tibble tibble
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom dplyr mutate filter
#' @importFrom ggplot2 ggplot aes xlab element_blank facet_wrap scale_x_date
#'  scale_alpha scale_y_continuous geom_ribbon geom_line theme geom_rug
#'  geom_point geom_vline annotate coord_cartesian geom_hline
#' @importFrom cowplot theme_cowplot panel_border
#'
#' @return simulation plot
#' @export
plot_timeseries_sims <- function (
    filename,
    simulations = NULL,
    type = c("notification", "infection", "reff"),
    dates,
    states,
    base_colour = grey(0.4),
    start_date = max(dates) - lubridate::dmonths(1),
    end_date = max(dates),
    case_validation_data = NULL,
    infection_nowcast = TRUE,
    case_forecast = FALSE,
    valid_mat = NULL,
    nowcast_start = NULL,
    dim_sim = c("1", "2"),
    reff_ylim = c(0, 2)) {

    # colours for plotting
    blue <- "steelblue3"
    purple <- "#C3A0E8"
    green <- RColorBrewer::brewer.pal(8, "Set2")[1]
    yellow <- RColorBrewer::brewer.pal(8, "Set2")[6]
    blue_green <- colorRampPalette(c("blue", green))(10)[8]
    yellow_green <- colorRampPalette(c("yellow", green))(10)[8]
    orange <- RColorBrewer::brewer.pal(8, "Set2")[2]
    pink <- RColorBrewer::brewer.pal(8, "Set2")[4]
    fifo <- "#A8EB12"

    #check type to plot
    type <- match.arg(type)
    #check dimension of simulation
    dim_sim <- match.arg(dim_sim)

    if (type == "notification") {
        cat("plotting case counts by notification date!")
        ylab_name <- "cases by notification date"
        ribbon_colour <- orange
    } else if (type == "infection") {
        cat("plotting infection counts by infection date!")
        ylab_name <- "infections by infection date"
        ribbon_colour <- pink
    } else if (type == "reff") {
        cat("plotting instanteneous reproduction number by infection date!")
        ylab_name <- expression(R["eff"]~from~"locally-acquired"~cases)
        ribbon_colour <- green
    }

    #calculate mean and CI values for ribbon plot
    mean <- apply(simulations, 2:3, FUN = "mean",na.rm = TRUE)
    ci_90_lo <- apply(simulations, 2:3, quantile, c(0.05),na.rm = TRUE)
    ci_90_hi <- apply(simulations, 2:3, quantile, c(0.95),na.rm = TRUE)
    ci_50_hi <- apply(simulations, 2:3, quantile, c(0.75),na.rm = TRUE)
    ci_50_lo <- apply(simulations, 2:3, quantile, c(0.25),na.rm = TRUE)

    #if forecast version of case timeseries is used, the calulated values
    #should already have extra forecasting days, so we simply annotate the
    #dates
    if (case_forecast) {
        extra_days <- nrow(mean) - length(dates)
        extra_dates <- dates[length(dates)] + seq(1:extra_days)
        dates <- c(dates,extra_dates)

        projection_at <- extra_dates[1]
    } else {extra_dates <- NULL} #no extra dates when no forecast



    #when using date by state dim for simulations
    if (dim_sim == "2") {

        #check if date sequence length and calculated values match
        if (length(dates) != nrow(mean)) {
            stop("Error: number of days in timeseries not equal to number of date labels provided! Did you mis-specify if forecasting is required?")
        }

        #construct a tibble of data state labels and the calculated vals
        vals <- rbind(
            mean,
            ci_90_lo,
            ci_90_hi,
            ci_50_hi,
            ci_50_lo
        )

        full_dates <- rep(dates,5)
        full_label <- rep(c("mean",
                            "ci_90_lo",
                            "ci_90_hi",
                            "ci_50_hi",
                            "ci_50_lo"),
                          each = length(dates)
        )

        df <- tibble::tibble(date = full_dates,
                             label = full_label
        )

        df <- cbind(df,vals)
        colnames(df) <- c(colnames(df)[1:2],states)

        df <- df |>
            tidyr::pivot_longer(cols = 3:ncol(df),
                                values_to = "value",
                                names_to = "state")

        df <- df |>
            tidyr::pivot_wider(names_from = label,
                               values_from = value)
    } else {
        #if dimension of simulation is a vector
        #build a df with combination of date and state as row
        df <- tibble::tibble(date = rep(dates,length(states)),
                             state = rep(states, each = length(dates))
                             )

        #subset to valid dates
        #when all dates are valid don't subset
        if (!is.null(valid_mat)) {
            df <- df[as.logical(valid_mat),]
        }
        #this should give us the same number of rows as the simulations

        vals <- cbind(
            mean,
            ci_90_lo,
            ci_90_hi,
            ci_50_hi,
            ci_50_lo
        )

        df <- cbind(df,vals)

        colnames(df)[3:7] <- c("mean",
                               "ci_90_lo",
                               "ci_90_hi",
                               "ci_50_hi",
                               "ci_50_lo")
    }



    df <- df |>
        dplyr::mutate(type = ifelse(date %in% extra_dates, "forecast","estimate"))

    df <- df |> dplyr::filter(date >= start_date, date <= end_date)


    #dynamic date label and breaks
    if (length(unique(df$date)) >= 180){
        date_breaks <- "3 month"
        date_minor_breaks <- "1 month" # minor breaks don't actually show on cowplots??
        date_labels <- "%b%y"
        x_text_angle <- 0 #legacy code for consistency with current plot styles - to discuss
        x_text_size <- 9
        x_text_hjust <- 0.5
        x_text_vjust <- 0.5
    } else if(length(unique(df$date)) < 90){
        date_breaks <- "1 week"
        date_minor_breaks <- "1 day"
        date_labels <- "%d%b"
        x_text_angle <- 0
        x_text_size <- 9
        x_text_hjust <- 0.5
        x_text_vjust <- 0.5
    } else {
        date_breaks <- "1 month"
        date_minor_breaks <- "2 weeks"
        date_labels <- "%b%y"
        x_text_angle <- 0
        x_text_size <- 9
        x_text_hjust <- 0.5
        x_text_vjust <- 0.5
    }




    #make the plot
    p <- ggplot2::ggplot(df) +

        ggplot2::aes(date, mean) +
        ggplot2::xlab(ggplot2::element_blank()) +
        ggplot2::facet_wrap(~ state, ncol = 2, scales = "free") +
        ggplot2::scale_x_date(date_breaks = date_breaks,
                     date_minor_breaks = date_minor_breaks,
                     date_labels = date_labels,
                     limits = c(as.Date(start_date),end_date)) +
        ggplot2::scale_alpha(range = c(0, 0.5)) +
        ggplot2::scale_y_continuous(name = ylab_name) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = ci_90_lo,
                        ymax = ci_90_hi),
                    fill = ribbon_colour,
                    alpha = 0.2) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = ci_50_lo,
                        ymax = ci_50_hi),
                    fill = ribbon_colour,
                    alpha = 0.5) +
        ggplot2::geom_line(ggplot2::aes(y = ci_90_lo),
                  colour = base_colour,
                  alpha = 0.8) +
        ggplot2::geom_line(ggplot2::aes(y = ci_90_hi),
                  colour = base_colour,
                  alpha = 0.8) +
        cowplot::theme_cowplot() +
        cowplot::panel_border(remove = TRUE) +
        ggplot2::theme(legend.position = "none",
              strip.background = ggplot2::element_blank(),
              strip.text = ggplot2::element_text(hjust = 0, face = "bold"),
              axis.title.y.right = ggplot2::element_text(vjust = 0.5, angle = 90),
              panel.spacing = ggplot2::unit(1.2, "lines"),
              axis.text.x = ggplot2::element_text(size = x_text_size,
                                         angle = x_text_angle,
                                         hjust = x_text_hjust,
                                         vjust = x_text_vjust)
        )
    #grey box for nowcast
    if (infection_nowcast) {
        p <- p +   ggplot2::geom_vline(xintercept = nowcast_start, linetype = "dashed", colour = "grey60") +
            ggplot2::annotate("rect",
                     xmin = nowcast_start,
                     xmax = max(df$date),
                     ymin = -Inf,
                     ymax = Inf,
                     fill = grey(0.5), alpha = 0.1)
    }

    #optional date rug for shorter plots
    if (length(unique(df$date)) < 90) {
        p <- p +
            ggplot2::geom_rug(
                ggplot2::aes(date),
                sides = "b",
                alpha = 1,
                size = 0.5,
                colour = grey(0.5)
            )
    }

    #add validation data as dots
    if (!is.null(case_validation_data) & type == 'reff') {
        stop("Error: cannot overlay case counts over reff, different units!")
    }

    if (!is.null(case_validation_data) & type != 'reff') {
        p <- p + ggplot2::geom_point(data = case_validation_data |>
                                dplyr::filter(date >= start_date,
                                              date <= end_date),
                                ggplot2::aes(x = date, y = count),
                            inherit.aes = TRUE,
                            size = 0.3)
    }


    #fix reff plot ylim
    if (type == 'reff') {
        # ylim <- c(min(df$ci_90_lo), max(df$ci_90_hi))
        p <- p + ggplot2::coord_cartesian(ylim = reff_ylim) +
            ggplot2::geom_hline(yintercept = 1, linetype = 'dotted', colour = "black")

    }

    pdf(filename, width = 10, height = 12)
    print(p)
    dev.off()

    return(filename)
}
