#' Define a ggplot function to plot the data on an x-y, y-z, and z-x views
#' Tailored for the the detection outlier aglo
#'
#' @param ms_df microseismic dataframe with x,y,z,m0
#' @param comp_df completions dataframe
#' @param surv_df survey dataframe
#' @param well well label
#' @param stage stage label
#' @param circle_data shapiro circle data
#' @param title boolean for title
#' @param legend_bool boolean for legend
#' @export
outlier_ggplot <- function(ms_df, comp_df, surv_df, circle_data, well = "1",
                           stage = "1", title = TRUE, legend_bool = FALSE) {

  ms_df$bool <- as.numeric(ms_df$bool)

  if (title) {
    title <- paste("(a)", "Well", well, "Stage", stage, "(Hydraulic Fracture)")
  } else {
    title <- "(a)"
  }

  xyplot <- ggplot(ms_df) +
    geom_point(aes(x = x, y = y, size = m0,
                   colour = factor(bool), fill = factor(bool)), alpha = 0.75) +
    geom_path(data = circle_data, aes(x = x, y = y)) +
    scale_size_continuous(name = "Magnitude") +
    scale_colour_npg(name = "Outlier") +
    scale_fill_npg(name = "Outlier")

  for (well_i in seq(length(unique(comp_df$well)))){
    xyplot <- xyplot +
      geom_path(data = surv_df[surv_df$well == well_i,], aes(x=x,y=y)) +
      geom_point(data = comp_df[comp_df$well == well_i,], aes(x=x,y=y),
                 col ='red', size = 1)
  }

  xyplot <- xyplot +
    ggtitle(title) +
    xlab('Easting (m)') +
    ylab('Northing (m)') +
    theme_minimal() +
    ylim(500,6500) +
    xlim(500,6500) +
    theme(legend.position="none")

  xzplot <- ggplot(ms_df) +
    geom_point(aes(x = x, y = z, size = m0,
                   colour = factor(bool), fill = factor(bool)), alpha = 0.75) +
    geom_path(data = circle_data, aes(x = x, y = z)) +
    scale_size_continuous(name = "Magnitude") +
    scale_colour_npg(name = "Outlier") +
    scale_fill_npg(name = "Outlier")

  for (well_i in seq(length(unique(comp_df$well)))){
    xzplot <- xzplot +
      geom_path(data = surv_df[surv_df$well == well_i,], aes(x=x,y=z)) +
      geom_point(data = comp_df[comp_df$well == well_i,], aes(x=x,y=z),
                 col ='red', size = 1)
  }

  xzplot <- xzplot +
    ggtitle('(c)') +
    xlab('Easting (m)') +
    ylab('Elevation (m)') +
    theme_minimal() +
    ylim(-4500,-2000) +
    xlim(500,6500) +
    theme(legend.position="none")

  zyplot <- ggplot(ms_df) +
    geom_point(aes(x = z, y = y, size = m0,
                   colour = factor(bool), fill = factor(bool)), alpha = 0.75) +
    geom_path(data = circle_data, aes(x = z, y = yz)) +
    scale_size_continuous(name = "Magnitude") +
    scale_colour_npg(name = "Outlier") +
    scale_fill_npg(name = "Outlier")

  for (well_i in seq(length(unique(comp_df$well)))){
    zyplot <- zyplot +
      geom_path(data = surv_df[surv_df$well == well_i,], aes(x=z,y=y)) +
      geom_point(data = comp_df[comp_df$well == well_i,], aes(x=z,y=y),
                 col ='red', size = 1)
  }

  zyplot <- zyplot +
    ggtitle('(b)') +
    xlab('Elevation (m)') +
    ylab('Northing (m)') +
    theme_minimal() +
    xlim(-4500,-2000) +
    ylim(500,6500) +
    theme(legend.position="none")

  legend_plt <- ggplot(ms_df) +
    geom_point(aes(x = x, y = y, colour = factor(bool))) +
    scale_colour_npg(name = "Outlier")

  if (legend_bool) {
    legend <- get_legend(legend_plt)
    g = arrangeGrob(xyplot,zyplot,xzplot,legend, ncol = 2, nrow = 2, widths = c(3,1), heights = c(3,1))
  } else {
    g = arrangeGrob(xyplot,zyplot,xzplot, ncol = 2, nrow = 2, widths = c(3,1), heights = c(3,1))
  }

  ggsave(file = paste("W",well,"_","S",stage,"_DET.jpg",sep=""),
           g, width = 20, height = 20, units = 'cm', dpi = 320)
}
