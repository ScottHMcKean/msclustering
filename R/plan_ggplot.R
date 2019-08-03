#' Define a ggplot function to plot the data on an x-y views only
#' With an inset legend, for publication (no title)
#'
#' @param ms_df microseismic dataframe with x,y,z,m0
#' @param comp_df completions dataframe
#' @param surv_df survey dataframe
#' @param well well label
#' @param stage stage label
#' @param zoom boolean for zoomed plot
#' @param title char for title string
#' @param label char to file name label
#' @export
plan_ggplot <- function(ms_df, comp_df, surv_df, well = "3", stage = "2",
                        title = '(a)', zoom = FALSE, label = ""){

  ms_df$bool <- as.numeric(ms_df$bool)

  xyplot <- ggplot(ms_df) +
    geom_point(aes(x = x, y = y, size = m0,
                   colour = factor(cluster), fill = factor(cluster)), alpha = 0.75) +
    stat_ellipse(aes(x = x, y = y,
                     colour = factor(cluster)), type = 'norm') +
    stat_ellipse(aes(x = x, y = y,
                     fill = factor(cluster)), geom = "polygon", type = 'norm', alpha = 0.25) +
    scale_size_continuous(name = "Magnitude", guide = 'none') +
    scale_colour_npg(name = "Cluster") +
    scale_fill_npg(name = "Cluster")

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
    theme(legend.position="bottom")

  if (zoom == TRUE) {
    xyplot <- xyplot +
      ylim(500,6500) +
      xlim(500,6500)
  } else {
    xyplot <- xyplot +
      ylim(500,6500) +
      xlim(500,6500)
  }

  ggsave(file = paste(label,"W",well,"_","S",stage,"_plan.jpg",sep=""),
         xyplot, width = 15, height = 15, units = 'cm', dpi = 320)
}
