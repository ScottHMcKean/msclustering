#' Define a ggplot function to plot the data on an x-y, y-z, and z-x views
#'
#' @param ms_df microseismic dataframe with x,y,z
#' @param well well label
#' @param stage stage label
#' @param outlier outlier flag
#' @param order order label
#' @param title title boolean
#' @param legend_bool legend plot boolean
#' @param label optional label for filenames
#' @export
gg_cluster_plot <- function(ms_df, comp_df, surv_df, well = "1",
                            stage = "1",outlier = 0, order = '0',
                            title = TRUE, legend_bool = FALSE,
                            label = "", output_path = "../output/"){

  if (title) {
    if (outlier == 0){
      title <- paste("(a)", "Well", well, "Stage", stage, "(Hydraulic Fracture)")
    } else {
      title <- paste("(a)", "Well", well, "Stage", stage, "(Outliers)")
    }
  } else {
    title <- "(a)"
  }

  xyplot <- ggplot(ms_df) +
    geom_point(aes(x = x, y = y,
                   colour = factor(cluster), fill = factor(cluster)), alpha = 0.75) +
    stat_ellipse(aes(x = x, y = y,
                     colour = factor(cluster)), type = 'norm') +
    stat_ellipse(aes(x = x, y = y,
                     fill = factor(cluster)), geom = "polygon", type = 'norm', alpha = 0.25) +
    scale_colour_discrete(name = "Cluster") +
    scale_fill_brewer(name = "Cluster")

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
    geom_point(aes(x = x, y = z,
                   colour = factor(cluster), fill = factor(cluster)), alpha = 0.75) +
    stat_ellipse(aes(x = x, y = z,
                     colour = factor(cluster)), type = 'norm') +
    stat_ellipse(aes(x = x, y = z,
                     fill = factor(cluster)), geom = "polygon", type = 'norm', alpha = 0.25) +
    scale_colour_ucscgb(name = "Cluster") +
    scale_fill_ucscgb(name = "Cluster")

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
    geom_point(aes(x = z, y = y,
                   colour = factor(cluster), fill = factor(cluster)), alpha = 0.75) +
    stat_ellipse(aes(x = z, y = y,
                     colour = factor(cluster)), type = 'norm') +
    stat_ellipse(aes(x = z, y = y,
                     fill = factor(cluster)), geom = "polygon", type = 'norm', alpha = 0.25) +
    scale_colour_ucscgb(name = "Cluster") +
    scale_fill_ucscgb(name = "Cluster")

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
    xlim(-2000,-4500) +
    ylim(500,6500) +
    theme(legend.position="none")

  legend_plt <- ggplot(ms_df) +
    geom_point(aes(x = x, y = y, colour = factor(cluster))) +
    scale_colour_ucscgb(name = "Cluster")

  if (legend_bool) {
    legend <- get_legend(legend_plt)
    g = arrangeGrob(xyplot,zyplot,xzplot,legend, ncol = 2, nrow = 2, widths = c(3,1), heights = c(3,1))
  } else {
    g = arrangeGrob(xyplot,zyplot,xzplot, ncol = 2, nrow = 2, widths = c(3,1), heights = c(3,1))
  }

  if (outlier == 0){
    ggsave(file = paste(output_path,order,label,"_W",well,"_","S",stage,"_HF.png",sep=""),
           g, width = 20, height = 20, units = 'cm', dpi = 320)
  } else {
    ggsave(file = paste(output_path,order,label,"_W",well,"_","S",stage,"_OUT.png",sep=""),
           g, width = 20, height = 20, units = 'cm', dpi = 320)
  }
}
