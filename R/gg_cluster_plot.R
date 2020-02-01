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
gg_stage_plot <- function(ms_df, comp_df, surv_df, well = 1,
                            stage = 1, outlier = 0,
                            title = TRUE, legend_bool = FALSE,
                            label = "", output_path = "../output/",
                            elev_lim = c(-3000,-2000),
                            north_lim = NULL,
                            east_lim = NULL) {

  if (title) {
    if (outlier == 0){
      title <- paste("Well", well, "Stage", stage, "(HF-Created Seismicity)")
    } else {
      title <- paste("Well", well, "Stage", stage, "(Induced Seismicity)")
    }
  }

  xyplot <- ggplot(ms_df) +
    geom_point(aes(x = x, y = y,
                   colour = factor(cluster), fill = factor(cluster)), alpha = 0.75) +
    stat_ellipse(aes(x = x, y = y,
                     colour = factor(cluster)), type = 'norm', level = 0.68) +
    stat_ellipse(aes(x = x, y = y,
                     fill = factor(cluster)), geom = "polygon", type = 'norm', 
                 alpha = 0.25, level = 0.68) +
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
    xlim(east_lim[1],east_lim[2]) + 
    ylim(north_lim[1],north_lim[2]) + 
    theme_minimal() +
    theme(legend.position="none") +
    coord_equal()

  xzplot <- ggplot(ms_df) +
    geom_point(aes(x = x, y = z,
                   colour = factor(cluster), fill = factor(cluster)), alpha = 0.75) +
    stat_ellipse(aes(x = x, y = z,
                     colour = factor(cluster)), type = 'norm', level = 0.68) +
    stat_ellipse(aes(x = x, y = z,
                     fill = factor(cluster)), geom = "polygon", type = 'norm', 
                 alpha = 0.25, level = 0.68) +
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
    xlim(east_lim[1],east_lim[2]) +
    ylim(elev_lim[1],elev_lim[2]) +
    theme(legend.position="none")

  zyplot <- ggplot(ms_df) +
    geom_point(aes(x = z, y = y,
                   colour = factor(cluster), fill = factor(cluster)), alpha = 0.75) +
    stat_ellipse(aes(x = z, y = y,
                     colour = factor(cluster)), type = 'norm', level = 0.68) +
    stat_ellipse(aes(x = z, y = y,
                     fill = factor(cluster)), geom = "polygon", type = 'norm', 
                 alpha = 0.25, level = 0.68) +
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
    xlim(elev_lim[1],elev_lim[2]) + 
    ylim(north_lim[1],north_lim[2]) + 
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
    ggsave(file = paste(output_path,label,"_W",well,"_","_HF.tiff",sep=""),
           g, width = 20, height = 20, units = 'cm', dpi = 320)
    ggsave(file = paste(output_path,label,"_W",well,"_","_HF.jpg",sep=""),
           g, width = 20, height = 20, units = 'cm', dpi = 320)
  } else {
    ggsave(file = paste(output_path,label,"_W",well,"_","_OUT.tiff",sep=""),
           g, width = 20, height = 20, units = 'cm', dpi = 320)
    ggsave(file = paste(output_path,label,"_W",well,"_","_OUT.jpg",sep=""),
           g, width = 20, height = 20, units = 'cm', dpi = 320)
  }
}


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
gg_well_plot_tif <- function(ms_df, comp_df, surv_df, well = 1,
                         outlier = 0, title_bool = TRUE, label = "", 
                         output_path = "../output/", elev_lim = NULL,
                         north_lim = NULL, east_lim = NULL){
  
  if (title_bool) {
    if (outlier == 0){
      title <- paste("Well", well, "(HF-Created Seismicity)")
    } else {
      title <- paste("Well", well, "(Induced Seismicity)")
    }
  } else {
    title <- ""
  }
  
  xyplot <- ggplot(ms_df) +
    geom_point(aes(x = x, y = y,
                   colour = factor(cluster), fill = factor(cluster)), alpha = 0.75) +
    stat_ellipse(aes(x = x, y = y,
                     colour = factor(cluster)), type = 'norm', level = 0.685) +
    stat_ellipse(aes(x = x, y = y,
                     fill = factor(cluster)), geom = "polygon", type = 'norm', 
                 alpha = 0.25, level = 0.68) +
    scale_color_viridis_d(name = "Cluster") +
    scale_fill_viridis_d(name = "Cluster")
  
  for (well_i in seq(length(unique(comp_df$well)))){
    xyplot <- xyplot +
      geom_path(data = surv_df[surv_df$well == well_i,], aes(x=x,y=y)) +
      geom_point(data = comp_df[comp_df$well == well_i,], aes(x=x,y=y),
                 col ='red', size = 1)
  }
  
  if (!is.null(east_lim)) {
    xyplot <- xyplot + xlim(east_lim[1],east_lim[2])
  }
  
  if (!is.null(east_lim)) {
    xyplot <- xyplot + ylim(north_lim[1],north_lim[2])
  }
  
  xyplot <- xyplot +
    ggtitle(title) +
    xlab('Easting (m)') +
    ylab('Northing (m)') +
    theme_minimal() +
    theme(legend.position="none")
  
  xzplot <- ggplot(ms_df) +
    geom_point(aes(x = x, y = z,
                   colour = factor(cluster), fill = factor(cluster)), alpha = 0.75) +
    stat_ellipse(aes(x = x, y = z,
                     colour = factor(cluster)), type = 'norm', level = 0.685) +
    stat_ellipse(aes(x = x, y = z,
                     fill = factor(cluster)), geom = "polygon", type = 'norm', 
                 alpha = 0.25, level = 0.68) +
    scale_color_viridis_d(name = "Cluster") +
    scale_fill_viridis_d(name = "Cluster")
  
  for (well_i in seq(length(unique(comp_df$well)))){
    xzplot <- xzplot +
      geom_path(data = surv_df[surv_df$well == well_i,], aes(x=x,y=z)) +
      geom_point(data = comp_df[comp_df$well == well_i,], aes(x=x,y=z),
                 col ='red', size = 1)
  }
  
  if (!is.null(east_lim)) {
    xzplot <- xzplot + xlim(east_lim[1],east_lim[2])
  }
  
  if (!is.null(elev_lim)) {
    xzplot <- xzplot + ylim(elev_lim[1],elev_lim[2])
  }
  
  xzplot <- xzplot +
    ggtitle('(c)') +
    xlab('Easting (m)') +
    ylab('Elevation (m)') +
    theme_minimal() +
    theme(legend.position="none")
  
  zyplot <- ggplot(ms_df) +
    geom_point(aes(x = z, y = y,
                   colour = factor(cluster), fill = factor(cluster)), alpha = 0.75) +
    stat_ellipse(aes(x = z, y = y,
                     colour = factor(cluster)), type = 'norm', level = 0.685) +
    stat_ellipse(aes(x = z, y = y,
                     fill = factor(cluster)), geom = "polygon", type = 'norm', 
                 alpha = 0.25, level = 0.68) +
    scale_color_viridis_d(name = "Cluster") +
    scale_fill_viridis_d(name = "Cluster")
  
  for (well_i in seq(length(unique(comp_df$well)))){
    zyplot <- zyplot +
      geom_path(data = surv_df[surv_df$well == well_i,], aes(x=z,y=y)) +
      geom_point(data = comp_df[comp_df$well == well_i,], aes(x=z,y=y),
                 col ='red', size = 1)
  }
  
  if (!is.null(elev_lim)) {
    zyplot <- zyplot + xlim(elev_lim[1],elev_lim[2])
  }
  
  if (!is.null(north_lim)) {
    zyplot <- zyplot + ylim(north_lim[1],north_lim[2])
  }
  
  zyplot <- zyplot +
    ggtitle('(b)') +
    xlab('Elevation (m)') +
    ylab('Northing (m)') +
    theme_minimal() +
    theme(legend.position="none")
  
  g = arrangeGrob(xyplot, zyplot, xzplot, ncol = 2, nrow = 2, 
                  widths = c(3,1), heights = c(3,1))

  if (outlier == 0){
    ggsave(file = paste(output_path,label,"_W",well,"_","_HF.tiff",sep=""),
           g, width = 20, height = 20, units = 'cm', dpi = 320)
    ggsave(file = paste(output_path,label,"_W",well,"_","_HF.jpg",sep=""),
           g, width = 20, height = 20, units = 'cm', dpi = 320)
  } else {
    ggsave(file = paste(output_path,label,"_W",well,"_","_OUT.tiff",sep=""),
           g, width = 20, height = 20, units = 'cm', dpi = 320)
    ggsave(file = paste(output_path,label,"_W",well,"_","_OUT.jpg",sep=""),
           g, width = 20, height = 20, units = 'cm', dpi = 320)
  }
}

gg_well_plot_eps <- function(ms_df, comp_df, surv_df, well = 1,
                             outlier = 0, label = "", 
                             output_path = "../output/", elev_lim = NULL,
                             north_lim = NULL, east_lim = NULL){
  
  xyplot <- ggplot(ms_df) +
    geom_point(aes(x = x, y = y,
                   colour = factor(cluster), fill = factor(cluster))) +
    stat_ellipse(aes(x = x, y = y,
                     colour = factor(cluster)), type = 'norm', level = 0.685) +
    scale_color_viridis_d(name = "Cluster") +
    scale_fill_viridis_d(name = "Cluster")
  
  for (well_i in seq(length(unique(comp_df$well)))){
    xyplot <- xyplot +
      geom_path(data = surv_df[surv_df$well == well_i,], aes(x=x,y=y)) +
      geom_point(data = comp_df[comp_df$well == well_i,], aes(x=x,y=y),
                 col ='red', size = 1)
  }
  
  if (!is.null(east_lim)) {
    xyplot <- xyplot + xlim(east_lim[1],east_lim[2])
  }
  
  if (!is.null(east_lim)) {
    xyplot <- xyplot + ylim(north_lim[1],north_lim[2])
  }
  
  xyplot <- xyplot +
    ggtitle('') +
    xlab('Easting (m)') +
    ylab('Northing (m)') +
    theme_minimal() +
    theme(legend.position="none")
  
  xzplot <- ggplot(ms_df) +
    geom_point(aes(x = x, y = z,
                   colour = factor(cluster), fill = factor(cluster))) +
    stat_ellipse(aes(x = x, y = z,
                     colour = factor(cluster)), type = 'norm', level = 0.685) +
    scale_color_viridis_d(name = "Cluster") +
    scale_fill_viridis_d(name = "Cluster")
  
  for (well_i in seq(length(unique(comp_df$well)))){
    xzplot <- xzplot +
      geom_path(data = surv_df[surv_df$well == well_i,], aes(x=x,y=z)) +
      geom_point(data = comp_df[comp_df$well == well_i,], aes(x=x,y=z),
                 col ='red', size = 1)
  }
  
  if (!is.null(east_lim)) {
    xzplot <- xzplot + xlim(east_lim[1],east_lim[2])
  }
  
  if (!is.null(elev_lim)) {
    xzplot <- xzplot + ylim(elev_lim[1],elev_lim[2])
  }
  
  xzplot <- xzplot +
    ggtitle('(c)') +
    xlab('Easting (m)') +
    ylab('Elevation (m)') +
    theme_minimal() +
    theme(legend.position="none")
  
  zyplot <- ggplot(ms_df) +
    geom_point(aes(x = z, y = y,
                   colour = factor(cluster), fill = factor(cluster))) +
    stat_ellipse(aes(x = z, y = y,
                     colour = factor(cluster)), type = 'norm', level = 0.685) +
    scale_color_viridis_d(name = "Cluster") +
    scale_fill_viridis_d(name = "Cluster")
  
  for (well_i in seq(length(unique(comp_df$well)))){
    zyplot <- zyplot +
      geom_path(data = surv_df[surv_df$well == well_i,], aes(x=z,y=y)) +
      geom_point(data = comp_df[comp_df$well == well_i,], aes(x=z,y=y),
                 col ='red', size = 1)
  }
  
  if (!is.null(elev_lim)) {
    zyplot <- zyplot + xlim(elev_lim[1],elev_lim[2])
  }
  
  if (!is.null(north_lim)) {
    zyplot <- zyplot + ylim(north_lim[1],north_lim[2])
  }
  
  zyplot <- zyplot +
    ggtitle('(b)') +
    xlab('Elevation (m)') +
    ylab('Northing (m)') +
    theme_minimal() +
    theme(legend.position="none")
  
  g = arrangeGrob(xyplot, zyplot, xzplot, ncol = 2, nrow = 2, 
                  widths = c(3,1), heights = c(3,1))
  
  if (outlier == 0){
    ggsave(file = paste(output_path,label,"_W",well,"_","_HF.eps",sep=""),
           g, width = 20, height = 20, units = 'cm')
  } else {
    ggsave(file = paste(output_path,label,"_W",well,"_","_OUT.eps",sep=""),
           g, width = 20, height = 20, units = 'cm')
  }
}