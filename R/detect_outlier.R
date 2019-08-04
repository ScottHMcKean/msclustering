#' Outlier & Induced Event Detection
#' @description
#' Now we investigate each stage using an R-T plot to see if we can remove outliers.
#' This approach follows Shapiro (2008) and Shapiro (2015), where an apparent
#' diffusivity is estimated based on the approximation to the PKN model of fracture
#' growth. See McKean et al. (2019).
#'
#' @param ms_df microseismic dataframe
#' @param comp_df completions dataframe
#' @param surv_df survey dataframe
#' @param plot boolean to activate plots or not
#' @export
detect_outlier <- function(ms_df, comp_df, plot = FALSE) {

  if (!('bool' %in% colnames(ms_df))) {
    ms_df$bool <- 0
  }

  for (well_num_i in unique(comp_df$well_num)) {
    stages = comp_df %>%
      filter(well_num == well_num_i) %>%
      pull(stage_num) %>%
      unique()

    for (stage_num_i in stages){
      stage = comp_df %>%
        filter(well_num == well_num_i & stage_num == stage_num_i)

      stage_ms = ms_df %>%
        filter(class_well == well_num_i & class_stage == stage_num_i)
      
      if (nrow(stage_ms) == 0){
        next
      }

      ms_dist = sqrt((stage_ms$x - stage$x)**2 + 
                       (stage_ms$y - stage$y)**2 +
                       (stage_ms$z - stage$z)**2)

      bool = (ms_dist - stage$r_ap) >= 0

      ms_df[(ms_df$class_well == well_num_i) & (ms_df$class_stage == stage_num_i),'bool'] <- bool

      # if (plot & (nrow(stage_ms) >= 5)) {
      #   # filename
      #   f_name = paste('outlier_W',well_num_i,'_S',stage_num_i,'.jpeg',sep = "")
      # 
      #   outlier_ggplot(stage_ms, comp_df, surv_df, circle_data = dat, well = well_num_i,
      #                  stage = stage_num_i, title = FALSE, legend_bool = TRUE)
      #  }
    }
  }
  return(ms_df)
}

#' create a circle around the stage center
#'
#' @param center vector
#' @param r radius
#' @param npoints int, number of points
#' @return dataframe
#' @export
circle3d <- function(center = c(0,0,0),r = 1, npoints = 100){
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  zz <- center[3] + r * sin(tt)
  yz <- center[2] + r * cos(tt)
  return(data.frame(x = xx, y = yy, z = zz, yz = yz))
}
