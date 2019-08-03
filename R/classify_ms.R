#' Classification algorithm for microseismic events
#' @description
#' Classify events based on time (events must occur after stage start)
#' to start of stage and proximity to center of stage
#' Events have to have the same column names, corresponding to x_col,y_col, and z_col
#' Times should be in posixct
#'
#' @param ms_df microseismic data in a dataframe
#' @param comp_df completions data
#' @param comp_time_col string to identify col of df with start day
#' @param ms_time_col string to identify col of df with start day
#' @export
classify_ms <- function(ms_df, comp_df,
                        comp_tcol = 't',
                        ms_tcol = 't') {
  ms_df$class_well = 0
  ms_df$class_stage = 0

  t_sd <- ms_df %>% pull(get(ms_tcol)) %>% sd

  ms_df <- ms_df %>%
    mutate(t_scale = (get(ms_tcol) - min(get(ms_tcol))) / t_sd)

  comp_df <- comp_df %>%
    mutate(t_scale = (get(comp_tcol) - min(get(comp_tcol))) / t_sd)

  ms_sd <- ms_df %>%
    summarize(x = sd(x), y = sd(y), z = sd(z))

  ms_mean <- ms_df %>%
    summarize(x = mean(x), y = mean(y), z = mean(z))

  ms_df <- ms_df %>%
    mutate(x_scale = (x - ms_mean$x)/ms_sd$x) %>%
    mutate(y_scale = (y - ms_mean$y)/ms_sd$y) %>%
    mutate(z_scale = (z - ms_mean$z)/ms_sd$z)

  comp_df <- comp_df %>%
    mutate(x_scale = (x - ms_mean$x)/ms_sd$x) %>%
    mutate(y_scale = (y - ms_mean$y)/ms_sd$y) %>%
    mutate(z_scale = (z - ms_mean$z)/ms_sd$z)

  # pick first point
  for (i in 1:nrow(ms_df)) {
    event = ms_df[i,]

    # positive time after stage start in minutes
    time_dist = (event$t - comp_df$t)*1440

    # events that occured before the start of a stage are excluded by imposing a large distance
    time_dist[time_dist < 0] <- 1E9

    # euclidean spatial distance
    euc_dist = sqrt((event$x - comp_df$x)**2 +
                      (event$y - comp_df$y)**2 +
                      (event$z - comp_df$z)**2)
    
    # add distances and find minimum to classify each stage
    net_dist = sqrt(time_dist**2 + euc_dist**2)

    if(sum(time_dist <= 1E5) == 0){
      ms_df[i,'class_well'] = 0
      ms_df[i,'class_stage'] = 0
    } else {
      ms_df[i,'class_well'] = comp_df[net_dist == min(net_dist),]$well_num
      ms_df[i,'class_stage'] = comp_df[net_dist == min(net_dist),]$stage_num
    }
  }
  return(ms_df)
}

#' Confusion matrix plot
#' @description
#' Create and plot a confusion matrix
#'
#' @param ms_df microseismic data in a dataframe
#' @param comp_df completions data
#' @param comp_start_col string to identify col of df with start day
#' @export
confusion_matrix <- function(ms_df) {
  ms_df$op_class = ms_df$well * 100 + ms_df$stage
  ms_df$ml_class = ms_df$class_well * 100 + ms_df$class_stage

  reclass = sum(ms_df$op_class != ms_df$ml_class)

  confusion_matrix <- as.data.frame(table(factor(ms_df$ml_class), factor(ms_df$op_class)))

  # plot confusion matrix
  ggplot(data = confusion_matrix, aes(x = Var1,y = Var2)) +
    geom_tile(aes(fill = Freq)) +
    scale_fill_gradient(low = 'blue',
                        high = 'red',
                        trans = 'log',
                        name = 'Log Frequency') +
    ylab('Operator Classification') +
    xlab('Algorithm Classification') +
    theme_minimal() +
    coord_fixed() +
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank()) +
    ggsave('classification.jpeg', height = 15, width = 15, units = 'cm', dpi = 320) +
    ggsave('classification.eps', height = 15, width = 15, units = 'cm')

  # plot reclassified points
  ggplot(ms_df) +
    geom_point(aes(x = x, y = y, colour = factor((op_class != ml_class))), size = 0.2, shape = 1) +
    scale_color_discrete(h.start = 150, name = 'reclassified', direction = -1) +
    xlim(0,5000) +
    ylim(0,5000) +
    coord_fixed() +
    xlab('Easting') +
    ylab('Northing') +
    theme_minimal() +
    ggsave('point_reclass.jpeg', height = 15, width = 15, units = 'cm', dpi = 320) +
    ggsave('point_reclass.eps', height = 15, width = 15, units = 'cm')
}
