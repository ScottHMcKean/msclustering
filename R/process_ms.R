#' Function to process raw microseismic data into a clean data.frame for subsequent analysis
#' @description function takes a raw ms dataframe with projected coordinates (e.g. UTM NAD 83)
#' and processes an analysi
#' @param ms_df raw microseismic dataframe
#' @param mag_complete moment magnitude of completeness of the dataset
#' @param ms_feat a list of features to carry into the cluster analysis
#' @param crs projection of the data
#' @param names of required columns - character vector of time, x, y, z, and magnitude
#' @return data.frame with t,x,y,z,and any other features
process_ms <- function(ms_df, mag_completness = -2,
                       loc_cols = c('x','y','z'),
                       date_time_cols = c('date','time'),
                       mag_col = c('mw'),
                       feat_cols = NA,
                       crs = 26911,
                       tz = 'UTC') {

  loc_cols = c('x_nad83','y_nad83','z')
  # get proper date - time
  if (length(date_time_cols == 1)){
    date = ms_df %>%
      pull(date_time_cols[1]) %>%
      lubridate::as_datetime(., tz = tz)
  } else {
    date = ms_df %>%
      mutate(date = date_time_cols[1] %>% ymd(tz = tz)) %>%
      mutate(time = date_time_cols[2] %>% ymd_hms(foo$start.time)())
  }

  # # assign locations to an sf dataframe with controlled column names
  # coordinates <- ms_df %>%
  #   dplyr::select(x = loc_cols[1], y = loc_cols[2], z = loc_cols[3])
  #
  # centroid <- colMeans(coordinates)
  #
  # rot = function(df, a, coords = c('x','y')){
  #   df[, coords] = df[, coords] * matrix(c(cos(a), sin(a), -sin(a), cos(a)), 2)
  #   df
  # }
  #
  # rotate_xy = function(a, x, y){
  #
  # }
  #
  # coordinates
  # centroid
  #
  # centered <- map2(.x = coordinates, .y = centroid, .f = function(.x,.y){.x - .y}) %>%
  #   as.data.frame %>%
  #   st_as_sf(coords = c('x','y','z'), remove = FALSE, crs = 26911) %>%
  #   mutate(c(x,y) = c(x,y) * rot_mat(pi/4))
  #
  #
  # centered

#
#
#   centered
#   %>%
#     rot(df = ., -pi/2, coords = c('x','y')) %>%
#     st_as_sf(coords = c('x','y','z'), remove = FALSE, crs = 26911) %>%
#     pull(geometry) %>%
#     plot()
#
#   plot(rotated$geometry)
#   # rotate and scale
#   out_df <- (out_df$geometry - st_centroid(out_df$geometry)) * rot(pi/4)
#
#   centroid <- colMeans(out_df %>% st_coordinates())

out_df

  # convert dates to POSIX format

}
