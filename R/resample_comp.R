#' Function to resample location of hydraulic fracturing
#' @description resamples completions dataframe with uniform distribution
#' @param survey_df well survey
#' @param comp_df completions dataframe
#' @return resampled completions data frame
#' @export
resample_comp_df <- function(comp_df, survey_df) {
  well_numbers <- survey_df %>%
    pull(well_num) %>%
    unique()

  # load and resample completions
  comp_df <- comp_df %>%
    mutate(md_start = md - perf_length) %>%
    mutate(resample_md = runif(nrow(.), min = md_start, max = md)) %>%
    mutate(unc_x = runif(nrow(.), min = -unc_md, max = unc_md)) %>%
    mutate(unc_y = runif(nrow(.), min = -unc_md, max = unc_md)) %>%
    mutate(unc_z = runif(nrow(.), min = -unc_md, max = unc_md))

  for (current_well in well_numbers) {
    well_survey_df <- survey_df %>%
      filter(well_num == current_well)

    spline_x <- splinefun(well_survey_df$md, well_survey_df$x)
    spline_y <- splinefun(well_survey_df$md, well_survey_df$y)
    spline_z <- splinefun(well_survey_df$md, well_survey_df$z)

    well_comp_df <- comp_df %>%
      filter(well_num == current_well) %>%
      mutate(x = spline_x(resample_md) + unc_x) %>%
      mutate(y = spline_y(resample_md) + unc_y) %>%
      mutate(z = spline_z(resample_md) + unc_z)

    if (current_well == well_numbers[1]){
      resampled_comp_df <- well_comp_df
    } else {
      resampled_comp_df <- rbind(resampled_comp_df, well_comp_df)
    }
  }
  resampled_comp_df
}
