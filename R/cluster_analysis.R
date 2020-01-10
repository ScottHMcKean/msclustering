#' Clean cluster names to match gmm df
#' @param df cluster or covariance df
#' @param col column to fix
#' @return df with cluster column
#' @export
clean_cluster_name <- function(df, col){
  df$cluster = df %>% pull(col) %>%
    str_remove_all(., 'W') %>%
    str_remove_all(., 'S') %>%
    str_remove_all(., '_C')
  
  df
}

#' Calculate the distance from points in a cluster to its plane, to be mapped
#' @param this_cluster the current cluster
#' @param ms_gmm_results ms_gmm_results
#' @param covariance the covariance dataframe
#' @return vector of distances between each point in a cluster and its plane
#' @export
point_to_plane_distance_cluster <- function(this_cluster, ms_gmm_results, covariance){
  cluster_points = as.matrix(
    ms_gmm_results[ms_gmm_results$cluster == this_cluster, c('x','y','z')]
  )
  
  this_cov = covariance[covariance$cluster == this_cluster,]
  
  centre = as.numeric(this_cov[1, c('x_coord','y_coord','z_coord')])
  
  cov_mat = matrix(
    as.matrix(this_cov[1, c('c11','c12','c13','c21','c22','c23','c31','c32','c33')]),
    nrow = 3
  )
  
  scale = as.numeric(this_cov[1, c('scale_x', 'scale_y', 'scale_z')])
  
  scaled_cov_mat = cov_mat*scale
  
  normal = eigen(scaled_cov_mat)[[2]][,3]
  
  dist = abs((rowSums(sweep(cluster_points,2,normal, FUN = "*")) - sum(centre*normal)))/Norm(normal)
}

#' Calculate the log likelihood of all points in a cluster
#' @param this_cluster the current cluster
#' @param ms_gmm_results ms_gmm_results
#' @param covariance the covariance dataframe
#' @return vector of log likelihood of each point in a cluster
#' @export
cluster_log_likelihood <- function(this_cluster, ms_gmm_results, covariance){
  cluster_points = as.matrix(
    ms_gmm_results[ms_gmm_results$cluster == this_cluster, c('x','y','z')]
  )
  
  this_cov = covariance[covariance$cluster == this_cluster,]
  
  centre = as.numeric(this_cov[1, c('x_coord','y_coord','z_coord')])
  
  cov_mat = matrix(
    as.matrix(this_cov[1, c('c11','c12','c13','c21','c22','c23','c31','c32','c33')]),
    nrow = 3
  )
  
  scale = as.numeric(this_cov[1, c('scale_x', 'scale_y', 'scale_z')])
  
  scaled_cov_mat = cov_mat*scale
  
  -0.5 * (log(det(scaled_cov_mat)) + mahalanobis(cluster_points,centre,scaled_cov_mat) + nrow(scaled_cov_mat)*log(2*pi))
}

#' Analyze clusters for size, flattening, likelihood, and point-to-plane distances
#' @param clusters cluster dataframe
#' @param ms_gmm_results ms_gmm_results dataframe
#' @param covariance the covariance dataframe
#' @return vector of log likelihood of each point in a cluster
#' @export
analyze_clusters <- function(clusters, ms_gmm_results, covariance){
  clusters <- ms_gmm_results %>% 
    dplyr::group_by(cluster) %>% 
    dplyr::summarise(n = n()) %>%
    dplyr::left_join(clusters, ., by = 'cluster') %>%
    dplyr::mutate(
      out_bool = str_detect(clusters$name,'OUT'),
      plane_area = (plane_l1*plane_l2),
      ellipsoid_volume = 4/3*pi*plane_l1*plane_l2*plane_l3,
      mean_axis_vs_n = (plane_l1+plane_l2+plane_l3) /3 / n,
      avg_radius_flattening = ((plane_l1+plane_l2)/2 - plane_l3)/((plane_l1+plane_l2)/2),
      l3_l1_ratio = plane_l3/plane_l1
    )
  
  cluster_distances <- map(.x = clusters$cluster, 
                           .f = point_to_plane_distance_cluster, 
                           ms_gmm_results = ms_gmm_results,
                           covariance = covariance)
  
  cluster_likelihood <- map(.x = clusters$cluster,
                            .f = cluster_log_likelihood, 
                            ms_gmm_results = ms_gmm_results,
                            covariance = covariance)
  
  clusters <- clusters %>%
    dplyr::mutate(
      rmse_point_plane_rmse = sapply(cluster_distances, function(x){sqrt(sum(x^2)/length(x))}),
      mean_log_likelihood = sapply(cluster_likelihood, mean)
    )
  
  clusters$dbscan <- clusters %>%
    dplyr::select(ellipsoid_volume, mean_axis_vs_n,
                  avg_radius_flattening, rmse_point_plane_rmse, 
                  mean_log_likelihood) %>%
    scale() %>%
    dbscan(., eps = 1, minPts = 8) %>% 
    .$cluster
  
  clusters
}
