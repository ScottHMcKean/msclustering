#' Run the Gaussian Mixture Model function
#'
#' @description
#' Run the MClust package with up to 30 clusters and a VVV model
#'
#' @param ms_df microseismic data in a dataframe
#' @param comp_df completions data
#' @param surv_df survey day
#' @param run run number
#' @param comp_tcol string to identify col of df with start day
#' @param ms_feat a list of columns for the clustering
#' @param ms_tcol string to identify time column of ms data
#' @param output_path output path to save csvs to
#' @param plot bool for calling stage based plot function
#' @export
run_gmm_parallel <- function(ms_df, comp_df, surv_df, run = 1,
                    comp_tcol = 't',
                    ms_feat = c('x','y','z','t'),
                    ms_tcol = 't',
                    output_path = '../output/',
                    plot = FALSE,
                    ...
                    ) {

  ellipse_df <- data.frame("name" = character(), "x0" = numeric(),"y0" = numeric(),
                           "z0" = numeric(), "pole_easting" = numeric(),
                           "pole_northing" = numeric(), "pole_up" = numeric(),
                           "plane_strike" = numeric(),"plane_dip" = numeric(),
                           "plane_l1" = numeric(),"plane_l2" = numeric(),
                           "plane_l3" = numeric(), "well_num_i" = numeric(),
                           "stage_num_i" = numeric(), "hf" = logical())

  cov_df <- data.frame('c11'= numeric(), 'c12'= numeric(), 'c13'= numeric(),
                       'c14'= numeric(), 'c21'= numeric(), 'c22'= numeric(),
                       'c23'= numeric(), 'c24'= numeric(), 'c31'= numeric(),
                       'c32'= numeric(), 'c33'= numeric(), 'c34'= numeric(),
                       'c41'= numeric(), 'c42'= numeric(), 'c43'= numeric(),
                       'c44'= numeric(), 'x_coord'= numeric(), 'y_coord'= numeric(),
                       'z_coord'= numeric(), 'time'= numeric(), 'name' = character(),
                       'scale_x'= numeric(), 'scale_y'= numeric(), 'scale_z'= numeric(),
                       'scale_t'= numeric())

  ms_df$cluster = "0"

  #sort comp_df by start time for image indexing
  comp_sort <- arrange(comp_df, get(comp_tcol))
  comp_sort$row <- row.names(comp_sort)

  for (well_num_i in unique(ms_df$class_well)){
    # Subset data to a single well and stage and select attributes in the data
    ## Process and subset data for plotting and analysis
    well_ms <- subset(ms_df, class_well == well_num_i)

    print(paste("well",well_num_i))
    print(table(ms_df$class_stage))

    #for (stage_num_i in c(2)){
    for (stage_num_i in unique(well_ms$class_stage)){
      print(paste("well",well_num_i,"stage",stage_num_i))

      ## Set a prefix to label all outputs and allow multiple data analyses
      prefix <- paste('GMM_',well_num_i,'_',stage_num_i,sep="")
      order <- as.character(comp_sort %>% filter(well_num == well_num_i & stage_num == stage_num_i) %>% select(row))

      ## Subset to a single stage
      stage_ms <- well_ms[well_ms$class_stage == stage_num_i,]

      ## Divide data into hydraulic fracture and outliers
      stage_hf_ms = stage_ms[stage_ms$bool == 0,ms_feat]
      stage_out_ms = stage_ms[stage_ms$bool == 1,ms_feat]

      # CLUSTER HF DATA
      if (nrow(stage_hf_ms) > 5){

        ## Select attributes
        stage_hf_ms_scaled <- scale(stage_hf_ms)

        ## Get scaling factors
        data_centers <- attr(stage_hf_ms_scaled,"scaled:center")
        data_scale <- attr(stage_hf_ms_scaled,"scaled:scale")

        # Run the gmm analysis
        gmm_res <- Mclust(stage_hf_ms_scaled, G = 1:25, modelNames = c("VVV"))

        # Assign clusters to overall results
        stage_ms[stage_ms$bool == 0,'cluster'] <- paste("HF",as.character(gmm_res$classification),sep = "")

        ms_df[(ms_df$class_well == well_num_i) & (ms_df$class_stage == stage_num_i) &
                (ms_df$bool == 0),'cluster'] <- paste(well_num_i,
                                                    "_",stage_num_i,
                                                    "_HF",
                                                    as.character(gmm_res$classification),sep = "")

        if (plot == TRUE) {
        gg_cluster_plot(stage_ms[stage_ms$bool == 0,], comp_df, surv_df,
                        well = well_num_i, stage = stage_num_i, outlier = 0, order,
                        title = title_bool, legend_bool = legend_bool)
        }

        # Assign means and covariances
        gmm_means <- gmm_res$parameters$mean
        gmm_covs <- gmm_res$parameters$variance$sigma

        # Set i
        for (i in 1:gmm_res$G){
          # name cluster
          name <- paste("W",well_num_i,"_S",stage_num_i,"_HF_C",i,sep="")
          hf_bool <- TRUE

          # Example on first mean and covariance
          # Extract 3D coordinates of cluster
          cov_i_scaled <- gmm_covs[,,i]
          mean_i_scaled <- gmm_means[,i]

          # Calculate unscaled mean
          mean_i <- mean_i_scaled * data_scale + data_centers

          # Calculate unscaled covariance matrix
          cov_i <- cov_i_scaled * data_scale
          cov_i_array <- as.numeric(cov_i)
          eigen_i <- eigen(cov_i)

          # Eigenvectors and covariance are unit values (direction cosines)
          eigvec <- eigen_i$vectors[1:3,1:3]
          eigval <-  2 * sqrt(9.488) * eigen_i$values[1:3]

          # indices largest and smallest eigenvector based on eigenvalue
          i_max <- which(eigval == max(eigval))
          i_min <- which(eigval == min(eigval))
          i_mid <- c(1,2,3)[c(-i_max,-i_min)]

          # Size of ellipsoid axes
          l1 <- sqrt(sum(eigvec[,i_max]^2)) * eigval[i_max]
          l2 <- sqrt(sum(eigvec[,i_mid]^2)) * eigval[i_mid]
          l3 <- sqrt(sum(eigvec[,i_min]^2)) * eigval[i_min]

          u <- eigvec[,i_max] / sqrt(sum(eigvec[,i_max]^2)) # unit vector of largest eigenvector
          v <- eigvec[,i_mid] / sqrt(sum(eigvec[,i_mid]^2)) # unit vector of second eigenvector
          pole <- c((u[2]*v[3] - u[3]*v[2]),(u[3]*v[1] - u[1]*v[3]),(u[1]*v[2]- u[2]*v[1]))

          # calculate strike and dip of plane from pole
          ## calculate the trend and plunge of the pole
          pole_tp <- cart2sph(enu2ned(pole))  #convert from enu to ned and spherical coordinates
          plane_sd <- pole2plane(pole_tp[1],pole_tp[2]) # convert to strike/dip of plane

          plane_sd * 180 / pi

          ## create an array to bind to an existing datafame or database
          obs <- data.frame("name" = name, "x0" = mean_i[1],"y0" = mean_i[2],
                            "z0" = mean_i[3], "pole_easting" = pole[1],
                            "pole_northing" = pole[2], "pole_up" = pole[3],
                            "plane_strike" = plane_sd[1],"plane_dip" = plane_sd[2],
                            "plane_l1" = l1,"plane_l2" = l2,"plane_l3" = l3,
                            "well_num_i" = well_num_i, "stage_num_i" = stage_num_i,
                            "cluster_num" = i, 'hf' = hf_bool)

          cov_obs <- c(cov_i_scaled, mean_i, data_scale)
          names(cov_obs) = c('c11', 'c12', 'c13', 'c14', 'c21', 'c22', 'c23', 'c24',
                             'c31', 'c32', 'c33', 'c34', 'c41', 'c42', 'c43', 'c44',
                             'x_coord', 'y_coord', 'z_coord', 'time', 'scale_x',
                             'scale_y', 'scale_z', 'scale_t')
          cov_obs <- data.frame(t(cov_obs))
          cov_obs['name'] = name

          cov_df <- rbind(cov_df,cov_obs)
          ellipse_df <- rbind(ellipse_df,obs)
        }
      } else {print("skipping hydraulic fracture data, insufficient number of points")}

      # CLUSTER OUTLIER DATA
      if (nrow(stage_out_ms) > 5){

        ## Select attributes
        stage_out_ms_scaled <- scale(stage_out_ms)

        ## Get scaling factors
        data_centers <- attr(stage_out_ms_scaled,"scaled:center")
        data_scale <- attr(stage_out_ms_scaled,"scaled:scale")

        # Run the gmm analysis
        gmm_res <- Mclust(stage_out_ms_scaled, G = 1:25, modelNames = c("VVV"))

        # Assign clusters to overall results
        stage_ms[stage_ms$bool == 1,'cluster'] <- paste("OUT",as.character(gmm_res$classification),sep = "")

        ms_df[(ms_df$class_well == well_num_i) & (ms_df$class_stage == stage_num_i) &
                ms_df$bool == 1,'cluster'] <- paste(well_num_i,
                                                        "_",stage_num_i,
                                                        "_OUT",
                                                        as.character(gmm_res$classification),sep = "")

        if (plot == TRUE) {
        gg_cluster_plot(stage_ms[stage_ms$bool == 1,], comp_df, surv_df,
                        well = well_num_i, stage = stage_num_i, outlier = 1, order,
                        title = title_bool, legend_bool = legend_bool)
        }

        # Assign means and covariances
        gmm_means <- gmm_res$parameters$mean
        gmm_covs <- gmm_res$parameters$variance$sigma

        # Set i
        for (i in 1:gmm_res$G){
          # name cluster
          name <- paste("W",well_num_i,"_S",stage_num_i,"_OUT_C",i,sep="")
          hf_bool <- FALSE 
          # Example on first mean and covariance
          # Extract 3D coordinates of cluster
          cov_i_scaled <- gmm_covs[,,i]
          mean_i_scaled <- gmm_means[,i]

          # Calculate unscaled mean
          mean_i <- mean_i_scaled * data_scale + data_centers

          # Calculate unscaled covariance matrix
          cov_i <- cov_i_scaled * data_scale
          cov_i_array <- as.numeric(cov_i)
          eigen_i <- eigen(cov_i)

          # Eigenvectors and covariance are unit values (direction cosines)
          eigvec <- eigen_i$vectors[1:3,1:3]
          eigval <-  2 * sqrt(9.488) * eigen_i$values[1:3]

          # indices largest and smallest eigenvector based on eigenvalue
          i_max <- which(eigval == max(eigval))
          i_min <- which(eigval == min(eigval))
          i_mid <- c(1,2,3)[c(-i_max,-i_min)]

          # Size of ellipsoid axes
          l1 <- sqrt(sum(eigvec[,i_max]^2)) * eigval[i_max]
          l2 <- sqrt(sum(eigvec[,i_mid]^2)) * eigval[i_mid]
          l3 <- sqrt(sum(eigvec[,i_min]^2)) * eigval[i_min]

          u <- eigvec[,i_max] / sqrt(sum(eigvec[,i_max]^2)) # unit vector of largest eigenvector
          v <- eigvec[,i_mid] / sqrt(sum(eigvec[,i_mid]^2)) # unit vector of second eigenvector
          pole <- c((u[2]*v[3] - u[3]*v[2]),(u[3]*v[1] - u[1]*v[3]),(u[1]*v[2]- u[2]*v[1]))

          # calculate strike and dip of plane from pole
          ## calculate the trend and plunge of the pole
          pole_tp <- cart2sph(enu2ned(pole))  #convert from enu to ned and spherical coordinates
          plane_sd <- pole2plane(pole_tp[1],pole_tp[2]) # convert to strike/dip of plane

          plane_sd * 180 / pi

          ## create an array to bind to an existing datafame or database
          obs <- data.frame("name" = name, "x0" = mean_i[1],"y0" = mean_i[2],
                            "z0" = mean_i[3], "pole_easting" = pole[1],
                            "pole_northing" = pole[2], "pole_up" = pole[3],
                            "plane_strike" = plane_sd[1],"plane_dip" = plane_sd[2],
                            "plane_l1" = l1,"plane_l2" = l2,"plane_l3" = l3,
                            "well_num_i" = well_num_i, "stage_num_i" = stage_num_i,
                            "cluster_num" = i, 'hf' = hf_bool)

          cov_obs <- c(cov_i_scaled, mean_i, data_scale)
          names(cov_obs) = c('c11', 'c12', 'c13', 'c14', 'c21', 'c22', 'c23', 'c24',
                             'c31', 'c32', 'c33', 'c34', 'c41', 'c42', 'c43', 'c44',
                             'x_coord', 'y_coord', 'z_coord', 'time', 'scale_x',
                             'scale_y', 'scale_z', 'scale_t')
          cov_obs <- data.frame(t(cov_obs))
          cov_obs['name'] = name

          cov_df <- rbind(cov_df,cov_obs)
          ellipse_df <- rbind(ellipse_df,obs)
        }
      } else {print("skipping outlier data, insufficient number of points")}
    }
  }

  write.csv(ellipse_df,paste0(output_path,"run_",run,"_clusters.csv"), row.names = FALSE)
  write.csv(cov_df,paste0(output_path,"run_",run,"_covariance.csv"),row.names = FALSE)
  write.csv(ms_df,paste0(output_path,"run_",run,"_gmm_results.csv"),row.names = FALSE)

  return(list(ms_df, ellipse_df, cov_df))
}

