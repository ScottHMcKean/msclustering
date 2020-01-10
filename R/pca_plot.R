#' Calculate and Plot a PCA unit circle
#' @param num_ml_df clean, numeric only dataframe of ml predictors
#' @param output_path output path for plot
#' @return ggplot
#' @export
pca_unit_circle_plot <- function(num_df, output_path, colour_vec = NULL){
  
  pca <- prcomp(scale(num_df))
  
  correlations = as.data.frame(cor(num_df,pca$x))
  
  # draw unit circle
  tt = seq(0, 2 * pi, length = 100)
  circle <- data.frame(x= 1 * cos(tt), y = 1 * sin(tt))
  
  # draw PCA arrows
  arrows <- data.frame(x1 = rep(0,nrow(correlations)), 
                       y1 = rep(0,nrow(correlations)),
                       x2 = correlations$PC1,
                       y2 = correlations$PC2)
  
  # scale PCA results to +/- 1 to fit on unit circle plot
  range <- apply(pca$x, 2, range)
  
  # pull coordinates of PCA
  pca_results <- as.data.frame(scale(pca$x, center = TRUE, scale = abs(range[1,])+abs(range[2,]))) 
  
  # custom ggplot of PCA results and unit circle
  plot <- ggplot() +
    geom_hline(yintercept = 0, colour = 'gray') +
    geom_vline(xintercept = 0, colour = 'gray')
  
  if(is.null(colour_vec)){
    plot <- plot +
      geom_point(data = pca_results, aes(x = PC1, y = PC2), alpha = 0.5)
  } else {
    plot <- plot +
      geom_point(data = pca_results, 
                 aes(x = PC1, y = PC2, colour = as.factor(colour_vec)), 
                 alpha = 0.5) +
      scale_color_brewer(palette = 'Set1', name = 'cluster')
  }
  
  plot +
    geom_path(data = circle, aes(x = x, y = y), colour = "gray65") +
    geom_segment(data = arrows, 
                 aes(x = x1, y = y1, xend = x2, yend = y2), colour = "gray65") +
    geom_text(data = correlations, 
              aes(x = PC1, y = PC2, label = rownames(correlations)), 
              colour = 'black', size = 2) +
    xlim(-1.1, 1.1) + 
    ylim(-1.1, 1.1) +
    coord_fixed() +
    ggtitle("PCA Correlation Circle") +
    theme_minimal() +
    ggsave(paste(output_path,"pca_circle.jpeg",sep=""), width = 6, height = 4, units = 'in')
}

