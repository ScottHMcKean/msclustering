#' Function to determine bin numbers for catalog determination
#' 
#' @description The ECDF function requires binning earthquakes into discrete categories. This
#' complicates our modelling process by turning a continuous distribution into a 
#' discrete categorical series of bins. There are There are several methods for determining
#' this binning: 
#' Sturges, H. (1926) The choice of a class-interval. J. Amer. Statist. Assoc., 21, 65-66.
#' Scott, D.W. (1979) On optimal and data-based histograms. Biometrika, 66, 605-610.
#' Freedman, D. and Diaconis, P. (1981) On this histogram as a density estimator: L2 theory. Zeit. Wahr. ver. Geb., 57, 453-476.
#' Wand, M. P. (1997). Data-based choice of histogram bin width. The American Statistician, 51(1), 59-64.
#' We use three common methods to calculate the deterministic b-value of 
#' the catalogue. These functions are also used in the bootstrap estimation
#' of uncertainity.
#' @param mw_vector a mw_vector of magnitudes
#' @param method a character of methods: 'rice', 'sturges', 'freedman'
#' @return an integer number of bins
#' @export
k_bins <- function(mw_vector, method = 'rice') {
  if (method == 'sturges') {
    bins <- log2(length(mw_vector)) + 1
  } else if (method == 'freedman') {
    bins <- (max(mw_vector) - min(mw_vector)) / (IQR(mw_vector) / length(mw_vector)^(1/3))
  } else if (method == 'rice') {
    bins <- 2 * length(mw_vector) ^ (1/3)
  } else {
    warning('No binning method specified, defaulting to Rice (1951)')
    bins <- 2 * length(mw_vector) ^ (1/3)
  }
  bins %>% ceiling()
}

#' ECDF Calculation using Freedman & Diaconis (1981) approach for probability 
#' mass function for an earthquake catalogue.
#'
#' @param cat earthquake catalogue
#' @param mwcol string denoting catalogue magnitude column
#' @export
calc_ecdf <- function(cat, mwcol = 'mw', k_bins = 10){
  # Create a sorted column of magnitudes
  cat$sort <- sort(cat$mw)
  
  # Calculate bins & centres using Freedman - Diaconis Rule
  bins <- seq(min(cat$mw),max(cat$mw),length.out = k_bins)
  centres <- (bins[2:length(bins)]-bins[1:length(bins)-1])/2+bins[1:length(bins)-1]
  
  # Calculate bin ECDF
  ecdf <- data.frame(centres)
  ecdf$cdf <- (1-ecdf(cat$mw)(centres))*length(cat$mw)
  ecdf$log_cdf <- log10(ecdf$cdf)
  # Calculate number of earthquakes per bin
  ecdf$count <- rowSums(table(cut(cat$mw,bins),cat$mw))
  
  ecdf
}

#' GGPlot function to plot an ECDF with minimal formatting
#'
#' @param ecdf ecdf dataframe with centres, cdf, and log_cdf
#' @export
plot_ecdf <- function(ecdf, out_file = 'ecdf'){
  ggplot(ecdf) +
    geom_point(aes(x = centres, y = cdf), size = 1.5) +
    geom_path(aes(x = centres, y = cdf)) +
    scale_y_log10() +
    ggtitle('ECDF Plot') +
    ylab('Occurences') +
    xlab('Magnitude') +
    theme_minimal() +
    ggsave(paste0(out_file,'.jpg'), width = 6, height = 6, dpi = 600) +
    ggsave(paste0(out_file,'.eps'), width = 6, height = 6)
}

#' Maximum curvature method linear model
#' 
#' @description Function to calculate the maximum curvature Mc, a value (intercept),
#' and b-value using a smoothing spline and predictor method. 
#' 
#' @param ecdf dataframe with ecdf
#' @return dataframe with model output
#' @export
train_maxc_model <- function(ecdf){
  spl_pred <- smooth.spline(ecdf$centres,ecdf$cdf) %>% predict(.,deriv=2)
  mmax_c <- spl_pred$x[which(abs(spl_pred$y) == max(abs(spl_pred$y)))]
  ecdf_filt_maxc <- ecdf[ecdf$centres > mmax_c,]
  
  fit <- lm(data = ecdf_filt_maxc, log_cdf ~ centres)
  
  maxc_res <- data.frame('method' = 'maxc',
                         'mmax_c'= mmax_c, 
                         'a'= fit$coefficients[1],
                         'b'= fit$coefficients[2],
                         'rmse' = yardstick::rmse_vec(fit$fitted.values, ecdf_filt_maxc$log_cdf),
                         'mae' = yardstick::mae_vec(fit$fitted.values, ecdf_filt_maxc$log_cdf),
                         'bin_width' = diff(ecdf$centres) %>% mean())
  
  rownames(maxc_res) <- 'max_curvature'
  
  maxc_res
}

#' Aki (1965) method
#' 
#' @description Function to calculate the maximum curvature Mc, 'a' value (intercept),
#' and b-value using the Aki (1965) method. 
#' 
#' @param ecdf dataframe with ecdf
#' @return dataframe with model output
#' @export
train_aki_model <- function(ecdf, cat){
  # Maximum of histogram method, maximum likelihood method
  Maki <- ecdf[min(which(ecdf$count == max(ecdf$count))),'centres']
  Mbar <- mean(cat$mw[cat$mw > Maki])
  b <- log10(exp(1)) / (Maki - Mbar)
  
  ecdf_filt <- ecdf[ecdf$centres >= Maki,]
  
  loss_func <- function(a,b,ecdf){
    yardstick::rmse_vec(b*ecdf$centres + a, ecdf$log_cdf)
  }
  
  # Use a vector of a values and determine best RMSE / MAE value
  a <- optimize(loss_func, interval = c(0,10), b = b, ecdf = ecdf_filt)$minimum
  
  # store aki results
  aki_res <- data.frame('method' = 'aki',
                        'mmax_c' = Maki,
                        'a' = a, 
                        'b' = b,
                        'rmse' = yardstick::rmse_vec(b*ecdf_filt$centres + a, ecdf_filt$log_cdf),
                        'mae' = yardstick::mae_vec(b*ecdf_filt$centres + a, ecdf_filt$log_cdf),
                        'bin_width' = diff(ecdf_filt$centres) %>% mean())
  
  rownames(aki_res) <- 'aki'
  
  aki_res
}

#' Gutenberg Richter Probability
#' 
#' @description Returns probability of earthquakes above magnitude zero. Can
#' accept a single value of M or a vector.
#' 
#' @param model dataframe with b value
#' @return vector of probabilities
#' @export
gr_prob <- function(model, M = 2){
  10^(model$b*M)
}


#' Goodness of Fit Method
#' 
#' @description Optimization algorithm to fit three dimensional goodness of 
#' fit b-value.
#' 
#' @param ecdf dataframe with ecdf
#' @return dataframe with model output
#' @export
train_gof_model <- function(ecdf){
  loss_func <- function(par,ecdf){
    ecdf_filt <- ecdf[ecdf$centres >= par[1],]
    yardstick::rmse_vec(par[2]*ecdf_filt$centres + par[3], ecdf_filt$log_cdf)
  }
  
  out <- optim(par = c(-1.5,-1,3),
               lower = c(-3, -2, 2),
               upper = c(-1, -0.5, 4),
               fn = loss_func,
               ecdf = ecdf,
               method = 'L-BFGS-B')
  
  ecdf_filt <- ecdf[ecdf$centres >= out$par[1],]
  
  gof_res <- data.frame('method' = 'gof',
                        'mmax_c' = out$par[1],
                        'a' = out$par[3], 
                        'b' = out$par[2],
                        'rmse' = yardstick::rmse_vec(out$par[2]*ecdf_filt$centres + out$par[3], ecdf_filt$log_cdf),
                        'mae' = yardstick::mae_vec(out$par[2]*ecdf_filt$centres + out$par[3], ecdf_filt$log_cdf),
                        'bin_width' = diff(ecdf_filt$centres) %>% mean())
  
  rownames(gof_res) <- 'gof'
  
  gof_res
}

#' Run three models to fit distribution, output parameters, and plot results
#' 
#' @param ecdf dataframe with ecdf
#' @param cat dataframe with earthquake catalogue
#' @return dataframe with model output
#' @export
run_plot_models <- function(ecdf, cat, plot = FALSE){
  maxc_mod <- train_maxc_model(ecdf)
  aki_mod <- train_aki_model(ecdf, cat)
  gof_mod <- train_gof_model(ecdf)
  
  if (plot) {
    M = seq(0,8,0.1)
    res_df <- data.frame(mag = M, maxc = gr_prob(maxc_mod, M), 
                         aki = gr_prob(aki_mod, M), gof = gr_prob(aki_mod, M))
    
    plt <- ggplot(res_df %>% gather(key = model, value = probability, -mag)) +
      geom_line(aes(x = mag, y = probability, colour = model), size = 1) +
      scale_y_log10() +
      ylab('Log Probability') +
      xlab('Magnitude') +
      scale_color_discrete(name = 'Model', labels = c('Aki','Goodness of Fit', 'Max Curvature')) +
      ggtitle('Gutenberg Richter Probability') +
      ggsave('gr_prob_plot.jpeg')
    
    print(plt)
  }
  
 out_df <- rbind(maxc_mod, aki_mod, gof_mod)
 
 out_df = out_df %>%
   dplyr::mutate(weight = rmse / mean(out_df$rmse)) %>%
   dplyr::mutate(norm_weight = weight / sum(weight))
}

#' Create a plot from the b-value analysis
#' 
#' @param ecdf dataframe with ecdf
#' @param cat dataframe with earthquake catalogue
#' @param out_df dataframe with model results
#' @return dataframe with model output
#' @export
b_value_plot <- function(ecdf,cat,out_df){
  # Import a discrete colour vector for lines
  col_vect <- brewer.pal(5,'Set1')
  aki_text <- paste('ML - Mc: ',sprintf("%.1f", round(out_df['aki','mmax_c'],2)),
                    ' b: ',sprintf("%.2f", -round(out_df['aki','b'],2)),
                    ' rmse: ',sprintf("%.2f", round(out_df['aki','rmse'],2),sep=''))
  gof_text <- paste('GOF - Mc: ',sprintf("%.1f", round(out_df['gof','mmax_c'],2)),
                    ' b: ',sprintf("%.2f", -round(out_df['gof','b'],2)),
                    ' rmse: ',sprintf("%.2f", round(out_df['gof','rmse'],2),sep=''))
  maxc_text <- paste('MC - Mc: ',sprintf("%.1f", round(out_df['max_curvature','mmax_c'],2)),
                     ' b: ',sprintf("%.2f", -round(out_df['max_curvature','b'],2)),
                     ' rmse: ',sprintf("%.2f", round(out_df['max_curvature','rmse'],2),sep=''))
  
  xrng = range(ecdf$centres)
  yrng = range(ecdf$cdf)
  
  # B-value plot 
  ggplot(cat, aes(mw)) +
    geom_histogram(aes(mw), binwidth = out_df['gof',"bin_width"]*3) +
    geom_abline(aes(intercept = out_df['aki','a'], slope =  out_df['aki','b'], col = 'Maximum\nLikelihood (ML)\n'), size = 1) +
    geom_text(aes(x, y, label = aki_text), data = data.frame(x = xrng[2], y = yrng[2]), 
              hjust = 1, vjust = 0, size = 3.25, col = col_vect[1]) +
    geom_abline(aes(intercept = out_df['gof','a'], slope = out_df['gof','b'], col = 'Goodness\nof Fit (GOF)\n'), size = 1) +
    geom_text(aes(x, y, label = gof_text), data = data.frame(x = xrng[2], y = yrng[2]), 
              hjust = 1, vjust = 2, size = 3.25, col = col_vect[2]) +
    geom_abline(aes(intercept = out_df['max_curvature','a'], slope = out_df['max_curvature','b'], col = 'Maximum\nCurvature (MC)\n'), size = 1) +
    geom_text(aes(x, y, label = maxc_text), data = data.frame(x = xrng[2], y = yrng[2]), 
              hjust = 1, vjust = 4, size = 3.25, col = col_vect[3]) +
    geom_point(aes(x=centres, y=cdf), data = ecdf) +
    ylab('Cumulative / Probability Density') +
    scale_fill_manual(name='My Lines', values=c("black", "blue")) +
    scale_y_log10(breaks=c(1,10,100,1000,10000)) +
    scale_color_manual('B-Value Models', values = c('Maximum\nLikelihood (ML)\n' = col_vect[1],
                                                    'Goodness\nof Fit (GOF)\n' = col_vect[2], 
                                                    'Maximum\nCurvature (MC)\n' = col_vect[3])) +
    xlab('Moment Magnitude') +
    theme_minimal() +
    ggsave(paste('b_value_plot.jpeg',sep=""), scale = 1.75, width = 90, height = 90, units = "mm")
}
