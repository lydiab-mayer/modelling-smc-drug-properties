#####################
# plot GP training results
####################

library(hetGP)
library(dplyr)

test_GP_plot = function(real_values, predicted_values, plot_file) {
    df_plot = cbind.data.frame(real_values, predicted_values)
    colnames(df_plot) = c("real_values", "predicted_values")
    cor_coeff = cor(real_values, predicted_values)
    p = ggplot(df_plot, aes(x=real_values, y=predicted_values)) + geom_point(aes(color = '#E69F00', alpha = .1)) +#geom_hex() + 
        scale_color_manual(values = c('#E69F00')) + 
        theme_bw(base_size=16) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
         labs( x = "Real values", y = "Predicted values", title = paste("r^2 =", formatC(cor_coeff, digits = 3, format = "f")))
    plot(p)
    ggsave(plot_file, p)
}

plot_figures = function(points_df, plot_file_1, plot_file_2, plot_file_3) {
    ggplot(points_df, aes(x=as.factor(points_df$Biting_pattern), y=points_df$Error)) + geom_violin(color = "#d95f0e") +
        theme_bw(base_size=16) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        facet_wrap(~Seasonality) + labs( x = "Biting pattern", y = "Absolute error") +
        theme(strip.background = element_rect(colour="white", fill="white")) + ylim(0, 1)
    ggsave(plot_file_1, width = 10, height = 5)
    
    p = ggplot(points_df, aes(x=True, y=Predicted)) + geom_point(color = "#d95f0e", fill = "#d95f0e", shape = 21, alpha = 0.5, size = 3) + #stat_binhex() + 
        theme_bw(base_size=14) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        facet_wrap(~Seasonality + Biting_pattern, ncol = 3) + labs( x = "True value", y = "Predicted value")+
        theme(strip.background = element_rect(colour="black", fill="white")) + xlim(0, 100) + ylim(0, 100) +
        geom_text(aes(label = paste("r^2 ==", round(points_df$Corr, digits = 2))), x=20, y=90, check_overlap = TRUE, parse = TRUE)
    ggsave(plot_file_2, plot = p, width = 6, height = 5)
    
    p = ggplot(points_df, aes(x=True, y=Predicted)) + stat_binhex() + scale_fill_viridis_c() +
        theme_bw(base_size=14) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        facet_wrap(~Seasonality + Biting_pattern, ncol = 3) + labs( x = "True value", y = "Predicted value")+
        theme(strip.background = element_rect(colour="black", fill="white")) + xlim(0, 100) + ylim(0, 100) +
        geom_text(aes(label = paste("r^2 ==", round(points_df$Corr, digits = 2))), x=20, y=90, check_overlap = TRUE, parse = TRUE)
    ggsave(plot_file_3, plot = p, width = 6, height = 5)
}

plot_performance_GP = function(train_dir, plot_dir, plot_title, exp_name, model_pattern) {
    file.names = dir(train_dir, pattern = ".RData", full.names = TRUE)
    error_df = NULL
    points_df = NULL
    for(i in 1:length(file.names)){
        load(file.names[i])
        error_df = rbind.data.frame(error_df, cbind.data.frame(cv_result$seasonality, cv_result$biting_pattern, cv_result$test_errors))
        points_df = rbind.data.frame(points_df, cbind.data.frame(cv_result$seasonality, 
                                                                 cv_result$biting_pattern, 
                                                                 cv_result$test_points[,1], 
                                                                 cv_result$test_points[,2], 
                                                                 abs(cv_result$test_points[,1] - 
                                                                     cv_result$test_points[,2]),
                                                                 cor(cv_result$test_points[,1], 
                                                                     cv_result$test_points[,2])))
    }
    colnames(error_df) = c("Seasonality", "Biting_pattern", "Test_error")
    colnames(points_df) = c("Seasonality", "Biting_pattern", "True", "Predicted", "Error", "Corr")
    error_df$Seasonality = factor(error_df$Seasonality, levels = c("perennial", "seasonal"))
    error_df$Biting_pattern = factor(error_df$Biting_pattern, levels = c("Low_indoor", "Mid_indoor", "High_indoor"))
    
    # Plot the results
    plot1 = paste(plot_dir, model_pattern, "_cv_error_", exp_name, ".pdf", sep="")
    plot2 = paste(plot_dir, model_pattern, "_cv_corr_", exp_name, ".pdf", sep="")
    plot3 = paste(plot_dir, model_pattern, "_cv_corr_hex_", exp_name, ".pdf", sep="")
    plot_figures(points_df, plot1, plot2, plot3)
}

plot_performance_GP_test = function(gp_dir, test_dir, ranges_file, plot_dir, 
                                    plot_title, exp_name, model_pattern) {
    file.names = dir(gp_dir, pattern = ".RData", full.names = TRUE)
    points_df = NULL
    ranges = load(ranges_file)
    for(i in 1:length(file.names)){
        gp_result_name = load(file.names[i])
        gp_result = get(gp_result_name)
        
        test_set_name = paste0(model_pattern, "_", 
                               gp_result$seasonality, "_", gp_result$biting_pattern, ".txt")
        test_set = read.table(paste0(test_dir, test_set_name), header = TRUE, as.is = TRUE, stringsAsFactors = FALSE)
        # remove the points where the initial prevalence was 0
        if (length(which(test_set$initial_prev == 0)) > 0){
            test_set = test_set[-which(test_set$initial_prev == 0), ]
        }
        test_data = test_set[, c(rownames(param_ranges), "prev_red")]
        
        D = ncol(test_data) - 1
        param_col = c(1:D)
        response_col = D + 1
        # apply the model on the test data and calculate the error
        prediction_train = predict(x = as.matrix(test_data[, param_col]), object = gp_result$GP_model)
        predicted_values = prediction_train$mean
        predicted_values[which(predicted_values < 0)] = 0
        predicted_values[which(predicted_values > 100)] = 100
        
        points_data_frame = cbind.data.frame(test_data[, param_col], 
                                             test_data[, response_col], predicted_values)
        colnames(points_data_frame) = c(rownames(param_ranges), "true_prev_red", "predicted_prev_red")
        agg_points_df = points_data_frame %>% group_by_at(rownames(param_ranges)) %>% summarise(true_prev_red = mean(true_prev_red, na.rm=TRUE),
                                                                       predicted_prev_red = mean(predicted_prev_red, na.rm=TRUE))
        
        MSE_train = sum((agg_points_df$predicted_prev_red - agg_points_df$true_prev_red)^2)/(100*nrow(agg_points_df))
        abs_error = abs(agg_points_df$predicted_prev_red - agg_points_df$true_prev_red)
        # calculate the correlation between true and predicted values
        corr = cor(agg_points_df$predicted_prev_red, agg_points_df$true_prev_red)
        points_df = rbind.data.frame(points_df, cbind.data.frame(gp_result$seasonality, 
                                                                 gp_result$biting_pattern, 
                                                                 agg_points_df$true_prev_red, 
                                                                 agg_points_df$predicted_prev_red, 
                                                                 abs_error, corr))
    }
    colnames(points_df) = c("Seasonality", "Biting_pattern", "True", "Predicted", "Error", "Corr")
    points_df$Seasonality = factor(points_df$Seasonality, levels = c("perennial", "seasonal"))
    points_df$Biting_pattern = factor(points_df$Biting_pattern, levels = c("Low_indoor", "Mid_indoor", "High_indoor"))
    
    # Plot the results
    plot1 = paste(plot_dir, model_pattern, "_error_", exp_name, ".pdf", sep="")
    plot2 = paste(plot_dir, model_pattern, "_corr_", exp_name, ".pdf", sep="")
    plot3 = paste(plot_dir, model_pattern, "_corr_hex_", exp_name, ".pdf", sep="")
    plot_figures(points_df, plot1, plot2, plot3)
}

# # Plot cross validation results
# plot_dir = "~/MMC/TPP/figures/gp_performance/MAB_once_3_years/"
# plot_performance_GP("~/MMC/TPP/simulations/MAB_once_3years_avg_prev/gp_4/trained/",
#                     plot_dir, "MAB once 3 years", "4", "seeds_MAB_once_3_years")
# plot_performance_GP("~/MMC/TPP/simulations/MAB_once_3years_avg_prev/gp_6/trained/",
#                     plot_dir, "MAB once 3 years", "6", "seeds_MAB_once_3_years")

# plot_dir = "~/MMC/TPP/figures/gp_performance/MAB_once_3_years/"
# plot_performance_GP("~/MMC/TPP/simulations/MAB_once_3years_cont_EIR/gp_4/trained/",
#                                         plot_dir, "MAB once 3 years", "4", "seeds_MAB_once_3_years")
# plot_performance_GP("~/MMC/TPP/simulations/MAB_once_3years_cont_EIR/gp_6/trained/",
#                     plot_dir, "MAB once 3 years", "6", "seeds_MAB_once_3_years")
# 
# # Plot the test set results
# plot_performance_GP_test("~/MMC/TPP/simulations/MAB_once_3years_cont_EIR/gp_4/trained/",
#                          "~/MMC/TPP/simulations/test_MAB_once_3years/postprocessing_4/",
#                          "~/MMC/TPP/simulations/MAB_once_3years_cont_EIR/param_ranges.RData",
#                          plot_dir, "MAB once 3 years", "test_4", "seeds_MAB_once_3_years")
# plot_performance_GP_test("~/MMC/TPP/simulations/MAB_once_3years_cont_EIR/gp_6/trained/",
#                          "~/MMC/TPP/simulations/test_MAB_once_3years/postprocessing_6/",
#                          "~/MMC/TPP/simulations/MAB_once_3years_cont_EIR/param_ranges.RData",
#                          plot_dir, "MAB once 3 years", "test_6", "seeds_MAB_once_3_years")
# 
# # Plot the adaptive sampling results on test set
# plot_performance_GP_test("~/MMC/TPP/simulations/MAB_once_3years_cont_EIR/as_4/",
#                          "~/MMC/TPP/simulations/test_MAB_once_3years/postprocessing_4/",
#                          "~/MMC/TPP/simulations/MAB_once_3years_cont_EIR/param_ranges.RData",
#                          plot_dir, "MAB once 3 years", "test_as_4", "seeds_MAB_once_3_years")
# plot_performance_GP_test("~/MMC/TPP/simulations/MAB_once_3years_cont_EIR/as_6/",
#                          "~/MMC/TPP/simulations/test_MAB_once_3years/postprocessing_6/",
#                          "~/MMC/TPP/simulations/MAB_once_3years_cont_EIR/param_ranges.RData",
#                          plot_dir, "MAB once 3 years", "test_as_6", "seeds_MAB_once_3_years")
# 
# # Plot the adaptive sampling results on training set
# plot_performance_GP_test("~/MMC/TPP/simulations/MAB_once_3years_cont_EIR/as_4/",
#                          "~/MMC/TPP/simulations/MAB_once_3years_cont_EIR/postprocessing_4/",
#                          "~/MMC/TPP/simulations/MAB_once_3years_cont_EIR/param_ranges.RData",
#                          plot_dir, "MAB once 3 years", "train_as_4", "seeds_MAB_once_3_years")
# plot_performance_GP_test("~/MMC/TPP/simulations/MAB_once_3years_cont_EIR/as_6/",
#                          "~/MMC/TPP/simulations/MAB_once_3years_cont_EIR/postprocessing_6/",
#                          "~/MMC/TPP/simulations/MAB_once_3years_cont_EIR/param_ranges.RData",
#                          plot_dir, "MAB once 3 years", "train_as_6", "seeds_MAB_once_3_years")
# 
