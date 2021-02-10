#########################################
# plot GP sensitivity analysis results
#########################################
library(ggplot2)
library(dplyr)

create_color_palette_sens = function() {
    p = c("Halflife", "Efficacy", "Coverage", "Access", "Drug_halflife", "Drug_efficacy")
    categories = factor(p, levels = p)
    library(RColorBrewer)
    # myColors = brewer.pal(length(categories),"Set1")
    myColors = c("#6CC3B9", "#3792BF", "#0D539C", "#fec44f", "#efedf5", "#bcbddc")
    names(myColors) = levels(categories)
    return(myColors)
    # colScale = scale_colour_manual(name = "grp",values = myColors)
}

plot_sens = function(plot_df, plot_val, plot_type, plot_title, plot_file) {
    if(plot_type == "barplot") {
        plotType = geom_bar(stat = "identity", alpha = 0.9) 
    } else {
        plotType = geom_area(alpha = 0.9)
    }
    
    if(plot_val == "main_effects") {
        y_val = "Main effects"
    } else {
        y_val = "Total effects"
    }
    intervention_colors = create_color_palette_sens()
    ggplot(plot_df, aes_string(x = factor(plot_df$EIR), y = plot_val, group = "parameter", fill = "parameter")) +
        theme_bw(base_size=14) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        plotType + theme(strip.background = element_rect(colour="white", fill="white")) +
        facet_wrap(~Seasonality + Biting_pattern) + scale_fill_manual(name = "parameter", values = intervention_colors) +
        labs( x = "EIR", y = y_val, title = plot_title) + theme(legend.position = "none") + 
        theme(legend.position="top") + guides(fill=guide_legend(title="")) 
    ggsave(plot_file, width = 10, height = 5)
}

plot_sens_GP = function(train_dir, plot_dir, param_ranges_file, plot_title, exp_name) {
    # train_dir = "~/MMC/elimination_modeling/results/GP/sensitivity_analysis/prev_red_E14_RCD_before/"
    load(param_ranges_file)
    file.names = dir(train_dir, pattern ="_sidx.RData", full.names = TRUE)
    final_plot_df = NULL
    points_df = NULL
    for(i in 1:length(file.names)){
        load(file.names[i])
        params = factor(rownames(param_ranges), levels = rownames(param_ranges))
        sens_df = cbind.data.frame(as.double(sobol_idx_list$EIR), sobol_idx_list$seasonality, sobol_idx_list$biting_pattern, params, sobol_idx_list$S_eff, sobol_idx_list$T_eff)
        final_plot_df = rbind.data.frame(final_plot_df, sens_df)
    }
    final_plot_df = final_plot_df[-which(final_plot_df$params == "EIR"),]
    colnames(final_plot_df) = c("EIR", "Seasonality", "Biting_pattern", "parameter", "main_effects", "total_effects")
    final_plot_df$parameter = factor(final_plot_df$parameter, levels = c("Coverage", "Efficacy", "Halflife", "Access"))
    # Normalize effects for area plot
    final_plot_df_avg = (final_plot_df %>% group_by(EIR, Seasonality, Biting_pattern) %>% 
                             mutate(main_effects = main_effects/sum(main_effects), 
                                    total_effects = total_effects/sum(total_effects))) 
    
    # Plotting Main effects 
    plot_sens(final_plot_df, "main_effects", "barplot", plot_title,
              paste(plot_dir, "main_effects_barplot_", exp_name, ".pdf", sep="")) 
    plot_sens(final_plot_df_avg, "main_effects", "area", plot_title,
              paste(plot_dir, "main_effects_area_", exp_name, ".pdf", sep="")) 
    # Plotting Total effects
    plot_sens(final_plot_df, "total_effects", "barplot", plot_title,
              paste(plot_dir, "total_effects_barplot_", exp_name, ".pdf", sep="")) 
    plot_sens(final_plot_df_avg, "total_effects", "area", plot_title,
              paste(plot_dir, "total_effects_area_", exp_name, ".pdf", sep="")) 
}

# # Only for testing:
# plot_dir = "~/MMC/TPP/figures/sensitivity_analysis/MAB_once_3_years/"
# 
# plot_sens_GP("~/MMC/TPP/simulations/MAB_once_3years_avg_prev/gp_6/sensitivity/", 
#                 plot_dir, "~/MMC/TPP/simulations/MAB_once_3years_avg_prev/param_ranges.RData",
#                 "", "MAB_once_3years_avg_prev_6")

