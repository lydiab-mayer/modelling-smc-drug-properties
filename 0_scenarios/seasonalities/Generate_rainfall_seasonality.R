###############################################################
###
### Adapting rainfall data to generate seasonal profiles
### Used for MOCK example
### NN
### 15.06.2021
###
###############################################################

# User 
user = strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]

# Code from:

############################################################
# SEASONALITY
#
# Determine national-level seasonality profile for each country
# based on average rainfall data over a number of years.
#
# IPM project by G.Yang, A.J.Shattock, and M.A.Penny
# Code by A.J.Shattock - andrewjames.shattock@unibas.ch
############################################################

rm(list = ls())
library(randomcoloR)
library(ggplot2)
library(dplyr)
library(tidyr)
library(directlabels)
library(ggrepel)

# ---------------------------------------------------------
# Compute seasonality patterns for each country.
# ---------------------------------------------------------

# Load rainfall data from file
rainfall_data = read.csv(paste0("/scicore/home/penny/",user,"/M3TPP/analysisworkflow/0_scenarios/seasonalities/rainfall.csv"))

# Group and average rainfall by country and month
rainfall_mean = rainfall_data %>%
  group_by(country, month) %>%
  summarise_at(vars(rainfall), list(~mean(., na.rm = TRUE)))

# Convert to dataframe
rainfall_df = as.data.frame(rainfall_mean)

# shift one month delay, i.e. Jan value shifted to Feb
rainfall_df$month = rainfall_df$month + 1
rainfall_df[rainfall_df == 13] = 1

# Initiate seasonality table
season_tab = NULL

# Country list
country_list = as.character(unique(rainfall_data$country))

# Iterate through countries
for (country in country_list) {
  
  # Extract rainful data for country and sort by month
  country_df = rainfall_df[rainfall_df$country == country, ]
  country_df = country_df[order(country_df$month), ]
  
  # Normalise to get seasonality
  seasonality = country_df$rainfall / max(country_df$rainfall)
  
  # Append this to other seasonality vectors
  season_tab = rbind(season_tab, seasonality)
}

# Rename columns and remove row names
# colnames(season_tab) = paste0("seasonality_", 1 : 12)
colnames(season_tab) = c("Jan","Feb","Mar","Apr","May","Jun",
                      "Jul","Aug","Sep","Oct","Nov","Dec")
rownames(season_tab) = NULL

# Bind a column for country name and convert to dataframe
country   = as.character(country_list)
season_df = as.data.frame(cbind(country, season_tab), 
                          stringsAsFactors = FALSE)

# Transform data
seasonality_plot_data_clean = season_df %>% 
  dplyr::group_by(country) %>% 
  tidyr::gather(key = "month", value = "value", -country) %>% 
  dplyr::ungroup() %>% 
  mutate(month = factor(month, levels = c("Jan","Feb","Mar","Apr","May","Jun",
                                          "Jul","Aug","Sep","Oct","Nov","Dec")),
         value = as.numeric(value))

# Add country name
library(readr)
all_country_codes = read_delim(paste0("/scicore/home/penny/",user,"/M3TPP/analysisworkflow/0_scenarios/seasonalities/all_country_codes.csv"), 
                               ";", escape_double = FALSE, trim_ws = TRUE)
seasonality_plot_data_clean = merge(seasonality_plot_data_clean,
                                    all_country_codes[,c("alpha_3","name","sub_region")], 
                                    by.x = "country", by.y = "alpha_3",
                                    all.x = T)

# Plot
ggplot(seasonality_plot_data_clean) +
  geom_line(aes(x=month, y=value, group = name, color = name))

# Re-order data to set max transmission in June
seasonality_plot_data_clean_reorder = seasonality_plot_data_clean %>% 
  dplyr::group_by(country, name) %>% 
  dplyr::arrange(country, desc(value)) %>% 
  dplyr::mutate(rank = ifelse(value == max(value), 1, NA)) %>% 
  dplyr::arrange(country, month, rank) %>% 
  dplyr::mutate(n_row = dplyr::row_number(),
         top_row = ifelse(rank == 1, n_row, n_row),
         top_row_rep = max(top_row, na.rm = T),
         diff = n_row - top_row_rep,
         diff_pos = ifelse(diff < 0, diff+12, diff)) %>% 
  dplyr::arrange(country, diff_pos) %>% 
  dplyr::mutate(month_adj = c("Jun","Jul","Aug","Sep","Oct","Nov",
                       "Dec","Jan","Feb","Mar","Apr","May"),
         month_adj = factor(month_adj, levels = c("Jan","Feb","Mar","Apr","May","Jun",
                                                  "Jul","Aug","Sep","Oct","Nov","Dec"))) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(country, name, sub_region, month, month_adj, value)

# Plot
ggplot(seasonality_plot_data_clean_reorder) +
  geom_line(aes(x=month_adj, y=value, group = name, color = name)) +
  labs(x = "month (peak centered to June)", y = "rainfall (normalized)") +
  theme_minimal() +
  # theme_dark() +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "none",
        panel.grid.minor.y = element_blank())

# # Save
# dev.copy(png, filename="/scicore/home/penny/nekkab0000/M3TPP/Experiments/Seasonality_analysis/plots/Rainfall_24.png", width = 750, height = 450)
# dev.off ()

# Group seasonal profiles
seasonality_plot_data_clean_reorder = seasonality_plot_data_clean_reorder %>% 
  dplyr::group_by(country) %>% 
  dplyr::mutate(type0 = NA,
         type0 = ifelse(month_adj == "Oct" & value > 0.25, "perennial",NA),
         type0 = ifelse(length(unique(type0)) > 1, "perennial","seasonal"),
         type1 = NA,
         type1 = ifelse(month_adj == "Feb" & value < 0.25, "seasonal",type1),
         type1 = ifelse(length(unique(type1)) > 1, "seasonal",type0),) %>% 
  dplyr::ungroup()

# Plot
ggplot(seasonality_plot_data_clean_reorder) +
  geom_line(aes(x=month_adj, y=value, group = name, color = name)) +
  facet_wrap(.~type1)

# Group by peaks & early or late tail
seasonality_plot_data_clean_reorder = seasonality_plot_data_clean_reorder %>% 
  dplyr::group_by(country) %>% 
  dplyr::mutate(type2 = NA,
         type2 = ifelse(month_adj == "Sep" & value > 0.25 & type1=="seasonal", "wide seasonal", NA),
         type2 = ifelse(month_adj == "Mar" & value > 0.25 & type1=="seasonal", "wide seasonal", type2),
         type2 = ifelse(length(unique(type2)) > 1, "wide seasonal",type2),
         type2 = ifelse(type1=="seasonal" & is.na(type2), "sharp seasonal", type2),
         type2 = ifelse(is.na(type2),type1,type2),
         type3 = ifelse(type2=="wide seasonal" & month_adj == "Sep" & value > 0.5, "wide tail seasonal",NA),
         type3 = ifelse(length(unique(type3)) > 1, "wide tail seasonal",type3), 
         type4 = ifelse(type2=="wide seasonal" & month_adj == "Mar" & value > 0.5, "wide head seasonal",NA),
         type4 = ifelse(length(unique(type4)) > 1, "wide head seasonal",type4),
         type5 = type2,
         type5 = ifelse(is.na(type3),type5,type3),
         type5 = ifelse(is.na(type4),type5,type4),
         type6 = NA,
         type6 = ifelse(month_adj == "Aug" & value > 0.5 & type1=="perennial","wide tail perennial",NA),
         type6 = ifelse(length(unique(type6)) > 1, "wide tail perennial",type6),
         type7 = ifelse(month_adj == "Jan" & value > 0.75 & type1=="perennial","2-peak perennial",NA),
         type7 = ifelse(length(unique(type7)) > 1, "2-peak perennial",type7),
         type8 = type5,
         type8 = ifelse(is.na(type6),type5,type6),
         type8 = ifelse(is.na(type7),type5,type7),
         type8 = ifelse(type8 == "2-peak perennial","two_peak_perennial",type8),
         type8 = gsub(" ","_",type8),
         type8 = factor(type8, levels = c("perennial", "two_peak_perennial","sharp_seasonal",
                                          "wide_seasonal","wide_head_seasonal","wide_tail_seasonal"))) %>% 
  dplyr::ungroup()

# Colors
n <- length(unique(seasonality_plot_data_clean_reorder$country))
# palette <- distinctColorPalette(n)
palette = c("#DD4BDE","#DAAE99","#DF4A9D","#AAECB9","#B4E8E4","#DE7985","#64E8D4",
            "#AC9BDD","#6E74DB","#7548E7","#CC7AD4","#63E192","#E06E43","#D3CEE1",
            "#83E64F","#E0B158","#E3E2C4","#92B17D","#E1E041","#D4E483","#74A3D1",
            "#7D746B","#5DC3D2","#DEA6CB")

# Plot
ggplot(seasonality_plot_data_clean_reorder) +
  geom_line(aes(x=month_adj, y=value, group = name, color = name)) +
  facet_wrap(.~type8) +
  scale_color_manual(values=palette) +
  geom_dl(aes(x=month_adj, y=value, label = name, color=name), method = list(dl.combine("last.bumpup"), cex = 0.7)) +
  # coord_cartesian(expand = c(1,2)) +
  labs(x = "month (peak centered to June)", y = "rainfall (normalized)") +
  expand_limits(x= c(0, 15)) +
  theme_minimal() +
  # theme_dark() +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "none",
        panel.grid.minor.y = element_blank())

# # Save
# dev.copy(png, filename="/scicore/home/penny/nekkab0000/M3TPP/Experiments/Seasonality_analysis/plots/Rainfall_24_grouped_6_cat.png", width = 750, height = 450)
# dev.off ()

# Transform data for Fourier conversion
Seasonality = seasonality_plot_data_clean_reorder[,c("country","month_adj","value")]
colnames(Seasonality) = c("Seasonality","Month","Value")
Seasonality = Seasonality %>% 
  group_by(Seasonality) %>% 
  tidyr::spread(Month, Value) %>% 
  ungroup() %>% 
  mutate(Seasonality = factor(Seasonality))

# # Save
# write.table(Seasonality, file = "/scicore/home/penny/nekkab0000/M3TPP/Experiments/Seasonality_analysis/Seasonality_24_rainfall_centered_June.txt",
#             sep="\t", quote = F, row.names = F)

# # Transform
# Seasonality_F = seasonality_plot_data_clean_reorder[,c("country","month_adj","value")] %>% 
#   group_by(country) %>% 
#   arrange(country,month_adj)
# 
# # Transform
# Seasonality_F_data = unlist(c(Seasonality_F[1,1],Seasonality_F[1:12,3],
#                               Seasonality_F[(1+12*1),1],Seasonality_F[(1+12*1):(12+12*1),3],
#                               Seasonality_F[(1+12*2),1],Seasonality_F[(1+12*2):(12+12*2),3],
#                               Seasonality_F[(1+12*3),1],Seasonality_F[(1+12*3):(12+12*3),3],
#                               Seasonality_F[(1+12*4),1],Seasonality_F[(1+12*4):(12+12*4),3],
#                               Seasonality_F[(1+12*5),1],Seasonality_F[(1+12*5):(12+12*5),3],
#                               Seasonality_F[(1+12*6),1],Seasonality_F[(1+12*6):(12+12*6),3],
#                               Seasonality_F[(1+12*7),1],Seasonality_F[(1+12*7):(12+12*7),3],
#                               Seasonality_F[(1+12*8),1],Seasonality_F[(1+12*8):(12+12*8),3],
#                               Seasonality_F[(1+12*9),1],Seasonality_F[(1+12*9):(12+12*9),3],
#                               Seasonality_F[(1+12*10),1],Seasonality_F[(1+12*10):(12+12*10),3],
#                               Seasonality_F[(1+12*11),1],Seasonality_F[(1+12*11):(12+12*11),3],
#                               Seasonality_F[(1+12*12),1],Seasonality_F[(1+12*12):(12+12*12),3],
#                               Seasonality_F[(1+12*13),1],Seasonality_F[(1+12*13):(12+12*13),3],
#                               Seasonality_F[(1+12*14),1],Seasonality_F[(1+12*14):(12+12*14),3],
#                               Seasonality_F[(1+12*15),1],Seasonality_F[(1+12*15):(12+12*15),3],
#                               Seasonality_F[(1+12*16),1],Seasonality_F[(1+12*16):(12+12*16),3],
#                               Seasonality_F[(1+12*17),1],Seasonality_F[(1+12*17):(12+12*17),3],
#                               Seasonality_F[(1+12*18),1],Seasonality_F[(1+12*18):(12+12*18),3],
#                               Seasonality_F[(1+12*19),1],Seasonality_F[(1+12*19):(12+12*19),3],
#                               Seasonality_F[(1+12*20),1],Seasonality_F[(1+12*20):(12+12*20),3],
#                               Seasonality_F[(1+12*21),1],Seasonality_F[(1+12*21):(12+12*21),3],
#                               Seasonality_F[(1+12*22),1],Seasonality_F[(1+12*22):(12+12*22),3],
#                               Seasonality_F[(1+12*23),1],Seasonality_F[(1+12*23):(12+12*23),3],
#                               Seasonality_F[(1+12*24),1],Seasonality_F[(1+12*24):(12+12*24),3]))
# 
# # Save
# write.table(Seasonality_F_data, file = "/scicore/home/penny/nekkab0000/M3TPP/Experiments/Seasonality_analysis/Seasonality_24_list.txt",
#             sep="\t", quote = F, row.names = F)



#####################################
# Average smoothing by seasonal group
#####################################

# Group and take mean of monthly rainfall
seasonal_group_avg = seasonality_plot_data_clean_reorder %>% 
  dplyr::select(type8, month_adj, value) %>% 
  dplyr::group_by(type8, month_adj) %>% 
  dplyr::summarise(value_avg = mean(value, na.rm = T)) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(type8, month_adj, value_avg)

# Colors
n <- length(unique(seasonal_group_avg$type8))
palette <- distinctColorPalette(n)

# Plot
ggplot(seasonal_group_avg) +
  geom_point(aes(x=month_adj, y=value_avg, group = type8), color = "gray50") +
  geom_line(aes(x=month_adj, y=value_avg, group = type8, color = type8)) +
  facet_wrap(.~type8) +
  scale_color_manual(values=palette) +
  # geom_dl(aes(x=month_adj, y=value_avg, label = type8, color=type8), method = list(dl.combine("last.bumpup"), cex = 0.7)) +
  labs(x = "month (peak centered to June)", y = "average rainfall (normalized)") +
  theme_minimal() +
  # theme_dark() +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "none",
        panel.grid.minor.y = element_blank())

# # Save
# dev.copy(png, filename="/scicore/home/penny/nekkab0000/M3TPP/Experiments/Seasonality_analysis/plots/Rainfall_24_grouped_6_cat_averaged.png", width = 750, height = 450)
# dev.off ()

# Transform
Seasonality_avg_F = seasonal_group_avg[,c("type8","month_adj","value_avg")] %>% 
  dplyr::arrange(type8,month_adj) %>% 
  dplyr::ungroup()

# Transform
n_groups = length(unique(Seasonality_avg_F$type8))
Seasonality_avg_F_data = unlist(c(Seasonality_avg_F[1,1],Seasonality_avg_F[1:12,3],
                              Seasonality_avg_F[(1+12*1),1],Seasonality_avg_F[(1+12*1):(12+12*1),3],
                              Seasonality_avg_F[(1+12*2),1],Seasonality_avg_F[(1+12*2):(12+12*2),3],
                              Seasonality_avg_F[(1+12*3),1],Seasonality_avg_F[(1+12*3):(12+12*3),3],
                              Seasonality_avg_F[(1+12*4),1],Seasonality_avg_F[(1+12*4):(12+12*4),3],
                              Seasonality_avg_F[(1+12*5),1],Seasonality_avg_F[(1+12*5):(12+12*5),3]))

# # Save
# write.table(Seasonality_avg_F_data, file = "/scicore/home/penny/nekkab0000/M3TPP/Experiments/Seasonality_analysis/Seasonality_6_groups_avg.txt",
#             sep="\t", quote = F, row.names = F)

########################################## SHIFTED (manually) 4 coefficient Fourier!

# Using data from Seasonality_6_groups.txt and shifting peaks to get new center with Fourier transformation
# perennial: stays the same
# two_peak_perennial: stays the same
# wide_head_seasonal: stays the same
# wide_tail_seasonal: shift peak forward 1 month
# wide_seasonal: shift peak back 1 month

# Fourier transformation with 3 components
Fourier_4_coef = NULL

# List of groups
unique(Seasonality_avg_F$type8)

# perennial
a1="-0.3372632662455241"; b1="-0.12224390109380086"
a2="-0.06829960644245148"; b2="-0.24529226620992026"
a3="-0.028720036149024963"; b3="-0.03452832500139872"
Fourier_4_coef = data.frame(Seasonality = "perennial",
                            a1=a1, b1=b1,
                            a2=a2, b2=b2,
                            a3=a3, b3=b3)

# two_peak_perennial
a1="0.046462630232175194"; b1="0.27645963430404663"
a2="0.39200230439503986"; b2="-0.32943689823150635"
a3="0.0372484028339386"; b3="0.008385255932807922"
Fourier_4_coef = rbind(Fourier_4_coef,
                       data.frame(Seasonality = "twopeakperennial",
                                  a1=a1, b1=b1,
                                  a2=a2, b2=b2,
                                  a3=a3, b3=b3))

# wide_head_seasonal
a1="-0.5257936318715414"; b1="1.1914607683817546"
a2="0.38757216930389404"; b2="0.1180146833260854"
a3="-0.03967851400375366"; b3="-0.05106572310129801"
Fourier_4_coef = rbind(Fourier_4_coef,
                       data.frame(Seasonality = "wideheadseasonal",
                                  a1=a1, b1=b1,
                                  a2=a2, b2=b2,
                                  a3=a3, b3=b3))

# wide_tail_seasonal (modified manually)
# 0.0506786881530791
# 0.0760570843540103
# 0.0897228014679389
# 0.185721178557586
# 0.513330668461677
# 1
# 0.95
# 0.6
# 0.6
# 0.6
# 0.6
# 0.2
a1="-1.2338287830352783"; b1="-0.5501441558202108"
a2="-0.16329852739969888"; b2="-0.37494564056396484"
a3="-0.14171481132507324"; b3="0.07149893045425415"
Fourier_4_coef = rbind(Fourier_4_coef,
                       data.frame(Seasonality = "widetailseasonal",
                                  a1=a1, b1=b1,
                                  a2=a2, b2=b2,
                                  a3=a3, b3=b3))

# wide_seasonal 
a1="-1.3715341885884602"; b1="1.196649392445882"
a2="0.15281548102696738"; b2="0.16501147548357645"
a3="-0.0321511427561442"; b3="-0.07909560203552246"
Fourier_4_coef = rbind(Fourier_4_coef,
                       data.frame(Seasonality = "wideseasonal",
                                  a1=a1, b1=b1,
                                  a2=a2, b2=b2,
                                  a3=a3, b3=b3))

# sharp_seasonal 
a1="-2.3171836535135903"; b1="2.114041487375895"
a2="0.384947935740153"; b2="0.3511592149734497"
a3="0.18097086747487387"; b3="-0.1918795108795166"
Fourier_4_coef = rbind(Fourier_4_coef,
                       data.frame(Seasonality = "sharpseasonal",
                                  a1=a1, b1=b1,
                                  a2=a2, b2=b2,
                                  a3=a3, b3=b3))

# Save
# write.table(Fourier_4_coef, file = "/scicore/home/penny/nekkab0000/M3TPP/Experiments/Seasonality_analysis/seasonality_centered_June_Fourier_4_coeff_6_groups_shifted5.txt",
#             sep="\t", quote = F, row.names = F)

########################################## MOVING AVERAGE

# Packages
library(smooth)

# Get rolling average
Seasonality_avg_F_rolling = Seasonality_avg_F %>%
  dplyr::group_by(type8) %>% 
  dplyr::arrange(month_adj, type8) %>%
  dplyr::mutate(bi_weekly = zoo::rollmean(value_avg, k = 2, fill = NA)) %>% 
  dplyr::ungroup()

# Merge
list_per_month = NULL
list_per_type = NULL
for(i in Seasonality_avg_F_rolling$type8){
  type = Seasonality_avg_F_rolling[which(Seasonality_avg_F_rolling$type8 == i),]
  type$n_month_adj = as.numeric(type$month_adj)
  for(j in type$n_month_adj){
    list_per_month[[j]] = rbind(cbind(i,j,as.numeric(type[j,"value_avg"])),
                           cbind(i,j+0.5,as.numeric(type[j,"bi_weekly"])))
  }
  list_per_month_all = do.call(rbind, list_per_month)
  list_per_type[[i]] = list_per_month_all
}
Seasonality_bimonthly = do.call(rbind, list_per_type)
Seasonality_bimonthly = data.frame(Seasonality_bimonthly, stringsAsFactors = F)
colnames(Seasonality_bimonthly) = c("type","month_num","value")
Seasonality_bimonthly = Seasonality_bimonthly[which(Seasonality_bimonthly$month_num != "12.5"),]
Seasonality_bimonthly$month_num = as.numeric(Seasonality_bimonthly$month_num)
Seasonality_bimonthly$value = as.numeric(Seasonality_bimonthly$value)

# Plot
ggplot(Seasonality_bimonthly) +
  geom_point(aes(x=month_num, y=value, group = type), color = "gray50") +
  geom_line(aes(x=month_num, y=value, group = type, color = type)) +
  facet_wrap(.~type) +
  # scale_color_manual(values=palette) +
  # geom_dl(aes(x=month_adj, y=value_avg, label = type8, color=type8), method = list(dl.combine("last.bumpup"), cex = 0.7)) +
  labs(x = "month (peak centered to June)", y = "average rainfall (normalized)") +
  theme_minimal() +
  # theme_dark() +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "none",
        panel.grid.minor.y = element_blank())

# Save data and move the values manually to adjust
# write.table(Seasonality_bimonthly, file = "/scicore/home/penny/nekkab0000/M3TPP/Experiments/Seasonality_analysis/seasonality_bimonthly_6_groups.csv",
#             sep=",", quote = F, row.names = F)


########################################## SHIFTED (bimonthly) 4 coefficient Fourier!

# Using data from Seasonality_6_groups.txt and shifting peaks to get new center with Fourier transformation
# move all 2 weeks back & widetailseasonal 1 month back (manually)

# Fourier transformation with 3 components
Fourier_4_coef = NULL

# perennial (move back 1 month)
a1="-0.34348832006039826"; b1="0.0375149068625077"
a2="-0.22421119524085004"; b2="-0.09370369496552841"
a3="-0.06013856763425081"; b3="0.028569091921267303"
Fourier_4_coef = data.frame(Seasonality = "perennial",
                            a1=a1, b1=b1,
                            a2=a2, b2=b2,
                            a3=a3, b3=b3)

# two_peak_perennial (moved back 2 weeks; peak lowered for lower EIR)
a1="0.027392159337582794"; b1="0.292134430097497"
a2="0.25072215951007343"; b2="-0.3785313316013502"
a3="0.028589533722918968"; b3="-0.024979117123976997"
Fourier_4_coef = rbind(Fourier_4_coef,
                       data.frame(Seasonality = "twopeakperennial",
                                  a1=a1, b1=b1,
                                  a2=a2, b2=b2,
                                  a3=a3, b3=b3))

# wide_head_seasonal (moved back 2 weeks; peak moved up for higher EIR; move back 2 weeks)
a1="-0.38356619295866595"; b1="1.2481913359268852"
a2="0.4174445193746816"; b2="0.04990679803101913"
a3="0.02517598349115123"; b3="0.0034630256502524667"
Fourier_4_coef = rbind(Fourier_4_coef,
                       data.frame(Seasonality = "wideheadseasonal",
                                  a1=a1, b1=b1,
                                  a2=a2, b2=b2,
                                  a3=a3, b3=b3))

# wide_tail_seasonal (1.5 month back + modified peaks; higher 1st peak; move forward 1 month)
a1="-1.3636315387228262"; b1="0.25123193989629333"
a2="-0.4411371479863706"; b2="-0.007745087794635607"
a3="0.012321574532467386"; b3="0.16161805650462274"
Fourier_4_coef = rbind(Fourier_4_coef,
                       data.frame(Seasonality = "widetailseasonal",
                                  a1=a1, b1=b1,
                                  a2=a2, b2=b2,
                                  a3=a3, b3=b3))

# wide_seasonal (moved back 1 month; increase peak for higher EIR; 2 weeks back)
a1="-0.7749619691268258"; b1="1.6528371727984885"
a2="0.26055636613265326"; b2="0.008244851361150328"
a3="0.061031947965207306"; b3="-0.0016192380824814672"
Fourier_4_coef = rbind(Fourier_4_coef,
                       data.frame(Seasonality = "wideseasonal",
                                  a1=a1, b1=b1,
                                  a2=a2, b2=b2,
                                  a3=a3, b3=b3))

# sharp_seasonal (moved 2 weeks back; increase peak for higher EIR)
a1="-1.992391005806301"; b1="2.378197379734205"
a2="0.42162575929061225"; b2="0.3168624795001486"
a3="0.2561972659567128"; b3="-0.0415584367254506"
Fourier_4_coef = rbind(Fourier_4_coef,
                       data.frame(Seasonality = "sharpseasonal",
                                  a1=a1, b1=b1,
                                  a2=a2, b2=b2,
                                  a3=a3, b3=b3))

# # Save
# write.table(Fourier_4_coef, file = "/scicore/home/penny/nekkab0000/M3TPP/Experiments/Seasonality_analysis/seasonality_centered_June_Fourier_4_coeff_6_groups_shifted9_bimonthly.txt",
#             sep="\t", quote = F, row.names = F)








# ########################################## ORIGINAL
# 
# # Fourier transformation with 3 components
# Fourier_3_coef = NULL
# 
# # List of groups
# unique(Seasonality_avg_F$type8)
# 
# # perennial
# a1="-0.3372632662455241"; b1="-0.12224390109380086"
# a2="-0.06829960644245148"; b2="-0.24529226620992026"
# Fourier_3_coef = data.frame(Seasonality = "perennial",
#                             a1=a1, b1=b1,
#                             a2=a2, b2=b2)
# 
# # two_peak_perennial
# a1="0.046462630232175194"; b1="0.27645963430404663"
# a2="0.39200230439503986"; b2="-0.32943689823150635"
# Fourier_3_coef = rbind(Fourier_3_coef,
#                        data.frame(Seasonality = "twopeakperennial",
#                             a1=a1, b1=b1,
#                             a2=a2, b2=b2))
# 
# # wide_head_seasonal
# a1="-0.5257936318715414"; b1="1.1914607683817546"
# a2="0.38757216930389404"; b2="0.1180146833260854"
# Fourier_3_coef = rbind(Fourier_3_coef,
#                        data.frame(Seasonality = "wideheadseasonal",
#                                   a1=a1, b1=b1,
#                                   a2=a2, b2=b2))
# 
# # wide_tail_seasonal
# a1="-1.4505152702331543"; a2="-0.030694733063379925"
# a2="-0.441175897916158"; b2="0.011104324211676916"
# Fourier_3_coef = rbind(Fourier_3_coef,
#                        data.frame(Seasonality = "widetailseasonal",
#                                   a1=a1, b1=b1,
#                                   a2=a2, b2=b2))
# 
# # wide_seasonal
# a1="-1.3715341885884602"; b1="1.196649392445882"
# a2="0.15281548102696738"; b2="0.16501147548357645"
# Fourier_3_coef = rbind(Fourier_3_coef,
#                        data.frame(Seasonality = "wideseasonal",
#                                   a1=a1, b1=b1,
#                                   a2=a2, b2=b2))
# 
# # sharp_seasonal
# a1="-2.3171836535135903"; b1="2.114041487375895"
# a2="0.384947935740153"; b2="0.3511592149734497"
# Fourier_3_coef = rbind(Fourier_3_coef,
#                        data.frame(Seasonality = "sharpseasonal",
#                                   a1=a1, b1=b1,
#                                   a2=a2, b2=b2))
# 
# # # Save
# # write.table(Fourier_3_coef, file = "/scicore/home/penny/nekkab0000/M3TPP/Experiments/Seasonality_analysis/seasonality_centered_June_Fourier_3_coeff_6_groups.txt",
# #             sep="\t", quote = F, row.names = F)

#########################################

# Fourier series function (https://stackoverflow.com/questions/41435777/perform-fourier-analysis-to-a-time-series-in-r/41465250):
nff = function(x = NULL, n = NULL, up = 10L, plot = TRUE, add = FALSE, main = NULL, ...){
  #The direct transformation
  #The first frequency is DC, the rest are duplicated
  dff = fft(x)
  #The time
  t = seq(from = 1, to = length(x))
  #Upsampled time
  nt = seq(from = 1, to = length(x)+1-1/up, by = 1/up)
  #New spectrum
  ndff = array(data = 0, dim = c(length(nt), 1L))
  ndff[1] = dff[1] #Always, it's the DC component
  if(n != 0){
    ndff[2:(n+1)] = dff[2:(n+1)] #The positive frequencies always come first
    #The negative ones are trickier
    ndff[length(ndff):(length(ndff) - n + 1)] = dff[length(x):(length(x) - n + 1)]
  }
  #The inverses
  indff = fft(ndff/12, inverse = TRUE)
  idff = fft(dff/12, inverse = TRUE)
  if(plot){
    if(!add){
      plot(x = t, y = x, pch = 16L, xlab = "Time", ylab = "Measurement",
           main = ifelse(is.null(main), paste(n, "harmonics"), main))
      lines(y = Mod(idff), x = t, col = adjustcolor(1L, alpha = 0.5))
    }
    lines(y = Mod(indff), x = nt, ...)
  }
  ret = data.frame(time = nt, y = Mod(indff))
  return(ret)
}

# Loop function through all 6 groups
series_data = NULL
points_data = NULL
for(i in unique(Seasonality_avg_F$type8)){
  # Get data per type
  y = as.vector(Seasonality_avg_F[which(Seasonality_avg_F$type8 == i),3])$value_avg
  
  # Time 
  t = 1:12
  
  # Range
  rg = diff(range(y))
  
  # Run function (3 components)
  series = nff(x = y, n = 3, up = 100L, col = 2L)
  
  # Add group
  series$type8 = i
  
  # Points
  points = data.frame(month = t, y = y, type8 = i)
  
  # Save
  series_data = rbind(series_data,series)
  points_data = rbind(points_data,points)
}

# # Unlist
# series_data = do.call(rbind, series_list)
# points_data = do.call(rbind, points_list)

# Factor groups
series_data$type8 = factor(series_data$type8, levels = c("perennial", "two_peak_perennial","sharp_seasonal",
"wide_seasonal","wide_head_seasonal","wide_tail_seasonal"))
points_data$type8 = factor(points_data$type8, levels = c("perennial", "two_peak_perennial","sharp_seasonal",
                                                           "wide_seasonal","wide_head_seasonal","wide_tail_seasonal"))

# Plot
ggplot() +
  geom_point(data = points_data, aes(x= month, y = y, group = type8), color = "gray50") +
  scale_x_continuous(breaks = seq(1,12,1), labels = c("Jan","Feb","Mar","Apr","May","Jun",
                                                      "Jul","Aug","Sep","Oct","Nov","Dec")) +
  geom_line(data = series_data, aes(x = time, y = y, group = type8, color = type8)) +
  facet_wrap(.~type8) +
  scale_color_manual(values=palette) +
  labs(x = "month (peak centered to June)", y = "average rainfall (normalized)") +
  theme_minimal() +
  # theme_dark() +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "none",
        panel.grid.minor.y = element_blank())

# Save
# dev.copy(png, filename="/scicore/home/penny/nekkab0000/M3TPP/Experiments/Seasonality_analysis/plots/Rainfall_24_grouped_6_cat_averaged_Fourier_3_coeff.png", width = 750, height = 450)
# dev.off ()



# ##########################################
# 
# # Fourier transformation with 6 components
# Fourier_6_coef = NULL
# 
# # AGO
# a1="-2.7895676294962564"; b1="0.14098519086837769"
# a2="-1.2612583637237549"; b2="0.31252602736155194"
# a3="-0.2934271494547526"; b3="0.39860864480336505"
# a4="-0.0471192995707194"; b4="0.1984338959058126"
# a5="-0.005415679886937141"; b5="0.10911770661671956"
# Fourier_6_coef = data.frame(Seasonality = "AGO",
#                             a1=a1, b1=b1,
#                             a2=a2, b2=b2,
#                             a3=a3, b3=b3,
#                             a4=a4, b4=b4,
#                             a5=a5, b5=b5)
# 
# xs <- seq(-2*pi,2*pi,pi/100)
# wave.1 <- sin(3*xs)
# wave.2 <- sin(10*xs)
# par(mfrow = c(1, 2))
# plot(xs,wave.1,type="l",ylim=c(-1,1)); abline(h=0,lty=3)
# plot(xs,wave.2,type="l",ylim=c(-1,1)); abline(h=0,lty=3)
# wave.3 <- 0.5 * wave.1 + 0.25 * wave.2
# plot(xs,wave.3,type="l"); title("Eg complex wave"); abline(h=0,lty=3)
# 
# # BEN
# a1="-0.9387714862823486"; b1="1.836237907409668"
# a2="0.6739360491434733"; b2="0.28585392236709595"
# a3="-0.06444287300109863"; b3="-0.18385382493336996"
# a4="-0.018822749455769856"; b4="0.06869750221570332"
# a5="0.07103974123795827"; b5="-0.06421609719594319"
# Fourier_6_coef = rbind(Fourier_6_coef, 
#                        data.frame(Seasonality = "BEN",
#                             a1=a1, b1=b1,
#                             a2=a2, b2=b2,
#                             a3=a3, b3=b3,
#                             a4=a4, b4=b4,
#                             a5=a5, b5=b5))
# 
# # BFA
# a1="-1.6241013209025066"; b1="2.742680549621582"
# a2="0.8970824877421061"; b2="0.4601362148920695"
# a3="0.06008005142211914"; b3="-0.18529963493347168"
# a4="-0.16285578409830728"; b4="0.049499248464902244"
# a5="0.014558481673399607"; b5="0.06791184842586517"
# Fourier_6_coef = rbind(Fourier_6_coef, 
#                        data.frame(Seasonality = "BFA",
#                                   a1=a1, b1=b1,
#                                   a2=a2, b2=b2,
#                                   a3=a3, b3=b3,
#                                   a4=a4, b4=b4,
#                                   a5=a5, b5=b5))
# 
# # CMR
# a1="-0.13401097059249878"; b1="1.214416742324829"
# a2="0.6186872720718384"; b2="-0.1749900778134664"
# a3="-0.10516266028086345"; b3="-0.06953330834706624"
# a4="0.0026548306147257486"; b4="0.0011580013670027256"
# a5="-0.05570864677429199"; b5="-0.0582513560851415"
# Fourier_6_coef = rbind(Fourier_6_coef, 
#                        data.frame(Seasonality = "CMR",
#                                   a1=a1, b1=b1,
#                                   a2=a2, b2=b2,
#                                   a3=a3, b3=b3,
#                                   a4=a4, b4=b4,
#                                   a5=a5, b5=b5))
# 
# # COD
# a1="-0.35849658648173016"; b1="-0.08349822958310445"
# a2="-0.1635348598162333"; b2="-0.25900622208913165"
# a3="-0.08793869614601135"; b3="0.0022387405236562095"
# a4="-0.0037685136000315347"; b4="0.010190263390541077"
# a5="-0.007104553282260895"; b5="-0.05007663865884145"
# Fourier_6_coef = rbind(Fourier_6_coef, 
#                        data.frame(Seasonality = "COD",
#                                   a1=a1, b1=b1,
#                                   a2=a2, b2=b2,
#                                   a3=a3, b3=b3,
#                                   a4=a4, b4=b4,
#                                   a5=a5, b5=b5))
# 
# # Save
# write.table(Fourier_6_coef, file = "/scicore/home/penny/nekkab0000/M3TPP/Experiments/Seasonality_analysis/Seasonality_centered_June_Fourier_6_coeff_examples.txt",
#             sep="\t", quote = F, row.names = F)