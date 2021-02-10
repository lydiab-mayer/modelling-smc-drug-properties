#############################
# Script for improving GP accuracy through adaptive sampling
# 
# 18.06.2019
# monica.golumbeanu@unibas.ch
#############################

GP_updated = GP_model

while(ok) {
    ok = FALSE
    Xcand = lhs(10, param_ranges)
    pred_cand = predict(x = Xcand, object = GP_updated$model)
    res_pred = cbind(Xcand, pred_cand$sd2)
    avg_var = mean(pred_cand$sd2)
    plot_gradient(plot_df)
    print("checking")
    if(n_iter<10) {
        ok = TRUE
        print("Resampling ...")
        ordered_sd_pred = res_pred[order(res_pred$sd2, decreasing = TRUE),]
        #new_SIS_data = call_SIS_sim(XX[1:10,1], XX[1:10,2]) call OM script
        train_data = rbind(train_data, new_OM_data)
        GP_updated = update_GP(GP_updated, new_OM_data))
    }
}
