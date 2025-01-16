filenames={'../raw_data/RMSE_val/Pi_for_rmse_calc.xlsx','../raw_data/RMSE_val/PCr_for_rmse_calc.xlsx','../raw_data/RMSE_val/ADP_for_rmse_calc.xlsx','../raw_data/RMSE_val/pH_for_rmse_calc.xlsx'};
m=length(filenames);
rmse=zeros(m+1,1);
r_labels={'Pi (mM)','PCr (mM)','ADP (mM)','pH','Force (N)'};
for i=1:m
    data_temp = readtable(filenames{i});
    rmse(i)=rms(data_temp{:,2}-data_temp{:,3});
end
rmse(m+1) = rms(62.88525156-39.30328223);
t_rmse=table(rmse,'RowNames',r_labels');
writetable(t_rmse,fullfile('Table_3_column3','RMSD_table3_column3.xlsx'),'WriteRowNames',true);