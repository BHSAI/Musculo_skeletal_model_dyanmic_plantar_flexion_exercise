% power
power = readtable('../raw_data/power_for_fitting_DPF_2.xlsx'); 
t1 = power{1:240,1}; power_exp = power{1:240,2}; power_err = power{1:240,3};
% iemg
iemg = readtable('../raw_data/Emg_for_fitting_DPF.xlsx');
t2 = iemg{1:38,1}; iemg_exp = iemg{1:38,2}; iemg_err = iemg{1:38,3};  
% phosphocreatine
PCR = readtable("../raw_data/Pcr_for_fitting_DPF.xlsx");
t3 = PCR{:,1}; PCr_exp = PCR{:,2};PCr_err= PCR{:,3}; 
% phosphate
phos = readtable('../raw_data/pi_for_fitting_DPF.xlsx'); 
t4 = phos{:,1}; Pi_exp = phos{:,2};Pi_err= phos{:,3};
% ADP
ADP = readtable("../raw_data/ADP_for_fitting_DPF.xlsx");
t5 = ADP{:,1}; p4 = ADP{:,2};err4= ADP{:,3}; 
ADP_exp=p4*10^-3;ADP_err=err4*10^-3;
% pH
PH = readtable("../raw_data/pH_for_fitting_DPF.xlsx");
t6 = PH{:,1}; p5 = PH{:,2};err5= PH{:,3};
H_exp = p5; H_err=err5;
% Compile datasets
cycle_index = [t2 t3 t4 t5 t6];
exp_data = [iemg_exp PCr_exp Pi_exp ADP_exp H_exp];
exp_err = [iemg_err PCr_err Pi_err ADP_err H_err];
y_labels={'Normalized EMG (%)','PCr (mM)','Pi (mM)','ADP (mM)','pH'};
exp_flag = [3 1 1 1 2];
q=length(y_labels);
filename={'iEMG_exp.pdf','PCr_exp.pdf','Pi_exp.pdf','ADP_exp.pdf','pH_exp.pdf'};
% plot all except power
for i=1:q
    figure(i);clf;
    errorbar(cycle_index(:,i),exp_data(:,i),exp_err(:,i), ...
    'LineStyle','none','Marker','.', 'MarkerSize',8,'Color',[0 0 0],'CapSize',3)
    xlim([0 250]);
    ylim([0 inf]);
    xlabel('Cycle Index');
    ylabel(y_labels(i));
    if exp_flag(i)==3
       ylim([30 70]);
    end
    if exp_flag(i)==2
       ylim([6.5 7.2]);
    end
    set(gca,'Unit','Inches')
    p = get(gca,'Position');
    set(gca,'Unit','Inches','Position',[p(1) p(2) 1.75 1.25]);
    exportgraphics(figure(i),fullfile('figure_3_subplots',filename{i}),'BackgroundColor','w','Resolution',300,'ContentType','vector');   
end
%plot power
figure(q+1);clf;
errorbar(t1(1:5:240),power_exp(1:5:240),power_err(1:5:240), ...
    'LineStyle','none','Marker','.', 'MarkerSize',8,'Color',[0 0 0],'CapSize',3)
xlim([0 250]);
ylim([0 inf]);
xlabel('Cycle Index');
ylabel('Power (W)');
set(gca,'Unit','Inches')
p = get(gca,'Position');
set(gca,'Unit','Inches','Position',[p(1) p(2) 1.75 1.25]);
exportgraphics(figure(q+1),fullfile('figure_3_subplots','Power_exp.pdf'),'BackgroundColor','w','Resolution',300,'ContentType','vector');   