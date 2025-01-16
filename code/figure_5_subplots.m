% use the parameters to simulate the state variables and compare against
% iemg 
clear variables
clear; close all;
% fit iemg to hill equation and get the hill parameters
data_resting=readtable('../raw_data/val_dataset/Initial_state_Avg.xlsx','Sheet','Summary'); % resting levels of state variables
table_emg=readtable("../raw_data/val_dataset/Emg_for_fitting_Avg.xlsx");
x = table_emg{:,1}; y1 = table_emg{:,2}; y2 = table_emg{:,2}-table_emg{:,3}; y3=table_emg{:,2} + table_emg{:,3};
[coeff1,~,~]=poly4(x,y1);
[coeff2,~,~]=poly4(x,y2);
[coeff3,~,~]=poly4(x,y3);
coeff_p=[coeff1 coeff2 coeff3];
[~,o]=size(coeff_p);
%load the estimated parameters
param_table=readtable('params/params.xlsx');
params = param_table.estimate;
% available experimental data
AxisFontSize = 12; LabelFontSize = 14;
% Set temperature and initial SL
TmpC = 37; 
% Set metabolite concentrations
MgATP = 8.2; 
MgADP = data_resting{1,2}*10^-3; 
Pi = data_resting{4,2}; % Experimentally estimated resting levels by Umass team
Pcr = data_resting{2,2};% Experimentally estimated resting levels by Umass team
SL0 = 3.23;
%pH = data_resting{3,2}; % Experimentally estimated resting levels by Umass team
pH=7.2;
H = 1e3*10^-pH; % mM
N0 = 1;
data_Pcr  = readtable('../raw_data/val_dataset/PCr_for_fitting_Avg.xlsx'); % Pcr force
cycle_index_exp=data_Pcr{:,1};
cycles=1:1:max(cycle_index_exp);
cycle_no = mean([7.277777778 7.130434783]); % No of cycles per 10s in US009, US010 and US011
cycle_time=10/cycle_no;
tspan = 0:0.1:cycle_time;
n=length(tspan);
m=length(cycles);
sim_Ftotal = zeros(m,o);
dCK_final = zeros(m,1);
dK3_final = zeros(m,1);
dCK_r_final = zeros(m,1);
sim_Ftotal_cycles=zeros(n,m);
dCK_cycle=zeros(n,m);
dK3_cycle=zeros(n,m);
dCK_r_cycle=zeros(n,m);
pi_p=zeros(m,o);
ADP_p=zeros(m,o);
Pcr_p=zeros(m,o);
H_p=zeros(m,o);
ATP_p=zeros(m,o);
for k=1:o
    init = [zeros(1,9),N0,SL0, Pi,MgADP, Pcr,H,MgATP]; % Initial conditions for the model
    for i=1:m
        SL_set=3.23;
        coeff=coeff_p(:,k);
        %dSL_set=-0.7732; % takes into account the contraction in G.medialis alone in US011
        dSL_set=-1*mean([1.163852393 0.331395203]); % takes into account the contraction in G.medialis alone in US011, US010 and US009
        iemg= ((coeff(1)*(cycles(i)^4))+(coeff(2)*(cycles(i)^3))+(coeff(3)*(cycles(i)^2))+(coeff(4)*(cycles(i)^1))+coeff(5))/100;
        options = odeset('RelTol',1e-3,'AbsTol',1e-6,'MaxStep',5e-3);
        [T, Y] = ode15s(@Model_XB_human_QC,tspan,init,options,TmpC,SL_set,params,iemg,dSL_set,Pcr,H);
        init(10)=Y(n,10);%N
        init(12)=Y(n,12);%Pi
        init(13)=Y(n,13);%ADP
        init(14)=Y(n,14);%Pcr
        init(15)=Y(n,15);%H
        init(16)=Y(n,16);%ATP
        pi_p(i,k)=Y(n,12);
        ADP_p(i,k)=Y(n,13);
        Pcr_p(i,k)=Y(n,14);
        H_p(i,k)=Y(n,15);
        ATP_p(i,k)=Y(n,16);
        for j=1:n
            [~, sim_Ftotal_cycles(j,i),~,~,~,~,~,~,dCK_cycle(j,i),dK3_cycle(j,i),dCK_r_cycle(j,i),~,~,~,~] = Model_XB_human_QC(T(j),Y(j,:),TmpC,SL_set,params,iemg,dSL_set,Pcr,H);
        end
        sim_Ftotal(i,k) = sim_Ftotal_cycles(n,i);
        dCK_final(i) = dCK_cycle(n,i);
        dK3_final(i) = dK3_cycle(n,i);
        dCK_r_final(i) = dCK_r_cycle(n,i);
    end
end
% convert simulated force into power
dispt=readtable("../raw_data/val_dataset/dsdt_for_fitting_Avg.xlsx");
dispt=dispt{1:453,:};
sim_power = dispt(:,2).*sim_Ftotal(:,1);
% proton concentration unit conversion or pH vs [H] interconversion
%H_fig = H_p*10^6; % unit conversion to nano molar
H_fig = -log10(H_p*10^-3); %pH vs [H] interconversion
% Compile simulations
simulations1=[pi_p(:,1) Pcr_p(:,1) ADP_p(:,1) H_fig(:,1)];
simulations2=[pi_p(:,2) Pcr_p(:,2) ADP_p(:,2) H_fig(:,2)];
simulations3=[pi_p(:,3) Pcr_p(:,3) ADP_p(:,3) H_fig(:,3)];
%y_labels={'Pi (mM)','PCr (mM)','ADP (mM)','[H] (nM)','ATP (mM)'};
y_labels={'Pi (mM)','PCr (mM)','ADP (mM)','pH','Power (W)'};
legend={'Mean iEMG','Min_iEMG'};
exp_flag = [1 1 1 2];
% compile experimental data
power = readtable('../raw_data/val_dataset/power_for_fitting_Avg.xlsx'); 
t6 = power{1:453,1}; power_exp = power{1:453,2}; power_error = power{1:453,3};
phos = readtable('../raw_data/val_dataset/Pi_for_fitting_Avg.xlsx'); 
t1 = phos{:,1}; Pi_exp = phos{:,2};Pi_err= phos{:,3};
PCR = readtable("../raw_data/val_dataset/Pcr_for_fitting_Avg.xlsx");
t3 = PCR{:,1}; PCr_exp = PCR{:,2};PCr_err= PCR{:,3}; 
ADP = readtable("../raw_data/val_dataset/ADP_for_fitting_Avg.xlsx");
t4 = ADP{:,1}; p4 = ADP{:,2};err4= ADP{:,3}; 
ADP_exp=p4*10^-3;ADP_err=err4*10^-3;
PH = readtable("../raw_data/val_dataset/pH_for_fitting_Avg.xlsx");
t5 = PH{:,1}; p5 = PH{:,2};err5= PH{:,3};
%H_exp=(10.^-p5)*1e9; p6d=(10.^-(p5+err5))*1e9; H_err=p6d-H_exp;
H_exp = p5; H_err=err5;
cycle_index=[t1 t3 t4 t5];
exp_data = [Pi_exp PCr_exp ADP_exp H_exp];
exp_err = [Pi_err PCr_err ADP_err H_err];
% Plot the figure for metabolites
[~,q]=size(simulations1);
%filename={'Pi.pdf','PCr.pdf','ADP.pdf','H.pdf','ATP.pdf'};
filename={'Pi_mod1.pdf','PCr_mod1.pdf','ADP_mod1.pdf','pH_mod1.pdf'};
filename2={'Pi_mod.xlsx','PCr_mod.xlsx','ADP_mod.xlsx','H_mod.xlsx'};
rmsd = zeros(length(y_labels),1);
fraction=zeros(length(y_labels),1);
for i=1:q
    figure(i);clf;
    %fill([cycles, fliplr(cycles)],[simulations2(:,i)', fliplr(simulations3(:,i)')],[0.7 0.7 0.7],EdgeColor='none'); hold on
    plot(cycles,simulations1(:,i),'linewidth',2,'Color','k'); hold on;
    %plot(cycles,simulations2(:,i),'--','linewidth',1.5,'Color','k'); hold on;
    %plot(cycles,simulations3(:,i),'--','linewidth',1.5,'Color','k');
    xlim([0 450]);
    ylim([0 inf]);
    xlabel('Cycle Index');
    ylabel(y_labels(i));
    x_txt=0.99;
    y_txt=0.05;
    if exp_flag(i)>0
        hold on;
        errorbar(cycle_index(:,i),exp_data(:,i),exp_err(:,i), ...
            'LineStyle','none','Marker','.', 'MarkerSize',8,'Color',[1 1 1]*0.5,'CapSize',3)
    end
    if exp_flag(i)==0
       ylim([0 12]);
    end
    if exp_flag(i)==2
       ylim([6.5 7.2]);
       y_txt = 0.9;
    end
    if exp_flag(i)>0
        cycle_index_rd=round(cycle_index(:,i));
        tf_temp = ismember(cycles,cycle_index_rd);
        tf_temp(m) = 1;
        simulations_temp=simulations1(tf_temp,i);
        table_temp=array2table([cycle_index(:,i),simulations_temp,exp_data(:,i),exp_err(:,i)],'VariableNames',{'cycle_index','simulation','data','error'});
        diff = simulations_temp-exp_data(:,i);
        rmsd(i)=rms(diff);
        %rmsd_check(i)= (sum((simulations_temp-exp_data(:,i)).^2)/length(exp_data))^0.5;
        w = length(diff);
        flag_temp = zeros(w,1);
        for j1=1:w
            if abs(diff(j1))>exp_err(j1,i)
                flag_temp(j1)=1;
            end
        end
        fraction(i)=(sum(~flag_temp)/w)*100;
        writetable(table_temp,fullfile(pwd,'figure_5_subplots',filename2{i}));
    end
    if exp_flag(i)>0
        rmse_string = sprintf('RMSE: %1.3f (%2.0f%%)',rmsd(i),round(fraction(i)));
        text(x_txt,y_txt,rmse_string,"Units","normalized","HorizontalAlignment","right","VerticalAlignment","baseline","FontSize",8);
    end
    set(gca,'Unit','Inches')
    p = get(gca,'Position');
    set(gca,'Unit','Inches','Position',[p(1) p(2) 1.75 1.25]);
    exportgraphics(figure(i),fullfile(pwd,'figure_5_subplots',filename{i}),'BackgroundColor','w','Resolution',300,'ContentType','vector');
%     close(i)
end
%calculation of RMSD and fraction for power
diff = power_exp-sim_power;
rmsd(5) = rms(diff);
w=length(diff);
flag_temp = zeros(w,1);
for j1=1:w
    if abs(diff(j1))>power_error(j1)
        flag_temp(j1)=1; 
    end
end
fraction(5) = (sum(~flag_temp)/w)*100;
%compile the RMSD
t_rmsd=table(rmsd,'RowNames',y_labels');
frac_table=table(fraction,'RowNames',y_labels);
writetable(t_rmsd,fullfile(pwd,'figure_5_subplots','RMSD_fig1.xlsx'),'WriteRowNames',true);
writetable(frac_table,fullfile(pwd,'figure_5_subplots','frac_fig1.xlsx'),'WriteRowNames',true);
% Plot the figure for force
figure(q+1);clf;
figure(q+1);clf;
errorbar(t6(1:5:453),power_exp(1:5:453),power_error(1:5:453), ...
            'LineStyle','none','Marker','.', 'MarkerSize',8,'Color',[1 1 1]*0.5,'CapSize',3);hold on;
plot(cycles,sim_power,'linewidth',2,'Color','k','DisplayName','Simulation'); 
rmse_string = sprintf('RMSE: %1.3f (%2.0f%%)',rmsd(q+1),round(fraction(q+1)));
text(1,0.9,rmse_string,"Units","normalized","HorizontalAlignment","right","VerticalAlignment","baseline","FontSize",8);
xlabel('Cycle Index');
ylabel('Power (W)');
set(gca,'Unit','Inches')
p = get(gca,'Position');
set(gca,'Unit','Inches','Position',[p(1) p(2) 1.75 1.25]);
exportgraphics(figure(q+1),fullfile(pwd,'figure_5_subplots','Force_mod.pdf'),'BackgroundColor','w','Resolution',300,'ContentType','vector');