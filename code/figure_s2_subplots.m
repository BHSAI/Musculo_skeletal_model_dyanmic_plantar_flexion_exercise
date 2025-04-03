% use the parameters to simulate the state variables and compare against
% iemg 
clear variables
clear; close all;
% fit iemg to hill equation and get the hill parameters
data_resting=readtable('../raw_data/Initial_state.xlsx','Sheet','Summary'); % resting levels of state variables
table_emg=readtable("../raw_data/Emg_for_fitting_DPF.xlsx");
x = table_emg{:,1}; y = table_emg{:,2};
[hill_a,hill_b,hill_c]=classic_hill(x,y);
%load the estimated parameters
load('params/solutions_SI.mat');
params = solutions(1,3).X;
% available experimental data
AxisFontSize = 12; LabelFontSize = 14;
% Set temperature and initial SL
TmpC = 37; 
% Set metabolite concentrations
MgATP = 8.2; 
MgADP = data_resting{1,2}*10^-3; 
%MgADP = 0; 
Pi = data_resting{4,2}; % Experimentally estimated resting levels by Umass team
%Pi = 0; % Experimentally estimated resting levels by Umass team
Pcr = data_resting{2,2};% Experimentally estimated resting levels by Umass team
%Pcr = 40; % Experimentally estimated resting levels by Umass team
SL0 = 3.23;%2.2; % um [114 mm = 1.3758 um; 116 mm = 1.5869; 118mm = 1.8114; 120mm = 2.0403; 128 mm = 2.8346; 130mm = 2.9728; 132mm= 3.0980]
pH = data_resting{3,2}; %Experimentally estimated resting levels by Umass team
H = 1e3*10^-pH; % mM
N0 = 1;
init = [zeros(1,9),N0,SL0, Pi,MgADP, Pcr,H,MgATP]; % Initial conditions for the model
data_Pcr  = readtable('../raw_data/Pcr_for_fitting_DPF.xlsx'); % Pcr force
cycle_index_exp=data_Pcr{:,1};
cycles=1:1:max(cycle_index_exp);
cycle_time=10/6.33;% average of 6.33 cycles for 10 secs calculated from the data
tspan = 0:0.1:cycle_time;
n=length(tspan);
m=length(cycles);
sim_Ftotal = zeros(m,1);
dCK_final = zeros(m,1);
dK3_final = zeros(m,1);
dCK_r_final = zeros(m,1);
sim_Ftotal_cycles=zeros(n,m);
dCK_cycle=zeros(n,m);
dK3_cycle=zeros(n,m);
dCK_r_cycle=zeros(n,m);
pi_p=zeros(m,1);
ADP_p=zeros(m,1);
Pcr_p=zeros(m,1);
H_p=zeros(m,1);
ATP_p=zeros(m,1);
for i=1:m
    SL_set=3.23;
    dSL_set=-0.68; %takes into account the contraction in G.medialis alone
    iemg=((hill_a*(cycles(i)^hill_b))/(cycles(i)+hill_c))/100;
    options = odeset('RelTol',1e-3,'AbsTol',1e-6,'MaxStep',5e-3);
    [T, Y] = ode15s(@Model_XB_human_QC_SI,tspan,init,options,TmpC,SL_set,params,iemg,dSL_set,Pcr,H);
    init(10)=Y(n,10);%N
    init(12)=Y(n,12);%Pi
    init(13)=Y(n,13);%ADP
    init(14)=Y(n,14);%Pcr
    init(15)=Y(n,15);%H
    init(16)=Y(n,16);%ATP
    pi_p(i)=Y(n,12);
    ADP_p(i)=Y(n,13);
    Pcr_p(i)=Y(n,14);
    H_p(i)=Y(n,15);
    ATP_p(i)=Y(n,16);
    for j=1:n
        [~, sim_Ftotal_cycles(j,i),~,~,~,~,~,~,dCK_cycle(j,i),dK3_cycle(j,i),dCK_r_cycle(j,i),~,~,~,~] = Model_XB_human_QC_SI(T(j),Y(j,:),TmpC,SL_set,params,iemg,dSL_set,Pcr,H);
    end
    sim_Ftotal(i) = sim_Ftotal_cycles(n,i);
    dCK_final(i) = dCK_cycle(n,i);
    dK3_final(i) = dK3_cycle(n,i);
    dCK_r_final(i) = dCK_r_cycle(n,i);
end
% convert simulated force into power
dispt=readtable("../raw_data/dsdt_for_fitting_DPF_2.xlsx");
dispt=dispt{1:240,:};
sim_power = dispt(:,2).*sim_Ftotal;
% Compile simulations
%H_fig = H_p*10^6; % [H] unit conversion from mM to nM
H_fig = -log10(H_p*10^-3); %pH vs [H] interconversion
simulations=[pi_p Pcr_p ADP_p H_fig];
%y_labels={'Pi (mM)','PCr (mM)','ADP (mM)','[H] (nM)','ATP (mM)'};
y_labels={'Pi (mM)','PCr (mM)','ADP (mM)','pH','Power (W)'};
% compile experimental data
exp_flag = [1 1 1 2];
power = readtable('../raw_data/power_for_fitting_DPF_2.xlsx'); 
t6 = power{1:240,1}; power_exp = power{1:240,2}; power_error = power{1:240,3};
phos = readtable('../raw_data/pi_for_fitting_DPF.xlsx'); 
t1 = phos{:,1}; Pi_exp = phos{:,2};Pi_err= phos{:,3};
PCR = readtable("../raw_data/Pcr_for_fitting_DPF.xlsx");
t3 = PCR{:,1}; PCr_exp = PCR{:,2};PCr_err= PCR{:,3}; 
ADP = readtable("../raw_data/ADP_for_fitting_DPF.xlsx");
t4 = ADP{:,1}; p4 = ADP{:,2};err4= ADP{:,3}; 
ADP_exp=p4*10^-3;ADP_err=err4*10^-3;
PH = readtable("../raw_data/pH_for_fitting_DPF.xlsx");
t5 = PH{:,1}; p5 = PH{:,2};err5= PH{:,3};
%H_exp=(10.^-p5)*1e9; p6d=(10.^-(p5+err5))*1e9; H_err=p6d-H_exp;
H_exp = p5; H_err=err5;
cycle_index=[t1 t3 t4 t5];
exp_data = [Pi_exp PCr_exp ADP_exp H_exp];
exp_err = [Pi_err PCr_err ADP_err H_err];
% Plot the figure for metabolites
[~,q]=size(simulations);
filename={'Pi_mod.pdf','PCr_mod.pdf','ADP_mod.pdf','H_mod.pdf'};
filename2={'Pi_mod.xlsx','PCr_mod.xlsx','ADP_mod.xlsx','H_mod.xlsx'};
rmsd = zeros(length(y_labels),1);
fraction=zeros(length(y_labels),1);
for i=1:q
    figure(i);clf;
    plot(cycles,simulations(:,i),'linewidth',2,'Color','k');
    xlim([0 250]);
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
        simulations_temp=simulations(tf_temp,i);
        table_temp=[cycle_index(:,i),simulations_temp,exp_data(:,i)];
        diff = simulations_temp-exp_data(:,i);
        rmsd(i)=rms(diff);
        %rmsd_check(i)= (sum((simulations_temp-exp_data(:,i)).^2)/length(exp_data))^0.5;
        w=length(diff);
        flag_temp = zeros(w,1);
        for j1=1:w
            if abs(diff(j1))>exp_err(j1,i)
                flag_temp(j1)=1;
            end
        end
        fraction(i)=(sum(~flag_temp)/w)*100;
        writematrix(table_temp,fullfile(pwd,'figure_s2_subplots',filename2{i}));
    end
    if exp_flag(i)>0
        rmse_string = sprintf('RMSE: %1.3f (%2.0f%%)',rmsd(i),round(fraction(i)));
        text(x_txt,y_txt,rmse_string,"Units","normalized","HorizontalAlignment","right","VerticalAlignment","baseline","FontSize",8);
    end
    set(gca,'Unit','Inches')
    p = get(gca,'Position');
    set(gca,'Unit','Inches','Position',[p(1) p(2) 1.75 1.25]);
    exportgraphics(figure(i),fullfile('figure_s2_subplots',filename{i}),'BackgroundColor','w','Resolution',300,'ContentType','vector');   
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
writetable(t_rmsd,fullfile(pwd,'figure_s2_subplots','RMSD_fig1.xlsx'),'WriteRowNames',true);
writetable(frac_table,fullfile(pwd,'figure_s2_subplots','frac_fig1.xlsx'),'WriteRowNames',true);
%plot power data
figure(q+1);clf;
errorbar(t6(1:5:240),power_exp(1:5:240),power_error(1:5:240), ...
            'LineStyle','none','Marker','.', 'MarkerSize',8,'Color',[1 1 1]*0.5,'CapSize',3);hold on;
plot(cycles,sim_power,'linewidth',2,'Color','k','DisplayName','Simulation'); 
rmse_string = sprintf('RMSE: %1.3f (%2.0f%%)',rmsd(q+1),round(fraction(q+1)));
text(1,0.9,rmse_string,"Units","normalized","HorizontalAlignment","right","VerticalAlignment","baseline","FontSize",8);
xlabel('Cycle Index');
ylabel('Power (W)');
set(gca,'Unit','Inches')
p = get(gca,'Position');
set(gca,'Unit','Inches','Position',[p(1) p(2) 1.75 1.25]);
exportgraphics(figure(q+1),fullfile(pwd,'figure_s2_subplots','power_mod.pdf'),'BackgroundColor','w','Resolution',300,'ContentType','vector');