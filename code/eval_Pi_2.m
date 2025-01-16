function [Y] = eval_Pi_2(params, iemg_prof, cycles, cut_off, Pi_set, dpidt_set, dHdt_set, dMgADPdt_set, dPCrdt_set)
data_resting=readtable('../raw_data/Initial_state.xlsx','Sheet','Summary'); % resting levels of state variables
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
cycle_time=10/6.33;% average of 6.33 cycles for 10 secs calculated from the data
tspan = 0:0.1:cycle_time;
n=length(tspan);
m=cycles;
sim_Ftotal = zeros(m,1);
sim_Ftotal_cycles=zeros(n,m);
pi_p=zeros(m,1);
ADP_p=zeros(m,1);
Pcr_p=zeros(m,1);
H_p=zeros(m,1);
ATP_p=zeros(m,1);
for i=1:m
    SL_set=3.23;
    dSL_set=-0.68; %takes into account the contraction in G.medialis alone
    iemg=iemg_prof(i);
    options = odeset('RelTol',1e-3,'AbsTol',1e-6,'MaxStep',5e-3);
    if i<=cut_off 
        [T, Y] = ode15s(@Model_XB_human_QC,tspan,init,options,TmpC,SL_set,params,iemg,dSL_set,Pcr,H);
        %init(13)=Y(n,13);%ADP
        init(12)=Y(n,12);%P
    else
        init(12) = Pi_set*pi_p(cut_off);
        [T, Y] = ode15s(@Model_XB_human_QC_metdyn_set,tspan,init,options,TmpC,SL_set,params,iemg,dSL_set,Pcr,H,dpidt_set,dHdt_set,dMgADPdt_set,dPCrdt_set);
    end
    init(10)=Y(n,10);%N
    %init(12)=Y(n,12);%P
    init(14)=Y(n,14);%Pcr
    init(15)=Y(n,15);%H
    init(16)=Y(n,16);%ATP
    init(13)=Y(n,13);%ADP
    pi_p(i)=Y(n,12);
    ADP_p(i)=Y(n,13);
    Pcr_p(i)=Y(n,14);
    H_p(i)=Y(n,15);
    ATP_p(i)=Y(n,16);
    for j=1:n
        [~, sim_Ftotal_cycles(j,i),~,~,~,~,~,~,~,~,~,~,~,~,~] = Model_XB_human_QC(T(j),Y(j,:),TmpC,SL_set,params,iemg,dSL_set,Pcr,H);
    end
    sim_Ftotal(i) = sim_Ftotal_cycles(n,i);
end
Y = [pi_p ADP_p Pcr_p H_p ATP_p sim_Ftotal];
end