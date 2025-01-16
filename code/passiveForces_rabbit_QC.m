function [FPEE,FSEE] = passiveForces_rabbit_QC(SL, SLset)
%% Passive Force formulation from Rockenfeller et al 2020
% Function to calculate the passive forces of the muscle
% The passive force is a function of the sarcomere length
lm = 3.8462e+04*SLset*1e-6; % m (optimal fiber lenth 10 cm /2.6 um sarcomere length (10e-2/2.6e-6 = 3.8462e+04)
% lm = 7.7225e4*SLset*1e-6; % m
% lm =0.114;
SLopt = 2.6; %um
% lCEopt=0.0198;  %[m] Optimal length of muscle fibres
% lCEopt=5e3*SLopt*1e-6; 
lCEopt=3.8462e+04*SLopt*1e-6;
% lCEopt=8.7225e3*SLopt*1e-6; 
LLPEE0=1*3.02;
vPEE=1*3.16;
dWdes=0.330;
FFPEE=0.135;
% Fmax=3400; % N % Soleus
% Fmax = 675; %N % Tibialis Anterior
Fmax = 6000;% N % vasta lateralis (Quadriceps)
lCE_temp = 3.8462e+04*SL*1e-6;
% lCE_temp = 8.7225e3*SL*1e-6;
KPEE=FFPEE.*Fmax./((lCEopt.*(dWdes+1-LLPEE0)).^vPEE);
lPEE0=LLPEE0.*lCEopt;
% lPEE0=LLPEE0.*SLopt;
% FPEE=0;
        if lCE_temp<lPEE0
            FPEE=0;
        else
            FPEE= KPEE.*(lCE_temp-lPEE0).^vPEE;
        end
        
% lSEE0=0.1017/120;
lSEE0=1*0.1017/80;
dUSEEnll=0.0785 ;
dUSEEl=1*0.0482;
% dFSEE0= 2700; % N Sol
% dFSEE0= 250; % N TA
dFSEE0= 2500; % N VL
lSEEnll=(1+dUSEEnll).*lSEE0;
vSEE=dUSEEnll./dUSEEl;
KSEEnl=dFSEE0./(dUSEEnll.*lSEE0).^vSEE;
KSEEl=dFSEE0./(dUSEEl.*lSEE0);

        if (lm-lCE_temp)<lSEE0
            FSEE=0;
        elseif (lm-lCE_temp)<lSEEnll
            FSEE= KSEEnl.*(lm-lCE_temp-lSEE0).^vSEE;
        else
            FSEE=(dFSEE0+KSEEl.*(lm-lCE_temp-lSEEnll));
        end
        
end