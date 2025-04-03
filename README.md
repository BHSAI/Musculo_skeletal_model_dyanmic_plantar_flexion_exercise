# A cross-bridge cycling model for understanding the role of metabolite accumulation in skeletal muscle fatigue 
# Abstract
Skeletal muscle fatigue is accompanied by the accumulation of metabolites, such as adenosine diphosphate (ADP), inorganic phosphate (Pi), and protons (H+). However, we lack a comprehensive understanding of the contribution of these metabolic changes to the development of muscle fatigue during intense exercise and the underlying mechanisms. To address this gap, we collected data from young adults performing a dynamic (0.75 Hz) plantar flexion exercise to task failure (642±104 s), including in vivo concentrations of metabolites and H+ measured by 31P magnetic resonance spectroscopy as well as muscle activation signals obtained via electromyography. Using these data, we developed and validated a human skeletal muscle model. Our model-based simulations suggested that to continue the plantar flexion exercise at the required power output, muscle activation should progressively increase. In the absence of this increased activation, we observed a reduction in force-generating capacity due to metabolite-mediated inhibition of actin-myosin cross-bridge cycling. Our simulations also showed that Pi reduced force production by 30% when we increased it 50% above the concentrations measured experimentally. A parameter sensitivity analysis suggested that force generation is strongly dependent on the rate of Pi release from the actin-myosin complex, and Pi inhibits force by increasing the rate of actin-myosin detachment. In addition, we proposed an alternate mechanism through which H+ might reduce muscle force generation during exercise. In contrast, elevated ADP levels did not significantly affect force generation. This study provides insight into the impact of metabolite accumulation on force generation and muscle fatigue development.
# Contents
This repository contains the data, scripts, and models used in this study. They are organized into two folders 1) raw_data and 2) codes.
- Contents of folder /raw_data are as follows:
    - Power, EMG, Phosphocreatine concentration, Pi concentration , ADP concentration and pH dataset used for parameterization are in the excel files power_for_fitting_DPF_2.xlsx, Emg_for_fitting_DPF.xlsx, Pcr_for_fitting_DPF.xlsx, pi_for_fitting_DPF.xlsx, ADP_for_fitting_DPF.xlsx, and pH_for_fitting_DPF.xlsx, respectively.
    - Resting state concentrations and the sarcomere shortening velocity data used for parameterization are provided in the excel files Initial_state.xlsx and dsdt_for_fitting_DPF_2.xlsx, respectively.
    - The data set used for validation is provided in the folder val_dataset
- Contents of folder /codes are as follows:
    - This folder contains the MATLAB codes (MATLAB/R2022a) that represent the musculoskeletal model and the scripts used to simulate the model results shown in this manuscript
    - Model_XB_human_QC.m encodes the 4-state crossbridge model and the kinetic model discussed in the manuscript and is used to simulate the force and metabolite dynamics
    - params.xlsx in folder /params contains the model parameters estimated using our parameterization routine and subsequently used to simulate the model
    - figure_4_subplots.m simulates the model in Model_XB_human_QC.m  to generate the subplots of Figure 4
    - figure_5_subplots.m simulates the model in Model_XB_human_QC.m  to generate the subplots of Figure 5
    - figure_6_a.mlx, figure_6_b.mlx, figure_6_c.mlx and figure_6_d.mlx simulates the model in Model_XB_human_QC.m to generate the subplot A, B, C and D in Figure 6
    - figure_8_b.mlx, figure_8_c.mlx, figure_8_d.mlx and figure_8_e.mlx simulates the model in Model_XB_human_QC_SI.m to generate the subplot B, C, D and E in Figure 8. Model_XB_human_QC_SI.m is the model that incorporates the alternate proton inhibition hypothesis depicted in Figure 8A.
    - figure_s2_sublots.m simulates the model in Model_XB_human_QC_SI.m  to generate the subplots of Figure S2.

# Set up
- The codes were written and tested in MATLAB/R2022a and they can be directly executed from MATLAB command prompt or from the editor
    - Example: To run the code figure_4_subplots.m, directly type 'figure_4_subplots' in the MATLAB command prompt or open the editor and run it by pressing F5 or clicking the 'Run' icon in the toolbar of MATLAB Editor.
