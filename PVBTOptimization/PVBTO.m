% ======================================================== % 
% ---- PV-Battery Tool (PVBT) Sizing Optimization  ------- %
% Version: 1.2  (2/2021) --------------------------------- %
% ======================================================== % 
% This code was developed by Ahmed A.Raouf Mohamed: ------ %
% ------ AARaoufM@gmail.com / amohamed06@qub.ac.uk ------- %
% https://pure.qub.ac.uk/en/persons/ahmed-mohamed -------- %
% Copyright @2021 ---------------------------------------- %
% --- https://github.com/ARa2of/PV-BESS-Analysis-Tool  --- %
% ======================================================== % 
% ---Details are  available in the following publication-- %
% A. A. R. Mohamed, R. J. Best, X. Liu, and D. J. Morrow,  %
%‘A Comprehensive Robust Techno-Economic Analysis and      %
% Sizing Tool for the Small-Scale PV and BESS’ in IEEE --- %
%-------Transactions on Energy Conversion, 2021.-----------%
%----- doi: https://doi.org/10.1109/TEC.2021.3107103 ------%
% ======================================================== % 
% ---  This work is part of INTERREG VA SPIRE2 Project --- %
% EPIC Research Cluster ------- Queen's University Belfast %
% ======================================================== % 
% Please check PVBTGuide.pdf for details on how to run the %
% code. -------------------------------------------------- %
% The measurements should be entered in inputs.csv ------- %
% -------------------------------------------------------- %
% ======================================================== %      

clc; close all; clear;
tic;
format long g
global Select_Optimizer Tariff kPF DataRes IRR OPTTY PVSize BESS BESSP DOD SOCMAX SOCMIN SOCI RE PVCost InvCost InvSize Lifinv Lifpv PVdeg PVOM PRP EXP EX YearI IR er Dr PCN TLS TLE PTHC PTHD LTY SOHM BP SC SaveR BPP ERRa FB DIAA Profile
%% 
Tariff="DT"; % DT for dual tariff (Economy 7 in this code), FT for flat tariff, TT for triple tariff (TIDE tariff in this code)
DataRes=30; %Data resolution 10 for 10 minutes reso, 30 for 30 minutes reso, 60 for 60 minutes(1 hour) reso and so on...
OPTTY="BPV"; %Determine the optimization type : BPV: find the BESS and PV sizes together, B: Optimize the BESS size only, PV: Optimize the PV size only. 
Select_Optimizer="GreyWolfOptimizer"; % Select from the following options 1) NOMAD, 2) GreyWolfOptimizer 3)PatternSearch

if OPTTY=="PV"
BESS=1E-10; %Actual BESS Capacity [kWh] : If you want to size the PV in presence of a BESS, add the BESS specs. If not, leave it as small value (not zero). 
BESSP=1E-10; %BESS Power rating [kW]
LPV = 1; %PV Lower limit of the search space 
UPV = 6; %PV Upper limit of the search space 
end
if OPTTY=="B"
PVSize=3.3;  % Add the size of the PV to calcualte the BESS accordingly   
LBS=2.4; %BESS Lower limit of the search space 
UBS=14;  %BESS Upper limit of the search space 
end
if OPTTY=="BPV"
LPV=1;   %PV Lower limit of the search space 
UPV=6;   %PV Upper limit of the search space 
LBS=2.4; %BESS Lower limit of the search space 
UBS=20;  %BESS Upper limit of the search space 
end
%% Load Profiles
Profile=readmatrix('InputsHalfHourly.csv');
D=Profile(:,1); %Demand
PV=Profile(:,2);
EV=Profile(:,3);
T=length(Profile(:,1));
ND=round((DataRes/60)*(T/24)); %Number of days
TD=T/ND; %Length of one day
tau=TD/(24); % {Time interval=1/tau}
%% BESS Inputs
DOD=0.95; %MAX DOD
SOCMAX=1; %Max SoC
LTY=15; % lifetime in years
SOHM=0.6; % Minimum State of Health 
costbess=500;
FB=costbess*15/100; %Fixed price of the BESS price that doesn't decline 
BP=costbess-FB; % BESS Price  £/kWh%%%%%%%%%%%%%
RE=0.95*0.95; %= 0.95(BESS) * 0.95(Inverter)
SOCMIN=SOCMAX-DOD; %Min SoC 
SOCI=SOCMIN; %Initial SOC that the simulations will start with. 
%% PV Inputs
PVCost=1742; %Cost in £/kW
InvCost=100; %inverter Cost in £/kW
InvSize=3.68; %inverter Size in kW
Lifinv=15; %Lifetime of inverter in years
Lifpv=30; %Lifetime of PV  in years
PVdeg=0.5/100; %PV annual degradation rate
PVOM=1/100; %PV annual O&M cost as a percentage of the capital cost 
%% Utility Inputs - Tariffs are in pence/kWh or cent/kWh
EXP=3.68; %Export Power Limit =3.68kW
% Export Tariff e.g. PowerNI Microgen Tariff (https://powerni.co.uk/products--services/renewableenergy/sell-electricity/)
EX=4.59; % Export Tariff = 4.59p/kWh 
if Tariff=="DT"
%Economy 7 Dual tariff.
SC=13.66; % standing charge £/day
HR=35.89; %Day Rate 8am-1am 
LR=17.16; % Night Rate 1am-8am 
TLS=1;  %Night rate start time
TLE=7; %Night Rate end time
TPR=[LR LR LR LR LR LR LR HR HR HR HR HR HR HR HR HR HR HR HR HR HR HR HR HR];%Tariff Profile - Economy 7
FCN=1;
end
if Tariff=="FT"
% Flat Tariff Data 
SC=13.66; % standing charge £/day
SR = 31; % Taken from PowerNI Website - 16.77p/kWh via Online Billing, Monthly Direct Debit incl VAT
TPR=[SR SR SR SR SR SR SR SR SR SR SR SR SR SR SR SR SR SR SR SR SR SR SR SR];%Tariff Profile - Flat Tariff
FCN=0;
end
if Tariff=="TT"
% TIDE Tariff  - Triple tariff
% (https://www.greenenergyuk.com/TariffInfoLabel.aspx?TARIFF_ID=4&IS_TWO_RATE=False&IS_DUAL_FUEL=True&GAS=False&ELECTRICITY=True&GSP_GROUP=_A&POSTCODE=AL1%203EZ)
SC=22; % standing charge £/day
HR=32.55; %high Rate 8am-1am 
MR=16.27; %Medium rate (Shoulder)
LR=7.91; % low Rate 1am-8am 
TLS=1;  %Low rate start time
TLE=6; %Night Rate end time
TPR=[LR LR LR LR LR LR MR MR MR MR MR MR MR MR MR MR HR HR HR MR MR MR MR LR];%Tariff Profile - TIDE
FCN=1;
end
PRP = repelem(TPR,tau);
%% Cost Benefit Analysis 
YearI=2025; %Year of installation
IR=4.5/100; %Interest Rate
er=2/100; %electricity annual growth rate
Dr=12/100;% Annual declining rate in BESS prices
%% BESS Control Strategy Inputs: Threshold Rule-based 
%% BESS Control Method Inputs: Threshold Rule-based 
PTHD=0; % Specify the Upper threshold for BESS Discharge
PTHC=0; % Specify the Lower threshold for BESS Charge  
PCNS=[0.7, 0.23, 0, 0.45]*FCN; %A percentage of the BESS capacity to be charged during each season ((Winter, Spring, Summer, Autumn)
DIA=0; %= 1: if you want the algorithm to start discharging after the end of low tariff period, =0: to discharge whenever the demand exceeds the generation at any time of the day
%using low tariff rate to maximize the energy arbitrage. Set all values to zero if you don’t want to use this feature 

MAINCODE0

