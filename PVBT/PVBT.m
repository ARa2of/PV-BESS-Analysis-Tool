% ======================================================== % 
% --------------   PV-Battery Tool (PVBT)   -------------- %
% Version: 1.1  (1/2021) --------------------------------- %
% ======================================================== % 
% This code was developed by Ahmed A.Raouf Mohamed: ------ %
% ------ Ra2ooof@gmail.com / amohamed06@qub.ac.uk -------- %
% https://pure.qub.ac.uk/en/persons/ahmed-mohamed -------- %
% Copyright @2021 ---------------------------------------- %
% --- https://github.com/ARa2of/PV-BESS-Analysis-Tool  --- %
% ======================================================== % 
% Details will be available in the following publication-- %
% A. A. R. Mohamed, R. J. Best, X. Liu, and D. J. Morrow,  %
%‘A Comprehensive Robust Techno-Economic Analysis and      %
% Sizing Tool for the Small-Scale PV and BESS’------------ %
% EPIC Research Cluster ------- Queen's University Belfast %
% ======================================================== % 
% ---  This work is part of INTERREG VA SPIRE2 Project --- %
% ======================================================== % 
% Please check PVBTGuide.pdf for details on how to run the %
% code. -------------------------------------------------- %
% The measurements should be entered in inputs.csv ------- %
% -------------------------------------------------------- %
% ======================================================== %    

%%
clc; close all; clear;
tic;
format long g
%% Main Inputs 
SaveR=0; %if 1=save results in excel files, other values=don't save : (saving results will reduce excution time) 
DataRes=30; %Data resolution 10 for 10 minutes reso, 30 for 30 minutes reso, 60 for 60 minutes(1 hour) reso and so on...
Tariff="TT"; % DT for dual tariff (Economy 7 in this code), FT for flat tariff, TT for triple tariff (TIDE tariff in this code)
%% Load Profiles
Profile=readmatrix('Inputs.csv'); % Insert the simulation data as explained in the Inputs.csv file
D=Profile(:,1); %Demand
PV=Profile(:,2);
EV=Profile(:,3);
T=length(Profile(:,1));
ND=round((DataRes/60)*(T/24)); %Number of days
TD=T/ND; %Length of one day
tau=TD/(24); % {Time interval=1/tau}
%% BESS Inputs
BESS=6.5; %Actual BESS Capacity [kWh]
DOD=0.945; %MAX DOD
SOCMAX=1; %Max SoC
BESSP=3.3; %BESS Power rating [kW]
RE=0.95*0.95; %= 0.95(BESS) * 0.95(Inverter)
SOHM=60/100; % Minimum State of Health e.g. 60% 
FB=70.875; %Fixed price of the BESS price that doesn't decline 
BP=436.59; % BESS Price  £/kWh
LTY=10; % lifetime in years
%%%%%%%%%%%%%
SOCMIN=SOCMAX-DOD; %Min SoC 
SOCI=SOCMIN; %Initial SOC that the simulations will start with. 
%% PV Inputs
PVSize=3.3; %Size in kWp
PV=PV.*PVSize;
PVCost=1400; %Cost in £/kW
Lifpv=30; %Lifetime of PV  in years
InvCost=100; %inverter Cost in £/kW
InvSize=3.68; %inverter Size in kW
Lifinv=15; %Lifetime of inverter in years
PVdeg=0.5/100; %PV annual degradation rate
PVOM=1/100; %PV annual O&M cost as a percentage of the capital cost 
%% Utility Inputs - Tariffs are in pence/kWh or cent/kWh
EXP=3.68; %Export Power Limit =3.68kW
% Export Tariff 
EX=5; % Export Tariff = 5 p/kWh 
if Tariff=="FT"
% Flat Tariff Data
SC=22; % standing charge £/day
SR = 16.6; % Taken from PowerNI Website - 16.77p/kWh via Online Billing, Monthly Direct Debit incl VAT
TPR=[SR SR SR SR SR SR SR SR SR SR SR SR SR SR SR SR SR SR SR SR SR SR SR SR];%Tariff Profile - Flat Tariff
FCN=0;
end
if Tariff=="DT"
%Economy 7 Dual tariff
SC=22; % standing charge £/day
HR=17.4; %Day Rate 8am-1am 
LR=7.91; % Night Rate 1am-8am 
TLS=1;  %Night rate start time
TLE=7; %Night Rate end time
TPR=[LR LR LR LR LR LR LR HR HR HR HR HR HR HR HR HR HR HR HR HR HR HR HR HR];%Tariff Profile - Economy 7
FCN=1;
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
YearI=2021; %Year of installation
IR=3.5/100; %Interest Rate
er=2/100; %electricity price annual growth rate
Dr=12/100;% Annual declining rate in BESS prices
%% BESS Control Method Inputs: Threshold Rule-based 
PTHD=0; % Specify the Upper threshold for BESS Discharge
PTHC=0; % Specify the Lower threshold for BESS Charge 
PCNS=[0.7, 0.23, 0, 0.45]*FCN; %A percentage of the BESS capacity to be charged during each season ((Winter, Spring, Summer, Autumn)
%using low tariff rate to maximize the energy arbitrage. Set all values to zero if you don’t want to use this feature 

MAINCODE
