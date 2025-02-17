if DIA==1
    DIAA=tau;
else
    DIAA=1;
end

PCN=zeros;
for g=1:ND
   if g>=1 && g<60 || g>335 && g<=365
       PCN(g)=PCNS(1);
   elseif g>=60 && g<150
        PCN(g)=PCNS(2);
           elseif g>=150 && g<240
        PCN(g)=PCNS(3);
           elseif g>=240 && g<=335
        PCN(g)=PCNS(4);
   end
end
lttyy=LTY;
NuB=ceil(Lifpv/lttyy);
BPF=zeros;

%%%%+++%%%%%
for i=1:LTY*(NuB)
   BPF(i)=((1-(Dr))^i)*BP;
end
if NuB==1
    BPP=BP;
end
if NuB==2
BPP=[BP BPF(LTY-1)];
end
if NuB==3
    BPP=[BP BPF(LTY-1) BPF((LTY*2)-1)];
end
if NuB==4
    BPP=[BP BPF(LTY-1) BPF((LTY*2)-1) BPF((LTY*3)-1)];
end
if NuB==5
    BPP=[BP BPF(LTY-1) BPF((LTY*2)-1) BPF((LTY*3)-1) BPF((LTY*4)-1)];
end

for i=1:(Lifpv*3) 
   ERRa(i)=(1+er)^(i-1); 
   IRR(i)=((1+IR)^(i-1));
end
%%%%+++%%%%%

if Lifpv/LTY <NuB
    Res=NuB-Lifpv/LTY;
elseif Lifpv/LTY == NuB
    Res=0;
end
TSM=ones(1,NuB);
TSM(end)=1-Res;
LYY=lttyy*TSM;
% EFGO=cumsum(LYY);
% GF=[1 EFGO(1:end-1)];
GF=zeros(1,NuB);
for TTO=1:ceil(Lifpv/LTY)
   YearIn(TTO)= YearI+TTO*LTY;
end
YearIns=[YearI YearIn(1:end-1)];



%% Optimization Part
if OPTTY=="BPV"
    OPTTY="PV";
fun=@MAINCODE;
kPF=1;
BESS=1E-10; %Actual BESS Capacity [kWh]
BESSP=1E-10; %BESS Power rating [kW]
lb = [LPV];
ub = [UPV];
x0 = [LPV];
A = [];
b = [];
Aeq = [];
beq = [];
nvars=length(lb);

if Select_Optimizer=="NOMAD"
	opts = optiset('solver','NOMAD','display','iter');
	Opt = opti('fun',fun,'bounds',lb,ub,'options',opts,'ndec',nvars,'x0',x0);
	[x,objec] = solve(Opt,x0);
end	
if Select_Optimizer=="Fmincon"
	x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub);
end	
if Select_Optimizer=="PatternSearch"
    x = patternsearch(fun,x0,A,b,Aeq,beq,lb,ub);
end	
	
PVSize=round(x(1),1); %Size in kW
    OPTTY="B";
fun=@MAINCODE;
lb2 = [LBS];
ub2 = [UBS];
x0 = [LBS];
A = [];
b = [];
Aeq = [];
beq = [];
nvars=length(lb);
for i=1:NuB
kPF=GF(i);
BP=BPP(i);
LTY=single(LYY(i));
if Select_Optimizer=="NOMAD"
	opts = optiset('solver','NOMAD','display','iter');
	Opt = opti('fun',fun,'bounds',lb,ub,'options',opts,'ndec',nvars,'x0',x0);
	[x,objec] = solve(Opt,x0);
end	
if Select_Optimizer=="Fmincon"
	x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub);
end	
if Select_Optimizer=="PatternSearch"
    x = patternsearch(fun,x0,A,b,Aeq,beq,lb,ub);
end	

BG(i)=x;
end
BES=round((BG),2); %Actual BESS Capacity [kWh]
BESP=round(1.245+0.304*(BG),2); %BESS Power rating [kW]
    OPTTY="BPV";
end

if OPTTY=="B"
fun=@MAINCODE;
lb = [LBS];
ub = [UBS];
x0 = [LBS];
A = [];
b = [];
Aeq = [];
beq = [];
nvars=length(lb);
for i=1:NuB
kPF=GF(i);
BP=BPP(i);
LTY=LYY(i);
if Select_Optimizer=="NOMAD"
	opts = optiset('solver','NOMAD','display','iter');
	Opt = opti('fun',fun,'bounds',lb,ub,'options',opts,'ndec',nvars,'x0',x0);
	[x,objec] = solve(Opt,x0);
end	
if Select_Optimizer=="Fmincon"
	x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub);
end	
if Select_Optimizer=="PatternSearch"
    x = patternsearch(fun,x0,A,b,Aeq,beq,lb,ub);
end	

BG(i)=x;
end
BES=round((BG),2); %Actual BESS Capacity [kWh]
BESP=round(1.245+0.304*(BG),2); %BESS Power rating [kW]
end

if OPTTY=="PV"

fun=@MAINCODE;
lb = [LPV];
ub = [UPV];
x0 = [LPV];
kPF=1;
A = [];
b = [];
Aeq = [];
beq = [];
nvars=length(lb);
if Select_Optimizer=="NOMAD"
	opts = optiset('solver','NOMAD','display','iter');
	Opt = opti('fun',fun,'bounds',lb,ub,'options',opts,'ndec',nvars,'x0',x0);
	[x,objec] = solve(Opt,x0);
end	
if Select_Optimizer=="Fmincon"
	x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub);
end	
if Select_Optimizer=="PatternSearch"
    x = patternsearch(fun,x0,A,b,Aeq,beq,lb,ub);
end	
end
%%

PV=PVSize.*PV;
if Tariff=="FT"
TLE=1;
TLS=1;
end
%% Repeat the simulations with the optimal results
TSS=zeros(1,NuB);NPVVV=zeros(1,NuB);

if OPTTY=="PV"
    k=1;
BCM
close
end

if OPTTY=="B"
    BI=zeros(1,NuB);BII=zeros(1,NuB);BESS=zeros(1,NuB);BESSP=zeros(1,NuB);BP=zeros(1,NuB);BILL1=zeros;BILL2=zeros;
    for jjy=1:NuB
BESS=BES(jjy);
BESSP=BESP(jjy);
BP=BPP(jjy);
kPF=GF(jjy);
LTY=single(LYY(jjy));
    BCM
    BI=BILL1;
    BII(jjy)=BILL2;
    TSS(jjy)=sum(TS1);
    NPVVV(jjy)=NPV1;
    aroii(jjy)=AROI1;
    SOHHY(jjy)=SOHY(LTY)*100;
    PBY(jjy)=(PB1);
    SBY(jjy)=SB1;
    LTYNY(jjy)=LTYN;
        SCCC(jjy)=SC2;
    SSS(jjy)=SS2;
    EERRR(jjy)=sum(EErr*-1);
    PCUAA(jjy)=PCUA;
    end
end

if OPTTY=="BPV"
%     BCM
    BI=zeros(1,NuB);BII=zeros(1,NuB);BESS=zeros(1,NuB);BESSP=zeros(1,NuB);BP=zeros(1,NuB);BILL1=zeros;BILL2=zeros;SOHHY=zeros;
    for jjy=1:NuB
BESS=BES(jjy);
BESSP=BESP(jjy);
BP=BPP(jjy);
kPF=GF(jjy);
LTY=single(LYY(jjy));
BCM
    BI=BILL1;
    BII(jjy)=BILL2;
    TSS(jjy)=sum(TS1);
    NPVVV(jjy)=NPV1;
    aroii(jjy)=AROI1;
    SOHHY(jjy)=SOHY(LTY)*100;
    PBY(jjy)=(PB1);
    SBY(jjy)=SB1;
    LTYNY(jjy)=LTYN;
    SCCC(jjy)=SC2;
    SSS(jjy)=SS2;
    EERRR(jjy)=sum(EErr*-1);
    PCUAA(jjy)=PCUA;
    end
end


%% Print
disp('      ')
disp('=============================================================================================')
disp('********************                Building Annual Data               **********************')
disp('=============================================================================================')
disp(['Total Annual Consumption:              ',num2str(round((sum(D)+sum(EV))*1/tau,3)),' kWh/year'])
disp(['Total PV Annual Generation:            ',num2str(round(sum(PV)*1/tau,3)),' kWh/year'])
disp(['Tariff Structure:                      ',num2str(Tariff)])

if OPTTY=="PV"
disp('=============================================================================================')
disp('******************** SIMULATION RESULTS OF THE PVBT (PV sizing w/o BESS) ********************')
disp('=============================================================================================')
disp(['Optimal PV size:               ',num2str(PVSize),' kWp'])
disp('-------------------------------------- One-Year Results -------------------------------------')
disp('---------------------------------------------------------------------------------------------')
disp(['Electricity Bill w/o PV w/o BESS:       ',num2str(round(BILL0,3)),' £'])
disp(['Electricity Bill w/ PV w/o BESS:        ',num2str(round(BILL1,3)),' £'])
disp(['Curtailed PV power w/o BESS:            ',num2str(round(PCUB,3)),' kWh'])
disp(['Exported PV power w/o BESS:             ',num2str(round(sum(EEr*-1),3)),' kWh'])
disp(['PV Self-Consumption w/o BESS:           ',num2str(round(SC1,3)),' %'])
disp(['Premises Self-Sufficiency w/o BESS:     ',num2str(round(SS1,3)),' %'])
disp('=============================================================================================')
disp(['++++++++++ PV Cost Benefit Analysis for a Lifetime of ',num2str(Lifpv),' Years  +++++++++++++'])
disp('---------------------------------------------------------------------------------------------')
disp(['Annual Saving:                    ',num2str(round(sum(BILL0-BILL1),3)),' £'])
disp(['Total Savings:                    ',num2str(round(sum(TS11),3)),' £'])
disp(['Net Present Value:                ',num2str(round(NPV11,3)), ' £'])
disp(['Annual Return on Investment:      ',num2str(round(AROI11,3)),' %'])
disp(['Discounted Payback Period:        ',num2str(PB11(1)),' Years'])
disp(['Subsidy Required:                 ',num2str((SB11)), ' £/kW'])
disp('=============================================================================================')
end

if OPTTY=="B"
disp('=============================================================================================')
disp('************ SIMULATION RESULTS OF THE PVBT (BESS sizing in the presence of PV) *************')
disp('=============================================================================================')
disp(['Optimal BESS sizes:               ',num2str(BES),' kWh'])
disp(['BESS ratings :                    ',num2str(BESP),' kW'])
disp(['Years of installations:           ',num2str(YearIns)])
disp(['BESS price for each year:         ',num2str(round(BPP+FB)), ' £/kWh'])
if NuB>=1
disp('=============================================================================================')
disp(['+++++++++++++++++++++++++++ First BESS Cost Benefit Analysis  +++++++++++++++++++++++++++++++'])
disp('---------------------------------------------------------------------------------------------')
disp(['Annual Savings:                   ',num2str(round((BI-BII(1)),3)),' £/year'])
disp('---------------------------------------------------------------------------------------------')
disp(['--------------------  Cost Benefit Analysis for a operation of ',num2str(LYY(1)),' Years  ---------------------'])
disp('---------------------------------------------------------------------------------------------')
disp(['Curtailed PV power:               ',num2str(round(PCUAA(1),3)),' kWh'])
disp(['Exported PV power:                ',num2str(round(EERRR(1),3)),' kWh'])
disp(['PV Self-Consumption:              ',num2str(round(SCCC(1),3)),' %'])
disp(['Premises Self-Sufficiency:        ',num2str(round(SSS(1),3)),' %'])
disp(['Total Savings:                    ',num2str(round(TSS(1),3)),' £'])
disp(['Net Present Value:                ',num2str(round(NPVVV(1),3)), ' £'])
disp(['Annual Return on Investment:      ',num2str(round(aroii(1),3)),' %'])
disp(['SoH at the end of Lifetime:       ',num2str(round(SOHHY(1),3)),' %'])
disp(['Operation years to minimum SoH:   ',num2str(LTYNY(1)),' Years'])
disp(['Discounted Payback Period:        ',num2str(PBY(1)),' Years'])
disp(['Subsidy Required:                 ',num2str((SBY(1))), ' £/kWh'])
end
if  NuB>=2
disp('=============================================================================================')
disp(['+++++++++++++++++++++++++++ Second BESS Cost Benefit Analysis  ++++++++++++++++++++++++++++++'])
disp('---------------------------------------------------------------------------------------------')
disp(['Annual Savings:                   ',num2str(round((BI-BII(2)),3)),' £/year'])
disp('---------------------------------------------------------------------------------------------')
disp(['-------------------- Cost Benefit Analysis for a operation of  ',num2str(LYY(2)),' Years  ---------------------'])
disp('---------------------------------------------------------------------------------------------')
disp(['Curtailed PV power:               ',num2str(round(PCUAA(2),3)),' kWh'])
disp(['Exported PV power:                ',num2str(round(EERRR(2),3)),' kWh'])
disp(['PV Self-Consumption:              ',num2str(round(SCCC(2),3)),' %'])
disp(['Premises Self-Sufficiency:        ',num2str(round(SSS(2),3)),' %'])
disp(['Total Savings:                    ',num2str(round(TSS(2),3)),' £'])
disp(['Net Present Value:                ',num2str(round(NPVVV(2),3)), ' £'])
disp(['Annual Return on Investment:      ',num2str(round(aroii(2),3)),' %'])
disp(['SoH at the end of Lifetime:       ',num2str(round(SOHHY(2),3)),' %'])
disp(['Operation years to minimum SoH:   ',num2str(LTYNY(2)),' Years'])
disp(['Discounted Payback Period:        ',num2str(PBY(2)),' Years'])
disp(['Subsidy Required:                 ',num2str((SBY(2))), ' £/kWh'])
    end
    if  NuB>=3
disp('=============================================================================================')
disp(['+++++++++++++++++++++++++++ Third BESS Cost Benefit Analysis  +++++++++++++++++++++++++++++++'])
disp('---------------------------------------------------------------------------------------------')
disp(['Annual Savings:                   ',num2str(round((BI-BII(3)),3)),' £/year'])
disp('---------------------------------------------------------------------------------------------')
disp(['-------------------- Cost Benefit Analysis for a operation of ',num2str(LYY(1)),' Years  ---------------------'])
disp('---------------------------------------------------------------------------------------------')
disp(['Curtailed PV power:               ',num2str(round(PCUAA(3),3)),' kWh'])
disp(['Exported PV power:                ',num2str(round(EERRR(3),3)),' kWh'])
disp(['PV Self-Consumption:              ',num2str(round(SCCC(3),3)),' %'])
disp(['Premises Self-Sufficiency:        ',num2str(round(SSS(3),3)),' %'])
disp(['Total Savings:                    ',num2str(round(TSS(3),3)),' £'])
disp(['Net Present Value:                ',num2str(round(NPVVV(3),3)), ' £'])
disp(['Annual Return on Investment:      ',num2str(round(aroii(3),3)),' %'])
disp(['SoH at the end of Lifetime:       ',num2str(round(SOHHY(3),3)),' %'])
disp(['Operation years to minimum SoH:   ',num2str(LTYNY(3)),' Years'])
disp(['Discounted Payback Period:        ',num2str(PBY(3)),' Years'])
disp(['Subsidy Required:                 ',num2str((SBY(3))), ' £/kWh'])
        end
end




if OPTTY=="BPV"
disp('=============================================================================================')
disp('**************SIMULATION RESULTS OF THE PVBT Sizing Optimization (PV and BESS)***************')
disp('=============================================================================================')
disp(['Optimal PV size:               ',num2str(PVSize),' kWp'])
disp('-------------------------------------- One-Year Results -------------------------------------')
disp('---------------------------------------------------------------------------------------------')
disp(['Electricity Bill w/o PV w/o BESS:       ',num2str(round(BILL0,3)),' £'])
disp(['Electricity Bill w/ PV w/o BESS:        ',num2str(round(BILL1,3)),' £'])
disp(['Curtailed PV power w/o BESS:            ',num2str(round(PCUB,3)),' kWh'])
disp(['Exported PV power w/o BESS:             ',num2str(round(sum(EEr*-1),3)),' kWh'])
disp(['PV Self-Consumption w/o BESS:           ',num2str(round(SC1,3)),' %'])
disp(['Premises Self-Sufficiency w/o BESS:     ',num2str(round(SS1,3)),' %'])
disp('=============================================================================================')
disp(['++++++++++++++++++ PV Cost Benefit Analysis for a Lifetime of ',num2str(Lifpv),' Years  +++++++++++++++++++++'])
disp('---------------------------------------------------------------------------------------------')
disp(['Annual Saving:                    ',num2str(round(sum(BILL0-BILL1),3)),' £'])
disp(['Total Savings:                    ',num2str(round(sum(TS11),3)),' £'])
disp(['Net Present Value:                ',num2str(round(NPV11,3)), ' £'])
disp(['Annual Return on Investment:      ',num2str(round(AROI11,3)),' %'])
disp(['Discounted Payback Period:        ',num2str(PB11(1)),' Years'])
disp(['Subsidy Required:                 ',num2str((SB11)), ' £/kW'])
disp('=============================================================================================')
disp('=============================================================================================')
disp(['Optimal BESS sizes:               ',num2str(BES),' kWh'])
disp(['BESS ratings :                    ',num2str(BESP),' kW'])
disp(['Years of installations:           ',num2str(YearIns)])
disp(['BESS price for each year:         ',num2str(round(BPP+FB)), ' £/kWh'])
if NuB>=1
disp('=============================================================================================')
disp(['+++++++++++++++++++++++++++ First BESS Cost Benefit Analysis  +++++++++++++++++++++++++++++++'])
disp('---------------------------------------------------------------------------------------------')
disp(['Annual Savings:                   ',num2str(round((BI-BII(1)),3)),' £/year'])
disp('---------------------------------------------------------------------------------------------')
disp(['--------------------  Cost Benefit Analysis for a operation of ',num2str(LYY(1)),' Years  --------------------'])
disp('---------------------------------------------------------------------------------------------')
disp(['Curtailed PV power:               ',num2str(round(PCUAA(1),3)),' kWh'])
disp(['Exported PV power:                ',num2str(round(EERRR(1),3)),' kWh'])
disp(['PV Self-Consumption:              ',num2str(round(SCCC(1),3)),' %'])
disp(['Premises Self-Sufficiency:        ',num2str(round(SSS(1),3)),' %'])
disp(['Total Savings:                    ',num2str(round(TSS(1),3)),' £'])
disp(['Net Present Value:                ',num2str(round(NPVVV(1),3)), ' £'])
disp(['Annual Return on Investment:      ',num2str(round(aroii(1),3)),' %'])
disp(['SoH at the end of Lifetime:       ',num2str(round(SOHHY(1),3)),' %'])
disp(['Operation years to minimum SoH:   ',num2str(LTYNY(1)),' Years'])
disp(['Discounted Payback Period:        ',num2str(PBY(1)),' Years'])
disp(['Subsidy Required:                 ',num2str((SBY(1))), ' £/kWh'])
end
if  NuB>=2
disp('=============================================================================================')
disp(['+++++++++++++++++++++++++++ Second BESS Cost Benefit Analysis  ++++++++++++++++++++++++++++++'])
disp('---------------------------------------------------------------------------------------------')
disp(['Annual Savings:                   ',num2str(round((BI-BII(2)),3)),' £/year'])
disp('---------------------------------------------------------------------------------------------')
disp(['-------------------- Cost Benefit Analysis for a operation of  ',num2str(LYY(2)),' Years  --------------------'])
disp('---------------------------------------------------------------------------------------------')
disp(['Curtailed PV power:               ',num2str(round(PCUAA(2),3)),' kWh'])
disp(['Exported PV power:                ',num2str(round(EERRR(2),3)),' kWh'])
disp(['PV Self-Consumption:              ',num2str(round(SCCC(2),3)),' %'])
disp(['Premises Self-Sufficiency:        ',num2str(round(SSS(2),3)),' %'])
disp(['Total Savings:                    ',num2str(round(TSS(2),3)),' £'])
disp(['Net Present Value:                ',num2str(round(NPVVV(2),3)), ' £'])
disp(['Annual Return on Investment:      ',num2str(round(aroii(2),3)),' %'])
disp(['SoH at the end of Lifetime:       ',num2str(round(SOHHY(2),3)),' %'])
disp(['Operation years to minimum SoH:   ',num2str(LTYNY(2)),' Years'])
disp(['Discounted Payback Period:        ',num2str(PBY(2)),' Years'])
disp(['Subsidy Required:                 ',num2str((SBY(2))), ' £/kWh'])
    end
    if  NuB>=3
disp('=============================================================================================')
disp(['+++++++++++++++++++++++++++ Third BESS Cost Benefit Analysis  +++++++++++++++++++++++++++++++'])
disp('---------------------------------------------------------------------------------------------')
disp(['Annual Savings:                   ',num2str(round((BI-BII(3)),3)),' £/year'])
disp('---------------------------------------------------------------------------------------------')
disp(['-------------------- Cost Benefit Analysis for a operation of ',num2str(LYY(3)),' Years  ---------------------'])
disp('---------------------------------------------------------------------------------------------')
disp(['Curtailed PV power:               ',num2str(round(PCUAA(3),3)),' kWh'])
disp(['Exported PV power:                ',num2str(round(EERRR(3),3)),' kWh'])
disp(['PV Self-Consumption:              ',num2str(round(SCCC(3),3)),' %'])
disp(['Premises Self-Sufficiency:        ',num2str(round(SSS(3),3)),' %'])
disp(['Total Savings:                    ',num2str(round(TSS(3),3)),' £'])
disp(['Net Present Value:                ',num2str(round(NPVVV(3),3)), ' £'])
disp(['Annual Return on Investment:      ',num2str(round(aroii(3),3)),' %'])
disp(['SoH at the end of Lifetime:       ',num2str(round(SOHHY(3),3)),' %'])
disp(['Operation years to minimum SoH:   ',num2str(LTYNY(3)),' Years'])
disp(['Discounted Payback Period:        ',num2str(PBY(3)),' Years'])
disp(['Subsidy Required:                 ',num2str((SBY(3))), ' £/kWh'])
        end
disp('=============================================================================================')
disp(['+++++++++++++++++++++++++++++   PV+BESS Cost Benefit Analysis  +++++++++++++++++++++++++++++'])
disp('---------------------------------------------------------------------------------------------')
disp(['---  PV+BESS Cost Benefit Analysis for a Lifetime of ',num2str(Lifpv),' Years (PV) and ',num2str(lttyy),' Years of ',num2str(NuB),' BESS ---'])
disp('---------------------------------------------------------------------------------------------')
disp(['Total Savings:                    ',num2str(round((sum(TS11)+sum(TSS)),3)),' £'])
disp(['Net Present Value:                ',num2str(round(NPV11+sum(NPVVV),3)), ' £'])
AROIIIIA=100*(NPV11+sum(NPVVV))/(Lifpv*CapexT1);
disp(['Annual Return on Investment:      ',num2str(round(AROIIIIA,3)),' %'])
disp(['Residual years in the last BESS:  ',num2str(round(LTY-LYY(end))), ' Years'])
disp(['Project Operation Period:         ',num2str(Lifpv),' Years'])
disp('---------------------------------------------------------------------------------------------')
disp('=============================================================================================')
disp('      ')
end
toc;