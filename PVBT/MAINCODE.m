%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PCN=0; %Percentage of the BESS capacity to charge at night with low tariff (The ToU Tariff data is being used), set it to 0 if you don't want to use this option. 
PCN=PCN*ones(1,ND);
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
if Tariff=="FT"
TLE=1;
TLS=1;
end
% if Tariff=="TT"
% TLE=THS;
% end
if DIA==1
    DIAA=tau;
else
    DIAA=1;
end
%% Program 1: Conventional Rule-based Algorithm (CRBA)
Pnet=zeros(1,T); %Net 
BESSC=zeros(1,T); %Charge
BESSCC=zeros(1,T); %Charge
BESSD=zeros(1,T); %Discharge
BESSR=zeros(1,T); %Residual
BESSRR=zeros(1,ND);
SOC=ones(1,T)*SOCMIN; %SOC
SOC(1)=SOCI;
FT=zeros;FCAL=zeros;FDOD=zeros;FSOC=zeros;FC=zeros;FCYC=zeros;FD=zeros;L=zeros; %Degradation Model parameters
g=2;f=1;gf=0;k=1;
for i=1:T
   Pnet(i)=D(i)+EV(i)-PV(i); 
   %Charge
   if Tariff=="DT" || Tariff=="TT"
%    Charge at night low ToU  
      if i==((TLS)*tau+gf)
    BESSRR(k)=(SOCMAX-SOC(i))*BESS;
   end
   if (i>=((TLS*tau)+gf)&& i<=(TLE*tau+gf)) 
         BESSC(i)=(PCN(k)*BESSRR(k)/RE)/(TLE-TLS);
   end
   end
   %Charge in normal time according to a Threshold
     if Pnet(i)<PTHC && SOC(i)<SOCMAX
      BESSR(i)=(SOCMAX-SOC(i))*BESS;
     if (PTHC-Pnet(i))*(1/tau)<=BESSR(i) && (PTHC-Pnet(i))<=BESSP 
      BESSC(i)=(PTHC-Pnet(i));
      elseif (PTHC-Pnet(i))>=BESSP  && BESSP*(1/tau)<=BESSR(i)
      BESSC(i)=BESSP;
     else 
      BESSC(i)=BESSR(i)*(tau);
     end
     end
   %Discharge
    if Pnet(i)>PTHD && SOC(i)>SOCMIN &&  (i>=(TLE*DIAA+gf))
    BESSR(i)=(SOC(i)-SOCMIN)*BESS;
    if (Pnet(i)-PTHD)*(1/tau)<=BESSR(i) && (Pnet(i)-PTHD)<=BESSP
    BESSD(i)=(Pnet(i)-PTHD);
    elseif (Pnet(i)-PTHD)>=BESSP && BESSP*(1/tau)<=BESSR(i)
    BESSD(i)=BESSP; 
    else
    BESSD(i)=BESSR(i)*(tau);    
    end
    end
   SOC(g)=SOC(i)+(((BESSC(i)*(1/tau)*RE))/BESS)-(((BESSD(i)*(1/tau))/RE)/BESS);  
    
   if SOC(g)<SOCMIN
      SOC(g)=SOCMIN; 
   end
            if SOC(g)>SOCMAX
      SOC(g)=SOCMAX; 
   end
   %Update days
   if i==TD*k
      gf=TD*k;
      k=k+1;
   end 

g=g+1;
end

%Count Cycles (Rainflow Count Algorithm)
c=rainflow(SOC);
cyc=c(:,1);
S=c(:,4);
E=c(:,5);
DODD=c(:,2);
SOCD=c(:,3);
LL1=S(1:end-1);
LL2=E(1:end-1);
S1=SOC(LL1);
S2=SOC(LL2);
NCYC=round(abs(sum(S1-(S2)))/(DOD)); %Number of full cycles at DoD
NCYCPF=sum(cyc); %Number of all cycles partial and full from the rainflow algorithm 

% Degradation Model (Cycle Aging)
for f=1:length(S1)
    FT(f)=(exp((6.93*10^-2)*(20-25)*(25/20)));
    FC(f)=(exp(0.263*((DODD(f)/(abs(LL1(f)-LL2(f)))*1/tau)-1)));
    FDOD(f)=(((8.95*10^4)*(DODD(f)^(-0.486)))-7.28*10^4)^-1;
    FSOC(f)=(exp(1.04*(SOCD(f)-(0.5))));
    FCYC(f)=FDOD(f)*FSOC(f)*FC(f)*FT(f)*(abs(S1(f)-S2(f))/(DOD));
end
FCAL=(4.14*(10^-10))*ND*24*60*(exp(1.04*(mean(SOCD)-0.5)))*(exp((6.93*10^-2)*(20-25)*(25/20)));
FD=sum(FCYC)+FCAL;
L=1-((0.0575*exp(FD*-121))+((1-0.0575)*exp(-FD)));
L2=FCYC(1:length(S1));
L2 = cumsum(L2);
LFL=(L2/max(L2));
LFL2=LFL*(BESS*L);
BESSn=[BESS BESS-LFL2];  %BESS residual Capcity 
BESSPer=100*BESSn/BESS;
SOH=1-L;
%% Metrices Calculations
BESSNP=BESSC-BESSD;
%Curtailed Power
PNETC1=zeros(1,T);PNETC2=zeros(1,T);PCU1=zeros;PCU2=zeros;PSC1=zeros;PSC2=zeros;p=1;PDD1=zeros(1,T);PDD2=zeros(1,T);n=1;m=1;
for ip=1:T
 PDD1(ip)=D(ip)+EV(ip);
 PNETC1(ip)=D(ip)+EV(ip)-PV(ip)+BESSNP(ip);
 PNETC2(ip)=D(ip)+EV(ip)-PV(ip);
 if PNETC1(ip)<0 && PNETC1(ip)*-1>EXP
   PCU1(ip)= (PNETC1(ip)*-1)-EXP;
 end
 if PNETC2(ip)<0 && PNETC2(ip)*-1>EXP
   PCU2(ip)= (PNETC2(ip)*-1)-EXP;
 end
end
PCUA=sum(PCU1)*1/tau;%power curtailed after BESS
PCUB=sum(PCU2)*1/tau;%power curtailed after BESS


%% Electricity Bill Calculations 
totold=zeros(1,ND);
totnew=zeros(1,ND);
totoldb=zeros(1,ND);
P=zeros(1,ND);
for j=1:ND
Pnetbb=zeros(1,TD); %Net
Pneto=zeros(1,TD); %Net
PNETNEW=zeros(1,TD); %Net
g=1;
for k=(1+TD*(j-1)):(TD*(j))
    Pnetbb(g)=D(k)+EV(k); 
    Pneto(g)=D(k)+EV(k)-PV(k); 
    PNETNEW(g)=D(k)+EV(k)-PV(k)+BESSNP(k);
    g=g+1;
end
Pneto(Pneto<-EXP)=0;
PNETNEW(PNETNEW<-EXP)=0;
% w/o BESS
PNETDO=max(0,Pneto);
PNETPVO=min(0,Pneto);
PNETPVO(PNETPVO<-EXP)=0;
EBO=PNETDO.*PRP+PNETPVO*EX;
totold(j)=((sum(EBO)*(1/tau))/100)+(SC/100);
EEr(j)=sum(PNETPVO)*1/tau;

% w/ BESS
PNETD=max(0,PNETNEW);
PNETPV=min(0,PNETNEW);
PNETPV(PNETPV<-EXP)=0;
EBN=PNETD.*PRP+PNETPV*EX;
totnew(j)=((sum(EBN)*(1/tau))/100)+(SC/100);
EErr(j)=sum(PNETPV)*1/tau;

EBBB=Pnetbb.*PRP;
totoldb(j)=((sum(EBBB)*(1/tau))/100)+(SC/100); %without PV or BESS

%save in file
if SaveR==1
range10 = strcat('A',num2str(j+1));
range11 = strcat('B',num2str(j+1));
range9 = strcat('C',num2str(j+1));
range1 = strcat('D',num2str(j+1));

writematrix(j,'Results.xls','Sheet',2,'Range', range10)
writematrix(totoldb(j),'Results.xls','Sheet',2,'Range', range11)
writematrix(totold(j),'Results.xls','Sheet',2,'Range', range9)
writematrix(totnew(j),'Results.xls','Sheet',2,'Range', range1)
end
end
RESBESSC=sum(BESSC)*1/tau+(sum(EEr)-sum(EErr));
%Self-Consumption(SC) and Self-Sufficiency (SS)
SC1=(((sum(PV)*1/tau)+sum(EEr))/(sum(PV)*1/tau))*100; %SC before BESS
SS1=(((sum(PV)*1/tau)+sum(EEr))/(sum(PDD1)*1/tau))*100; %SS before BESS

SC2=((((sum(PV)*1/tau)+sum(EErr)))/(((sum(PV)*1/tau))))*100; %SC AFTER BESS
SS2=((((sum(PV)*1/tau)+sum(EErr))+(sum(BESSD)*1/tau))/((sum(PDD1)*1/tau)+(sum(BESSC)*1/tau)))*100; %SS after BESS

BILL1=sum(round(totold,2)); %Elec bill w/o BESS
BILL2=sum(round(totnew,2)); %Elec bill w BESS
BILL0=sum(round(totoldb,2));%Elec bill w/o PV and BESS
%% Cost Benefit Analysis 
%PV
DEG=[0 cumsum(PVdeg*ones(1,50-1))];
CAPEX1=(PVCost*PVSize)+(InvSize*InvCost*Lifpv/Lifinv);
TS11=0;NPV11=0;AROI11=0;
for i=1:Lifpv
   TS11(i)=(((BILL0-BILL1)*((1+er)^(i-1))*(1-DEG(i)))-((PVCost*PVSize)*PVOM))/((1+IR)^(i-1));
end
NPV11=sum(TS11)-CAPEX1;
AROI11=100*NPV11/(((CAPEX1))*Lifpv);

hg=1;TS111=0;NPV111=0;PB11=0;
for j=1:50
  for i=1:j
   TS111(i)=(((BILL0-BILL1)*((1+er)^(i-1))*(1-DEG(i)))-((PVCost*PVSize)*PVOM))/((1+IR)^(i-1));
end  
NPV111(j) =sum(TS111)-CAPEX1;
 if NPV111(j)>=0
     PB11(hg)=j;
    hg=hg+1; 
 else 
     PB11(hg)="N/A";
 end
end

if NPV11<0
SB11=round(abs(NPV11)/PVSize,3);
else
    SB11=0;
end

%Battery
SOHY=[1 1-cumsum(L(end)*ones(1,LTY-1))];
TS1=0;NPV1=0;AROI1=0;PB1=0;
for i=1:LTY
   TS1(i)=((BILL1-BILL2)*((1+er)^(i-1))*SOHY(i))/((1+IR)^(i-1));
end
NPV1=sum(TS1)-(((BP+FB)*BESS));
AROI1=100*NPV1/((((BP+FB)*BESS))*LTY);

%Determine the cost benfit analysis at 60% SoH
SOHYY=[1 1-cumsum(L(end)*ones(1,100-1))];
f=1;
for lk=1:length(SOHYY)
    if SOHYY(lk)<=SOHM
        h(f)=lk;
        f=f+1;
    end
end
LTYN=h(1); %lIFETIME AT 60% SOH
for i=1:LTYN
   TS2(i)=((BILL1-BILL2)*((1+er)^(i-1))*SOHYY(i))/((1+IR)^(i-1));
end
NPV2=sum(TS2)-(((BP+FB)*BESS));
AROI2=100*NPV2/((((BP+FB)*BESS))*LTYN);
%Determine valid Payback period
hg=1;PB3=zeros;
for j=1:50
  for i=1:j
   TS3(i)=((BILL1-BILL2)*((1+er)^(i-1))*SOHYY(i))/((1+IR)^(i-1));
end  
   NPV3(j)=sum(TS3)-(((BP+FB)*BESS));
 if NPV3(j)>=0
     PB3(hg)=j;
     hg=hg+1;
 end
end
if sum(PB3)==0
PB3(1)=inf;
end

if PB3(1)>LTYN
    PB2='N/A';
end
if PB3(1)>LTY
    PB1='N/A';
end
if PB3(1)<=LTY
    PB1=PB3(1);
end
if PB3(1)<=LTYN
    PB2=PB3(1);
end

if NPV1<0
SB1=round(abs(NPV1)/BESS,3);
else
    SB1=0;
end
if NPV2<0
SB2=round(abs(NPV2)/BESS,3);
else
    SB2=0;
end

% PV+BESS For BESS Lifetime
NuB=ceil(Lifpv/LTY);
for i=1:(Lifpv*3)
   ERRa(i)=(1+er)^(i-1); 
   IRR(i)=((1+IR)^(i-1));
end

if Lifpv/LTY <NuB
    Res=NuB-Lifpv/LTY;
elseif Lifpv/LTY == NuB
    Res=0;
end
TSM=ones(1,NuB);
TSM(end)=1-Res;
LYY=TSM*LTY;
% EFGO=cumsum(LYY);
% GF=[1 EFGO(1:end-1)];
GF=zeros(1,NuB);


BPF=zeros;
for i=1:LTY*(NuB-1)
   BPF(i)=((1-(Dr))^i)*BP;
end
BPF=[BP BPF(1:end-1)]+FB;

jjp=0;
go=1;
for oi=1:NuB
BPRi(oi)= BPF(go+jjp*(LTY)); 
go=0;
jjp=jjp+1;
end

for lp=1:NuB
    TS1A=0;
    kPF(lp)=GF(lp);
for i=1:LYY(lp)
   TS1A(i)=((BILL1-BILL2)*(ERRa(kPF(lp)+i))*SOHYY(i))/(IRR(kPF(lp)+i));
end
TST(lp)=sum(TS1A);
end

% TOTS1=sum(TS11)+((Lifpv/LTY)*sum(TS1));
TOTS1=sum(TS11)+sum(TST);
CapexT1=CAPEX1+sum(BPRi*BESS);
NPVT1=TOTS1-CapexT1;
AROIT1=100*NPVT1/(((CapexT1))*Lifpv);

if NPVT1<0
SB4=round(abs(NPVT1)/BESS,3);
else
    SB4=0;
end

for TT=1:ceil(Lifpv/LTY)
   YearIn(TT)= YearI+TT*LTY;
end
YearIns=[YearI YearIn(1:end-1)];

TS1B=repmat(TS1,ceil(Lifpv/LTY));
TS1B = TS1B';
TS1B = TS1B(:)';
hg=1;
for j=1:Lifpv
for b=1:j
   TOTB1(b)=TS11(b)+TS1B(1);
end
NPVT1PB(j)=sum(TOTB1)-CapexT1;

 if NPVT1PB(j)>=0
     PB4(hg)=j;
    hg=hg+1; 
 else 
     PB4(hg)="N/A";
 end
end

% PV+BESS For BESS SoH 
NuB2=ceil(Lifpv/LTYN);
BPF2=zeros;
for i=1:LTYN*(NuB2-1)
   BPF2(i)=((1-(Dr))^i)*BP;
end
BPF2=[BP BPF2(1:end-1)]+FB;

jjp=0;
go=1;
for oi=1:NuB2
BPRii(oi)= BPF2(go+jjp*(LTYN)); 
go=0;
jjp=jjp+1;
end

if Lifpv/LTYN <NuB2
    resp=Lifpv-LTYN;
    if resp>=LTYN
            Res2=(Lifpv/LTYN)-floor(Lifpv/LTYN);
    end
    if resp<LTYN
        Res2=resp/LTYN;
    end
end
if Lifpv/LTYN == NuB2
    Res2=1;
end

TSM2=ones(1,NuB2);
TSM2(end)=Res2;
LYY2=TSM2*LTYN;
EFGO2=cumsum(LYY2);
GF2=[1 EFGO2(1:end-1)];


for lp=1:NuB2
    TS1A=0;
    kPF2(lp)=GF2(lp);
for i=1:LYY2(lp)
   TS1A(i)=((BILL1-BILL2)*(ERRa(kPF2(lp)+i))*SOHYY(i))/(IRR(kPF2(lp)+i));
end
TST2(lp)=sum(TS1A);
end

TOTS12=sum(TS11)+(sum(TST2));
CapexT12=CAPEX1+sum(BPRii*BESS);
NPVT12=TOTS12-CapexT12;
AROIT12=100*NPVT12/(((CapexT12))*Lifpv);
if NPVT12<0
SB5=round(abs(NPVT12)/BESS,3);
else
    SB5=0;
end

for TT=1:ceil(Lifpv/LTYN)
   YearIn2(TT)= YearI+TT*LTYN;
end
YearIns2=[YearI YearIn2(1:end-1)];

TS2B=repmat(TS2,ceil(Lifpv/LTYN));
TS2B = TS2B';
TS2B = TS2B(:)';
hg=1;
for j=1:Lifpv
for b=1:j
   TOTB2(b)=TS11(b)+TS2B(1);
end
NPVT2PB(j)=sum(TOTB2)-CapexT1;

 if NPVT2PB(j)>=0
     PB5(hg)=j;
    hg=hg+1; 
 else 
     PB5(hg)="N/A";
 end
end


%% Figures and Results
Pnet=D+EV-PV;
PNETW=D+EV-PV+BESSNP';

if SaveR==1
RES=[(1:T)' Pnet PNETW BESSNP' SOC(1:T)'];
range101 = strcat('A',num2str(2));
writematrix(RES,'Results.xls','Sheet',1,'Range', range101)
range101A = strcat('A',num2str(2));
writematrix(BESSPer(end),'Results.xls','Sheet',3,'Range', range101A)
range101B = strcat('B',num2str(2));
writematrix(BESSn(end),'Results.xls','Sheet',3,'Range', range101B)
end

figure
subplot (1,2,1)
yyaxis left
ylabel('Power [kW]')
hold on
plot(PNETW,'r-')
hold on
plot(Pnet,'k')
yyaxis right
ylabel('SOC [%]')
stairs(SOC,'b--')
legend({'w/ BESS','w/o BESS','BESS SOC'},'Location','northeast','AutoUpdate','off')
axis([1 T -inf inf])
set(gca,'FontSize',14,'FontName', 'Times New Roman')
xlabel('Time')

subplot (1,2,2)
newcolors = [0, 0.4470, 0.7410
0.8500, 0.3250, 0.0980];  
colororder(newcolors)
x=linspace(0,ND,(length(BESSn)));
yyaxis left
ylabel('BESS Capacity [kWh]')
hold on
Y1=plot(x,BESSn,'-');
hold on
yyaxis right
ylabel('SoH [%]')
Y2=plot(x,BESSPer,'--');
legend({'BESS Capacity', 'BESS SoH'},'Location','northeast','AutoUpdate','off')
xlim([0 ND])
set(gca,'FontSize',14,'FontName', 'Times New Roman')
xlabel('Days')
grid on;

LTY
%% Print
disp('      ')
disp('=============================================================================================')
disp('********************                  Building  Data                   **********************')
disp('=============================================================================================')
disp(['Total Annual Consumption:              ',num2str(round((sum(D)+sum(EV))*1/tau,3)),' kWh/year'])
disp(['Total PV Annual Generation:            ',num2str(round(sum(PV)*1/tau,3)),' kWh/year'])
disp(['Tariff Structure:                      ',num2str(Tariff)])
disp('=============================================================================================')
disp('********************           SIMULATION RESULTS OF THE PVBT          **********************')
disp('=============================================================================================')
disp('-------------------------------------- One-Year Results -------------------------------------')
disp('---------------------------------------------------------------------------------------------')
disp(['Electricity Bill w/o PV w/o BESS:       ',num2str(round(BILL0,3)),' £'])
disp(['Electricity Bill w/ PV w/o BESS:        ',num2str(round(BILL1,3)),' £'])
disp(['Electricity Bill w/ PV w/ BESS:         ',num2str(round(BILL2,3)),' £'])
disp(['Curtailed PV power w/o BESS:            ',num2str(round(PCUB,3)),' kWh'])
disp(['Curtailed PV power w/ BESS:             ',num2str(round(PCUA,3)),' kWh'])
disp(['Exported PV power w/o BESS:             ',num2str(round(sum(EEr*-1),3)),' kWh'])
disp(['Exported PV power w/ BESS:              ',num2str(round(sum(EErr*-1),3)),' kWh'])
disp(['PV Self-Consumption w/o BESS:           ',num2str(round(SC1,3)),' %'])
disp(['PV Self-Consumption w/ BESS:            ',num2str(round(SC2,3)),' %'])
disp(['Premises Self-Sufficiency w/o BESS:     ',num2str(round(SS1,3)),' %'])
disp(['Premises Self-Sufficiency w/ BESS:      ',num2str(round(SS2,3)),' %'])
disp(['BESS State of Health (SoH):             ',num2str(round(SOH*100,3)),' %'])
disp('=============================================================================================')
disp(['+++++++++++++++++++  PV Cost Benefit Analysis for a Lifetime of ',num2str(Lifpv),' Years  ++++++++++++++++++++'])
disp('---------------------------------------------------------------------------------------------')
disp(['Annual Saving:                    ',num2str(round(sum(BILL0-BILL1),3)),' £'])
disp(['Total Savings:                    ',num2str(round(sum(TS11),3)),' £'])
disp(['Net Present Value:                ',num2str(round(NPV11,3)), ' £'])
disp(['Annual Return on Investment:      ',num2str(round(AROI11,3)),' %'])
disp(['Discounted Payback Period:        ',num2str(PB11(1)),' Years'])
disp(['Subsidy Required:                 ',num2str((SB11)), ' £/kW'])
disp('=============================================================================================')
disp(['+++++++++++++++++++++++++++++++   BESS Cost Benefit Analysis  +++++++++++++++++++++++++++++++'])
disp('---------------------------------------------------------------------------------------------')
disp(['Annual Savings:                   ',num2str(round(sum(BILL1-BILL2),3)),' £/year'])
disp('---------------------------------------------------------------------------------------------')
disp(['--------------------  Cost Benefit Analysis for a Lifetime of ',num2str(LTY),' Years  ---------------------'])
disp('---------------------------------------------------------------------------------------------')
disp(['Total Savings:                    ',num2str(round(sum(TS1),3)),' £'])
disp(['Net Present Value:                ',num2str(round(NPV1,3)), ' £'])
disp(['Annual Return on Investment:      ',num2str(round(AROI1,3)),' %'])
disp(['SoH at the end of Lifetime:       ',num2str(round(SOHY(LTY)*100,3)),' %'])
disp(['Discounted Payback Period:        ',num2str(PB1),' Years'])
disp(['Subsidy Required:                 ',num2str((SB1)), ' £/kWh'])
disp('---------------------------------------------------------------------------------------------')
disp(['--------------------  Cost Benefit Analysis for a Lifetime of ',num2str(SOHM*100),'% SoH  ----------------------'])
disp('---------------------------------------------------------------------------------------------')
disp(['Total Savings:                    ',num2str(round(sum(TS2),3)),' £'])
disp(['Net Present Value:                ',num2str(round(NPV2,3)), ' £'])
disp(['Annual Return on Investment:      ',num2str(round(AROI2,3)),' %'])
disp(['Years of Operation:               ',num2str(LTYN),' Years'])
disp(['Discounted Payback Period:        ',num2str(PB2),' Years'])
disp(['Subsidy Required:                 ',num2str((SB2)), ' £/kWh'])
disp('=============================================================================================')
disp(['+++++++++++++++++++++++++++++   PV+BESS Cost Benefit Analysis  +++++++++++++++++++++++++++++'])
disp('---------------------------------------------------------------------------------------------')
disp(['---  PV+BESS Cost Benefit Analysis for a Lifetime of ',num2str(Lifpv),' Years (PV) and ',num2str(LTY),' Years of ',num2str(NuB),' BESS ---'])
disp('---------------------------------------------------------------------------------------------')
disp(['Total Savings:                    ',num2str(round(sum(TOTS1),3)),' £'])
disp(['Net Present Value:                ',num2str(round(NPVT1,3)), ' £'])
disp(['Annual Return on Investment:      ',num2str(round(AROIT1,3)),' %'])
disp(['Number of BESS:                   ',num2str(round(NuB))])
disp(['Years of installations:           ',num2str(YearIns)])
disp(['BESS price for each year:         ',num2str(round(BPRi)), ' £/kWh'])
disp(['Residual years in the last BESS:  ',num2str(round((LTY*NuB)-Lifpv)), ' Years'])
disp(['Project Operation Period:         ',num2str(Lifpv),' Years'])
disp(['Discounted Payback Period:        ',num2str(PB4(1)),' Years'])
disp(['Subsidy Required:                 ',num2str((SB4)), ' £/kWh'])
disp('---------------------------------------------------------------------------------------------')
disp(['---  PV+BESS Cost Benefit Analysis for a Lifetime of ',num2str(Lifpv),' Years (PV) and ',num2str(SOHM*100),'%  SoH of ',num2str(NuB2),' BESS ---'])
disp('---------------------------------------------------------------------------------------------')
disp(['Total Savings:                    ',num2str(round(sum(TOTS12),3)),' £'])
disp(['Net Present Value:                ',num2str(round(NPVT12,3)), ' £'])
disp(['Annual Return on Investment:      ',num2str(round(AROIT12,3)),' %'])
disp(['Number of BESS:                   ',num2str(round(NuB2))])
disp(['Years of installations:           ',num2str(YearIns2)])
disp(['BESS price for each year:         ',num2str(round(BPRii)), ' £/kWh'])
disp(['Residual years in the last BESS:  ',num2str(round((LTYN*NuB2)-Lifpv)), ' Years'])
disp(['Project Operation Period:         ',num2str(Lifpv),' Years'])
disp(['Discounted Payback Period:        ',num2str(PB5(1)),' Years'])
disp(['Subsidy Required:                 ',num2str((SB5)), ' £/kWh'])
disp('---------------------------------------------------------------------------------------------')
disp('=============================================================================================')
disp('      ')


if SaveR==1
R1 = strcat('B',num2str(1));R2 = strcat('B',num2str(2));R3 = strcat('B',num2str(3));R4 = strcat('B',num2str(4));R5 = strcat('B',num2str(5));R6 = strcat('B',num2str(6));
R7 = strcat('B',num2str(7));R8 = strcat('B',num2str(8));R9 = strcat('B',num2str(9));R10 = strcat('B',num2str(10));R11 = strcat('B',num2str(11));R12 = strcat('B',num2str(12));
R13 = strcat('B',num2str(13));R15 = strcat('B',num2str(15));R17 = strcat('B',num2str(17));R19 = strcat('B',num2str(19));R21 = strcat('B',num2str(21));R23 = strcat('B',num2str(23));
R14 = strcat('B',num2str(14));R16 = strcat('B',num2str(16));R18 = strcat('B',num2str(18));R20 = strcat('B',num2str(20));R22 = strcat('B',num2str(22));R24 = strcat('B',num2str(24));
R25 = strcat('B',num2str(25));R27 = strcat('B',num2str(27));R29 = strcat('B',num2str(29));R31 = strcat('B',num2str(31));R33 = strcat('B',num2str(33));R35 = strcat('B',num2str(35));
R26 = strcat('B',num2str(26));R28 = strcat('B',num2str(28));R30 = strcat('B',num2str(30));R32 = strcat('B',num2str(32));R34 = strcat('B',num2str(34));R36 = strcat('B',num2str(36));
R37 = strcat('B',num2str(37));R39 = strcat('B',num2str(39));R41 = strcat('B',num2str(41));R43 = strcat('B',num2str(43));R45 = strcat('B',num2str(45));R47 = strcat('B',num2str(47));
R38 = strcat('B',num2str(38));R40 = strcat('B',num2str(40));R42 = strcat('B',num2str(42));R44 = strcat('B',num2str(44));R46 = strcat('B',num2str(46));R48 = strcat('B',num2str(48));
R50 = strcat('B',num2str(50));R51 = strcat('B',num2str(51));R53 = strcat('B',num2str(53));R55 = strcat('B',num2str(55));R57 = strcat('B',num2str(57));R59 = strcat('B',num2str(59));
R49 = strcat('B',num2str(49));R52 = strcat('B',num2str(52));R54 = strcat('B',num2str(54));R56 = strcat('B',num2str(56));R58 = strcat('B',num2str(58));R60 = strcat('B',num2str(60));R61 = strcat('B',num2str(61));
R63 = strcat('B',num2str(63));

writematrix((round((sum(D)+sum(EV))*1/tau,3)),'Results.xls','Sheet',4,'Range', R1)
writematrix((round(sum(PV)*1/tau,3)),'Results.xls','Sheet',4,'Range', R2)
writematrix((round(BILL0,3)),'Results.xls','Sheet',4,'Range', R4)
writematrix((round(BILL1,3)),'Results.xls','Sheet',4,'Range', R5)
writematrix((round(BILL2,3)),'Results.xls','Sheet',4,'Range', R6)
writematrix((round(PCUB,3)),'Results.xls','Sheet',4,'Range', R7)
writematrix((round(PCUA,3)),'Results.xls','Sheet',4,'Range', R8)
writematrix((round(sum(EEr*-1),3)),'Results.xls','Sheet',4,'Range', R9)
writematrix((round(sum(EErr*-1),3)),'Results.xls','Sheet',4,'Range', R10)
writematrix((round(SC1,3)),'Results.xls','Sheet',4,'Range', R11)
writematrix((round(SC2,3)),'Results.xls','Sheet',4,'Range', R12)
writematrix((round(SS1,3)),'Results.xls','Sheet',4,'Range', R13)
writematrix((round(SS2,3)),'Results.xls','Sheet',4,'Range', R14)
writematrix((round(SOH*100,3)),'Results.xls','Sheet',4,'Range', R15)
writematrix((sum(BILL0-BILL1)),'Results.xls','Sheet',4,'Range', R17)
writematrix((round(sum(TS11),3)),'Results.xls','Sheet',4,'Range', R18)
writematrix((round(NPV11,3)),'Results.xls','Sheet',4,'Range', R19)
writematrix((round(AROI11,3)),'Results.xls','Sheet',4,'Range', R20)
writematrix((PB11(1)),'Results.xls','Sheet',4,'Range', R21)
writematrix(((SB11)),'Results.xls','Sheet',4,'Range', R22)
writematrix((round(sum(BILL1-BILL2),3)),'Results.xls','Sheet',4,'Range', R24)
writematrix((round(sum(TS1),3)),'Results.xls','Sheet',4,'Range', R26)
writematrix((round(NPV1,3)),'Results.xls','Sheet',4,'Range', R27)
writematrix((round(AROI1,3)),'Results.xls','Sheet',4,'Range', R28)
writematrix((round(SOHY(LTY)*100,3)),'Results.xls','Sheet',4,'Range', R29)
writematrix((PB1),'Results.xls','Sheet',4,'Range', R30)
writematrix((SB1),'Results.xls','Sheet',4,'Range', R31)
writematrix((round(sum(TS2),3)),'Results.xls','Sheet',4,'Range', R33)
writematrix((round(NPV2,3)),'Results.xls','Sheet',4,'Range', R34)
writematrix((round(AROI2,3)),'Results.xls','Sheet',4,'Range', R35)
writematrix(LTYN,'Results.xls','Sheet',4,'Range', R36)
writematrix(PB2,'Results.xls','Sheet',4,'Range', R37)
writematrix(SB2,'Results.xls','Sheet',4,'Range', R38)
writematrix(round(sum(TOTS1),3),'Results.xls','Sheet',4,'Range', R41)
writematrix((round(NPVT1,3)),'Results.xls','Sheet',4,'Range', R42)
writematrix(round(AROIT1,3),'Results.xls','Sheet',4,'Range', R43)
writematrix(round(NuB),'Results.xls','Sheet',4,'Range', R44)
writematrix(YearIns,'Results.xls','Sheet',4,'Range', R45)
writematrix(round(BPRi),'Results.xls','Sheet',4,'Range', R46)
writematrix(round((LTY*NuB)-Lifpv),'Results.xls','Sheet',4,'Range', R47)
writematrix(Lifpv,'Results.xls','Sheet',4,'Range', R48)
writematrix(PB4(1),'Results.xls','Sheet',4,'Range', R49)
writematrix((SB4),'Results.xls','Sheet',4,'Range', R50)
writematrix(round(sum(TOTS12),3),'Results.xls','Sheet',4,'Range', R52)
writematrix(round(NPVT12,3),'Results.xls','Sheet',4,'Range', R53)
writematrix(round(AROIT12,3),'Results.xls','Sheet',4,'Range', R54)
writematrix(round(NuB2),'Results.xls','Sheet',4,'Range', R55)
writematrix(YearIns2,'Results.xls','Sheet',4,'Range', R56)
writematrix(round(BPRii),'Results.xls','Sheet',4,'Range', R57)
writematrix(round((LTYN*NuB2)-Lifpv),'Results.xls','Sheet',4,'Range', R58)
writematrix(Lifpv,'Results.xls','Sheet',4,'Range', R59)
writematrix(PB5(1),'Results.xls','Sheet',4,'Range', R60)
writematrix(SB4,'Results.xls','Sheet',4,'Range', R61)
writematrix((Tariff),'Results.xls','Sheet',4,'Range', R63)

end



toc;

