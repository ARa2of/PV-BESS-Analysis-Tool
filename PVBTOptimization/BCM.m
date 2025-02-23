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
    if Pnet(i)>PTHD && SOC(i)>SOCMIN && (i>=(TLE*DIAA+gf))
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
        if k > length(PCN)
            k = length(PCN); % Ensure k stays within valid bounds
        end
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
PCUB=sum(PCU2)*1/tau;%power curtailed before BESS

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
totold(j)=(sum(EBO)*(1/tau))/100;
EEr(j)=sum(PNETPVO)*1/tau;

% w/ BESS
PNETD=max(0,PNETNEW);
PNETPV=min(0,PNETNEW);
PNETPV(PNETPV<-EXP)=0;
EBN=PNETD.*PRP+PNETPV*EX;
totnew(j)=(sum(EBN)*(1/tau))/100;
EErr(j)=sum(PNETPV)*1/tau;

EBBB=Pnetbb.*PRP;
totoldb(j)=(sum(EBBB)*(1/tau))/100; %without PV or BESS
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
k=k+1;
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
SOHY=[1 1-cumsum(L(end)*ones(1,(60)))];%%%%+++%%%%%
TS1=0;NPV1=0;AROI1=0;PB1=0;
for i=1:LTY
   TS1(i)=((BILL1-BILL2)*(ERRa(kPF+i))*SOHY(i))/(IRR(kPF+i));
end
NPV1=sum(TS1)-(((BP+FB)*BESS));
AROI1=100*NPV1/((((BP+FB)*BESS))*LTY);

%Determine the cost benfit analysis at 60% SoH
SOHYY=[1 1-cumsum(L(end)*ones(1,60))];
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
hg=1;PB3=zeros;TS3=zeros;NPV3=zeros;
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
    PB1=inf;
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

NuB=ceil(Lifpv/lttyy);
BPF=zeros;
for i=1:lttyy*(NuB-1)
   BPF(i)=((1-(Dr))^i)*BP;
end
BPF=[BP BPF(1:end-1)]+FB;

jjp=0;
go=1;
for oi=1:NuB
BPRi(oi)= BPF(go+jjp*(lttyy)); 
go=0;
jjp=jjp+1;
end

TOTS1=sum(TS11)+((Lifpv/lttyy)*sum(TS1));
CapexT1=CAPEX1+sum(BPRi*BESS);
NPVT1=TOTS1-CapexT1;


AROIT1=100*NPVT1/(((CapexT1))*Lifpv);

if NPVT1<0
SB4=round(abs(NPVT1)/BESS,3);
else
    SB4=0;
end

TS1B=repmat(TS1,ceil(Lifpv/lttyy));
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
TOTS12=sum(TS11)+((Lifpv/LTYN)*sum(TS2));
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