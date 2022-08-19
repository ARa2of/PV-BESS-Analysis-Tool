%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function objec=MAINCODE(x)
global ERRa kPF IRR Tariff DataRes OPTTY PVSize BESS BESSP DOD SOCMAX SOCMIN SOCI RE PVCost InvCost InvSize Lifinv Lifpv PVdeg PVOM PRP EXP EX YearI IR er Dr PCN TLS TLE PTHC PTHD LTY SOHM BP SC FB DIAA
%% Load Profiles
Profile=readmatrix('Inputs.csv');
D=Profile(:,1); %Demand
PV=Profile(:,2);
EV=Profile(:,3);
T=length(Profile(:,1));
ND=round((DataRes/60)*(T/24)); %Number of days
TD=T/ND; %Length of one day
tau=TD/(24); % {Time interval=1/tau}

if OPTTY=="B"
BESS=round((x(1)),2); %Actual BESS Capacity [kWh]
BESSP=1.245+0.304*(BESS); %BESS Power rating [kW]
end
if OPTTY=="PV"
PVSize=round(x(1),1); %Size in kW
end
% 
PV=PVSize.*PV;

if Tariff=="FT"
TLE=1;
TLS=1;
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
L=0;
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

end
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
     PB11(hg)=inf;
 end
end

if NPV11<0
SB11=round(abs(NPV11)/PVSize,3);
else
    SB11=0;
end

%Battery
SOHY=[1 1-cumsum(L(end)*ones(1,(60)))];%%%%+++%%%%%

TS1=0;NPV1=0;AROI1=0;PB1=0;h=zeros;
for i=1:LTY
   TS1(i)=((BILL1-BILL2)*(ERRa(kPF+i))*SOHY(i))/(IRR(kPF+i));
end
NPV1=sum(TS1)-(((BP+FB)*BESS));
AROI1=100*NPV1/((((BP+FB)*BESS))*LTY);

if OPTTY=="B"
objec=NPV1*-1;
end
if OPTTY=="PV"
objec=NPV11*-1;
end


end
