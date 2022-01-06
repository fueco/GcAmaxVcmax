%%%%%Calculate the Gc by inverting the Penman-Monteith equation

clear all
str='.../fluxnetdata/';%set fluxnet data path
cd(str);
d=dir('*.csv');

%input flux data
i=1; %site id number
name1=d(i).name;
name=name1(1:length(name1)-4);

fLoad_ABACUS_csv

YEAR=floor(TIMESTAMP_START/1e8);
MONTH=floor((TIMESTAMP_START-YEAR*1e8)/1e6);
DAY=floor((TIMESTAMP_START-YEAR*1e8-MONTH*1e6)/1e4);
HOUR=floor((TIMESTAMP_START-YEAR*1e8-MONTH*1e6-DAY*1e4)/100)+rem((TIMESTAMP_START-YEAR*1e8-MONTH*1e6-DAY*1e4)/100,1)*5/3;
LEAP=zeros(length(YEAR),1);LEAP(rem(YEAR/4,1)==0)=1;
dMONTH=[0 31 28 31 30 31 30 31 31 30 31 30 31];
dMONTH_L=[0 31 29 31 30 31 30 31 31 30 31 30 31];
csMONTH=cumsum(dMONTH);
csMONTH_L=cumsum(dMONTH_L);
DOY=DAY+csMONTH(MONTH)';
CSLEAP=csMONTH_L(MONTH);
DOY(LEAP==1)=DAY(LEAP==1)+CSLEAP(LEAP==1)';

for kk=1:length(YEAR)
    if NEE_VUT_REF(kk)==-9999
        GPP_NT_VUT_REF(kk)=nan;
    end
    if NEE_VUT_REF_QC(kk)>=2
        GPP_NT_VUT_REF(kk)=nan;
    end
    if GPP_NT_VUT_REF(kk)<0
        GPP_NT_VUT_REF(kk)=nan;
    end
    if SWC_F_MDS_1(kk)<0
        SWC_F_MDS_1(kk)=nan;
    end
end

%Calculation of Gc
dat=[USTAR,WS_F,TA_F,PA_F,NETRAD,G_F_MDS,LE_F_MDS,VPD_F,SWC_F_MDS_1,SW_IN_F,GPP_NT_VUT_REF];

k=0.4;       %von K constant 
u8=dat(:,1); %ustar
ws=dat(:,2); %wind speed
ra1=log(exp(k*ws./u8)-0.7).*log(exp(k*ws./u8)-0.7);
ra=ra1./(ws*k*k);ga=1./ra; %[m s-1]
cp=1006;%[J kg-1 K-1]
ta=dat(:,3);
lv=2500-2.36*ta; %[kJ kg-1]
Pa=dat(:,4);     %Atmospheric pressure
psy=(cp*Pa)./(0.622*lv*1000);%Psychrometric constant [kPa K-1]
% D=(4098*(0.6108*exp((17.27*ta)./(ta+237.3)))./((ta+237.3).*(ta+237.3)));%[kPa c-1],
D1=(2508*(exp((17.27*ta)./(ta+237.3)))./((ta+237.3).*(ta+237.3)));%[kPa k-1],
%Rn-G wm-2
Rn=dat(:,5); %net radiation
G=dat(:,6);  %soil heat flux
Q=Rn-G;
%Density of dry air
da=(101.3*1000)./(287.04*(ta+273.15)); %[kg m-3]
lvE=dat(:,7);
vpd=0.1*dat(:,8); %kPa
fzi=ga.*psy;
fmu=(D1.*Q+da*cp.*ga.*vpd)./lvE-(D1+psy);
gc=fzi./fmu;  %m/s
gc(gc<0)=nan;
gc(gc>0.06)=nan;%m/s   %[0,2000 mmol m-2 s-1]
