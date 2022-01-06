%%%%%Calculate the Vcmax

clear all
str='.../fluxnetdata/';%set flux data path
cd(str);
d=dir('*.csv');
data(1,1:13)=nan;

%input flux data
i=1;%site id number
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

str='.../fullsetflux2015/data_REddyProcOutput/';%This is the path of Amax data
cd(str);dd=dir('*.csv');
name1=dd(i).name;
A = readtable(name1);beta2=A.beta;k2=A.k;
c=find(A.beta>0);
for j=1:length(c)-1
    if c(j+1)-c(j)<673
        beta2(c(j)+1:c(j+1)-1,1)=A.beta(c(j));
        k2(c(j)+1:c(j+1)-1,1)=A.k(c(j));
    else
        beta2(c(j)+1:c(j)-1+672,1)=A.beta(c(j));
        k2(c(j)+1:c(j)-1+672,1)=A.k(c(j));
    end
end
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
    if TA_F(kk)<-100
        TA_F(kk)=nan;
    end
end
beta=beta2;
for kk=1:length(YEAR)
    if VPD_F(kk)>10
        beta(kk)=beta2(kk).*exp(-k2(kk).*(VPD_F(kk)-10));
    end
    if beta2(kk)<0
        beta(kk)=nan;
    end
    if VPD_F(kk)<-100
        beta(kk)=nan;
        VPD_F(kk)=nan;
    end
end

for year=min(YEAR):max(YEAR)
    for k=1:max(DOY(YEAR==year))
        data0(k,1)=nanmean(USTAR(beta>0&NETRAD>0&CO2_F_MDS_QC<2&CO2_F_MDS>0&GPP_NT_VUT_REF>=0&NEE_VUT_REF_QC<2&VPD_F>-1000&TA_F>-100&TA_F_QC<2&VPD_F_QC<2&SW_IN_F_QC<2&SW_IN_F>=500&SWC_F_MDS_1>=0&SWC_F_MDS_1_QC<2&LE_F_MDS>=0&G_F_MDS>=-1000&LE_F_MDS_QC<2&G_F_MDS_QC<2&USTAR>=0.2&WS_F>=0&P_F==0&HOUR<=14&HOUR>=11&DOY==k&YEAR==year));
        data0(k,2)=nanmean(WS_F(beta>0&NETRAD>0&CO2_F_MDS_QC<2&CO2_F_MDS>0&GPP_NT_VUT_REF>=0&NEE_VUT_REF_QC<2&VPD_F>-1000&TA_F>-100&TA_F_QC<2&VPD_F_QC<2&SW_IN_F_QC<2&SW_IN_F>=500&SWC_F_MDS_1>=0&SWC_F_MDS_1_QC<2&LE_F_MDS>=0&G_F_MDS>=-1000&LE_F_MDS_QC<2&G_F_MDS_QC<2&USTAR>=0.2&WS_F>=0&P_F==0&HOUR<=14&HOUR>=11&DOY==k&YEAR==year));
        data0(k,3)=nanmean(TA_F(beta>0&NETRAD>0&CO2_F_MDS_QC<2&CO2_F_MDS>0&GPP_NT_VUT_REF>=0&NEE_VUT_REF_QC<2&VPD_F>-1000&TA_F>-100&TA_F_QC<2&VPD_F_QC<2&SW_IN_F_QC<2&SW_IN_F>=500&SWC_F_MDS_1>=0&SWC_F_MDS_1_QC<2&LE_F_MDS>=0&G_F_MDS>=-1000&LE_F_MDS_QC<2&G_F_MDS_QC<2&USTAR>=0.2&WS_F>=0&P_F==0&HOUR<=14&HOUR>=11&DOY==k&YEAR==year));
        data0(k,4)=nanmean(PA_F(beta>0&NETRAD>0&CO2_F_MDS_QC<2&CO2_F_MDS>0&GPP_NT_VUT_REF>=0&NEE_VUT_REF_QC<2&VPD_F>-1000&TA_F>-100&TA_F_QC<2&VPD_F_QC<2&SW_IN_F_QC<2&SW_IN_F>=500&SWC_F_MDS_1>=0&SWC_F_MDS_1_QC<2&LE_F_MDS>=0&G_F_MDS>=-1000&LE_F_MDS_QC<2&G_F_MDS_QC<2&USTAR>=0.2&WS_F>=0&P_F==0&HOUR<=14&HOUR>=11&DOY==k&YEAR==year));
        data0(k,5)=nanmean(NETRAD(beta>0&NETRAD>0&CO2_F_MDS_QC<2&CO2_F_MDS>0&GPP_NT_VUT_REF>=0&NEE_VUT_REF_QC<2&VPD_F>-1000&TA_F>-100&TA_F_QC<2&VPD_F_QC<2&SW_IN_F_QC<2&SW_IN_F>=500&SWC_F_MDS_1>=0&SWC_F_MDS_1_QC<2&LE_F_MDS>=0&G_F_MDS>=-1000&LE_F_MDS_QC<2&G_F_MDS_QC<2&USTAR>=0.2&WS_F>=0&P_F==0&HOUR<=14&HOUR>=11&DOY==k&YEAR==year));
        data0(k,6)=nanmean(G_F_MDS(beta>0&NETRAD>0&CO2_F_MDS_QC<2&CO2_F_MDS>0&GPP_NT_VUT_REF>=0&NEE_VUT_REF_QC<2&VPD_F>-1000&TA_F>-100&TA_F_QC<2&VPD_F_QC<2&SW_IN_F_QC<2&SW_IN_F>=500&SWC_F_MDS_1>=0&SWC_F_MDS_1_QC<2&LE_F_MDS>=0&G_F_MDS>=-1000&LE_F_MDS_QC<2&G_F_MDS_QC<2&USTAR>=0.2&WS_F>=0&P_F==0&HOUR<=14&HOUR>=11&DOY==k&YEAR==year));
        data0(k,7)=nanmean(LE_F_MDS(beta>0&NETRAD>0&CO2_F_MDS_QC<2&CO2_F_MDS>0&GPP_NT_VUT_REF>=0&NEE_VUT_REF_QC<2&VPD_F>-1000&TA_F>-100&TA_F_QC<2&VPD_F_QC<2&SW_IN_F_QC<2&SW_IN_F>=500&SWC_F_MDS_1>=0&SWC_F_MDS_1_QC<2&LE_F_MDS>=0&G_F_MDS>=-1000&LE_F_MDS_QC<2&G_F_MDS_QC<2&USTAR>=0.2&WS_F>=0&P_F==0&HOUR<=14&HOUR>=11&DOY==k&YEAR==year));
        data0(k,8)=nanmean(VPD_F(beta>0&NETRAD>0&CO2_F_MDS_QC<2&CO2_F_MDS>0&GPP_NT_VUT_REF>=0&NEE_VUT_REF_QC<2&VPD_F>-1000&TA_F>-100&TA_F_QC<2&VPD_F_QC<2&SW_IN_F_QC<2&SW_IN_F>=500&SWC_F_MDS_1>=0&SWC_F_MDS_1_QC<2&LE_F_MDS>=0&G_F_MDS>=-1000&LE_F_MDS_QC<2&G_F_MDS_QC<2&USTAR>=0.2&WS_F>=0&P_F==0&HOUR<=14&HOUR>=11&DOY==k&YEAR==year));
        data0(k,9)=nanmean(SWC_F_MDS_1(beta>0&NETRAD>0&CO2_F_MDS_QC<2&CO2_F_MDS>0&GPP_NT_VUT_REF>=0&NEE_VUT_REF_QC<2&VPD_F>-1000&TA_F>-100&TA_F_QC<2&VPD_F_QC<2&SW_IN_F_QC<2&SW_IN_F>=500&SWC_F_MDS_1>=0&SWC_F_MDS_1_QC<2&LE_F_MDS>=0&G_F_MDS>=-1000&LE_F_MDS_QC<2&G_F_MDS_QC<2&USTAR>=0.2&WS_F>=0&P_F==0&HOUR<=14&HOUR>=11&DOY==k&YEAR==year));
        data0(k,10)=nanmean(SW_IN_F(beta>0&NETRAD>0&CO2_F_MDS_QC<2&CO2_F_MDS>0&GPP_NT_VUT_REF>=0&NEE_VUT_REF_QC<2&VPD_F>-1000&TA_F>-100&TA_F_QC<2&VPD_F_QC<2&SW_IN_F_QC<2&SW_IN_F>=500&SWC_F_MDS_1>=0&SWC_F_MDS_1_QC<2&LE_F_MDS>=0&G_F_MDS>=-1000&LE_F_MDS_QC<2&G_F_MDS_QC<2&USTAR>=0.2&WS_F>=0&P_F==0&HOUR<=14&HOUR>=11&DOY==k&YEAR==year));
        data0(k,11)=nanmean(beta(beta>0&NETRAD>0&CO2_F_MDS_QC<2&CO2_F_MDS>0&GPP_NT_VUT_REF>=0&NEE_VUT_REF_QC<2&VPD_F>-1000&TA_F>-100&TA_F_QC<2&VPD_F_QC<2&SW_IN_F_QC<2&SW_IN_F>=500&SWC_F_MDS_1>=0&SWC_F_MDS_1_QC<2&LE_F_MDS>=0&G_F_MDS>=-1000&LE_F_MDS_QC<2&G_F_MDS_QC<2&USTAR>=0.2&WS_F>=0&P_F==0&HOUR<=14&HOUR>=11&DOY==k&YEAR==year));
        data0(k,12)=nanmean(CO2_F_MDS(beta>0&NETRAD>0&CO2_F_MDS_QC<2&CO2_F_MDS>0&GPP_NT_VUT_REF>=0&NEE_VUT_REF_QC<2&VPD_F>-1000&TA_F>-100&TA_F_QC<2&VPD_F_QC<2&SW_IN_F_QC<2&SW_IN_F>=500&SWC_F_MDS_1>=0&SWC_F_MDS_1_QC<2&LE_F_MDS>=0&G_F_MDS>=-1000&LE_F_MDS_QC<2&G_F_MDS_QC<2&USTAR>=0.2&WS_F>=0&P_F==0&HOUR<=14&HOUR>=11&DOY==k&YEAR==year));
        data0(k,13)=nanmean(beta(beta>0&NETRAD>0&CO2_F_MDS_QC<2&CO2_F_MDS>0&GPP_NT_VUT_REF>=0&NEE_VUT_REF_QC<2&VPD_F>-1000&TA_F>-100&TA_F_QC<2&VPD_F_QC<2&SW_IN_F_QC<2&SW_IN_F>=500&SWC_F_MDS_1>=0&SWC_F_MDS_1_QC<2&LE_F_MDS>=0&G_F_MDS>=-1000&LE_F_MDS_QC<2&G_F_MDS_QC<2&USTAR>=0.2&WS_F>=0&P_F==0&HOUR<=14&HOUR>=11&DOY==k&YEAR==year));
    end
    data=[data;data0];
    clear data0
end
data= rmmissing(data,1);
dat=data(data(:,10)>250&data(:,3)>15&data(:,8)>5,:);
k=0.4; %von K constant 
u8=dat(:,1); %ustar
ws=dat(:,2); %wind speed
ra1=log(exp(k*ws./u8)-0.7).*log(exp(k*ws./u8)-0.7);
ra=ra1./(ws*k*k);ga=1./ra; %[m s-1]
cp=1006;%[J kg-1 K-1]
ta=dat(:,3);
lv=2500-2.36*ta; %[kJ kg-1]
Pa=dat(:,4); %Atmospheric pressure
psy=(cp*Pa)./(0.622*lv*1000);%Psychrometric constant [kPa K-1]
% D=(4098*(0.6108*exp((17.27*ta)./(ta+237.3)))./((ta+237.3).*(ta+237.3)));%[kPa c-1],
D1=(2508*(exp((17.27*ta)./(ta+237.3)))./((ta+237.3).*(ta+237.3)));%[kPa k-1],
%Rn-G wm-2
Rn=dat(:,5); %net radiation
G=dat(:,6);  %soil heat flux
Q=Rn-G;
%Density of dry air
da=(101.3*1000)./(287.04*(ta+273.15));%[kg m-3]
lvE=dat(:,7);
vpd=0.1*dat(:,8);
fzi=ga.*psy;
fmu=(D1.*Q+da*cp.*ga.*vpd)./lvE-(D1+psy);
gc=fzi./fmu;  %m/s
gc(gc<0)=nan;
gc(gc>0.06)=nan;%m/s   %[0,2000 mmol m-2 s-1]

rc=1./gc;
ci=dat(:,12)-(dat(:,11).*(ra+1.6*rc)/44.6);%molCO2 mol-1%https://seer.sct.embrapa.br/index.php/agrometeoros/article/view/26527
%km calcualte the Michaelis-Menton coefficient 
patm = Pa*1000; %pa
rat = Pa/101.325;
R = 8.314;
O2 = 2.09476e5;
Kc25 = 41.03 * rat;
Ko25 = 28210 * rat;
Hkc = 79430;
Hko = 36380;
temp=ta;
temp_k = 273.15 + temp;
Kc_pa =Kc25.*exp(Hkc*((temp_k - 298.15)./ (298.15 * R * temp_k)));
Ko_pa =Ko25.*exp(Hko * ((temp_k - 298.15)./(298.15 * R * temp_k)));
O2_pa = O2 * (1e-6) * patm;
Km_pa = Kc_pa .* (1 + O2_pa./Ko_pa);
%gammastar: CO2 compensation point (Pa)
gammastar25 = 4.332 * rat;  % Pa
Hgm=37830; % J mol-1
R = 8.314;        % J K-1 mol-1
O2 = 2.09476e5; % ppm
O2_0 = O2 * 1e-6 * 101325;
O2_z = O2 * 1e-6 * patm;

gStar_pa = gammastar25.*(exp((Hgm/R)*(1/298.15-1./temp_k)));
Amax=dat(:,13);
vcmax=Amax.*(ci*101325+Km_pa)./(ci*101325-gStar_pa);
Ha= 55000;
Hd= 200000;
adelS= 663.1;
trefK=25+273.15;
R=8.314;
kbeg=exp(Ha*(temp_k-trefK)./(trefK*R*temp_k));
kend=((1+exp((trefK*(adelS)-Hd)/(trefK*R)))./(1+exp((temp_k*(adelS)-Hd)./(temp_k*R))));
vcmax25=vcmax./(kbeg.*kend);
