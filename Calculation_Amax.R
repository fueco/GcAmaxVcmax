
#The method for calculating the Amax is following Luo X, Keenan TF. Global evidence for the acclimation of ecosystem photosynthesis to light. Nat Ecol Evol, (2020). https://doi.org/10.1038/s41559-020-1258-7.
#The REddyProc R package (https://github.com/bgctw/REddyProc). 

rm(list=ls())

#please install packages if not yet installed
library(sirad)
library(chron)
library(REddyProc)
library(plyr)
library(mlegp)

for (xxx in 1:sitenumber) {
  setwd('.../fullset/fullsetflux2015/')  #set path
  
  # set the location of the data
  dataInHome<-'.../fullset/fullsetflux2015/' 
  siteList<-list.files(dataInHome, pattern='FULLSET', all.files=FALSE,
                       full.names=FALSE)
  
  siteName=siteList[xxx]
  cSiteShort<-substr(siteName,5,10)
  
  print(paste('###################### Loading data for ', cSiteShort))
  
  tmp<-read.table(siteName, header = TRUE, sep = ",")
  
  TIMESTAMP_START <- tmp$TIMESTAMP_START
  TIMESTAMP_END <- tmp$TIMESTAMP_END
  VPD_F <- tmp$VPD_F
  VPD_F_QC <- tmp$VPD_F_QC
  TA_F <- tmp$TA_F
  TA_F_QC <- tmp$TA_F_QC
  TS_F <- tmp$TS_F_MDS_1
  TS_F_QC <- tmp$TS_F_MDS_1_QC
  NIGHT <- tmp$NIGHT
  SW_IN_F <- tmp$SW_IN_F
  SW_IN_F_QC <- tmp$SW_IN_F_QC
  SW_IN_F_MDS <- tmp$SW_IN_F_MDS
  SW_IN_F_MDS_QC <- tmp$SW_IN_F_MDS_QC
  NEE_VUT_USTAR50 <- tmp$NEE_VUT_USTAR50
  NEE_VUT_USTAR50_QC <- tmp$NEE_VUT_USTAR50_QC
  NEE_VUT_USTAR50_RANDUNC <- tmp$NEE_VUT_USTAR50_RANDUNC
  
  tmpx<-NEE_VUT_USTAR50
  tmpy<-tmp$NEE_VUT_USTAR50_RANDUNC
  tmpy[tmpy<=-2000]<-NA # remove -9999s
  indX<-tmpx>=0
  
  myline.fit <- lm(tmpy[indX] ~ tmpx[indX])
  
  a<-myline.fit$coefficients[1]
  b<-myline.fit$coefficients[2]
  print(paste('###################### uncertainty parameters: ', a, b))
  
  NEE_VUT_USTAR50_RANDUNC[NEE_VUT_USTAR50>=0]<-a+NEE_VUT_USTAR50[NEE_VUT_USTAR50>=0]*b
  
  tmpx<-NEE_VUT_USTAR50
  tmpy <- tmp$NEE_VUT_USTAR50_RANDUNC
  tmpy[tmpy<=-2000]<-NA
  indX<-tmpx<0
  
  myline.fit <- lm(tmpy[indX] ~ tmpx[indX])
  
  a<-myline.fit$coefficients[1]
  b<-myline.fit$coefficients[2]
  print(paste('###################### uncertainty parameters: ', a, b))
  NEE_VUT_USTAR50_RANDUNC[NEE_VUT_USTAR50<0]<-a-abs(NEE_VUT_USTAR50[NEE_VUT_USTAR50<0])*b
  
  SW_IN_POT <- tmp$SW_IN_POT
  
  # Creating DAY Variable
  DAY<-1-NIGHT
  
  # Creating the data frame with the proper values of julian, day, month, year etc from the TIMESTAMP_START variable in FLUXNET dataframe
  timeStampChar <- as.character(TIMESTAMP_START)
  Data.F <- data.frame(
    year = as.integer(substr(timeStampChar,1,4))
    ,month = as.integer(substr(timeStampChar,5,6))
    ,day = as.integer(substr(timeStampChar,7,8))
    ,hour = as.integer(substr(timeStampChar,9,10)) +
      as.integer(substr(timeStampChar,11,12))/60)
  
  #Removing NA 
  Data.F<-Data.F[!is.na(as.integer(substr(timeStampChar,1,4))),]  
  Data.F$julian<-julian(Data.F$month,Data.F$day,Data.F$year)-julian(1,1,na.omit(unique(Data.F$year))[1])
  
  #Creating data frame with the variable of interest
  Data.F <- cbind(Data.F, 
                  data.frame(NEE=NEE_VUT_USTAR50,
                             NEE_QC=NEE_VUT_USTAR50_QC,
                             NEE_SE=NEE_VUT_USTAR50_RANDUNC,
                             SW_IN_F=SW_IN_F,
                             SW_IN_F_QC=SW_IN_F_QC,
                             SW_IN_POT=SW_IN_POT,
                             TA_F=TA_F,
                             TA_F_QC=TA_F_QC,
                             VPD_F=VPD_F,
                             VPD_F_QC=VPD_F_QC,
                             NIGHT=NIGHT,
                             DAY=DAY))
  
  #Convertime time in posix (not necessary but useful for the plots with REddyProc)
  dfall_posix  <- fConvertTimeToPosix(Data.F, 'YMDH', Year = 'year', Month='month', Day = 'day', Hour = 'hour')
  #dfall_posix  <- fConvertTimeToPosix(Data.F, format = "%Y/%m/%d %H:%M", Year.s = 'year', Month.s='month', Day.s = 'day', Hour.s = 'hour')
  
  hourz=unique(dfall_posix$hour)
  nRecInDay.i=length(hourz)
  
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # START - RUN THE REddyProc partitioning
  
  # Run the REddyProc NT partitioning
  #EddyProc.C <- sEddyProc(cSiteShort,dfall_posix, c('NEE','NEE_QC','SW_IN_F','TA_F','TA_F_QC','VPD_F'),ColPOSIXTime='DateTime',DTS=nRecInDay.i)
  #EddyProc.C$sSetLocationInfo(LatDeg=42.5, LongDeg=-72.19, TimeZoneHour=-5)  #Location of DE-Tharandt
  #EddyProc.C$sMRFluxPartition(FluxVar='NEE'      ##<< Variable name of column with original and filled net ecosystem fluxes (NEE)
  #                           ,QFFluxVar='NEE_QC' ##<< Quality flag of NEE variable
  #                            ,QFFluxValue=0       ##<< Value of quality flag for _good_ (original) data
  #                          ,TempVar='TA_F'    ##<< Filled air- or soil temperature variable (degC)
  #                            ,QFTempVar='TA_F_QC'##<< Quality flag of filled temperature variable
  #                            ,QFTempValue=0       ##<< Value of temperature quality flag for _good_ (original) data
  #                            ,RadVar='SW_IN_F'  )	# night time partitioning -> Reco, GPP
  #FilledEddyDataNT.F <- EddyProc.C$sExportResults()
  
  # Run the REddyProc DT partitioning at zero threshold
  print('###################### ')
  print('###################### ')
  print(cSiteShort)
  print('###################### Running Day-time partitioning')
  
  # Need to run on each year independently as the smoothing temperature sensitivity estimates code can not handle large arrays
  yearsUnique<-unique(Data.F$year)
  
  #yearsUnique<-yearsUnique[-1]
  for (ii in yearsUnique){
    print(c('running year',ii))
    indX<- Data.F$year==ii
    Data.FcYear<-Data.F[indX,]
    
    tmpPart<-NULL
    
    try(tmpPart <- partitionNEEGL(Data.FcYear,NEEVar="NEE",QFNEEVar="NEE_QC",QFNEEValue = 0,NEESdVar="NEE_SE",
                                  TempVar="TA_F",QFTempVar="TA_F_QC",VPDVar="VPD_F",QFVPDVar="VPD_F_QC",
                                  RadVar="SW_IN_F",QFRadVar="SW_IN_F_QC",PotRadVar="DAY",suffix="", nRecInDay= nRecInDay.i,
                                  controlGLPart=partGLControl(nBootUncertainty=10L, isAssociateParmsToMeanOfValids=FALSE, 
                                                              isLasslopPriorsApplied=TRUE, isBoundLowerNEEUncertainty = FALSE,
                                                              smoothTempSensEstimateAcrossTime=TRUE)), silent = FALSE)
    
    if(is.null(tmpPart)){
      # if the year returned no results, then create an empty result array 'emptyTmp'
      # need to do this to concatenate later
      RecoDTVar.s <- 'Reco_DT'
      GPPDTVar.s <- 'GPP_DT'
      RecoDTSdVar.s <- paste0(RecoDTVar.s,"_SD")
      GPPDTSdVar.s <- paste0(GPPDTVar.s,"_SD")
      
      emptyTmp <- data.frame(
        FP_VARnight=rep(NA_real_,nrow(Data.FcYear))	##<< NEE filtered for nighttime records (others NA)
        ,FP_VARday=NA_real_		##<< NEE filtered for daytime recores (others NA)
        ,NEW_FP_Temp=NA_real_	##<< temperature after filtering for quality flag degree Celsius
        ,NEW_FP_VPD=NA_real_	##<< vapour pressure deficit after filtering for quality flag, hPa
        ,FP_RRef_Night=NA_real_		##<< basal respiration estimated from LRC of daytime window  (W/m2)
        ,FP_qc =NA_integer_
        ,FP_dRecPar=NA_integer_
        ,FP_errorcode=NA_real_
        ,FP_GPP2000 =NA_real_   
        ,FP_OPT_VPD=NA_real_  
        ,FP_OPT_NoVPD=NA_real_
        ,FP_k=NA_real_
        ,FP_beta=NA_real_	
        ,FP_alpha=NA_real_
        ,FP_RRef=NA_real_		##<< basal respiration estimated from LRC of daytime window  (W/m2)
        ,FP_E0=NA_real_			##<< temperature sensitivity estimated from nighttime NEE window  in Kelvin (degK)
        ,FP_k_sd=NA_real_       
        ,FP_beta_sd=NA_real_    
        ,FP_alpha_sd=NA_real_   
        ,FP_RRef_sd=NA_real_    
        ,FP_E0_sd=NA_real_
        ,Reco_DT=NA_real_		##<< NEE filtered for daytime recores (others NA)
        ,GPP_DT=NA_real_		##<< NEE filtered for daytime recores (others NA)
        ,Reco_DT_SD=NA_real_		##<< NEE filtered for daytime recores (others NA)
        ,GPP_DT_SD=NA_real_		##<< NEE filtered for daytime recores (others NA)
      )
      attr(emptyTmp$FP_VARnight, 'varnames') <- paste('NEE.NEE_QC_0_night', sep='')
      attr(emptyTmp$FP_VARday, 'varnames') <- paste('NEE.NEE_QC_0_day', sep='')
      attr(emptyTmp[[RecoDTVar.s]], 'varnames') <- RecoDTVar.s
      attr(emptyTmp[[GPPDTVar.s]], 'varnames') <- GPPDTVar.s
      attr(emptyTmp[[RecoDTSdVar.s]], 'varnames') <- RecoDTSdVar.s
      attr(emptyTmp[[GPPDTSdVar.s]], 'varnames') <- GPPDTSdVar.s
      attr(emptyTmp$FP_VARnight, 'units') <- NULL
      attr(emptyTmp$FP_VARday, 'units') <- NULL
      attr(emptyTmp[[RecoDTVar.s]], 'units') <- NULL
      attr(emptyTmp[[GPPDTVar.s]], 'units') <- NULL
      attr(emptyTmp[[RecoDTSdVar.s]], 'units') <- NULL
      attr(emptyTmp[[GPPDTSdVar.s]], 'units') <- NULL
      
      tmpPart <- emptyTmp
    }
    
    if (ii==yearsUnique[1]){
      df.REddy<-tmpPart
    }else{
      df.REddy<-rbind(df.REddy,tmpPart)
    }
  }
  
  # replace NA with -9999
  Data.F[is.na(Data.F)] <- -9999
  
  #FilledEddyDataNT.F[is.na(FilledEddyDataNT.F)]<- -9999
  df.REddy[is.na(df.REddy)] <- -9999
  
  # write the night-time partitioned data to file
  #FilledEddyDataNT.F[is.na(FilledEddyDataNT.F)] <- -9999
  #write.table(FilledEddyDataNT.F, file = paste("./data_REddyProcOutput/",cSiteShort,"REddyProc_NT_","VUT_USTAR50",".csv",sep=""), append = FALSE, sep = ",")
  
  # write the multiple day-time partitioned data to file
  GPP_DT<-as.data.frame(cbind(Data.F$year,df.REddy$GPP_DT))
  names(GPP_DT)<-c("year","GPP_DT")
  Reco_DT<-as.data.frame(cbind(Data.F$year,df.REddy$Reco_DT))
  names(Reco_DT)<-c("year","Reco_DT")
  
  #write.table(GPP_DT, file = paste("./data_REddyProcOutput/",cSiteShort,"_GPP_DT_","VUT_USTAR50",".csv",sep=""), append = FALSE, sep = ",",row.names=F)
  #write.table(Reco_DT, file = paste("./data_REddyProcOutput/",cSiteShort,"_Reco_DT_","VUT_USTAR50",".csv",sep=""), append = FALSE, sep = ",",row.names=F)
  
  # write out the parameters from each paritioning experiment
  parameters<-cbind(Data.F$year,Data.F$month,Data.F$julian,Data.F$hour,Data.F$DAY,
                    df.REddy$FP_RRef_Night,df.REddy$FP_RRef,df.REddy$FP_E0,df.REddy$FP_alpha,df.REddy$FP_beta,df.REddy$FP_k,df.REddy$FP_qc)
  colnames(parameters) <- c("Year", "Month", "Day","Hour","DAYNIGHT","R_night","R_ref","E0","alpha","beta","k","qc")
  write.csv(parameters, file = paste(getwd(),"/data_REddyProcOutput/",cSiteShort,"REddyProc_parameters",".csv",sep=""),  row.names  = FALSE)
  rm(list=setdiff(ls(), "xxx"))
}

