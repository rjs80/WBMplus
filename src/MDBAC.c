/******************************************************************************
 
GHAAS Water Balance/Transport Model V3.0
Global Hydrologic Archive and Analysis System
Copyright 1994-2011, UNH - CCNY/CUNY

tao7753@gmail.com  

*******************************************************************************/
#include <stdio.h>
#include <cm.h>
#include <MF.h>
#include <MD.h>
#include <math.h>

// input
            
static int _MDInDischargeID            = MFUnset;
static int _MDInDischMeanID            = MFUnset;
static int _MDInwaterStorageID         = MFUnset;
static int _MDInContributingAreaID     = MFUnset;
static int _MDInDischarge0ID           = MFUnset;
static int _MDInRunoffVolumeID         = MFUnset;
static int _MDInRunoffID               = MFUnset;
static int _MDInRunoffVolID	       = MFUnset;
static int _MDInwaterTotalVolumeID     = MFUnset;
static int _MDInRiverStorageChgID      = MFUnset;
static int _MDInRiverStorageID         = MFUnset;
//static int _MDInLocalLoad_BacteriaID   = MFUnset;
static int _MDInRiverDepthID	       = MFUnset;
static int _MDInRiverWidthID           = MFUnset;
static int _MDInRiverbedWidthID        = MFUnset;
static int _MDInWidthAdjustID          = MFUnset;
static int _MDInRiverOrderID           = MFUnset;
static int _MDInQxT_WaterTempID        = MFUnset;
static int _MDInSolarRadiationID       = MFUnset;
//static int _MDInKsID		       = MFUnset;
static int _MDPrecipitationID          = MFUnset;
//static int _MDAPrecipitationID = MFUnset;
static int _MDAirTemperatureID         = MFUnset;
//static int _MDInBacteriaTotalInID             = MFUnset;
//static int _MDInBacteriaTotalInMixingID       = MFUnset;
static int _MDInLandUseSubID           = MFUnset;              // RJS 022714
static int _MDInLandUseForID           = MFUnset;
//static int _MDInWetlandsID             = MFUnset;  
//static int _MDInCellLengthID           = MFUnset;
//static int _MDCellLengthID           = MFUnset;
//static int _MDInPropROGroundWaterID    = MFUnset;
static int _MDAlphaHTSID               = MFUnset; //alpha exchange coefficient

// output

//static int _MDOutcross_AID                  = MFUnset;
static int _MDOutLocalLoad_BacteriaID       = MFUnset;
static int _MDFluxMixing_BacteriaID         = MFUnset;
static int _MDFlux_BacteriaID               = MFUnset;
static int _MDFluxMC_BacteriaID             = MFUnset;
//static int _MDFluxOnlyMC_BacteriaID         = MFUnset;
static int _MDpreFlux_BacteriaID            = MFUnset;
static int _MDOutpreFlux_BacteriaID         = MFUnset;
static int _MDOutpreFluxMC_BacteriaID       = MFUnset;
static int _MDOutConcMixing_BacteriaID      = MFUnset;
static int _MDStoreWaterMixing_BacteriaID   = MFUnset;
static int _MDStoreWater_BacteriaID         = MFUnset;

//static int _MDStoreWaterMixing_MUSLEID      = MFUnset;
//static int _MDStoreWater_MUSLEID            = MFUnset;

static int _MDpostConcMixing_BacteriaID     = MFUnset;
static int _MDOutpostConcMixing_BacteriaID  = MFUnset;
static int _MDOutpostConc_BacteriaID        = MFUnset;
static int _MDOutpostConcMC_BacteriaID      = MFUnset;
//static int _MDOutpostConcOnlyMC_BacteriaID  = MFUnset;
//static int _MDpostConcOnlyMC_BacteriaID  = MFUnset;
static int _MDOutpostConcHTS_BacteriaID     = MFUnset;
static int _MDRunoffConc_BacteriaID         = MFUnset;
static int _MDExport_BacteriaID         = MFUnset;
static int _MDInRiverVelocityID          = MFUnset;
static int _MDtotalMassRemoved_BacteriaID       = MFUnset;
static int _MDMCMassRemoved_BacteriaID          = MFUnset;
//static int _MDOnlyMCMassRemoved_BacteriaID      = MFUnset;

//static int _MDKsID                              = MFUnset;
//static int _MDtResiID                           = MFUnset;
//static int _MDOuttResiID                                     = MFUnset;
static int _MDflowPathRemoval_BacteriaID                     = MFUnset;
static int _MDflowPathRemovalMixing_BacteriaID               = MFUnset;
static int _MDflowPathRemovalMC_BacteriaID               = MFUnset;
static int _MDflowPathRemovalHTS_BacteriaID               = MFUnset;
static int _MDOutBacteriaTotalInID                           = MFUnset;
static int _MDOutBacteriaTotalInMixingID                     = MFUnset;
static int _MDHTSRemoval_BacteriaID                          = MFUnset;
//static int _MDTEHTSID                           = MFUnset; //The fractions entering HTS // which one?
//static int _MDOutTEHTSID                        = MFUnset; //The fractions entering HTS 
//static int _MDStoreWater_Bacteria2ID		= MFUnset;
static int _MDStoreWaterMC_BacteriaID		= MFUnset;
static int _MDOutpreConc_BacteriaID             = MFUnset; 
static int _MDOutpreConcMC_BacteriaID           = MFUnset;	// RJS 120408

static void _MDBacteriav2 (int itemID) {
    float HCWA                          = 0.0;
    //    float aA			        = 0.0;		// main channel area (m2)
    //float cross_A			= 0.0;		// main channel area (m2)
   // float wetlands                      = 0.0;
  //  float wetlandArea                   = 0.0;
    float luSub                         = 0.0;
    float luFor                         = 0.0;
    float RunoffVolume                  = 0.0;  // m3/s
    float RunoffVol                     = 0.0;  // m3/s
    float Runoff                        = 0.0;
    float waterStorage                  = 0.0;  
    float waterStorageChange                 = 0.0;  
    float LocalLoad_Bacteria                 = 0.0;  // CFU/d
    float preFluxMixing_Bacteria             = 0.0;  // CFU/d
    float preFluxMC_Bacteria                 = 0.0;  // CFU/d
    float postFluxMixing_Bacteria            = 0.0;  // CFU/d
    float storeWaterMixing_Bacteria          = 0.0;  // CFU/d
    float postStoreWaterMixing_Bacteria      = 0.0;  // CFU/d
    float discharge                          = 0.0;  // m3/s
    //float CellLength                          = 0.0;  // m3/s    
    //float PropROGroundWater             = 0.0;  
    float ContributingArea              = 0.0;  //?
    float dischargePre                       = 0.0;  // m3/s
    float waterStoragePrev                   = 0.0;  // m3/s  
    float BacteriaTotalInMixing              = 0.0;  // CFU/d
    float waterTotalVolume                   = 0.0;	// m3/d
    float flowPathRemovalMixing_Bacteria     = 0.0;  // CFU/d
    //float flowPathRemovalHTS_Bacteria        = 0.0;  // CFU/d
    float postConcMixing_Bacteria            = 0.0;  // CFU/100ml
    float massBalanceMixing_Bacteria         = 0.0;  // CFU/d
    float massBalanceMC_Bacteria             = 0.0; 
    float preConcMixing_Bacteria             = 0.0;  // CFU/100ml
    float preFlux_Bacteria                   = 0.0;  // CFU/d
    float storeWater_Bacteria                = 0.0;  // CFU/d
    float storeWaterMC_Bacteria                = 0.0;  // CFU/d 
    float BacteriaTotalIn                    = 0.0;  // CFU/d
    float BacteriaTotalInMC                    = 0.0;  // CFU/d    
    float preConc_Bacteria                   = 0.0;  // 
    float preConcMC_Bacteria                   = 0.0;  
    float HL                                 = 0.0;  // m/d
    float width                              = 0.0;  // m
    float widthAdjust                        = 0.0;

    float removal                            = 0.0;  // proportional removal
    float totalMassRemoved_Bacteria          = 0.0;  // CFU/d
    float MCMassRemoved_Bacteria             = 0.0;  // CFU/d
    //float OnlyMCMassRemoved_Bacteria         = 0.0;  // CFU/d
    float flowPathRemoval_Bacteria           = 0.0;  // CFU/d
    float flowPathRemovalMC_Bacteria           = 0.0;  // CFU/d
    float postConc_Bacteria                  = 0.0;  // CFU/100ml
    float postConcMC_Bacteria                = 0.0;  // CFU/100ml
    //float postConcOnlyMC_Bacteria            = 0.0;  // CFU/100ml
    float postConcHTS_Bacteria               = 0.0;  // CFU/100ml
    float postFlux_Bacteria                  = 0.0;  // CFU/d
    float postFluxMC_Bacteria                = 0.0;  // CFU/d
    //float postFluxOnlyMC_Bacteria            = 0.0;  // CFU/d
    float postStoreWater_Bacteria            = 0.0;  // CFU/d
    float postStoreWaterMC_Bacteria          = 0.0;  // CFU/d
    //float postStoreWaterOnlyMC_Bacteria      = 0.0;  // CFU/d
    float massBalance_Bacteria               = 0.0;  // CFU/d
    float massBalance_Q                      = 0.0;
    float order                              = 0.0;  // river order
    //float length                             = 0.0;  // m
    float cell_area                     = 0.0;  // km2
    float RunoffConc_Bacteria                = 0.0;  //CFU/100ml
    float tResi                 = 0.0; // Resident time, day
    float dL			= 0.0;		// cell length (m)
    float Velocity		= 0.0; 
    float QxT_WaterTemp         =0.0; 
    float A =0.0;  // Temperature adjustment factor 
    float SolarRadiation =0.0;   
    //float Ks             =0.0;
    float Precipitation       = 0.0;
    //float APrecipitation = 0.0;
    float TEHTS               = 0.0; 
    float AirTemperature      = 0.0; 
    float AlphaHTS            = 0.0; 
    float HTSRemoval_Bacteria = 0.0; 
    float tStorHZ	      = 0.0;		// time of storage or residence time in HZ (days)
    
    
 //   wetlands             = MFVarGetFloat (_MDInWetlandsID,           itemID, 0.0);   // proportion wetlands
    luSub                = MFVarGetFloat (_MDInLandUseSubID,         itemID, 0.0) / 100;
    luFor                = MFVarGetFloat (_MDInLandUseForID,         itemID, 0.0);
    //CellLength           = MFVarGetFloat (_MDInCellLengthID,         itemID, 0.0);
    dL                   = MFModelGetLength (itemID);               // km converted to m
    discharge            = MFVarGetFloat (_MDInDischargeID,          itemID, 0.0); // m3/sec, discharge leaving the grid cell, after routing!          
    //PropROGroundWater    = MFVarGetFloat (_MDInPropROGroundWaterID,  itemID, 0.0);
    ContributingArea     = MFVarGetFloat (_MDInContributingAreaID,          itemID, 0.0); // m3/sec, discharge leaving the grid cell, after routing!
    dischargePre	 = MFVarGetFloat (_MDInDischarge0ID,         itemID, 0.0); // m3/sec, discharge from upstream PLUS local runoff, before routing!
    RunoffVolume         = MFVarGetFloat (_MDInRunoffVolumeID,       itemID, 0.0); // m3/s
    RunoffVol            = MFVarGetFloat (_MDInRunoffVolID,      itemID, 0.0); // m3/sec
               Runoff          = MFVarGetFloat (_MDInRunoffID,      itemID, 0.0); 
                waterStorageChange   = MFVarGetFloat (_MDInRiverStorageChgID,    itemID, 0.0); // m3/s
		waterStorage         = MFVarGetFloat (_MDInRiverStorageID,       itemID, 0.0); // m3/s           
                preFlux_Bacteria           = MFVarGetFloat (_MDFlux_BacteriaID,             itemID, 0.0);	// CFU/day RJS 091108
                preFluxMC_Bacteria         = MFVarGetFloat (_MDFluxMC_BacteriaID,           itemID, 0.0);
                preFluxMixing_Bacteria     = MFVarGetFloat (_MDFluxMixing_BacteriaID,       itemID, 0.0); // CFU/day 
                storeWater_Bacteria        = MFVarGetFloat (_MDStoreWater_BacteriaID,       itemID, 0.0);	// CFU/day RJS 091108
                storeWaterMC_Bacteria        = MFVarGetFloat (_MDStoreWaterMC_BacteriaID,       itemID, 0.0);
                storeWaterMixing_Bacteria  = MFVarGetFloat (_MDStoreWaterMixing_BacteriaID, itemID, 0.0); // CFU/day
                order                      = MFVarGetFloat (_MDInRiverOrderID,           itemID, 0.0); // strahler order
                cell_area            = MFModelGetArea(itemID); // (unit?)
                waterStoragePrev     = waterStorage - waterStorageChange;                                    // m3/sec     
                waterTotalVolume         = (discharge + waterStorage) * 86400;                                   // m3/d
              //  HL                       = width > 0.00001 ? discharge / (width * length) * 86400 : 0.0;      // m/d
               Precipitation		 = MFVarGetFloat(_MDPrecipitationID,              itemID, 0.0);  
               //APrecipitation		 = MFVarGetFloat(_MDAPrecipitationID,              itemID, 0.0);
               AirTemperature		 = MFVarGetFloat(_MDAirTemperatureID,              itemID, 0.0);
               Velocity                  = MFVarGetFloat(_MDInRiverVelocityID,            itemID,0.0);          
          //RunoffConc_Bacteria=pow(10, 1.2339+0.0538*Precipitation+0.0343*AirTemperature+0.0129*luSub*100);
          //RunoffConc_Bacteria=pow(10, -1.0044+0.0539*Precipitation+0.0343*AirTemperature+0.0346*luSub*100+0.0372*(luDec+luCon+luMix)*100);
          
         // Can we have a citation for this relationship?          
         Precipitation=Precipitation > 25.0 ? 25.0 : Precipitation; 
         RunoffConc_Bacteria=pow(10, 0.8676+0.0488*Precipitation+0.0462*AirTemperature+0.0139*luSub*100+0.0052*(luFor)*100);
          
         // RunoffConc_Bacteria=pow(10, 0.8676+0.0488*Precipitation+0.0462*AirTemperature+0.0139*luSub*100+0.0052*(1-luSub));
          
          //RunoffConc_Bacteria=1;
                
                //RunoffConc_Bacteria=pow(10,0.76+0.02*impAreaFrac*100+0.04*Precipitation+0.01*APrecipitation+0.01*AirTemperature+2.59*USLE_K);
          LocalLoad_Bacteria        = RunoffConc_Bacteria /100* 1000000* RunoffVol *86400; 
                //LocalLoad_Bacteria        = RunoffVol *86400; 
                //BacteriaTotalIn           = LocalLoad_Bacteria + preFlux_Bacteria + storeWater_Bacteria;                    // CFU/d  
                BacteriaTotalIn           = LocalLoad_Bacteria + preFlux_Bacteria + storeWater_Bacteria;                    // CFU/d                  
                BacteriaTotalInMC         = LocalLoad_Bacteria + preFluxMC_Bacteria + storeWaterMC_Bacteria;                 
                BacteriaTotalInMixing     = LocalLoad_Bacteria + preFluxMixing_Bacteria + storeWaterMixing_Bacteria;        // CFU/d                 
 
                	

               // cross_A              =4.6883*pow(discharge,0.611);
                  //Velocity	     = discharge/cross_A; 
                
                tResi                = (dL) / Velocity / 86400;             // days,                
              //  tResi                = (dL) / Velocity / 86400;
                
                if (waterTotalVolume > 0.0001) {
              //  wetlandArea          = ((MFModelGetArea (itemID) * wetlands) - (dL * width)) * HCWA;                    // m2 
               // wetlandArea          = wetlandArea < 0.0 ? MFModelGetArea (itemID) * wetlands * HCWA : wetlandArea;     // m2
                    
                  preConc_Bacteria          = BacteriaTotalIn / waterTotalVolume / 1000/1000*100;                          // total_in: CFU/d , water total volume: m3/d, CFU/100ml          
                  preConcMC_Bacteria        = BacteriaTotalInMC / waterTotalVolume / 1000/1000*100;
                  preConcMixing_Bacteria    = BacteriaTotalInMixing / waterTotalVolume / 10000;                    // CFU/100ml      

                       // if (order > 2.0) {

                  // Can we have a citation for this?
                         A=pow(1*1.07, (QxT_WaterTemp-20)); //Merr_Bac1225: 1.18, default: 1.07, Merr0117: 0.8                       
                            //Ks=SolarRadiation*1; // Cho et al., 2012, solar intensity coefficient [m2/MJ/day] =1.                            
                           // removal              = 1.0 - pow(2.718281828, -(1.0+Ks) * (tResi)*A);  // proportional removal, k20=1, Ɵ=1.02-1.08 (mattews 2007)
                            //removal              = 1.0 - pow(2.718281828, -1*0.8*(tResi)*A);// J Iudicello, DA Chin, 2014, ASSUME TrESI=0.1
                            removal              = 1.0 - pow(2.718281828, -1*0.8*(tResi)*A);
                            //removal   =0;
                                //totalMassRemoved_Bacteria = waterStorage/waterTotalVolume*BacteriaTotalIn+removal * BacteriaTotalIn*(1-waterStorage/waterTotalVolume);                            // CFU/d
                            //RiverWidthLmp=18.8*pow(discharge,0.1684);
                            //TEHTS=9.5/1000000*9.56*0.45*pow(discharge,0.82)*dL/discharge;//AlphaHTS=9.5/1000000
                           // TEHTS=9.5/1000000*aA*dL/discharge;//AlphaHTS=9.5/1000000
                            
                            //TEHTS=9.5/1000000*width*depth*dL/discharge;//AlphaHTS=9.5/1000000
                            
                            //TEHTS=1*9.5/1000000*4.6883*pow(discharge,0.611)*dL/discharge > 1.0 ? 1.0 :1*9.5/1000000*4.6883*pow(discharge,0.611)*dL/discharge;//oyster
                            //TEHTS=1*42.8/1000000*4.6883*pow(discharge,0.611)*dL/discharge > 1.0 ? 1.0 :1*42.8/1000000*4.6883*pow(discharge,0.611)*dL/discharge;
                            
                            TEHTS=9.5/1000000*tResi*86400;//AlphaHTS=9.5/1000000
                            //TEHTS=0.49/1000000*5*pow(discharge,0.5358)*dL/discharge;//AlphaHTS=9.5/1000000
                            MCMassRemoved_Bacteria=removal * BacteriaTotalInMC;
                            //MCMassRemoved_Bacteria=0;
                            //OnlyMCMassRemoved_Bacteria=removal * BacteriaTotalIn;
                            
                            //tStorHZ    = (0.35 / (9.5/1000000) ) / 86400;
                            
                            //HTSRemoval_Bacteria    = TEHTS > 1? BacteriaTotalIn : TEHTS* BacteriaTotalIn;

                            HTSRemoval_Bacteria    = TEHTS* BacteriaTotalIn;

                            totalMassRemoved_Bacteria = removal*BacteriaTotalIn*(1-TEHTS)+HTSRemoval_Bacteria ; 

                     //   }
                  
                   postConcMixing_Bacteria   = (BacteriaTotalInMixing-flowPathRemovalMixing_Bacteria)/ waterTotalVolume / 10000;                   // CFU/100ml      
                }
                
                else {
                  
                  flowPathRemoval_Bacteria       = BacteriaTotalIn;                       // CFU/d       
                  flowPathRemovalMC_Bacteria= BacteriaTotalInMC;
                  flowPathRemovalMixing_Bacteria = BacteriaTotalInMixing;                 // CFU/d
            
                  postConc_Bacteria              = 0.0;                              // CFU/100ml                               
                  postConcMC_Bacteria            = 0.0;                              // CFU/100ml
                  //postConcOnlyMC_Bacteria        = 0.0;
                  postConcHTS_Bacteria           = 0.0;                              // CFU/100ml
                  postConcMixing_Bacteria        = 0.0;                              // CFU/100ml
    
                }
   
                postConc_Bacteria             = waterTotalVolume > 0.0001 ? (BacteriaTotalIn - totalMassRemoved_Bacteria - flowPathRemoval_Bacteria) / waterTotalVolume / 10000 : 0.0;     // CFU/100ml
                postConcMC_Bacteria           = waterTotalVolume > 0.0001 ? (BacteriaTotalInMC - MCMassRemoved_Bacteria- flowPathRemovalMC_Bacteria)/ waterTotalVolume  / 10000 : 0.0;     // CFU/100ml
                //postConcOnlyMC_Bacteria       = waterTotalVolume > 0.0001 ? (BacteriaTotalIn - OnlyMCMassRemoved_Bacteria)/ waterTotalVolume  / 10000 : 0.0;     // CFU/100ml
                //postConcHTS_Bacteria          = waterTotalVolume > 0.0001 ? (BacteriaTotalIn - HTSRemoval_Bacteria- flowPathRemovalHTS_Bacteria)/ waterTotalVolume / 10000 : 0.0;
                postFlux_Bacteria 	      = (discharge * MDConst_m3PerSecTOm3PerDay) * postConc_Bacteria * 1000000 / 100;                                                        // CFU/day
                postFluxMC_Bacteria 	      = (discharge * MDConst_m3PerSecTOm3PerDay) * postConcMC_Bacteria * 1000000 / 100;
                //postFluxOnlyMC_Bacteria       = (discharge * MDConst_m3PerSecTOm3PerDay) * postConcOnlyMC_Bacteria * 1000000 / 100;
                postFluxMixing_Bacteria       = (discharge * MDConst_m3PerSecTOm3PerDay) * postConcMixing_Bacteria * 1000000 / 100;                                                  // CFU/day
                //postStoreWater_Bacteria       = waterStorage * postConc_Bacteria *1000000/ 100;            // CFU/day
                //postStoreWaterMC_Bacteria     = waterStorage * postConcMC_Bacteria *1000000/ 100;            // CFU/day
                //postStoreWaterOnlyMC_Bacteria = waterStorage * postConcOnlyMC_Bacteria *1000000/ 100;  

                //postStoreWaterMixing_Bacteria = waterStorage * postConcMixing_Bacteria *1000000/ 100;	 // CFU/day
                postStoreWater_Bacteria       = (waterStorage * MDConst_m3PerSecTOm3PerDay) * postConc_Bacteria *1000000/ 100;            // CFU/day
                postStoreWaterMC_Bacteria     = (waterStorage * MDConst_m3PerSecTOm3PerDay) * postConcMC_Bacteria *1000000/ 100;         // CFU/day
                postStoreWaterMixing_Bacteria = (waterStorage * MDConst_m3PerSecTOm3PerDay) * postConcMixing_Bacteria *1000000/ 100;	 // CFU/day
                
                //massBalance_Q = Precipitation - (etp + runoff + grdWaterChg + snowPackChg + soilMoistChg );
                
                //massBalance_Bacteria       = BacteriaTotalIn > 0.00001 ? (BacteriaTotalIn - (postFlux_Bacteria + postStoreWater_Bacteria + flowPathRemoval_Bacteria + totalMassRemoved_Bacteria)) : 0.0;    // proportion of total CFU in
                //massBalanceMixing_Bacteria = (BacteriaTotalInMixing - (postFluxMixing_Bacteria + postStoreWaterMixing_Bacteria + flowPathRemovalMixing_Bacteria)) / BacteriaTotalInMixing;                          // proportion of total CFU in

                massBalance_Bacteria       = BacteriaTotalIn > 0.1 ? (BacteriaTotalIn - (postFlux_Bacteria + postStoreWater_Bacteria + totalMassRemoved_Bacteria + flowPathRemoval_Bacteria)) / BacteriaTotalIn : 0.0;    // proportion of total CFU in
                massBalanceMC_Bacteria       = BacteriaTotalInMC > 0.1 ? (BacteriaTotalInMC - (postFluxMC_Bacteria + postStoreWaterMC_Bacteria + MCMassRemoved_Bacteria + flowPathRemovalMC_Bacteria)) / BacteriaTotalInMC : 0.0;    // proportion of total CFU in
                massBalanceMixing_Bacteria = (BacteriaTotalInMixing - (postFluxMixing_Bacteria + postStoreWaterMixing_Bacteria + flowPathRemovalMixing_Bacteria)) / BacteriaTotalInMixing;                          // proportion of total CFU in
                
                if (massBalance_Bacteria > 0.0001) {
                    printf("itemID = %d, %d-%d-%d, MB_Bacteria = %f\n",itemID, MFDateGetCurrentYear(), MFDateGetCurrentMonth(), MFDateGetCurrentDay(), massBalance_Bacteria);                  
                    printf("Q = %f, BacteriaTotalIn=%f, postFluxMixing_Bacteria=%f, postFlux_Bacteria=%f, postStoreWater_Bacteria=%f, flowPathRemoval_Bacteria=%f,totalMassRemoved_Bacteria=%f,removal=%f,TEHTS=%f,dL=%f,luSub=%f\n",discharge,BacteriaTotalIn,postFluxMixing_Bacteria, postFlux_Bacteria,postStoreWater_Bacteria,flowPathRemoval_Bacteria, totalMassRemoved_Bacteria, removal,TEHTS,dL,luSub);
                   // printf("LocalLoad_Bacteria=%f ,postStoreWaterMixing_Bacteria=%f, storeWater_Bacteria=%f, preCpnc_Bacteria, postConc_Bacteria=%f, waterStoragePrev=%f, waterStorage=%f, waterStorageChange=%f\n",LocalLoad_Bacteria ,postStoreWaterMixing_Bacteria, storeWater_Bacteria, preConc_Bacteria, postConc_Bacteria, waterStoragePrev, waterStorage, waterStorageChange );
                }              
                
               // printf("itemID = %d....", itemID);
               
                 //if (itemID == 102) printf("*** itemID = %d, %d-%d-%d, BacteriaTotalInMixing=%f, preFluxMixing_Bacteria = %f, postFluxMixing_Bacteria = %f, LocalLoad_Bacteria = %f, storeMixingWater_Bacteria = %f\n",itemID, MFDateGetCurrentYear(), MFDateGetCurrentMonth(), MFDateGetCurrentDay(), BacteriaTotalInMixing, preFluxMixing_Bacteria, postFluxMixing_Bacteria, LocalLoad_Bacteria, storeWaterMixing_Bacteria);
                
                //if (itemID == 101) {
                    //printf("*** itemID = %d, %d-%d-%d, BacteriaTotalInMixing=%f, preFluxMixing_Bacteria = %f, postFluxMixing_Bacteria = %f, LocalLoad_Bacteria = %f, storeMixingWater_Bacteria = %f\n",itemID, MFDateGetCurrentYear(), MFDateGetCurrentMonth(), MFDateGetCurrentDay(), BacteriaTotalInMixing, preFluxMixing_Bacteria, postFluxMixing_Bacteria, LocalLoad_Bacteria, storeWaterMixing_Bacteria);
                    //printf("totalMassRemoved_Bacteria = %f, flowPathRemoval_Bacteria = %f, postConc_Bacteria=%f, postFlux_Bacteria = %f, postStoreWater_Bacteria = %f\n", totalMassRemoved_Bacteria, flowPathRemoval_Bacteria, postConc_Bacteria, postFlux_Bacteria, postStoreWater_Bacteria);
                    //printf("dischargePre = %f, discharge = %f, waterStorage = %f, waterStoragePrev = %f\n", dischargePre, discharge, waterStorage, waterStoragePrev);
                    //printf("Q = %f, BacteriaTotalIn=%f, postFlux_Bacteria=%f, postStoreWater_Bacteria=%f, flowPathRemoval_Bacteria=%f,totalMassRemoved_Bacteria=%f\n",discharge,BacteriaTotalIn,postFlux_Bacteria,postStoreWater_Bacteria,flowPathRemoval_Bacteria, totalMassRemoved_Bacteria);
		    //printf("waterTotalVolume = %f, HTSRemoval_Bacteria = %f, TEHTS = %f, MCMassRemoved_Bacteria = %f\n",waterTotalVolume, HTSRemoval_Bacteria, TEHTS, MCMassRemoved_Bacteria);
                    //printf("removal = %f, tResi = %f, A = %f, MB_Bacteria = %f\n", removal, tResi, A, massBalance_Bacteria);
                    //printf("preFlux_Bacteria = %f, storeWater_Bacteria = %f\n", preFlux_Bacteria, storeWater_Bacteria);
                    //printf("dL = %f\n", dL);
                              
                //}

               // MFVarSetFloat (_MDOutpreFluxMixing_BacteriaID,       itemID, preFluxMixing_Bacteria);                
                MFVarSetFloat (_MDFluxMixing_BacteriaID,             itemID, postFluxMixing_Bacteria);
                      
                MFVarSetFloat (_MDFlux_BacteriaID,                   itemID, postFlux_Bacteria);
                               
                MFVarSetFloat (_MDOutpreFluxMC_BacteriaID,           itemID, preFluxMC_Bacteria);
                MFVarSetFloat (_MDMCMassRemoved_BacteriaID,          itemID, MCMassRemoved_Bacteria); 
                
                MFVarSetFloat (_MDFluxMC_BacteriaID,                 itemID, postFluxMC_Bacteria);
                          
                //MFVarSetFloat (_MDFluxOnlyMC_BacteriaID,             itemID, postFluxOnlyMC_Bacteria);
                MFVarSetFloat (_MDOutLocalLoad_BacteriaID,           itemID, LocalLoad_Bacteria);
                MFVarSetFloat (_MDOutConcMixing_BacteriaID,          itemID, postConcMixing_Bacteria);
                MFVarSetFloat (_MDOutpostConc_BacteriaID,            itemID, postConc_Bacteria);
                MFVarSetFloat (_MDOutpostConcMC_BacteriaID,          itemID, postConcMC_Bacteria);
                //MFVarSetFloat (_MDOutpostConcOnlyMC_BacteriaID,      itemID, postConcOnlyMC_Bacteria);
                //MFVarSetFloat (_MDpostConcOnlyMC_BacteriaID,         itemID, postConcOnlyMC_Bacteria);
                MFVarSetFloat (_MDOutpostConcHTS_BacteriaID,         itemID, postConcHTS_Bacteria);
                MFVarSetFloat (_MDOutpostConcMixing_BacteriaID,      itemID, postConcMixing_Bacteria);
                MFVarSetFloat (_MDStoreWaterMixing_BacteriaID,       itemID, postStoreWaterMixing_Bacteria);
                MFVarSetFloat (_MDStoreWater_BacteriaID,             itemID, postStoreWater_Bacteria);
                MFVarSetFloat (_MDStoreWaterMC_BacteriaID,             itemID, postStoreWaterMC_Bacteria);
		//MFVarSetFloat (_MDStoreWater_Bacteria2ID,            itemID, storeWater_Bacteria);
                MFVarSetFloat (_MDRunoffConc_BacteriaID,             itemID, RunoffConc_Bacteria);
               // MFVarSetFloat (_MDMCMassRemoved_BacteriaID,          itemID, MCMassRemoved_Bacteria);
                //MFVarSetFloat (_MDOnlyMCMassRemoved_BacteriaID,      itemID, OnlyMCMassRemoved_Bacteria);
                //                MFVarSetFloat (_MDOutpreFlux_BacteriaID,             itemID, preFlux_Bacteria);
                MFVarSetFloat (_MDtotalMassRemoved_BacteriaID,       itemID, totalMassRemoved_Bacteria);   
              //  MFVarSetFloat (_MDOutpreFluxHTS_BacteriaID,           itemID, preFluxHTS_Bacteria);
               // MFVarSetFloat (_MDHTSRemoval_BacteriaID,             itemID, HTSRemoval_Bacteria);
               // MFVarSetFloat (_MDKsID,             itemID, Ks); 
                MFVarSetFloat (_MDflowPathRemoval_BacteriaID,             itemID, flowPathRemoval_Bacteria); 
                MFVarSetFloat (_MDflowPathRemovalMC_BacteriaID,             itemID, flowPathRemovalMC_Bacteria); 
                MFVarSetFloat (_MDflowPathRemovalMixing_BacteriaID,       itemID, flowPathRemovalMixing_Bacteria); 
                //MFVarSetFloat (_MDTEHTSID,             itemID, TEHTS); 
                //              MFVarSetFloat (_MDOutTEHTSID,             itemID, TEHTS); 
                MFVarSetFloat (_MDOutBacteriaTotalInID,             itemID, BacteriaTotalIn); 
//                MFVarSetFloat (_MDOutBacteriaTotalInMixingID,             itemID, BacteriaTotalInMixing);
               // MFVarSetFloat (_MDOutcross_AID,                      itemID, cross_A);                
                //MFVarSetFloat (_MDInwaterStorageID,             itemID, waterStorage); 

}

int MDBACDef () {
	MFDefEntering ("Bacteriav2 Routing");
	
        if (
             //((_MDInWetlandsID               = MFVarGetID (MDVarFracWetlandArea,                             "-",    MFInput,  MFState, MFBoundary))   == CMfailed) ||           // RJS 112513
                ((_MDInLandUseSubID                 = MFVarGetID (MDVarLandUseSpatialSub,                           "-",    MFInput,  MFState, MFBoundary))   == CMfailed) ||			//RJS 112211
                ((_MDInLandUseForID                 = MFVarGetID (MDVarFracForestedArea,                           "-",    MFInput,  MFState, MFBoundary))   == CMfailed) ||
              //  ((_MDInCellLengthID                 = MFVarGetID (MDVarCellLength,                           "m",    MFInput,  MFState, MFBoundary))   == CMfailed) ||
                    
                ((_MDInDischarge0ID                 = MFVarGetID (MDVarDischarge0,                               "m3/s",    MFInput,  MFState, MFBoundary))   == CMfailed) ||	
                
                ((_MDInContributingAreaID           = MFVarGetID (MDVarContributingArea,                       "m2",      MFInput,  MFState, MFBoundary))   == CMfailed) ||	
                ((_MDInDischargeID                  = MFVarGetID (MDVarDischarge,                              "m3/s",    MFInput,  MFState, MFBoundary))   == CMfailed) ||
                ((_MDInDischMeanID                  = MFVarGetID (MDVarDischMean,                              "m3/s",    MFInput,  MFState, MFBoundary))   == CMfailed) ||
                ((_MDInRiverVelocityID              = MFVarGetID (MDVarRiverVelocity,                            "-",       MFInput,  MFState, MFBoundary))   == CMfailed) ||          
                ((_MDInRiverStorageID               = MFVarGetID (MDVarRiverStorage,                           "m3/day",    MFInput,  MFState, MFInitial))    == CMfailed) ||	
                ((_MDInRiverStorageChgID            = MFVarGetID (MDVarRiverStorageChg,                        "m3/day",    MFInput,  MFState, MFBoundary))   == CMfailed) ||    
                ((_MDInRunoffVolumeID               = MFVarGetID (MDVarRunoffVolume, 		               "m3/s",      MFInput,  MFState, MFBoundary))   == CMfailed) ||       
                ((_MDInRunoffVolID                  = MFVarGetID (MDVarRunoffVolume, 		               "m3/s",      MFInput,  MFState, MFBoundary))   == CMfailed) ||
                ((_MDInRunoffID                     = MFVarGetID (MDVarRunoff, 		                       "mm/day",    MFInput,  MFState, MFBoundary))   == CMfailed) ||
                ((_MDInRiverOrderID                 = MFVarGetID (MDVarRiverOrder,                                "-",    MFInput,  MFState, MFBoundary))   == CMfailed) ||
                ((_MDInRiverDepthID                 = MFVarGetID (MDVarRiverDepth,                                "-",    MFInput,  MFState, MFBoundary))   == CMfailed) ||
                ((_MDRunoffConc_BacteriaID          = MFVarGetID (MDVarRunoffConc_Bacteria,               "CFU/100mL",     MFOutput,  MFState, MFInitial))    == CMfailed) ||	               
                ((_MDOutLocalLoad_BacteriaID        = MFVarGetID (MDVarLocalLoad_Bacteria,               "CFU/d",   MFOutput,  MFFlux,  MFBoundary))   == CMfailed) ||
                
                //((_MDOutLocalLoad_BacteriaID        = MFVarGetID (MDVarLocalLoad_Bacteria,               "CFU/d",   MFOutput,  MFFlux,  MFBoundary))   == CMfailed) ||
                
                //((_MDOutcross_AID                     	= MFVarGetID (MDVarcross_A,                            "m2",   MFOutput,  MFState, MFBoundary))   == CMfailed) || 
                
                ((_MDOutConcMixing_BacteriaID         	= MFVarGetID (MDVarConcMixing_Bacteria,    	        "CFU/100ml",   MFOutput,  MFState, MFBoundary))   == CMfailed) ||	  
                ((_MDStoreWaterMixing_BacteriaID      	= MFVarGetID (MDVarStoreWaterMixingBacteria,          "CFU/day",   MFOutput,  MFState, MFInitial))    == CMfailed) ||	
                ((_MDFluxMixing_BacteriaID         	= MFVarGetID (MDVarFluxMixing_Bacteria,    	      "CFU/day",   MFRoute,   MFFlux,  MFBoundary))   == CMfailed) ||	 
                ((_MDOutpreFluxMC_BacteriaID         	= MFVarGetID (MDVarpreFluxMC_Bacteria,    	      "CFU/day",   MFRoute,   MFFlux,  MFBoundary))   == CMfailed) ||	                 
                ((_MDOutpostConc_BacteriaID         	= MFVarGetID (MDVarpostConc_Bacteria,    	        "CFU/100mL",   MFOutput,  MFState, MFBoundary))   == CMfailed) ||	  
                ((_MDOutpostConcMC_BacteriaID         	= MFVarGetID (MDVarpostConcMC_Bacteria,    	        "CFU/100mL",   MFOutput,  MFState, MFBoundary))   == CMfailed) ||	  
                //((_MDOutpostConcOnlyMC_BacteriaID       = MFVarGetID (MDVarpostConcOnlyMC_Bacteria,    	        "CFU/100mL",   MFOutput,  MFState, MFBoundary))   == CMfailed) ||
                ((_MDOutpostConcHTS_BacteriaID         	= MFVarGetID (MDVarpostConcHTS_Bacteria,    	        "CFU/100mL",   MFOutput,  MFState, MFBoundary))   == CMfailed) ||
                ((_MDOutpostConcMixing_BacteriaID        = MFVarGetID (MDVarpostConcMixing_Bacteria,    	"CFU/100mL",   MFOutput,  MFState, MFBoundary))   == CMfailed) ||	  
                
                ((_MDOutBacteriaTotalInID        = MFVarGetID (MDVarBacteriaTotalIn,    	"CFU/day",   MFOutput,  MFState, MFBoundary))   == CMfailed) ||	  
                
                ((_MDStoreWater_BacteriaID               = MFVarGetID (MDVarStoreWater_Bacteria,             "CFU/day",   MFOutput,  MFState, MFInitial))    == CMfailed) ||	
                ((_MDStoreWaterMC_BacteriaID              = MFVarGetID (MDVarStoreWaterMC_Bacteria,             "CFU/day",   MFOutput,  MFState, MFInitial))    == CMfailed) ||	
                ((_MDFlux_BacteriaID                     = MFVarGetID (MDVarFlux_Bacteria,                     "CFU/day",   MFRoute,   MFFlux,  MFBoundary))   == CMfailed) ||	
                ((_MDFluxMC_BacteriaID                   = MFVarGetID (MDVarFluxMC_Bacteria,                   "CFU/day",   MFRoute,   MFFlux,  MFBoundary))   == CMfailed) ||
                //((_MDFluxOnlyMC_BacteriaID               = MFVarGetID (MDVarFluxOnlyMC_Bacteria,                   "CFU/day",   MFRoute,   MFFlux,  MFBoundary))   == CMfailed) ||
                ((_MDMCMassRemoved_BacteriaID            = MFVarGetID (MDVarMCMassRemoved_Bacteria,            "CFU/day",   MFOutput,   MFFlux,  MFBoundary))   == CMfailed) ||
                //((_MDOnlyMCMassRemoved_BacteriaID        = MFVarGetID (MDVarOnlyMCMassRemoved_Bacteria,               "CFU/day",   MFOutput,   MFFlux,  MFBoundary))   == CMfailed) ||
                ((_MDtotalMassRemoved_BacteriaID         = MFVarGetID (MDVartotalMassRemoved_Bacteria,                "CFU/day",   MFOutput,   MFFlux,  MFBoundary))   == CMfailed) ||	 
                //((_MDKsID                           = MFVarGetID (MDVarKs,                   "1/day",   MFOutput,   MFFlux,  MFBoundary))   == CMfailed) ||               
                ((_MDPrecipitationID                = MFVarGetID (MDVarPrecipitation,                   "mm",   MFInput,   MFFlux,  MFBoundary))   == CMfailed) ||          
                //((_MDAPrecipitationID                = MFVarGetID (MDVarAPrecipitation,                   "mm",   MFInput,   MFFlux,  MFBoundary))   == CMfailed) ||
                ((_MDAirTemperatureID                  = MFVarGetID (MDVarAirTemperature,                   "degC",   MFInput,   MFState,  MFBoundary))   == CMfailed) ||
                ((_MDInQxT_WaterTempID                 = MFVarGetID (MDVarWTemp_QxT,                   "degC",   MFInput,   MFState,  MFBoundary))   == CMfailed) ||
                //                ((_MDOuttResiID                        = MFVarGetID (MDVartResi,                   "d",   MFOutput,   MFState,  MFBoundary))   == CMfailed) ||	 
                ((_MDflowPathRemoval_BacteriaID        = MFVarGetID (MDVarflowPathRemoval_Bacteria,                   "CFU/day",   MFOutput,   MFFlux,  MFBoundary))   == CMfailed) ||	
                ((_MDflowPathRemovalMC_BacteriaID      = MFVarGetID (MDVarflowPathRemovalMC_Bacteria,                   "CFU/day",   MFOutput,   MFFlux,  MFBoundary))   == CMfailed) ||	
                ((_MDflowPathRemovalMixing_BacteriaID  = MFVarGetID (MDVarflowPathRemovalMixing_Bacteria,             "CFU/day",   MFOutput,   MFFlux,  MFBoundary))   == CMfailed) ||	
                ((_MDHTSRemoval_BacteriaID             = MFVarGetID (MDVarHTSRemoval_Bacteria,                        "CFU/day",   MFOutput,   MFFlux,  MFBoundary))   == CMfailed) ||	
                //                ((_MDOutTEHTSID                        = MFVarGetID (MDVarTEHTS,                                      "-",   MFOutput,  MFState,  MFBoundary))   == CMfailed) ||
                
      (MFModelAddFunction (_MDBacteriav2) == CMfailed)) return (CMfailed);
        
      MFDefLeaving ("Bacteriav2 Routing");
	
	return (_MDFlux_BacteriaID); 
}


