/******************************************************************************

GHAAS Water Balance/Transport Model V3.0
Global Hydrologic Archive and Analysis System
Copyright 1994-2011, UNH - CCNY/CUNY

MDSpecCond.c  - Input and Routing of Specific Conductance

shan.zuidema@unh.edu  
 * 
 * Used in conjunction with chloride.  Interprets salinazation of freshwater networks.
 * 

*******************************************************************************/
#include <stdio.h>
#include <cm.h>
#include <MF.h>
#include <MD.h>
#include <math.h>

const char* MDVarQuantiles[103] = {"Quantile0001","Quantile0010","Quantile0100","Quantile0200","Quantile0300","Quantile0400","Quantile0500","Quantile0600","Quantile0700","Quantile0800","Quantile0900","Quantile1000","Quantile1100","Quantile1200","Quantile1300","Quantile1400","Quantile1500","Quantile1600","Quantile1700","Quantile1800","Quantile1900","Quantile2000","Quantile2100","Quantile2200","Quantile2300","Quantile2400","Quantile2500","Quantile2600","Quantile2700","Quantile2800","Quantile2900","Quantile3000","Quantile3100","Quantile3200","Quantile3300","Quantile3400","Quantile3500","Quantile3600","Quantile3700","Quantile3800","Quantile3900","Quantile4000","Quantile4100","Quantile4200","Quantile4300","Quantile4400","Quantile4500","Quantile4600","Quantile4700","Quantile4800","Quantile4900","Quantile5000","Quantile5100","Quantile5200","Quantile5300","Quantile5400","Quantile5500","Quantile5600","Quantile5700","Quantile5800","Quantile5900","Quantile6000","Quantile6100","Quantile6200","Quantile6300","Quantile6400","Quantile6500","Quantile6600","Quantile6700","Quantile6800","Quantile6900","Quantile7000","Quantile7100","Quantile7200","Quantile7300","Quantile7400","Quantile7500","Quantile7600","Quantile7700","Quantile7800","Quantile7900","Quantile8000","Quantile8100","Quantile8200","Quantile8300","Quantile8400","Quantile8500","Quantile8600","Quantile8700","Quantile8800","Quantile8900","Quantile9000","Quantile9100","Quantile9200","Quantile9300","Quantile9400","Quantile9500","Quantile9600","Quantile9700","Quantile9800","Quantile9900","Quantile9990","Quantile9999"};

// input
static int _MDInTotalSurfRunoff        = MFUnset;
static int _MDInDischargeID            = MFUnset;
static int _MDInDischarge0ID           = MFUnset;
static int _MDInRunoffVolumeID         = MFUnset;
static int _MDInRunoffPoolReleaseID    = MFUnset;
static int _MDInStormRunoffTotalID    = MFUnset;
static int _MDInRiverStorageChgID      = MFUnset;
static int _MDInRiverStorageID         = MFUnset;
static int _MDInRiverWidthID           = MFUnset;
static int _MDInRiverOrderID           = MFUnset;
static int _MDInBaseFlowID             = MFUnset;

//Quantiles
static int _MDRunoffQuantileID[103] = {MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset, MFUnset};

// Loading calculation
static int _MDInSubFractionID          = MFUnset;
static int _MDInSSURGO_PERMID           = MFUnset;
static int _MDInSSURGO_AWCID            = MFUnset;
//static int _MDInLocalElevationID       = MFUnset; // unused SZ061014

// output  - SC = Specific Conductance
static int _MDOutLocalLoadSCID         = MFUnset;

static int _MDFluxMixing_SCID          = MFUnset; 
static int _MDFlux_SCID                = MFUnset; 
static int _MDOutConcMixing_SCID       = MFUnset;
static int _MDStoreWaterMixing_SCID    = MFUnset;
static int _MDStoreWater_SCID          = MFUnset;
static int _MDOutPostConc_SCID         = MFUnset;

static float _InterpQuantile (float *Q_fdc, float q_today) {
    float P[103];
    int i;
    P[0]=0.0001;P[1]=0.001;
    for (i=2;i<101;i++){
        P[i] = (i-1)/100.;
    }
    P[101]=0.999;P[102]=0.9999;
    i=0;
    while (Q_fdc[i] <= q_today){i++;}
    if (i == 0){
        return 0.0; // Baseflow 
    } else if (i == 102) {
        return 1.0; // Extreme flow
    } else {
        return (P[i]-P[i-1])/(Q_fdc[i]-Q_fdc[i-1])*(q_today-Q_fdc[i-1])+P[i-1] ;
    }
}

static void _MDSpecCond (int itemID) {
    
    float surfaceRunoff                 = 0.0;
    float runoffVol                     = 0.0;
    float stormflowVol                  = 0.0;
    float baseflowVol                   = 0.0;
    float runoffpoolreleaseVol          = 0.0;
    float waterStorage                  = 0.0;
    float waterStorageChange            = 0.0;
    float discharge                     = 0.0;
    float dischargePre                  = 0.0;
    float waterStoragePrev              = 0.0;
    float waterTotalVolume              = 0.0;	
    
    // New Variables // 
    
    float Developed                           = 0.0;
    float preFluxMixing_SC              = 0.0;
    float postFluxMixing_SC             = 0.0;
    float storeWaterMixing_SC           = 0.0;
    float postStoreWaterMixing_SC       = 0.0;
    float SCTotalInMixing               = 0.0;
    float flowPathRemovalMixing_SC      = 0.0;
    float postConcMixing_SC             = 0.0;
    float massBalanceMixing_SC          = 0.0;  
    float preConcMixing_SC              = 0.0;
    float preFlux_SC                    = 0.0;
    float storeWater_SC                 = 0.0;
    float baseflowConc_SC               = 0.0; 
    float stormflowConc_SC              = 0.0;
    float SCTotalIn                     = 0.0;
    float preConc_SC                    = 0.0;
    float Vf                            = 0.0; // Zero - no removal for ionic content - maybe slight evaporative enrichment: neglected 35.0; // m/yr
    float HL                            = 0.0;
    float width                         = 0.0;
    float removal                       = 0.0;
    float totalMassRemoved_SC           = 0.0;
    float flowPathRemoval_SC            = 0.0;
    float postConc_SC                   = 0.0;
    float postFlux_SC                   = 0.0;
    float postStoreWater_SC             = 0.0;
    float massBalance_SC                = 0.0;
    float localLoad_SC                  = 0.0;
    float order                         = 0.0;
    float length                        = 0.0;
    float SCTotalInMixing_Conc          = 0.0;
    float surfaceRunoffVol              = 0.0;
    float surfflowConc_SC               = 0.0;
    float clean_intrcp                  = 0.0;
    float SurfRO_Quantiles[103];
    float stormProbability;
    
    float ssurgo_PERM                   = 4.0;
    float ssurgo_AWC                    = 0.1;
    
    float dilutionRate                  = 0.0;
    float dilutionIntcp                 = 0.0;
    
    int i;
    
    //for (i=0;i<103;i++){
    //    SurfRO_Quantiles[i] = MFVarGetFloat(_MDRunoffQuantileID[i],    itemID, 0.0);
    //}
    
    surfaceRunoff        = MFVarGetFloat (_MDInTotalSurfRunoff,      itemID, 0.0); // mm/d
    surfaceRunoffVol     = surfaceRunoff * MFModelGetArea(itemID) / (MFModelGet_dt () * 1000.0);
    runoffVol            = MFVarGetFloat (_MDInRunoffVolumeID,       itemID, 0.0); // m3/s
    baseflowVol          = MFVarGetFloat (_MDInBaseFlowID,            itemID,0.0) * MFModelGetArea (itemID) / (MFModelGet_dt () * 1000.0); // TODO: Either need to define BaseFlowVolumeDef or calculate internally here.
    stormflowVol         = MFVarGetFloat (_MDInStormRunoffTotalID,   itemID,0.0) * MFModelGetArea (itemID) / (MFModelGet_dt () * 1000.0); 
    runoffpoolreleaseVol  = MFVarGetFloat (_MDInRunoffPoolReleaseID, itemID,0.0) * MFModelGetArea (itemID) / (MFModelGet_dt () * 1000.0);
    waterStorageChange   = MFVarGetFloat (_MDInRiverStorageChgID,    itemID, 0.0); //m3/sec
    waterStorage         = MFVarGetFloat (_MDInRiverStorageID,       itemID, 0.0); //m3/sec
    discharge            = MFVarGetFloat (_MDInDischargeID,          itemID, 0.0); // m3/sec, discharge leaving the grid cell, after routing!
    dischargePre	 = MFVarGetFloat (_MDInDischarge0ID,         itemID, 0.0); // m3/sec, discharge from upstream PLUS local runoff, before routing!
    order                = MFVarGetFloat (_MDInRiverOrderID,         itemID, 0.0);

    // New Variables //
    Developed            = MFVarGetFloat (_MDInSubFractionID,        itemID, 0.0);  // proportion developed land
    ssurgo_PERM          = MDMaximum(MFVarGetFloat (_MDInSSURGO_PERMID,        itemID, 0.0),0.0001);  // Permeability derived from SSURGO data m/d
    ssurgo_AWC           = MDMaximum(0.001 * MFVarGetFloat(_MDInSSURGO_AWCID,         itemID, 0.0),0.0001);  // AWC derived from SSURGO data (as opposed to the Reynold's method used by FrAMES?)
    
    //stormProbability  = _InterpQuantile(&SurfRO_Quantiles,surfaceRunoff);  // This is the cumulative flow probability
    // Exceedance probability is the corresponding survival function (e.g. 1 - StormProbability).
    // Per my analysis (e.g. ipynb/TestingDilutionCurves) I'm using the exceedance probability - 
  
// DEFINE IONIC CONTENT:  ic = m3*(uS/cm)
    // Assumes that ionic strength behaves linearly and conservatively.  Strictly this is only valid when
    // ionic strength is controlled by a few conservative ions.  Therefore, this is valid for salt impacted
    // streams; however, this breaks down at low ionic strength waters.  Therefore - this deep in the chloride analysis
    // we are operating with the understanding that we are only really looking at chloride impacts - not chloride itself. - SZ 6/12/2014
   
    preFluxMixing_SC     = MFVarGetFloat (_MDFluxMixing_SCID,        itemID, 0.0); // ic/day 
    storeWaterMixing_SC  = MFVarGetFloat (_MDStoreWaterMixing_SCID,  itemID, 0.0); // ic/day       
    preFlux_SC           = MFVarGetFloat (_MDFlux_SCID,              itemID, 0.0);	 // ic/day
    storeWater_SC        = MFVarGetFloat (_MDStoreWater_SCID,        itemID, 0.0);	 // ic/day
    width	         = MFVarGetFloat (_MDInRiverWidthID,    	 itemID, 0.0);	 // m			// moved here 031209
    length               = MFModelGetLength(itemID) / 1000;

    waterStoragePrev     = (waterStorage - waterStorageChange)*86400.;                              // m3/day 
    waterTotalVolume     = discharge * 86400;                       // m3/d
    HL                   = discharge > 0.0001 ? discharge / (width * length) * 86400 : 0.0;                // m/d

    /// BASEFLOW LOADING OF SPECIFIC CONDUCTANCE  ///
    ///  attributed to the groundwater and surface runoff pools - impervious runof is attributed with clean water
    /// shan.zuidema@unh.edu for details on loading estimation ///
    /// WINTER TIME runoff is not implemented yet ///
    //baseflowConc_SC       = 25.94 + 13.32 * (Developed); // Baseflow_SC relation (0.9 Quantile based) to Developed area.  Developed should already be in PER
    baseflowConc_SC = clean_intrcp + 9.602 * (Developed); // Baseflow SC relation (regression to zero storm flow) to Devloped area.
    clean_intrcp = 21.72;//25.94; //21.72;
    /// SUMMER-TIME DILIUTION ///
    // Convert to log
    ssurgo_PERM = log(ssurgo_PERM); ssurgo_AWC = log(ssurgo_AWC); Developed = log(MDMaximum(Developed,0.0001));
    // Constrain predictive parameters to the observed bounds
    ssurgo_PERM = MDMaximum(ssurgo_PERM,-0.077833); ssurgo_PERM = MDMinimum(ssurgo_PERM,1.986004);
    ssurgo_AWC = MDMaximum(ssurgo_AWC,-2.317663); ssurgo_AWC = MDMinimum(ssurgo_AWC,-1.649231);
   
    // Calculate dilution 
    dilutionRate = -0.0938*ssurgo_PERM - 0.1158 * ssurgo_AWC;
    dilutionIntcp = -0.0938*ssurgo_PERM - 0.6207 * ssurgo_AWC + 0.0117 * (Developed / 100.);
    // Constrain dilution to observed bounds
    // dilutionRate = MDMaximum(dilutionRate,-0.034711); dilutionRate = MDMinimum(dilutionRate,0.534144);
    //dilutionIntcp = MDMaximum(dilutionIntcp,0.974424); dilutionIntcp = MDMinimum(dilutionIntcp,1.252073);
    dilutionIntcp = 1.00;
    dilutionRate = 0.15;
    //stormflowConc_SC      = (dilutionIntcp + log10(MDMaximum(1.0 - stormProbability,0.0))*dilutionRate)*baseflowConc_SC; // SC dependence on catchment soil/development attributes and storm flow probability
    //stormflowConc_SC      = stormflowConc_SC < 21.72 ? 21.72 : stormflowConc_SC;
   surfflowConc_SC = MDMaximum(0.64 * baseflowConc_SC,clean_intrcp); // 0.8244
     // Average fraction of baseflow SC exhibited during storm events (at daily flow exceedance of 1% P1)
    // TODO: if this can be determined to be the result of some land-scape metric, then this could at least be related to flow probability?
   stormflowConc_SC = MDMaximum(0.20 * baseflowConc_SC,0.8*clean_intrcp);
    // Calculate input loadings
    localLoad_SC         = (baseflowConc_SC * baseflowVol + surfflowConc_SC * runoffpoolreleaseVol + stormflowVol*stormflowConc_SC ) * 86400 ; // ic/day
     
    SCTotalInMixing      = localLoad_SC + preFluxMixing_SC  + storeWaterMixing_SC;         // ic/day  
    SCTotalIn            = localLoad_SC + preFlux_SC  + storeWater_SC;          // ic/day  
    
    // Calculate mixed concentration
    SCTotalInMixing_Conc = discharge > 0.0000001 ? SCTotalInMixing / waterTotalVolume : SCTotalInMixing; // uS/cm
   
    if (discharge > 0.0000001) {
        // Calculate pre-tranformation concentration - for cells with accumulated discharge
        preConcMixing_SC     = SCTotalInMixing / waterTotalVolume;      // uS/cm  
        preConc_SC     = SCTotalIn  / waterTotalVolume ; // uS/cm
        // Initializes post tranformation concentration ?
        postConcMixing_SC    = SCTotalInMixing / waterTotalVolume; // uS/cm

        if (order > 2.0) {
        // Calculate tranformation only for cells with any accumulated flow
                //removal = 1.0 - pow(2.718281828, -1.0 * (Vf / 365) / HL);
                removal = 0.0;
                totalMassRemoved_SC  = removal * SCTotalIn;
        } 
    } else {
        // Force concentrations to local values for minimal accumulation / dessicated cells
        flowPathRemoval_SC          = SCTotalIn;                               
        flowPathRemovalMixing_SC    = SCTotalInMixing;
        postConc_SC                 = baseflowConc_SC; //  
        postConcMixing_SC           = baseflowConc_SC; //
    }
    // Calculates first-order uptake like a reactive solute.  (For comparison purposes only)  // Sets to local SC for minimal accumulation/dessicated cells
    postConc_SC    = discharge > 0.0000001 ? (SCTotalIn - totalMassRemoved_SC - flowPathRemoval_SC) / waterTotalVolume : baseflowConc_SC;   //uS/cm             
    postFlux_SC        = (discharge * MDConst_m3PerSecTOm3PerDay) * postConc_SC ;         // ic/day
    postStoreWater_SC  = (waterStorage ) * postConc_SC ;      // ic/day
    
    //Calculates the fluxes/storage for the mixing (appropriate for speccond) case
    postFluxMixing_SC  = (discharge * MDConst_m3PerSecTOm3PerDay) * postConcMixing_SC ;         // ic/day
    postStoreWaterMixing_SC  = (waterStorage ) * postConcMixing_SC ;      // ic/day
    
    // Calculates the mass balances
    massBalance_SC  = SCTotalIn > 0.00001 ? (SCTotalIn - (postFlux_SC + postStoreWater_SC + flowPathRemoval_SC + totalMassRemoved_SC)) / SCTotalIn : 0.0;
    massBalanceMixing_SC  = (SCTotalInMixing - (postFluxMixing_SC + postStoreWaterMixing_SC + flowPathRemovalMixing_SC)) / SCTotalInMixing;
      
    // Prints mass balance errors
    
    if (MFDateGetCurrentYear() > 0){
        if ( (massBalanceMixing_SC > 0.001)) {
           printf("itemID = %d, %d-%d-%d, MB_SC = %f\n",itemID, MFDateGetCurrentYear(), MFDateGetCurrentMonth(), MFDateGetCurrentDay(), massBalanceMixing_SC);
           printf("\tTotalIn: %f, Local: %f , UpFlux: %f, InStore: %f\n",SCTotalInMixing,localLoad_SC,preFluxMixing_SC,storeWaterMixing_SC);
           printf("\tDownFlux: %f, OutStore: %f , flowPathRemove: %f\n",postFluxMixing_SC,postStoreWaterMixing_SC,flowPathRemovalMixing_SC);
        }
        
        //if ((itemID == 8310)) {
            //printf("itemID = %d, order %f, Store_SC %f, Store %f Store2 %f StoreDelta %f StorePrev %f StoreSCconc %f \n ",itemID, order, storeWater_SC, waterStorage*86400.,(discharge - dischargePre)*86400.,waterStorageChange*86400.,waterStoragePrev, storeWater_SC / waterStoragePrev ) ;
        //}
        //if (postConcMixing_SC < 21.70){
            //printf("itemID = %d, order %f, Store_SC %f, Store %f Store2 %f StoreDelta %f StorePrev %f StoreSCconc %f \n ",itemID, order, storeWater_SC, waterStorage*86400.,(discharge - dischargePre)*86400.,waterStorageChange*86400.,waterStoragePrev, storeWater_SC / waterStoragePrev ) ;
          //  printf("itemID %d ro %f bf %f srf %f str %f bfSC %f srfSC %f strSC %f \n",itemID,runoffVol,baseflowVol,runoffpoolreleaseVol,stormflowVol,baseflowConc_SC,stormflowConc_SC,clean_intrcp);
        //}
    }
    
     
    // Set Output
    MFVarSetFloat (_MDFluxMixing_SCID,              itemID, postFluxMixing_SC);
    MFVarSetFloat (_MDFlux_SCID,                    itemID, postFlux_SC);

    MFVarSetFloat (_MDOutLocalLoadSCID,             itemID, localLoad_SC);

    MFVarSetFloat (_MDOutConcMixing_SCID,           itemID, postConcMixing_SC);
    MFVarSetFloat (_MDOutPostConc_SCID,             itemID, postConc_SC);

    MFVarSetFloat (_MDStoreWaterMixing_SCID,        itemID, postStoreWaterMixing_SC);
    MFVarSetFloat (_MDStoreWater_SCID,              itemID, postStoreWater_SC);
}

enum {MDcalculate, MDinput, MDinput2, MDnone};

int MDSpecCondDef () {
    int i;
    int  optID = MFUnset;											//SZ 08212014
    const char *optStr, *optName = MDOptSpecConductance;								//SZ 08212014
    const char *options [] = { MDCalculateStr, MDInputStr, MDInput2Str, MDNoneStr, (char *) NULL };		//SZ 08212014

    MFDefEntering ("Specific Conductance Routing");

    if ((optStr = MFOptionGet (optName)) != (char *) NULL) optID = CMoptLookup (options, optStr, true);  //SZ 08212014

    switch (optID) {	//SZ 08212014

        case MDcalculate:
            
            //for(i=0;i<103;i++) {
              //  if ((_MDRunoffQuantileID[i] = MFVarGetID(MDVarQuantiles[i],                              "mm",   MFInput,    MFState,    MFBoundary)) == CMfailed) return (CMfailed);                
            //}

            // May add case options in future.  I 

            if (
                ((_MDInRiverWidthID                 = MDDischargeDef ())     == CMfailed) ||
                ((_MDInTotalSurfRunoff              = MFVarGetID (MDVarTotalSurfRunoff,                           "mm",    MFOutput, MFFlux, MFBoundary))  == CMfailed) ||
                ((_MDInDischarge0ID                 = MFVarGetID (MDVarDischarge0,                               "m3/s",    MFInput,  MFState, MFBoundary))   == CMfailed) ||	
                ((_MDInDischargeID                  = MFVarGetID (MDVarDischarge,                                "m3/s",    MFInput,  MFState, MFBoundary))   == CMfailed) ||
                ((_MDInRiverStorageID               = MFVarGetID (MDVarRiverStorage,                           "m3/day",    MFInput,  MFState, MFInitial))    == CMfailed) ||	
                ((_MDInRiverStorageChgID            = MFVarGetID (MDVarRiverStorageChg,                        "m3/day",    MFInput,  MFState, MFBoundary))   == CMfailed) ||    
                ((_MDInRunoffVolumeID               = MFVarGetID (MDVarRunoffVolume, 		             "m3/s",    MFInput,  MFState, MFBoundary))   == CMfailed) ||
                ((_MDInStormRunoffTotalID           = MFVarGetID (MDVarStormRunoffTotal,                          "mm",   MFOutput, MFFlux, MFBoundary))  == CMfailed) ||
                ((_MDInRunoffPoolReleaseID          = MFVarGetID (MDVarRunoffPoolRelease,                        "mm",  MFOutput, MFFlux, MFBoundary))  == CMfailed) ||
                ((_MDInBaseFlowID                   = MFVarGetID (MDVarBaseFlow,                                  "mm",    MFInput,  MFFlux,MFBoundary))    == CMfailed) ||
                ((_MDInRiverOrderID                 = MFVarGetID (MDVarRiverOrder,                                "-",    MFInput,  MFState, MFBoundary))   == CMfailed) ||
                ((_MDInSubFractionID                = MFVarGetID (MDVarLandUseSpatialSub,                         "-",    MFInput,  MFState, MFBoundary))   == CMfailed) ||
                ((_MDInSSURGO_AWCID                 = MFVarGetID (MDVarSoilAvailWaterCapInput,                   "mm",    MFInput,  MFState, MFBoundary))   == CMfailed) || // Note that I'm redefining the AWC input to SSURGO data...
                ((_MDInSSURGO_PERMID                = MFVarGetID (MDVarSoilPermeability,                         "m/d",    MFInput,  MFState, MFBoundary))   == CMfailed) ||  
                ((_MDOutLocalLoadSCID               = MFVarGetID (MDVarLocalLoadSC,               "ic/d",   MFOutput,  MFFlux,  MFBoundary))   == CMfailed) ||
                ((_MDOutConcMixing_SCID             = MFVarGetID (MDVarConcMixingSC,    	      "uS/cm",   MFOutput,  MFFlux, MFBoundary))   == CMfailed) ||	               
                ((_MDStoreWaterMixing_SCID          = MFVarGetID (MDVarStoreWaterMixingSC,        "ic/day",   MFOutput,  MFState, MFInitial))    == CMfailed) ||	
                ((_MDFluxMixing_SCID                = MFVarGetID (MDVarFluxMixingSC,    	      "ic/day",   MFRoute,   MFFlux,  MFBoundary))   == CMfailed) ||	  

                ((_MDOutPostConc_SCID         	= MFVarGetID (MDVarPostConcSC,    	      "uS/cm",   MFOutput,  MFFlux, MFBoundary))   == CMfailed) ||	               
                ((_MDStoreWater_SCID                = MFVarGetID (MDVarStoreWaterSC,              "ic/day",   MFOutput,  MFState, MFInitial))    == CMfailed) ||	
                ((_MDFlux_SCID                      = MFVarGetID (MDVarFluxSC,                    "ic/day",   MFRoute,   MFFlux,  MFBoundary))   == CMfailed) ||	  


            (MFModelAddFunction (_MDSpecCond) == CMfailed)) return (CMfailed);
            break;

        default: 
            MFOptionMessage (optName, optStr, options); return (CMfailed);
    }

    MFDefLeaving ("Specific Conductance Routing");
    return (_MDFluxMixing_SCID);
}
