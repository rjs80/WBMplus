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

// input
static int _MDInDischargeID            = MFUnset;
static int _MDInDischarge0ID           = MFUnset;
static int _MDInRunoffVolumeID         = MFUnset;
static int _MDInRunoffID               = MFUnset;
static int _MDInRunoffPoolReleaseID    = MFUnset;
static int _MDInStormRunoffTotalID    = MFUnset;
static int _MDInRiverStorageChgID      = MFUnset;
static int _MDInRiverStorageID         = MFUnset;
static int _MDInRiverWidthID           = MFUnset;
static int _MDInRiverOrderID           = MFUnset;
static int _MDInBaseFlowID             = MFUnset;
static int _MDInRunoffCorrID          = MFUnset;

static int _MDInWTempRiverID           =MFUnset;

// Loading calculation
static int _MDInSubFractionID          = MFUnset;
//static int _MDInLocalElevationID       = MFUnset; // unused SZ061014

// output  - SC = Specific Conductance
static int _MDOutLocalLoadSCID         = MFUnset;

static int _MDFluxMixing_SCID          = MFUnset; 
static int _MDFlux_SCID                = MFUnset; 
static int _MDOutConcMixing_SCID       = MFUnset;
static int _MDStoreWaterMixing_SCID    = MFUnset;
static int _MDStoreWater_SCID          = MFUnset;
static int _MDOutPostConc_SCID         = MFUnset;

static void _MDSpecCond (int itemID) {

    float runoffVol                     = 0.0;
    float stormflowVol                  = 0.0;
    float baseflowVol                   = 0.0;
    float runoffpoolrelaseVol            = 0.0;
    float waterStorage                  = 0.0;
    float waterStorageChange            = 0.0;
    float discharge                     = 0.0;
    float dischargePre                  = 0.0;
    float waterStoragePrev              = 0.0;
    float waterTotalVolume              = 0.0;	
    
    // New Variables // 
    
    float dev                           = 0.0;
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
    
    runoffVol            = MFVarGetFloat (_MDInRunoffVolumeID,       itemID, 0.0); // m3/s
    baseflowVol          = MFVarGetFloat (_MDInBaseFlowID,            itemID,0.0) * MFModelGetArea (itemID) / (MFModelGet_dt () * 1000.0); // TODO: Either need to define BaseFlowVolumeDef or calculate internally here.
    stormflowVol         = MFVarGetFloat (_MDInStormRunoffTotalID,   itemID,0.0) * MFModelGetArea (itemID) / (MFModelGet_dt () * 1000.0); 
    runoffpoolrelaseVol  = MFVarGetFloat (_MDInRunoffPoolReleaseID, itemID,0.0) * MFModelGetArea (itemID) / (MFModelGet_dt () * 1000.0);
    waterStorageChange   = MFVarGetFloat (_MDInRiverStorageChgID,    itemID, 0.0);
    waterStorage         = MFVarGetFloat (_MDInRiverStorageID,       itemID, 0.0);
    discharge            = MFVarGetFloat (_MDInDischargeID,          itemID, 0.0); // m3/sec, discharge leaving the grid cell, after routing!
    dischargePre	 = MFVarGetFloat (_MDInDischarge0ID,         itemID, 0.0); // m3/sec, discharge from upstream PLUS local runoff, before routing!
    order                = MFVarGetFloat (_MDInRiverOrderID,         itemID, 0.0);

    // New Variables //

    dev                  = MFVarGetFloat (_MDInSubFractionID,        itemID, 0.0);  // proportion developed land

// DEFINE IONIC CONTENT:  ic = m3*(uS/cm)
    // Assumes that ionic strength behaves linearly and conservatively.  Strictly this is only valid when
    // ionic strength is controlled by a few conservative ions.  Therefore, this is valid for salt impacted
    // streams; however, this breaks down at low ionic strength waters.  Therefore - this deep in the chloride analysis
    // we are operating with the understanding that we are only really looking at chloride impacts - not chloride itself. - SZ 6/12/2014
   
    preFluxMixing_SC     = MFVarGetFloat (_MDFluxMixing_SCID,        itemID, 0.0); // ic/day 
    storeWaterMixing_SC  = MFVarGetFloat (_MDStoreWaterMixing_SCID,  itemID, 0.0); // ic/day       
    preFlux_SC           = MFVarGetFloat (_MDFlux_SCID,              itemID, 0.0);	 // ic/day RJS 091108
    storeWater_SC        = MFVarGetFloat (_MDStoreWater_SCID,        itemID, 0.0);	 // ic/day RJS 091108
    width	         = MFVarGetFloat (_MDInRiverWidthID,    	 itemID, 0.0);	 // m			// moved here 031209
    length               = MFModelGetLength(itemID) / 1000;



    waterStoragePrev     = waterStorage - waterStorageChange;                              // m3/sec     
    waterTotalVolume     = discharge * 86400;                       // m3/d
    HL                   = discharge > 0.0001 ? discharge / (width * length) * 86400 : 0.0;                // m/d

    /// BASEFLOW LOADING OF SPECIFIC CONDUCTANCE  ///
    ///  attributed to the groundwater and surface runoff pools - impervious runof is attributed with clean water
    /// shan.zuidema@unh.edu for details on loading estimation ///
    /// WINTER TIME runoff is not implemented yet ///
    baseflowConc_SC       = 53.18 + 9.09 * (dev); // Baseflow_SC relation to Developed area.  Developed should already be in PER
    float stormfrac_SC = 0.488 ; // Average fraction of baseflow SC exhibited during storm events (at daily flow exceedance of 1% P1)
    // TODO: if this can be determined to be the result of some land-scape metric, then this could at least be related to flow probability?
    
    // Calculate input loadings
    localLoad_SC         = (baseflowConc_SC * baseflowVol + stormfrac_SC * baseflowConc_SC * (stormflowVol + runoffpoolrelaseVol) ) * 86400 ; // ic/day
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
    postStoreWater_SC  = (waterStorage * MDConst_m3PerSecTOm3PerDay) * postConc_SC ;      // ic/day
    
    //Calculates the fluxes/storage for the mixing (appropriate for speccond) case
    postFluxMixing_SC  = (discharge * MDConst_m3PerSecTOm3PerDay) * postConcMixing_SC ;         // ic/day
    postStoreWaterMixing_SC  = (waterStorage * MDConst_m3PerSecTOm3PerDay) * postConcMixing_SC ;      // ic/day
    
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

int MDSpecCondDef () {
	MFDefEntering ("Specific Conductance Routing");
	
        if (
	    ((_MDInWTempRiverID                 = MDWTempRiverRouteDef ()) == CMfailed) ||	
            ((_MDInRiverWidthID                 = MDRiverWidthDef ())     == CMfailed) ||
//          ((_MDInLitterFall_POCID             = MDLitterFallDef ()) == CMfailed) ||
//          ((_MDInLocalLoad_DOCID              = MFVarGetID (MDVarLocalLoadDOC,                             "kg/d",    MFInput,  MFFlux,  MFBoundary))   == CMfailed) ||
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
            ((_MDOutLocalLoadSCID               = MFVarGetID (MDVarLocalLoadSC,               "ic/d",   MFOutput,  MFFlux,  MFBoundary))   == CMfailed) ||
            ((_MDOutConcMixing_SCID         	= MFVarGetID (MDVarConcMixingSC,    	      "uS/cm",   MFOutput,  MFState, MFBoundary))   == CMfailed) ||	               
	    ((_MDStoreWaterMixing_SCID      	= MFVarGetID (MDVarStoreWaterMixingSC,        "ic/day",   MFOutput,  MFState, MFInitial))    == CMfailed) ||	
            ((_MDFluxMixing_SCID         	= MFVarGetID (MDVarFluxMixingSC,    	      "ic/day",   MFRoute,   MFFlux,  MFBoundary))   == CMfailed) ||	  

            ((_MDOutPostConc_SCID         	= MFVarGetID (MDVarPostConcSC,    	      "uS/cm",   MFOutput,  MFState, MFBoundary))   == CMfailed) ||	               
	    ((_MDStoreWater_SCID                = MFVarGetID (MDVarStoreWaterSC,              "ic/day",   MFOutput,  MFState, MFInitial))    == CMfailed) ||	
            ((_MDFlux_SCID                      = MFVarGetID (MDVarFluxSC,                    "ic/day",   MFRoute,   MFFlux,  MFBoundary))   == CMfailed) ||	  
    

	(MFModelAddFunction (_MDSpecCond) == CMfailed)) return (CMfailed);
        
	MFDefLeaving ("Specific Conductance Routing");
	return (_MDFluxMixing_SCID);
}
