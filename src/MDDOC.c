/******************************************************************************

GHAAS Water Balance/Transport Model V3.0
Global Hydrologic Archive and Analysis System
Copyright 1994-2011, UNH - CCNY/CUNY

MDDOC.c  - Input and Routing of DOC 

rob.stewart@unh.edu  

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
static int _MDInRiverStorageChgID      = MFUnset;
static int _MDInRiverStorageID         = MFUnset;
static int _MDInLocalLoad_DOCID        = MFUnset;
static int _MDInLocalLoad_LTRID        = MFUnset;
static int _MDFlux_DIN_denitID         = MFUnset;
static int _MDInTotalEvaporationID     = MFUnset;

// output

static int _MDFluxMixing_DOCID         = MFUnset;
static int _MDFluxMixing_LTRID         = MFUnset;
static int _MDOutConcMixing_DOCID      = MFUnset;
static int _MDOutConcMixing_LTRID      = MFUnset;
static int _MDStoreWaterMixing_DOCID   = MFUnset;
static int _MDStoreWaterMixing_LTRID   = MFUnset;
static int _MDDeltaStorageMixing_DOCID = MFUnset;
static int _MDDeltaStorageMixing_LTRID = MFUnset;
static int _MDPostConcMixing_DOCID     = MFUnset;
static int _MDPostConcMixing_LTRID     = MFUnset;

static void _MDDOC (int itemID) {

    float runoffVol                     = 0.0;
    float waterStorage                  = 0.0;
    float waterStorageChange            = 0.0;
    float localLoad_DOC                 = 0.0;
    float localLoad_LTR                 = 0.0;
    float preFluxMixing_DOC             = 0.0;
    float preFluxMixing_LTR             = 0.0;
    float postFluxMixing_DOC            = 0.0;
    float postFluxMixing_LTR            = 0.0;
    float storeWaterMixing_DOC          = 0.0;
    float storeWaterMixing_LTR          = 0.0;
    float postStoreWaterMixing_DOC      = 0.0;
    float postStoreWaterMixing_LTR      = 0.0;
    float discharge                     = 0.0;
    float dischargePre                  = 0.0;
    float waterStoragePrev              = 0.0;
    float DOCTotalInMixing              = 0.0;
    float LTRTotalInMixing              = 0.0;
    float runoffConc_DOC                = 0.0;
    float runoffConc_LTR                = 0.0;
    float waterTotalVolume              = 0.0;	
    float DIN                           = 0.0;
    float flowPathRemovalMixing_DOC     = 0.0;
    float flowPathRemovalMixing_LTR     = 0.0;
    float postConcMixing_DOC            = 0.0;
    float postConcMixing_LTR            = 0.0;
    float massBalanceMixing_DOC         = 0.0;
    float massBalanceMixing_LTR         = 0.0;  
    float preConcMixing_DOC             = 0.0;
    float preConcMixing_LTR             = 0.0;
    float totalEvap                     = 0.0;

    
                runoffVol            = MFVarGetFloat (_MDInRunoffVolumeID,       itemID, 0.0); // m3/s
		waterStorageChange   = MFVarGetFloat (_MDInRiverStorageChgID,    itemID, 0.0);
		waterStorage         = MFVarGetFloat (_MDInRiverStorageID,       itemID, 0.0);
                localLoad_DOC  	     = MFVarGetFloat (_MDInLocalLoad_DOCID,      itemID, 0.0); // kg/day
                localLoad_LTR        = MFVarGetFloat (_MDInLocalLoad_LTRID,      itemID, 0.0); // kg/day
                preFluxMixing_DOC    = MFVarGetFloat (_MDFluxMixing_DOCID,       itemID, 0.0); // kg/day 
                preFluxMixing_LTR    = MFVarGetFloat (_MDFluxMixing_LTRID,       itemID, 0.0); // kg/day 
                storeWaterMixing_DOC = MFVarGetFloat (_MDStoreWaterMixing_DOCID, itemID, 0.0); // kg/day 
                storeWaterMixing_LTR = MFVarGetFloat (_MDStoreWaterMixing_LTRID, itemID, 0.0); // kg/day 
                discharge            = MFVarGetFloat (_MDInDischargeID,          itemID, 0.0); // m3/sec, discharge leaving the grid cell, after routing!
                dischargePre	     = MFVarGetFloat (_MDInDischarge0ID,         itemID, 0.0); // m3/sec, discharge from upstream PLUS local runoff, before routing!
                DIN                  = MFVarGetFloat (_MDFlux_DIN_denitID,       itemID, 0.0); // to start DIN model
                totalEvap            = MFVarGetFloat (_MDInTotalEvaporationID,   itemID, 0.0); // m3/day
                
                waterStoragePrev     = waterStorage - waterStorageChange;                               // m3/sec     
                waterTotalVolume     = (waterStoragePrev + dischargePre) * 86400;                       // m3/d
                     
                
                DOCTotalInMixing     = localLoad_DOC + preFluxMixing_DOC + storeWaterMixing_DOC;        // kg/day  
                LTRTotalInMixing     = localLoad_LTR + preFluxMixing_LTR + storeWaterMixing_LTR;        // kg/day
                runoffConc_DOC	     = runoffVol > 0.0 ? (localLoad_DOC * 1000000) / (runoffVol * 86400 * 1000) : 0.0;	
                runoffConc_LTR       = runoffVol > 0.0 ? (localLoad_LTR * 1000000) / (runoffVol * 86400 * 1000) : 0.0;

                preConcMixing_DOC    = DOCTotalInMixing / waterTotalVolume * 1000;
                preConcMixing_LTR    = LTRTotalInMixing / waterTotalVolume * 1000;
                              
                if (discharge <= 0.000001) {
                    flowPathRemovalMixing_DOC = DOCTotalInMixing;
                    flowPathRemovalMixing_LTR = LTRTotalInMixing;
                }
                
                postConcMixing_DOC      = (DOCTotalInMixing - flowPathRemovalMixing_DOC) / (waterTotalVolume - totalEvap) * 1000;				    // mg/L
                postConcMixing_LTR      = (LTRTotalInMixing - flowPathRemovalMixing_LTR) / (waterTotalVolume - totalEvap) * 1000;
                
                postFluxMixing_DOC 	 = (discharge * MDConst_m3PerSecTOm3PerDay) * postConcMixing_DOC / 1000;		// kg/day
                postFluxMixing_LTR       = (discharge * MDConst_m3PerSecTOm3PerDay) * postConcMixing_LTR / 1000;         // kg/day
                postStoreWaterMixing_DOC = (waterStorage * MDConst_m3PerSecTOm3PerDay) * postConcMixing_DOC / 1000;	// kg/day
                postStoreWaterMixing_LTR = (waterStorage * MDConst_m3PerSecTOm3PerDay) * postConcMixing_LTR / 1000;      // kg/day

                massBalanceMixing_DOC = ((localLoad_DOC + preFluxMixing_DOC + storeWaterMixing_DOC) - (postFluxMixing_DOC + postStoreWaterMixing_DOC + flowPathRemovalMixing_DOC)) / (localLoad_DOC + storeWaterMixing_DOC + preFluxMixing_DOC);
                massBalanceMixing_LTR = ((localLoad_LTR + preFluxMixing_LTR + storeWaterMixing_LTR) - (postFluxMixing_LTR + postStoreWaterMixing_LTR + flowPathRemovalMixing_LTR)) / (localLoad_LTR + storeWaterMixing_LTR + preFluxMixing_LTR);
 
 //               if ((itemID == 2890) || (itemID == 1599)) {
                if (((massBalanceMixing_DOC > 0.0003) || (massBalanceMixing_LTR > 0.0003)) && (localLoad_DOC + storeWaterMixing_DOC + preFluxMixing_DOC > 0.000001)) {
                    printf("******itemID = %d, y = %d, m = %d, d = %d, massBalance_DOC = %f, massBalance_LTR = %f, preConc_DOC = %f, preConc_LTR = %f\n", itemID, MFDateGetCurrentYear(), MFDateGetCurrentMonth(), MFDateGetCurrentDay(), massBalanceMixing_DOC, massBalanceMixing_LTR, preConcMixing_DOC, preConcMixing_LTR);
                    printf("Evap = %f, dischargePre = %f, discharge = %f, runoffVol = %f, waterStoragePrev = %f, waterStorage = %f, waterTotalVol = %f\n", totalEvap, dischargePre, discharge, runoffVol, waterStoragePrev, waterStorage, waterTotalVolume);
                    printf("totalIn_DOC = %f, totalIn_LTR = %f, localLoad_DOC = %f, localLoad_LTR = %f, preFlux_DOC = %f, preFlux_LTR = %f, storeWater_DOC = %f, storeWater_LTR = %f\n", DOCTotalInMixing, LTRTotalInMixing, localLoad_DOC, localLoad_LTR, preFluxMixing_DOC, preFluxMixing_LTR, storeWaterMixing_DOC, storeWaterMixing_LTR);
                    printf("postConc_DOC = %f, postConc_LTR = %f, postFlux_DOC = %f, postFlux_LTR = %f, postStore_DOC = %f, postStore_LTR = %f, flowPathRem_DOC = %f, flowPathRem_LTR = %f\n", postConcMixing_DOC, postConcMixing_LTR, postFluxMixing_DOC, postFluxMixing_LTR, postStoreWaterMixing_DOC, postStoreWaterMixing_LTR, flowPathRemovalMixing_DOC, flowPathRemovalMixing_LTR);
                }
  //              }
               
                MFVarSetFloat (_MDFluxMixing_DOCID,             itemID, postFluxMixing_DOC);
                MFVarSetFloat (_MDFluxMixing_LTRID,             itemID, postFluxMixing_LTR);
                MFVarSetFloat (_MDFluxMixing_DOCID,             itemID, postConcMixing_DOC);
                MFVarSetFloat (_MDFluxMixing_LTRID,             itemID, postConcMixing_LTR);
                MFVarSetFloat (_MDStoreWaterMixing_DOCID,       itemID, postStoreWaterMixing_DOC);
                MFVarSetFloat (_MDStoreWaterMixing_LTRID,       itemID, postStoreWaterMixing_LTR);

   
}

int MDDOCDef () {

	MFDefEntering ("DOC Routing");
	
        if (
	    ((_MDFlux_DIN_denitID               = MDDINDef ()) == CMfailed) ||	
            ((_MDInLocalLoad_DOCID              = MFVarGetID (MDVarLocalLoadDOC,                             "kg/d",    MFInput,  MFFlux,  MFBoundary))   == CMfailed) ||
            ((_MDInLocalLoad_LTRID              = MFVarGetID (MDVarLocalLoadLTR,                             "kg/d",    MFInput,  MFFlux,  MFBoundary))   == CMfailed) ||
            ((_MDInDischarge0ID                 = MFVarGetID (MDVarDischarge0,                               "m3/s",    MFInput,  MFState, MFBoundary))   == CMfailed) ||	
	    ((_MDInDischargeID                  = MFVarGetID (MDVarDischarge,                                "m3/s",    MFInput,  MFState, MFBoundary))   == CMfailed) ||
            ((_MDInRiverStorageID               = MFVarGetID (MDVarRiverStorage,                           "m3/day",    MFInput,  MFState, MFInitial))    == CMfailed) ||	
	    ((_MDInRiverStorageChgID            = MFVarGetID (MDVarRiverStorageChg,                        "m3/day",    MFInput,  MFState, MFBoundary))   == CMfailed) ||    
            ((_MDInRunoffVolumeID               = MFVarGetID (MDVarRunoffVolume, 		             "m3/s",    MFInput,  MFState, MFBoundary))   == CMfailed) ||
            ((_MDInTotalEvaporationID           = MFVarGetID (MDVarTotalEvaporation,                           "m3",    MFInput,  MFFlux,  MFBoundary)) == CMfailed) ||
   
            ((_MDOutConcMixing_DOCID         	= MFVarGetID (MDVarConcMixingDOC,    	        "mg/L",   MFOutput,  MFState, MFBoundary))   == CMfailed) ||	  
            ((_MDOutConcMixing_LTRID         	= MFVarGetID (MDVarConcMixingLTR,    	        "mg/L",   MFOutput,  MFState, MFBoundary))   == CMfailed) ||	               
	    ((_MDStoreWaterMixing_DOCID      	= MFVarGetID (MDVarStoreWaterMixingDOC,       "kg/day",   MFOutput,  MFState, MFInitial))    == CMfailed) ||	
	    ((_MDStoreWaterMixing_LTRID      	= MFVarGetID (MDVarStoreWaterMixingLTR,       "kg/day",   MFOutput,  MFState, MFInitial))    == CMfailed) ||	
            ((_MDFluxMixing_DOCID         	= MFVarGetID (MDVarFluxMixingDOC,    	      "kg/day",   MFRoute,   MFFlux,  MFBoundary))   == CMfailed) ||	 
            ((_MDFluxMixing_LTRID         	= MFVarGetID (MDVarFluxMixingLTR,    	      "kg/day",   MFRoute,   MFFlux,  MFBoundary))   == CMfailed) ||	  


	(MFModelAddFunction (_MDDOC) == CMfailed)) return (CMfailed);
        
	MFDefLeaving ("DOC Routing");
	return (_MDFluxMixing_DOCID); 
}
