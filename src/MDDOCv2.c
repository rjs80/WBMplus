/******************************************************************************

GHAAS Water Balance/Transport Model V3.0
Global Hydrologic Archive and Analysis System
Copyright 1994-2011, UNH - CCNY/CUNY

MDDOCv2.c - DOC Modeling Class

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
static int _MDInWetlandsID             = MFUnset;
static int _MDInRiverWidthID           = MFUnset;
static int _MDInRiverOrderID           = MFUnset;
static int _MDInDOCmID                 = MFUnset;
static int _MDInDOCbID                 = MFUnset;
static int _MDInVfID                   = MFUnset;

// output

static int _MDOutLocalLoadDOCID        = MFUnset;
static int _MDFluxMixing_DOCID         = MFUnset;
static int _MDFlux_DOCID               = MFUnset; 
static int _MDOutConcMixing_DOCID      = MFUnset;
static int _MDStoreWaterMixing_DOCID   = MFUnset;
static int _MDStoreWater_DOCID         = MFUnset;
static int _MDDeltaStorageMixing_DOCID = MFUnset;
static int _MDPostConcMixing_DOCID     = MFUnset;
static int _MDOutPostConc_DOCID        = MFUnset;
static int _MDRunoffConc_DOCID         = MFUnset;

static void _MDDOCv2 (int itemID) {

    float runoffVol                     = 0.0;  // m3/s
    float waterStorage                  = 0.0;  // m3/s
    float waterStorageChange            = 0.0;  // m3/s
    float localLoad_DOC                 = 0.0;  // kg/d
    float preFluxMixing_DOC             = 0.0;  // kg/d
    float postFluxMixing_DOC            = 0.0;  // kg/d
    float storeWaterMixing_DOC          = 0.0;  // kg/d
    float postStoreWaterMixing_DOC      = 0.0;  // kg/d
    float discharge                     = 0.0;  // m3/s
    float dischargePre                  = 0.0;  // m3/s
    float waterStoragePrev              = 0.0;  // m3/s
    float DOCTotalInMixing              = 0.0;  // kg/d
    float runoffConc_DOC                = 0.0;  // mg/L
    float waterTotalVolume              = 0.0;	// m3/d
    float flowPathRemovalMixing_DOC     = 0.0;  // kg/d
    float postConcMixing_DOC            = 0.0;  // mg/L
    float massBalanceMixing_DOC         = 0.0;  // kg/d
    float preConcMixing_DOC             = 0.0;  // mg/L
    float wetlands                      = 0.0;  // proportion wetlands  
    float preFlux_DOC                   = 0.0;  // kg/d
    float storeWater_DOC                = 0.0;  // kg/d
    float DOCTotalIn                    = 0.0;  // kg/d
    float preConc_DOC                   = 0.0;  // mg/
    float Vf                            = 35.0; // m/yr
    float HL                            = 0.0;  // m/d
    float width                         = 0.0;  // m
    float removal                       = 0.0;  // proportional removal
    float totalMassRemoved_DOC          = 0.0;  // kg/d
    float flowPathRemoval_DOC           = 0.0;  // kg/d
    float postConc_DOC                  = 0.0;  // mg/L
    float postFlux_DOC                  = 0.0;  // kg/d
    float postStoreWater_DOC            = 0.0;  // kg/d
    float massBalance_DOC               = 0.0;  // kg/d
    float order                         = 0.0;  // river order
    float length                        = 0.0;  // m
    float DOC_m                         = 0.0;  // slope DOC
    float DOC_b                         = 0.0;  // intercept DOC
    
    
                discharge            = MFVarGetFloat (_MDInDischargeID,          itemID, 0.0); // m3/sec, discharge leaving the grid cell, after routing!
                dischargePre	     = MFVarGetFloat (_MDInDischarge0ID,         itemID, 0.0); // m3/sec, discharge from upstream PLUS local runoff, before routing!
                runoffVol            = MFVarGetFloat (_MDInRunoffVolumeID,       itemID, 0.0); // m3/s
		waterStorageChange   = MFVarGetFloat (_MDInRiverStorageChgID,    itemID, 0.0); // m3/s
		waterStorage         = MFVarGetFloat (_MDInRiverStorageID,       itemID, 0.0); // m3/s 
              
                preFlux_DOC          = MFVarGetFloat (_MDFlux_DOCID,             itemID, 0.0);	// kg/day RJS 091108
                preFluxMixing_DOC    = MFVarGetFloat (_MDFluxMixing_DOCID,       itemID, 0.0); // kg/day 
                storeWater_DOC       = MFVarGetFloat (_MDStoreWater_DOCID,       itemID, 0.0);	// kg/day RJS 091108
                storeWaterMixing_DOC = MFVarGetFloat (_MDStoreWaterMixing_DOCID, itemID, 0.0); // kg/day 
                wetlands             = MFVarGetFloat (_MDInWetlandsID,           itemID, 0.0); // proportion wetlands
                order                = MFVarGetFloat (_MDInRiverOrderID,         itemID, 0.0); // strahler order
                width	             = MFVarGetFloat (_MDInRiverWidthID,    	 itemID, 0.0);	// m			// moved here 031209
                length               = MFModelGetLength(itemID) / 1000;
                     
//                DOC_m                = MFVarGetFloat (_MDInDOCmID,               itemID, 0.0); // slope DOC
//                DOC_b                = MFVarGetFloat (_MDInDOCbID,               itemID, 0.0); // intercept DOC
                Vf                   = MFVarGetFloat (_MDInVfID,                 itemID, 0.0); // m/yr
                
   
                waterStoragePrev     = waterStorage - waterStorageChange;                                    // m3/sec     
                waterTotalVolume     = (discharge + waterStorage) * 86400;                                   // m3/d
                HL                   = discharge > 0.0001 ? discharge / (width * length) * 86400 : 0.0;      // m/d

                
                runoffConc_DOC       = 0.065 + 0.671 * (wetlands * 100);                                // mg/L
                localLoad_DOC        = runoffConc_DOC * runoffVol * 86400 / 1000;                       // kg/d
        
                DOCTotalIn           = localLoad_DOC + preFlux_DOC + storeWater_DOC;                    // kg/d  
                DOCTotalInMixing     = localLoad_DOC + preFluxMixing_DOC + storeWaterMixing_DOC;        // kg/d                 
               
                if (discharge > 0.0000001) {
                
                  preConc_DOC          = DOCTotalIn / waterTotalVolume * 1000;                          // mg/L          
                  preConcMixing_DOC    = DOCTotalInMixing / waterTotalVolume * 1000;                    // mg/L      

                        if (order > 2.0) {
                  
                                removal              = 1.0 - pow(2.718281828, -1.0 * (Vf / 365) / HL);  // proportional removal
                                totalMassRemoved_DOC = removal * DOCTotalIn;                            // kg/d
 
                        }
                  
                   postConcMixing_DOC   = DOCTotalInMixing / waterTotalVolume * 1000;                   // mg/L
                            
                }
                
                else {
                  
                  flowPathRemoval_DOC       = DOCTotalIn;                       // kg/d       
                  flowPathRemovalMixing_DOC = DOCTotalInMixing;                 // kg/d
            
                  postConc_DOC              = 0.0;                              // mg/L                               
                  postConcMixing_DOC        = 0.0;                              // mg/L
    
                }
   
                postConc_DOC       = discharge > 0.0000001 ? (DOCTotalIn - totalMassRemoved_DOC - flowPathRemoval_DOC) / waterTotalVolume * 1000 : 0.0;     // mg/L
                postFlux_DOC 	   = (discharge * MDConst_m3PerSecTOm3PerDay) * postConc_DOC / 1000;                                                        // kg/day
                postFluxMixing_DOC = (discharge * MDConst_m3PerSecTOm3PerDay) * postConcMixing_DOC / 1000;                                                  // kg/day

                postStoreWater_DOC       = (waterStorage * MDConst_m3PerSecTOm3PerDay) * postConc_DOC / 1000;            // kg/day
                postStoreWaterMixing_DOC = (waterStorage * MDConst_m3PerSecTOm3PerDay) * postConcMixing_DOC / 1000;	 // kg/day

                massBalance_DOC       = DOCTotalIn > 0.00001 ? (DOCTotalIn - (postFlux_DOC + postStoreWater_DOC + flowPathRemoval_DOC + totalMassRemoved_DOC)) : 0.0;    // proportion of total kg in
                massBalanceMixing_DOC = (DOCTotalInMixing - (postFluxMixing_DOC + postStoreWaterMixing_DOC + flowPathRemovalMixing_DOC)) / DOCTotalInMixing;                          // proportion of total kg in

                if (massBalance_DOC > 0.00001) {
                    printf("itemID = %d, %d-%d-%d, MB_DOC = %f\n",itemID, MFDateGetCurrentYear(), MFDateGetCurrentMonth(), MFDateGetCurrentDay(), massBalance_DOC);                  
                    printf("Q = %f, DOCTotalIn=%f, postFlux_DOC=%f, postStoreWater_DOC=%f, flowPathRemoval_DOC=%f\n",discharge,DOCTotalIn,postFlux_DOC,postStoreWater_DOC,flowPathRemoval_DOC);
                }              
                
  
  //              if ((itemID == 20) || (itemID == 9)) {
  //                  printf("*** itemID = %d, %d-%d-%d, preFlux_DOC = %f, localLoad_DOC = %f, storeWater_DOC = %f\n",itemID, MFDateGetCurrentYear(), MFDateGetCurrentMonth(), MFDateGetCurrentDay(), preFlux_DOC, localLoad_DOC, storeWater_DOC);
  //                  printf("totalMassRemoved_DOC = %f, flowPathRemoval_DOC = %f, postFlux_DOC = %f, postStoreWater_DOC = %f\n", totalMassRemoved_DOC, flowPathRemoval_DOC, postFlux_DOC, postStoreWater_DOC);
  //                  printf("dischargePre = %f, discharge = %f, waterStorage = %f, waterStoragePrev = %f\n", dischargePre, discharge, waterStorage, waterStoragePrev);
  //                  printf("MB_DOC = %f\n", massBalance_DOC);
  //              }
                
                MFVarSetFloat (_MDFluxMixing_DOCID,             itemID, postFluxMixing_DOC);
                MFVarSetFloat (_MDFlux_DOCID,                   itemID, postFlux_DOC);
                MFVarSetFloat (_MDOutLocalLoadDOCID,            itemID, localLoad_DOC);
                MFVarSetFloat (_MDOutConcMixing_DOCID,          itemID, postConcMixing_DOC);
                MFVarSetFloat (_MDOutPostConc_DOCID,            itemID, postConc_DOC);
                MFVarSetFloat (_MDStoreWaterMixing_DOCID,       itemID, postStoreWaterMixing_DOC);
                MFVarSetFloat (_MDStoreWater_DOCID,             itemID, postStoreWater_DOC);
 
                
}

int MDDOCv2Def () {

	MFDefEntering ("DOCv2 Routing");
	
        if (
            ((_MDInRiverWidthID                 = MDRiverWidthDef ())     == CMfailed) ||
            ((_MDInDischarge0ID                 = MFVarGetID (MDVarDischarge0,                               "m3/s",    MFInput,  MFState, MFBoundary))   == CMfailed) ||	
	    ((_MDInDischargeID                  = MFVarGetID (MDVarDischarge,                                "m3/s",    MFInput,  MFState, MFBoundary))   == CMfailed) ||
            ((_MDInRiverStorageID               = MFVarGetID (MDVarRiverStorage,                           "m3/day",    MFInput,  MFState, MFInitial))    == CMfailed) ||	
	    ((_MDInRiverStorageChgID            = MFVarGetID (MDVarRiverStorageChg,                        "m3/day",    MFInput,  MFState, MFBoundary))   == CMfailed) ||    
            ((_MDInRunoffVolumeID               = MFVarGetID (MDVarRunoffVolume, 		             "m3/s",    MFInput,  MFState, MFBoundary))   == CMfailed) ||
            ((_MDInWetlandsID                   = MFVarGetID (MDVarFracWetlandArea,                             "-",    MFInput,  MFState, MFBoundary))   == CMfailed) ||           // RJS 112513
            ((_MDInRiverOrderID                 = MFVarGetID (MDVarRiverOrder,                                "-",    MFInput,  MFState, MFBoundary))   == CMfailed) ||
            ((_MDRunoffConc_DOCID       	= MFVarGetID (MDVarRunoffConcDOC,               "mg/L",     MFOutput,  MFState, MFInitial))    == CMfailed) ||	
            ((_MDInDOCmID                       = MFVarGetID (MDVarDOCm,                                        "-",    MFInput,  MFState, MFBoundary))   == CMfailed) ||
            ((_MDInDOCbID                       = MFVarGetID (MDVarDOCb,                                        "-",    MFInput,  MFState, MFBoundary))   == CMfailed) ||
            ((_MDInVfID                         = MFVarGetID (MDVarVf,                                       "m/yr",    MFInput,  MFState, MFBoundary))   == CMfailed) ||
            ((_MDOutLocalLoadDOCID              = MFVarGetID (MDVarLocalLoadDOC2,               "kg/d",   MFOutput,  MFFlux,  MFBoundary))   == CMfailed) ||
            ((_MDOutConcMixing_DOCID         	= MFVarGetID (MDVarConcMixingDOC,    	        "mg/L",   MFOutput,  MFState, MFBoundary))   == CMfailed) ||	  
	    ((_MDStoreWaterMixing_DOCID      	= MFVarGetID (MDVarStoreWaterMixingDOC,       "kg/day",   MFOutput,  MFState, MFInitial))    == CMfailed) ||	
            ((_MDFluxMixing_DOCID         	= MFVarGetID (MDVarFluxMixingDOC,    	      "kg/day",   MFRoute,   MFFlux,  MFBoundary))   == CMfailed) ||	 
            ((_MDOutPostConc_DOCID         	= MFVarGetID (MDVarPostConcDOC,    	        "mg/L",   MFOutput,  MFState, MFBoundary))   == CMfailed) ||	  
	    ((_MDStoreWater_DOCID               = MFVarGetID (MDVarStoreWaterDOC,             "kg/day",   MFOutput,  MFState, MFInitial))    == CMfailed) ||	
            ((_MDFlux_DOCID                     = MFVarGetID (MDVarFluxDOC,                   "kg/day",   MFRoute,   MFFlux,  MFBoundary))   == CMfailed) ||	 
  

	(MFModelAddFunction (_MDDOCv2) == CMfailed)) return (CMfailed);
        
	MFDefLeaving ("DOCv2 Routing");
	return (_MDFlux_DOCID); 
}

