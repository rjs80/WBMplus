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
#include <stdlib.h>

// input

static int _MDInDischargeID            = MFUnset;
static int _MDInDischarge0ID           = MFUnset;
static int _MDInRunoffVolumeID         = MFUnset;
static int _MDInRunoffID               = MFUnset;
static int _MDInRiverStorageChgID      = MFUnset;
static int _MDInRiverStorageID         = MFUnset;
static int _MDInLocalLoad_DOCID        = MFUnset;
static int _MDInLitterFall_POCID       = MFUnset;
static int _MDFlux_DIN_denitID         = MFUnset;
static int _MDInTotalEvaporationID     = MFUnset;
static int _MDInWetlandsID             = MFUnset;
static int _MDInImpFractionID          = MFUnset;
static int _MDInHumanLandUseID         = MFUnset;
static int _MDInRiverWidthID           = MFUnset;
static int _MDInRiverOrderID           = MFUnset;
static int _MDInDOCmID                 = MFUnset;
static int _MDInDOCbID                 = MFUnset;
static int _MDInClmID                  = MFUnset;
static int _MDInClbID                  = MFUnset;
static int _MDInDINmID                 = MFUnset;
static int _MDInDINbID                 = MFUnset;
static int _MDInVfID                   = MFUnset;
static int _MDInError_DOCID            = MFUnset;
static int _MDInError_DOC2ID           = MFUnset;
static int _MDInError_NO3ID            = MFUnset;
static int _MDInError_ClID             = MFUnset;
static int _MDInError_HPOA1ID          = MFUnset;
static int _MDInError_HPOA2ID          = MFUnset;
static int _MDInError_HPOA3ID          = MFUnset;
static int _MDInError_nHPOAID          = MFUnset;
static int _MDInError_TPIAID           = MFUnset;
static int _MDInError_HPIID            = MFUnset;
static int _MDInLoadAdjustID           = MFUnset;

// output

static int _MDOutLocalLoadDOCID        = MFUnset;
static int _MDOutLocalLoadDOC2ID       = MFUnset;
static int _MDOutLocalLoadDINID        = MFUnset;
static int _MDOutLocalLoadClID         = MFUnset;
static int _MDOutLocalLoadHPOA1ID      = MFUnset;
static int _MDOutLocalLoadHPOA2ID      = MFUnset;
static int _MDOutLocalLoadHPOA3ID      = MFUnset;
static int _MDOutLocalLoadnHPOAID      = MFUnset;
static int _MDOutLocalLoadTPIAID       = MFUnset;
static int _MDOutLocalLoadHPIID        = MFUnset;

static int _MDFluxMixing_DOCID         = MFUnset;
static int _MDFluxMixing_DOC2ID        = MFUnset;
static int _MDFluxMixing_POCID         = MFUnset;
static int _MDFluxMixing_DINID         = MFUnset; 
static int _MDFluxMixing_ClID          = MFUnset; 
static int _MDFluxMixing_HPOA1ID       = MFUnset; 
static int _MDFluxMixing_HPOA2ID       = MFUnset; 
static int _MDFluxMixing_HPOA3ID       = MFUnset; 
static int _MDFluxMixing_nHPOAID       = MFUnset; 
static int _MDFluxMixing_TPIAID        = MFUnset; 
static int _MDFluxMixing_HPIID         = MFUnset; 

static int _MDFlux_DINID               = MFUnset; 
static int _MDFlux_ClID                = MFUnset; 
static int _MDFlux_DOCID               = MFUnset; 
static int _MDFlux_DOC2ID              = MFUnset; 
static int _MDFlux_HPOA1ID             = MFUnset; 
static int _MDFlux_HPOA2ID             = MFUnset; 
static int _MDFlux_HPOA3ID             = MFUnset; 
static int _MDFlux_nHPOAID             = MFUnset; 
static int _MDFlux_TPIAID              = MFUnset; 
static int _MDFlux_HPIID               = MFUnset; 

static int _MDOutConcMixing_DOCID      = MFUnset;
static int _MDOutConcMixing_DOC2ID     = MFUnset;
static int _MDOutConcMixing_POCID      = MFUnset;
static int _MDOutConcMixing_DINID      = MFUnset;
static int _MDOutConcMixing_ClID       = MFUnset;
static int _MDOutConcMixing_HPOA1ID    = MFUnset;
static int _MDOutConcMixing_HPOA2ID    = MFUnset;
static int _MDOutConcMixing_HPOA3ID    = MFUnset;
static int _MDOutConcMixing_nHPOAID    = MFUnset;
static int _MDOutConcMixing_TPIAID     = MFUnset;
static int _MDOutConcMixing_HPIID      = MFUnset;

static int _MDStoreWaterMixing_DOCID   = MFUnset;
static int _MDStoreWaterMixing_DOC2ID  = MFUnset;
static int _MDStoreWaterMixing_POCID   = MFUnset;
static int _MDStoreWaterMixing_DINID   = MFUnset;
static int _MDStoreWaterMixing_ClID    = MFUnset;
static int _MDStoreWaterMixing_HPOA1ID = MFUnset;
static int _MDStoreWaterMixing_HPOA2ID = MFUnset;
static int _MDStoreWaterMixing_HPOA3ID = MFUnset;
static int _MDStoreWaterMixing_nHPOAID = MFUnset;
static int _MDStoreWaterMixing_TPIAID  = MFUnset;
static int _MDStoreWaterMixing_HPIID   = MFUnset;

static int _MDStoreWater_DOCID         = MFUnset;
static int _MDStoreWater_DOC2ID        = MFUnset;
static int _MDStoreWater_DINID         = MFUnset;
static int _MDStoreWater_ClID          = MFUnset;
static int _MDStoreWater_HPOA1ID       = MFUnset;
static int _MDStoreWater_HPOA2ID       = MFUnset;
static int _MDStoreWater_HPOA3ID       = MFUnset;
static int _MDStoreWater_nHPOAID       = MFUnset;
static int _MDStoreWater_TPIAID        = MFUnset;
static int _MDStoreWater_HPIID         = MFUnset;

static int _MDDeltaStorageMixing_DOCID  = MFUnset;
static int _MDDeltaStorageMixing_DOC2ID = MFUnset;
static int _MDDeltaStorageMixing_POCID  = MFUnset;
static int _MDPostConcMixing_DOCID      = MFUnset;
static int _MDPostConcMixing_DOC2ID     = MFUnset;
static int _MDPostConcMixing_POCID      = MFUnset;
static int _MDOutPostConc_DOCID         = MFUnset;
static int _MDOutPostConc_DOC2ID        = MFUnset;
static int _MDOutPostConc_DINID         = MFUnset;
static int _MDOutPostConc_ClID          = MFUnset;
static int _MDOutPostConc_HPOA1ID       = MFUnset;
static int _MDOutPostConc_HPOA2ID       = MFUnset;
static int _MDOutPostConc_HPOA3ID       = MFUnset;
static int _MDOutPostConc_nHPOAID       = MFUnset;
static int _MDOutPostConc_TPIAID        = MFUnset;
static int _MDOutPostConc_HPIID         = MFUnset;
static int _MDRunoffConc_DOCID          = MFUnset;
static int _MDRunoffConc_DOC2ID         = MFUnset;
static int _MDRunoffConc_DINID          = MFUnset;
static int _MDRunoffConc_ClID           = MFUnset;
static int _MDRunoffConc_HPOA1ID        = MFUnset;
static int _MDRunoffConc_HPOA2ID        = MFUnset;
static int _MDRunoffConc_HPOA3ID        = MFUnset;
static int _MDRunoffConc_nHPOAID        = MFUnset;
static int _MDRunoffConc_TPIAID         = MFUnset;
static int _MDRunoffConc_HPIID          = MFUnset;

static void _MDDOC (int itemID) {

    float runoffVol                     = 0.0;  // m3/s
    float waterStorage                  = 0.0;  // m3/s
    float waterStorageChange            = 0.0;  // m3/s
    float localLoad_DOC                 = 0.0;  // kg/d
    float localLoad_POC                 = 0.0;  // kg/d
    float preFluxMixing_DOC             = 0.0;  // kg/d
    float preFluxMixing_POC             = 0.0;  // kg/d
    float postFluxMixing_DOC            = 0.0;  // kg/d
    float postFluxMixing_POC            = 0.0;  // kg/d
    float storeWaterMixing_DOC          = 0.0;  // kg/d
    float storeWaterMixing_POC          = 0.0;  // kg/d
    float postStoreWaterMixing_DOC      = 0.0;  // kg/d
    float postStoreWaterMixing_POC      = 0.0;  // kg/d
    float discharge                     = 0.0;  // m3/s
    float dischargePre                  = 0.0;  // m3/s
    float waterStoragePrev              = 0.0;  // m3/s
    float DOCTotalInMixing              = 0.0;  // kg/d
    float POCTotalInMixing              = 0.0;  // kg/d
    float runoffConc_DOC                = 0.0;  // mg/L
    float waterTotalVolume              = 0.0;	// m3/d
    float DIN                           = 0.0;  // kg/d
    float flowPathRemovalMixing_DOC     = 0.0;  // kg/d
    float flowPathRemovalMixing_POC     = 0.0;  // kg/d
    float postConcMixing_DOC            = 0.0;  // mg/L
    float postConcMixing_POC            = 0.0;  // mg/L
    float massBalanceMixing_DOC         = 0.0;  // kg/d
    float massBalanceMixing_POC         = 0.0;  // kg/d
    float preConcMixing_DOC             = 0.0;  // mg/L
    float preConcMixing_POC             = 0.0;  // mg/L
    float litterFall_POC                = 0.0;  // kg/d
    float wetlands                      = 0.0;  // proportion wetlands
    
    // New Variables // 
    
    float imp                           = 0.0;  // proportion impervious
    float human                         = 0.0;  // proportion human land use
    float preFluxMixing_DIN             = 0.0;  // kg/d
    float preFluxMixing_DOC2            = 0.0;  // kg/d
    float preFluxMixing_Cl              = 0.0;  // kg/d
    float preFluxMixing_HPOA1           = 0.0;  // kg/d
    float preFluxMixing_HPOA2           = 0.0;  // kg/d
    float preFluxMixing_HPOA3           = 0.0;  // kg/d
    float preFluxMixing_nHPOA           = 0.0;  // kg/d
    float preFluxMixing_TPIA            = 0.0;  // kg/d
    float preFluxMixing_HPI             = 0.0;  // kg/d
    
    float postFluxMixing_DIN            = 0.0;  // kg/d
    float postFluxMixing_DOC2           = 0.0;  // kg/d
    float postFluxMixing_Cl             = 0.0;  // kg/d
    float postFluxMixing_HPOA1          = 0.0;  // kg/d
    float postFluxMixing_HPOA2          = 0.0;  // kg/d
    float postFluxMixing_HPOA3          = 0.0;  // kg/d
    float postFluxMixing_nHPOA          = 0.0;  // kg/d
    float postFluxMixing_TPIA           = 0.0;  // kg/d
    float postFluxMixing_HPI            = 0.0;  // kg/d
    
    float storeWaterMixing_DIN          = 0.0;  // kg/d
    float storeWaterMixing_DOC2         = 0.0;  // kg/d
    float storeWaterMixing_Cl           = 0.0;  // kg/d
    float storeWaterMixing_HPOA1        = 0.0;  // kg/d
    float storeWaterMixing_HPOA2        = 0.0;  // kg/d
    float storeWaterMixing_HPOA3        = 0.0;  // kg/d
    float storeWaterMixing_nHPOA        = 0.0;  // kg/d
    float storeWaterMixing_TPIA         = 0.0;  // kg/d
    float storeWaterMixing_HPI          = 0.0;  // kg/d
   
    float postStoreWaterMixing_DIN      = 0.0;  // kg/d
    float postStoreWaterMixing_DOC2     = 0.0;  // kg/d
    float postStoreWaterMixing_Cl       = 0.0;  // kg/d
    float postStoreWaterMixing_HPOA1    = 0.0;  // kg/d
    float postStoreWaterMixing_HPOA2    = 0.0;  // kg/d
    float postStoreWaterMixing_HPOA3    = 0.0;  // kg/d
    float postStoreWaterMixing_nHPOA    = 0.0;  // kg/d
    float postStoreWaterMixing_TPIA     = 0.0;  // kg/d
    float postStoreWaterMixing_HPI      = 0.0;  // kg/d

    float DINTotalInMixing              = 0.0;  // kg/d
    float DOC2TotalInMixing             = 0.0;  // kg/d
    float ClTotalInMixing               = 0.0;  // kg/d
    float HPOA1TotalInMixing            = 0.0;  // kg/d
    float HPOA2TotalInMixing            = 0.0;  // kg/d
    float HPOA3TotalInMixing            = 0.0;  // kg/d
    float nHPOATotalInMixing            = 0.0;  // kg/d
    float TPIATotalInMixing             = 0.0;  // kg/d
    float HPITotalInMixing              = 0.0;  // kg/d

    float flowPathRemovalMixing_DIN     = 0.0;  // kg/d
    float flowPathRemovalMixing_DOC2    = 0.0;  // kg/d
    float flowPathRemovalMixing_Cl      = 0.0;  // kg/d
    float flowPathRemovalMixing_HPOA1   = 0.0;  // kg/d
    float flowPathRemovalMixing_HPOA2   = 0.0;  // kg/d
    float flowPathRemovalMixing_HPOA3   = 0.0;  // kg/d
    float flowPathRemovalMixing_nHPOA   = 0.0;  // kg/d
    float flowPathRemovalMixing_TPIA    = 0.0;  // kg/d
    float flowPathRemovalMixing_HPI     = 0.0;  // kg/d
   
    float postConcMixing_DIN            = 0.0;  // kg/d
    float postConcMixing_DOC2           = 0.0;  // kg/d
    float postConcMixing_Cl             = 0.0;  // kg/d
    float postConcMixing_HPOA1          = 0.0;  // kg/d
    float postConcMixing_HPOA2          = 0.0;  // kg/d
    float postConcMixing_HPOA3          = 0.0;  // kg/d
    float postConcMixing_nHPOA          = 0.0;  // kg/d
    float postConcMixing_TPIA           = 0.0;  // kg/d
    float postConcMixing_HPI            = 0.0;  // kg/d
   
    float massBalanceMixing_DIN         = 0.0;  // kg/d
    float massBalanceMixing_DOC2        = 0.0;  // kg/d
    float massBalanceMixing_Cl          = 0.0;  // kg/d
    float massBalanceMixing_HPOA1       = 0.0;  // kg/d
    float massBalanceMixing_HPOA2       = 0.0;  // kg/d
    float massBalanceMixing_HPOA3       = 0.0;  // kg/d
    float massBalanceMixing_nHPOA       = 0.0;  // kg/d
    float massBalanceMixing_TPIA        = 0.0;  // kg/d
    float massBalanceMixing_HPI         = 0.0;  // kg/d
    
    float preConcMixing_DIN             = 0.0;  // kg/d
    float preConcMixing_DOC2            = 0.0;  // kg/d
    float preConcMixing_Cl              = 0.0;  // kg/d
    float preConcMixing_HPOA1           = 0.0;  // kg/d
    float preConcMixing_HPOA2           = 0.0;  // kg/d
    float preConcMixing_HPOA3           = 0.0;  // kg/d
    float preConcMixing_nHPOA           = 0.0;  // kg/d
    float preConcMixing_TPIA            = 0.0;  // kg/d
    float preConcMixing_HPI             = 0.0;  // kg/d
  
    float preFlux_DIN                   = 0.0;  // kg/d
    float storeWater_DIN                = 0.0;  // kg/d
    float preFlux_DOC2                  = 0.0;  // kg/d
    float storeWater_DOC2               = 0.0;  // kg/d
    float preFlux_Cl                    = 0.0;  // kg/d
    float storeWater_Cl                 = 0.0;  // kg/d
    float preFlux_DOC                   = 0.0;  // kg/d
    float storeWater_DOC                = 0.0;  // kg/d
    float preFlux_HPOA1                 = 0.0;  // kg/d
    float storeWater_HPOA1              = 0.0;  // kg/d
    float preFlux_HPOA2                 = 0.0;  // kg/d
    float storeWater_HPOA2              = 0.0;  // kg/d
    float preFlux_HPOA3                 = 0.0;  // kg/d
    float storeWater_HPOA3              = 0.0;  // kg/d
    float preFlux_nHPOA                 = 0.0;  // kg/d
    float storeWater_nHPOA              = 0.0;  // kg/d
    float preFlux_TPIA                  = 0.0;  // kg/d
    float storeWater_TPIA               = 0.0;  // kg/d
    float preFlux_HPI                   = 0.0;  // kg/d
    float storeWater_HPI                = 0.0;  // kg/d
                           
    float runoffConc_DIN                = 0.0;  // mg/L
    float runoffConc_DOC2               = 0.0;  // mg/L
    float runoffConc_Cl                 = 0.0;  // mg/L
    float runoffConc_HPOA1              = 0.0;  // mg/L
    float runoffConc_HPOA2              = 0.0;  // mg/L
    float runoffConc_HPOA3              = 0.0;  // mg/L
    float runoffConc_nHPOA              = 0.0;  // mg/L
    float runoffConc_TPIA               = 0.0;  // mg/L
    float runoffConc_HPI                = 0.0;  // mg/L
    
    float DOCTotalIn                    = 0.0;  // kg/d
    float DOC2TotalIn                   = 0.0;  // kg/d
    float DINTotalIn                    = 0.0;  // kg/d
    float ClTotalIn                     = 0.0;  // kg/d
    float HPOA1TotalIn                  = 0.0;  // kg/d
    float HPOA2TotalIn                  = 0.0;  // kg/d
    float HPOA3TotalIn                  = 0.0;  // kg/d
    float nHPOATotalIn                  = 0.0;  // kg/d
    float TPIATotalIn                   = 0.0;  // kg/d
    float HPITotalIn                    = 0.0;  // kg/d
    
    float preConc_DIN                   = 0.0;  // mg/L
    float preConc_DOC2                  = 0.0;  // mg/L
    float preConc_DOC                   = 0.0;  // mg/L
    float preConc_Cl                    = 0.0;  // mg/L
    float preConc_HPOA1                 = 0.0;  // mg/L
    float preConc_HPOA2                 = 0.0;  // mg/L
    float preConc_HPOA3                 = 0.0;  // mg/L
    float preConc_nHPOA                 = 0.0;  // mg/L
    float preConc_TPIA                  = 0.0;  // mg/L
    float preConc_HPI                   = 0.0;  // mg/L

    float Vf                            = 35.0; // m/yr
    float HL                            = 0.0;  // m/d
    float width                         = 0.0;  // m
    float removal                       = 0.0;  // proportional removal
    
    float totalMassRemoved_DIN          = 0.0;  // kg/d
    float totalMassRemoved_DOC2         = 0.0;  // kg/d
    float totalMassRemoved_DOC          = 0.0;  // kg/d
    float totalMassRemoved_Cl           = 0.0;  // kg/d
    float totalMassRemoved_HPOA1        = 0.0;  // kg/d
    float totalMassRemoved_HPOA2        = 0.0;  // kg/d
    float totalMassRemoved_HPOA3        = 0.0;  // kg/d
    float totalMassRemoved_nHPOA        = 0.0;  // kg/d
    float totalMassRemoved_TPIA         = 0.0;  // kg/d
    float totalMassRemoved_HPI          = 0.0;  // kg/d
  
    float flowPathRemoval_DOC           = 0.0;  // kg/d
    float flowPathRemoval_DOC2          = 0.0;  // kg/d
    float flowPathRemoval_DIN           = 0.0;  // kg/d
    float flowPathRemoval_Cl            = 0.0;  // kg/d
    float flowPathRemoval_HPOA1         = 0.0;  // kg/d
    float flowPathRemoval_HPOA2         = 0.0;  // kg/d
    float flowPathRemoval_HPOA3         = 0.0;  // kg/d
    float flowPathRemoval_nHPOA         = 0.0;  // kg/d
    float flowPathRemoval_TPIA          = 0.0;  // kg/d
    float flowPathRemoval_HPI           = 0.0;  // kg/d
 
    float postConc_DIN                  = 0.0;  // mg/L
    float postConc_DOC                  = 0.0;  // mg/L
    float postConc_DOC2                 = 0.0;  // mg/L
    float postConc_Cl                   = 0.0;  // mg/L
    float postConc_HPOA1                = 0.0;  // mg/L
    float postConc_HPOA2                = 0.0;  // mg/L
    float postConc_HPOA3                = 0.0;  // mg/L
    float postConc_nHPOA                = 0.0;  // mg/L
    float postConc_TPIA                 = 0.0;  // mg/L
    float postConc_HPI                  = 0.0;  // mg/L
  
    float postFlux_DIN                  = 0.0;  // kg/d
    float postFlux_DOC                  = 0.0;  // kg/d
    float postFlux_DOC2                 = 0.0;  // kg/d
    float postFlux_Cl                   = 0.0;  // kg/d
    float postFlux_HPOA1                = 0.0;  // kg/d
    float postFlux_HPOA2                = 0.0;  // kg/d
    float postFlux_HPOA3                = 0.0;  // kg/d
    float postFlux_nHPOA                = 0.0;  // kg/d
    float postFlux_TPIA                 = 0.0;  // kg/d
    float postFlux_HPI                  = 0.0;  // kg/d
    
    float postStoreWater_DOC            = 0.0;  // kg/d
    float postStoreWater_DOC2           = 0.0;  // kg/d
    float postStoreWater_DIN            = 0.0;  // kg/d
    float postStoreWater_Cl             = 0.0;  // kg/d
    float postStoreWater_HPOA1          = 0.0;  // kg/d
    float postStoreWater_HPOA2          = 0.0;  // kg/d
    float postStoreWater_HPOA3          = 0.0;  // kg/d
    float postStoreWater_nHPOA          = 0.0;  // kg/d
    float postStoreWater_TPIA           = 0.0;  // kg/d
    float postStoreWater_HPI            = 0.0;  // kg/d
  
    float massBalance_DOC               = 0.0;  // kg/d
    float massBalance_DOC2              = 0.0;  // kg/d
    float massBalance_DIN               = 0.0;  // kg/d
    float massBalance_Cl                = 0.0;  // kg/d
    float massBalance_HPOA1             = 0.0;  // kg/d
    float massBalance_HPOA2             = 0.0;  // kg/d
    float massBalance_HPOA3             = 0.0;  // kg/d
    float massBalance_nHPOA             = 0.0;  // kg/d
    float massBalance_TPIA              = 0.0;  // kg/d
    float massBalance_HPI               = 0.0;  // kg/d
 
    float localLoad_DIN                 = 0.0;  // kg/d
    float localLoad_DOC2                = 0.0;  // kg/d
    float localLoad_Cl                  = 0.0;  // kg/d
    float localLoad_HPOA1               = 0.0;  // kg/d
    float localLoad_HPOA2               = 0.0;  // kg/d
    float localLoad_HPOA3               = 0.0;  // kg/d
    float localLoad_nHPOA               = 0.0;  // kg/d
    float localLoad_TPIA               = 0.0;  // kg/d
    float localLoad_HPI               = 0.0;  // kg/d
    
    float order                         = 0.0;  // river order
    float length                        = 0.0;  // m
    float ClTotalInMixing_Conc          = 0.0;  // kg/d
    
    float DIN_m                         = 0.0;  // slope DIN
    float DIN_b                         = 0.0;  // intercept DIN
    float Cl_m                          = 0.0;  // slope chloride
    float Cl_b                          = 0.0;  // intercept chloride
    float DOC_m                         = 0.0;  // slope DOC
    float DOC_b                         = 0.0;  // intercept DOC
    
    float x1                            = 0.0;  //
    float x2                            = 0.0;  //
    float w                             = 0.0;
    float y1                            = 0.0;
    float SEE_DOC                       = 1.0483;       // (mg/L) Standard Error of Estimate or Residual Standard Deviation
    float SEE_NO3                       = 0.4325;       // (mg/L) Standard Error of Estimate or Residual Standard Deviation
    float SEE_Cl                        = 53.6182;      // (mg/L) Standard Error of Estimate or Residual Standard Deviation
    float error_DOC                     = 0.0;          // mg/L
    float error_DOC2                    = 0.0;          // mg/L
    float error_NO3                     = 0.0;          // mg/L
    float error_Cl                      = 0.0;          // mg/L
    float error_HPOA1                   = 0.0;          // mg/L HPOA linear relationship
    float error_HPOA2                   = 0.0;          // mg/L HPOA polynomial relationship
    float error_HPOA3                   = 0.0;          // mg/L HPOA polynomial relationship
    float error_nHPOA                   = 0.0;          // mg/L nonHPOA linear relationship
    float error_TPIA                    = 0.0;          // mg/L linear
    float error_HPI                     = 0.0;          // mg/L linear
    float u1                            = 0.0;  
    float u2                            = 0.0;  
    float loadAdjust                    = 0.0;

 //             localLoad_DOC  	     = MFVarGetFloat (_MDInLocalLoad_DOCID,      itemID, 0.0); // kg/day TEM inputs
 //             localLoad_POC        = MFVarGetFloat (_MDInLitterFall_POCID,     itemID, 0.0); // kg/day TEM inputs
 //             preFluxMixing_POC    = MFVarGetFloat (_MDFluxMixing_POCID,       itemID, 0.0); // kg/day TEM inputs
 //             storeWaterMixing_POC = MFVarGetFloat (_MDStoreWaterMixing_POCID, itemID, 0.0); // kg/day TEM inputs
 
                runoffVol             = MFVarGetFloat (_MDInRunoffVolumeID,       itemID, 0.0); // m3/s
		waterStorageChange    = MFVarGetFloat (_MDInRiverStorageChgID,    itemID, 0.0); // m3/s
		waterStorage          = MFVarGetFloat (_MDInRiverStorageID,       itemID, 0.0); // m3/s 
                preFluxMixing_DOC     = MFVarGetFloat (_MDFluxMixing_DOCID,       itemID, 0.0); // kg/day 
                storeWaterMixing_DOC  = MFVarGetFloat (_MDStoreWaterMixing_DOCID, itemID, 0.0); // kg/day 
                preFluxMixing_DOC2    = MFVarGetFloat (_MDFluxMixing_DOC2ID,       itemID, 0.0); // kg/day 
                storeWaterMixing_DOC2 = MFVarGetFloat (_MDStoreWaterMixing_DOC2ID, itemID, 0.0); // kg/day 
                discharge             = MFVarGetFloat (_MDInDischargeID,          itemID, 0.0); // m3/sec, discharge leaving the grid cell, after routing!
                dischargePre	      = MFVarGetFloat (_MDInDischarge0ID,         itemID, 0.0); // m3/sec, discharge from upstream PLUS local runoff, before routing!
                DIN                   = MFVarGetFloat (_MDFlux_DIN_denitID,       itemID, 0.0); // to start DIN model 
                wetlands              = MFVarGetFloat (_MDInWetlandsID,           itemID, 0.0); // proportion wetlands
                order                 = MFVarGetFloat (_MDInRiverOrderID,         itemID, 0.0); // strahler order
                DIN_m                 = MFVarGetFloat (_MDInDINmID,               itemID, 0.0); // slope DIN
                DIN_b                 = MFVarGetFloat (_MDInDINbID,               itemID, 0.0); // intercept DIN
                DOC_m                 = MFVarGetFloat (_MDInDOCmID,               itemID, 0.0); // slope DOC
                DOC_b                 = MFVarGetFloat (_MDInDOCbID,               itemID, 0.0); // intercept DOC
                Cl_m                  = MFVarGetFloat (_MDInClmID,                itemID, 0.0); // slope chloride
                Cl_b                  = MFVarGetFloat (_MDInClbID,                itemID, 0.0); // intercept chloride
                Vf                    = MFVarGetFloat (_MDInVfID,                 itemID, 0.0); // m/yr
                        
                // New Variables //
                
                human                = MFVarGetFloat (_MDInHumanLandUseID,       itemID, 0.0);  // percent human land use
                imp                  = MFVarGetFloat (_MDInImpFractionID,        itemID, 0.0);  // proportion impervious surface
                preFluxMixing_DIN    = MFVarGetFloat (_MDFluxMixing_DINID,       itemID, 0.0);  // kg/day 
                preFluxMixing_Cl     = MFVarGetFloat (_MDFluxMixing_ClID,        itemID, 0.0);  // kg/day 
                storeWaterMixing_DIN = MFVarGetFloat (_MDStoreWaterMixing_DINID, itemID, 0.0);  // kg/day 
                storeWaterMixing_Cl  = MFVarGetFloat (_MDStoreWaterMixing_ClID,  itemID, 0.0);  // kg/day       
                preFlux_DIN          = MFVarGetFloat (_MDFlux_DINID,             itemID, 0.0);	// kg/day RJS 091108
                storeWater_DIN       = MFVarGetFloat (_MDStoreWater_DINID,       itemID, 0.0);	// kg/day RJS 091108
                preFlux_Cl           = MFVarGetFloat (_MDFlux_ClID,              itemID, 0.0);	// kg/day RJS 091108
                storeWater_Cl        = MFVarGetFloat (_MDStoreWater_ClID,        itemID, 0.0);	// kg/day RJS 091108
                preFlux_DOC          = MFVarGetFloat (_MDFlux_DOCID,             itemID, 0.0);	// kg/day RJS 091108
                storeWater_DOC       = MFVarGetFloat (_MDStoreWater_DOCID,       itemID, 0.0);	// kg/day RJS 091108
                preFlux_DOC2         = MFVarGetFloat (_MDFlux_DOC2ID,            itemID, 0.0);	// kg/day RJS 091108
                storeWater_DOC2      = MFVarGetFloat (_MDStoreWater_DOC2ID,      itemID, 0.0);	// kg/day RJS 091108
                width	             = MFVarGetFloat (_MDInRiverWidthID,    	 itemID, 0.0);	// m			// moved here 031209
                length               = MFModelGetLength(itemID) / 1000;
                loadAdjust           = MFVarGetFloat (_MDInLoadAdjustID,         itemID, 0.0);  // adjusts local load
               
                preFluxMixing_HPOA1     = MFVarGetFloat (_MDFluxMixing_HPOA1ID,        itemID, 0.0);  // kg/day 
                storeWaterMixing_HPOA1  = MFVarGetFloat (_MDStoreWaterMixing_HPOA1ID,  itemID, 0.0);  // kg/day       
  		preFlux_HPOA1           = MFVarGetFloat (_MDFlux_HPOA1ID,              itemID, 0.0);	// kg/day RJS 091108
                storeWater_HPOA1        = MFVarGetFloat (_MDStoreWater_HPOA1ID,        itemID, 0.0);	// kg/day RJS 091108
               
                preFluxMixing_HPOA2     = MFVarGetFloat (_MDFluxMixing_HPOA2ID,        itemID, 0.0);  // kg/day 
                storeWaterMixing_HPOA2  = MFVarGetFloat (_MDStoreWaterMixing_HPOA2ID,  itemID, 0.0);  // kg/day       
  		preFlux_HPOA2           = MFVarGetFloat (_MDFlux_HPOA2ID,              itemID, 0.0);	// kg/day RJS 091108
                storeWater_HPOA2        = MFVarGetFloat (_MDStoreWater_HPOA2ID,        itemID, 0.0);	// kg/day RJS 091108
                
                preFluxMixing_HPOA3     = MFVarGetFloat (_MDFluxMixing_HPOA3ID,        itemID, 0.0);  // kg/day 
                storeWaterMixing_HPOA3  = MFVarGetFloat (_MDStoreWaterMixing_HPOA3ID,  itemID, 0.0);  // kg/day       
  		preFlux_HPOA3           = MFVarGetFloat (_MDFlux_HPOA3ID,              itemID, 0.0);	// kg/day RJS 091108
                storeWater_HPOA3        = MFVarGetFloat (_MDStoreWater_HPOA3ID,        itemID, 0.0);	// kg/day RJS 091108
                              
	        preFluxMixing_nHPOA     = MFVarGetFloat (_MDFluxMixing_nHPOAID,        itemID, 0.0);  // kg/day 
                storeWaterMixing_nHPOA  = MFVarGetFloat (_MDStoreWaterMixing_nHPOAID,  itemID, 0.0);  // kg/day       
  		preFlux_nHPOA           = MFVarGetFloat (_MDFlux_nHPOAID,              itemID, 0.0);	// kg/day RJS 091108
                storeWater_nHPOA        = MFVarGetFloat (_MDStoreWater_nHPOAID,        itemID, 0.0);	// kg/day RJS 091108
               
		preFluxMixing_TPIA     = MFVarGetFloat (_MDFluxMixing_TPIAID,        itemID, 0.0);  // kg/day 
                storeWaterMixing_TPIA  = MFVarGetFloat (_MDStoreWaterMixing_TPIAID,  itemID, 0.0);  // kg/day       
  		preFlux_TPIA           = MFVarGetFloat (_MDFlux_TPIAID,              itemID, 0.0);	// kg/day RJS 091108
                storeWater_TPIA        = MFVarGetFloat (_MDStoreWater_TPIAID,        itemID, 0.0);	// kg/day RJS 091108
               
		preFluxMixing_HPI     = MFVarGetFloat (_MDFluxMixing_HPIID,        itemID, 0.0);  // kg/day 
                storeWaterMixing_HPI  = MFVarGetFloat (_MDStoreWaterMixing_HPIID,  itemID, 0.0);  // kg/day       
  		preFlux_HPI           = MFVarGetFloat (_MDFlux_HPIID,              itemID, 0.0);	// kg/day RJS 091108
                storeWater_HPI        = MFVarGetFloat (_MDStoreWater_HPIID,        itemID, 0.0);	// kg/day RJS 091108
               
		
                waterStoragePrev     = waterStorage - waterStorageChange;                                    // m3/sec     
                waterTotalVolume     = discharge * 86400;                                                    // m3/d
                HL                   = discharge > 0.0001 ? discharge / (width * length) * 86400 : 0.0;      // m/d
 //               runoffConc_DOC	     = runoffVol > 0.0 ? (localLoad_DOC * 1000000) / (runoffVol * 86400 * 1000) : 0.0;	// TEM INPUTS
 //               DOCTotalInMixing     = localLoad_DOC + preFluxMixing_DOC + storeWaterMixing_DOC;        // kg/day                TEM INPUTS
  //              if (itemID == 100) printf("HL = %f, width = %f, discharge = %f, length = %f\n", HL, width, discharge, length);
  //              runoffConc_DOC       = (5.4 + (60.0 * wetlands * 100)) / 1000000 * 12 * 1000;               
  //              DOCTotalInMixing     = (runoffConc_DOC * runoffVol * 86400 / 1000) + preFluxMixing_DOC + storeWaterMixing_DOC;        // kg/day  
  //              POCTotalInMixing     = localLoad_POC + preFluxMixing_POC + storeWaterMixing_POC;        // kg/day
 
 //               runoffConc_DOC         = pow(10, (DOC_m * (wetlands * 100)) + DOC_b) > 0.0 ? pow(10, (DOC_m * (wetlands * 100)) + DOC_b) : 0.0; // mg/L
 //               runoffConc_DIN         = pow(10, (DIN_m * (human) + DIN_b)) > 0.0 ? pow(10, (DIN_m * (human) + DIN_b)) : 0.0;                   // mg/L
 //               runoffConc_Cl          = (Cl_m * (imp * 100)) + Cl_b;                                                                           // mg/L
                
                
                /********** Read-In Error Days **********/              

                error_DOC     = MFVarGetFloat (_MDInError_DOCID,       itemID, 0.0);	// mg/L
                error_DOC2    = MFVarGetFloat (_MDInError_DOC2ID,      itemID, 0.0);	// mg/L
                error_NO3     = MFVarGetFloat (_MDInError_NO3ID,       itemID, 0.0);	// mg/L
                error_Cl      = MFVarGetFloat (_MDInError_ClID,        itemID, 0.0);	// mg/L
                error_HPOA1   = MFVarGetFloat (_MDInError_HPOA1ID,     itemID, 0.0);	// mg/L
                error_HPOA2   = MFVarGetFloat (_MDInError_HPOA2ID,     itemID, 0.0);	// mg/L
                error_HPOA3   = MFVarGetFloat (_MDInError_HPOA3ID,     itemID, 0.0);	// mg/L
                error_nHPOA   = MFVarGetFloat (_MDInError_nHPOAID,     itemID, 0.0);	// mg/L
                error_TPIA    = MFVarGetFloat (_MDInError_TPIAID,      itemID, 0.0);	// mg/L
                error_HPI     = MFVarGetFloat (_MDInError_HPIID,       itemID, 0.0);	// mg/L

                /********** Monte Carlo **********/              
         /*      
                do {
                    u1 = rand();                   
                    u1 = u1 / (RAND_MAX);
                    u2 = rand();
                    u2 = u2 / (RAND_MAX);
                    
                    x1 = 2.0 * u1 - 1.0;
                    x2 = 2.0 * u2 - 1.0;
                    w  = x1 * x1 + x2 * x2;    
                //    printf("ID = %d, %d-%d-%d, u1 = %f, u2 = %f, x1 = %f, x2 = %f, w = %f\n",itemID, MFDateGetCurrentYear(), MFDateGetCurrentMonth(), MFDateGetCurrentDay(), u1, u2, x1,x2,w);
                } while (w >= 1.0);
                
                w         = sqrt((-2.0 * log(w)) / w);
                y1        = x1 * w;
                
                error_DOC = y1 * SEE_DOC + 0.0;         // solving for value given Z (y1), standard deviation (SEE_DOC), and sample mean (0)
                error_NO3 = y1 * SEE_NO3 + 0.0;         // solving for value given Z (y1), standard deviation (SEE_NO3), and sample mean (0)
                error_Cl  = y1 * SEE_Cl + 0.0;          // solving for value given Z (y1), standard deviation (SEE_Cl), and sample mean (0)
         */       
                /*********************************/        
                
                runoffConc_DOC         = ((0.065 + 0.671 * (wetlands * 100)) + error_DOC) > 0.0 ? (0.065 + 0.671 * (wetlands * 100)) + error_DOC : 0.0;
                runoffConc_DOC2         = ((2.87 + 0.314 * (wetlands * 100) + 0.00882 * pow((wetlands * 100),2)) + error_DOC2) > 0.0 ? ((2.87 + 0.314 * (wetlands * 100) + 0.00882 * pow((wetlands * 100),2)) + error_DOC2) : 0.0;
                runoffConc_DIN         = ((-0.016 + 0.013 * human) + error_NO3) > 0.0 ? (-0.016 + 0.013 * human) + error_NO3 : 0.0;
                runoffConc_Cl          = ((3.83 + 8.9 * (imp * 100)) + error_Cl) > 0.0 ? (3.83 + 8.9 * (imp * 100)) + error_Cl : 0.0;
                runoffConc_HPOA1       = ((-2.13 + 0.546 * (wetlands * 100)) + error_HPOA1) > 0.0 ? (-2.13 + 0.546 * (wetlands * 100)) + error_HPOA1 : 0.0;
                runoffConc_HPOA2       = ((1.51 + 0.046 * (wetlands * 100) + 0.0116 * pow((wetlands * 100),2)) + error_HPOA2) > 0.0 ? (1.51 + 0.046 * (wetlands * 100) + 0.0116 * pow((wetlands * 100),2)) + error_HPOA2 : 0.0;
                runoffConc_HPOA3       = ((0.341 + 0.328 * (wetlands * 100)) + error_HPOA3) > 0.0 ? (0.341 + 0.328 * (wetlands * 100)) + error_HPOA3 : 0.0;
                runoffConc_nHPOA       = ((2.26 + 0.149 * (wetlands * 100)) + error_nHPOA) > 0.0 ? (2.26 + 0.149 * (wetlands * 100)) + error_nHPOA : 0.0;
                runoffConc_TPIA        = ((0.24 + 0.09 * (wetlands * 100)) + error_TPIA) > 0.0 ? (0.24 + 0.09 * (wetlands * 100)) + error_TPIA : 0.0;
                runoffConc_HPI         = ((0.554 + 0.0714 * (wetlands * 100)) + error_HPI) > 0.0 ? (0.554 + 0.0714 * (wetlands * 100)) + error_HPI : 0.0;

                /*********************************/        
             
             //   if (MFDateGetCurrentYear() > 2000) printf("DOC = %f, DOC2 = %f, DIN = %f, Cl = %f, HPOA1 = %f, HPOA2 = %f, nHPOA = %f, TPIA = %f, HPI = %f\n", error_DOC, error_DOC2, error_NO3, error_Cl, error_HPOA1, error_HPOA2, error_nHPOA, error_TPIA, error_HPI);
                
    //            runoffConc_DOC         = 0.065 + 0.671 * (wetlands * 100);
    //            runoffConc_DIN         = (-0.016 + 0.013 * human) > 0.0 ? -0.016 + 0.013 * human : 0.0;
    //            runoffConc_Cl          = 3.83 + 8.9 * (imp * 100);
              
                localLoad_DOC        = runoffConc_DOC * runoffVol * 86400 / 1000 * loadAdjust;                       // kg/d
                localLoad_DOC2       = runoffConc_DOC2 * runoffVol * 86400 / 1000 * loadAdjust;                      // kg/d
                localLoad_DIN        = runoffConc_DIN * runoffVol * 86400 / 1000 * loadAdjust;                       // kg/d
                localLoad_Cl         = runoffConc_Cl * runoffVol * 86400 / 1000 * loadAdjust;                        // kg/d
                localLoad_HPOA1      = runoffConc_HPOA1 * runoffVol * 86400 / 1000 * loadAdjust;                     // kg/d
                localLoad_HPOA2      = runoffConc_HPOA2 * runoffVol * 86400 / 1000 * loadAdjust;                     // kg/d
                localLoad_HPOA3      = runoffConc_HPOA3 * runoffVol * 86400 / 1000 * loadAdjust;                     // kg/d
                localLoad_nHPOA      = runoffConc_nHPOA * runoffVol * 86400 / 1000 * loadAdjust;                     // kg/d
                localLoad_TPIA       = runoffConc_TPIA * runoffVol * 86400 / 1000 * loadAdjust;                      // kg/d
                localLoad_HPI        = runoffConc_HPI * runoffVol * 86400 / 1000 * loadAdjust;                       // kg/d

                DOCTotalInMixing     = localLoad_DOC + preFluxMixing_DOC + storeWaterMixing_DOC;        // kg/d  
                DOC2TotalInMixing    = localLoad_DOC2 + preFluxMixing_DOC2 + storeWaterMixing_DOC2;        // kg/d  
                DINTotalInMixing     = localLoad_DIN + preFluxMixing_DIN + storeWaterMixing_DIN;        // kg/d  
                ClTotalInMixing      = localLoad_Cl + preFluxMixing_Cl  + storeWaterMixing_Cl;          // kg/d  
                HPOA1TotalInMixing   = localLoad_HPOA1 + preFluxMixing_HPOA1  + storeWaterMixing_HPOA1;          // kg/d  
                HPOA2TotalInMixing   = localLoad_HPOA2 + preFluxMixing_HPOA2  + storeWaterMixing_HPOA2;          // kg/d  
                HPOA3TotalInMixing   = localLoad_HPOA3 + preFluxMixing_HPOA3  + storeWaterMixing_HPOA3;          // kg/d  
                nHPOATotalInMixing   = localLoad_nHPOA + preFluxMixing_nHPOA  + storeWaterMixing_nHPOA;          // kg/d  
                TPIATotalInMixing    = localLoad_TPIA + preFluxMixing_TPIA  + storeWaterMixing_TPIA;             // kg/d  
                HPITotalInMixing     = localLoad_HPI + preFluxMixing_HPI  + storeWaterMixing_HPI;                // kg/d  
            
                ClTotalInMixing_Conc = ClTotalInMixing / waterTotalVolume * 1000;                       // mg/L
                
                DOCTotalIn           = localLoad_DOC + preFlux_DOC + storeWater_DOC;                    // kg/day  
                DOC2TotalIn          = localLoad_DOC2 + preFlux_DOC2 + storeWater_DOC2;                 // kg/day  
                DINTotalIn           = localLoad_DIN + preFlux_DIN + storeWater_DIN;                    // kg/day  
                ClTotalIn            = localLoad_Cl + preFlux_Cl  + storeWater_Cl;                      // kg/day  
                HPOA1TotalIn         = localLoad_HPOA1 + preFlux_HPOA1  + storeWater_HPOA1;             // kg/day  
                HPOA2TotalIn         = localLoad_HPOA2 + preFlux_HPOA2  + storeWater_HPOA2;             // kg/day  
                HPOA3TotalIn         = localLoad_HPOA3 + preFlux_HPOA3  + storeWater_HPOA3;             // kg/day  
                nHPOATotalIn         = localLoad_nHPOA + preFlux_nHPOA  + storeWater_nHPOA;             // kg/day  
                TPIATotalIn          = localLoad_TPIA + preFlux_TPIA  + storeWater_TPIA;                // kg/day  
                HPITotalIn           = localLoad_HPI + preFlux_HPI  + storeWater_HPI;                   // kg/day  

    //            printf("waterStoragePrev = %f, waterStorage = %f, waterStorageChange = %f\n", waterStoragePrev, waterStorage, waterStorageChange);
                
                if (discharge > 0.0000001) {
                
                  preConcMixing_DOC     = DOCTotalInMixing / waterTotalVolume * 1000;                    // mg/L      
                  preConcMixing_DOC2    = DOC2TotalInMixing / waterTotalVolume * 1000;                   // mg/L      
                  preConcMixing_DIN     = DINTotalInMixing / waterTotalVolume * 1000;                    // mg/L                        
                  preConcMixing_Cl      = ClTotalInMixing / waterTotalVolume * 1000;                     // mg/L
                  preConcMixing_HPOA1   = HPOA1TotalInMixing / waterTotalVolume * 1000;                  // mg/L
                  preConcMixing_HPOA2   = HPOA2TotalInMixing / waterTotalVolume * 1000;                  // mg/L
                  preConcMixing_HPOA3   = HPOA3TotalInMixing / waterTotalVolume * 1000;                  // mg/L
                  preConcMixing_nHPOA   = nHPOATotalInMixing / waterTotalVolume * 1000;                  // mg/L
                  preConcMixing_TPIA    = TPIATotalInMixing / waterTotalVolume * 1000;                   // mg/L
                  preConcMixing_HPI     = HPITotalInMixing / waterTotalVolume * 1000;                    // mg/L

                  preConc_DOC     = DOCTotalIn / waterTotalVolume * 1000;                                // mg/L          
                  preConc_DOC2    = DOC2TotalIn / waterTotalVolume * 1000;                               // mg/L          
                  preConc_DIN     = DINTotalIn / waterTotalVolume * 1000;                                // mg/L                                             
                  preConc_Cl      = ClTotalIn  / waterTotalVolume * 1000;                                // mg/L
                  preConc_HPOA1   = HPOA1TotalIn  / waterTotalVolume * 1000;                             // mg/L
                  preConc_HPOA2   = HPOA2TotalIn  / waterTotalVolume * 1000;                             // mg/L
                  preConc_HPOA3   = HPOA3TotalIn  / waterTotalVolume * 1000;                             // mg/L
                  preConc_nHPOA   = nHPOATotalIn  / waterTotalVolume * 1000;                             // mg/L
                  preConc_TPIA    = TPIATotalIn  / waterTotalVolume * 1000;                              // mg/L
                  preConc_HPI     = HPITotalIn  / waterTotalVolume * 1000;                               // mg/L

                        if (order > 2.0) {
                  
                                removal = 1.0 - pow(2.718281828, -1.0 * (Vf / 365) / HL);               // proportional removal

                                totalMassRemoved_DIN     = removal * DINTotalIn;                            // kg/d
                                totalMassRemoved_DOC     = removal * DOCTotalIn;                            // kg/d
                                totalMassRemoved_DOC2    = removal * DOC2TotalIn;                           // kg/d
                                totalMassRemoved_Cl      = removal * ClTotalIn;                             // kg/d
                                totalMassRemoved_HPOA1   = removal * HPOA1TotalIn;                          // kg/d
                                totalMassRemoved_HPOA2   = removal * HPOA2TotalIn;                          // kg/d
                                totalMassRemoved_HPOA3   = removal * HPOA3TotalIn;                          // kg/d
                                totalMassRemoved_nHPOA   = removal * nHPOATotalIn;                          // kg/d
                                totalMassRemoved_TPIA    = removal * TPIATotalIn;                           // kg/d
                                totalMassRemoved_HPI     = removal * HPITotalIn;                            // kg/d
                                   
 //                             preConcMixing_DOC    = DOCTotalInMixing / waterTotalVolume * 1000;            // TEM version
 //                             preConcMixing_POC    = POCTotalInMixing / waterTotalVolume * 1000;            // TEM version                                                                     
 //                             postConcMixing_DOC   = DOCTotalInMixing / (waterTotalVolume - totalEvap) * 1000;	// TEM version		    // mg/L
 //                             postConcMixing_POC   = POCTotalInMixing / (waterTotalVolume - totalEvap) * 1000;      // TEM version
          
                        }
                  
                   postConcMixing_DOC     = DOCTotalInMixing / waterTotalVolume * 1000;                   // mg/L
                   postConcMixing_DOC2    = DOC2TotalInMixing / waterTotalVolume * 1000;                  // mg/L
                   postConcMixing_DIN     = DINTotalInMixing / waterTotalVolume * 1000;                   // mg/L
                   postConcMixing_Cl      = ClTotalInMixing / waterTotalVolume * 1000;                    // mg/L
                   postConcMixing_HPOA1   = HPOA1TotalInMixing / waterTotalVolume * 1000;                 // mg/L
                   postConcMixing_HPOA2   = HPOA2TotalInMixing / waterTotalVolume * 1000;                 // mg/L
                   postConcMixing_HPOA3   = HPOA3TotalInMixing / waterTotalVolume * 1000;                 // mg/L
                   postConcMixing_nHPOA   = nHPOATotalInMixing / waterTotalVolume * 1000;                 // mg/L
                   postConcMixing_TPIA    = TPIATotalInMixing / waterTotalVolume * 1000;                  // mg/L
                   postConcMixing_HPI     = HPITotalInMixing / waterTotalVolume * 1000;                   // mg/L
                               
                }
                
                else {
                  
                  flowPathRemoval_DOC       = DOCTotalIn;                       // kg/d
                  flowPathRemoval_DOC2      = DOC2TotalIn;                      // kg/d
                  flowPathRemoval_DIN       = DINTotalIn;                       // kg/d
                  flowPathRemoval_Cl        = ClTotalIn;                        // kg/d
                  flowPathRemoval_HPOA1     = HPOA1TotalIn;                     // kg/d
                  flowPathRemoval_HPOA2     = HPOA2TotalIn;                     // kg/d
                  flowPathRemoval_HPOA3     = HPOA3TotalIn;                     // kg/d
                  flowPathRemoval_nHPOA     = nHPOATotalIn;                     // kg/d
                  flowPathRemoval_TPIA      = TPIATotalIn;                      // kg/d
                  flowPathRemoval_HPI       = HPITotalIn;                       // kg/d
                 
                  flowPathRemovalMixing_DOC   = DOCTotalInMixing;                 // kg/d
                  flowPathRemovalMixing_DOC2  = DOC2TotalInMixing;                // kg/d
                  flowPathRemovalMixing_DIN   = DINTotalInMixing;                 // kg/d
                  flowPathRemovalMixing_Cl    = ClTotalInMixing;                  // kg/d
                  flowPathRemovalMixing_HPOA1 = HPOA1TotalInMixing;               // kg/d
                  flowPathRemovalMixing_HPOA2 = HPOA2TotalInMixing;               // kg/d
                  flowPathRemovalMixing_HPOA3 = HPOA3TotalInMixing;               // kg/d
                  flowPathRemovalMixing_nHPOA = nHPOATotalInMixing;               // kg/d
                  flowPathRemovalMixing_TPIA  = TPIATotalInMixing;                // kg/d
                  flowPathRemovalMixing_HPI   = HPITotalInMixing;                 // kg/d
               
                  postConc_DOC        = 0.0;                                    // mg/L                            
                  postConc_DOC2       = 0.0;                                    // mg/L                            
                  postConc_DIN        = 0.0;                                    // mg/L
                  postConc_Cl         = 0.0;                                    // mg/L
                  postConc_HPOA1      = 0.0;                                    // mg/L
                  postConc_HPOA2      = 0.0;                                    // mg/L
                  postConc_HPOA3      = 0.0;                                    // mg/L
                  postConc_nHPOA      = 0.0;                                    // mg/L
                  postConc_TPIA       = 0.0;                                    // mg/L
                  postConc_HPI        = 0.0;                                    // mg/L
                 
                  postConcMixing_DOC        = 0.0;                              // mg/L
                  postConcMixing_DOC2       = 0.0;                              // mg/L
                  postConcMixing_DIN        = 0.0;                              // mg/L
                  postConcMixing_Cl         = 0.0;                              // mg/L
                  postConcMixing_HPOA1      = 0.0;                              // mg/L
                  postConcMixing_HPOA2      = 0.0;                              // mg/L
                  postConcMixing_HPOA3      = 0.0;                              // mg/L
                  postConcMixing_nHPOA      = 0.0;                              // mg/L
                  postConcMixing_TPIA       = 0.0;                              // mg/L
                  postConcMixing_HPI        = 0.0;                              // mg/L
                  
//                flowPathRemovalMixing_DOC = DOCTotalInMixing;         // TEM version
//                flowPathRemovalMixing_POC = POCTotalInMixing;         // TEM version
                
//                postConcMixing_DOC        = 0.0;                      // TEM version
//                postConcMixing_POC        = 0.0;                      // TEM version
                
                }
   
                postConc_DOC   = discharge > 0.0000001 ? (DOCTotalIn - totalMassRemoved_DOC - flowPathRemoval_DOC) / waterTotalVolume * 1000 : 0.0;     // mg/L
                postConc_DOC2  = discharge > 0.0000001 ? (DOC2TotalIn - totalMassRemoved_DOC2 - flowPathRemoval_DOC2) / waterTotalVolume * 1000 : 0.0;  // mg/L
                postConc_DIN   = discharge > 0.0000001 ? (DINTotalIn - totalMassRemoved_DIN - flowPathRemoval_DIN) / waterTotalVolume * 1000 : 0.0;     // mg/L
                postConc_Cl    = discharge > 0.0000001 ? (ClTotalIn - totalMassRemoved_Cl - flowPathRemoval_Cl) / waterTotalVolume * 1000 : 0.0;        // mg/L 
                postConc_HPOA1 = discharge > 0.0000001 ? (HPOA1TotalIn - totalMassRemoved_HPOA1 - flowPathRemoval_HPOA1) / waterTotalVolume * 1000 : 0.0;        // mg/L 
                postConc_HPOA2 = discharge > 0.0000001 ? (HPOA2TotalIn - totalMassRemoved_HPOA2 - flowPathRemoval_HPOA2) / waterTotalVolume * 1000 : 0.0;        // mg/L 
                postConc_HPOA3 = discharge > 0.0000001 ? (HPOA3TotalIn - totalMassRemoved_HPOA3 - flowPathRemoval_HPOA3) / waterTotalVolume * 1000 : 0.0;        // mg/L 
                postConc_nHPOA = discharge > 0.0000001 ? (nHPOATotalIn - totalMassRemoved_nHPOA - flowPathRemoval_nHPOA) / waterTotalVolume * 1000 : 0.0;        // mg/L 
                postConc_TPIA  = discharge > 0.0000001 ? (TPIATotalIn - totalMassRemoved_TPIA - flowPathRemoval_TPIA) / waterTotalVolume * 1000 : 0.0;           // mg/L 
                postConc_HPI   = discharge > 0.0000001 ? (HPITotalIn - totalMassRemoved_HPI - flowPathRemoval_HPI) / waterTotalVolume * 1000 : 0.0;              // mg/L 
                
                postFlux_DOC 	   = (discharge * MDConst_m3PerSecTOm3PerDay) * postConc_DOC / 1000;        // kg/day
                postFlux_DOC2	   = (discharge * MDConst_m3PerSecTOm3PerDay) * postConc_DOC2 / 1000;       // kg/day
                postFlux_DIN       = (discharge * MDConst_m3PerSecTOm3PerDay) * postConc_DIN / 1000;        // kg/day
                postFlux_Cl        = (discharge * MDConst_m3PerSecTOm3PerDay) * postConc_Cl / 1000;         // kg/day
                postFlux_HPOA1     = (discharge * MDConst_m3PerSecTOm3PerDay) * postConc_HPOA1 / 1000;      // kg/day
                postFlux_HPOA2     = (discharge * MDConst_m3PerSecTOm3PerDay) * postConc_HPOA2 / 1000;      // kg/day
                postFlux_HPOA3     = (discharge * MDConst_m3PerSecTOm3PerDay) * postConc_HPOA3 / 1000;      // kg/day
                postFlux_nHPOA     = (discharge * MDConst_m3PerSecTOm3PerDay) * postConc_nHPOA / 1000;      // kg/day
                postFlux_TPIA      = (discharge * MDConst_m3PerSecTOm3PerDay) * postConc_TPIA / 1000;       // kg/day
                postFlux_HPI       = (discharge * MDConst_m3PerSecTOm3PerDay) * postConc_HPI / 1000;        // kg/day

                postStoreWater_DOC   = (waterStorage * MDConst_m3PerSecTOm3PerDay) * postConc_DOC / 1000;     // kg/day
                postStoreWater_DOC2  = (waterStorage * MDConst_m3PerSecTOm3PerDay) * postConc_DOC2 / 1000;    // kg/day
                postStoreWater_DIN   = (waterStorage * MDConst_m3PerSecTOm3PerDay) * postConc_DIN / 1000;     // kg/day
                postStoreWater_Cl    = (waterStorage * MDConst_m3PerSecTOm3PerDay) * postConc_Cl / 1000;      // kg/day
                postStoreWater_HPOA1 = (waterStorage * MDConst_m3PerSecTOm3PerDay) * postConc_HPOA1 / 1000;   // kg/day
                postStoreWater_HPOA2 = (waterStorage * MDConst_m3PerSecTOm3PerDay) * postConc_HPOA2 / 1000;   // kg/day
                postStoreWater_HPOA3 = (waterStorage * MDConst_m3PerSecTOm3PerDay) * postConc_HPOA3 / 1000;   // kg/day
                postStoreWater_nHPOA = (waterStorage * MDConst_m3PerSecTOm3PerDay) * postConc_nHPOA / 1000;   // kg/day
                postStoreWater_TPIA  = (waterStorage * MDConst_m3PerSecTOm3PerDay) * postConc_TPIA / 1000;    // kg/day
                postStoreWater_HPI   = (waterStorage * MDConst_m3PerSecTOm3PerDay) * postConc_HPI / 1000;     // kg/day

                postFluxMixing_DOC 	 = (discharge * MDConst_m3PerSecTOm3PerDay) * postConcMixing_DOC / 1000;	 // kg/day
                postFluxMixing_DOC2 	 = (discharge * MDConst_m3PerSecTOm3PerDay) * postConcMixing_DOC2 / 1000;	 // kg/day
                postFluxMixing_DIN       = (discharge * MDConst_m3PerSecTOm3PerDay) * postConcMixing_DIN / 1000;         // kg/day
                postFluxMixing_Cl        = (discharge * MDConst_m3PerSecTOm3PerDay) * postConcMixing_Cl / 1000;          // kg/day
                postFluxMixing_HPOA1     = (discharge * MDConst_m3PerSecTOm3PerDay) * postConcMixing_HPOA1 / 1000;       // kg/day
                postFluxMixing_HPOA2     = (discharge * MDConst_m3PerSecTOm3PerDay) * postConcMixing_HPOA2 / 1000;       // kg/day
                postFluxMixing_HPOA3     = (discharge * MDConst_m3PerSecTOm3PerDay) * postConcMixing_HPOA3 / 1000;       // kg/day
                postFluxMixing_nHPOA     = (discharge * MDConst_m3PerSecTOm3PerDay) * postConcMixing_nHPOA / 1000;       // kg/day
                postFluxMixing_TPIA      = (discharge * MDConst_m3PerSecTOm3PerDay) * postConcMixing_TPIA / 1000;        // kg/day
                postFluxMixing_HPI       = (discharge * MDConst_m3PerSecTOm3PerDay) * postConcMixing_HPI / 1000;         // kg/day

                postStoreWaterMixing_DOC   = (waterStorage * MDConst_m3PerSecTOm3PerDay) * postConcMixing_DOC / 1000;	 // kg/day
                postStoreWaterMixing_DOC2  = (waterStorage * MDConst_m3PerSecTOm3PerDay) * postConcMixing_DOC2 / 1000;	 // kg/day
                postStoreWaterMixing_DIN   = (waterStorage * MDConst_m3PerSecTOm3PerDay) * postConcMixing_DIN / 1000;      // kg/day
                postStoreWaterMixing_Cl    = (waterStorage * MDConst_m3PerSecTOm3PerDay) * postConcMixing_Cl / 1000;       // kg/day
                postStoreWaterMixing_HPOA1 = (waterStorage * MDConst_m3PerSecTOm3PerDay) * postConcMixing_HPOA1 / 1000; // kg/day
                postStoreWaterMixing_HPOA2 = (waterStorage * MDConst_m3PerSecTOm3PerDay) * postConcMixing_HPOA2 / 1000; // kg/day
                postStoreWaterMixing_HPOA3 = (waterStorage * MDConst_m3PerSecTOm3PerDay) * postConcMixing_HPOA3 / 1000; // kg/day
                postStoreWaterMixing_nHPOA = (waterStorage * MDConst_m3PerSecTOm3PerDay) * postConcMixing_nHPOA / 1000; // kg/day
                postStoreWaterMixing_TPIA  = (waterStorage * MDConst_m3PerSecTOm3PerDay) * postConcMixing_TPIA / 1000;  // kg/day
                postStoreWaterMixing_HPI   = (waterStorage * MDConst_m3PerSecTOm3PerDay) * postConcMixing_HPI / 1000;   // kg/day

                massBalanceMixing_DOC   = (DOCTotalInMixing - (postFluxMixing_DOC + postStoreWaterMixing_DOC + flowPathRemovalMixing_DOC)) / DOCTotalInMixing;    // proportion of total kg in
                massBalanceMixing_DOC2  = (DOC2TotalInMixing - (postFluxMixing_DOC2 + postStoreWaterMixing_DOC2 + flowPathRemovalMixing_DOC2)) / DOC2TotalInMixing;    // proportion of total kg in
                massBalanceMixing_DIN   = (DINTotalInMixing - (postFluxMixing_DIN + postStoreWaterMixing_DIN + flowPathRemovalMixing_DIN)) / DINTotalInMixing;    // proportion of total kg in
                massBalanceMixing_Cl    = (ClTotalInMixing - (postFluxMixing_Cl + postStoreWaterMixing_Cl + flowPathRemovalMixing_Cl));                           // proportion of total kg in
                massBalanceMixing_HPOA1 = (HPOA1TotalInMixing - (postFluxMixing_HPOA1 + postStoreWaterMixing_HPOA1 + flowPathRemovalMixing_HPOA1));                           // proportion of total kg in
                massBalanceMixing_HPOA2 = (HPOA2TotalInMixing - (postFluxMixing_HPOA2 + postStoreWaterMixing_HPOA2 + flowPathRemovalMixing_HPOA2));                           // proportion of total kg in
                massBalanceMixing_HPOA3 = (HPOA3TotalInMixing - (postFluxMixing_HPOA3 + postStoreWaterMixing_HPOA3 + flowPathRemovalMixing_HPOA3));                           // proportion of total kg in
                massBalanceMixing_nHPOA = (nHPOATotalInMixing - (postFluxMixing_nHPOA + postStoreWaterMixing_nHPOA + flowPathRemovalMixing_nHPOA));                           // proportion of total kg in
                massBalanceMixing_TPIA  = (TPIATotalInMixing - (postFluxMixing_TPIA + postStoreWaterMixing_TPIA + flowPathRemovalMixing_TPIA));                               // proportion of total kg in
                massBalanceMixing_HPI   = (HPITotalInMixing - (postFluxMixing_HPI + postStoreWaterMixing_HPI + flowPathRemovalMixing_HPI));                                   // proportion of total kg in

                massBalance_DOC   = DOCTotalIn > 0.00001 ? (DOCTotalIn - (postFlux_DOC + postStoreWater_DOC + flowPathRemoval_DOC + totalMassRemoved_DOC)) / DOCTotalIn : 0.0;    // proportion of total kg in
                massBalance_DOC2  = DOC2TotalIn > 0.00001 ? (DOC2TotalIn - (postFlux_DOC2 + postStoreWater_DOC2 + flowPathRemoval_DOC2 + totalMassRemoved_DOC2)) / DOC2TotalIn : 0.0;    // proportion of total kg in
                massBalance_DIN   = DINTotalIn > 0.00001 ? (DINTotalIn - (postFlux_DIN + postStoreWater_DIN + flowPathRemoval_DIN + totalMassRemoved_DIN)) / DINTotalIn : 0.0;    // proportion of total kg in
                massBalance_Cl    = ClTotalIn > 0.00001 ? (ClTotalIn - (postFlux_Cl + postStoreWater_Cl + flowPathRemoval_Cl + totalMassRemoved_Cl)) : 0.0;                       // proportion of total kg in
                massBalance_HPOA1 = HPOA1TotalIn > 0.00001 ? (HPOA1TotalIn - (postFlux_HPOA1 + postStoreWater_HPOA1 + flowPathRemoval_HPOA1 + totalMassRemoved_HPOA1)) : 0.0;                       // proportion of total kg in
                massBalance_HPOA2 = HPOA2TotalIn > 0.00001 ? (HPOA2TotalIn - (postFlux_HPOA2 + postStoreWater_HPOA2 + flowPathRemoval_HPOA2 + totalMassRemoved_HPOA2)) : 0.0;                       // proportion of total kg in
                massBalance_HPOA3 = HPOA3TotalIn > 0.00001 ? (HPOA3TotalIn - (postFlux_HPOA3 + postStoreWater_HPOA3 + flowPathRemoval_HPOA3 + totalMassRemoved_HPOA3)) : 0.0;                       // proportion of total kg in
                massBalance_nHPOA = nHPOATotalIn > 0.00001 ? (nHPOATotalIn - (postFlux_nHPOA + postStoreWater_nHPOA + flowPathRemoval_nHPOA + totalMassRemoved_nHPOA)) : 0.0;                       // proportion of total kg in
                massBalance_TPIA  = TPIATotalIn > 0.00001 ? (TPIATotalIn - (postFlux_TPIA + postStoreWater_TPIA + flowPathRemoval_TPIA + totalMassRemoved_TPIA)) : 0.0;                       // proportion of total kg in
                massBalance_HPI   = HPITotalIn > 0.00001 ? (HPITotalIn - (postFlux_HPI + postStoreWater_HPI + flowPathRemoval_HPI + totalMassRemoved_HPI)) : 0.0;                             // proportion of total kg in

 //               if ((massBalance_DOC > 0.003) || (massBalance_DIN > 0.003) || (massBalance_Cl > 0.003)) {
 //                   printf("itemID = %d, %d-%d-%d, MB_DOC = %f, MB_DIN = %f, MB_Cl = %f\n",itemID, MFDateGetCurrentYear(), MFDateGetCurrentMonth(), MFDateGetCurrentDay(), massBalance_DOC, massBalance_DIN, massBalance_Cl);                  
 //                   printf("Q = %f, DINTotalIn=%f, postFlux_DIN=%f, postStoreWater_DIN=%f, flowPathRemoval_DIN=%f\n",discharge,DINTotalIn,postFlux_DIN,postStoreWater_DIN,flowPathRemoval_DIN);
  //              }              
  //              if ((massBalanceMixing_DOC > 0.00003) || (massBalanceMixing_DIN > 0.00003) || (massBalanceMixing_Cl > 0.00003)) {
  //                  printf("itemID = %d, %d-%d-%d, MIXING_DOC = %f, MIXING_DIN = %f, MIXING_Cl = %f\n",itemID, MFDateGetCurrentYear(), MFDateGetCurrentMonth(), MFDateGetCurrentDay(), massBalanceMixing_DOC, massBalanceMixing_DIN, massBalanceMixing_Cl);                  
  //                  printf("Q = %f, DOCTotalInMix = %f, postFluxMix_DOC = %f, postStoreWaterMix_DOC = %f, flowPathRemovalMix_DOC = %f\n", discharge, DOCTotalInMixing, postFluxMixing_DOC, postStoreWaterMixing_DOC, flowPathRemovalMixing_DOC);
  //              }    
                
 //               if ((itemID == 1576) || (itemID == 1568)) {
 //                printf("ID = %d, %d-%d-%d Q = %f, QPre = %f, ROVol = %f, runoffConc_Cl = %f, localLoad_Cl = %f, diff1 = %f, diff2 = %f\n",itemID, MFDateGetCurrentYear(), MFDateGetCurrentMonth(), MFDateGetCurrentDay(),discharge, dischargePre, runoffVol, runoffConc_Cl, localLoad_Cl, postConcMixing_Cl - ClTotalInMixing_Conc, postFluxMixing_Cl - ClTotalInMixing);
 //                printf("preFluxMixing_Cl = %f, ClTotalInMixing = %f, preConcMixing_Cl = %f, postFluxMixing_Cl = %f, postConcMixing_Cl = %f, MB = %f\n", preFluxMixing_Cl, ClTotalInMixing, preConcMixing_Cl, postFluxMixing_Cl, postConcMixing_Cl,massBalanceMixing_Cl);
                
 //               }
                
 //               postFluxMixing_DOC 	 = (discharge * MDConst_m3PerSecTOm3PerDay) * postConcMixing_DOC / 1000;		// kg/day
 //               postFluxMixing_POC       = (discharge * MDConst_m3PerSecTOm3PerDay) * postConcMixing_POC / 1000;         // kg/day
 //               postStoreWaterMixing_DOC = (waterStorage * MDConst_m3PerSecTOm3PerDay) * postConcMixing_DOC / 1000;	// kg/day
 //               postStoreWaterMixing_POC = (waterStorage * MDConst_m3PerSecTOm3PerDay) * postConcMixing_POC / 1000;      // kg/day

 //               massBalanceMixing_DOC = ((localLoad_DOC + preFluxMixing_DOC + storeWaterMixing_DOC) - (postFluxMixing_DOC + postStoreWaterMixing_DOC + flowPathRemovalMixing_DOC)) / (localLoad_DOC + storeWaterMixing_DOC + preFluxMixing_DOC);
 //               massBalanceMixing_POC = ((localLoad_POC + preFluxMixing_POC + storeWaterMixing_POC) - (postFluxMixing_POC + postStoreWaterMixing_POC + flowPathRemovalMixing_POC)) / (localLoad_POC + storeWaterMixing_POC + preFluxMixing_POC);
 
 //               if ((itemID == 31) || (itemID == 32)) {
 //               if (((massBalanceMixing_DOC > 0.0003) || (massBalanceMixing_POC > 0.0003)) && (localLoad_DOC + storeWaterMixing_DOC + preFluxMixing_DOC > 0.000001)) {
 //                   printf("******itemID = %d, y = %d, m = %d, d = %d, massBalance_DOC = %f, massBalance_POC = %f, preConc_DOC = %f, preConc_POC = %f\n", itemID, MFDateGetCurrentYear(), MFDateGetCurrentMonth(), MFDateGetCurrentDay(), massBalanceMixing_DOC, massBalanceMixing_POC, preConcMixing_DOC, preConcMixing_POC);
 //                   printf("Evap = %f, dischargePre = %f, discharge = %f, runoffVol = %f, waterStoragePrev = %f, waterStorage = %f, waterTotalVol = %f\n", totalEvap, dischargePre, discharge, runoffVol, waterStoragePrev, waterStorage, waterTotalVolume);
 //                   printf("totalIn_DOC = %f, totalIn_POC = %f, localLoad_DOC = %f, localLoad_POC = %f, preFlux_DOC = %f, preFlux_POC = %f, storeWater_DOC = %f, storeWater_POC = %f\n", DOCTotalInMixing, POCTotalInMixing, localLoad_DOC, localLoad_POC, preFluxMixing_DOC, preFluxMixing_POC, storeWaterMixing_DOC, storeWaterMixing_POC);
 //                   printf("postConc_DOC = %f, postConc_POC = %f, postFlux_DOC = %f, postFlux_POC = %f, postStore_DOC = %f, postStore_POC = %f, flowPathRem_DOC = %f, flowPathRem_POC = %f\n", postConcMixing_DOC, postConcMixing_POC, postFluxMixing_DOC, postFluxMixing_POC, postStoreWaterMixing_DOC, postStoreWaterMixing_POC, flowPathRemovalMixing_DOC, flowPathRemovalMixing_POC);
 //               }
  //              }
               
  //              MFVarSetFloat (_MDFluxMixing_DOCID,             itemID, postFluxMixing_DOC);
  //              MFVarSetFloat (_MDFluxMixing_POCID,             itemID, postFluxMixing_POC);
  //              MFVarSetFloat (_MDOutConcMixing_DOCID,          itemID, postConcMixing_DOC);
  //              MFVarSetFloat (_MDOutConcMixing_POCID,          itemID, postConcMixing_POC);
  //             MFVarSetFloat (_MDStoreWaterMixing_DOCID,       itemID, postStoreWaterMixing_DOC);
  //              MFVarSetFloat (_MDStoreWaterMixing_POCID,       itemID, postStoreWaterMixing_POC);
  //              MFVarSetFloat (_MDRunoffConc_DOCID,             itemID, runoffConc_DOC);
  
                MFVarSetFloat (_MDFluxMixing_DOCID,             itemID, postFluxMixing_DOC);
                MFVarSetFloat (_MDFluxMixing_DOC2ID,            itemID, postFluxMixing_DOC2);
                MFVarSetFloat (_MDFluxMixing_DINID,             itemID, postFluxMixing_DIN);
                MFVarSetFloat (_MDFluxMixing_ClID,              itemID, postFluxMixing_Cl);
                MFVarSetFloat (_MDFluxMixing_HPOA1ID,           itemID, postFluxMixing_HPOA1);
                MFVarSetFloat (_MDFluxMixing_HPOA2ID,           itemID, postFluxMixing_HPOA2);
                MFVarSetFloat (_MDFluxMixing_HPOA3ID,           itemID, postFluxMixing_HPOA3);
                MFVarSetFloat (_MDFluxMixing_nHPOAID,           itemID, postFluxMixing_nHPOA);
                MFVarSetFloat (_MDFluxMixing_TPIAID,            itemID, postFluxMixing_TPIA);
                MFVarSetFloat (_MDFluxMixing_HPIID,             itemID, postFluxMixing_HPI);
              
                MFVarSetFloat (_MDFlux_DOCID,                   itemID, postFlux_DOC);
                MFVarSetFloat (_MDFlux_DOC2ID,                  itemID, postFlux_DOC2);
                MFVarSetFloat (_MDFlux_DINID,                   itemID, postFlux_DIN);
                MFVarSetFloat (_MDFlux_ClID,                    itemID, postFlux_Cl);
                MFVarSetFloat (_MDFlux_HPOA1ID,                 itemID, postFlux_HPOA1);
                MFVarSetFloat (_MDFlux_HPOA2ID,                 itemID, postFlux_HPOA2);
                MFVarSetFloat (_MDFlux_HPOA3ID,                 itemID, postFlux_HPOA3);
                MFVarSetFloat (_MDFlux_nHPOAID,                 itemID, postFlux_nHPOA);
                MFVarSetFloat (_MDFlux_TPIAID,                  itemID, postFlux_TPIA);
                MFVarSetFloat (_MDFlux_HPIID,                   itemID, postFlux_HPI);

                MFVarSetFloat (_MDOutLocalLoadDOCID,           itemID, localLoad_DOC);
                MFVarSetFloat (_MDOutLocalLoadDOC2ID,          itemID, localLoad_DOC2);
                MFVarSetFloat (_MDOutLocalLoadDINID,           itemID, localLoad_DIN);
                MFVarSetFloat (_MDOutLocalLoadClID,            itemID, localLoad_Cl);
                MFVarSetFloat (_MDOutLocalLoadHPOA1ID,         itemID, localLoad_HPOA1);
                MFVarSetFloat (_MDOutLocalLoadHPOA2ID,         itemID, localLoad_HPOA2);
                MFVarSetFloat (_MDOutLocalLoadHPOA3ID,         itemID, localLoad_HPOA3);
                MFVarSetFloat (_MDOutLocalLoadnHPOAID,         itemID, localLoad_nHPOA);
                MFVarSetFloat (_MDOutLocalLoadTPIAID,         itemID, localLoad_TPIA);
                MFVarSetFloat (_MDOutLocalLoadHPIID,          itemID, localLoad_HPI);
              
                MFVarSetFloat (_MDOutConcMixing_DOCID,          itemID, postConcMixing_DOC);
                MFVarSetFloat (_MDOutConcMixing_DOC2ID,         itemID, postConcMixing_DOC2);
                MFVarSetFloat (_MDOutConcMixing_DINID,          itemID, postConcMixing_DIN);
                MFVarSetFloat (_MDOutConcMixing_ClID,           itemID, postConcMixing_Cl);
                MFVarSetFloat (_MDOutConcMixing_HPOA1ID,        itemID, postConcMixing_HPOA1);
                MFVarSetFloat (_MDOutConcMixing_HPOA2ID,        itemID, postConcMixing_HPOA2);
                MFVarSetFloat (_MDOutConcMixing_HPOA3ID,        itemID, postConcMixing_HPOA3);
                MFVarSetFloat (_MDOutConcMixing_nHPOAID,        itemID, postConcMixing_nHPOA);
                MFVarSetFloat (_MDOutConcMixing_TPIAID,         itemID, postConcMixing_TPIA);
                MFVarSetFloat (_MDOutConcMixing_HPIID,          itemID, postConcMixing_HPI);
              
                MFVarSetFloat (_MDOutPostConc_DOCID,            itemID, postConc_DOC);
                MFVarSetFloat (_MDOutPostConc_DOC2ID,           itemID, postConc_DOC2);
                MFVarSetFloat (_MDOutPostConc_DINID,            itemID, postConc_DIN);
                MFVarSetFloat (_MDOutPostConc_ClID,             itemID, postConc_Cl);
                MFVarSetFloat (_MDOutPostConc_HPOA1ID,          itemID, postConc_HPOA1);
                MFVarSetFloat (_MDOutPostConc_HPOA2ID,          itemID, postConc_HPOA2);
                MFVarSetFloat (_MDOutPostConc_HPOA3ID,          itemID, postConc_HPOA3);
                MFVarSetFloat (_MDOutPostConc_nHPOAID,          itemID, postConc_nHPOA);
                MFVarSetFloat (_MDOutPostConc_TPIAID,           itemID, postConc_TPIA);
                MFVarSetFloat (_MDOutPostConc_HPIID,            itemID, postConc_HPI);
               
                MFVarSetFloat (_MDStoreWaterMixing_DOCID,       itemID, postStoreWaterMixing_DOC);
                MFVarSetFloat (_MDStoreWaterMixing_DOC2ID,      itemID, postStoreWaterMixing_DOC2);
                MFVarSetFloat (_MDStoreWaterMixing_DINID,       itemID, postStoreWaterMixing_DIN);
                MFVarSetFloat (_MDStoreWaterMixing_ClID,        itemID, postStoreWaterMixing_Cl);
                MFVarSetFloat (_MDStoreWaterMixing_HPOA1ID,     itemID, postStoreWaterMixing_HPOA1);
                MFVarSetFloat (_MDStoreWaterMixing_HPOA2ID,     itemID, postStoreWaterMixing_HPOA2);
                MFVarSetFloat (_MDStoreWaterMixing_HPOA3ID,     itemID, postStoreWaterMixing_HPOA3);
                MFVarSetFloat (_MDStoreWaterMixing_nHPOAID,     itemID, postStoreWaterMixing_nHPOA);
                MFVarSetFloat (_MDStoreWaterMixing_TPIAID,      itemID, postStoreWaterMixing_TPIA);
                MFVarSetFloat (_MDStoreWaterMixing_HPIID,       itemID, postStoreWaterMixing_HPI);
              
                MFVarSetFloat (_MDStoreWater_DOCID,             itemID, postStoreWater_DOC);
                MFVarSetFloat (_MDStoreWater_DOC2ID,            itemID, postStoreWater_DOC2);
                MFVarSetFloat (_MDStoreWater_DINID,             itemID, postStoreWater_DIN);
                MFVarSetFloat (_MDStoreWater_ClID,              itemID, postStoreWater_Cl);
                MFVarSetFloat (_MDStoreWater_HPOA1ID,           itemID, postStoreWater_HPOA1);
                MFVarSetFloat (_MDStoreWater_HPOA2ID,           itemID, postStoreWater_HPOA2);
                MFVarSetFloat (_MDStoreWater_HPOA3ID,           itemID, postStoreWater_HPOA3);
                MFVarSetFloat (_MDStoreWater_nHPOAID,           itemID, postStoreWater_nHPOA);
                MFVarSetFloat (_MDStoreWater_TPIAID,            itemID, postStoreWater_TPIA);
                MFVarSetFloat (_MDStoreWater_HPIID,             itemID, postStoreWater_HPI);
               
                MFVarSetFloat (_MDRunoffConc_DOCID,              itemID, runoffConc_DOC);
                MFVarSetFloat (_MDRunoffConc_DOC2ID,             itemID, runoffConc_DOC2);
                MFVarSetFloat (_MDRunoffConc_DINID,              itemID, runoffConc_DIN);
                MFVarSetFloat (_MDRunoffConc_ClID,               itemID, runoffConc_Cl);
                MFVarSetFloat (_MDRunoffConc_HPOA1ID,            itemID, runoffConc_HPOA1);
                MFVarSetFloat (_MDRunoffConc_HPOA2ID,            itemID, runoffConc_HPOA2);
                MFVarSetFloat (_MDRunoffConc_HPOA3ID,            itemID, runoffConc_HPOA3);
                MFVarSetFloat (_MDRunoffConc_nHPOAID,            itemID, runoffConc_nHPOA);
                MFVarSetFloat (_MDRunoffConc_TPIAID,             itemID, runoffConc_TPIA);
                MFVarSetFloat (_MDRunoffConc_HPIID,              itemID, runoffConc_HPI);  
                
}

int MDDOCDef () {

	MFDefEntering ("DOC Routing");
	
        if (
	    ((_MDFlux_DIN_denitID               = MDDINDef ()) == CMfailed) ||	
            ((_MDInRiverWidthID                 = MDRiverWidthDef ())     == CMfailed) ||
//          ((_MDInLitterFall_POCID             = MDLitterFallDef ()) == CMfailed) ||
//          ((_MDInLocalLoad_DOCID              = MFVarGetID (MDVarLocalLoadDOC,                             "kg/d",    MFInput,  MFFlux,  MFBoundary))   == CMfailed) ||
            ((_MDInDischarge0ID                 = MFVarGetID (MDVarDischarge0,                               "m3/s",    MFInput,  MFState, MFBoundary))   == CMfailed) ||	
	    ((_MDInDischargeID                  = MFVarGetID (MDVarDischarge,                                "m3/s",    MFInput,  MFState, MFBoundary))   == CMfailed) ||
            ((_MDInRiverStorageID               = MFVarGetID (MDVarRiverStorage,                           "m3/day",    MFInput,  MFState, MFInitial))    == CMfailed) ||	
	    ((_MDInRiverStorageChgID            = MFVarGetID (MDVarRiverStorageChg,                        "m3/day",    MFInput,  MFState, MFBoundary))   == CMfailed) ||    
            ((_MDInRunoffVolumeID               = MFVarGetID (MDVarRunoffVolume, 		             "m3/s",    MFInput,  MFState, MFBoundary))   == CMfailed) ||
            ((_MDInTotalEvaporationID           = MFVarGetID (MDVarTotalEvaporation,                           "m3",    MFInput,  MFFlux,  MFBoundary)) == CMfailed) ||
            ((_MDInWetlandsID                   = MFVarGetID (MDVarFracWetlandArea,                             "-",    MFInput,  MFState, MFBoundary))   == CMfailed) ||           // RJS 112513
            ((_MDInHumanLandUseID               = MFVarGetID (MDVarHumanLandUse,                                "-",    MFInput,  MFState, MFBoundary))   == CMfailed) ||
            ((_MDInRiverOrderID                 = MFVarGetID (MDVarRiverOrder,                                "-",    MFInput,  MFState, MFBoundary))   == CMfailed) ||
            ((_MDRunoffConc_DOCID       	= MFVarGetID (MDVarRunoffConcDOC,               "mg/L",     MFOutput,  MFState, MFInitial))    == CMfailed) ||	
            ((_MDRunoffConc_DOC2ID       	= MFVarGetID (MDVarRunoffConcDOC2,              "mg/L",     MFOutput,  MFState, MFInitial))    == CMfailed) ||	
            ((_MDRunoffConc_DINID       	= MFVarGetID (MDVarRunoffConcDIN2,               "mg/L",     MFOutput,  MFState, MFInitial))    == CMfailed) ||	
            ((_MDRunoffConc_ClID                = MFVarGetID (MDVarRunoffConcCl,               "mg/L",     MFOutput,  MFState, MFInitial))    == CMfailed) ||	
            ((_MDRunoffConc_HPOA1ID             = MFVarGetID (MDVarRunoffConcHPOA1,            "mg/L",     MFOutput,  MFState, MFInitial))    == CMfailed) ||	
            ((_MDRunoffConc_HPOA2ID             = MFVarGetID (MDVarRunoffConcHPOA2,            "mg/L",     MFOutput,  MFState, MFInitial))    == CMfailed) ||	
            ((_MDRunoffConc_HPOA3ID             = MFVarGetID (MDVarRunoffConcHPOA3,            "mg/L",     MFOutput,  MFState, MFInitial))    == CMfailed) ||	
            ((_MDRunoffConc_nHPOAID             = MFVarGetID (MDVarRunoffConcnHPOA,            "mg/L",     MFOutput,  MFState, MFInitial))    == CMfailed) ||	
            ((_MDRunoffConc_TPIAID              = MFVarGetID (MDVarRunoffConcTPIA,             "mg/L",     MFOutput,  MFState, MFInitial))    == CMfailed) ||	
            ((_MDRunoffConc_HPIID               = MFVarGetID (MDVarRunoffConcHPI,              "mg/L",     MFOutput,  MFState, MFInitial))    == CMfailed) ||	
            ((_MDInImpFractionID                = MFVarGetID (MDVarImpFracSpatial,                              "-",    MFInput,  MFState, MFBoundary))   == CMfailed) ||
            ((_MDInDOCmID                       = MFVarGetID (MDVarDOCm,                                        "-",    MFInput,  MFState, MFBoundary))   == CMfailed) ||
            ((_MDInDOCbID                       = MFVarGetID (MDVarDOCb,                                        "-",    MFInput,  MFState, MFBoundary))   == CMfailed) ||
            ((_MDInDINmID                       = MFVarGetID (MDVarDINm,                                        "-",    MFInput,  MFState, MFBoundary))   == CMfailed) ||
            ((_MDInDINbID                       = MFVarGetID (MDVarDINb,                                        "-",    MFInput,  MFState, MFBoundary))   == CMfailed) ||
            ((_MDInClmID                        = MFVarGetID (MDVarClm,                                         "-",    MFInput,  MFState, MFBoundary))   == CMfailed) ||
            ((_MDInClbID                        = MFVarGetID (MDVarClb,                                         "-",    MFInput,  MFState, MFBoundary))   == CMfailed) ||
            ((_MDInVfID                         = MFVarGetID (MDVarVf,                                       "m/yr",    MFInput,  MFState, MFBoundary))   == CMfailed) ||
            ((_MDInError_DOCID                  = MFVarGetID (MDVarErrorDOC,                                 "mg/L",    MFInput,  MFState, MFBoundary))   == CMfailed) ||
            ((_MDInError_DOC2ID                 = MFVarGetID (MDVarErrorDOC2,                                 "mg/L",    MFInput,  MFState, MFBoundary))   == CMfailed) ||
            ((_MDInError_NO3ID                  = MFVarGetID (MDVarErrorNO3,                                 "mg/L",    MFInput,  MFState, MFBoundary))   == CMfailed) ||
            ((_MDInError_ClID                   = MFVarGetID (MDVarErrorCl,                                 "mg/L",    MFInput,  MFState, MFBoundary))   == CMfailed) ||
            ((_MDInError_HPOA1ID                = MFVarGetID (MDVarErrorHPOA1,                              "mg/L",    MFInput,  MFState, MFBoundary))   == CMfailed) ||
            ((_MDInError_HPOA2ID                = MFVarGetID (MDVarErrorHPOA2,                              "mg/L",    MFInput,  MFState, MFBoundary))   == CMfailed) ||
            ((_MDInError_HPOA3ID                = MFVarGetID (MDVarErrorHPOA3,                              "mg/L",    MFInput,  MFState, MFBoundary))   == CMfailed) ||
            ((_MDInError_nHPOAID                = MFVarGetID (MDVarErrornHPOA,                              "mg/L",    MFInput,  MFState, MFBoundary))   == CMfailed) ||
            ((_MDInError_TPIAID                 = MFVarGetID (MDVarErrorTPIA,                               "mg/L",    MFInput,  MFState, MFBoundary))   == CMfailed) ||
            ((_MDInError_HPIID                  = MFVarGetID (MDVarErrorHPI,                                "mg/L",    MFInput,  MFState, MFBoundary))   == CMfailed) ||
            ((_MDInLoadAdjustID                 = MFVarGetID (MDVarLoadAdjust,                              "-",    MFInput,  MFState, MFBoundary))   == CMfailed) ||          
            ((_MDOutLocalLoadDOCID              = MFVarGetID (MDVarLocalLoadDOC,               "kg/d",   MFOutput,  MFFlux,  MFBoundary))   == CMfailed) ||
            ((_MDOutLocalLoadDOC2ID             = MFVarGetID (MDVarLocalLoadDOC2,               "kg/d",   MFOutput,  MFFlux,  MFBoundary))   == CMfailed) ||
            ((_MDOutLocalLoadDINID              = MFVarGetID (MDVarLocalLoadDIN2,               "kg/d",   MFOutput,  MFFlux,  MFBoundary))   == CMfailed) ||
            ((_MDOutLocalLoadClID               = MFVarGetID (MDVarLocalLoadCl2,                "kg/d",   MFOutput,  MFFlux,  MFBoundary))   == CMfailed) ||
            ((_MDOutLocalLoadHPOA1ID            = MFVarGetID (MDVarLocalLoadHPOA1,              "kg/d",   MFOutput,  MFFlux,  MFBoundary))   == CMfailed) ||
            ((_MDOutLocalLoadHPOA2ID            = MFVarGetID (MDVarLocalLoadHPOA2,              "kg/d",   MFOutput,  MFFlux,  MFBoundary))   == CMfailed) ||
            ((_MDOutLocalLoadHPOA3ID            = MFVarGetID (MDVarLocalLoadHPOA3,              "kg/d",   MFOutput,  MFFlux,  MFBoundary))   == CMfailed) ||
            ((_MDOutLocalLoadnHPOAID            = MFVarGetID (MDVarLocalLoadnHPOA,              "kg/d",   MFOutput,  MFFlux,  MFBoundary))   == CMfailed) ||
            ((_MDOutLocalLoadTPIAID             = MFVarGetID (MDVarLocalLoadTPIA,               "kg/d",   MFOutput,  MFFlux,  MFBoundary))   == CMfailed) ||
            ((_MDOutLocalLoadHPIID              = MFVarGetID (MDVarLocalLoadHPI,                "kg/d",   MFOutput,  MFFlux,  MFBoundary))   == CMfailed) ||
           ((_MDOutConcMixing_DOCID         	= MFVarGetID (MDVarConcMixingDOC,    	        "mg/L",   MFOutput,  MFState, MFBoundary))   == CMfailed) ||	  
           ((_MDOutConcMixing_DOC2ID         	= MFVarGetID (MDVarConcMixingDOC2,    	        "mg/L",   MFOutput,  MFState, MFBoundary))   == CMfailed) ||	  
//          ((_MDOutConcMixing_POCID         	= MFVarGetID (MDVarConcMixingPOC,    	        "mg/L",   MFOutput,  MFState, MFBoundary))   == CMfailed) ||	               
            ((_MDOutConcMixing_DINID         	= MFVarGetID (MDVarConcMixingDIN2,    	        "mg/L",   MFOutput,  MFState, MFBoundary))   == CMfailed) ||	               
            ((_MDOutConcMixing_ClID         	= MFVarGetID (MDVarConcMixingCl,    	        "mg/L",   MFOutput,  MFState, MFBoundary))   == CMfailed) ||	               
            ((_MDOutConcMixing_HPOA1ID         	= MFVarGetID (MDVarConcMixingHPOA1,    	        "mg/L",   MFOutput,  MFState, MFBoundary))   == CMfailed) ||	               
            ((_MDOutConcMixing_HPOA2ID         	= MFVarGetID (MDVarConcMixingHPOA2,    	        "mg/L",   MFOutput,  MFState, MFBoundary))   == CMfailed) ||	               
            ((_MDOutConcMixing_HPOA3ID         	= MFVarGetID (MDVarConcMixingHPOA3,    	        "mg/L",   MFOutput,  MFState, MFBoundary))   == CMfailed) ||	               
            ((_MDOutConcMixing_nHPOAID         	= MFVarGetID (MDVarConcMixingnHPOA,    	        "mg/L",   MFOutput,  MFState, MFBoundary))   == CMfailed) ||	               
            ((_MDOutConcMixing_TPIAID         	= MFVarGetID (MDVarConcMixingTPIA,    	        "mg/L",   MFOutput,  MFState, MFBoundary))   == CMfailed) ||	               
            ((_MDOutConcMixing_HPIID         	= MFVarGetID (MDVarConcMixingHPI,    	        "mg/L",   MFOutput,  MFState, MFBoundary))   == CMfailed) ||	               
	    ((_MDStoreWaterMixing_DOCID      	= MFVarGetID (MDVarStoreWaterMixingDOC,       "kg/day",   MFOutput,  MFState, MFInitial))    == CMfailed) ||	
	    ((_MDStoreWaterMixing_DOC2ID      	= MFVarGetID (MDVarStoreWaterMixingDOC2,       "kg/day",   MFOutput,  MFState, MFInitial))    == CMfailed) ||	
//	    ((_MDStoreWaterMixing_POCID      	= MFVarGetID (MDVarStoreWaterMixingPOC,       "kg/day",   MFOutput,  MFState, MFInitial))    == CMfailed) ||	
	    ((_MDStoreWaterMixing_DINID      	= MFVarGetID (MDVarStoreWaterMixingDIN2,      "kg/day",   MFOutput,  MFState, MFInitial))    == CMfailed) ||	
	    ((_MDStoreWaterMixing_ClID      	= MFVarGetID (MDVarStoreWaterMixingCl,        "kg/day",   MFOutput,  MFState, MFInitial))    == CMfailed) ||	
	    ((_MDStoreWaterMixing_HPOA1ID      	= MFVarGetID (MDVarStoreWaterMixingHPOA1,     "kg/day",   MFOutput,  MFState, MFInitial))    == CMfailed) ||	
	    ((_MDStoreWaterMixing_HPOA2ID      	= MFVarGetID (MDVarStoreWaterMixingHPOA2,     "kg/day",   MFOutput,  MFState, MFInitial))    == CMfailed) ||	
	    ((_MDStoreWaterMixing_HPOA3ID      	= MFVarGetID (MDVarStoreWaterMixingHPOA3,     "kg/day",   MFOutput,  MFState, MFInitial))    == CMfailed) ||	
 	    ((_MDStoreWaterMixing_nHPOAID      	= MFVarGetID (MDVarStoreWaterMixingnHPOA,     "kg/day",   MFOutput,  MFState, MFInitial))    == CMfailed) ||	
  	    ((_MDStoreWaterMixing_TPIAID      	= MFVarGetID (MDVarStoreWaterMixingTPIA,      "kg/day",   MFOutput,  MFState, MFInitial))    == CMfailed) ||	
 	    ((_MDStoreWaterMixing_HPIID      	= MFVarGetID (MDVarStoreWaterMixingHPI,       "kg/day",   MFOutput,  MFState, MFInitial))    == CMfailed) ||	
          ((_MDFluxMixing_DOCID         	= MFVarGetID (MDVarFluxMixingDOC,    	      "kg/day",   MFRoute,   MFFlux,  MFBoundary))   == CMfailed) ||	 
           ((_MDFluxMixing_DOC2ID         	= MFVarGetID (MDVarFluxMixingDOC2,    	      "kg/day",   MFRoute,   MFFlux,  MFBoundary))   == CMfailed) ||	 
//          ((_MDFluxMixing_POCID         	= MFVarGetID (MDVarFluxMixingPOC,    	      "kg/day",   MFRoute,   MFFlux,  MFBoundary))   == CMfailed) ||	  
            ((_MDFluxMixing_DINID         	= MFVarGetID (MDVarFluxMixingDIN2,    	      "kg/day",   MFRoute,   MFFlux,  MFBoundary))   == CMfailed) ||	  
            ((_MDFluxMixing_ClID         	= MFVarGetID (MDVarFluxMixingCl,    	      "kg/day",   MFRoute,   MFFlux,  MFBoundary))   == CMfailed) ||	  
            ((_MDFluxMixing_HPOA1ID         	= MFVarGetID (MDVarFluxMixingHPOA1,    	      "kg/day",   MFRoute,   MFFlux,  MFBoundary))   == CMfailed) ||	  
            ((_MDFluxMixing_HPOA2ID         	= MFVarGetID (MDVarFluxMixingHPOA2,    	      "kg/day",   MFRoute,   MFFlux,  MFBoundary))   == CMfailed) ||	  
            ((_MDFluxMixing_HPOA3ID         	= MFVarGetID (MDVarFluxMixingHPOA3,    	      "kg/day",   MFRoute,   MFFlux,  MFBoundary))   == CMfailed) ||	  
            ((_MDFluxMixing_nHPOAID         	= MFVarGetID (MDVarFluxMixingnHPOA,    	      "kg/day",   MFRoute,   MFFlux,  MFBoundary))   == CMfailed) ||	  
            ((_MDFluxMixing_TPIAID         	= MFVarGetID (MDVarFluxMixingTPIA,    	      "kg/day",   MFRoute,   MFFlux,  MFBoundary))   == CMfailed) ||	  
            ((_MDFluxMixing_HPIID         	= MFVarGetID (MDVarFluxMixingHPI,    	      "kg/day",   MFRoute,   MFFlux,  MFBoundary))   == CMfailed) ||	  

            ((_MDOutPostConc_DOCID         	= MFVarGetID (MDVarPostConcDOC,    	        "mg/L",   MFOutput,  MFState, MFBoundary))   == CMfailed) ||	  
            ((_MDOutPostConc_DOC2ID         	= MFVarGetID (MDVarPostConcDOC2,    	        "mg/L",   MFOutput,  MFState, MFBoundary))   == CMfailed) ||	  
            ((_MDOutPostConc_DINID         	= MFVarGetID (MDVarPostConcDIN2,    	        "mg/L",   MFOutput,  MFState, MFBoundary))   == CMfailed) ||	               
            ((_MDOutPostConc_ClID         	= MFVarGetID (MDVarPostConcCl,    	        "mg/L",   MFOutput,  MFState, MFBoundary))   == CMfailed) ||	               
            ((_MDOutPostConc_HPOA1ID         	= MFVarGetID (MDVarPostConcHPOA1,    	        "mg/L",   MFOutput,  MFState, MFBoundary))   == CMfailed) ||	               
            ((_MDOutPostConc_HPOA2ID         	= MFVarGetID (MDVarPostConcHPOA2,    	        "mg/L",   MFOutput,  MFState, MFBoundary))   == CMfailed) ||	               
            ((_MDOutPostConc_HPOA3ID         	= MFVarGetID (MDVarPostConcHPOA3,    	        "mg/L",   MFOutput,  MFState, MFBoundary))   == CMfailed) ||	               
            ((_MDOutPostConc_nHPOAID         	= MFVarGetID (MDVarPostConcnHPOA,    	        "mg/L",   MFOutput,  MFState, MFBoundary))   == CMfailed) ||	               
            ((_MDOutPostConc_TPIAID         	= MFVarGetID (MDVarPostConcTPIA,    	        "mg/L",   MFOutput,  MFState, MFBoundary))   == CMfailed) ||	               
            ((_MDOutPostConc_HPIID         	= MFVarGetID (MDVarPostConcHPI,    	        "mg/L",   MFOutput,  MFState, MFBoundary))   == CMfailed) ||	               
	    ((_MDStoreWater_DOCID               = MFVarGetID (MDVarStoreWaterDOC,             "kg/day",   MFOutput,  MFState, MFInitial))    == CMfailed) ||	
	    ((_MDStoreWater_DOC2ID              = MFVarGetID (MDVarStoreWaterDOC2,             "kg/day",   MFOutput,  MFState, MFInitial))    == CMfailed) ||	
	    ((_MDStoreWater_DINID               = MFVarGetID (MDVarStoreWaterDIN2,            "kg/day",   MFOutput,  MFState, MFInitial))    == CMfailed) ||	
	    ((_MDStoreWater_ClID                = MFVarGetID (MDVarStoreWaterCl,              "kg/day",   MFOutput,  MFState, MFInitial))    == CMfailed) ||	
	    ((_MDStoreWater_HPOA1ID             = MFVarGetID (MDVarStoreWaterHPOA1,              "kg/day",   MFOutput,  MFState, MFInitial))    == CMfailed) ||	
	    ((_MDStoreWater_HPOA2ID             = MFVarGetID (MDVarStoreWaterHPOA2,              "kg/day",   MFOutput,  MFState, MFInitial))    == CMfailed) ||	
	    ((_MDStoreWater_HPOA3ID             = MFVarGetID (MDVarStoreWaterHPOA3,              "kg/day",   MFOutput,  MFState, MFInitial))    == CMfailed) ||	
	    ((_MDStoreWater_nHPOAID             = MFVarGetID (MDVarStoreWaternHPOA,              "kg/day",   MFOutput,  MFState, MFInitial))    == CMfailed) ||	
	    ((_MDStoreWater_TPIAID              = MFVarGetID (MDVarStoreWaterTPIA,              "kg/day",   MFOutput,  MFState, MFInitial))    == CMfailed) ||	
	    ((_MDStoreWater_HPIID               = MFVarGetID (MDVarStoreWaterHPI,              "kg/day",   MFOutput,  MFState, MFInitial))    == CMfailed) ||	
            ((_MDFlux_DOCID                     = MFVarGetID (MDVarFluxDOC,                   "kg/day",   MFRoute,   MFFlux,  MFBoundary))   == CMfailed) ||	 
            ((_MDFlux_DOC2ID                    = MFVarGetID (MDVarFluxDOC2,                   "kg/day",   MFRoute,   MFFlux,  MFBoundary))   == CMfailed) ||	 
            ((_MDFlux_DINID                     = MFVarGetID (MDVarFluxDIN2,    	      "kg/day",   MFRoute,   MFFlux,  MFBoundary))   == CMfailed) ||	  
            ((_MDFlux_ClID                      = MFVarGetID (MDVarFluxCl,                    "kg/day",   MFRoute,   MFFlux,  MFBoundary))   == CMfailed) ||	  
            ((_MDFlux_HPOA1ID                   = MFVarGetID (MDVarFluxHPOA1,                 "kg/day",   MFRoute,   MFFlux,  MFBoundary))   == CMfailed) ||	  
            ((_MDFlux_HPOA2ID                   = MFVarGetID (MDVarFluxHPOA2,                 "kg/day",   MFRoute,   MFFlux,  MFBoundary))   == CMfailed) ||	  
            ((_MDFlux_HPOA3ID                   = MFVarGetID (MDVarFluxHPOA3,                 "kg/day",   MFRoute,   MFFlux,  MFBoundary))   == CMfailed) ||	  
            ((_MDFlux_nHPOAID                   = MFVarGetID (MDVarFluxnHPOA,                 "kg/day",   MFRoute,   MFFlux,  MFBoundary))   == CMfailed) ||	  
            ((_MDFlux_TPIAID                    = MFVarGetID (MDVarFluxTPIA,                  "kg/day",   MFRoute,   MFFlux,  MFBoundary))   == CMfailed) ||	  
            ((_MDFlux_HPIID                     = MFVarGetID (MDVarFluxHPI,                   "kg/day",   MFRoute,   MFFlux,  MFBoundary))   == CMfailed) ||	  
   

	(MFModelAddFunction (_MDDOC) == CMfailed)) return (CMfailed);
        
	MFDefLeaving ("DOC Routing");
	return (_MDFluxMixing_DOCID); 
}
