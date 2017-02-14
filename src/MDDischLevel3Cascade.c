/******************************************************************************

GHAAS Water Balance/Transport Model V3.0
Global Hydrologic Archive and Analysis System
Copyright 1994-2011, UNH - CCNY/CUNY

MDDichLevel3Cascade.c

balazs.fekete@unh.edu
 * (Except balazs never did it - so Shan adapted Alex P's interpretation). 12-Sep-2014

*******************************************************************************/

#include <stdio.h>
#include <cm.h>
#include <MF.h>
#include <MD.h>
// Input
static int _MDInCascadeC0ID          = MFUnset; // SZ 091214
static int _MDInRunoffVolumeID       = MFUnset;
static int _MDInDischargeID          = MFUnset;
static int _MDInPropROStormWaterID   = MFUnset;    // RJS 100313
static int _MDInPropROSurfaceWaterID = MFUnset;    // RJS 100313
static int _MDInPropROGroundWaterID  = MFUnset;    // RJS 100313
static int _MDQFlagID                = MFUnset;    // RJS 012214
static int _MDInPropROFlagID         = MFUnset;    // RJS 012214
static int _MDPropRSFlagID         = MFUnset;    // RJS 012214

static int _MDQStormWaterID    = MFUnset;    // RJS 100313
static int _MDQSurfaceWaterID  = MFUnset;    // RJS 100313
static int _MDQGroundWaterID   = MFUnset;    // RJS 100313
static int _MDPropRSStormWaterID   = MFUnset;    // RJS 100313
static int _MDPropRSSurfaceWaterID = MFUnset;    // RJS 100313
static int _MDPropRSGroundWaterID  = MFUnset;    // RJS 100313
static int _MDRunoffFlagID	   = MFUnset;	 // RJS 011714
// Output
static int _MDOutDischAux0ID    = MFUnset;
static int _MDOutDischAux1ID    = MFUnset;
static int _MDOutDischLevel3ID  = MFUnset;
static int _MDOutRiverStorChgID = MFUnset;
static int _MDOutRiverStorageID = MFUnset;
static int _MDOutQPreID         = MFUnset;
static int _MDOutQCurID         = MFUnset;
static int _MDOutQOutID         = MFUnset;
static int _MDOutDischRJSID     = MFUnset;

static int _MDOutQinStormWaterID    = MFUnset;    // RJS 100313     This is simply to keep documentation of the prop coming from upstream
static int _MDOutQinSurfaceWaterID  = MFUnset;    // RJS 100313     This is simply to keep documentation of the prop coming from upstream
static int _MDOutQinGroundWaterID   = MFUnset;    // RJS 100313     This is simply to keep documentation of the prop coming from upstream
static int _MDOutPropRSinStormWaterID   = MFUnset;    // RJS 100313     This is simply to keep documentation of the prop remaining from prev storage
static int _MDOutPropRSinSurfaceWaterID = MFUnset;    // RJS 100313     This is simply to keep documentation of the prop remaining from prev storage
static int _MDOutPropRSinGroundWaterID  = MFUnset;    // RJS 100313     This is simply to keep documentation of the prop remaining from prev storage
static int _MDOutQinFlagID              = MFUnset;    // RJS 012214     This is simply to keep documentation of the prop coming from upstream
static int _MDOutPropRSinFlagID         = MFUnset;    // RJS 012214     This is simply to keep documentation of the prop remaining from prev storage


static void _MDDischLevel3Cascade(int itemID) {
    // Input
	float C0;              // Linear Routing Storage parameter (dimensionless)
	float runoff;          // Runoff [mm/dt]
// Output
	float inDischCurrent;  // Upstream discharge at the current time step [m3/s]
	float outDisch;        // Downstream discharge [m3/s]
// Local
	//float inDischPrevious; // Upstream discharge at the previous time step [m3/s] // Not needed for cascade routing
	float storChg;         // River Storage Change [m3]
	float storage;         // River Storage [m3]
        float storagePre;      // River Storage [m3] at beginning of timestep]
        float outDisch_RJS    = 0.0;       // RJS 100213 simply to keep track of initial outDisch value in timestep
        float propStW_RO      = 0.0;       // RJS 100313 proportion of runoff that is Stormwater
        float propSuW_RO      = 0.0;       // RJS 100313 proportion of runoff that is surface water
        float propGrW_RO      = 0.0;       // RJS 100313 proportion of runoff that is ground water
        
        float StW_Qin       = 0.0;       // RJS 100313 volume of incoming discharge that is Stormwater
        float SuW_Qin       = 0.0;       // RJS 100313 volume of incoming discharge that is surface water
        float GrW_Qin       = 0.0;       // RJS 100313 volume of incoming discharge that is ground water
        
        float StW_Qout       = 0.0;       // RJS 100313 volume of outgoing discharge that is Stormwater
        float SuW_Qout       = 0.0;       // RJS 100313 volume of outgoing discharge that is surface water
        float GrW_Qout       = 0.0;       // RJS 100313 volume of outgoing discharge that is ground water
        
        float propStW_out       = 0.0;       // RJS 100313 outgoing proportion that is Stormwater
        float propSuW_out       = 0.0;       // RJS 100313 outgoing proportion that is surface water
        float propGrW_out       = 0.0;       // RJS 100313 outgoing proportion that is ground water
        
        float propStW_RSin       = 0.0;       // RJS 100313 proportion of 'incoming' river storage that is Stormwater
        float propSuW_RSin       = 0.0;       // RJS 100313 proportion of 'incoming' river storage that is surface water
        float propGrW_RSin       = 0.0;       // RJS 100313 proportion of 'incoming' river storage that is ground water
        
        float StW_RO            = 0.0;
        float SuW_RO            = 0.0;
        float GrW_RO            = 0.0;
  
        float StW_RSin          = 0.0;
        float SuW_RSin          = 0.0;
        float GrW_RSin          = 0.0;    
        
        float balance1           = 0.0;
        float balance2           = 0.0;

	//Note: Flagged (source mapping) discharge is not implemented in Cascade (linear reservoir) routing

	C0 = MFVarGetFloat (_MDInCascadeC0ID,   itemID, 1.0);

	runoff          = MFVarGetFloat (_MDInRunoffVolumeID,     itemID, 0.0);
 	//inDischPrevious = MFVarGetFloat (_MDOutDischAux0ID,       itemID, 0.0);
	outDisch        = MFVarGetFloat (_MDOutDischAux1ID,       itemID, 0.0);
	inDischCurrent  = MFVarGetFloat (_MDInDischargeID,        itemID, 0.0) + 0.5 * runoff;
	storage         = MFVarGetFloat (_MDOutRiverStorageID,    itemID, 0.0);
        
        propStW_RO      = MFVarGetFloat (_MDInPropROStormWaterID,   itemID, 0.0);
        propSuW_RO      = MFVarGetFloat (_MDInPropROSurfaceWaterID, itemID, 0.0);
        propGrW_RO      = MFVarGetFloat (_MDInPropROGroundWaterID,  itemID, 0.0);     
        
        StW_Qin     = MFVarGetFloat (_MDQStormWaterID,      itemID, 0.0);
        SuW_Qin     = MFVarGetFloat (_MDQSurfaceWaterID,    itemID, 0.0);
        GrW_Qin     = MFVarGetFloat (_MDQGroundWaterID,     itemID, 0.0);
        
        propStW_RSin    = MFVarGetFloat (_MDPropRSStormWaterID,     itemID, 0.0);
        propSuW_RSin    = MFVarGetFloat (_MDPropRSSurfaceWaterID,   itemID, 0.0);
        propGrW_RSin    = MFVarGetFloat (_MDPropRSGroundWaterID,    itemID, 0.0);

        StW_RO = 0.5*runoff * propStW_RO;                           // m3/s 
        SuW_RO = 0.5*runoff * propSuW_RO;                           // m3/s
        GrW_RO = 0.5*runoff * propGrW_RO;                           // m3/s
            
        StW_RSin = storage * propStW_RSin;                      // m3/s 
        SuW_RSin = storage * propSuW_RSin;                      // m3/s
        GrW_RSin = storage * propGrW_RSin;                      // m3/s
        
	balance1 = (StW_RO + SuW_RO + GrW_RO + StW_Qin + SuW_Qin + GrW_Qin + StW_RSin + SuW_RSin + GrW_RSin) - (inDischCurrent + storage );
        
        propStW_out = (inDischCurrent + storage) > 0.0 ? (StW_RO + StW_Qin + StW_RSin) / (inDischCurrent + storage ) : 0.33333;
        propSuW_out = (inDischCurrent + storage) > 0.0 ? (SuW_RO + SuW_Qin + SuW_RSin) / (inDischCurrent + storage ) : 0.33333;
        propGrW_out = (inDischCurrent + storage) > 0.0 ? (GrW_RO + GrW_Qin + GrW_RSin) / (inDischCurrent + storage ) : 0.33333;
        
       	
        outDisch_RJS = outDisch;   // RJS 100213
        
        // Routing Calc
	outDisch = C0 * (inDischCurrent+storage) + 0.5 * runoff; 
        outDisch = outDisch < 0.0 ? inDischCurrent + 0.5 * runoff: outDisch;                  // RJS 100213 needed for storm water routing (impervious + lakes)
	storagePre = storage; 
        storage  = (inDischCurrent+storage)*(1-C0);
        if (storage < 0.0){
            storChg = -storagePre;
            storage = 0.0; }
        else { 
            storChg = storage - storagePre;
        }
        
       StW_Qout = propStW_out * (outDisch - 0.5*runoff) + StW_RO; // Add The other half of runoff disaggregated by flowpath
        SuW_Qout = propSuW_out * (outDisch - 0.5*runoff) + SuW_RO;
        GrW_Qout = propGrW_out * (outDisch - 0.5*runoff) + GrW_RO;

        balance2 = (StW_Qout + SuW_Qout + GrW_Qout) - outDisch;
     
        
        if (((inDischCurrent > 0.00001)&&(fabs(balance1/inDischCurrent) > 0.00001)) || ((outDisch >0.00001)&&(fabs(balance2/outDisch) > 0.00001))) {
            printf("\n\n itemID: %i [%i-%i-%i] %f:  %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f \n",itemID,MFDateGetCurrentYear(),MFDateGetCurrentMonth(),MFDateGetCurrentDay(),balance1, inDischCurrent, runoff, storagePre, StW_RO, SuW_RO, GrW_RO, StW_Qin, SuW_Qin, GrW_Qin, StW_RSin, SuW_RSin, GrW_RSin);
            printf("%f, %f, %f, %f, %f,  \n", balance2, StW_Qout, SuW_Qout, GrW_Qout, outDisch );
        }


        
//	if (itemID == 0 || itemID == 9245 || itemID == 11090) printf("**Discharge Cascade** itemID=%d, day = %d, outDisch = %f, inDischCurrent = %f, C0 = %f, storage = %f, storagePre = %f, storChg = %f\n", itemID, MFDateGetCurrentDay(), outDisch, inDischCurrent, C0,storage,storagePre,storChg);
//	if (itemID == 33 || itemID == 32) printf("C0 = %f, C1 = %f, C2 = %f, storage = %f, storChg = %f\n", C0, C1, C2, storage, storChg);

	MFVarSetFloat (_MDOutDischAux0ID,            itemID, inDischCurrent);
//	MFVarSetFloat (_MDOutDischAux0ID,            itemID, outDisch);		// 030113 RJS
	MFVarSetFloat (_MDOutDischAux1ID,            itemID, outDisch);
	MFVarSetFloat (_MDOutDischLevel3ID,          itemID, outDisch);
	MFVarSetFloat (_MDOutRiverStorChgID,         itemID, storChg);
	MFVarSetFloat (_MDOutRiverStorageID,         itemID, storage);
        MFVarSetFloat (_MDOutQOutID,                 itemID, outDisch);                  // RJS 100213 documentation for water balance
        MFVarSetFloat (_MDOutQCurID,                 itemID, inDischCurrent);            // RJS 100213 documentation for water balance
        MFVarSetFloat (_MDOutDischRJSID,             itemID, outDisch_RJS);              // RJS 100213 documentation for water balance
        MFVarSetFloat (_MDQStormWaterID,             itemID, StW_Qout);           // RJS 100313 Volume routed downstream
        MFVarSetFloat (_MDQSurfaceWaterID,           itemID, SuW_Qout);           // RJS 100313 Volume routed downstream
        MFVarSetFloat (_MDQGroundWaterID,            itemID, GrW_Qout);           // RJS 100313 Volume routed downstream
        MFVarSetFloat (_MDPropRSStormWaterID,        itemID, propStW_out);           // RJS 100313 Prop remaining in storage
        MFVarSetFloat (_MDPropRSSurfaceWaterID,      itemID, propSuW_out);           // RJS 100313 Prop remaining in storage
        MFVarSetFloat (_MDPropRSGroundWaterID,       itemID, propGrW_out);           // RJS 100313 Prop remaining in storage
        MFVarSetFloat (_MDOutQinStormWaterID,        itemID, StW_Qin);      // RJS 100313 Documentation of incoming prop from upstream
        MFVarSetFloat (_MDOutQinSurfaceWaterID,      itemID, SuW_Qin);      // RJS 100313 Documentation of incoming prop from upstream
        MFVarSetFloat (_MDOutQinGroundWaterID,       itemID, GrW_Qin);      // RJS 100313 Documentation of incoming prop from upstream
        MFVarSetFloat (_MDOutPropRSinStormWaterID,   itemID, propStW_out);     // RJS 100313 Documentation of initial prop for storage
        MFVarSetFloat (_MDOutPropRSinSurfaceWaterID, itemID, propSuW_out);     // RJS 100313 Documentation of initial prop for storage
        MFVarSetFloat (_MDOutPropRSinGroundWaterID,  itemID, propGrW_out);     // RJS 100313 Documentation of initial prop for storage    
}

int MDDischLevel3CascadeDef () {
    
    if (_MDOutDischLevel3ID != MFUnset) return (_MDOutDischLevel3ID);

	MFDefEntering ("Discharge Cascading Reservoirs");

	if (((_MDInRunoffVolumeID          = MDRunoffVolumeDef ()) == CMfailed) ||
	    ((_MDInCascadeC0ID             = MDDischLevel3CascadeCoeffDef ()) == CMfailed) ||
	    ((_MDInDischargeID             = MFVarGetID (MDVarDischarge,        "m3/s",   MFInput,  MFState, MFBoundary)) == CMfailed) ||
	    ((_MDInPropROStormWaterID      = MFVarGetID (MDVarPropROStormWater,    "-",   MFInput, MFState, MFBoundary)) == CMfailed) ||       // RJS 100313
            ((_MDInPropROSurfaceWaterID    = MFVarGetID (MDVarPropROSurfaceWater,  "-",   MFInput, MFState, MFBoundary)) == CMfailed) ||       // RJS 100313 
            ((_MDInPropROGroundWaterID     = MFVarGetID (MDVarPropROGroundWater,   "-",   MFInput, MFState, MFBoundary)) == CMfailed) ||       // RJS 100313                            
            ((_MDInPropROFlagID            = MFVarGetID (MDVarPropROFlag,          "-",   MFInput, MFState, MFBoundary)) == CMfailed) ||
            ((_MDQFlagID                   = MFVarGetID (MDVarQFlag,            "m3/s",   MFRoute, MFState, MFBoundary)) == CMfailed) ||
            ((_MDPropRSFlagID              = MFVarGetID (MDVarPropRSFlag,       "-",      MFOutput, MFState, MFInitial)) == CMfailed) || 
            ((_MDOutQinFlagID              = MFVarGetID (MDVarPropQinFlag,      "-",      MFOutput, MFState, MFBoundary))  == CMfailed) ||
            ((_MDOutPropRSinFlagID         = MFVarGetID (MDVarPropRSinFlag,     "-",      MFOutput, MFState, MFBoundary))  == CMfailed) ||
            ((_MDOutDischAux0ID            = MFVarGetID (MDVarDischarge0,       "m3/s",   MFOutput, MFState, MFInitial))  == CMfailed) ||
	    ((_MDOutDischAux1ID            = MFVarGetID (MDVarDischarge1,       "m3/s",   MFOutput, MFState, MFInitial))  == CMfailed) ||
	    ((_MDOutDischLevel3ID          = MFVarGetID ("__DischLevel3",       "m3/s",   MFOutput, MFState, MFBoundary)) == CMfailed) ||
	    ((_MDOutRiverStorChgID         = MFVarGetID (MDVarRiverStorageChg,  "m3",     MFOutput, MFFlux,  MFBoundary)) == CMfailed) ||
	    ((_MDOutRiverStorageID         = MFVarGetID (MDVarRiverStorage,     "m3",     MFOutput, MFState, MFInitial))  == CMfailed) ||
            ((_MDOutQPreID                 = MFVarGetID (MDVarQPre,             "m3/s",   MFOutput, MFState, MFInitial))  == CMfailed) ||
            ((_MDOutQOutID                 = MFVarGetID (MDVarQOut,             "m3/s",   MFOutput, MFState, MFInitial))  == CMfailed) ||
            ((_MDOutQCurID                 = MFVarGetID (MDVarQCur,             "m3/s",   MFOutput, MFState, MFInitial))  == CMfailed) ||           
            ((_MDQStormWaterID             = MFVarGetID (MDVarQStormWater,    "m3/s",    MFRoute, MFState, MFBoundary)) == CMfailed) ||       // RJS 100313
            ((_MDQSurfaceWaterID           = MFVarGetID (MDVarQSurfaceWater,  "m3/s",    MFRoute, MFState, MFBoundary)) == CMfailed) ||       // RJS 100313 
            ((_MDQGroundWaterID            = MFVarGetID (MDVarQGroundWater,   "m3/s",    MFRoute, MFState, MFBoundary)) == CMfailed) ||       // RJS 100313  
	    ((_MDPropRSStormWaterID        = MFVarGetID (MDVarPropRSStormWater,   "-",   MFOutput, MFState, MFInitial)) == CMfailed) ||       // RJS 100313
            ((_MDPropRSSurfaceWaterID      = MFVarGetID (MDVarPropRSSurfaceWater, "-",   MFOutput, MFState, MFInitial)) == CMfailed) ||       // RJS 100313 
            ((_MDPropRSGroundWaterID       = MFVarGetID (MDVarPropRSGroundWater,  "-",   MFOutput, MFState, MFInitial)) == CMfailed) ||       // RJS 100313  	    
            ((_MDOutQinStormWaterID    = MFVarGetID (MDVarPropQinStormWater,   "-",   MFOutput, MFState, MFBoundary)) == CMfailed) ||       // RJS 100313
            ((_MDOutQinSurfaceWaterID  = MFVarGetID (MDVarPropQinSurfaceWater,   "-", MFOutput, MFState, MFBoundary)) == CMfailed) ||       // RJS 100313
            ((_MDOutQinGroundWaterID   = MFVarGetID (MDVarPropQinGroundWater,   "-",  MFOutput, MFState, MFBoundary)) == CMfailed) ||       // RJS 100313
            ((_MDOutPropRSinStormWaterID   = MFVarGetID (MDVarPropRSinStormWater,   "-",  MFOutput, MFState, MFBoundary)) == CMfailed) ||       // RJS 100313
            ((_MDOutPropRSinSurfaceWaterID = MFVarGetID (MDVarPropRSinSurfaceWater,  "-", MFOutput, MFState, MFBoundary)) == CMfailed) ||       // RJS 100313
            ((_MDOutPropRSinGroundWaterID  = MFVarGetID (MDVarPropRSinGroundWater,   "-", MFOutput, MFState, MFBoundary)) == CMfailed) ||       // RJS 100313     
            ((_MDOutDischRJSID          = MFVarGetID (MDVarDischRJS,         "m3/s",   MFOutput, MFState, MFInitial))  == CMfailed) ||
            (MFModelAddFunction(_MDDischLevel3Cascade) == CMfailed)) return (CMfailed);
	MFDefLeaving ("Discharge Cascading Reservoirs");
	return (_MDOutDischLevel3ID);
	
}
