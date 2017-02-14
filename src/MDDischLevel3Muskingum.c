/******************************************************************************

GHAAS Water Balance/Transport Model V3.0
Global Hydrologic Archive and Analysis System
Copyright 1994-2011, UNH - CCNY/CUNY

MDDichLevel3Muskingum.c

balazs.fekete@unh.edu

*******************************************************************************/

#include <stdio.h>
#include <math.h>
#include <cm.h>
#include <MF.h>
#include <MD.h>

// Input
static int _MDInMuskingumC0ID        = MFUnset;
static int _MDInMuskingumC1ID        = MFUnset;
static int _MDInMuskingumC2ID        = MFUnset;
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


static void _MDDischLevel3Muskingum (int itemID) {
// Input
	float C0;              // Muskingum C0 coefficient (current inflow)
	float C1;              // Muskingum C1 coefficient (previous inflow)
	float C2;              // MUskingum C2 coefficient (previous outflow) 
	float runoff;          // Runoff [mm/dt]
// Output
	float inDischCurrent;  // Upstream discharge at the current time step [m3/s]
	float outDisch;        // Downstream discharge [m3/s]
// Local
	float inDischPrevious; // Upstream discharge at the previous time step [m3/s]
	float storChg;         // River Storage Change [m3]
	float storage;         // River Storage [m3]
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

        float propFlag_RO        = 0.0;       // RJS 012214
        float propFlag_RSin      = 0.0;       // RJS 012214
        float propFlag_out       = 0.0;       // RJS 012214
        
        float Flag_Qin           = 0.0;       // RJS 012214
        float Flag_Qout          = 0.0;       // RJS 012214
        float Flag_RO            = 0.0;       // RJS 012214
        float Flag_RSin          = 0.0;       // RJS 012214
              
	C0 = MFVarGetFloat (_MDInMuskingumC0ID,   itemID, 1.0);
	C1 = MFVarGetFloat (_MDInMuskingumC1ID,   itemID, 0.0);
	C2 = MFVarGetFloat (_MDInMuskingumC2ID,   itemID, 0.0);

	runoff          = MFVarGetFloat (_MDInRunoffVolumeID,     itemID, 0.0);
 	inDischPrevious = MFVarGetFloat (_MDOutDischAux0ID,       itemID, 0.0);
	outDisch        = MFVarGetFloat (_MDOutDischAux1ID,       itemID, 0.0);
	inDischCurrent  = MFVarGetFloat (_MDInDischargeID,        itemID, 0.0) + runoff;
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

        Flag_Qin	= MFVarGetFloat (_MDQFlagID,            itemID, 0.0);	// RJS 011714
	propFlag_RO	= MFVarGetFloat (_MDInPropROFlagID,     itemID, 0.0);	// RJS 011714
	propFlag_RSin	= MFVarGetFloat (_MDPropRSFlagID,       itemID, 0.0);   // RJS 011714

	Flag_RO   = runoff * propFlag_RO; 			// m3/s
	Flag_RSin = storage * propFlag_RSin;			// m3/s
	
	propFlag_out = (inDischCurrent + storage) > 0.0 ? (Flag_RO + Flag_Qin + Flag_RSin) / (inDischCurrent + storage) : 0.5;
                
        StW_RO = runoff * propStW_RO;                           // m3/s 
        SuW_RO = runoff * propSuW_RO;                           // m3/s
        GrW_RO = runoff * propGrW_RO;                           // m3/s
            
        StW_RSin = storage * propStW_RSin;                      // m3/s 
        SuW_RSin = storage * propSuW_RSin;                      // m3/s
        GrW_RSin = storage * propGrW_RSin;                      // m3/s
        
        propStW_out = (inDischCurrent + storage) > 0.0 ? (StW_RO + StW_Qin + StW_RSin) / (inDischCurrent + storage) : 0.33333;
        propSuW_out = (inDischCurrent + storage) > 0.0 ? (SuW_RO + SuW_Qin + SuW_RSin) / (inDischCurrent + storage) : 0.33333;
        propGrW_out = (inDischCurrent + storage) > 0.0 ? (GrW_RO + GrW_Qin + GrW_RSin) / (inDischCurrent + storage) : 0.33333;
        
        balance1 = (StW_RO + SuW_RO + GrW_RO + StW_Qin + SuW_Qin + GrW_Qin + StW_RSin + SuW_RSin + GrW_RSin) - (inDischCurrent + storage);
        
        outDisch_RJS = outDisch;   // RJS 100213
        
// original
	outDisch = C0 * inDischCurrent + C1 * inDischPrevious + C2 * outDisch;
        outDisch = outDisch < 0.0 ? inDischCurrent : outDisch;                  // RJS 100213 needed for storm water routing (impervious + lakes)
	storChg  = inDischCurrent - outDisch;
	storage  = storage + storChg > 0.0 ? storage + storChg : 0.0;
	
        StW_Qout = propStW_out * outDisch;
        SuW_Qout = propSuW_out * outDisch;
        GrW_Qout = propGrW_out * outDisch;

	Flag_Qout = propFlag_out * outDisch;
        
        balance2 = (StW_Qout + SuW_Qout + GrW_Qout) - outDisch;
     
        
//if (MFDateGetCurrentYear() >= 2000) {        
//        if ((itemID == 453) || (itemID == 233)) {
//            printf("%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f,",inDischCurrent, runoff, storage, StW_RO, SuW_RO, GrW_RO, StW_Qin, SuW_Qin, GrW_Qin, StW_RSin, SuW_RSin, GrW_RSin);
//            printf("%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, ", StW_Qin, SuW_Qin, GrW_Qin, propStW_RSin, propSuW_RSin, propGrW_RSin, propStW_out, propSuW_out, propGrW_out, balance1, StW_Qout, SuW_Qout, GrW_Qout, outDisch, balance2, propFlag_RO, Flag_Qin, propFlag_RSin, Flag_Qout, propFlag_out);
//        }
//}

       
//	if (itemID == 33 || itemID == 32) printf("**Discharge Musk** itemID=%d, day = %d, outDisch = %f, inDischCurrent = %f, inDischPrevious = %f\n", itemID, MFDateGetCurrentDay(), outDisch, inDischCurrent, inDischPrevious);
//	if (itemID == 33 || itemID == 32) printf("C0 = %f, C1 = %f, C2 = %f, storage = %f, storChg = %f\n", C0, C1, C2, storage, storChg);

	MFVarSetFloat (_MDOutDischAux0ID,            itemID, inDischCurrent);
//	MFVarSetFloat (_MDOutDischAux0ID,            itemID, outDisch);		// 030113 RJS
	MFVarSetFloat (_MDOutDischAux1ID,            itemID, outDisch);
	MFVarSetFloat (_MDOutDischLevel3ID,          itemID, outDisch);
	MFVarSetFloat (_MDOutRiverStorChgID,         itemID, storChg);
	MFVarSetFloat (_MDOutRiverStorageID,         itemID, storage);
        MFVarSetFloat (_MDQFlagID,                   itemID, Flag_Qout);                 // RJS 012214   Volume of Flag routed downstream
        MFVarSetFloat (_MDPropRSFlagID,              itemID, propFlag_out);              // RJS 012214   Proportion of Flag remaining in storage
        MFVarSetFloat (_MDOutQinFlagID,              itemID, Flag_Qin);                  // RJS 012214   Documentation of incoming prop from upstream
        MFVarSetFloat (_MDOutPropRSinFlagID,         itemID, propFlag_RSin);             // RJS 012214   Documentation of intial prop for storage
        MFVarSetFloat (_MDOutQPreID,                 itemID, inDischPrevious);           // RJS 100213 documentation for water balance
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

int MDDischLevel3MuskingumDef () {

	if (_MDOutDischLevel3ID != MFUnset) return (_MDOutDischLevel3ID);

	MFDefEntering ("Discharge Muskingum");

	if (((_MDInRunoffVolumeID          = MDRunoffVolumeDef ()) == CMfailed) ||
	    ((_MDInMuskingumC0ID           = MDDischLevel3MuskingumCoeffDef ()) == CMfailed) ||
	    ((_MDInMuskingumC1ID           = MFVarGetID (MDVarMuskingumC1,      MFNoUnit, MFInput,  MFState, MFBoundary)) == CMfailed) ||
	    ((_MDInMuskingumC2ID           = MFVarGetID (MDVarMuskingumC2,      MFNoUnit, MFInput,  MFState, MFBoundary)) == CMfailed) ||
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
            (MFModelAddFunction(_MDDischLevel3Muskingum) == CMfailed)) return (CMfailed);
	MFDefLeaving ("Discharge Muskingum");
	return (_MDOutDischLevel3ID);
}

