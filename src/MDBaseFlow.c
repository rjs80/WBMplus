/******************************************************************************

GHAAS Water Balance/Transport Model V3.0
Global Hydrologic Archive and Analysis System
Copyright 1994-2011, UNH - CCNY/CUNY

MDBaseFlow.c

balazs.fekete@unh.edu

*******************************************************************************/

#include <stdio.h>
#include <math.h>
#include <cm.h>
#include <MF.h>
#include <MD.h>

// Input
static int _MDInRechargeID           = MFUnset;
static int _MDInIrrGrossDemandID     = MFUnset;
static int _MDInIrrReturnFlowID      = MFUnset;
static int _MDInIrrAreaFracID        = MFUnset;
static int _MDInSmallResReleaseID    = MFUnset;
static int _MDInGroundWatBETAID      = MFUnset;
// Output
static int _MDOutGrdWatID            = MFUnset;
static int _MDOutGrdWatChgID         = MFUnset;
static int _MDOutGrdWatRechargeID    = MFUnset;
static int _MDOutGrdWatUptakeID      = MFUnset;
static int _MDOutBaseFlowID          = MFUnset;
static int _MDOutIrrUptakeGrdWaterID = MFUnset;
static int _MDOutIrrUptakeExternalID = MFUnset;

static float _MDGroundWatBETA = 0.016666667;

static void _MDBaseFlow (int itemID) {
// Input
	float irrDemand;               // Irrigation demand [mm/dt]
	float irrReturnFlow;           // Irrigational return flow [mm/dt]
	float irrAreaFraction;         // Irrigated area fraction
// Output
	float grdWater;                // Groundwater size   [mm]
	float grdWaterChg;             // Groundwater change [mm/dt]
	float grdWaterRecharge;        // Groundwater recharge [mm/dt]
	float grdWaterUptake;          // Groundwater uptake [mm/dt]
	float baseFlow          = 0.0; // Base flow from groundwater [mm/dt]
	float irrUptakeGrdWater = 0.0; // Irrigational water uptake from shallow groundwater [mm/dt]
	float irrUptakeExt      = 0.0; // Unmet irrigational water demand [mm/dt]
// Local
                     
	grdWaterChg = grdWater = MFVarGetFloat (_MDOutGrdWatID,  itemID, 0.0);
	if (grdWater < 0.0) grdWaterChg = grdWater = 0.0;			//RJS 071511
	grdWaterRecharge = MFVarGetFloat (_MDInRechargeID, itemID, 0.0);
	grdWater = grdWater + grdWaterRecharge;

	if ((_MDInIrrGrossDemandID != MFUnset) &&
	    (_MDInIrrReturnFlowID  != MFUnset) &&
	    (_MDInIrrAreaFracID    != MFUnset) &&
		((irrAreaFraction   = MFVarGetFloat (_MDInIrrAreaFracID,   itemID, 0.0)) > 0.0)) {

		irrReturnFlow = MFVarGetFloat (_MDInIrrReturnFlowID,  itemID, 0.0);
		irrDemand     = MFVarGetFloat (_MDInIrrGrossDemandID, itemID, 0.0);

		grdWater         = grdWater         + irrReturnFlow;
		grdWaterRecharge = grdWaterRecharge + irrReturnFlow;
		//if (itemID==9)printf("Item %i, returnflow %f, demand %f, grdWater %f area %f \n",itemID,irrReturnFlow,irrDemand,grdWater,irrAreaFraction);
		if (_MDInSmallResReleaseID    != MFUnset) irrDemand = irrDemand - MFVarGetFloat(_MDInSmallResReleaseID,itemID,0.0);
		if (_MDOutIrrUptakeGrdWaterID != MFUnset) {
			if (irrDemand < grdWater) {
				// Irrigation demand is satisfied from groundwater storage 
				irrUptakeGrdWater = irrDemand;
				grdWater = grdWater - irrUptakeGrdWater;
			}
			else {
				// Irrigation demand needs external source
				irrUptakeGrdWater = grdWater;
				irrUptakeExt      = irrDemand - irrUptakeGrdWater;
				grdWater = 0.0;
			}
			MFVarSetFloat (_MDOutIrrUptakeGrdWaterID, itemID, irrUptakeGrdWater);
		}
		else irrUptakeExt = irrDemand;
		MFVarSetFloat (_MDOutIrrUptakeExternalID, itemID, irrUptakeExt);
	}
        _MDGroundWatBETA = _MDInGroundWatBETAID == MFUnset ? _MDGroundWatBETA : MFVarGetFloat(_MDInGroundWatBETAID, itemID, _MDGroundWatBETA);
	baseFlow    = grdWater * _MDGroundWatBETA;
	grdWater    = grdWater - baseFlow;
	grdWaterChg = grdWater - grdWaterChg;

	//if ((itemID == 486)) printf("y = %d, m = %d, d = %d, baseFlow = %f, grdWater = %f, grdWaterChg = %f, grdWaterRecharge = %f\n", MFDateGetCurrentYear(), MFDateGetCurrentMonth(), MFDateGetCurrentDay(), baseFlow, grdWater, grdWaterChg, grdWaterRecharge);	//RJS 071511

  //      if (itemID == 1) printf("Beta = %f\n", _MDGroundWatBETA);
        
//in= irrReturnFlow+grdWaterRecharge;
//out = irrUptakeGrdWater + grdWaterChg;
//float balance;
//balance = in -out;
//if (balance > 0.1)printf ("Balance %f in %f out %f rech%f \n",balance,in,out,grdWaterRecharge);
	grdWaterUptake = baseFlow + irrUptakeGrdWater;

	MFVarSetFloat (_MDOutGrdWatID,         itemID, grdWater);
    MFVarSetFloat (_MDOutGrdWatChgID,      itemID, grdWaterChg);
    MFVarSetFloat (_MDOutGrdWatRechargeID, itemID, grdWaterRecharge);
    MFVarSetFloat (_MDOutGrdWatUptakeID,   itemID, grdWaterUptake);
	MFVarSetFloat (_MDOutBaseFlowID,       itemID, baseFlow);
}

static void _MDBaseFlow2 (int itemID) {
// Input

// Output
	float grdWater;                // Groundwater size   [mm]
	float grdWaterChg;             // Groundwater change [mm/dt]
	float grdWaterRecharge;        // Groundwater recharge [mm/dt]
	float grdWaterUptake;          // Groundwater uptake [mm/dt]
	float baseFlow          = 0.0; // Base flow from groundwater [mm/dt]

// Local
                     
	grdWaterChg = grdWater = MFVarGetFloat (_MDOutGrdWatID,  itemID, 0.0);
	if (grdWater < 0.0) grdWaterChg = grdWater = 0.0;			//RJS 071511
	grdWaterRecharge = MFVarGetFloat (_MDInRechargeID, itemID, 0.0);
	grdWater = grdWater + grdWaterRecharge;

	baseFlow    = grdWater * _MDGroundWatBETA;
	grdWater    = grdWater - baseFlow;
	grdWaterChg = grdWater - grdWaterChg;

	
	MFVarSetFloat (_MDOutGrdWatID,         itemID, grdWater);
        MFVarSetFloat (_MDOutGrdWatChgID,      itemID, grdWaterChg);
        MFVarSetFloat (_MDOutGrdWatRechargeID, itemID, grdWaterRecharge);
	MFVarSetFloat (_MDOutBaseFlowID,       itemID, baseFlow);
}

enum { MDcalculate, MDinput2 , MDspatial};         // RJS 060214 // SZ 100114

int MDBaseFlowDef () {
	float par;
	const char *optStr;
        
        int  optID = MFUnset;                                                                                   // RJS 060214
	const char *optName = MDVarRunoff;                                                                      // RJS 060214
	const char *options [] = { MDCalculateStr, MDInput2Str, "spatially", (char *) NULL };                                // RJS 060214

	if (_MDOutBaseFlowID != MFUnset) return (_MDOutBaseFlowID);

	MFDefEntering ("Base flow");

	if ((optStr  = MFOptionGet (optName)) != (char *) NULL) optID = CMoptLookup (options, optStr, true);    // RJS 060214

        switch (optID) {
            case MDspatial:
                if ((_MDInGroundWatBETAID     = MFVarGetID(MDParGroundWatBETA,       "1/d", MFInput, MFState, MFBoundary)) == CMfailed) return (CMfailed);
             case MDcalculate: 	
                 if (_MDInGroundWatBETAID == MFUnset) {
                    if (((optStr = MFOptionGet (MDParGroundWatBETA))  != (char *) NULL) && (sscanf (optStr,"%f",&par) == 1)) _MDGroundWatBETA = par;
                 }
                if (((_MDInRechargeID       = MDRainInfiltrationDef ()) == CMfailed) ||
                    ((_MDInIrrGrossDemandID = MDIrrGrossDemandDef   ()) == CMfailed)) return (CMfailed);

                if ( _MDInIrrGrossDemandID != MFUnset) {
                        if (((_MDInSmallResReleaseID    = MDSmallReservoirReleaseDef ()) == CMfailed) ||
                            ((_MDInIrrAreaFracID        = MDIrrigatedAreaDef         ()) ==  CMfailed) ||
                            ((_MDInIrrReturnFlowID      = MFVarGetID (MDVarIrrReturnFlow,     "mm", MFInput,  MFFlux,  MFBoundary)) == CMfailed) ||
                            ((_MDOutIrrUptakeExternalID = MFVarGetID (MDVarIrrUptakeExternal, "mm", MFOutput, MFFlux,  MFBoundary)) == CMfailed) ||
                            ((_MDOutIrrUptakeGrdWaterID = MDIrrUptakeGrdWaterDef     ()) == CMfailed))
                            return CMfailed;
                 }
                if (((_MDOutGrdWatID                = MFVarGetID (MDVarGroundWater,         "mm", MFOutput, MFState, MFInitial))  == CMfailed) ||
                    ((_MDOutGrdWatChgID             = MFVarGetID (MDVarGroundWaterChange,   "mm", MFOutput, MFFlux,  MFBoundary)) == CMfailed) ||
                    ((_MDOutGrdWatRechargeID        = MFVarGetID (MDVarGroundWaterRecharge, "mm", MFOutput, MFFlux,  MFBoundary)) == CMfailed) ||
                    ((_MDOutGrdWatUptakeID          = MFVarGetID (MDVarGroundWaterUptake,   "mm", MFOutput, MFFlux,  MFBoundary)) == CMfailed) ||
                    ((_MDOutBaseFlowID              = MFVarGetID (MDVarBaseFlow,            "mm", MFOutput, MFFlux,  MFBoundary)) == CMfailed) ||
                     (MFModelAddFunction (_MDBaseFlow) == CMfailed)) return (CMfailed);
                break;
            case MDinput2:
                 if (((_MDInRechargeID               = MDRainInfiltrationDef ()) == CMfailed) ||
                     ((_MDOutGrdWatID                = MFVarGetID (MDVarGroundWater,         "mm", MFOutput, MFState, MFInitial))  == CMfailed) ||
                     ((_MDOutGrdWatChgID             = MFVarGetID (MDVarGroundWaterChange,   "mm", MFOutput, MFFlux,  MFBoundary)) == CMfailed) ||
                     ((_MDOutGrdWatRechargeID        = MFVarGetID (MDVarGroundWaterRecharge, "mm", MFOutput, MFFlux,  MFBoundary)) == CMfailed) ||
                     ((_MDOutGrdWatUptakeID          = MFVarGetID (MDVarGroundWaterUptake,   "mm", MFOutput, MFFlux,  MFBoundary)) == CMfailed) ||
                     ((_MDOutBaseFlowID              = MFVarGetID (MDVarBaseFlow,            "mm", MFOutput, MFFlux,  MFBoundary)) == CMfailed) ||
                     (MFModelAddFunction (_MDBaseFlow2) == CMfailed)) return (CMfailed);														// RJS 061312
                break;	      
        }
	MFDefLeaving ("Base flow ");
	return (_MDOutBaseFlowID);
}


 



