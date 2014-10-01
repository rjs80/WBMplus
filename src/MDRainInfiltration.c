/******************************************************************************

GHAAS Water Balance/Transport Model V3.0
Global Hydrologic Archive and Analysis System
Copyright 1994-2011, UNH - CCNY/CUNY

MDInfiltration.c

balazs.fekete@unh.edu

*******************************************************************************/

#include <stdio.h>
#include <string.h>
#include <cm.h>
#include <MF.h>
#include <MD.h>

// Input
static int _MDInRainWaterSurplusID      = MFUnset;
static int _MDInRainWaterSurplusCID     = MFUnset;
static int _MDInRainWaterSurplusDID     = MFUnset;
static int _MDInRainWaterSurplusMID     = MFUnset;

// Output
static int _MDOutRainSurfRunoffID       = MFUnset;
static int _MDOutRainInfiltrationID     = MFUnset;
static float _MDInfiltrationFrac = 0.5;
static int  _MDInfiltrationFractionID   = MFUnset;
static int _MDInSaturationExcessRunoffID = MFUnset;
static int _MDInRainInfiltrationID = MFUnset;
static void _MDRainInfiltrationSimple (int itemID) {

	float surplus;
	float surfRunoff;
	float infiltration;

	if (_MDInfiltrationFractionID != MFUnset)
		_MDInfiltrationFrac = MFVarGetFloat(_MDInfiltrationFractionID,itemID,0.0);

	surplus = MFVarGetFloat(_MDInRainWaterSurplusID, itemID, 0.0);
	surfRunoff   = surplus * (1.0 - _MDInfiltrationFrac);
	infiltration = surplus *_MDInfiltrationFrac;
	MFVarSetFloat (_MDOutRainSurfRunoffID,       itemID, surfRunoff);
	MFVarSetFloat (_MDOutRainInfiltrationID,     itemID, infiltration);
   //     printf("Gamma = %f, Infiltraction %f surfRunoff %f \n",_MDInfiltrationFrac, infiltration,surfRunoff);
   //     if (itemID == 1) printf("Gamma = %f\n", _MDInfiltrationFrac);
}

static void _MDRainInfiltrationSimple2 (int itemID) {

	float surplusC;
        float surplusD;
        float surplusM;
	float surfRunoff;
	float infiltration;

	if (_MDInfiltrationFractionID != MFUnset)
		_MDInfiltrationFrac = MFVarGetFloat(_MDInfiltrationFractionID,itemID,0.0);

	surplusC = MFVarGetFloat(_MDInRainWaterSurplusCID, itemID, 0.0);
	surplusD = MFVarGetFloat(_MDInRainWaterSurplusDID, itemID, 0.0);
        surplusM = MFVarGetFloat(_MDInRainWaterSurplusMID, itemID, 0.0);
	surfRunoff   = (surplusC + surplusD + surplusM) * (1.0 - _MDInfiltrationFrac);
	infiltration = (surplusC + surplusD + surplusM) *_MDInfiltrationFrac;
	MFVarSetFloat (_MDOutRainSurfRunoffID,       itemID, surfRunoff);
	MFVarSetFloat (_MDOutRainInfiltrationID,     itemID, infiltration);
 //  if ((itemID == 8640) || (itemID == 4596))   printf("id = %d, Gamma = %f, Infiltration %f surfRunoff %f \n",itemID, _MDInfiltrationFrac, infiltration,surfRunoff);
   //     if (itemID == 1) printf("Gamma = %f\n", _MDInfiltrationFrac);
}

static void _MDRainInfiltrationSaturation (int itemID){
		MFVarSetFloat (_MDOutRainSurfRunoffID,       itemID, MFVarGetFloat(_MDInSaturationExcessRunoffID, itemID,0.0));
		MFVarSetFloat (_MDOutRainInfiltrationID,     itemID, MFVarGetFloat(_MDOutRainInfiltrationID, itemID,0.0));
}

enum { MDinput, MDinput2, MDsimple, MDvarying, MDspatial };

int MDRainInfiltrationDef () {
	int  optID = MFUnset;
	int ret =0;
	const char *optStr, *optName = MDVarRainInfiltration;
	const char *options [] = { MDInputStr, MDInput2Str, "simple", "varying","spatially" ,(char *) NULL };
	float par;
	//printf ("THE framework = greatest time sink ever invented\n");
	if (_MDOutRainInfiltrationID != MFUnset) return (_MDOutRainInfiltrationID);

	if (((optStr = MFOptionGet (MDParInfiltrationFrac))  != (char *) NULL) && (sscanf (optStr,"%f",&par) == 1)) _MDInfiltrationFrac = par;		//RJS 082812, Gamma wasn't read in until this edit
	
	const char *soilMoistureOptions [] = { "bucket", "layers", (char *) NULL };
		int soilMoistureOptionID;
		 //TODO Add baseflow from layered SM to infiltration!
			if (((optStr = MFOptionGet (MDOptSoilMoisture))  == (char *) NULL) || ((soilMoistureOptionID = CMoptLookup (soilMoistureOptions, optStr, true)) == CMfailed)) {
						CMmsgPrint(CMmsgUsrError," Soil Moisture mode not specifed! Options = 'bucket' or 'layers'\n");
						return CMfailed;
					}
		
	
	MFDefEntering ("Rainfed Infiltration");
	 
	if ((optStr = MFOptionGet (optName)) != (char *) NULL) optID = CMoptLookup (options, optStr, true);
	
	//if ((ret = MDPermafrostDef())                == CMfailed) return CMfailed;
	
	if (soilMoistureOptionID ==1){ //layer is the soil moisture option infiltration will be calculated differently. 
		//if ((_MDInRainWaterSurplusID = MDRainWaterSurplusDef ()) == CMfailed) return (CMfailed);
		
		
		if ((ret = MDRainSMoistChgLayeredSoilDef()) == CMfailed) return CMfailed;
		if ((_MDOutRainSurfRunoffID       = MFVarGetID (MDVarRainSurfRunoff,       "mm", MFOutput, MFFlux, MFBoundary)) == CMfailed)return CMfailed; 
		if ((_MDOutRainInfiltrationID     = MFVarGetID (MDVarRainInfiltration,           "mm", MFOutput, MFFlux, MFBoundary)) == CMfailed) return CMfailed;
		if ((_MDInSaturationExcessRunoffID       = MFVarGetID (MDVarSaturationExcessflow,       "mm", MFInput, MFFlux, MFBoundary)) == CMfailed)return CMfailed; 
		if ((_MDInRainInfiltrationID     = MFVarGetID (MDVarRainInfiltration,           "mm", MFInput, MFFlux, MFBoundary)) == CMfailed) return CMfailed;
		
		
		if (MFModelAddFunction (_MDRainInfiltrationSaturation) == CMfailed) return CMfailed;
		MFDefLeaving  ("Rainfed Infiltration");
		return (_MDOutRainInfiltrationID);

	}
	
	switch (optID) {
		case MDinput:
			_MDOutRainInfiltrationID = MFVarGetID (MDVarRainInfiltration, "mm", MFInput, MFState, MFBoundary);
			break;
                case MDinput2:
                    if (_MDInfiltrationFractionID != MFUnset) {
				if (((optStr = MFOptionGet (MDParInfiltrationFrac))  != (char *) NULL) &&
				    (sscanf (optStr,"%f",&par) == 1))
					_MDInfiltrationFrac = par;
				else goto Stop;
			}
         //           if (((_MDInRainWaterSurplusID          = MFVarGetID (MDVarRainWaterSurplus,     "mm", MFInput,  MFFlux, MFBoundary)) == CMfailed) ||
                        if (((_MDInRainWaterSurplusDID     = MFVarGetID (MDVarRainWaterSurplusD,    "mm", MFInput,  MFFlux, MFBoundary)) == CMfailed) ||
			    ((_MDInRainWaterSurplusCID     = MFVarGetID (MDVarRainWaterSurplusC,    "mm", MFInput,  MFFlux, MFBoundary)) == CMfailed) ||
     			    ((_MDInRainWaterSurplusMID     = MFVarGetID (MDVarRainWaterSurplusM,    "mm", MFInput,  MFFlux, MFBoundary)) == CMfailed) ||
                            ((_MDOutRainSurfRunoffID       = MFVarGetID (MDVarRainSurfRunoff,       "mm", MFOutput, MFFlux, MFBoundary)) == CMfailed) ||
			    ((_MDOutRainInfiltrationID     = MFVarGetID (MDVarRainInfiltration,     "mm", MFOutput, MFFlux, MFBoundary)) == CMfailed) ||
                            (MFModelAddFunction (_MDRainInfiltrationSimple2) == CMfailed)) return (CMfailed);	
                        break;
		case MDspatial:
			if (( _MDInfiltrationFractionID = MFVarGetID (MDParInfiltrationFrac, "mm", MFInput, MFState, MFBoundary)) == CMfailed) return (CMfailed);
			//break;		// RJS 082812 // SZ commented out 09292014
		case MDsimple:
		case MDvarying:
			if ((_MDInRainWaterSurplusID = MDRainWaterSurplusDef ()) == CMfailed) return (CMfailed);
			if (_MDInfiltrationFractionID == MFUnset) {
				if (((optStr = MFOptionGet (MDParInfiltrationFrac))  != (char *) NULL) &&
				    (sscanf (optStr,"%f",&par) == 1))
					_MDInfiltrationFrac = par;
				else goto Stop;
			}
			if ((_MDOutRainSurfRunoffID       = MFVarGetID (MDVarRainSurfRunoff,       "mm", MFOutput, MFFlux, MFBoundary)) == CMfailed) return CMfailed;
			if ((_MDOutRainInfiltrationID     = MFVarGetID (MDVarRainInfiltration,     "mm", MFOutput, MFFlux, MFBoundary)) == CMfailed) return CMfailed;
			if  (MFModelAddFunction (_MDRainInfiltrationSimple) == CMfailed) return (CMfailed);
			break;
		default: MFOptionMessage (optName, optStr, options); return (CMfailed);
	}
	MFDefLeaving  ("Rainfed Infiltration");
	return (_MDOutRainInfiltrationID);
Stop:
	MFOptionMessage (optName, optStr, options);
	return (CMfailed);
}
