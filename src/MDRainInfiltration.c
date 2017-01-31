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
static int _MDInRainWaterSurplusConID     = MFUnset;
static int _MDInRainWaterSurplusDecID     = MFUnset;
static int _MDInRainWaterSurplusMixID     = MFUnset;
static int _MDInRainWaterSurplusLawnID    = MFUnset;
static int _MDInRainWaterSurplusAgID      = MFUnset;
static int _MDInRainWaterSurplusWetID     = MFUnset;
static int _MDInLandUseSpatialConID         = MFUnset;
static int _MDInLandUseSpatialDecID         = MFUnset;
static int _MDInLandUseSpatialMixID         = MFUnset;
static int _MDInLandUseSpatialAgID          = MFUnset;
static int _MDInLandUseSpatialLawnID        = MFUnset;
static int _MDInLandUseSpatialWetID         = MFUnset;
static int _MDInLandUseSpatialImpID         = MFUnset;
static int _MDInLandUseSpatialH2OID         = MFUnset;

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

	surplusC = MFVarGetFloat(_MDInRainWaterSurplusConID, itemID, 0.0);
	surplusD = MFVarGetFloat(_MDInRainWaterSurplusDecID, itemID, 0.0);
        surplusM = MFVarGetFloat(_MDInRainWaterSurplusMixID, itemID, 0.0);
	surfRunoff   = (surplusC + surplusD + surplusM) * (1.0 - _MDInfiltrationFrac);
	infiltration = (surplusC + surplusD + surplusM) *_MDInfiltrationFrac;
	MFVarSetFloat (_MDOutRainSurfRunoffID,       itemID, surfRunoff);
	MFVarSetFloat (_MDOutRainInfiltrationID,     itemID, infiltration);
 //  if ((itemID == 8640) || (itemID == 4596))   printf("id = %d, Gamma = %f, Infiltration %f surfRunoff %f \n",itemID, _MDInfiltrationFrac, infiltration,surfRunoff);
   //     if (itemID == 1) printf("Gamma = %f\n", _MDInfiltrationFrac);
}
static void _MDRainInfiltrationPnET (int itemID) {

	float surplusConif;
        float surplusDecid;
        float surplusMixed;
        float surplusLawn;
        float surplusAg;
        float surplusWet;
        float fracConif;
        float fracMixed;
        float fracDecid;
        float fracLawn;
        float fracAg;
        float fracWet;
        float fracImp;
        float fracH2O;
	float surfRunoff;
	float infiltration;
    
	if (_MDInfiltrationFractionID != MFUnset)
		_MDInfiltrationFrac = MFVarGetFloat(_MDInfiltrationFractionID,itemID,0.0);

	surplusConif = MFVarGetFloat(_MDInRainWaterSurplusConID, itemID, 0.0);
	surplusDecid = MFVarGetFloat(_MDInRainWaterSurplusDecID, itemID, 0.0);
        surplusMixed = MFVarGetFloat(_MDInRainWaterSurplusMixID, itemID, 0.0);
        surplusLawn = MFVarGetFloat(_MDInRainWaterSurplusLawnID, itemID, 0.0);
        surplusAg = MFVarGetFloat(_MDInRainWaterSurplusAgID,itemID, 0.0);
        surplusWet = MFVarGetFloat(_MDInRainWaterSurplusWetID,itemID,0.0);
        fracConif = MFVarGetFloat(_MDInLandUseSpatialConID, itemID,0.0);
        fracDecid = MFVarGetFloat(_MDInLandUseSpatialDecID, itemID, 0.0);
        fracMixed = MFVarGetFloat(_MDInLandUseSpatialMixID, itemID, 0.0);
        fracAg = MFVarGetFloat(_MDInLandUseSpatialAgID, itemID, 0.0)*0.01;
        fracLawn = MFVarGetFloat(_MDInLandUseSpatialLawnID, itemID, 0.0);
        fracWet = MFVarGetFloat(_MDInLandUseSpatialWetID, itemID, 0.0);
        fracImp = MFVarGetFloat(_MDInLandUseSpatialImpID, itemID, 0.0);
        fracH2O = MFVarGetFloat(_MDInLandUseSpatialH2OID, itemID, 0.0);
        float tmpLandBalance;
        // Temp store weighted runoff in output surfRunoff
        // Note: No water surplus from impervious or open water areas
        
	surfRunoff   = surplusConif*fracConif + surplusDecid*fracDecid + surplusMixed*fracMixed 
                        + surplusLawn*fracLawn + surplusAg*fracAg 
                        + surplusWet*fracWet;
	
        infiltration = surfRunoff *_MDInfiltrationFrac;
        surfRunoff = surfRunoff * (1.0 - _MDInfiltrationFrac); // Calculate output surface runoff
        if ((MFDateGetCurrentYear() > 0) && (infiltration != infiltration)) printf("itemID %d: conif %.2e decid %.2e mixed %.2e lawn %.2e ag %.2e wet %.2e\n",itemID,surplusConif, surplusDecid, surplusMixed, surplusLawn,surplusAg,surplusWet);
	MFVarSetFloat (_MDOutRainSurfRunoffID,       itemID, surfRunoff);
	MFVarSetFloat (_MDOutRainInfiltrationID,     itemID, infiltration);
 //  if ((itemID == 8640) || (itemID == 4596))   printf("id = %d, Gamma = %f, Infiltration %f surfRunoff %f \n",itemID, _MDInfiltrationFrac, infiltration,surfRunoff);
   //     if (itemID == 1) printf("Gamma = %f\n", _MDInfiltrationFrac);
}

static void _MDRainInfiltrationSaturation (int itemID){
		MFVarSetFloat (_MDOutRainSurfRunoffID,       itemID, MFVarGetFloat(_MDInSaturationExcessRunoffID, itemID,0.0));
		MFVarSetFloat (_MDOutRainInfiltrationID,     itemID, MFVarGetFloat(_MDOutRainInfiltrationID, itemID,0.0));
}

enum { MDinput, MDinput2, MDPnET, MDspatial, MDsimple, MDvarying};

int MDRainInfiltrationDef () {
	int  optID = MFUnset;
	int ret =0;
	const char *optStr, *optName = MDVarRainInfiltration;
	const char *options [] = { MDInputStr, MDInput2Str, MDPnETStr,"spatially", "simple", "varying",(char *) NULL };
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
                        if (((_MDInRainWaterSurplusDecID     = MFVarGetID (MDVarRainWaterSurplusDecid,    "mm", MFInput,  MFFlux, MFBoundary)) == CMfailed) ||
			    ((_MDInRainWaterSurplusConID     = MFVarGetID (MDVarRainWaterSurplusConif,    "mm", MFInput,  MFFlux, MFBoundary)) == CMfailed) ||
     			    ((_MDInRainWaterSurplusMixID     = MFVarGetID (MDVarRainWaterSurplusMixed,    "mm", MFInput,  MFFlux, MFBoundary)) == CMfailed) ||
                            ((_MDOutRainSurfRunoffID       = MFVarGetID (MDVarRainSurfRunoff,       "mm", MFOutput, MFFlux, MFBoundary)) == CMfailed) ||
			    ((_MDOutRainInfiltrationID     = MFVarGetID (MDVarRainInfiltration,     "mm", MFOutput, MFFlux, MFBoundary)) == CMfailed) ||
                            (MFModelAddFunction (_MDRainInfiltrationSimple2) == CMfailed)) return (CMfailed);	
                        break;
                case MDPnET:
                    if (_MDInfiltrationFractionID != MFUnset) {
                        if (((optStr = MFOptionGet (MDParInfiltrationFrac))  != (char *) NULL) &&
                                        (sscanf (optStr,"%f",&par) == 1))
                                            _MDInfiltrationFrac = par;
                                    else goto Stop;
                    }
                    if (((_MDInRainWaterSurplusDecID     = MFVarGetID (MDVarRainWaterSurplusDecid,    "mm", MFInput,  MFFlux, MFBoundary)) == CMfailed) ||
			    ((_MDInRainWaterSurplusConID     = MFVarGetID (MDVarRainWaterSurplusConif,    "mm", MFInput,  MFFlux, MFBoundary)) == CMfailed) ||
     			    ((_MDInRainWaterSurplusMixID     = MFVarGetID (MDVarRainWaterSurplusMixed,    "mm", MFInput,  MFFlux, MFBoundary)) == CMfailed) ||
                            ((_MDInRainWaterSurplusLawnID     = MFVarGetID (MDVarRainWaterSurplusLawn,    "mm", MFInput,  MFFlux, MFBoundary)) == CMfailed) ||
                            ((_MDInRainWaterSurplusAgID     = MFVarGetID (MDVarRainWaterSurplusAg,    "mm", MFInput,  MFFlux, MFBoundary)) == CMfailed) ||
                            ((_MDInRainWaterSurplusWetID        = MFVarGetID (MDVarRainWaterSurplusWet,     "mm", MFInput, MFFlux, MFBoundary)) == CMfailed) ||
                            ((_MDInLandUseSpatialConID        = MFVarGetID (MDVarLandUseSpatialCon,     "-",  MFInput, MFState, MFBoundary)) == CMfailed) ||
                            ((_MDInLandUseSpatialDecID        = MFVarGetID (MDVarLandUseSpatialDec,     "-",  MFInput, MFState, MFBoundary)) == CMfailed) ||
                            ((_MDInLandUseSpatialMixID        = MFVarGetID (MDVarLandUseSpatialMix,     "-",  MFInput, MFState, MFBoundary)) == CMfailed) ||
                            ((_MDInLandUseSpatialAgID        = MFVarGetID (MDVarLandUseSpatialAg,     "-",  MFInput, MFState, MFBoundary)) == CMfailed) ||
                            ((_MDInLandUseSpatialLawnID        = MFVarGetID (MDVarLawnFraction,     "-",  MFInput, MFState, MFBoundary)) == CMfailed) ||
                            ((_MDInLandUseSpatialWetID        = MFVarGetID (MDVarFracWetlandArea,     "-",  MFInput, MFState, MFBoundary)) == CMfailed) ||
                            ((_MDInLandUseSpatialImpID     = MFVarGetID (MDVarImpFracSpatial,     "-",  MFInput, MFState, MFBoundary)) == CMfailed) ||  
                            ((_MDInLandUseSpatialH2OID     = MFVarGetID (MDVarH2OFracSpatial,     "-",  MFInput, MFState, MFBoundary)) == CMfailed) ||  
                            ((_MDOutRainSurfRunoffID       = MFVarGetID (MDVarRainSurfRunoff,       "mm", MFOutput, MFFlux, MFBoundary)) == CMfailed) ||
			    ((_MDOutRainInfiltrationID     = MFVarGetID (MDVarRainInfiltration,     "mm", MFOutput, MFFlux, MFBoundary)) == CMfailed) ||
                            (MFModelAddFunction (_MDRainInfiltrationPnET) == CMfailed)) return (CMfailed);	
                        break;
                
		case MDspatial:
			if (( _MDInfiltrationFractionID = MFVarGetID (MDParInfiltrationFrac, "-", MFInput, MFState, MFBoundary)) == CMfailed) return (CMfailed);
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
