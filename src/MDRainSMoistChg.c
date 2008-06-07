/******************************************************************************

GHAAS Water Balance/Transport Model V3.0
Global Hydrologic Archive and Analysis System
Copyright 1994-2007, University of New Hampshire

MDRainSMoistChg.c

balazs.fekete@unh.edu

*******************************************************************************/

#include <stdio.h>
#include <math.h>
#include <cm.h>
#include <MF.h>
#include <MD.h>

static float _MDAWCap;
static float _MDSoilMoistALPHA = 5.0;

static float _MDDryingFunc (float sMoist) {
	float gm;
	gm = (1.0 - exp (- _MDSoilMoistALPHA * sMoist / _MDAWCap)) / (1.0 - exp (- _MDSoilMoistALPHA)); 
	return (gm);
}

// Input
static int _MDInAirTMeanID          = MFUnset;
static int _MDInPrecipID            = MFUnset;
static int _MDInPotETID             = MFUnset;
static int _MDInInterceptID         = MFUnset;
static int _MDInSPackChgID          = MFUnset;
static int _MDInSoilAvailWaterCapID = MFUnset;
static int _MDInIrrAreaFracID       = MFUnset;
 

// Output
static int _MDOutEvaptrsID          = MFUnset;
static int _MDOutSoilMoistCellID    = MFUnset;
static int _MDOutSoilMoistID        = MFUnset;
static int _MDOutSMoistChgID        = MFUnset;

static void _MDRainSMoistChg (int itemID) {	
// Input
	float airT;              // Air temperature [degreeC]
	float precip;            // Precipitation [mm/dt]
	float pet;               // Potential evapotranspiration [mm/dt]
	float intercept   = 0.0; // Interception (when the interception module is turned on) [mm/dt]
	float sPackChg;          // Snow pack change [mm/dt]
	float irrAreaFrac = 0.0; // Irrigated area fraction
//	TODO float impAreaFrac = 0.0; // Impervious area fraction RJS 01-17-08
//	TODO float H2OAreaFrac = 0.0; // water area fraction RJS 01-17-08
//	float runofftoPerv;      // runoff from impervious to pervious [mm/dt]  RJS 01-17-08
	float def;               // water deficit [mm/dt]
	float sMoist      = 0.0; // Soil moisture [mm/dt]
// Output
	float sMoistChg   = 0.0; // Soil moisture change [mm/dt]
	float transp      = 0.0; // Transpiration [mm]
	float excess      = 0.0;
	float evapotrans;
	float balance;
// Local
	float waterIn;

	airT         = MFVarGetFloat (_MDInAirTMeanID,          itemID, 0.0);
	precip       = MFVarGetFloat (_MDInPrecipID,            itemID, 0.0);
 	sPackChg     = MFVarGetFloat (_MDInSPackChgID,          itemID, 0.0);
	pet          = MFVarGetFloat (_MDInPotETID,             itemID, 0.0);
	_MDAWCap     = MFVarGetFloat (_MDInSoilAvailWaterCapID, itemID, 0.0);
	sMoist       = MFVarGetFloat (_MDOutSoilMoistCellID,    itemID, 0.0);
	intercept    = _MDInInterceptID   != MFUnset ? MFVarGetFloat (_MDInInterceptID,   itemID, 0.0) : 0.0;
	irrAreaFrac  = _MDInIrrAreaFracID != MFUnset ? MFVarGetFloat (_MDInIrrAreaFracID, itemID, 0.0) : 0.0;
	waterIn      = precip - intercept - sPackChg;

	pet = pet > intercept ? pet - intercept : 0.0;

	if ((airT > 0.0) && (_MDAWCap > 0.0)) {
		if (waterIn >= pet) {
			sMoistChg = waterIn - pet > _MDAWCap - sMoist ? waterIn - pet : _MDAWCap - sMoist;
		}
		else {
			sMoistChg = (-1.0 * _MDDryingFunc (sMoist)) * (pet - waterIn);
		}
		if (sMoistChg > _MDAWCap - sMoist) sMoistChg = _MDAWCap - sMoist;
		if (sMoistChg < 0.0) sMoistChg = 0.0;
		sMoist = sMoist + sMoistChg;
	}
	else  sMoistChg = 0.0;		
		
	if (precip - intercept - sPackChg - sMoistChg > pet) evapotrans = pet;
	else evapotrans = precip - intercept - sPackChg - sMoistChg;

	excess = precip - sPackChg - intercept - evapotrans - sMoistChg;

	balance = waterIn - intercept - evapotrans - sMoistChg - excess;
	if ((fabs (balance) > 0.001) && (_MDAWCap > 0.0)) printf ("balance = %f sMoist = %f, precip=%f sMoistChg = %f, transp = %f, waterIn = %f, pet = %f, excess = %f, def = %f itemID = %i AWC %f\n", balance, sMoist, precip, sMoistChg, transp, waterIn, pet, excess, def,itemID, _MDAWCap);

	MFVarSetFloat (_MDOutSoilMoistCellID, itemID, sMoist);
	MFVarSetFloat (_MDOutEvaptrsID,       itemID, evapotrans * (1.0 - irrAreaFrac)); //RJS 01-17-08 "- impAreaFrac - H2OAreaFrac"
	MFVarSetFloat (_MDOutSoilMoistID,     itemID, sMoist     * (1.0 - irrAreaFrac)); //RJS 01-17-08 "- impAreaFrac - H2OAreaFrac"
	MFVarSetFloat (_MDOutSMoistChgID,     itemID, sMoistChg  * (1.0 - irrAreaFrac)); //RJS 01-17-08 "- impAreaFrac - H2OAreaFrac"
}

int MDRainSMoistChgDef () {
	int ret = 0;
	float par;
	const char *optStr;
	if (_MDOutSMoistChgID != MFUnset) return (_MDOutSMoistChgID);

	MFDefEntering ("Rainfed Soil Moisture");
	if (((optStr = MFOptionGet (MDParSoilMoistALPHA))  != (char *) NULL) && (sscanf (optStr,"%f",&par) == 1)) _MDSoilMoistALPHA = par;

	if ((ret = MDIrrGrossDemandDef ()) == CMfailed) return (CMfailed);
	if ((ret != MFUnset)  &&
	  ((_MDInIrrAreaFracID         = MDIrrigatedAreaDef    ())==  CMfailed) )
		return (CMfailed);
	if (((_MDInPrecipID            = MDPrecipitationDef     ()) == CMfailed) ||
	    ((_MDInSPackChgID          = MDSPackChgDef          ()) == CMfailed) ||
	    ((_MDInPotETID             = MDRainPotETDef         ()) == CMfailed) ||
	    ((_MDInInterceptID         = MDRainInterceptDef     ()) == CMfailed) ||
	    ((_MDInSoilAvailWaterCapID = MDSoilAvailWaterCapDef ()) == CMfailed) ||
	    ((_MDInAirTMeanID          = MFVarGetID (MDVarAirTemperature,             "degC", MFInput,  MFState, MFBoundary)) == CMfailed) ||
 	    ((_MDOutEvaptrsID          = MFVarGetID (MDVarRainEvapotranspiration,     "mm",   MFOutput, MFFlux,  MFBoundary)) == CMfailed) ||
	    ((_MDOutSoilMoistCellID    = MFVarGetID (MDVarRainSoilMoistureCell,       "mm",   MFOutput, MFState, MFInitial))  == CMfailed) ||
	    ((_MDOutSoilMoistID        = MFVarGetID (MDVarRainSoilMoisture,           "mm",   MFOutput, MFState, MFBoundary)) == CMfailed) ||
        ((_MDOutSMoistChgID        = MFVarGetID (MDVarRainSoilMoistChange,        "mm",   MFOutput, MFState, MFBoundary)) == CMfailed) ||
        (MFModelAddFunction (_MDRainSMoistChg) == CMfailed)) return (CMfailed);
	MFDefLeaving ("Rainfed Soil Moisture");
	return (_MDOutSMoistChgID);
}
