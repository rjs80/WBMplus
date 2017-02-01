/******************************************************************************

GHAAS Water Balance/Transport Model V3.0
Global Hydrologic Archive and Analysis System
Copyright 1994-2011, UNH - CCNY/CUNY

MDDichRouteCascadeCoeff.c

shan.zuidema@unh.edu

*******************************************************************************/

#include <stdio.h>
#include <math.h>
#include <cm.h>
#include <MF.h>
#include <MD.h>

// Input
static int _MDInRiverbedVelocityMeanID  = MFUnset;
static int _MDInRiverbedShapeExponentID = MFUnset;
static int _MDInSinuosityID             = MFUnset;
// Output
static int _MDOutCascadeC0ID          = MFUnset;

static void _MDDischRouteCascadeCoeff (int itemID) {
// Input
	float vMean;            // Mean velocity
        float Si;               // Sinuosity
// Output
	float C0;               // Cascading reservoirs Storage C0 coefficient
// Local
	float dt;               // time step length [s]
	float dL;               // Cell length [m]
	
	dL        = MFModelGetLength (itemID);
        Si        = MFVarGetFloat (_MDInSinuosityID,            itemID, 0.0);	
	vMean     = 2.18 * (1.024 -0.077*log(1e-06*MFModelGetArea(itemID)/0.5))/3.6; // Wave velocity (m/s) from same source as Alex (wbm_trans).  Other numbers include 0.4 (m/s))
        // SZ added factor of 10 for testing ... 2016-12-15.  Cascade routing too strong of an effect - trying to reduce the dampening
	if (CMmathEqualValues (vMean,     0.0)) {
            // Can't have mean velocity defining a zero storage coefficient.
            // Need to default to something ... 1.0=accumulation, 0.95 is fairly reasonable
            //  but need to use an end-member. (C0=0 is alright for time-varying parameter definition)
		MFVarSetFloat (_MDOutCascadeC0ID, itemID, 1.0); 
		return;
	}
	if (CMmathEqualValues (dL,        0.0)) { 
	    // Falling back to flow-accumulation
		MFVarSetFloat (_MDOutCascadeC0ID, itemID, 1.0); 
		return;
	}
	dt = MFModelGet_dt ();
        
	C0 = (isnan(vMean) || vMean < 0.0001) ? 1. : 1. / ( 1. + (Si*dL/vMean/dt) ); // After Alex P's method.
	MFVarSetFloat (_MDOutCascadeC0ID, itemID, C0);
}

enum { MDinput, MDstatic };

int MDDischLevel3CascadeCoeffDef() {
	int  optID = MFUnset;
	const char *optStr, *optName = MDOptCascade;
	const char *options [] = { MDInputStr, "static", (char *) NULL };

	if (_MDOutCascadeC0ID != MFUnset) return (_MDOutCascadeC0ID);

	MFDefEntering ("Cascading Reservoirs Coefficients");
	if ((optStr = MFOptionGet (optName)) != (char *) NULL) optID = CMoptLookup (options, optStr, true);
	switch (optID) {
		case MDinput:
			if (((_MDOutCascadeC0ID = MFVarGetID (MDVarMuskingumC0, MFNoUnit,   MFInput,  MFState,  MFBoundary)) == CMfailed) )
				return (CMfailed); 
			break;
		case MDstatic:
			if (((_MDInRiverbedShapeExponentID  = MDRiverbedShapeExponentDef ()) == CMfailed) ||
			    ((_MDInSinuosityID              = MFVarGetID (MDVarSinuosity,  "m/m",  MFInput,    MFState,    MFBoundary )) == CMfailed ) ||
			    ((_MDInRiverbedVelocityMeanID   = MFVarGetID (MDVarRiverbedVelocityMean, "m/s",    MFInput,  MFState, MFBoundary)) == CMfailed) ||
			    ((_MDOutCascadeC0ID             = MFVarGetID (MDVarCascadeC0,          MFNoUnit, MFOutput, MFState, MFBoundary)) == CMfailed) ||
			    (MFModelAddFunction (_MDDischRouteCascadeCoeff) == CMfailed)) return (CMfailed);
			break;
		// TODO: add option to calculate coefficients based on time varying at-a-station relations.
                default: MFOptionMessage (optName, optStr, options); return (CMfailed);
	}
	MFDefLeaving ("Cascading Reservoirs Coefficients");
	return (_MDOutCascadeC0ID);
}
