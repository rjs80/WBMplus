/******************************************************************************

GHAAS Water Balance/Transport Model V3.0
Global Hydrologic Archive and Analysis System
Copyright 1994-2011, UNH - CCNY/CUNY

MDRiverWidth.c

balazs.fekete@unh.edu

*******************************************************************************/

#include <stdio.h>
#include <math.h>
#include <cm.h>
#include <MF.h>
#include <MD.h>

// Inputs
static int _MDInDischargeID             = MFUnset;
static int _MDInRiverbedShapeExponentID = MFUnset;
static int _MDInRiverbedAvgDepthMeanID  = MFUnset;
static int _MDInRiverbedWidthMeanID     = MFUnset;
static int _MDInRiverbedVelocityMeanID  = MFUnset;
static int _MDInAtaSiteWidthCoefID      = MFUnset;
static int _MDInAtaSiteDepthCoefID      = MFUnset;
static int _MDInAtaSiteVelocityCoefID   = MFUnset;
static int _MDInAtaSiteWidthExpID      = MFUnset;
static int _MDInAtaSiteDepthExpID      = MFUnset;
static int _MDInAtaSiteVelocityExpID   = MFUnset;

// Outputs
static int _MDOutRiverDepthID           = MFUnset;
static int _MDOutRiverWidthID           = MFUnset;
static int _MDOutRiverVelocityID        = MFUnset;

static void _MDRiverWidth (int itemID) {
// Input
	float discharge; // Discharge [m3/s]
	float shapeExp;  // Riverbed shape exponent
	float avgDepth;  // Average depth at mean discharge [m]
	float avgWidth;  // Average width at mean discharge [m]
	float velocity;  // Flow velocity [m/s]
// Output
	float depth;     // Flow depth at current discharge [m]
	float width;     // Flow width at current discharge [m]
// Local
	float alpha;     // Shape coefficient
	float area;      // Cross-section area [m2]
        float Width_Coef    = 0.0;
        float Depth_Coef    = 0.0;
        float Velocity_Coef = 0.0;
        float Width_Exp     = 0.0;
        float Depth_Exp     = 0.0;
        float Velocity_Exp  = 0.0;
        
	discharge = MFVarGetFloat (_MDInDischargeID,              itemID, 0.0);
	shapeExp  = MFVarGetFloat (_MDInRiverbedShapeExponentID,  itemID, 0.0);
	avgDepth  = MFVarGetFloat (_MDInRiverbedAvgDepthMeanID,   itemID, 0.0);
	avgWidth  = MFVarGetFloat (_MDInRiverbedWidthMeanID,      itemID, 0.0);
	velocity  = MFVarGetFloat (_MDInRiverbedVelocityMeanID,   itemID, 0.0);
        
        Width_Coef    = MFVarGetFloat (_MDInAtaSiteWidthCoefID,        itemID, 0.0);
        Depth_Coef    = MFVarGetFloat (_MDInAtaSiteWidthCoefID,        itemID, 0.0);
        Velocity_Coef = MFVarGetFloat (_MDInAtaSiteVelocityCoefID,     itemID, 0.0);
        Width_Exp     = MFVarGetFloat (_MDInAtaSiteWidthExpID,        itemID, 0.0);
        Depth_Exp     = MFVarGetFloat (_MDInAtaSiteWidthExpID,        itemID, 0.0);
        Velocity_Exp  = MFVarGetFloat (_MDInAtaSiteVelocityExpID,     itemID, 0.0);


	if (CMmathEqualValues (discharge, 0.0) ||
	    CMmathEqualValues (avgDepth,  0.0) ||
	    CMmathEqualValues (avgWidth,  0.0) ||
	    CMmathEqualValues (velocity,  0.0)) {
		width = depth = 0.0;
	}
	else	{
		
            /**************** Balaz's Approach*************/
        //      alpha = avgDepth / pow (avgWidth, shapeExp);
	//	area  = discharge / velocity;
	//	width = pow (((shapeExp + 1.0) * area) / (shapeExp * alpha), 1.0 / (shapeExp + 1));
	//	depth = alpha * pow (width, shapeExp);
                
            /**************** Wil's Approach **************/    
              
             width = Width_Coef * pow(discharge, Width_Exp);  
             depth = Depth_Coef * pow(discharge, Depth_Exp) / 100.;  
             velocity = Velocity_Coef * pow(discharge, Velocity_Exp);  
            
	}
       
	MFVarSetFloat (_MDOutRiverDepthID,   itemID, depth);
	MFVarSetFloat (_MDOutRiverWidthID,   itemID, width);
	MFVarSetFloat (_MDOutRiverVelocityID,   itemID, velocity);
}

int MDRiverWidthDef () {

	if (_MDOutRiverWidthID != MFUnset) return (_MDOutRiverWidthID);

	MFDefEntering ("River Geometry");

	if (((_MDInDischargeID             = MDDischargeDef ())             == CMfailed) ||
	    ((_MDInRiverbedShapeExponentID = MDRiverbedShapeExponentDef ()) == CMfailed) ||
	    ((_MDInRiverbedAvgDepthMeanID  = MFVarGetID (MDVarRiverbedAvgDepthMean,  "m",      MFInput,  MFState, MFBoundary)) == CMfailed) ||
            ((_MDInRiverbedWidthMeanID     = MFVarGetID (MDVarRiverbedWidthMean,     "m",      MFInput,  MFState, MFBoundary)) == CMfailed) ||
            ((_MDInRiverbedVelocityMeanID  = MFVarGetID (MDVarRiverbedVelocityMean,  "m/s",    MFInput,  MFState, MFBoundary)) == CMfailed) ||
            ((_MDInAtaSiteWidthCoefID      = MFVarGetID (MDVarAtaSiteWidthCoef,      "-",      MFInput,  MFState, MFBoundary)) == CMfailed) ||	    
            ((_MDInAtaSiteDepthCoefID      = MFVarGetID (MDVarAtaSiteDepthCoef,      "-",      MFInput,  MFState, MFBoundary)) == CMfailed) ||	    
            ((_MDInAtaSiteVelocityCoefID   = MFVarGetID (MDVarAtaSiteVelocityCoef,   "-",      MFInput,  MFState, MFBoundary)) == CMfailed) ||	    
            ((_MDInAtaSiteWidthExpID      = MFVarGetID (MDVarAtaSiteWidthExp,      "-",      MFInput,  MFState, MFBoundary)) == CMfailed) ||	    
            ((_MDInAtaSiteDepthExpID      = MFVarGetID (MDVarAtaSiteDepthExp,      "-",      MFInput,  MFState, MFBoundary)) == CMfailed) ||	    
            ((_MDInAtaSiteVelocityExpID   = MFVarGetID (MDVarAtaSiteVelocityExp,   "-",      MFInput,  MFState, MFBoundary)) == CMfailed) ||	    
            ((_MDOutRiverDepthID           = MFVarGetID (MDVarRiverDepth,            "m",      MFOutput, MFState, MFBoundary)) == CMfailed) ||
            ((_MDOutRiverWidthID           = MFVarGetID (MDVarRiverWidth,            "m",      MFOutput, MFState, MFBoundary)) == CMfailed) ||
            ((_MDOutRiverVelocityID        = MFVarGetID (MDVarRiverVelocity,         "m/s",    MFOutput, MFState, MFBoundary)) == CMfailed) ||
	    (MFModelAddFunction (_MDRiverWidth) == CMfailed)) return (CMfailed);
	MFDefLeaving ("River Geometry");

	return (_MDOutRiverWidthID);
}
