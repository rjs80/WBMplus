/******************************************************************************

GHAAS Water Balance/Transport Model V3.0
Global Hydrologic Archive and Analysis System
Copyright 1994-2011, UNH - CCNY/CUNY

MDSPackChg.c

balazs.fekete@unh.edu

 * Zuidema (9/2014) Updates to provide option of melt factor definition and include development.
 * TODO:
 *  Separate melt factor definition functions
 *  Provide options for: 1) Melt factor definition, 2) Including development snowpack
 *  Provide separate calculation for development snowpack.
*******************************************************************************/

#include <stdio.h>
#include <cm.h>
#include <math.h>
#include <MF.h>
#include <MD.h>

// Input
static int _MDInAtMeanID    = MFUnset;
static int _MDInPrecipID    = MFUnset;
static int _MDInWinterOnsetID = MFUnset;
static int _MDInForestCoverID = MFUnset;
static int _MDInImpFractionID = MFUnset;
static int _MDInHCIAID        =MFUnset;
// Output
static int _MDOutSnowPackID = MFUnset;
static int _MDOutSPackChgID = MFUnset;
static int _MDOutSnowMeltID = MFUnset;
static int _MDOutSnowFallID = MFUnset;
static int _MDOutSnowDensityID = MFUnset;
static int _MDOutSnowDepthID = MFUnset;
static int _MDCalculateSoilTempID = MFUnset;
static float _MDSnowMeltThreshold = 1.0;
static float _MDFallThreshold = -1.0;
static int _MDOutSnowPackDaysID = MFUnset;
static int _MDOutSnowMeltFactorID = MFUnset;
static int _MDOutSnowPackAgeID   = MFUnset;
static int _MDOutImpSnowFallROID = MFUnset;

static void _MDSPackChg (int itemID) {
// Input
	float airT;
	float precip;
	float winterOnsetDoY;
	float initialDensity = 150;
	float sDensitySlope = 3;
        float snowPackDays = 0.0;
	if( _MDCalculateSoilTempID==1)winterOnsetDoY = MFVarGetFloat(_MDInWinterOnsetID, itemID,1.0);
//	printf ("Anf SnowPackChange \n");
	// Local
	float sPack;
	float sPackChg = 0.0;
	float sDensity=0.0;
	float sDepth = 0.0;
	int snowAge=0;
        float packAge;
	float densityOfWater =1000;
	if (MFDateGetDayOfYear()  - winterOnsetDoY > 0){
	snowAge = MFDateGetDayOfYear() - winterOnsetDoY;
	}else{
		snowAge = 365 - winterOnsetDoY + MFDateGetDayOfYear();
	}
	
	sPack  = MFVarGetFloat (_MDOutSnowPackID, itemID, 0.0);
	if (MFVarTestMissingVal (_MDInAtMeanID,itemID) || MFVarTestMissingVal (_MDInPrecipID, itemID)) {
		MFVarSetFloat (_MDOutSnowFallID, itemID, 0.0);
		MFVarSetFloat (_MDOutSnowMeltID, itemID, 0.0);
		MFVarSetFloat (_MDOutSnowPackID, itemID, sPack);	
		MFVarSetFloat (_MDOutSPackChgID, itemID, 0.0);
		return; 
	}

//printf ("SnowMeltThreshold= %f SnoFall %f\n",_MDSnowMeltThreshold, _MDFallThreshold);
	airT   = MFVarGetFloat (_MDInAtMeanID,    itemID, 0.0);
	precip = MFVarGetFloat (_MDInPrecipID,    itemID, 0.0);

        packAge = MFVarGetFloat(_MDOutSnowPackAgeID,itemID,0.0);
        
	if (airT < _MDFallThreshold) {  /* Accumulating snow pack */
		MFVarSetFloat (_MDOutSnowFallID, itemID, precip);
		MFVarSetFloat (_MDOutSnowMeltID, itemID, 0.0);
		MFVarSetFloat (_MDOutSnowPackID, itemID, sPack + precip);
		MFVarSetFloat (_MDOutSPackChgID, itemID, precip);
                if (precip > 0.0) {
                    MFVarSetFloat(_MDOutSnowPackAgeID,itemID, 1.0); // New fallen Snow. 
                } else {
                    MFVarSetFloat(_MDOutSnowPackAgeID,itemID,packAge+1.0); // Age the snow.
	}
	}
	else if (airT > _MDSnowMeltThreshold) { /* Melting snow pack */
		sPackChg = MFVarGetFloat (_MDOutSnowMeltFactorID, itemID, 0.0); //2.63 + 2.55 * airT + 0.0912 * airT * precip;
		sPackChg = - (sPack < sPackChg ? sPack : sPackChg);
		MFVarSetFloat (_MDOutSnowFallID, itemID, 0.0);
		MFVarSetFloat (_MDOutSnowMeltID, itemID, fabs(sPackChg));
		MFVarSetFloat (_MDOutSPackChgID, itemID, sPackChg);
		MFVarSetFloat (_MDOutSnowPackID, itemID, sPack + sPackChg);
                if ((sPack + sPackChg) > 0.0) {
                    MFVarSetFloat(_MDOutSnowPackAgeID,itemID,packAge+1.0);
                } else {
                    MFVarSetFloat(_MDOutSnowPackAgeID,itemID,0.0);
	}
	}
	else { /* No change when air temperature is in [-1.0,1.0] range */
		MFVarSetFloat (_MDOutSnowFallID, itemID, 0.0);
		MFVarSetFloat (_MDOutSnowMeltID, itemID, 0.0);
		MFVarSetFloat (_MDOutSnowPackID, itemID, sPack);	
		MFVarSetFloat (_MDOutSPackChgID, itemID, 0.0);
                MFVarSetFloat (_MDOutSnowPackAgeID, itemID, packAge+1.0);
	}
	
		sDensity = (initialDensity + (snowAge * sDensitySlope));
		if (sPack > 0.0 ) sDepth = sPack  * densityOfWater / sDensity; //in mm
		//printf ("sAge %i sDens %f sPack %f sDepth %f \n", snowAge, sDensity, sPack, sDepth);
                if (sPack > 0.0) snowPackDays = 1;
                //if (itemID == 1) printf("SnowFall = %f, SnowMelt = %f\n", _MDFallThreshold, _MDSnowMeltThreshold);
                
		MFVarSetFloat(_MDOutSnowDensityID,itemID,sDensity);  
		MFVarSetFloat(_MDOutSnowDepthID,itemID, sDepth); 
		MFVarSetFloat(_MDOutSnowPackDaysID, itemID, snowPackDays); 
//	printf ("Ende SnowPackChange \n");

//		if (itemID == 54914) printf("-- m = %d, d = %d, SF = %f, SM = %f, airT = %f, precip = %f, sPack = %f, sPackChg = %f\n", MFDateGetCurrentMonth(), MFDateGetCurrentDay(), _MDFallThreshold, _MDSnowMeltThreshold, airT, precip, sPack, sPackChg);

}
void _MDDingmanMeltFactor (int itemID) {
    // Follows Male & Gray (1981),Gray & Prowse (1993) and Dingman (2002)
    float airT = MFVarGetFloat( _MDInAtMeanID,itemID,0.0); //2.63 + 2.55 * airT + 0.0912 * airT * precip;
    float precip = MFVarGetFloat(_MDInPrecipID,itemID,0.0);

    if (precip == 0) {
        // Male & Gray
        // TODO Should also incorporate slope factor! Assumed to be 1 at our resolution.
        float albedo;
        float F = MFVarGetFloat(_MDInForestCoverID,itemID,1.0);
        float age = MDMinimum(14.,MFVarGetFloat( _MDOutSnowPackAgeID,itemID,0.0));
        albedo = 0.7854-0.037946*age+0.0014*pow(age,2); // From Dingman 2002
        float fsl = 1.0;
        MFVarSetFloat(_MDOutSnowMeltFactorID,itemID,4.0*(1-albedo)*exp(-4.*F)*fsl*(airT-_MDSnowMeltThreshold));
    } else {
        MFVarSetFloat(_MDOutSnowMeltFactorID,itemID,(0.74 + 0.007 * precip )*(airT -_MDSnowMeltThreshold));
    }
    
}

void _MDWilmontMeltFactor (int itemID) {
    //Input 
    // TODO: Check citation (Wilmont) - should this be airT or (airT-meltT)?
    
    float airT = MFVarGetFloat( _MDInAtMeanID,itemID,0.0); //2.63 + 2.55 * airT + 0.0912 * airT * precip;
    float precip = MFVarGetFloat(_MDInPrecipID,itemID,0.0);
    MFVarSetFloat(_MDOutSnowMeltFactorID,itemID,2.63 + 2.55 * airT + 0.0912 * airT * precip);
}

void _MDDevelopedSnowFallRunoff (int itemID) {
    // Snow falling on impervious areas (or a percentage thereof?) is directed straight to runoff
    float hcia;
    
    float snowfall;
    float snowPackChange;
    float snowPack;
    float snowMelt;
    float impervSnowFallRunoff;
    hcia = MFVarGetFloat(_MDInImpFractionID,itemID,0.0) * MFVarGetFloat(_MDInHCIAID,itemID,0.0);
    snowfall = MFVarGetFloat(_MDOutSnowFallID,itemID,0.0);
    snowPack = MFVarGetFloat(_MDOutSnowPackID,itemID,0.0);
    snowPackChange = MFVarGetFloat(_MDOutSPackChgID,itemID,0.0);
    snowMelt = MFVarGetFloat(_MDOutSnowMeltID,itemID,0.0);
        
    // Calculate the volume of snow-falling on impervious area
    impervSnowFallRunoff = snowfall > 0.0 ? snowfall*hcia : 0.0;
    
    // Remove this portion from the (whole-cell-averaged) snow-pack SWE depth
    snowPack = MDMaximum(snowPack - impervSnowFallRunoff,0.0);
    
    // Add this volume to snowPackChange (as a whole-cell-averaged depth)
    snowPackChange = snowPackChange - impervSnowFallRunoff;
    //  Assign this value to snowMelt as well.
    snowMelt = fabs(snowPackChange);
    
    MFVarSetFloat(_MDOutSnowPackID,itemID,snowPack);
    MFVarSetFloat(_MDOutSPackChgID,itemID,snowPackChange);
    MFVarSetFloat(_MDOutSnowMeltID,itemID,snowMelt);
    MFVarSetFloat(_MDOutImpSnowFallROID,itemID,impervSnowFallRunoff);
}

int MDSPackChgDef () {

	if (_MDOutSPackChgID != MFUnset) return (_MDOutSPackChgID);
	MFDefEntering ("Snow Pack Change");
	const char *optStr;
	const char *soilTemperatureOptions [] = { "none", "calculate", (char *) NULL };
	int soilTemperatureID;
	float par;
	if (((optStr = MFOptionGet (MDParSnowMeltThreshold))  != (char *) NULL) && (sscanf (optStr,"%f",&par) == 1))
		_MDSnowMeltThreshold = par;
	
	if (((optStr = MFOptionGet (MDParFallThreshold))  != (char *) NULL) && (sscanf (optStr,"%f",&par) == 1))
		_MDFallThreshold= par;

	if (((optStr = MFOptionGet (MDOptSoilTemperature))  == (char *) NULL) || ((soilTemperatureID = CMoptLookup (soilTemperatureOptions, optStr, true)) == CMfailed)) {
					CMmsgPrint(CMmsgUsrError," Soil TemperatureOption not specifed! Options = 'none' or 'calculate'\n");
					return CMfailed;
				}

	if (soilTemperatureID == 1 ){
		_MDCalculateSoilTempID=1;
		if ((_MDInWinterOnsetID    = MFVarGetID (MDVarWinterOnsetDoy, "DoY", MFInput,  MFState, MFBoundary)) == CMfailed) return CMfailed;
	}

	if (((_MDInPrecipID       = MDPrecipitationDef ()) == CMfailed) ||
	    ((_MDInAtMeanID       = MFVarGetID (MDVarAirTemperature, "degC", MFInput,  MFState, MFBoundary)) == CMfailed) ||
	    ((_MDOutSnowFallID    = MFVarGetID (MDVarSnowFall,       "mm",   MFOutput, MFFlux,  MFBoundary)) == CMfailed) ||
	    ((_MDOutSnowMeltID    = MFVarGetID (MDVarSnowMelt,       "mm",   MFOutput, MFFlux,  MFBoundary)) == CMfailed) ||
	    ((_MDOutSnowDensityID = MFVarGetID (MDVarSnowDensity,    "mm",   MFOutput, MFState, MFBoundary)) == CMfailed) ||
	    ((_MDOutSnowDepthID   = MFVarGetID (MDVarSnowDepth,      "mm",   MFOutput, MFState, MFBoundary)) == CMfailed) ||
	    ((_MDOutSnowPackID    = MFVarGetID (MDVarSnowPack,       "mm",   MFOutput, MFState, MFInitial))  == CMfailed) ||
	    ((_MDOutSPackChgID    = MFVarGetID (MDVarSnowPackChange, "mm",   MFOutput, MFFlux,  MFBoundary)) == CMfailed) ||
	    ((_MDOutSnowPackDaysID= MFVarGetID (MDVarSnowPackDays, "-",      MFOutput, MFFlux,  MFBoundary)) == CMfailed) ||
            ((_MDOutSnowPackAgeID = MFVarGetID (MDVarSnowPackAge,     "d",      MFOutput,MFState,MFBoundary)) == CMfailed) ||
            ((_MDOutSnowMeltFactorID= MFVarGetID (MDVarSnowMeltFactor, "-",  MFOutput, MFState, MFBoundary)) == CMfailed)
            ) return (CMfailed);
        //Choose appropriate function for melt calculation.
        enum { MFcalculate, MFcalculate2 };
        const char *snowMeltCalculationOptions [] = {MDCalculateStr,MDCalculate2Str,(char *) NULL };
        int snowMeltCalcID;
        if ((optStr = MFOptionGet (MDOptSnowMeltCalculation) ) == (char *) NULL) {
            optStr = MDCalculateStr;            
            CMmsgPrint(CMmsgWarning," Snow Melt Factor Calculation method not specified - defaulting to Wilmont 1983.\n");
        } 
        if ((snowMeltCalcID = CMoptLookup (snowMeltCalculationOptions,optStr,true)) == CMfailed) {
            CMmsgPrint(CMmsgUsrError," Snow Melt Factor Calculation method incorrectly specified. Options are 'calculate' (Wilmont) or 'calculate2' (Dingman).\n");
            return (CMfailed);
        }
        
        switch (snowMeltCalcID) {
            case MFcalculate:
                if (MFModelAddFunction(_MDWilmontMeltFactor) == CMfailed) return (CMfailed);
                break;
            case MFcalculate2:
                if ( 
                    ((_MDInForestCoverID       = MFVarGetID (MDVarFracForestedArea, MFNoUnit, MFInput,MFState,MFBoundary)) == CMfailed) ||
                    (MFModelAddFunction(_MDDingmanMeltFactor) == CMfailed)
                    ) return (CMfailed);
                break;
            default: MFOptionMessage (MDOptSnowMeltCalculation, optStr, snowMeltCalculationOptions); return (CMfailed);
        }

        if (MFModelAddFunction (_MDSPackChg) == CMfailed) return (CMfailed);
        
        // Add in check for development ... then read-in Impervious and HCIA ... then add model function
        enum { MFnoImp, MFcalcImpRO };
        const char *snowImpervMeltOptions [] = { MDNoneStr, MDCalculateStr,(char *) NULL };
        int impervSnowMeltCalcID;
        if ((optStr = MFOptionGet (MDOptImperviousMeltCalc) ) == (char *) NULL)  {
            optStr = MDNoneStr;
            CMmsgPrint(CMmsgWarning," Impervious Snow Fall runoff method not specified - defaulting to none (SPackChgDef).\n");
        } 
        if ((impervSnowMeltCalcID = CMoptLookup (snowImpervMeltOptions,optStr,true)) == CMfailed) {
            CMmsgPrint(CMmsgUsrError," Impervious Snow Fall runoff method incorrectly specified.  Options are 'none' and 'calculate'.\n");
            return (CMfailed);
        }
        
        switch (impervSnowMeltCalcID) {
            case MFnoImp:
                // Nothing to be done
                break;
            case MFcalcImpRO:
                if ( 
                    ((_MDOutImpSnowFallROID    = MFVarGetID ( MDVarImpSnowFallRunoff,"mm",   MFOutput, MFFlux,  MFBoundary)) == CMfailed) ||
                    ((_MDInImpFractionID        = MFVarGetID (MDVarImpFracSpatial,  "-",    MFInput,  MFState, MFBoundary)) == CMfailed) ||
                    ((_MDInHCIAID               = MFVarGetID (MDVarHCIA,            "-",    MFInput,  MFState, MFBoundary)) == CMfailed) ||
                    (MFModelAddFunction(_MDDevelopedSnowFallRunoff) == CMfailed)
                    ) return (CMfailed);
                break;
            default: MFOptionMessage (MDOptImperviousMeltCalc, optStr, snowImpervMeltOptions); return (CMfailed);
        }
        
	MFDefLeaving ("Snow Pack Change");
	return (_MDOutSPackChgID);
}

int MDSPackMeltDef () {

	if (_MDOutSnowMeltID != MFUnset) return (_MDOutSnowMeltID);

	if ((MDSPackChgDef () == CMfailed) ||
	    ((_MDOutSnowMeltID   = MFVarGetID (MDVarSnowMelt, "mm", MFInput, MFFlux, MFBoundary)) == CMfailed))
		return (CMfailed);
	return (_MDOutSnowMeltID);
}

