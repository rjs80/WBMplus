/******************************************************************************

GHAAS Water Balance/Transport Model V3.0
Global Hydrologic Archive and Analysis System
Copyright 1994-2011, UNH - CCNY/CUNY

MDWaterBalance.c

dominik.wisser@unh.edu
This is meant to check the vertical water balance for each grid cell. It does 
NOT include any water that is flowing laterally and should not be used to call BCG....
*******************************************************************************/

#include <stdio.h>
#include <math.h>
#include <cm.h>
#include <MF.h>
#include <MD.h>

// Input
static int _MDInPrecipID             = MFUnset;
static int _MDInEvaptrsID            = MFUnset;
static int _MDInSnowPackChgID        = MFUnset;
static int _MDInSoilMoistChgID       = MFUnset;
static int _MDInGrdWatChgID          = MFUnset;
static int _MDInRunoffID             = MFUnset;
static int _MDInDischargeID          = MFUnset;

static int _MDInIrrEvapotranspID     = MFUnset;
static int _MDInIrrSoilMoistChgID    = MFUnset;
static int _MDInIrrAreaFracID        = MFUnset;
static int _MDInIrrGrossDemandID     = MFUnset;
static int _MDInIrrReturnFlowID      = MFUnset;
static int _MDInIrrUptakeRiverID     = MFUnset;
static int _MDInIrrUptakeGrdWaterID  = MFUnset;
static int _MDInIrrUptakeExcessID    = MFUnset;
static int _MDInSmallResReleaseID    = MFUnset;
static int _MDInSmallResStorageChgID = MFUnset;
static int _MDInSmallResEvapoID      = MFUnset;

// Surface Runoff Pool
static int _MDInRunoffPoolChgID		 = MFUnset;
static int _MDInRunoffPoolID		 = MFUnset;
static int _MDInRunoffPoolRechargeID = MFUnset;
static int _MDInRunoffPoolReleaseID	 = MFUnset;

// Storm Runoff
static int _MDInStormRunoffImpID 	 = MFUnset;
static int _MDInStormRunoffH2OID 	 = MFUnset;
static int _MDInPrecipPervID  		 = MFUnset;
static int _MDInStormRunoffTotalID   = MFUnset;
static int _MDInRunofftoPervID 	 	 = MFUnset;

// Soil Moisture
static int _MDInEvaptrsNotScaledID      = MFUnset;
static int _MDInSoilMoistNotScaledID    = MFUnset;
static int _MDInSoilMoistID        	    = MFUnset;
static int _MDInSoilAvailWaterCapID     = MFUnset;
static int _MDInImpFractionID           = MFUnset;
static int _MDInH2OFractionID		    = MFUnset;
static int _MDInExcessID			    = MFUnset;
static int _MDInExcessNotScaledID	    = MFUnset;
static int _MDInSoilMoistChgNotScaledID = MFUnset;

// Rain Infiltration
static int _MDInRainWaterSurplusID  = MFUnset;
static int _MDInRainSurfRunoffID   = MFUnset;
static int _MDInRainInfiltrationID = MFUnset;

// Groundwater
static int _MDInGrdWatID                = MFUnset;
static int _MDInGrdWatRechargeID        = MFUnset;
static int _MDInBaseFlowID              = MFUnset;

// Snowpack
static int _MDInSnowFallID    = MFUnset;
static int _MDInSnowMeltID    = MFUnset;
static int _MDInSnowPackID    = MFUnset;

// Runoff
static int _MDInTotalSurfRunoffID  = MFUnset;


// Output
static int _MDOutWaterBalanceID      = MFUnset;
static int _MDOutIrrUptakeBalanceID  = MFUnset;
static int _MDOutIrrWaterBalanceID   = MFUnset;

static void _MDWaterBalance(int itemID) {
// Input
	float precip        	 = MFVarGetFloat(_MDInPrecipID,         itemID, 0.0);
	float etp           	 = MFVarGetFloat(_MDInEvaptrsID,        itemID, 0.0);
	float snowPackChg   	 = MFVarGetFloat(_MDInSnowPackChgID,    itemID, 0.0);
	float soilMoistChg  	 = MFVarGetFloat(_MDInSoilMoistChgID,   itemID, 0.0);
	float grdWaterChg   	 = MFVarGetFloat(_MDInGrdWatChgID,      itemID, 0.0);
	float runoff        	 = MFVarGetFloat(_MDInRunoffID,         itemID, 0.0);

	float runoffPoolChg 	 = MFVarGetFloat(_MDInRunoffPoolChgID,  itemID, 0.0);
	float runoffPool    	 = MFVarGetFloat(_MDInRunoffPoolID,  itemID, 0.0);
	float runoffPoolRecharge = MFVarGetFloat(_MDInRunoffPoolRechargeID,  itemID, 0.0);
	float runoffPoolRelease  = MFVarGetFloat(_MDInRunoffPoolReleaseID,  itemID, 0.0);

	float awCap        	 	 = MFVarGetFloat(_MDInSoilAvailWaterCapID,        itemID, 0.0);
	float h2oAreaFrac		 = MFVarGetFloat(_MDInH2OFractionID,              itemID, 0.0);
	float impAreaFrac		 = MFVarGetFloat(_MDInImpFractionID,              itemID, 0.0);
	float snowfall			 = MFVarGetFloat(_MDInSnowFallID,              	  itemID, 0.0);
	float snowmelt			 = MFVarGetFloat(_MDInSnowMeltID,                 itemID, 0.0);
	float snowpack			 = MFVarGetFloat(_MDInSnowPackID,                 itemID, 0.0);

	float stormRunoffImp			 = MFVarGetFloat(_MDInStormRunoffImpID,           itemID, 0.0);
	float stormRunoffH2O			 = MFVarGetFloat(_MDInStormRunoffH2OID,           itemID, 0.0);
	float stormRunoffTotal			 = MFVarGetFloat(_MDInStormRunoffTotalID,         itemID, 0.0);
	float precipPerv				 = MFVarGetFloat(_MDInPrecipPervID,               itemID, 0.0);
	float runoffToPerv				 = MFVarGetFloat(_MDInRunofftoPervID,             itemID, 0.0);

	float soilMoist					 = MFVarGetFloat(_MDInSoilMoistID,          itemID, 0.0);
	float soilMoistNotScaled		 = MFVarGetFloat(_MDInSoilMoistNotScaledID, itemID, 0.0);
	float etpNotScaled				 = MFVarGetFloat(_MDInEvaptrsNotScaledID,   itemID, 0.0);
	float excess					 = MFVarGetFloat(_MDInExcessID,             itemID, 0.0);
	float excessNotScaled			 = MFVarGetFloat(_MDInExcessNotScaledID,    itemID, 0.0);
	float soilMoistChgNotScaled		 = MFVarGetFloat(_MDInSoilMoistChgNotScaledID, itemID, 0.0);

	float surplus 					 = MFVarGetFloat(_MDInRainWaterSurplusID,   itemID, 0.0);
	float surfRunoff				 = MFVarGetFloat(_MDInRainSurfRunoffID,     itemID, 0.0);
	float infiltration				 = MFVarGetFloat(_MDInRainInfiltrationID,   itemID, 0.0);

	float grdWater					 = MFVarGetFloat(_MDInGrdWatID,             itemID, 0.0);
	float grdWaterRecharge			 = MFVarGetFloat(_MDInGrdWatRechargeID,     itemID, 0.0);
	float baseFlow					 = MFVarGetFloat(_MDInBaseFlowID,           itemID, 0.0);

	float totalSurfRunoff			 = MFVarGetFloat(_MDInTotalSurfRunoffID,    itemID, 0.0);

	float irrAreaFrac        = 0.0;
	float irrGrossDemand     = 0.0;
	float irrReturnFlow      = 0.0;
	float irrEvapotransp     = 0.0;
	float irrSoilMoistChg    = 0.0;
	float irrUptakeGrdWater  = 0.0;
	float irrUptakeRiver     = 0.0;
	float irrUptakeExcess    = 0.0;
	float smallResStorageChg = 0.0;
	float smallResRelease    = 0.0;
	float smallResEvapo      = 0.0;


// Output
	float balance;
	float balance2;
	float balance2b;

	if (_MDInIrrGrossDemandID != MFUnset) { 
		irrAreaFrac       = MFVarGetFloat (_MDInIrrAreaFracID,            itemID, 0.0);
		irrGrossDemand    = MFVarGetFloat (_MDInIrrGrossDemandID,         itemID, 0.0);
		irrReturnFlow     = MFVarGetFloat (_MDInIrrReturnFlowID,          itemID, 0.0);
		irrEvapotransp    = MFVarGetFloat (_MDInIrrEvapotranspID,         itemID, 0.0);
		irrSoilMoistChg   = _MDInIrrSoilMoistChgID   != MFUnset ? MFVarGetFloat (_MDInIrrSoilMoistChgID,   itemID, 0.0) : 0.0;
		irrUptakeGrdWater = _MDInIrrUptakeGrdWaterID != MFUnset ? MFVarGetFloat (_MDInIrrUptakeGrdWaterID, itemID, 0.0) : 0.0;
		irrUptakeRiver    = _MDInIrrUptakeRiverID    != MFUnset ? MFVarGetFloat (_MDInIrrUptakeRiverID,    itemID, 0.0) : 0.0;
		irrUptakeExcess   = MFVarGetFloat (_MDInIrrUptakeExcessID,        itemID, 0.0);

		if (_MDInSmallResReleaseID != MFUnset) {
			smallResRelease    = MFVarGetFloat (_MDInSmallResReleaseID,    itemID, 0.0);
			smallResStorageChg = MFVarGetFloat (_MDInSmallResStorageChgID, itemID, 0.0);
			smallResEvapo      = MFVarGetFloat (_MDInSmallResEvapoID,      itemID, 0.0);
		}

		balance = (precip - snowPackChg) * irrAreaFrac + irrGrossDemand - irrEvapotransp - irrSoilMoistChg - irrReturnFlow;
		MFVarSetFloat (_MDOutIrrWaterBalanceID, itemID, balance);

		balance = irrGrossDemand - (irrUptakeGrdWater + irrUptakeRiver + irrUptakeExcess + smallResRelease);
		MFVarSetFloat (_MDOutIrrUptakeBalanceID, itemID, balance);
	}
	balance  = precip + irrUptakeRiver + irrUptakeExcess - (etp + runoff + grdWaterChg + snowPackChg + soilMoistChg + smallResStorageChg);
	balance2b = (precip + irrUptakeRiver + irrUptakeExcess - (etp + runoff + grdWaterChg + snowPackChg + soilMoistChg + smallResStorageChg + runoffPoolChg)) / (etp + runoff + grdWaterChg + snowPackChg + soilMoistChg + smallResStorageChg + runoffPoolChg);
	balance2 = precip + irrUptakeRiver + irrUptakeExcess - (etp + runoff + grdWaterChg + snowPackChg + soilMoistChg + smallResStorageChg + runoffPoolChg);
//	printf("d = %d, m = %d, y = %d, waterbalance = %f\n", MFDateGetCurrentDay(), MFDateGetCurrentMonth(), MFDateGetCurrentYear(), balance2);

//	if (itemID == 125) {
	if (fabs (balance2) > 0.0001 ) {
//		printf ("TIEM %i %d %d %d WaterBalance! %f precip %f etp = %f runoff = %f grdWaterChg = %f snowPackChg = %f soilMoistChg = %f runoffPoolChg %f\n runoffPool = %f, runoffPoolRecharge = %f, runoffPoolRelease = %f\n", itemID, MFDateGetCurrentMonth(), MFDateGetCurrentDay(), MFDateGetCurrentYear(), balance2,precip,etp,runoff,grdWaterChg,snowPackChg,soilMoistChg,runoffPoolChg,runoffPool, runoffPoolRecharge, runoffPoolRelease);
		printf ("%i, %d, %d, %d, %f, %f, %f, %f, %f, %f,", itemID, MFDateGetCurrentYear(), MFDateGetCurrentMonth(), MFDateGetCurrentDay(), balance2b, balance2, awCap, impAreaFrac, h2oAreaFrac, precip);
		printf ("%f, %f, %f, %f,", snowPackChg, snowfall, snowmelt, snowpack);
		printf ("%f, %f, %f, %f, %f,", stormRunoffImp, stormRunoffH2O, stormRunoffTotal, precipPerv, runoffToPerv);
		printf ("%f, %f, %f, %f, %f, %f, %f, %f,", soilMoistChg, soilMoistChgNotScaled, soilMoist, soilMoistNotScaled, etp, etpNotScaled, excess, excessNotScaled);
		printf ("%f, %f, %f,", surplus, surfRunoff, infiltration);
		printf ("%f, %f, %f, %f,", grdWater, grdWaterChg, grdWaterRecharge, baseFlow);
		printf ("%f, %f, %f, %f,", runoffPool, runoffPoolChg, runoffPoolRecharge, runoffPoolRelease);
		printf ("%f, %f\n", totalSurfRunoff, runoff);
	}

	MFVarSetFloat (_MDOutWaterBalanceID, itemID , balance);
}

int MDWaterBalanceDef() {
 
	MFDefEntering ("WaterBalance");
	if ((                                  MDAccumBalanceDef     ()  == CMfailed) ||
	    ((_MDInPrecipID                  = MDPrecipitationDef    ()) == CMfailed) ||
	    ((_MDInDischargeID               = MDDischargeDef        ()) == CMfailed) ||
	 
	    ((_MDInSnowPackChgID             = MDSPackChgDef         ()) == CMfailed) ||
	    ((_MDInSoilMoistChgID            = MDSoilMoistChgDef     ()) == CMfailed) ||
	    ((_MDInRunoffID                  = MDRunoffDef           ()) == CMfailed) ||
	    ((_MDInEvaptrsID                 = MFVarGetID (MDVarEvapotranspiration,      	"mm",   MFInput,  MFFlux,  MFBoundary)) == CMfailed) ||
	    ((_MDInGrdWatChgID               = MFVarGetID (MDVarGroundWaterChange,       	"mm",   MFInput,  MFFlux,  MFBoundary)) == CMfailed) ||
	    ((_MDInRunoffPoolChgID           = MFVarGetID (MDVarRunoffPoolChg,           	"mm",   MFInput,  MFFlux,  MFBoundary)) == CMfailed) ||
	    ((_MDInRunoffPoolID              = MFVarGetID (MDVarRunoffPool,         		"mm", MFInput, MFState, MFInitial))  == CMfailed) ||
	    ((_MDInRunoffPoolRechargeID      = MFVarGetID (MDVarRunoffPoolRecharge, 		"mm", MFInput, MFFlux,  MFBoundary)) == CMfailed) ||
	    ((_MDInRunoffPoolReleaseID       = MFVarGetID (MDVarRunoffPoolRelease,  		"mm", MFInput, MFFlux,  MFBoundary)) == CMfailed) ||
	    ((_MDInStormRunoffImpID    		 = MFVarGetID (MDVarStormRunoffImp,  			"mm",   MFInput, MFFlux,  MFBoundary)) == CMfailed) ||
	    ((_MDInStormRunoffH2OID    		 = MFVarGetID (MDVarStormRunoffH2O,  			"mm",   MFInput, MFFlux,  MFBoundary)) == CMfailed) ||
	    ((_MDInPrecipPervID        		 = MFVarGetID (MDVarPrecipPerv,      			"mm",   MFInput, MFFlux, MFBoundary)) == CMfailed) ||
	    ((_MDInStormRunoffTotalID  		 = MFVarGetID (MDVarStormRunoffTotal,			"mm",   MFInput, MFFlux, MFBoundary))  == CMfailed) ||
	    ((_MDInRunofftoPervID      		 = MFVarGetID (MDVarRunofftoPerv,    			"mm",   MFInput, MFFlux, MFBoundary))  == CMfailed) ||
 	    ((_MDInEvaptrsID         		 = MFVarGetID (MDVarRainEvapotranspiration,     "mm",   MFInput, MFFlux,  MFBoundary)) == CMfailed) ||
 	    ((_MDInEvaptrsNotScaledID		 = MFVarGetID (MDVarRainETnotScaled,            "mm",   MFInput, MFFlux,  MFBoundary)) == CMfailed) ||	// RJS 082812
	    ((_MDInSoilMoistNotScaledID		 = MFVarGetID (MDVarRainSoilMoistureNotScaled,  "mm",   MFInput, MFState, MFInitial))  == CMfailed) ||
	    ((_MDInSoilMoistChgNotScaledID	 = MFVarGetID (MDVarRainSoilMoistureChangeNotScaled,  "mm",   MFInput, MFState, MFInitial))  == CMfailed) ||
	    ((_MDInSoilMoistID       		 = MFVarGetID (MDVarRainSoilMoisture,           "mm",   MFInput, MFState, MFBoundary)) == CMfailed) ||
	    ((_MDInSoilAvailWaterCapID 	 	 = MFVarGetID (MDVarSoilAvailWaterCap,          "mm",   MFInput,  MFFlux, MFInitial)) == CMfailed) ||
 	    ((_MDInExcessID           		 = MFVarGetID (MDVarExcess,                     "mm",   MFInput, MFFlux,  MFBoundary)) == CMfailed) ||	// RJS 091813
 	    ((_MDInExcessNotScaledID    	 = MFVarGetID (MDVarExcessNotScaled,            "mm",   MFInput, MFFlux,  MFBoundary)) == CMfailed) ||	// RJS 091813
 	    ((_MDInImpFractionID       		 = MFVarGetID (MDVarImpFracSpatial,             "-",    MFInput,  MFState, MFBoundary)) == CMfailed) ||   // RJS 082812
	    ((_MDInH2OFractionID       	 	 = MFVarGetID (MDVarH2OFracSpatial,             "-",    MFInput,  MFState, MFBoundary)) == CMfailed) ||   // RJS 082812
		((_MDInRainWaterSurplusID  		 = MFVarGetID (MDVarRainWaterSurplus,           "mm",    MFInput, MFFlux,  MFBoundary)) == CMfailed) ||
		((_MDInRainSurfRunoffID          = MFVarGetID (MDVarRainSurfRunoff,       		"mm", MFInput, MFFlux, MFBoundary)) == CMfailed) ||
		((_MDInRainInfiltrationID        = MFVarGetID (MDVarRainInfiltration,     		"mm", MFInput, MFFlux, MFBoundary)) == CMfailed) ||
		((_MDInGrdWatID                	 = MFVarGetID (MDVarGroundWater,         		"mm", MFInput, MFState, MFInitial))  == CMfailed) ||
		((_MDInGrdWatRechargeID        	 = MFVarGetID (MDVarGroundWaterRecharge, 		"mm", MFInput, MFFlux,  MFBoundary)) == CMfailed) ||
		((_MDInBaseFlowID              	 = MFVarGetID (MDVarBaseFlow,            		"mm", MFInput, MFFlux,  MFBoundary)) == CMfailed) ||
	    ((_MDInSnowFallID    			 = MFVarGetID (MDVarSnowFall,       			"mm",   MFInput, MFFlux,  MFBoundary)) == CMfailed) ||
	    ((_MDInSnowMeltID   			 = MFVarGetID (MDVarSnowMelt,       			"mm",   MFInput, MFFlux,  MFBoundary)) == CMfailed) ||
	    ((_MDInSnowPackID   			 = MFVarGetID (MDVarSnowPack,       			"mm",   MFInput, MFState, MFInitial))  == CMfailed) ||
		((_MDInTotalSurfRunoffID 		 = MFVarGetID (MDVarTotalSurfRunoff, 			"mm",   MFInput, MFFlux,  MFBoundary)) == CMfailed) ||
	    ((_MDOutWaterBalanceID           = MFVarGetID (MDVarWaterBalance,       		"mm",   MFOutput, MFFlux,  MFBoundary)) == CMfailed) ||
	    (MFModelAddFunction(_MDWaterBalance) == CMfailed))
	    return (CMfailed);
	if ((_MDInIrrGrossDemandID           = MDIrrGrossDemandDef    ()) != MFUnset) {
		if ((_MDInIrrGrossDemandID == CMfailed) ||
	        ((_MDInIrrUptakeRiverID      = MDIrrUptakeRiverDef    ()) == CMfailed) ||
	        ((_MDInIrrUptakeGrdWaterID   = MDIrrUptakeGrdWaterDef ()) == CMfailed) ||
	        ((_MDInIrrSoilMoistChgID     = MDIrrSoilMoistChgDef   ()) == CMfailed) ||
	        ((_MDInIrrAreaFracID         = MDIrrigatedAreaDef    ())==  CMfailed) ||
	        ((_MDInIrrEvapotranspID      = MFVarGetID (MDVarIrrEvapotranspiration, "mm",   MFInput,  MFFlux,  MFBoundary)) == CMfailed) ||
	        ((_MDInIrrReturnFlowID       = MFVarGetID (MDVarIrrReturnFlow,         "mm",   MFInput,  MFFlux,  MFBoundary)) == CMfailed) || 
	        ((_MDInIrrUptakeExcessID     = MFVarGetID (MDVarIrrUptakeExcess,       "mm",   MFInput,  MFFlux,  MFBoundary)) == CMfailed) ||
		    ((_MDOutIrrUptakeBalanceID   = MFVarGetID (MDVarIrrUptakeBalance,      "mm",   MFOutput, MFFlux,  MFBoundary)) == CMfailed) ||
		    ((_MDOutIrrWaterBalanceID    = MFVarGetID (MDVarIrrWaterBalance,       "mm",   MFOutput, MFFlux,  MFBoundary)) == CMfailed))
	    	return (CMfailed);		
		if ((_MDInSmallResReleaseID        = MDSmallReservoirReleaseDef ()) != MFUnset) {
			if (( _MDInSmallResReleaseID == CMfailed) ||
			    ((_MDInSmallResEvapoID      = MFVarGetID (MDVarSmallResEvaporation,   "mm", MFInput, MFFlux,  MFBoundary)) == CMfailed) ||
			    ((_MDInSmallResStorageChgID = MFVarGetID (MDVarSmallResStorageChange, "mm", MFInput, MFState, MFInitial))  == CMfailed))
			    return (CMfailed);
		}
	}
	MFDefLeaving ("WaterBalance");
	return (_MDOutWaterBalanceID);	
}
