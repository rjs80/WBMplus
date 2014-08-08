/******************************************************************************

GHAAS Water Balance/Transport Model V3.0
Global Hydrologic Archive and Analysis System
Copyright 1994-2011, UNH - CCNY/CUNY

WBMMain.c

balazs.fekete@unh.edu

*******************************************************************************/
#include "wbm.h"

enum { MDpet, MDsurplus, MDinfiltration, MDrunoff, MDdischarge,  MDwatertemp, MDthermal, MDthermal2, MDthermal3, MDbalance, MDgeometry, MDbgc, MDbgc_DIN, MDbgc_DINPLUSBIOMASS, MDbgc_DOC, MDfecal, MDDO2, MDDIN, MDIrrigation, MDDOC,MDSpecCond,MDChloride};

int main (int argc,char *argv []) {
	int argNum;
	int  optID = MDbalance;
	const char *optStr, *optName = MDOptModel;
	const char *options [] = { "pet", "surplus", "infiltration", "runoff", "discharge",  "watertemp", "thermal", "thermal2", "thermal3", "balance", "geometry", "bgc", "bgc_DIN","bgc_DINPLUSBIOMASS", "bgc_DOC", "fecal", "DO2", "DIN", "irrigation", "DOC","SpecCond","chloride", (char *) NULL };

	argNum = MFOptionParse (argc,argv);

	if ((optStr = MFOptionGet (optName)) != (char *) NULL) optID = CMoptLookup (options, optStr, true);

	switch (optID) {
		case MDpet:                   return (MFModelRun (argc,argv,argNum,MDRainPotETDef));
		case MDsurplus:               return (MFModelRun (argc,argv,argNum,MDRainWaterSurplusDef));
		case MDinfiltration:          return (MFModelRun (argc,argv,argNum,MDRainInfiltrationDef));
		case MDrunoff:                return (MFModelRun (argc,argv,argNum,MDRunoffDef));
		case MDdischarge:             return (MFModelRun (argc,argv,argNum,MDDischargeDef));
		case MDbalance:               return (MFModelRun (argc,argv,argNum,MDWaterBalanceDef));
		case MDwatertemp:             return (MFModelRun (argc,argv,argNum,MDWTempRiverRouteDef));
		case MDthermal:               return (MFModelRun (argc,argv,argNum,MDThermalInputsDef));		// RJS 013112
		case MDthermal2:              return (MFModelRun (argc,argv,argNum,MDThermalInputs2Def));	// RJS 062012
		case MDthermal3:              return (MFModelRun (argc,argv,argNum,MDThermalInputs3Def));	// RJS 112712
		case MDgeometry:              return (MFModelRun (argc,argv,argNum,MDRiverWidthDef));
		case MDbgc:                   return (MFModelRun (argc,argv,argNum,MDBgcRoutingDef));
		case MDbgc_DOC:               return (MFModelRun (argc,argv,argNum,MDBgcDOCRoutingDef));
		case MDbgc_DIN:               return (MFModelRun (argc,argv,argNum,MDBgcDINRoutingDef));
		case MDbgc_DINPLUSBIOMASS:    return (MFModelRun (argc,argv,argNum,MDBgcDINPlusBiomassRoutingDef));
		case MDDO2:                   return (MFModelRun (argc,argv,argNum,MDDO2Def));		// RJS 111612
		case MDDIN:                   return (MFModelRun (argc,argv,argNum,MDDINDef));		// RJS 121013
                case MDIrrigation:            return (MFModelRun (argc,argv,argNum,MDIrrigationDef));   // RJS 121013
                case MDDOC:                   return (MFModelRun (argc,argv,argNum,MDDOCDef));          // RJS 011914
                case MDSpecCond:              return (MFModelRun (argc,argv,argNum,MDSpecCondDef));          // SZ 061014
                case MDChloride:              return (MFModelRun (argc,argv,argNum,MDChlorideDef));          // SZ 061614
		default: MFOptionMessage (optName, optStr, options); return (CMfailed);
	}
	return (CMfailed);
}
