/******************************************************************************
GHAAS Water Balance/Transport Model V3.0
Global Hydrologic Archive and Analysis System
Copyright 1994-2007, University of New Hampshire

MDNEP.c

Ken.r.sheehan@gmail.com

Module for GPP and R

*******************************************************************************/

#include <stdio.h>
#include <string.h>
#include <cm.h>
#include <MF.h>
#include <MD.h>
#include <math.h>

static int _MDInDischargeID                                     = MFUnset; // m3/s
static int _MDInRiverWidthID                                    = MFUnset; // m3/s
static int _MDInPARBenthicID					= MFUnset; // MJ/m2/d converted to mol/m2/day 
static int _MDInSCALERBasinID                                   = MFUnset;
static int _MDInRunoffVolumeID                                  = MFUnset;
static int _MDInCanopyCover2ID                                  = MFUnset;
static int _MDInOrderSwitchID                                   = MFUnset;
static int _MDInRiverOrderID                                    = MFUnset;
static int _MDInWTemp_QxTID                                     = MFUnset;
static int _MDInTrefID                                          = MFUnset;

static int _MDOutNEPID                                          = MFUnset; // kg/d NEP
static int _MDOutGPPID                                          = MFUnset;
static int _MDOutRkgdID                                         = MFUnset;
static int _MDOutGPPgm2dID                                      = MFUnset;
static int _MDOutRgm2dID                                        = MFUnset;
static int _MDOutDryCID                                         = MFUnset;
static int _MDOutDOCID                                          = MFUnset;
static int _MDOutTotalCID                                       = MFUnset;

static void _MDNEP (int itemID) {

    float GPP_gm2d           = 0.0;
    float GPP_gm2d_ref       = 0.0;
    float GPP_kgd            = 0.0;
    float GPP_mgL            = 0.0;
    float R_a_gm2d           = 0.0;
    float R_a_gm2d_ref       = 0.0;
    float R_h_gm2d           = 0.0;
    float R_total_gm2d       = 0.0;
    float R_total_kgd        = 0.0;
    float R_total_mgL        = 0.0;
    float NEP_gm2d           = 0.0;
    float NEP_kgd            = 0.0;
    float NEP_mgL            = 0.0;
    float Q_out              = 0.0;
    
    float PAR_benthic        = 0.0;
    float Width              = 0.0;
    float Length             = 0.0;
    float water_Riv_total_in = 0.0;
    float SCALER_BasinID     = 0.0;
    float C_Lit              = 0.0; // (g/m2/yr) Carbon litterfall 
    float C_Lat              = 0.0; // (g/m/yr) Carbon lateral movement/blow-in
    float DOC_mgL            = 0.0; // mg/L DOC entering via runoff
    float DOC_kgd            = 0.0;
    float totalC_kgd         = 0.0;
    float runoffVol          = 0.0;
    float dryC_kgd           = 0.0;
    float canopy_cover       = 0.0;
    float C_Lit_kgd          = 0.0;
    float C_Lat_kgd          = 0.0;
    float orderSwitch        = 0.0;
    float riverOrder         = 0.0;
    float waterT             = 0.0;
    float Tref               = 0.0;
    
    	Q_out			= MFVarGetFloat (_MDInDischargeID,       itemID, 0.0);                  // cubic meters/second m3/sec - Discharge leaving this grid cell after routing.
        water_Riv_total_in      = Q_out;
        Width                   = MFVarGetFloat (_MDInRiverWidthID, itemID, 0.0);			// meters
        PAR_benthic             = MFVarGetFloat (_MDInPARBenthicID, itemID, 0.0) * 2.0513 / 1000000;    // reading in MJ/m2/d and converting to moles/m2/d
        Length			= MFModelGetLength(itemID);									// meters
        SCALER_BasinID          = MFVarGetFloat (_MDInSCALERBasinID, itemID, 0.0);                      // BNZ = 1, CWT = 2, KNZ = 3, LUQ = 4, MER = 5, TLK = 6
        runoffVol               = MFVarGetFloat (_MDInRunoffVolumeID,    itemID, 0.0);                  // cubic meters per second m3/sec
        canopy_cover            = MFVarGetFloat (_MDInCanopyCover2ID, itemID, 0.0);
        orderSwitch             = MFVarGetFloat (_MDInOrderSwitchID, itemID, 0.0);
        riverOrder              = MFVarGetFloat (_MDInRiverOrderID, itemID, 0.0);
  	waterT  		= MFVarGetFloat (_MDInWTemp_QxTID,       itemID, 0.0); 		// degrees celsius 
  	Tref                    = MFVarGetFloat (_MDInTrefID,       itemID, 0.0); 		// degrees celsius 

 //       if (MFDateGetCurrentYear() > 2000) printf("Width = %f, Q = %f\n", Width, Q_out);
      
        if (riverOrder >= orderSwitch) {
      if (water_Riv_total_in >= 0.000001) {
          
        if (SCALER_BasinID == 1.0) {          // BNZ
        GPP_gm2d_ref = 0.04 * pow(PAR_benthic, 0.93);           
        R_a_gm2d_ref = 3.8;
        C_Lit = 165;            // g/m2/yr
        C_Lat = 5.0;            // g/m/yr
        DOC_mgL = 6.1;          // mg/L
        }
    
        else if (SCALER_BasinID == 2.0) {     // CWT
        GPP_gm2d_ref = 0.21 * pow(PAR_benthic, 0.097);           
        R_a_gm2d_ref = 3.6 * pow(PAR_benthic, 0.125);
        C_Lit = 350;            // g/m2/yr
        C_Lat = 35.5;           // g/m/yr
        DOC_mgL = 3.0;          // mg/L
        }
    
        else if (SCALER_BasinID == 3.0) {     // KNZ
        GPP_gm2d_ref = 0.02 * pow(PAR_benthic, 0.96);           // was 0.0589, to the 0.955
        R_a_gm2d_ref = 1.3                           ;          // was 0.13 to the 0.8
        C_Lit = 347.0;           // g/m2/yr
        C_Lat = 184.5;           // g/m/yr
        DOC_mgL = 2.5;           // mg/L
        }

        else if (SCALER_BasinID == 4.0) {     // LUQ
        GPP_gm2d_ref = 1.440 * pow(PAR_benthic, 0.72);           
        R_a_gm2d_ref = 4.700 * pow(PAR_benthic, 0.09);
        C_Lit = 350.0;           // g/m2/yr
        C_Lat = 35.5;            // g/m/yr
        DOC_mgL = 3.8;           // mg/L
        }
    
        else if (SCALER_BasinID == 5.0) {     // MER
        GPP_gm2d_ref = 0.1930 * pow(PAR_benthic, 0.68);           
        R_a_gm2d_ref = 5.7000 * pow(PAR_benthic, 0.059);
        C_Lit = 235.0;           // g/m2/yr
        C_Lat = 35.5;            // g/m/yr
        DOC_mgL = 4.3;           // mg/L
        }
    
        else if (SCALER_BasinID == 6.0) {     // TLK
        GPP_gm2d_ref = 0.43 * pow(PAR_benthic, 0.30);           
        R_a_gm2d_ref = 2.7;
        C_Lit = 0.0;             // g/m2/yr
        C_Lat = 2500.0;          // g/m/yr
        DOC_mgL = 6.4;           // mg/L       
        }
  
 //      printf("%d-%d-%d, itemID = %d, PAR_benthic = %f, GPP = %f, R = %f\n", MFDateGetCurrentYear(), MFDateGetCurrentMonth(), MFDateGetCurrentDay(), itemID, PAR_benthic, GPP_gm2d_ref, R_a_gm2d_ref);          
        
        C_Lit_kgd                = (C_Lit / 365) * Width * Length * canopy_cover / 1000;
        C_Lat_kgd                = (C_Lat / 365) * 2 * Length / 1000;
        
        DOC_kgd                  = DOC_mgL * runoffVol * 86400 / 1000;
        dryC_kgd                 = C_Lit_kgd + C_Lat_kgd;
        totalC_kgd               = dryC_kgd + DOC_kgd; 
       
	GPP_gm2d		 = GPP_gm2d_ref;
	R_a_gm2d		 = R_a_gm2d_ref; 

        
       if ((SCALER_BasinID == 1) || (SCALER_BasinID == 6) || (SCALER_BasinID == 3)) {                        // for TLK and BNZ and KNZ
 //     GPP_gm2d                = GPP_gm2d_ref * pow(2, ((waterT - Tref) / 10));
        R_a_gm2d                = R_a_gm2d_ref * pow(2, ((waterT - Tref) / 10));
       }
        
// Ken's GPP and R equations //
        
 //       GPP_gm2d     /*KNZ*/	= ((20.1641 * tanh((0.01887* (PAR_benthic))/20.1641))/1000 * 1440);    // g / m2 / day (chlorophyll_A * Photosynthesis_max + param);DO in milligrams per meter squared per day = now using julian et al, 2008
//        GPP_gm2d     /*CWT*/	= (0.61 * PAR_benthic) + 0.43;                                          //*/((20.1641 * tanh((0.01887* (PAR_benthic))/20.1641))/1000 * 1440); // (0.61 * PAR_benthic) + 0.43; //0.1;	// CWT uses second equation, KNZ the first. 					// g / m2 / day (chlorophyll_A * Photosynthesis_max + param);DO in milligrams per meter squared per day = now using julian et al, 2008
// 	R_a_gm2d     /*KNZ*/	= -1 * (-0.0014212 * (PAR_benthic) + 7.9221131)/1000 * 1440 ;       // KNZ
//	R_a_gm2d     /*CWT*/	= -1.34 * log(GPP_gm2d) -3.94; // -1 * (-0.0014212 * (PAR_benthic) + 7.9221131)/1000 * 1440 ;       // CWT
 
        GPP_kgd			= GPP_gm2d * (Length * Width) / 1000;                                     // kg / day oygen per grid cell of GPP based on Roberts equation
	GPP_mgL			= GPP_kgd / water_Riv_total_in /86400 * 1000 ;                          // mg / L
       
//	R_a_gm2d		= -1.34 * log(GPP_gm2d) -3.06;                                          // */-1 * (-0.0014212 * (PAR_benthic) + 7.9221131)/1000 * 1440 ;   // R from SCALER data, for CWT use : -1.34 * log(GPP_gm2d) -3.06;  CWT Walker Branch Equation **/Matt Tre g / m2 / day respirationtman n from autotrophs
	R_h_gm2d		= 0.0;									// g / m2 / day respiration from heterotrophs
	R_total_gm2d		= -1 * (R_a_gm2d + R_h_gm2d);									// g / m2 / day Total respiration from auto- and heterotrophs.
	R_total_kgd		= R_total_gm2d * (Length * Width)/ 1000; 				// kg /day total respiration per day per grid cell
	R_total_mgL		= R_total_kgd / water_Riv_total_in / 86400 * 1000;                      // mg / L
            
	NEP_gm2d 		= GPP_gm2d + R_total_gm2d;                                              // g / m2 / day
	NEP_kgd			= NEP_gm2d * (Length * Width)/ 1000;
	NEP_mgL			= NEP_kgd / water_Riv_total_in / 86400 * 1000;
 
      }
        }
        
 /*       
        if (itemID == 1) {
            printf("%d-%d-%d, GPP_gm2d = %f, GPP_kgd = %f, GPP_mgL = %f\n",MFDateGetCurrentYear(), MFDateGetCurrentMonth(), MFDateGetCurrentDay(), GPP_gm2d, GPP_kgd, GPP_mgL);
            printf("R_a_gm2d = %f, R_h_gm2d = %f, R_total_gm2d = %f, R_total_kgd = %f, R_total_mgL = %f\n",R_a_gm2d, R_h_gm2d, R_total_gm2d, R_total_kgd, R_total_mgL);
            printf("NEP_gm2d = %f, NEP_kgd = %f, NEP_mgL = %f\n", NEP_gm2d, NEP_kgd, NEP_mgL);
            printf("PAR_benthic = %f, order = %f, orderSwitch = %f, water_Riv_Total_in = %f\n", PAR_benthic, riverOrder, orderSwitch, water_Riv_total_in);

        }
   */    
        MFVarSetFloat (_MDOutNEPID, 	 itemID,   NEP_kgd);
        MFVarSetFloat (_MDOutGPPID, 	 itemID,   GPP_kgd);
        MFVarSetFloat (_MDOutRkgdID, 	 itemID,   R_total_kgd);
        MFVarSetFloat (_MDOutGPPgm2dID,	 itemID,   GPP_gm2d);
        MFVarSetFloat (_MDOutRgm2dID,	 itemID,   R_total_gm2d);
        MFVarSetFloat (_MDOutDryCID,	 itemID,   dryC_kgd);
        MFVarSetFloat (_MDOutDOCID,	 itemID,   DOC_kgd);
        MFVarSetFloat (_MDOutTotalCID,	 itemID,   totalC_kgd);

}


int MDNEPDef ()       {					//curly brackets indicate the start of the function (mddo.def)- everything within the curly brackets is part of the function.
MFDefEntering ("Net Ecosystem Production (NEP)");
if  (   ((_MDInPARBenthicID       = MFVarGetID (MDVarPARBenthic,             "MJ/m2/d",  MFInput, MFState, MFBoundary)) == CMfailed) ||
        ((_MDInDischargeID        = MFVarGetID (MDVarDischarge,                 "m3/s",  MFInput, MFState, MFBoundary))	== CMfailed) ||
     	((_MDInRiverWidthID	  = MFVarGetID (MDVarRiverWidth,                   "m",  MFInput, MFState, MFBoundary)) == CMfailed) ||
        ((_MDInSCALERBasinID      = MFVarGetID (MDVarSCALERBasinID,                "-",  MFInput, MFState, MFBoundary)) == CMfailed) ||
        ((_MDInRunoffVolumeID 	  = MFVarGetID (MDVarRunoffVolume,              "m3/s",  MFInput, MFState, MFBoundary)) == CMfailed) ||
        ((_MDInCanopyCover2ID 	  = MFVarGetID (MDVarCanopyCover2,                 "-",  MFInput, MFState, MFBoundary)) == CMfailed) ||
        ((_MDInOrderSwitchID 	  = MFVarGetID (MDVarOrderSwitch,                  "-",  MFInput, MFState, MFBoundary)) == CMfailed) ||
        ((_MDInRiverOrderID 	  = MFVarGetID (MDVarRiverOrder,                  "-",  MFInput, MFState, MFBoundary)) == CMfailed) ||
      	((_MDInWTemp_QxTID        = MFVarGetID (MDVarWTemp_QxT,                 "degC",  MFInput,  MFState, MFBoundary))   == CMfailed) ||	//RJS 013112
       	((_MDInTrefID             = MFVarGetID (MDVarTref,                       "degC",  MFInput,  MFState, MFBoundary))   == CMfailed) ||	//RJS 013112
        ((_MDOutNEPID   	  = MFVarGetID (MDVarNEP,                       "kg/d", MFOutput, MFFlux, MFBoundary)) == CMfailed) ||
     	((_MDOutGPPID   	  = MFVarGetID (MDVarGPP,                       "kg/d", MFOutput, MFFlux, MFBoundary)) == CMfailed) ||
     	((_MDOutGPPgm2dID   	  = MFVarGetID (MDVarGPPgm2d,                 "g/m2/d", MFOutput, MFState, MFBoundary)) == CMfailed) ||
     	((_MDOutRkgdID   	  = MFVarGetID (MDVarRkgd,                      "kg/d", MFOutput, MFFlux, MFBoundary)) == CMfailed) ||
     	((_MDOutRgm2dID   	  = MFVarGetID (MDVarRgm2d,                   "g/m2/d", MFOutput, MFState, MFBoundary)) == CMfailed) ||
     	((_MDOutDryCID   	  = MFVarGetID (MDVarDryC_kgd,                      "kg/d", MFOutput, MFFlux, MFBoundary)) == CMfailed) ||
       	((_MDOutDOCID   	  = MFVarGetID (MDVarDOC_kgd,                       "kg/d", MFOutput, MFFlux, MFBoundary)) == CMfailed) ||
       	((_MDOutTotalCID   	  = MFVarGetID (MDVarTotalC_kgd,                    "kg/d", MFOutput, MFFlux, MFBoundary)) == CMfailed) ||

        ((MFModelAddFunction (_MDNEP) == CMfailed))) return (CMfailed);
MFDefLeaving("Net Ecosystem Production (NEP)");
return (_MDOutNEPID);
}
