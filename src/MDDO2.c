/******************************************************************************
GHAAS Water Balance/Transport Model V3.0
Global Hydrologic Archive and Analysis System
Copyright 1994-2007, University of New Hampshire

MDDO2.c

Ken.r.sheehan@gmail.com

Module for dissolved oxygen (SCALER)

*******************************************************************************/

#include <stdio.h>
#include <string.h>
#include <cm.h>
#include <MF.h>
#include <MD.h>
#include <math.h>

// INPUT													Units
static int _MDInDischargeID					= MFUnset; // m3 / s
static int _MDInDischarge0ID					= MFUnset; // m3 / s
static int _MDInRunoffID					= MFUnset; // mm / day
static int _MDInRunoffVolumeID					= MFUnset; // m3 / s
static int _MDInBaseFlowID					= MFUnset; // mm / day

static int _MDInWTemp_QxTID					= MFUnset; // degrees celsius
static int _MDInRiverStorageID					= MFUnset; // m3 / s
static int _MDInRiverStorageChgID				= MFUnset; // m3 / s
static int _MDInRiverWidthID					= MFUnset; // m
static int _MDInRiverDepthID					= MFUnset; // m

static int _MDInRiverStorMassDO2ID				= MFUnset; // kg/time
static int _MDInRiverStorMassChgDO2ID                           = MFUnset; // kg/time
static int _MDInRiverMixingMassDO2ID                            = MFUnset; // kg/time
static int _MDInRiverStorMixingMassDO2ID                        = MFUnset; // kg/time
static int _MDInRiverStorMixingMassChgDO2ID                     = MFUnset; // kg/time
static int _MDInWTempSurfROID					= MFUnset; // degrees celsius
static int _MDInWTempGrdWaterID					= MFUnset; // degress celsius
static int _MDInPARBenthicID					= MFUnset; // MJ/m2/d converted to mol/m2/day 
static int _MDInNEPID                                           = MFUnset;  // NEP kg/d
static int _MDInGPPID                                           = MFUnset;  // GPP kg/d
static int _MDInRkgdID                                          = MFUnset;  // R kg/d
static int _MDInTauID                    = MFUnset;
static int _MDInPhiID                    = MFUnset;
static int _MDInOrderSwitchID            = MFUnset;
static int _MDInRiverOrderID            = MFUnset;
static int _MDInKFactorID                = MFUnset;
static int _MDInSlopeID                = MFUnset;
static int _MDInSWrunoffDOsatID         = MFUnset;
static int _MDInGWrunoffDOsatID         = MFUnset;
static int _MDInAerationApproachID	= MFUnset;


// METABOLISM INPUT

// OUTPUT
static int _MDOutRiverMassDO2ID					= MFUnset; // kg/day
static int _MDOutRiverStorMassDO2ID				= MFUnset; // kg/day
static int _MDOutRiverStorMassChgDO2ID                          = MFUnset; // kg/day
static int _MDOutRiverConcDO2ID					= MFUnset; // mg/L mean daily concentration per grid cell
static int _MDOutRiverMixingMassDO2ID                           = MFUnset; // kg/day
static int _MDOutRiverMixingStorMassDO2ID                       = MFUnset; // kg/day
static int _MDOutRiverMixingStorMassChgDO2ID                    = MFUnset; // kg/day
static int _MDOutRiverMixingConcDO2ID                           = MFUnset; // mg/L
static int _MDOutNetChangeDO2ID					= MFUnset; // mg/L
static int _MDOutKID						= MFUnset; // Aeration coefficient (calculated in each grid cell)
static int _MDOutKTempCorrectedID				= MFUnset; // Aeration coefficient corrected for temperature
static int _MDOutAerationID					= MFUnset; // Aeration in g / m2 / day
static int _MDOutBsaID						= MFUnset; // Basic Benthic surface area (Length * Width)
static int _MDOutNetDOID			    		= MFUnset; // Net DO output in g / m2 / day
static int _MDOutCsat2ID					= MFUnset; // DO at Sat for Calc of Aeration in the equation = K(t) * (Csat-Cmeas)
static int _MDOutAerationMassID					= MFUnset; // DO mass in kg/day of Aeration
static int _MDOutAerationConcID					= MFUnset; // mg / L : Concentration of DO in the grid cell from aeration
static int _MDOutTravelTimeID					= MFUnset; // Traveltime of water in each grid cell in
static int _MDOutDODiffID					= MFUnset; // Difference between DO concentration output from each cell and concentration at saturation.
static int _MDOutdDO2dtGm2dID					= MFUnset; // gm/m2/day of change in DO per grid cell
static int _MDOutdDO2dtKgdID					= MFUnset; // kg/day of change in DO per grid cell
static int _MDOutGWMassID					= MFUnset; // kg/day of DO in Groundwater per grid cell
static int _MDOutGWGm2dID					= MFUnset; // g/m2/day of DO in Groundwater per grid cell
static int _MDOutmassbalanceDO2ID                               = MFUnset;
static int _MDOutSWMassID                                       = MFUnset;
static int _MDOutMB2DO2ID                                       = MFUnset;
static int _MDOutAerationAdjustID                               = MFUnset;

static void _MDDO2 (int itemID) {
	// ESTABLISHING LOCAL VARIABLES
	int day                                 = 0;	// day of month
	int month                               = 0;	// month of year
	int year                                = 0;	// year
	//float Q_in                            = 0.0;	// incoming discharge (m3/s) includes runoff - - Note- Discharge from upstream plus (+) local runoff before routing. Muskingum method has not yet been applied.
	float Q_out                             = 0.0;	// outgoing discharge (m3/s) // cubic meters/second m3/sec - Discharge leaving this grid cell after routing.
	float Runoff                            = 0.0;	// local grid cell runoff (mm/d)
	float RunoffVolume                      = 0.0;	// local grid cell runoff (m3/s)
	float RiverTemp                         = 0.0;	// River temperature (deg C)
	float do2_SWrunoff_conc                 = 0.0;	// concentration of do2 in SWrunoff (mg/L)
	float do2_GWrunoff_conc                 = 0.0;	// concentration of do2 in GWrunoff (mg/L)
	float do2_SWrunoff_mass                 = 0.0;	// mass of do2 in SWrunoff (kg/d)
	float do2_GWrunoff_mass                 = 0.0; 	// mass of do2 in GWrunoff (kg/d)
	float do2_Riv_mass_in                   = 0.0;	// do2 mass in incoming river discharge (kg/d)
	float do2_Riv_mass_out                  = 0.0;	// o2 mass leaving in river grid cell (kg/d)
	float do2_RivStor_mass_pre              = 0.0; 	// mass of do2 remaining from yesterday in local river (kg/d)
	float do2_RivStor_mass_post             = 0.0;	// mass of do2 to remain in river until tomorrow kg/day
	float water_RivStor_pre		 	= 0.0; 	// volume of water remaining from yesterday in river grid cell (m3/d)
	float water_RivStor_post 	 	= 0.0;	// volume of water that remains in grid cell until tomorrow (m3/d)
	float water_RivStor_chg 	 	= 0.0; 	// change in volume of water stored in river over current time step (m3/day)
	float do2_RivStor_mass_chg		= 0.0; 	// (kg/day)
	float do2_RivStor_mass_chg2		= 0.0; 	// (kg/day)
	float do2_Riv_mass_total_pre 		= 0.0; 	// total mass of do2 in River before processing (kg/d)
	float water_Riv_total_in     		= 0.0;	// total water in river cubic meters/day (m3/sec)
	float do2_Riv_conc_total_pre 		= 0.0; 	// mg / L total concentration of do2 in River before processing
	float do2_Riv_conc_total_post		= 0.0; 	// mg / L total concentration of do2 in River after processing
	float do2_Riv_mixing_mass_in 		= 0.0; 	// (kg/day)
	float do2_RivStor_mixing_mass_post 	= 0.0; 	// (kg/day)
	float do2_RivStor_mixing_mass_chg 	= 0.0; 	// (kg/day)
	float do2_RivStor_mixing_mass_pre 	= 0.0; 	// (kg/day)
	float do2_RivStor_mixing_mass_chg2 	= 0.0; 	// (kg/day)
	float do2_Riv_mixing_mass_total 	= 0.0; 	// (kg/day)
	float do2_Riv_mixing_conc_total		= 0.0; 	// (mg/L)
	float do2_Riv_mass_total_post		= 0.0; 	// (kg/day)ned
	float do2_Riv_mixing_mass_out		= 0.0; 	// (kg/day)
	//float massbalanceDo2			= 0.0;

	// ENTERING LOCAL VARIABLES FOR METABOLISM- Based on general equestion of NEP = GPP - Respiration of Heterotrophs - Respiration of Autotrophs
	float GPP_kgd				= 0.0;  // kg / day GPP in the grid cell.
	float Width				= 0.0; 	// meters
	float Depth				= 0.0; 	// meters
	float Length				= 0.0; 	// meters
	float Slope				= 0.0; 	// Degrees calculated by rise in meters over run in meters
	float K					= 0.0; 	// a per unit time value which must match rest of equation, which is per daily time step. Though often reported in minute or  1/minute
	float Net_DO				= 0.0; 	// g / m2 / day - net change in DO per grid cell due to Reaeration, GPP and Respiration.
	float Aeration_gm2d			= 0.0; 	// g / m2 / day (+/-) oxygen as a function of K_temp corrected * (Cmeas - Csat)
	float Aeration_kgd			= 0.0; 	// kg / day Mass of oxygen in the grid cell added or removed as a result of Aeration.
	float Aeration_kgd_adjust		= 0.0; 	// kg / day Change in Mass of oxygen in the grid cell added or removed as a result of Aeration.
	float Aeration_mgL			= 0.0;  // mg / L : Concentration of oxgen in the grid cell, used during calculation of A, which contains Cmeas, all oxygen in the cell including GPP_gm2d, R_gm2d, GWRunoff_conc and SWRunoff_conc.
	float K_temp_corrected 			= 0.0; 	// K corrected for temperature based on Bott 2006 where K(t) = K(1.026^(temp - 20))
	float C_sat2				= 0.0; 	// Concentration at Saturation mg/L
		//float C_meas_mean_def_mass	= 0.0;  // Measure mean 24 hour deficit in mass of DO kg/day
	//float C_meas				= 0.0; 	//
	//float Total_DO			= 0.0; 	// Total DO in the grid cell (different concentrations)
	float C_meas_conc			= 0.0; 	// concentration of DO (mg/l)

	float Travel_Time 			= 0.0; 	// days
	float Velocity				= 0.0; 	// m / s :
	float Bsa				= 0.0; 	// Benthic surface area (meters squared)
	float DO_Diff				= 0.0; 	// proportional difference in DO vs DO at saturation (percent)
	float Baseflow_mmd			= 0.0; 	// mm  / day
	float GWrunoff_vol			= 0.0; 	// m3  / sec
	float SWrunoff_vol     			= 0.0;  // m3  / sec
	float SWrunoff_T			= 0.0;	// deg Celsius
	float GWrunoff_T			= 0.0;  // deg Celsius
	float NEP_kgd 				= 0.0;  // NEP in kg / day : Converted from g / m2 / day
        float do2_Riv_mass_total_wNEP           = 0.0;  // kg / day
        float do2_Riv_conc_total_wNEP           = 0.0;  // mg / L
	float dDO2dt_kgd			= 0.0;  // kg / day
	float dDO2dt_kgd_new			= 0.0;  // kg / day
	float dDO2dt_kgd_diff			= 0.0;  // kg / day
	float dDO2dt_gm2d			= 0.0;  // g / m2 / day
	float dDO2dt_mgL			= 0.0;  // mg / L
	float C_sat2_kgd			= 0.0;  // kg / day value for water at oxygen saturation
	//float RiverOrder			= 0.0;
	float Ecan				= 0.0;
	float do2_GWrunoff_gm2d			= 0.0; // Groundwater DO in g/m2/day
        float massbalanceDO2                    = 0.0;
        float R_kgd                             = 0.0;
        float foo                               = 0.0;
        float Aer0                              = 0.0;
        float MB2                               = 0.0;
        float tau                               = 0.0;
        float phi                               = 0.0;
        float orderSwitch                       = 0.0;
        float riverOrder                        = 0.0;
        float Kfactor                           = 0.0;
        float SWrunoff_DOsat                    = 0.0;
        float GWrunoff_DOsat                    = 0.0;
	float K600				= 0.0;
	float Aeration_Approach			= 0.0;
        float Width2                            = 0.0;
        int cell                                = 780;

	// READ-IN VALUES FROM THE GREATER MODEL
	day 			   		= MFDateGetCurrentDay();
	month			   		= MFDateGetCurrentMonth();
	year 			   		= MFDateGetCurrentYear();
//	Q_in			   		= MFVarGetFloat (_MDInDischarge0ID,      itemID, 0.0);		// NO routing in this model, so we removed this variable: this includes local runoff - Note- Discharge from upstream plus (+) local runoff before routing.																					//Muskingum method has not yet been applied.
	Q_out			   		= MFVarGetFloat (_MDInDischargeID,       itemID, 0.0);  	// cubic meters/second m3/sec - Discharge leaving this grid cell after routing.
	Runoff			   		= MFVarGetFloat (_MDInRunoffID,          itemID, 0.0);  	// mm/day
	RunoffVolume                            = MFVarGetFloat (_MDInRunoffVolumeID,    itemID, 0.0);  	// cubic meters per second m3/sec
	RiverTemp		   		= MFVarGetFloat (_MDInWTemp_QxTID,       itemID, 0.0); 		// degrees celsius 9.722
	water_RivStor_post                      = MFVarGetFloat (_MDInRiverStorageID,    itemID, 0.0);  	// kg/day
	water_RivStor_chg                       = MFVarGetFloat (_MDInRiverStorageChgID, itemID, 0.0);		// kg/day
//       Width                                   = MFVarGetFloat (_MDInRiverWidthID, itemID, 0.0);			// meters
        Depth                                   = MFVarGetFloat (_MDInRiverDepthID, itemID, 0.0);			// meters
        Width                                   = MFVarGetFloat (_MDInRiverWidthID, itemID, 0.0);			// meters
        Length			   		= MFModelGetLength(itemID);									// meters
        tau = MFVarGetFloat (_MDInTauID,   itemID, 0.0);
        phi = MFVarGetFloat (_MDInPhiID,   itemID, 0.0);
        
    orderSwitch                                 = MFVarGetFloat (_MDInOrderSwitchID, itemID, 0.0);              // order river in which calcs proceed    
    riverOrder                                  = MFVarGetFloat (_MDInRiverOrderID, itemID, 0.0);              // River Order    
    Baseflow_mmd                                = MFVarGetFloat (_MDInBaseFlowID, itemID, 0.0);				// Baseflow mm / day
    GWrunoff_T                                  = MFVarGetFloat (_MDInWTempGrdWaterID, itemID, 0.0);		// GWrunoff temperature (degC) for KNZ Runs, set this parameter to 9.722
    SWrunoff_T                                  = MFVarGetFloat (_MDInWTempSurfROID,   itemID, 0.0);		// SWrunoff temperature (degC) for KNZ Runs, set this parameter to 9.722
    Kfactor                                    = MFVarGetFloat (_MDInKFactorID,   itemID, 0.0);		//  K value correction term
    Slope                                      = MFVarGetFloat (_MDInSlopeID,   itemID, 0.0);		//  River Slope

	do2_Riv_mass_in       			= MFVarGetFloat (_MDOutRiverMassDO2ID,    itemID, 0.0); 			// mass of do2 coming from the upstream grid cell
//	do2_RivStor_mass_post 			= MFVarGetFloat (_MDOutRiverStorMassDO2ID, itemID, 0.0);			// is always 0, unless routing is activated
//	do2_RivStor_mass_chg  			= MFVarGetFloat (_MDOutRiverStorMassChgDO2ID, itemID, 0.0);			// is always 0, unless routing is activated
	do2_Riv_mixing_mass_in                  = MFVarGetFloat (_MDOutRiverMixingMassDO2ID, itemID, 0.0); 			// mass of do2 in mixing coming from the upstream grid cell
//	do2_RivStor_mixing_mass_post            = MFVarGetFloat (_MDOutRiverMixingStorMassDO2ID, itemID, 0.0);		// is always 0, unless routing is activated
//	do2_RivStor_mixing_mass_chg             = MFVarGetFloat (_MDOutRiverMixingStorMassChgDO2ID, itemID, 0.0);	// is always 0, unless routing activated
	NEP_kgd                                 = MFVarGetFloat (_MDInNEPID, itemID, 0.0); 			// NEP kgd
	GPP_kgd                                 = MFVarGetFloat (_MDInGPPID, itemID, 0.0); 			// NEP kgd
	R_kgd                                   = MFVarGetFloat (_MDInRkgdID, itemID, 0.0); 			// NEP kgd
        SWrunoff_DOsat                          = MFVarGetFloat (_MDInSWrunoffDOsatID, itemID, 0.0);
        GWrunoff_DOsat                          = MFVarGetFloat (_MDInGWrunoffDOsatID, itemID, 0.0);
        Aeration_Approach                       = MFVarGetFloat (_MDInAerationApproachID, itemID, 0.0);		// 1 = Ken, 2 = Raymond, 3 = Bowden

        
        Width2 = tau * pow (Q_out, phi);
	// LOADING CALCULATIONS
	GWrunoff_vol		   		= Runoff > 0.000001 ? (Baseflow_mmd / Runoff) * RunoffVolume : 0.0;		// m3 / s
	SWrunoff_vol		  	 	= Runoff > 0.000001 ? RunoffVolume - GWrunoff_vol : 0.0;					// m3 / s
	do2_SWrunoff_conc                       = (SWrunoff_DOsat / 100) /* % saturation */ * (14.62 - (0.3898 * SWrunoff_T) + (0.006969 * SWrunoff_T * SWrunoff_T) - (5.897 * 0.00001 * (SWrunoff_T * SWrunoff_T * SWrunoff_T)) * 0.99964); // mg / L Based on the Duke 1973 method.  Changed from 0.95 to 0.7 based on Wil's input 10/5/2016
	do2_GWrunoff_conc                       = (GWrunoff_DOsat / 100) /* % saturation */ * (14.62 - (0.3898 * GWrunoff_T) + (0.006969 * GWrunoff_T * GWrunoff_T) - (5.897 * 0.00001 * (GWrunoff_T * GWrunoff_T * GWrunoff_T)) * 0.99964); // mg / L Currently set at 30% Saturation.
	do2_SWrunoff_mass                       = ((do2_SWrunoff_conc * SWrunoff_vol)/ 1000) * 86400; 				// kg/d : converted from conc
	do2_GWrunoff_mass                       = ((do2_GWrunoff_conc * GWrunoff_vol)/ 1000) * 86400;				// kg/d : converted from conc


	water_RivStor_pre                       = water_RivStor_post - water_RivStor_chg;				// m3 / day  (0 unless routing activated)
	do2_RivStor_mass_pre                    = do2_RivStor_mass_post - do2_RivStor_mass_chg;				// kg/day (0 unless routing activated)

	do2_Riv_mass_total_pre                  = do2_Riv_mass_in + do2_SWrunoff_mass + do2_GWrunoff_mass ;  	// kg/day

	water_Riv_total_in                      = Q_out; // m3 / s
	do2_Riv_conc_total_pre                  = water_Riv_total_in > 0.00000001 ? (do2_Riv_mass_total_pre * 1000) / (water_Riv_total_in * 86400) : 0.0; //
        
   
	// Metabolism Calculations
	Velocity		= (Width > 0.01) && (Depth > 0.01) ? water_Riv_total_in / (Width * Depth) : 0.01; 								//m / s
//	Hydraulic_Radius	= Depth > 0.01 ? (Width * Depth) / ((Depth * Depth) + Width): 0.01;
//	n 			= 0.052; // n may be varied. Larger n slows mannings velocity, and smaller n speeds velocity up. n is also known as the roughness coefficient.
//	Velocity_Mannings	= (1 / n) * ((pow (Hydraulic_Radius, 0.666667) * pow (((Slope * 0.25 )/ 100), 0.5))); 		// Slope in the input data is a layer with values of percent (50 for example). Slope is divided by 100 to get it in the correct format for the equation
//	Depth2			= (Width > 0.01) && (Depth > 0.01) ? sqrt(water_Riv_total_in / ((Width/Depth) * Velocity_Mannings)) : 0.01; 	// m : depth corrected for mannings velocity
//	Width2			= Depth > 0.01 ? (Width / Depth) * Depth2 : 0.01; 										// m :  width corrected for mannings velocity
//	Hydraulic_Radius_Corr   = Depth > 0.01 ? (Width2 * Depth2) / ((Depth2 * Depth2) + Width2): 0.01;																								   	   	//example: n = (0.0028 + 0.005 + 0.010 + 0.02 + 0.025)1.15 =  from "methods in stream ecology" pg. 66;   m where m is degree of meandering (minor: 1.0, appreciable: 1.15, severe: 1.3)
//      Velocity_Mannings_Corr      = (1 / n) * ((pow (Hydraulic_Radius_Corr, 0.666667) * pow (((Slope * 0.25 )/ 100), 0.5)));

	do2_GWrunoff_gm2d	= do2_GWrunoff_mass * 1000 / (Length * Width);						// g/m2/day DO in GW converted from GW kg/day

	Travel_Time		= Velocity > 0.0001 ? (Length / Velocity) / 86400 : 0.0;					// Time it takes in days for a molecule of water to travel one grid cell. Value is a proportion of a day.
         
	do2_Riv_mass_total_wNEP = do2_Riv_mass_total_pre + NEP_kgd;
	do2_Riv_mass_total_wNEP = do2_Riv_mass_total_pre + NEP_kgd < 0? 0 : do2_Riv_mass_total_pre + NEP_kgd;

	do2_Riv_conc_total_wNEP	= do2_Riv_mass_total_wNEP / water_Riv_total_in / 86400 * 1000 ; // mg / L

//        if (itemID == cell) {
//            printf("****** do2_Riv_mass_total_pre = %f, conc_pre = %f, NEP_kgd = %f, GPP = %f, R = %f, do2_Riv_mass_total_wNEP = %f, conc_total_wNEP = %f\n", do2_Riv_mass_total_pre, do2_Riv_conc_total_pre, NEP_kgd, GPP_kgd, R_kgd, do2_Riv_mass_total_wNEP, do2_Riv_conc_total_wNEP);
//        }
      
        
	if (Aeration_Approach == 1) {
		K = Kfactor * (543 * (pow (((Slope * 1.0) /100), 0.6236) * (pow (Velocity,0.5325)) * pow (Depth,-0.7258))); 	// RJS 10/2/2016 added Kfactor per Ken's request						//((1740 * (pow (Velocity_Mannings, 0.46))) * (pow ((Slope/100), 0.79)) * (pow (Depth, 0.74)))/1440; // Jha et al. 2001 equation (5.79 * (pow (Velocity_Mannings, 0.5) * pow (Depth, -0.25)))/1440;//check this, but equation is for a day, so divided by 1440 to get it in a version per minutes to match rest of equation
	}
	
	if (Aeration_Approach == 2) {
		K600 = 1162 * pow((Slope / 100), 0.77) * pow(Velocity, 0.85);	// Raymond et al.
		K    = K600 / Depth;
	}
	
	if (Aeration_Approach == 3) {
		K    = 28300 * (Slope / 100) * Velocity;	// Bowden
	}
	
	K_temp_corrected   	= K * pow (1.025,(RiverTemp - 20));// K corrected to 20 C

	C_sat2			= 14.62 - (0.3898 * RiverTemp) + (0.006969 * RiverTemp * RiverTemp) - (5.897 * 0.00001 * (RiverTemp * RiverTemp * RiverTemp)) * 0.99964; //
   	C_sat2_kgd		= C_sat2 * water_Riv_total_in / 1000 * 86400; 																							// kg / day at saturation
	C_meas_conc 		= do2_Riv_conc_total_wNEP;

        if (riverOrder >= orderSwitch) {
	if (water_Riv_total_in >= 0.000001) {

        Bsa			= Length * Width;// m2 : depending upon which velocity is wanted, either use width or width2 if using mannings velocity.

        
   	Aeration_mgL		= (K_temp_corrected * (C_sat2 - C_meas_conc)) * Travel_Time; 					//(mg / L )Note: had C_sat2_kgd --> changed to C_sat2, which is mg/L ; Multiplied by Travel time (proportion of a day to get through the grid cell) to get the daily value.
  	Aeration_kgd		= Aeration_mgL * water_Riv_total_in / 1000 * 86400; 										// kg / day Aeration
 	Aeration_gm2d		= Aeration_kgd * 1000 / (Length * Width); 	

//        if (itemID == cell) {
//            printf("Initial: A_mgL = %f, A_kgd = %f, A_gm2d = %f, C_sat2 = %f, C_sat2_kgd = %f, C_meas_conc = %f\n", Aeration_mgL, Aeration_kgd, Aeration_gm2d, C_sat2, C_sat2_kgd, C_meas_conc);    
//        }
        
        
/*        if ((itemID == 49) || (itemID == 48))  {
            
         printf("::::: Initial Aeration :::::  Aeration_kgd = %f, Aeration_gm2d = %f, Aeration_mgL = %f, Aeration_kgd_adjust = %f, Aer0 = %f\n", Aeration_kgd, Aeration_gm2d, Aeration_mgL, Aeration_kgd_adjust, Aer0);
         
        }
*/        
//    dDO2dt_kgd			= dDO2dt_kgd > C_sat2_kgd ? C_sat2_kgd: GPP_kgd + R_total_kgd + Aeration_kgd;
//    dDO2dt_kgd_new			= dDO2dt_kgd < C_sat2_kgd ? dDO2dt_kgd: C_sat2_kgd;
//    dDO2dt_kgd_diff						= dDO2dt_kgd_new - dDO2dt_kgd;								// new 12/1/2015
////    dDO2dt_kgd			= dDO2dt_kgd < C_sat2_kgd ? dDO2dt_kgd: C_sat2_kgd;	// old 12/01/15

//    Aeration_kgd 					= Aeration_kgd + dDO2dt_kgd_diff;								// new 12/1/2015
//    Aeration_gm2d					= Aeration_kgd * 1000 / (Length * Width); 	 					// new 12/1/2015
//    Aeration_mgL					= Aeration_kgd / water_Riv_total_in * 1000 / 86400;	 			// new 12/1/2015

//    dDO2dt_gm2d			= dDO2dt_kgd_new * 1000 / Length * Width; 																								// GPP_gm2d - R_total_gm2d + Aeration_gm2d;
//    Aeration_mgL		= Aeration_mgL > C_sat2 ? C_sat2 * Travel_Time :  (K_temp_corrected * (C_sat2 - C_meas_conc)) * Travel_Time;
//    Aeration_kgd		= Aeration_mgL * water_Riv_total_in / 1000 * 86400; 										// kg / day Aeration
//    Aeration_gm2d		= Aeration_kgd * 1000 / (Length * Width); 													// g / m2 / day : This is the same unit in which GPP_gm2d and R_gm2d are initially entered as parameters in the model.
//    dDO2dt_mgL		= dDO2dt_kgd_new / water_Riv_total_in / 86400 * 1000; 																					//GPP_mgL - R_total_mgL + Aeration_mgL;

//    do2_Riv_mass_total_post 	= do2_Riv_mass_total_pre + dDO2dt_kgd_new > C_sat2_kgd? C_sat2_kgd: do2_Riv_mass_total_pre + dDO2dt_kgd_new;								//
//    do2_Riv_mass_total_post 	= do2_Riv_mass_total_post < C_sat2_kgd? do2_Riv_mass_total_post: C_sat2_kgd;
////  do2_Riv_mass_total_post   = do2_Riv_mass_total_post < 0 ? 0 : do2_Riv_mass_total_post;

/// START NEW ///
        foo                             = do2_Riv_mass_total_pre + Aeration_kgd + GPP_kgd + R_kgd;
        Aer0                            = Aeration_kgd;
    do2_Riv_mass_total_post 		= do2_Riv_mass_total_pre + Aeration_kgd + GPP_kgd + R_kgd > C_sat2_kgd ? C_sat2_kgd : do2_Riv_mass_total_pre + Aeration_kgd + GPP_kgd + R_kgd;	// new 12/2/2015
    do2_Riv_mass_total_post       	= do2_Riv_mass_total_post <= 0 ? 0.00001 : do2_Riv_mass_total_post;					// new 12/2/2015
    Aeration_kgd_adjust			= do2_Riv_mass_total_post - (do2_Riv_mass_total_pre + Aeration_kgd + GPP_kgd + R_kgd);		// new 12/2/2015														// new 12/2/2015
    Aeration_kgd			= Aeration_kgd + Aeration_kgd_adjust;									// new 12/2/2015
    Aeration_gm2d			= Aeration_kgd * 1000 / (Length * Width); 	 							// new 12/2/2015
    Aeration_mgL			= Aeration_kgd / water_Riv_total_in * 1000 / 86400;	 						// new 12/2/2015

//    if (itemID == cell) {
//        printf("foo = %f, do2_Riv_mass_total_post = %f, A_kgd_adjust = %f, A_kgd = %f, A_mgL = %f, A_gm2d = %f\n", foo, do2_Riv_mass_total_post, Aeration_kgd_adjust, Aeration_kgd, Aeration_mgL, Aeration_gm2d);
//    }
    
/// END NEW ///

    do2_Riv_conc_total_post 	= do2_Riv_mass_total_post / water_Riv_total_in / 86400 * 1000 < 0? 0: do2_Riv_mass_total_post / water_Riv_total_in / 86400 * 1000;

    DO_Diff			= do2_Riv_conc_total_post / C_sat2;

	}
        
	else {

	do2_Riv_conc_total_post = do2_Riv_conc_total_pre;
	do2_Riv_mass_total_post = do2_Riv_mass_total_pre;

	}
        }
        
       	else {

	do2_Riv_conc_total_post = do2_Riv_conc_total_pre;
	do2_Riv_mass_total_post = do2_Riv_mass_total_pre;

	} 
	do2_Riv_mass_out	= ((do2_Riv_conc_total_post * water_Riv_total_in) / 1000) * 86400; // kg/day attempt at fixing the above equation, which appears to be incorrect KS October 16, 2013

	do2_RivStor_mass_post  	= ((do2_Riv_conc_total_post * water_RivStor_post) / 1000) * 86400 ; //kg/day - // (routing, this will be 0 until routing is changed) kg/day attempt at fixing the above equation, which appears to be incorrect KS October 16, 2013
	do2_RivStor_mass_chg2  	= do2_RivStor_mass_post - do2_RivStor_mass_pre; //kg/day	(this will be 0 until routing is changed)

	//mixing calcs
	do2_RivStor_mixing_mass_pre 	= do2_RivStor_mixing_mass_post - do2_RivStor_mixing_mass_chg;	// (this will be 0 until routing is changed)
	do2_Riv_mixing_mass_total   	= do2_Riv_mixing_mass_in + do2_SWrunoff_mass + do2_GWrunoff_mass; 				// kg/d
	do2_Riv_mixing_conc_total   	= water_Riv_total_in > 0.0 ? (do2_Riv_mixing_mass_total / (water_Riv_total_in * 86400)) * 1000000 / 1000 : 0.0;	// mg/L
	do2_Riv_mixing_mass_out      	= do2_Riv_mixing_conc_total * Q_out / 1000000 * 1000 * 86400;
	do2_RivStor_mixing_mass_post 	= do2_Riv_mixing_conc_total * water_RivStor_post / 1000000 * 1000;		// (this will be 0 until routing is changed)
	do2_RivStor_mixing_mass_chg2 	= do2_RivStor_mixing_mass_post - do2_RivStor_mixing_mass_pre;		// (this will be 0 until routing is changed)
	// mass balance of d02
	massbalanceDO2			= (do2_Riv_mass_in + do2_SWrunoff_mass + do2_GWrunoff_mass) - do2_Riv_mass_out + NEP_kgd + Aeration_kgd; // kg/day
        MB2                             = (do2_Riv_mass_in + do2_SWrunoff_mass + do2_GWrunoff_mass) > 0.0 ? massbalanceDO2 / (do2_Riv_mass_in + do2_SWrunoff_mass + do2_GWrunoff_mass) : 0.0;
        
//	if (((itemID == 3760) || (itemID == 2107)) && (fabs(MB2) > 0.01)) {
/* 	if ((itemID == 1))  {

        	printf("***** %d-%d-%d, itemID = %d, MB = %f, MBprop= %f\n", MFDateGetCurrentYear(), MFDateGetCurrentMonth(), MFDateGetCurrentDay(), itemID, massbalanceDO2, massbalanceDO2 / ((do2_Riv_mass_in + do2_SWrunoff_mass + do2_GWrunoff_mass)));
                printf("********* %d-%d-%d, itemID = %d\n", MFDateGetCurrentYear(), MFDateGetCurrentMonth(), MFDateGetCurrentDay(), itemID);
		printf("Bsa = %f, Depth = %f, Width = %f, Length = %f, Slope = %f\n", Bsa, Depth, Width, Length, Slope);
		printf("RunoffVolume = %f, waterIn = %f, GWrunoff_vol = %f, SWrunoff_vol = %f, Runoff = %f, Baseflow_mmd = %f, \n", RunoffVolume, water_Riv_total_in, GWrunoff_vol, SWrunoff_vol, Runoff, Baseflow_mmd);
		printf("C_sat2 = %f, K = lines(BP$Distance, BP$PARBenthic, type="l",col="brown",lwd = 3)%f, KtempCorrected = %f, Kfactor = %f, C_sat2_kgd = %f, C_meas_conc = %f\n", C_sat2, K, K_temp_corrected, Kfactor, C_sat2_kgd, C_meas_conc);
		printf("do2_Riv_conc_total_post = %f\n", do2_Riv_conc_total_post);
		printf("DO_Diff = %f, RiverTemp = %f\n", DO_Diff, RiverTemp);
		printf("SWrunoff_conc = %f, GWrunoff_conc = %f, SWrunoff_mass = %f, GWrunoff_mass = %f\n",do2_SWrunoff_conc,do2_GWrunoff_conc,do2_SWrunoff_mass,do2_GWrunoff_mass);       
		printf("Riv_mass_in = %f, Riv_mass_out = %f, Riv_mass_pre = %f, Riv_mass_post = %f, foo = %f, \n", do2_Riv_mass_in, do2_Riv_mass_out,do2_Riv_mass_total_pre, do2_Riv_mass_total_post, foo);
		printf("Travel_Time = %f, velocity = %f\n", Travel_Time, Velocity);
		printf("Velocity = %f, Width = %f, Depth = %f\n", Velocity, Width, Depth);
		printf("NEP_kgd = %f, GPP_kgd = %f, R_kgd = %f\n", NEP_kgd, GPP_kgd, R_kgd);
		printf("Riv_mass_total_wNEP = %f, Riv_conc_total_wNEP = %f\n", do2_Riv_mass_total_wNEP, do2_Riv_conc_total_wNEP);
		printf("Aeration_kgd = %f, Aeration_gm2d = %f, Aeration_mgL = %f, Aeration_kgd_adjust = %f, Aer0 = %f\n", Aeration_kgd, Aeration_gm2d, Aeration_mgL, Aeration_kgd_adjust, Aer0);

	}
*/ 
 //       printf("w1=%f,w2=%f       ",Width,Width2);
        
	MFVarSetFloat (_MDOutRiverMassDO2ID, 		  itemID,               do2_Riv_mass_out);
	MFVarSetFloat (_MDOutRiverStorMassDO2ID,          itemID,          do2_RivStor_mass_post);
	//MFVarSetFloat (_MDOutRiverStorMassChgDO2ID,      itemID,         do2_RivStor_mass_chg2);
	MFVarSetFloat (_MDOutRiverConcDO2ID,              itemID,        do2_Riv_conc_total_post);
	MFVarSetFloat (_MDOutRiverMixingMassDO2ID,        itemID,        do2_Riv_mixing_mass_out);
	//MFVarSetFloat (_MDOutRiverMixingStorMassDO2ID,   itemID,  do2_RivStor_mixing_mass_post);
	//MFVarSetFloat (_MDOutRiverMixingStorMassChgDO2ID, itemID, do2_RivStor_mixing_mass_chg2);
	MFVarSetFloat (_MDOutRiverMixingConcDO2ID,        itemID,      do2_Riv_mixing_conc_total);
        MFVarSetFloat (_MDOutKID,                     	  itemID,                              K);
	MFVarSetFloat (_MDOutKTempCorrectedID,            itemID,               K_temp_corrected);
	MFVarSetFloat (_MDOutAerationID,                  itemID,                  Aeration_gm2d);
	MFVarSetFloat (_MDOutAerationConcID,              itemID,                   Aeration_mgL);
	MFVarSetFloat (_MDOutBsaID,                       itemID,                            Bsa);
	MFVarSetFloat (_MDOutNetDOID,                     itemID,                         Net_DO);
	MFVarSetFloat (_MDOutCsat2ID, 			  itemID,                         C_sat2);
	MFVarSetFloat (_MDOutAerationMassID,              itemID,                   Aeration_kgd);
	MFVarSetFloat (_MDOutTravelTimeID,             	  itemID,                    Travel_Time);
	MFVarSetFloat (_MDOutDODiffID, 			  itemID,                        DO_Diff);
	MFVarSetFloat (_MDOutdDO2dtKgdID, 		  itemID,                     dDO2dt_kgd_new);
	MFVarSetFloat (_MDOutdDO2dtGm2dID, 		  itemID, 		     dDO2dt_gm2d);
	MFVarSetFloat (_MDOutGWMassID, 			  itemID,              do2_GWrunoff_mass);
	MFVarSetFloat (_MDOutGWGm2dID, 			  itemID,   	       do2_GWrunoff_gm2d);
	MFVarSetFloat (_MDOutSWMassID, 			  itemID,              do2_SWrunoff_mass);
	MFVarSetFloat (_MDOutmassbalanceDO2ID, 		  itemID, 		  massbalanceDO2);
	MFVarSetFloat (_MDOutMB2DO2ID,                    itemID, 		  MB2);
        MFVarSetFloat (_MDOutAerationAdjustID,            itemID,                 Aeration_kgd_adjust);

        
//float balance = 																				//KRS 091912 defining the balance equation for DO- ask rob if we need a balance equation for this, because it's not directly a water function.
}
int MDDO2Def ()       {					//curly brackets indicate the start of the function (mddo.def)- everything within the curly brackets is part of the function.
MFDefEntering ("Dissolved Oxygen2");
if  (  // ((_MDInPARBenthicID                    = MDBgcRiverLightDef())                                                                                          == CMfailed)    ||
        ((_MDInNEPID                           = MDNEPDef())                                                                                          == CMfailed)    ||
	((_MDInPARBenthicID                    = MFVarGetID (MDVarPARBenthic,             "MJ/m2/d",  MFInput,  MFState, MFBoundary))   == CMfailed) ||	//RJS 013112
	((_MDInWTemp_QxTID                     = MFVarGetID (MDVarWTemp_QxT,                 "degC",  MFInput,  MFState, MFBoundary))   == CMfailed) ||	//RJS 013112
        ((_MDInDischarge0ID		       = MFVarGetID (MDVarDischarge0,                "m3/s",  MFInput,  MFState, MFBoundary)) 	== CMfailed) 	||
	((_MDInDischargeID                     = MFVarGetID (MDVarDischarge,                 "m3/s",  MFInput,  MFState, MFBoundary)) 	== CMfailed) 	||
	((_MDInRunoffID 		       = MFVarGetID (MDVarRunoff,                      "mm",  MFInput,   MFFlux, MFBoundary)) 	== CMfailed) 	||
	((_MDInRiverWidthID		       = MFVarGetID (MDVarRiverWidth,                   "m",  MFInput,  MFState, MFBoundary)) 	== CMfailed) 	||
	((_MDInRiverDepthID		       = MFVarGetID (MDVarRiverDepth,                   "m",  MFInput,  MFState, MFBoundary)) 	== CMfailed) 	||
	((_MDInRunoffVolumeID 		       = MFVarGetID (MDVarRunoffVolume,              "m3/s",   MFInput,  MFState, MFBoundary)) 	== CMfailed) 	||
	((_MDInRiverStorageID 		       = MFVarGetID (MDVarRiverStorage,              "m3/s",   MFInput,  MFState,  MFInitial)) 	== CMfailed) 	|| // initial means there needs to be an initial value in the cell to starte the model.will set it to zero to start the model.
	((_MDInRiverStorageChgID	       = MFVarGetID (MDVarRiverStorageChg,           "m3/d",   MFInput,  MFState,  MFInitial))  == CMfailed)    ||
        ((_MDInWTempSurfROID                   = MFVarGetID (MDVarWTempSurfRunoff,           "degC",   MFInput,  MFState, MFBoundary))  == CMfailed) ||
	((_MDInWTempGrdWaterID                 = MFVarGetID (MDVarWTempGrdWater,             "degC",   MFInput,  MFState,  MFBoundary))  == CMfailed) ||
	((_MDInSWrunoffDOsatID                 = MFVarGetID (MDVarSWrunoffDOsat,             "%",   MFInput,  MFState,  MFBoundary))  == CMfailed) ||
	((_MDInGWrunoffDOsatID                 = MFVarGetID (MDVarGWrunoffDOsat,             "%",   MFInput,  MFState,  MFBoundary))  == CMfailed) ||
	((_MDInBaseFlowID                      = MFVarGetID (MDVarBaseFlow,                    "mm",   MFInput,   MFFlux, MFBoundary))  == CMfailed) ||
     	((_MDInGPPID                           = MFVarGetID (MDVarGPP,                       "kg/d", MFInput, MFState, MFBoundary)) == CMfailed) ||
     	((_MDInRkgdID                          = MFVarGetID (MDVarRkgd,                      "kg/d", MFInput, MFState, MFBoundary)) == CMfailed) ||
        ((_MDInTauID                           = MFVarGetID (MDVarTau,                   "-",      MFInput,  MFState, MFBoundary)) == CMfailed) ||
        ((_MDInPhiID                           = MFVarGetID (MDVarPhi,                   "-",      MFInput,  MFState, MFBoundary)) == CMfailed) ||                          
        ((_MDInRiverOrderID                    = MFVarGetID (MDVarRiverOrder,                  "-",  MFInput, MFState, MFBoundary)) == CMfailed) ||
        ((_MDInOrderSwitchID                   = MFVarGetID (MDVarOrderSwitch,                  "-",  MFInput, MFState, MFBoundary)) == CMfailed) ||
        ((_MDInAerationApproachID              = MFVarGetID (MDVarAerationApproach,             "-",  MFInput, MFState, MFBoundary)) == CMfailed) ||
        ((_MDInKFactorID                       = MFVarGetID (MDVarKFactor,                      "-",  MFInput, MFState, MFBoundary)) == CMfailed) ||
        ((_MDInSlopeID                         = MFVarGetID (MDVarRiverbedSlope2,                      "-",  MFInput, MFState, MFBoundary)) == CMfailed) ||
        ((_MDOutRiverMassDO2ID 		       = MFVarGetID (MDVarRiverMassDO2,              "kg/d",   MFRoute,   MFFlux, MFBoundary)) 	== CMfailed) 	||
	((_MDOutRiverConcDO2ID 		       = MFVarGetID (MDVarRiverConcDO2,              "mg/l",  MFOutput,  MFState, MFBoundary)) 	== CMfailed) 	||
	((_MDOutRiverStorMassDO2ID             = MFVarGetID (MDVarRiverStorMassDO2,            "kg",  MFOutput,  MFState,  MFInitial)) 	== CMfailed) 	||
	((_MDOutRiverStorMassChgDO2ID          = MFVarGetID (MDVarRiverStorMassChgDO2,       "kg/d",  MFOutput,   MFFlux, MFBoundary)) 	== CMfailed) 	||
	((_MDOutRiverMixingMassDO2ID           = MFVarGetID (MDVarRiverMixingMassDO2,        "kg/d",   MFRoute,   MFFlux, MFBoundary)) 	== CMfailed) 	||
	((_MDOutRiverMixingConcDO2ID           = MFVarGetID (MDVarRiverMixingConcDO2,        "mg/l",  MFOutput,  MFState, MFBoundary)) 	== CMfailed) ||
	((_MDOutRiverMixingStorMassDO2ID       = MFVarGetID (MDVarRiverStorMixingMassDO2,      "kg",  MFOutput,  MFState,  MFInitial)) 	== CMfailed) 	||
	((_MDOutRiverMixingStorMassChgDO2ID    = MFVarGetID (MDVarRiverMixingStorMassChgDO2, "kg/d",  MFOutput,   MFFlux, MFBoundary)) 	== CMfailed) ||
	((_MDOutNetChangeDO2ID                 = MFVarGetID (MDVarNetChangeDO2,              "mg/l",  MFOutput,  MFState, MFBoundary)) 	== CMfailed) ||
	((_MDOutNetDOID   		       = MFVarGetID (MDVarNetDO,                     "mg/l",  MFOutput,  MFState, MFBoundary)) 	== CMfailed) ||
	((_MDOutKID                            = MFVarGetID (MDVarK,                           "-" ,  MFOutput,  MFState, MFBoundary))  == CMfailed) ||
	((_MDOutAerationID                     = MFVarGetID (MDVarAeration,             "g/m2/day" ,  MFOutput,  MFState, MFBoundary))  == CMfailed) ||
	((_MDOutAerationAdjustID               = MFVarGetID (MDVarAerationAdjust,          "kg/day",  MFOutput,  MFFlux, MFBoundary))  == CMfailed) ||
	((_MDOutAerationMassID                 = MFVarGetID (MDVarAerationMass,            "kg/day",  MFOutput,  MFFlux, MFBoundary))  == CMfailed) ||
	((_MDOutAerationConcID                 = MFVarGetID (MDVarAerationConc,              "mg/l",  MFOutput,  MFState, MFBoundary))  == CMfailed) ||
	((_MDOutBsaID                          = MFVarGetID (MDVarBsa,                         "m2",  MFOutput,  MFState, MFBoundary))  == CMfailed) ||
	((_MDOutKTempCorrectedID               = MFVarGetID (MDVarKTempCorrected,               "-",  MFOutput,  MFState, MFBoundary))  == CMfailed) ||
	((_MDOutCsat2ID                        = MFVarGetID (MDVarCsat2,                     "mg/l",  MFOutput,  MFState, MFBoundary))  == CMfailed) ||
	((_MDOutTravelTimeID                   = MFVarGetID (MDVarTravelTime,                   "s",  MFOutput,  MFState, MFBoundary))  == CMfailed) ||
	((_MDOutDODiffID               	       = MFVarGetID (MDVarDODiff,                       "%",  MFOutput,  MFState, MFBoundary))  == CMfailed) ||
	((_MDOutdDO2dtKgdID                    = MFVarGetID (MDVardDO2dtKgd,               "kg/day",  MFOutput,  MFFlux, MFBoundary))  == CMfailed) || 
	((_MDOutdDO2dtGm2dID                   = MFVarGetID (MDVardDO2dtGm2d,            "g/m2/day",  MFOutput,  MFState, MFBoundary))  == CMfailed) ||
	((_MDOutSWMassID                       = MFVarGetID (MDVarSWMass,                  "kg/day",  MFOutput,  MFFlux, MFBoundary))  == CMfailed) ||
	((_MDOutGWMassID                       = MFVarGetID (MDVarGWMass,                  "kg/day",  MFOutput,  MFFlux, MFBoundary))  == CMfailed) ||
	((_MDOutGWGm2dID                       = MFVarGetID (MDVarGWGm2d,                "g/m2/day",  MFOutput,  MFState, MFBoundary))  == CMfailed) ||
	((_MDOutmassbalanceDO2ID	       = MFVarGetID (MDVarmassbalanceDO2 ,         "kg/day" , MFOutput, MFFlux, MFBoundary))  == CMfailed) ||
	((_MDOutMB2DO2ID                       = MFVarGetID (MDVarMB2DO2 ,               "-" , MFOutput, MFState, MFBoundary))  == CMfailed) ||

	((MFModelAddFunction (_MDDO2) == CMfailed))) return (CMfailed);
MFDefLeaving("Dissolved Oxygen2");
return (_MDOutRiverConcDO2ID);
}












