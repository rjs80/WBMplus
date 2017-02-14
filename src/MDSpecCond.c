/******************************************************************************

GHAAS Water Balance/Transport Model V3.0
Global Hydrologic Archive and Analysis System
Copyright 1994-2011, UNH - CCNY/CUNY

MDSpecCond.c  - Input and Routing of Specific Conductance

shan.zuidema@unh.edu  
 * 
 * Interprets salinazation of freshwater networks.
 * 

*******************************************************************************/
#include <stdio.h>
#include <cm.h>
#include <MF.h>
#include <MD.h>
#include <math.h>

// input
static int _MDInPrecipitationID        = MFUnset;
static int _MDInBaseFlowID             = MFUnset;
static int _MDInRunoffPoolReleaseID    = MFUnset;
static int _MDInSurfaceRunoffPoolID    = MFUnset;
static int _MDInSurfaceRunoffID        = MFUnset;
static int _MDInStormRunoffTotalID     = MFUnset;
static int _MDInDischargeID            = MFUnset;
static int _MDInRiverStorageID         = MFUnset;
static int _MDInRiverStorageChgID      = MFUnset;
static int _MDInRiverbedVelocityMeanID = MFUnset;
static int _MDInSinuosityID	       = MFUnset;
static int _MDInAirTemperatureID       = MFUnset;
static int _MDInAvailableWaterCapacityID = MFUnset;
static int _MDInSnowPackID              = MFUnset;
static int _MDInSPackChgID              = MFUnset;
static int _MDInSoilMoistID             = MFUnset;
static int _MDInGrdWatRechargeID        = MFUnset;
static int _MDInSoilPercolationID       = MFUnset;
static int _MDInGrdWatID                = MFUnset;
// Upstream Merge input
static int _MDInDINFluxID               = MFUnset;
static int _MDInWTempRiverID            = MFUnset;
static int _MDInRiverWidthID            = MFUnset;
static int _MDInLitterFall_POCID        = MFUnset;
static int _MDInLocalLoad_DOCID         = MFUnset;

//Parameters
static float _MDparPassiveSoilStoreFactor = 0.1  ; // Factor.  Multiplied by AWC to define passive storage in soil pool 
static float _MDSnowFallThreshold      = -0.29; // deg C where snowfall occurs
static float _MDparDeicerCl           = 3000.0; // Chloride in de-icing treatment (mg/L) 
static float _MDInDeicerClID           = MFUnset; // Same (Cl in de-icing) but in case of spatial input
static float _MDparPrecipCl           = 0.108; // Median Cl in precip (NADP - NH02) // 3.0; // Cl in clean precipitation (mg/L) (Back of the envelope from NADP wet dep. 2 mg/L)
static float _MDInHCIAID               = MFUnset; // Fraction of impervious areas that are directly connected to streams
static float _MDInWinterHCIAID         = MFUnset; // Fraction of impervious areas that are directly connected to streams (during snowfall)
static float _MDparGeoChem_ClWx_Rate  = 2.74e-7; // Ground Water Weathering Rate (kg / m3 / day) (BGC1997)
static float _MDparGreaterGroundWaterMixing = 100.0; // Added depth of the groundwater pool controlling mixing (about a years runoff) (mm)
static float _MDparGreaterGroundWaterExchange = 0.015; // Exchange rate of deep groundwater pool with baseflow generating groundwater (d-1)
static float _MDparDevelop_ClWx_Rate   = 2.74e-2; // Development scaled Cl loading rate (kg / m2 / day / developed area fraction)
static float _MDparAg_ClWx_Rate        = 1.36619e-4; // Agriculture scaled Cl (Cl) loading rate (kg Cl / m2 /day / agriculture area fraction)
// For future sensitivity tests
static float _MDparCleanTransitionYear     = 9999.0; // If present - year in which a switch to a different deicer concentration occurs
static float _MDparCleanTransitionConc     = 0.0; //  If present - the new concentration of deicer.

// Loading calculation
static int _MDInSubFractionID          = MFUnset;
static int _MDInAgFractionID         = MFUnset;
static int _MDInImpFracSpatialID       = MFUnset;
static int _MDInH2OFracSpatialID       = MFUnset;
static int _MDInPopulationDensityID    = MFUnset;
static int _MDInAtmChlorideID           = MFUnset;
static int _MDInPnET_Imp_MeltID         = MFUnset;
static int _MDInPnET_Imp_StormflowID    = MFUnset;

// output  - SC = Specific Conductance (uS / cm)
static int _MDOutLocalLoadSCID      = MFUnset;
static int _MDOutFlux_SCID          = MFUnset; 
static int _MDOutStoreWater_SCID    = MFUnset;
static int _MDOutPostConc_SCID      = MFUnset;
static int _MDOutSurfaceRunoffPool_SCID = MFUnset;

// output - Cl - Chloride concentration (mg/L) or mass (kg)
static int _MDOutClSnowPackID       = MFUnset;
static int _MDOutClRootZoneID       = MFUnset;
static int _MDOutClGrdWatID         = MFUnset;
static int _MDOutClSurfaceRunoffPoolID  = MFUnset;
static int _MDOutClFluxID           = MFUnset;
static int _MDOutClStoreID          = MFUnset;
static int _MDOutLocalLoadClID      = MFUnset;
static int _MDOutPostConc_ClID      = MFUnset;
static int _MDOutConcClgwPreID       = MFUnset;
static int _MDOutConcClimmPreID      = MFUnset;
static int _MDOutCldeicerInputID     = MFUnset;
static int _MDOutCltotalInputID      = MFUnset;


static void _MDSpecCond (int itemID) {

    float stormflowVol                  = 0.0;
    float baseflowVol                   = 0.0;
    float surfaceRunoffVol              = 0.0;
    float runoffpoolreleaseVol          = 0.0;
    float surfaceRunoffPool             = 0.0;
    float discharge                     = 0.0;
    float waterStorage                  = 0.0;
    
    // New Variables // 
    float Developed               = 0.0;
    float preFlux_SC              = 0.0;
    float postFlux_SC             = 0.0;
    float storeWater_SC           = 0.0;
    float postStoreWater_SC       = 0.0;
    float SCTotalIn               = 0.0;
    float postConc_SC             = 0.0;
    float massBalance_SC          = 0.0;  
    
    float localLoad_SC            = 0.0;
    
    float baseflow_SC           = 0.0; // Ionic content of baseflow (Development dependent) (uS/cm)
    float stormflow_SC           = 0.0; // Ionic content of stormflow (bf/Ta dependent) (uS/cm)
    float soilflow_SC             = 0.0; // Ionic content of soil runoff to the SroPool from soil (uS/cm)
    float surfflowPool_SC         = 0.0; // Ionic content of RunoffPool (Accumulates loading from HCIA Snowfall) (uS/cm)
    float clean_intrcp            = 0.0; // Clean ionic content of precipitation 
    float dev_bf_SC_rate          = 0.0; // Rate of ionic content increase based on development (uS/cm / % dev)
    float airTemperature          = 0.0; // degC
    float snowTreatTemperature    = 0.0; // Mean daily air Temperature where treatment occurs
    float soilflowDiluteFactor    = 0.0; // Dilution factor (against bf ionic content) for surface runoff to detention pool
    float stormflowDiluteFactor   = 0.0; // Dilution factor (against bf ionic content) for stormflow pool
    float stormflowTreatFactor    = 0.0;
        
    baseflowVol          = MFVarGetFloat (_MDInBaseFlowID,            itemID,0.0) * MFModelGetArea (itemID) / (MFModelGet_dt () * 1000.0); // m3/sec
    runoffpoolreleaseVol = MFVarGetFloat (_MDInRunoffPoolReleaseID, itemID,0.0) * MFModelGetArea (itemID) / (MFModelGet_dt () * 1000.0); // m3/sec
    stormflowVol         = MFVarGetFloat (_MDInStormRunoffTotalID,   itemID,0.0) * MFModelGetArea (itemID) / (MFModelGet_dt () * 1000.0); // m3/sec
    surfaceRunoffVol     = MFVarGetFloat (_MDInSurfaceRunoffID,      itemID,0.0) * MFModelGetArea (itemID) / (MFModelGet_dt () * 1000.0); //m3/sec
    surfaceRunoffPool    = MFVarGetFloat (_MDInSurfaceRunoffPoolID, itemID, 0.0) * MFModelGetArea (itemID) / 1000.0; // m3
    surfflowPool_SC      = MFVarGetFloat (_MDOutSurfaceRunoffPool_SCID,itemID,0.0); // uS/cm in SROpool 

    discharge            = MFVarGetFloat (_MDInDischargeID,          itemID, 0.0); // m3/sec, discharge leaving the grid cell, after routing!
    waterStorage         = MFVarGetFloat (_MDInRiverStorageID,       itemID, 0.0); // m3/sec, storage rate remaining in grid cell at end of timestep - routed to retention
    airTemperature       = MFVarGetFloat (_MDInAirTemperatureID,     itemID, 0.0); // degrees C
    // New Variables //
    Developed            = MFVarGetFloat (_MDInSubFractionID,        itemID, 0.0);  // proportion developed land

// DEFINE IONIC CONTENT:  ic = m3*(uS/cm)
    // Assumes that ionic strength behaves linearly and conservatively.  Strictly this is only valid when
    // ionic strength is controlled by a few conservative ions.  Therefore, this is valid for salt impacted
    // streams; however, this breaks down at low ionic strength waters.  Therefore - this deep in the chloride analysis
    // we are operating with the understanding that we are only really looking at chloride impacts - not chloride itself. - SZ 6/12/2014
    preFlux_SC           = MFVarGetFloat (_MDOutFlux_SCID,        itemID, 0.0); // ic/day 
    storeWater_SC        = MFVarGetFloat (_MDOutStoreWater_SCID,  itemID, 0.0); // ic/day       
   
    /// BASEFLOW LOADING OF SPECIFIC CONDUCTANCE  ///
    ///  attributed to the groundwater and surface runoff pools - impervious runof is attributed with clean water
    //baseflowConc_SC       = 25.94 + 13.32 * (Developed); // Baseflow_SC relation (0.9 Quantile based) to Developed area.  Developed should already be in PER
    clean_intrcp = 21.72; //25.94; //;//; 
    dev_bf_SC_rate = 9.602; // 13.32; //
    snowTreatTemperature =  _MDSnowFallThreshold;
    soilflowDiluteFactor = 0.64;
    //stormflowDiluteFactor = 0.2;
    stormflowTreatFactor = 8.0;
    baseflow_SC = clean_intrcp + dev_bf_SC_rate * (Developed); // Baseflow SC relation (regression to zero storm flow) to Devloped area.
    
    stormflow_SC = (airTemperature < snowTreatTemperature) ? stormflowTreatFactor*Developed/100.*baseflow_SC : MDMaximum(soilflowDiluteFactor * baseflow_SC,0.8*clean_intrcp); // Winter-time salting.
    soilflow_SC = MDMaximum( soilflowDiluteFactor * baseflow_SC ,1.5*clean_intrcp);
    surfflowPool_SC = (surfflowPool_SC * surfaceRunoffPool + soilflow_SC * surfaceRunoffVol*MFModelGet_dt() + stormflow_SC * stormflowVol*MFModelGet_dt()) / (surfaceRunoffPool + surfaceRunoffVol*MFModelGet_dt() + stormflowVol*MFModelGet_dt());  
    
    // Calculate input loadings
    localLoad_SC   = (baseflow_SC * baseflowVol + surfflowPool_SC * (runoffpoolreleaseVol + stormflowVol) ) * 86400 ; // ic/day
    
    SCTotalIn      = localLoad_SC + preFlux_SC  + storeWater_SC;         // ic/day  
    SCTotalIn = discharge <= 0.0000001 ? 0.0 : SCTotalIn;
    postConc_SC = discharge <= 0.0000001 ? 0.0 : SCTotalIn / (discharge * 86400.); // uS/cm : eliminate ionic mass flux on negiglible flow days
    //Calculates the fluxes/storage for the mixing (appropriate for speccond) case
    postFlux_SC  = (discharge * MDConst_m3PerSecTOm3PerDay) * postConc_SC ;         // ic/day
    postStoreWater_SC  = (waterStorage ) * postConc_SC ;      // ic/day
    
    // Calculates the mass balances
    massBalance_SC  = (SCTotalIn - (postFlux_SC + postStoreWater_SC ))/MDMaximum(SCTotalIn,0.0000001*baseflow_SC);
    
    // Print mass balance errors
    if (MFDateGetCurrentYear() > 0){
        if ( (massBalance_SC > 0.001)) {
           printf("itemID = %d, %d-%d-%d, MB_SC = %f\n",itemID, MFDateGetCurrentYear(), MFDateGetCurrentMonth(), MFDateGetCurrentDay(), massBalance_SC);
           printf("\tTotalIn: %f, Local: %f , UpFlux: %f, InStore: %f\n",SCTotalIn,localLoad_SC,preFlux_SC,storeWater_SC);
           printf("\tDownFlux: %f, discharge; %f, OutStore: %f , storage: %f\n",postFlux_SC,discharge, postStoreWater_SC,waterStorage);
        }
    }
    
    // Set Output
    MFVarSetFloat (_MDOutFlux_SCID,             itemID, postFlux_SC);
    MFVarSetFloat (_MDOutSurfaceRunoffPool_SCID,       itemID, surfflowPool_SC);
    MFVarSetFloat (_MDOutLocalLoadSCID,         itemID, localLoad_SC);
    MFVarSetFloat (_MDOutPostConc_SCID,         itemID, postConc_SC);
    MFVarSetFloat (_MDOutStoreWater_SCID,       itemID, postStoreWater_SC);
}

static int _treat_month(int month) {
    int _treating_months[6] = { 1, 2, 3, 4, 11, 12 };
    int i = 0;
    for (i = 0;i<=5;i++ ){
        if (month == _treating_months[i]) return 1;
    }
    return 0;
}

static void _MDTotalDissIons(int itemID){
    // NOTE: NAMES ARE INAPPROPRIATE NOW - USING AS TRACKING CHLORIDE SPECIFICALLY
    // Water balance variables
    float precipVol                  = 0.0; // m3/s
    float stormflowVol               = 0.0; // m3/s
    float percolation                = 0.0; // mm/d
    float recharge                   = 0.0; // mm/d
    float groundWater                = 0.0; // mm
    float baseflow                   = 0.0; // mm/d
    float soilMoist                  = 0.0; // mm
    float infiltration               = 0.0; // mm/d
    float surfRunoff                 = 0.0; // mm/d
    float runoffpoolrelease          = 0.0; // mm/d
    float surfaceRunoffPool          = 0.0; // mm
    float snowPack                   = 0.0; // mm
    float snowPackChange             = 0.0; // mm/d
    float discharge                  = 0.0; // m3/s
    float waterStorage               = 0.0; //m3/d
    float waterStorageChg	     = 0.0; //m3/d
    float riverbedVelocity	     = 0.0; // m/s
    float sinuosity		     = 0.0; // m/m
    float immGWExchange_rate            = _MDparGreaterGroundWaterExchange; // Immobile zone exchange rate (d^-1)
    float immGWbeta                     = 0.00; // Relative size of immobile zone (L^3 / L^3)
    float immFlux                       = 0.0; // kg / day (net exchange with immobile zone)
    float betaTerm                  =   1.0; // Excess dilution resulting from immobile zone exchange
    float Conc_GW_last              = 0.0; // Concentration Cl (kg/m3) last timestep
    float Conc_imm_last             = -1.0; // Concentration Cl (kg/m3) last timestep in immobile zone)
    
    // Cl WBM Variables // 
    float PassiveSoilStore		= 10.0; // mm
    float ConcClPrecipitation    = 0.0; // Cl mg /L
    float postConc_SC             = 0.0; // uS/cm
    float ConcClTreat            = 0.0; // Cl mg/L
    float FluxImp2Runoff          = 0.0; // kg/d
    float snowFallFromImp2SPack   = 0.0; // mm/d
    float snowFallFromSky2SPack  = 0.0; // mm/d
    float FluxImp2SnowPack        = 0.0; // kg/d
    float FluxSky2SnowPack      = 0.0; // kg/d
    float MassClSnowPackPre    = 0.0; // kg Cl
    float MassClSnowPack       = 0.0; // kg Cl
    float ConcClSnowPack       = 0.0; // mg/L Cl
    float FluxMelt              = 0.0; // kg/d Cl
    float FluxPrecip            = 0.0; // kg/d Cl
    float MassClRootZonePre    = 0.0; // kg Cl
    float MassClRootZone       = 0.0; // kg Cl
    float ConcClRootZone       = 0.0; // mg/L Cl
    float FluxPercolation       = 0.0; // kg/d
    float FluxRecharge          = 0.0; // kg /d
    float FluxSurfFlow           = 0.0; // kg/d
    float FluxWeathering        = 0.0; // kg/d
    float FluxDevel              = 0.0; // kg/d
    float FluxAg                = 0.0; // kg/d
    float MassClGroundWaterPre =  0.0; // kg Cl
    float MassClGroundWater     = 0.0; // kg Cl
    float ConcClGroundWater     = 0.0; // mg/L
    float FluxBaseFlow          = 0.0; // kg/d
    float MassClSurfROPoolPre  = 0.0; // kg Cl
    float MassClSurfROPool     = 0.0; // kg Cl
    float ConcClSurfROPool     = 0.0; // mg/L
    float FluxSurfROPoolRelease = 0.0; // kg/d
    float FluxH2O               = 0.0; // kg/d
    // Cl WTM Variables //
    float ClTotalIn            = 0.0; // kg/d
    float local_load            = 0.0; // kg/d
    float preFlux_Cl           = 0.0; // kg/d
    float storeWater_Cl        = 0.0; // kg/d
    float postConc_Cl          = 0.0; // mg/L
    float postFlux_Cl          = 0.0; // kg/d
    float postStoreWater_Cl    = 0.0; // kg/d
    float massBalance_Cl       = 0.0; // kg/d
    float C0		       = 0.0; // Linear Routing parameter


    // Misc
    float airTemperature         = 0.0; // degC
    float imp                    = 0.0; // fraction of grid cell impervious cover
    float dev                    = 0.0; // fraction of grid cell developed cover
    float h2o                   = 0.0; // fraction of grid cell open water
    float pop                   = 0.0; // Population density (person/km2)
    float ag                    = 0.0; // fraction of grid cell agriculture cover
    float hcia                   = 0.0; // fraction of impervious - directly hydrologically connected
    float whcia                  = 0.0; /// fraction of impervious - directly hydrologically connected for snowfall
    float dt                     = 1.0; // Timestep is one 1.0 day
    float cleanYear              = 0.0; // Year to switch deicer loading concentration
    float cleanDeicer            = 0.0; // Deicer concentration to switch to
    int month                    = MFDateGetCurrentMonth(); 
    int year                     = MFDateGetCurrentYear(); 
    
    /********************* WBM:  Calculate Cell Ion Balance  *********************/
    precipVol                  = MFVarGetFloat(_MDInPrecipitationID,    itemID, 0.0) * MFModelGetArea(itemID) / (MFModelGet_dt() * 1000.0); // m3/s
    stormflowVol               = MFVarGetFloat(_MDInStormRunoffTotalID, itemID, 0.0) * MFModelGetArea(itemID) / (MFModelGet_dt() * 1000.0); // m3/s
    snowPack                   = MFVarGetFloat(_MDInSnowPackID,     itemID, 0.0); // mm
    snowPackChange             = MFVarGetFloat(_MDInSPackChgID,itemID,0.0); // mm/d
    soilMoist                  = MFVarGetFloat(_MDInSoilMoistID,    itemID, 0.0); // mm
    PassiveSoilStore           = MFVarGetFloat(_MDInAvailableWaterCapacityID, itemID, 0.0) * _MDparPassiveSoilStoreFactor; // mm
//    wiltMoist                  = wiltMoist > 10 ? wiltMoist : 10; // set minimum of 10 mm for wilting capacity
    recharge                   = MFVarGetFloat(_MDInGrdWatRechargeID,     itemID, 0.0); // mm/d
    groundWater                = MFVarGetFloat(_MDInGrdWatID,       itemID, 0.0) + _MDparGreaterGroundWaterMixing; // mm
    immGWbeta                  = _MDparGreaterGroundWaterMixing * 0.001 * MFModelGetArea(itemID);
    baseflow                   = MFVarGetFloat(_MDInBaseFlowID,     itemID, 0.0); // mm/d
    surfRunoff                 = MFVarGetFloat(_MDInSurfaceRunoffID,itemID,0.0); // mm/d
    runoffpoolrelease          = MFVarGetFloat(_MDInRunoffPoolReleaseID,itemID,0.0); // mm/d
    surfaceRunoffPool          = MFVarGetFloat(_MDInSurfaceRunoffPoolID,itemID,0.0); // mm
    ConcClPrecipitation        = MFVarGetFloat(_MDInAtmChlorideID,itemID,0.0); // mg/L
    MassClSnowPack            = MFVarGetFloat(_MDOutClSnowPackID, itemID,0.0); // kg Cl
    MassClRootZone            = MFVarGetFloat(_MDOutClRootZoneID, itemID,0.0); // kg Cl
    MassClGroundWater         = MFVarGetFloat(_MDOutClGrdWatID,   itemID,0.0); // kg Cl
    MassClSurfROPool          = MFVarGetFloat(_MDOutClSurfaceRunoffPoolID,itemID,0.0); // kg Cl
    
    Conc_GW_last               = MFVarGetFloat(_MDOutConcClgwPreID,  itemID,     0.0); // mg Cl / L
    Conc_imm_last              = MFVarGetFloat(_MDOutConcClimmPreID, itemID,     0.0); // mg Cl/ L
    
    airTemperature             = MFVarGetFloat (_MDInAirTemperatureID,     itemID, 0.0); // degrees C
    imp                        = MFVarGetFloat (_MDInImpFracSpatialID,        itemID, 0.0);  // proportion Impervious Cover
    dev                        = MFVarGetFloat (_MDInSubFractionID,       itemID, 0.0); // proportion developed cover
    h2o                        = MFVarGetFloat (_MDInH2OFracSpatialID,     itemID,0.0); // proportion open water
    pop                        = MFVarGetFloat (_MDInPopulationDensityID,   itemID, 0.0); // Population density
    ag                         = MFVarGetFloat (_MDInAgFractionID,       itemID, 0.0); // proportion developed cover
    hcia                       = MFVarGetFloat (_MDInHCIAID,               itemID, 0.0); // hydrologically connected impervious
    whcia                      = MFVarGetFloat (_MDInWinterHCIAID,          itemID, 0.0); /// hydrologically connected impervious in winter
    if (_MDInSoilPercolationID != MFUnset){
        percolation = MFVarGetFloat(_MDInSoilPercolationID,  itemID, 0.0); // mmd/
    } else {
        percolation = 0.0;
    }
    
    if (_MDInDeicerClID != MFUnset) _MDparDeicerCl = MFVarGetFloat(_MDInDeicerClID,   itemID, 1234.); // If Deicer is passed as datasource, load data

    if (_MDparCleanTransitionYear <= year) _MDparDeicerCl = _MDparCleanTransitionConc; 
    ConcClTreat = ((airTemperature <= _MDSnowFallThreshold) & (_treat_month(month))) ? _MDparDeicerCl : ConcClPrecipitation; // mg/L
    // Calculate direct impervious snowfall runoff  - On SnowFall days: only runoff is HCIA SnowFall.  On Non-snowffall days: ConcClTreat = 0.
    FluxImp2Runoff = ConcClTreat*(0.001)*(stormflowVol)*MFModelGet_dt(); // kg/day # This represents imp and h2o on warm days.
    
    
    /*Calculate accumulation of Cl in "snowpack" (really just a snow-packish pool) */
    // Snow Fall / Snowpack
    snowFallFromImp2SPack = snowPackChange > 0.0 ? snowPackChange : 0.0;
    snowFallFromImp2SPack = snowFallFromImp2SPack * MDMaximum((1-whcia-0.4),0.00)*imp; // Removing 40% of impervious area that are expected to be from roofs ((Southworth and Ben-Joseph 1996)
    snowFallFromSky2SPack = snowPackChange > 0.0 ? snowPackChange - snowFallFromImp2SPack : 0.0; // Remaining snowpack development (including on roofs) is clean precip
    snowFallFromImp2SPack = snowFallFromImp2SPack * MFModelGetArea (itemID) / (MFModelGet_dt () * 1000.0); // Water flux to snowpack from impervious areas
    snowFallFromSky2SPack = snowFallFromSky2SPack * MFModelGetArea (itemID) / (MFModelGet_dt () * 1000.0); // Water flux to snowpack directly from sky (and non HCIA impervious areas)
    FluxImp2SnowPack = ConcClTreat*(0.001)*(snowFallFromImp2SPack)*MFModelGet_dt(); // kg/day Cl to accumulate in snowpack deicing impervious areas
    FluxSky2SnowPack = ConcClPrecipitation*(0.001)*(snowFallFromSky2SPack)*MFModelGet_dt(); // kg/day Cl to accumulate in snowpack from clean precip
    MassClSnowPackPre = MassClSnowPack;
    MassClSnowPack = MassClSnowPack + FluxImp2SnowPack * dt + FluxSky2SnowPack * dt; // kg Cl 
    // Note: Do not need to calculate concentration using sum of snowPack+Melt, because Melt only occurs on non-snowfall days
    if (snowPackChange <= 0.0) { // melt day - calculate concentration pre-melt
        ConcClSnowPack = snowPack > 0.0 ? (MassClSnowPack ) / ((snowPack - snowPackChange) * MFModelGetArea(itemID) / 1000.) * (1000.) : 0.0; // mg/L Cl
    } else { // snowfall day - calculate concentration after snow-fall (unused)
        ConcClSnowPack = snowPack > 0.0 ? (MassClSnowPack ) / (snowPack * MFModelGetArea(itemID) / 1000.) * (1000.) : 0.0; // mg/L Cl
    }
    
    // Only pathway out of snow pack is melt to RootZone
    FluxMelt = snowPackChange < 0.0 ? ConcClSnowPack*(0.001)*fabs(snowPackChange)*MFModelGetArea(itemID) / ( 1000.0 ) : 0.0; // kg/day
    // Presumably due to rounding error .... Not all mass leaves snow-pack at the end of the melt season.
    // So flush all mass to FluxMelt if SnowPack == 0.0
    FluxMelt = snowPack == 0.0 ? MassClSnowPack : FluxMelt;
    // Update SnowPack Mass with melt.
    MassClSnowPack = MDMaximum(MassClSnowPack - FluxMelt * dt,0.0);
        
    // Clean precipitation flux
    infiltration = snowPackChange == 0 ? precipVol - stormflowVol : 0.0;
    FluxPrecip = ConcClPrecipitation * (0.001) * (infiltration)*MFModelGet_dt(); // kg/day
 
    // Split operator mass balance in root zone
    /* Calculate accumulation in Root Zone */
    MassClRootZonePre = MassClRootZone; // kg Cl
    MassClRootZone = MassClRootZone + FluxMelt*dt + FluxPrecip*dt; // kg Cl
    ConcClRootZone = (MassClRootZone ) / ((soilMoist + PassiveSoilStore+dt*(surfRunoff+recharge+percolation)) * MFModelGetArea(itemID) / 1000.) * (1000.); // mg/L Cl
    
    // Groundwater
    FluxRecharge = ConcClRootZone*(0.001)*recharge*MFModelGetArea(itemID) / ( 1000.) ; // kg/day Cl
    FluxPercolation = ConcClRootZone*(0.001)*percolation*MFModelGetArea(itemID) / (1000.) ; // kg/day Cl
    FluxSurfFlow = ConcClRootZone*(0.001)*surfRunoff*MFModelGetArea(itemID) / ( 1000.) ; // kg/day Cl
    // If Soil Moist goes to zero - then mass stays in the Root Zone.
    // Update RootZone Mass with fluxes out
    MassClRootZone = MDMaximum(MassClRootZone - FluxRecharge * dt - FluxSurfFlow * dt - FluxPercolation * dt, 0.0); // kg Cl

//    if ((Conc_GW_last == 0) & (MFDateGetCurrentYear() < 0)) {
//        // For first year of spin-up ignore immobile zone.  Calculate mean concentration across entire gw domain at end of first year.
//        // Then continue spinup and simulation with M-Im MT following Silva et al. 2008
//        MassClGroundWaterPre = MassClGroundWater; // kg  Cl
//        FluxWeathering = _MDparGeoChem_ClWx_Rate * (0.001)*groundWater*MFModelGetArea(itemID); // kg/day
//        FluxDevel = _MDparDevelop_ClWx_Rate*pop*MFModelGetArea(itemID)*1.0e-6; // kg / day
//        FluxAg    = _MDparAg_ClWx_Rate*ag/100.*MFModelGetArea(itemID); // kg /day
//        betaTerm = ((groundWater - _MDparGreaterGroundWaterMixing)*MFModelGetArea(itemID)*0.001);
//        MassClGroundWater = ( MassClGroundWaterPre + (FluxPercolation+FluxRecharge+FluxWeathering+FluxDevel+FluxAg) * dt ) ;  // Mass in GW (no ImmStor)
//        ConcClGroundWater = groundWater > 0.0 ? 1000.*MassClGroundWater / ((groundWater + baseflow*dt)*MFModelGetArea(itemID)*0.001) : 0.0; // Conc Cl in GW mg/L
//        FluxBaseFlow = 0.001 * (ConcClGroundWater)*(0.001)*baseflow*MFModelGetArea(itemID);
//        MassClGroundWater = MassClGroundWater - FluxBaseFlow * dt; 
//        // Last day of first year ... average concentration is applied throughout whole groundwater (mobile+immobile domain)
//        if ((MFDateGetCurrentMonth() == 12) & (MFDateGetCurrentDay() == 31)) {
//		Conc_imm_last = ConcClGroundWater;
//        } else {
//            ConcClGroundWater = 0.0;
//        }
   // } else {
        // Calculate accumulation in Ground Water
        MassClGroundWaterPre = MassClGroundWater; // kg  Cl
        FluxWeathering = _MDparGeoChem_ClWx_Rate * (0.001)*groundWater*MFModelGetArea(itemID); // kg/day
        FluxDevel = _MDparDevelop_ClWx_Rate*pop*MFModelGetArea(itemID)*1.0e-6; // kg / day  (population density given in persons/km2 )
        FluxAg    = _MDparAg_ClWx_Rate*ag/100.*MFModelGetArea(itemID); // kg /day
        // Calculate immobile zone fluxes (assume eulerian scheme)
        immFlux = immGWbeta * immGWExchange_rate *0.001* Conc_imm_last * exp(-immGWExchange_rate*dt) - immGWbeta * immGWExchange_rate * 0.001*Conc_GW_last * exp(-immGWExchange_rate * dt) ; 
        betaTerm = ((groundWater - _MDparGreaterGroundWaterMixing)*MFModelGetArea(itemID)*0.001)+immGWbeta*(1.-exp(-immGWExchange_rate * dt));
        ConcClGroundWater = 1000.0 * ( 0.001 * Conc_GW_last + (FluxRecharge+FluxWeathering+FluxDevel+immFlux)/betaTerm ) / (1 + (baseflow * MFModelGetArea(itemID)*0.001*dt) / betaTerm); // Conc Cl in gw mg/L
        Conc_imm_last = 1000.* (0.001*Conc_imm_last*exp(-immGWExchange_rate * dt) + 0.001*Conc_GW_last * (1. - exp (-immGWExchange_rate * dt)) + \
                                ((0.001*ConcClGroundWater - 0.001*Conc_GW_last )/dt) * ( 1. - (1./immGWExchange_rate)*(1. - exp(-immGWExchange_rate * dt)))); // mg /L
        // Only pathway out of groundwater is baseflow
        FluxBaseFlow = 0.001*ConcClGroundWater*(0.001)*baseflow*MFModelGetArea(itemID) ; // kg/day Cl
        MassClGroundWater = MDMaximum((0.001*ConcClGroundWater*0.001*(groundWater-_MDparGreaterGroundWaterMixing)+0.001*Conc_imm_last*0.001*_MDparGreaterGroundWaterMixing)*MFModelGetArea(itemID),0.0); // kg Cl
   // }
    // Calculate accumulation in Surface Runoff Pool
    MassClSurfROPoolPre = MassClSurfROPool; // kg Cl
    MassClSurfROPool = MassClSurfROPool + FluxSurfFlow * dt;
    ConcClSurfROPool = surfaceRunoffPool > 0.0 ? (MassClSurfROPool)/ ((surfaceRunoffPool + runoffpoolrelease) * MFModelGetArea(itemID) / 1000.) * (1000.0) : 0.0 ; // mg/L Cl
    // Only pathway out of surfROpool is SurfROPoolRelease
    FluxSurfROPoolRelease = ConcClSurfROPool * (0.001) * runoffpoolrelease*MFModelGetArea(itemID) / ( 1000.) ; // kg/day Cl
    FluxSurfROPoolRelease = surfaceRunoffPool == 0.0 ? MassClSurfROPool / dt : FluxSurfROPoolRelease;
    // Update Surface RO Pool mass with flux out
    MassClSurfROPool = MDMaximum(MassClSurfROPool - FluxSurfROPoolRelease * dt,0.0); // kg Cl

    MFVarSetFloat( _MDOutClSnowPackID , itemID, MassClSnowPack); 
    MFVarSetFloat( _MDOutClRootZoneID , itemID, MassClRootZone);
    MFVarSetFloat( _MDOutClGrdWatID , itemID, MassClGroundWater);
    MFVarSetFloat( _MDOutClSurfaceRunoffPoolID, itemID, MassClSurfROPool ); 
    MFVarSetFloat( _MDOutConcClgwPreID, itemID, ConcClGroundWater );
    MFVarSetFloat( _MDOutConcClimmPreID, itemID, Conc_imm_last ); 
    // Mass Balance sets
    MFVarSetFloat( _MDOutCldeicerInputID,itemID, FluxImp2Runoff+FluxImp2SnowPack);
    MFVarSetFloat( _MDOutCltotalInputID,itemID,FluxImp2Runoff+FluxImp2SnowPack+FluxSky2SnowPack+FluxPrecip+FluxWeathering+FluxDevel+FluxAg);
    
    /********************* WTM:  Calculate Network Ion Balance  *********************/
    discharge                  = MFVarGetFloat(_MDInDischargeID,    itemID,0.0); // m3/s
    waterStorage               = MFVarGetFloat(_MDInRiverStorageID, itemID,0.0); // m3/s
    waterStorageChg	      = MFVarGetFloat(_MDInRiverStorageChgID, itemID,0.0); //m3/s
    preFlux_Cl                = MFVarGetFloat(_MDOutClFluxID,     itemID,0.0); // kg/d
    storeWater_Cl             = MFVarGetFloat(_MDOutClStoreID,    itemID,0.0); // kg/d
        
    local_load = (FluxImp2Runoff + FluxBaseFlow + FluxSurfROPoolRelease);  // kg/day Cl
        
    ClTotalIn = local_load + preFlux_Cl  + storeWater_Cl;    // kg/day Cl
    // Updated 2016-12-10 to calculate cascade routing

    sinuosity = MFVarGetFloat(_MDInSinuosityID,itemID,0.0); // 
    float dL = MFModelGetLength(itemID); //
    riverbedVelocity = MFVarGetFloat(_MDInRiverbedVelocityMeanID,itemID,0.22); //

    C0 = riverbedVelocity != riverbedVelocity ? 1.0 : 1.0 / (1.0 + (sinuosity*dL/riverbedVelocity/86400.0)); //

    postFlux_Cl = discharge > 0.00 ? C0 * ClTotalIn : 0.0;
    postFlux_Cl = postFlux_Cl < 0.0 ? local_load + preFlux_Cl : postFlux_Cl;
    float storageCl_Chg  = 0.0; //
    postStoreWater_Cl = ClTotalIn * (1.0 - C0); 
    if (postStoreWater_Cl < 0.0) {
	storageCl_Chg = -storeWater_Cl;
	postStoreWater_Cl = 0.0; }
    else {
	    storageCl_Chg = postStoreWater_Cl - storeWater_Cl;
    }

   // ClTotalIn = discharge <= 0.0000001 ? 0.0 : ClTotalIn;   // kg/day Cl // eliminate ionic mass flux on negiglible flow days
    postConc_Cl = discharge > 0.0 ? (1000.) * postFlux_Cl / (discharge * 86400.) : 0.0; // mg/L Cl 
    //Calculates the fluxes/storage for the mixing (appropriate for speccond) case
    //postFlux_Cl  = (discharge * MDConst_m3PerSecTOm3PerDay) * (0.001) * postConc_Cl ;         // kg/day Cl
    //postStoreWater_Cl  = (waterStorage * MDConst_m3PerSecTOm3PerDay ) * (0.001) * postConc_Cl ;      // kg/day Cl
    
    // Calculates the mass balances
    massBalance_Cl  = (ClTotalIn - (postFlux_Cl + postStoreWater_Cl ));
  

    if ( ClTotalIn == 999.999e3 ) {
	printf("#### itemID = %d (%d-%d-%d) LocalLoad: %f preFlux_Cl: %f storeWaterCl: %f \n \
			MassClSnowPack: %f MassClRootZone: %f MassClGroundWater: %f MassClSurROPool: %f ConcCLGroundWater: %f ConcImmLast: %f \n \
			SnowPack: %f SnowPackChg: %f RiverbedVelocity: %f dL: %f C0: %f \n", itemID,MFDateGetCurrentYear(),MFDateGetCurrentMonth(),MFDateGetCurrentDay(),local_load,preFlux_Cl,storeWater_Cl, \
			MassClSnowPack,MassClRootZone,MassClGroundWater,MassClSurfROPool,ConcClGroundWater,Conc_imm_last,snowPack,snowPackChange,riverbedVelocity,dL,C0);
    }


    // Print mass balance errors
    if (MFDateGetCurrentYear() > 0) {
        if  ( (fabs(massBalance_Cl) > 0.001) && ((fabs(massBalance_Cl)/ClTotalIn) > 0.001) ) {
	   printf("itemID = %d, %d-%d-%d, MB_SC = %f\n",itemID, MFDateGetCurrentYear(), MFDateGetCurrentMonth(), MFDateGetCurrentDay(), massBalance_Cl);
           printf("\tTotalIn: %g, Local: %g , UpFlux: %g, InStore: %g CO: %f v: %f dL: %f \n",ClTotalIn,local_load,preFlux_Cl,storeWater_Cl,C0,riverbedVelocity,dL);
           printf("\tDownFlux: %g, discharge; %g, OutStore: %g , storage: %g\n",postFlux_Cl,discharge, postStoreWater_Cl,waterStorage);
        } 
   }
    
    // Set Output
    MFVarSetFloat (_MDOutClFluxID,        itemID, postFlux_Cl);
    MFVarSetFloat (_MDOutClStoreID,       itemID, postStoreWater_Cl);
    MFVarSetFloat (_MDOutLocalLoadClID,   itemID, local_load);
    MFVarSetFloat (_MDOutPostConc_ClID,   itemID, postConc_Cl);
    
    float b = 31.2; // EPSCoR statewide  Trowbridge: -22.0 ... Novotny: -37.25 (inverse)
    float m = 3.96; // EPSCoR statewide   Trowbridge: 0.307      Novotny: 0.25 (inverse)
    postConc_SC = MDMaximum(m * postConc_Cl + b,0.0);
    
    MFVarSetFloat (_MDOutPostConc_SCID,       itemID, postConc_SC); // uS/cm
}

static void _MDChloride_PnET(int itemID){
    // Water balance variables
    float stormflowVol               = 0.0; // m3/s
    float Imp_stormflow              = 0.0; // mm/d
    float recharge                   = 0.0; // mm/d
    float groundWater                = 0.0; // mm
    float baseflow                   = 0.0; // mm/d
    float surfRunoff                 = 0.0; // mm/d
    float runoffpoolrelease          = 0.0; // mm/d
    float surfaceRunoffPool          = 0.0; // mm
    float ImpMelt                   = 0.0; // mm/d
    float ImpRain                   = 0.0; // mm/d
    float discharge                  = 0.0; // m3/s
    float waterStorage               = 0.0; //m3/s
    
    float immGWExchange_rate            = _MDparGreaterGroundWaterExchange; // Immobile zone exchange rate (d^-1)
    float immGWbeta                     = 0.00; // Relative size of immobile zone (L^3 / L^3)
    float immFlux                       = 0.0; // kg / day (net exchange with immobile zone)
    float betaTerm                  =   1.0; // Excess dilution resulting from immobile zone exchange
    float Conc_GW_last              = 0.0; // Concentration Cl (kg/m3) last timestep
    float Conc_imm_last             = -1.0; // Concentration Cl (kg/m3) last timestep in immobile zone)
    
    // Cl WBM Variables // 
    float Precip_Chloride_mean    = 0.270; // mg Cl /L 
    float postConc_SC             = 0.0; // uS/cm
    float ConcClTreat            = 0.0; // Cl mg/L
    float FluxImp2Runoff          = 0.0; // kg/d
    float snowFallFromImp2SPack   = 0.0; // mm/d
    float snowFallFromSky2SPack  = 0.0; // mm/d
    float FluxImp2SnowPack        = 0.0; // kg/d
    float FluxSky2SnowPack      = 0.0; // kg/d
    float MassClSnowPackPre    = 0.0; // kg Cl
    float MassClSnowPack       = 0.0; // kg Cl
    float ConcClSnowPack       = 0.0; // mg/L Cl
    float FluxMelt              = 0.0; // kg/d Cl
    float FluxPrecip            = 0.0; // kg/d Cl
    float MassClRootZonePre    = 0.0; // kg Cl
    float MassClRootZone       = 0.0; // kg Cl
    float ConcClRootZone       = 0.0; // mg/L Cl
    float FluxPercolation       = 0.0; // kg/d
    float FluxRecharge          = 0.0; // kg /d
    float FluxSurfFlow           = 0.0; // kg/d
    float FluxWeathering        = 0.0; // kg/d
    float FluxDevel              = 0.0; // kg/d
    float FluxAg                = 0.0; // kg/d
    float MassClGroundWaterPre =  0.0; // kg Cl
    float MassClGroundWater     = 0.0; // kg Cl
    float ConcClGroundWater     = 0.0; // mg/L
    float FluxBaseFlow          = 0.0; // kg/d
    float MassClSurfROPoolPre  = 0.0; // kg Cl
    float MassClSurfROPool     = 0.0; // kg Cl
    float ConcClSurfROPool     = 0.0; // mg/L
    float FluxSurfROPoolRelease = 0.0; // kg/d
    float FluxH2O               = 0.0; // kg/d
    // Cl WTM Variables //
    float ClTotalIn            = 0.0; // kg/d
    float local_load            = 0.0; // kg/d
    float preFlux_Cl           = 0.0; // kg/d
    float storeWater_Cl        = 0.0; // kg/d
    float postConc_Cl          = 0.0; // mg/L
    float postFlux_Cl          = 0.0; // kg/d
    float postStoreWater_Cl    = 0.0; // kg/d
    float massBalance_Cl       = 0.0; // kg/d

    // Misc
    float airTemperature         = 0.0; // degC
    float awc                    = 0.0; // mm
    float imp                    = 0.0; // fraction of grid cell impervious cover
    float dev                    = 0.0; // fraction of grid cell developed cover
    float h2o                   = 0.0; // fraction of grid cell open water
    float pop                   = 0.0; // Population density (person/km2)
    float ag                    = 0.0; // fraction of grid cell agriculture cover
    float hcia                   = 0.0; // fraction of impervious - directly hydrologically connected
    float whcia                  = 0.0; /// fraction of impervious - directly hydrologically connected for snowfall
    float dt                     = 1.0; // Timestep is one 1.0 day
    float CarHabitat            = 0.6; // Fraction of impervious areas getting treatment - see Shan for sources
    int month                    = MFDateGetCurrentMonth(); 
    
    /********************* WBM:  Calculate Cell Ion Balance  *********************/
    stormflowVol               = MFVarGetFloat(_MDInStormRunoffTotalID, itemID, 0.0); // mm/d
    Imp_stormflow              = MFVarGetFloat(_MDInPnET_Imp_StormflowID,itemID,0.0); // mm/d
    ImpMelt                    = MFVarGetFloat(_MDInPnET_Imp_MeltID,        itemID,0.0); // mm/d (terms of direct impervious area only) // 20160504-Changed to whole imperv area
    awc                        = MDMaximum( MFVarGetFloat(_MDInAvailableWaterCapacityID, itemID, 0.0), 10.0); // mm Force to 100mm
    recharge                   = MFVarGetFloat(_MDInGrdWatRechargeID,     itemID, 0.0); // mm/d
    groundWater                = MFVarGetFloat(_MDInGrdWatID,       itemID, 0.0) + _MDparGreaterGroundWaterMixing; // mm
    immGWbeta                  = _MDparGreaterGroundWaterMixing * 0.001 * MFModelGetArea(itemID);
    baseflow                   = MFVarGetFloat(_MDInBaseFlowID,     itemID, 0.0); // mm/d
    surfRunoff                 = MFVarGetFloat(_MDInSurfaceRunoffID,itemID,0.0); // mm/d
    runoffpoolrelease          = MFVarGetFloat(_MDInRunoffPoolReleaseID,itemID,0.0); // mm/d
    surfaceRunoffPool          = MFVarGetFloat(_MDInSurfaceRunoffPoolID,itemID,0.0); // mm
    //MassClSnowPack            = MFVarGetFloat(_MDOutClSnowPackID, itemID,0.0); // kg Cl
    MassClRootZone            = MFVarGetFloat(_MDOutClRootZoneID, itemID,0.0); // kg Cl
    MassClGroundWater         = MFVarGetFloat(_MDOutClGrdWatID,   itemID,0.0); // kg Cl
    MassClSurfROPool          = MFVarGetFloat(_MDOutClSurfaceRunoffPoolID,itemID,0.0); // kg Cl
   
    Conc_GW_last               = MFVarGetFloat(_MDOutConcClgwPreID,  itemID,     0.0); // mg Cl / L
    Conc_imm_last              = MFVarGetFloat(_MDOutConcClimmPreID, itemID,     0.0); // mg Cl/ L
    
    airTemperature             = MFVarGetFloat (_MDInAirTemperatureID,     itemID, 0.0); // degrees C
    imp                        = MFVarGetFloat (_MDInImpFracSpatialID,        itemID, 0.0);  // proportion Impervious Cover
    hcia                       = MFVarGetFloat (_MDInHCIAID,        itemID, 0.0);  // proportion Impervious Cover hydrologically/directly connected to waterbodies
    whcia                      = MFVarGetFloat (_MDInWinterHCIAID,          itemID, 0.0); /// hydrologically connected impervious in winter
    dev                        = MFVarGetFloat (_MDInSubFractionID,       itemID, 0.0); // proportion developed cover
    h2o                        = MFVarGetFloat (_MDInH2OFracSpatialID,     itemID,0.0); // proportion open water
    pop                        = MFVarGetFloat (_MDInPopulationDensityID,   itemID, 0.0); // Population density
    ag                         = MFVarGetFloat (_MDInAgFractionID,       itemID, 0.0); // proportion developed cover

    
    if (_MDInDeicerClID != MFUnset) _MDparDeicerCl = MFVarGetFloat(_MDInDeicerClID,   itemID, 1234.); // If Deicer is passed as datasource, load data
   
    float Whole_ImpMelt = 0.0; // Melt rate (mm/d) across whole impervious area
    float Flux_WholeImp = 0.0; // Chloride flux (kg/d) from whole impervious area from deicer application
    float Fraction_Runoff = 0.0; // Fraction of total deicer chloride applied in melt that is directed to storm runoff (kg/d)
    
//    Whole_ImpMelt = hcia > 0.0 ? hcia > 0.15 ? ImpMelt / 0.2 : ImpMelt / 0.2 : 0.0; // mm /d : Zaixing's claim: The amount of melt scaled to the whole impervious area (PnET passes only melt on direct impervious)
// I do not think that is true.  I'm trying without dividing through by / HCIA
    Whole_ImpMelt = ImpMelt; // Requested that Zaixing redifine the output variable as in terms of whole impervious area.
    Flux_WholeImp = Whole_ImpMelt * MFModelGetArea(itemID) * imp * CarHabitat * _MDparDeicerCl * 1.0e-6;// convert from mg/d to kg/d
    Fraction_Runoff = MDMinimum(hcia,CarHabitat)/CarHabitat; // Fraction of treated melt water directly to runoff

    FluxImp2Runoff = Flux_WholeImp * Fraction_Runoff; // kg Cl / d Total storm runoff flux from HC impervious and open water
       
    // We only track mass of melt from non hydrologically connected CarHabitat impervious areas to root zone
    FluxMelt = Flux_WholeImp - FluxImp2Runoff; // kg/day Dirty melt to pervious
    
//    if ((itemID == 836)|(itemID == 282)) printf("%d hcia: %f Deicer %f ImpMelt %f Whole Melt Imp %f Flux WholeImp %f FracRunoff %f FluxImp2RO %f FluxMelt %f \n",itemID,hcia,_MDparDeicerCl,  ImpMelt, Whole_ImpMelt, Flux_WholeImp,Fraction_Runoff, FluxImp2Runoff,FluxMelt); //
    /* Calculate accumulation in Root Zone */
    MassClRootZonePre = MassClRootZone; // kg Cl
    MassClRootZone = MassClRootZone + FluxMelt*dt; // kg Cl
    ConcClRootZone = MassClRootZone / (0.001*(awc+surfRunoff+recharge)* MFModelGetArea(itemID)) * 1000. ; // kg / m3 to mg/L Cl (From just Deicer)
    
    // Groundwater (Remember this ... we're assuming that all water has an additional concentration from precipitation)
    FluxRecharge = (ConcClRootZone)*(0.001)*recharge*MFModelGetArea(itemID) / ( 1000.) ; // kg/day Cl
    FluxSurfFlow = (ConcClRootZone)*(0.001)*surfRunoff*MFModelGetArea(itemID) / ( 1000.) ; // kg/day Cl
    float FluxPrecip2GW = Precip_Chloride_mean*0.001*recharge*MFModelGetArea(itemID) / 1000.; // kg/day Cl 
    float FluxPrecip2SW = Precip_Chloride_mean*0.001*surfRunoff*MFModelGetArea(itemID) / 1000.; // kg/day Cl
    
    if ((Conc_GW_last == 0) & (MFDateGetCurrentYear() < 0)) {
        // For first year of spin-up ignore immobile zone.  Calculate mean concentration across entire gw domain at end of first year.
        // Then continue spinup and simulation with M-Im MT following Silva et al. 2008
        MassClGroundWaterPre = MassClGroundWater; // kg  Cl
        FluxWeathering = _MDparGeoChem_ClWx_Rate * (0.001)*groundWater*MFModelGetArea(itemID); // kg/day
        FluxDevel = _MDparDevelop_ClWx_Rate*pop*MFModelGetArea(itemID)*1.0e-6; // kg / day
        FluxAg    = _MDparAg_ClWx_Rate*ag/100.*MFModelGetArea(itemID); // kg /day
        betaTerm = ((groundWater - _MDparGreaterGroundWaterMixing)*MFModelGetArea(itemID)*0.001);
        MassClGroundWater = ( MassClGroundWaterPre + (FluxPrecip2GW+FluxRecharge+FluxWeathering+FluxDevel+FluxAg) * dt ) ;  // Mass in GW (no ImmStor)
        ConcClGroundWater = groundWater > 0.0 ? 1000.*MassClGroundWater / ((groundWater + baseflow*dt)*MFModelGetArea(itemID)*0.001) : 0.0; // Conc Cl in GW mg/L
        FluxBaseFlow = 0.001 * (ConcClGroundWater)*(0.001)*baseflow*MFModelGetArea(itemID);
        MassClGroundWater = MassClGroundWater - FluxBaseFlow * dt; 
        // Last day of first year ... average concentration is applied throughout whole groundwater (mobile+immobile domain)
        if ((MFDateGetCurrentMonth() == 12) & (MFDateGetCurrentDay() == 31)) {
		Conc_imm_last = ConcClGroundWater;
        } else {
            ConcClGroundWater = 0.0;
        }
    } else {
        // Calculate accumulation in Ground Water
        MassClGroundWaterPre = MassClGroundWater; // kg  Cl
        FluxWeathering = _MDparGeoChem_ClWx_Rate * (0.001)*groundWater*MFModelGetArea(itemID); // kg/day
        FluxDevel = _MDparDevelop_ClWx_Rate*pop*MFModelGetArea(itemID)*1.0e-6; // kg / day  (population density given in persons/km2 )
        FluxAg    = _MDparAg_ClWx_Rate*ag/100.*MFModelGetArea(itemID); // kg /day
        // Calculate immobile zone fluxes (assume eulerian scheme)
        immFlux = immGWbeta * immGWExchange_rate *0.001* Conc_imm_last * exp(-immGWExchange_rate*dt) - immGWbeta * immGWExchange_rate * 0.001*Conc_GW_last * exp(-immGWExchange_rate * dt) ; 
        betaTerm = ((groundWater - _MDparGreaterGroundWaterMixing)*MFModelGetArea(itemID)*0.001)+immGWbeta*(1.-exp(-immGWExchange_rate * dt));
        ConcClGroundWater = 1000.0 * ( 0.001 * Conc_GW_last + (FluxPrecip2GW+FluxRecharge+FluxWeathering+FluxDevel+immFlux)/betaTerm ) / (1 + (baseflow * MFModelGetArea(itemID)*0.001*dt) / betaTerm); // Conc Cl in gw mg/L
        Conc_imm_last = 1000.* (0.001*Conc_imm_last*exp(-immGWExchange_rate * dt) + 0.001*Conc_GW_last * (1. - exp (-immGWExchange_rate * dt)) + \
                                ((0.001*ConcClGroundWater - 0.001*Conc_GW_last )/dt) * ( 1. - (1./immGWExchange_rate)*(1. - exp(-immGWExchange_rate * dt)))); // mg /L
        // Only pathway out of groundwater is baseflow
        FluxBaseFlow = 0.001*ConcClGroundWater*(0.001)*baseflow*MFModelGetArea(itemID) ; // kg/day Cl
        MassClGroundWater = MDMaximum((0.001*ConcClGroundWater*0.001*(groundWater-_MDparGreaterGroundWaterMixing)+0.001*Conc_imm_last*0.001*_MDparGreaterGroundWaterMixing)*MFModelGetArea(itemID),0.0); // kg Cl
    }
    // If Soil Moist goes to zero - then mass stays in the Root Zone.
    // Update RootZone Mass with fluxes out
    MassClRootZone = MDMaximum(MassClRootZone - FluxRecharge * dt - FluxSurfFlow * dt, 0.0); // kg Cl
    // Calculate accumulation in Surface Runoff Pool
    MassClSurfROPoolPre = MassClSurfROPool; // kg Cl
    MassClSurfROPool = MassClSurfROPool + FluxSurfFlow * dt + FluxPrecip2SW * dt;
    ConcClSurfROPool = surfaceRunoffPool > 0.0 ? (MassClSurfROPool)/ ((surfaceRunoffPool + runoffpoolrelease) * MFModelGetArea(itemID) / 1000.) * (1000.0) : 0.0 ; // mg/L Cl
    // Only pathway out of surfROpool is SurfROPoolRelease
    FluxSurfROPoolRelease = ConcClSurfROPool * (0.001) * runoffpoolrelease*MFModelGetArea(itemID) / ( 1000.) ; // kg/day Cl
    FluxSurfROPoolRelease = surfaceRunoffPool == 0.0 ? MassClSurfROPool / dt : FluxSurfROPoolRelease;
    // Update Surface RO Pool mass with flux out
    MassClSurfROPool = MDMaximum(MassClSurfROPool - FluxSurfROPoolRelease * dt,0.0); // kg Cl

    MFVarSetFloat( _MDOutClRootZoneID , itemID, MassClRootZone);
    MFVarSetFloat( _MDOutClGrdWatID , itemID, MassClGroundWater);
    MFVarSetFloat( _MDOutClSurfaceRunoffPoolID, itemID, MassClSurfROPool ); 
    MFVarSetFloat( _MDOutConcClgwPreID, itemID, ConcClGroundWater );
    MFVarSetFloat( _MDOutConcClimmPreID, itemID, Conc_imm_last ); 
    // Mass Balance sets
    MFVarSetFloat( _MDOutCldeicerInputID,itemID, FluxImp2Runoff+FluxMelt);
    MFVarSetFloat( _MDOutCltotalInputID,itemID,FluxImp2Runoff+FluxMelt+FluxPrecip2GW+FluxPrecip2SW+FluxWeathering+FluxDevel+FluxAg);
    
    /********************* WTM:  Calculate Network Ion Balance  *********************/
    discharge                  = MFVarGetFloat(_MDInDischargeID,    itemID,0.0); // m3/s
    waterStorage               = MFVarGetFloat(_MDInRiverStorageID, itemID,0.0); // m3/s
    preFlux_Cl                = MFVarGetFloat(_MDOutClFluxID,     itemID,0.0); // kg/d
    storeWater_Cl             = MFVarGetFloat(_MDOutClStoreID,    itemID,0.0); // kg/d
        
    local_load = (FluxImp2Runoff + FluxBaseFlow + FluxSurfROPoolRelease);  // kg/day Cl
        
    ClTotalIn = local_load + preFlux_Cl  + storeWater_Cl;    // kg/day Cl
    ClTotalIn = discharge <= 0.0000001 ? 0.0 : ClTotalIn;   // kg/day Cl // eliminate ionic mass flux on negiglible flow days
    postConc_Cl = discharge <= 0.0000001 ? 0.0 : (1000.) * ClTotalIn / (discharge * 86400.); // mg/L Cl 
    //Calculates the fluxes/storage for the mixing (appropriate for speccond) case
    postFlux_Cl  = (discharge * MDConst_m3PerSecTOm3PerDay) * (0.001) * postConc_Cl ;         // kg/day Cl
    postStoreWater_Cl  = (waterStorage * MDConst_m3PerSecTOm3PerDay ) * (0.001) * postConc_Cl ;      // kg/day Cl
    
    // Calculates the mass balances
    massBalance_Cl  = (ClTotalIn - (postFlux_Cl + postStoreWater_Cl ))/MDMaximum(ClTotalIn,0.0000001);
    
    // Print mass balance errors
    if (MFDateGetCurrentYear() > 0){
        if ( (massBalance_Cl > 0.001)) {
           printf("itemID = %d, %d-%d-%d, MB_SC = %f\n",itemID, MFDateGetCurrentYear(), MFDateGetCurrentMonth(), MFDateGetCurrentDay(), massBalance_Cl);
           printf("\tTotalIn: %f, Local: %f , UpFlux: %f, InStore: %f\n",ClTotalIn,local_load,preFlux_Cl,storeWater_Cl);
           printf("\tDownFlux: %f, discharge; %f, OutStore: %f , storage: %f\n",postFlux_Cl,discharge, postStoreWater_Cl,waterStorage);
        }
    }
    
    // Set Output
    MFVarSetFloat (_MDOutClFluxID,        itemID, postFlux_Cl);
    MFVarSetFloat (_MDOutClStoreID,       itemID, postStoreWater_Cl);
    MFVarSetFloat (_MDOutLocalLoadClID,   itemID, local_load);
    MFVarSetFloat (_MDOutPostConc_ClID,   itemID, postConc_Cl);
    
    float b = 31.2; // Trowbridge: -22.0 ... Novotny: -37.25
    float m = 3.96; // Trowbridge: 0.307      Novotny: 0.25
    postConc_SC = m * postConc_Cl + b;
    
    MFVarSetFloat (_MDOutPostConc_SCID,       itemID, postConc_SC); // uS/cm
}

enum {MDcalculate, MDcalculate2, MDcalculate3, MDnone};

int MDSpecCondDef () {

    float par;
    const char *InputOptStr;
    if (((InputOptStr = MFOptionGet (MDParFallThreshold))  != (char *) NULL) && (sscanf (InputOptStr,"%f",&par) == 1)) _MDSnowFallThreshold = par;
    if (((InputOptStr = MFOptionGet (MDParClWxRate)) != (char *) NULL) && (sscanf(InputOptStr,"%f",&par) == 1)) _MDparGeoChem_ClWx_Rate = par;    
    if (((InputOptStr = MFOptionGet (MDParDevelopClWxRate)) != (char *) NULL) && ( sscanf(InputOptStr,"%f",&par) == 1)) _MDparDevelop_ClWx_Rate = par;
    if (((InputOptStr = MFOptionGet (MDParDeepGroundWater)) != (char *) NULL) && (sscanf(InputOptStr,"%f",&par) == 1)) _MDparGreaterGroundWaterMixing = par;
    if (((InputOptStr = MFOptionGet (MDParDeepGroundWaterALPHA)) != (char *) NULL) &&(sscanf(InputOptStr,"%f",&par) == 1)) _MDparGreaterGroundWaterExchange = par;
    if (((InputOptStr = MFOptionGet (MDParAgricultureClRate)) != (char *) NULL) && (sscanf(InputOptStr,"%f",&par) == 1)) _MDparAg_ClWx_Rate = par;
    if (((InputOptStr = MFOptionGet (MDParCleanTransitionYear)) != (char *) NULL) && (sscanf(InputOptStr,"%f",&par) == 1)) _MDparCleanTransitionYear = par;
    if (((InputOptStr = MFOptionGet (MDParCleanTransitionConc)) != (char *) NULL) && (sscanf(InputOptStr,"%f",&par) == 1)) _MDparCleanTransitionConc = par;
    if (((InputOptStr = MFOptionGet (MDParPassiveSoilStoreFactor)) != (char *) NULL) && (sscanf(InputOptStr,"%f",&par) == 1)) _MDparPassiveSoilStoreFactor = par;
        
        InputOptStr = MFOptionGet (MDParDeicerClconc);
        CMmsgPrint(CMmsgUsrError," DeicerOption: %s checks as spatially? %d \n", InputOptStr,strcmp(InputOptStr,"spatially")==0);
        
    if (((InputOptStr  = MFOptionGet (MDParDeicerClconc)) != (char *) NULL) && (strcmp(InputOptStr,"spatially")==0)){
            CMmsgPrint(CMmsgUsrError," SpecCondDef: identified as spatially");
        if ((_MDInDeicerClID     = MFVarGetID(MDParDeicerClconcInput,       "mg/L", MFInput, MFState, MFBoundary)) == CMfailed) return (CMfailed);
    }
    CMmsgPrint(CMmsgUsrError," SpecCondDef: InputOptStr: %s InDeicerClID %d \n",InputOptStr,_MDInDeicerClID);
    if (_MDInDeicerClID == MFUnset) {
       if (((InputOptStr = MFOptionGet (MDParDeicerClconc))  != (char *) NULL) && (sscanf (InputOptStr,"%f",&par) == 1)) _MDparDeicerCl = par;
    }
   
        // Test if the percolation pathway is active (if its not - PercolationBETA should not be in the Options)
        if ((InputOptStr = MFOptionGet(MDParSoilPercolationBETA)) != (char *) NULL) {        
            if ((_MDInSoilPercolationID = MFVarGetID (MDVarSoilPercolation,     "mm", MFInput, MFFlux,   MFBoundary)) == CMfailed) return (CMfailed); // SZ 10212014
        }
    int  optID = MFUnset;											//SZ 08212014
    const char *optStr, *optName = MDOptSpecConductance;								//SZ 08212014
    const char *options [] = { MDCalculateStr, MDCalculate2Str, MDCalculate3Str, MDNoneStr, (char *) NULL };		//SZ 08212014

	MFDefEntering ("Specific Conductance Routing");
	
    if ((optStr = MFOptionGet (optName)) != (char *) NULL) optID = CMoptLookup (options, optStr, true);  //SZ 08212014

    if (
  //         ((_MDInDischargeID                  = MDDischargeDef()) == CMfailed ) ||
           ((_MDInDINFluxID                    = MDDINDef ()) == CMfailed) ||     // Needed for merging with upstream
 //           ((_MDInWTempRiverID               =   MDWTempRiverRouteDef ()) == CMfailed)  ||
          ((_MDInWTempRiverID                 = MFVarGetID (MDVarWTemp_QxT,              "degC",      MFInput, MFState, MFBoundary)) == CMfailed)   ||
 //         ((_MDInRiverWidthID                 = MDRiverWidthDef ())     == CMfailed) ||
//          ((_MDInLitterFall_POCID             = MDLitterFallDef ()) == CMfailed) ||
//          ((_MDInLocalLoad_DOCID              = MFVarGetID (MDVarLocalLoadDOC,                             "kg/d",    MFInput,  MFFlux,  MFBoundary))   == CMfailed) ||
	    ((_MDInDischargeID                  = MFVarGetID (MDVarDischarge,                                "m3/s",    MFInput,  MFState, MFBoundary))   == CMfailed) ||
            ((_MDInRiverStorageID               = MFVarGetID (MDVarRiverStorage,                           "m3/s",    MFInput,  MFState, MFInitial))    == CMfailed) ||	
	    ((_MDInRiverStorageChgID		= MFVarGetID (MDVarRiverStorageChg,			   "m3/s",    MFInput,  MFState, MFInitial))    == CMfailed) ||		
	    ((_MDInSinuosityID                   = MFVarGetID (MDVarSinuosity,                                  "m/m",    MFInput,  MFState,MFBoundary))    == CMfailed) ||
	    ((_MDInRiverbedVelocityMeanID       = MFVarGetID (MDVarRiverbedVelocityMean,                       "m/s",    MFInput,  MFState,MFBoundary))    == CMfailed) ||
	    ((_MDInBaseFlowID                   = MFVarGetID (MDVarBaseFlow,                                  "mm",    MFInput,  MFFlux,MFBoundary))    == CMfailed) ||
            ((_MDInStormRunoffTotalID           = MFVarGetID (MDVarStormRunoffTotal,                          "mm",   MFInput, MFFlux, MFBoundary))  == CMfailed) ||
            ((_MDInRunoffPoolReleaseID          = MFVarGetID (MDVarRunoffPoolRelease,                        "mm",  MFInput, MFFlux, MFBoundary))  == CMfailed) ||
            ((_MDInSurfaceRunoffID              = MFVarGetID (MDVarRainSurfRunoff,                           "mm",   MFOutput, MFFlux,  MFBoundary))  == CMfailed) || // Runoff Pool doesn't define SurfRunoff only RainSurfRunoff
            ((_MDInSurfaceRunoffPoolID          = MFVarGetID (MDVarRunoffPool,                           "mm", MFOutput, MFState, MFInitial))  == CMfailed)     ||
            ((_MDInAirTemperatureID             = MFVarGetID (MDVarAirTemperature,         "degC",       MFInput,  MFState, MFBoundary)) == CMfailed) ||
            ((_MDInSubFractionID                = MFVarGetID (MDVarLandUseSpatialSub,                         "-",    MFInput,  MFState, MFBoundary))   == CMfailed) ||  
            ((_MDInHCIAID                       = MFVarGetID (MDVarHCIA,                                    "MFNoUnit", MFInput, MFState, MFBoundary))      == CMfailed) ||
            ((_MDInPopulationDensityID          = MFVarGetID (MDVarPopulationDensity,       "p/km2",       MFInput, MFState,  MFBoundary ))       == CMfailed)  ||
            ((_MDInAgFractionID                = MFVarGetID (MDVarLandUseSpatialAg,                         "-",    MFInput,  MFState, MFBoundary))   == CMfailed) ||    
            ((_MDOutPostConc_SCID               = MFVarGetID (MDVarPostSpecCond,    	      "uS/cm",   MFOutput,  MFState, MFBoundary))   == CMfailed)               

        ) return (CMfailed);
                
        
    switch (optID) {	//SZ 08212014

        case MDcalculate: // "calculate" <- Uses simple empirical baseflow relation with development
            
            if (
                ((_MDOutLocalLoadSCID               = MFVarGetID (MDVarLocalLoadSC,               "ic/d",   MFOutput,  MFFlux,  MFBoundary))   == CMfailed) ||
                ((_MDOutStoreWater_SCID             = MFVarGetID (MDVarStoreWaterSC,              "ic/day",   MFOutput,  MFState, MFInitial))    == CMfailed) ||	
                ((_MDOutFlux_SCID                   = MFVarGetID (MDVarFluxSC,                    "ic/day",   MFRoute,   MFFlux,  MFBoundary))   == CMfailed) ||
                ((_MDOutSurfaceRunoffPool_SCID      = MFVarGetID (MDVarSurfRunoffPoolSC,          "uS/cm",  MFOutput, MFState, MFInitial))      == CMfailed) ||
                ((MFModelAddFunction (_MDSpecCond) == CMfailed))
            ) return (CMfailed);
            break;
        case MDcalculate2: // "calculate2" <- Uses impervious snowfall generating Cl input to system
            if(
                ((_MDInPrecipitationID                = MFVarGetID (MDVarPrecipitation,      "mm", MFInput, MFFlux,  MFBoundary)) == CMfailed) ||
                ((_MDInWinterHCIAID                   = MFVarGetID (MDVarWinterHCIA,        "MFNoUnit", MFInput, MFState, MFBoundary))      == CMfailed) ||
                ((_MDInImpFracSpatialID                = MFVarGetID (MDVarImpFracSpatial,          "-",    MFInput,  MFState, MFBoundary))   == CMfailed) ||    
                ((_MDInH2OFracSpatialID                = MFVarGetID (MDVarH2OFracSpatial,          "-",    MFInput,  MFState, MFBoundary))   == CMfailed) ||    
                ((_MDInAtmChlorideID                    = MFVarGetID (MDVarAtmChloride,          "mg/L",    MFInput,  MFState, MFBoundary))   == CMfailed) ||    
                ((_MDInSnowPackID                      = MFVarGetID (MDVarSnowPack,                "mm",   MFInput, MFState, MFInitial))  == CMfailed) ||
                ((_MDInSPackChgID                      = MFVarGetID (MDVarSnowPackChange,           "mm",  MFInput, MFFlux,  MFBoundary)) == CMfailed) ||
                ((_MDInSoilMoistID                     = MFVarGetID (MDVarRainSoilMoisture,           "mm",MFInput, MFState, MFBoundary)) == CMfailed) ||    
                ((_MDInAvailableWaterCapacityID        = MFVarGetID (MDVarSoilAvailWaterCapInput,         "mm",     MFInput,  MFState, MFBoundary)) == CMfailed) ||
                ((_MDInGrdWatRechargeID                = MFVarGetID (MDVarGroundWaterRecharge,      "mm", MFInput, MFFlux,  MFBoundary)) == CMfailed) ||
                ((_MDInGrdWatID                        = MFVarGetID (MDVarGroundWater,              "mm", MFInput, MFState, MFInitial))  == CMfailed) ||
                ((_MDOutClSnowPackID                   = MFVarGetID(MDVarSnowPackCl,              "kg", MFOutput, MFState, MFInitial)) == CMfailed) ||
                ((_MDOutClRootZoneID                   = MFVarGetID(MDVarRootZoneCl,              "kg", MFOutput, MFState, MFInitial)) == CMfailed) ||
                ((_MDOutClGrdWatID                     = MFVarGetID(MDVarGrdWatCl,              "kg", MFOutput, MFState, MFInitial)) == CMfailed) ||
                ((_MDOutClSurfaceRunoffPoolID          = MFVarGetID(MDVarSurfROPoolCl,              "kg", MFOutput, MFState, MFInitial)) == CMfailed) ||
                ((_MDOutClFluxID                       = MFVarGetID(MDVarFluxCl,                  "kg/d", MFRoute, MFFlux, MFBoundary)) == CMfailed) ||
                ((_MDOutClStoreID                      = MFVarGetID(MDVarStoreWaterCl,              "kg/d", MFOutput, MFState, MFInitial)) == CMfailed) ||
                ((_MDOutLocalLoadClID                  = MFVarGetID(MDVarLocalLoadCl,              "kg/d", MFOutput, MFFlux, MFBoundary)) == CMfailed) ||
                ((_MDOutPostConc_ClID                  = MFVarGetID(MDVarPostConcCl,               "mg/L", MFOutput, MFState, MFBoundary)) == CMfailed) || 
                ((_MDOutConcClgwPreID                  = MFVarGetID(MDVarConcClgwPre,              "mg/L", MFOutput, MFState, MFInitial)) == CMfailed) || 
                ((_MDOutConcClimmPreID                  = MFVarGetID(MDVarConcClimmPre,            "mg/L", MFOutput, MFState, MFInitial)) == CMfailed) || 
                ((_MDOutCldeicerInputID                 = MFVarGetID(MDVarFluxClimpInput,         "kg/d",MFOutput,MFFlux, MFBoundary)) == CMfailed)  ||
                ((_MDOutCltotalInputID                 = MFVarGetID(MDVarFluxCltotalInput,         "kg/d",MFOutput,MFFlux, MFBoundary)) == CMfailed)  ||
                ((MFModelAddFunction (_MDTotalDissIons) == CMfailed))
                    ) return (CMfailed) ;
            break;
        case MDcalculate3: // "calculate3" <- Chloride integrated with PnET - tracks deicer at surfaces - mixes other sources in Groundwater
            if(
                ((_MDInWinterHCIAID                   = MFVarGetID (MDVarWinterHCIA,        "MFNoUnit", MFInput, MFState, MFBoundary))      == CMfailed) ||
                ((_MDInImpFracSpatialID                = MFVarGetID (MDVarImpFracSpatial,          "-",    MFInput,  MFState, MFBoundary))   == CMfailed) ||    
                ((_MDInH2OFracSpatialID                = MFVarGetID (MDVarH2OFracSpatial,          "-",    MFInput,  MFState, MFBoundary))   == CMfailed) ||    
                ((_MDInAvailableWaterCapacityID        = MFVarGetID (MDVarSoilAvailWaterCapInput,         "mm",     MFInput,  MFState, MFBoundary)) == CMfailed) ||
                ((_MDInPnET_Imp_MeltID                 = MFVarGetID (MDVarImperviousSnowMelt,           "mm",   MFInput,    MFFlux,    MFBoundary)) == CMfailed) ||
                ((_MDInPnET_Imp_StormflowID            = MFVarGetID (MDVarStormRunoffImp,           "mm",   MFInput,    MFFlux,    MFBoundary)) == CMfailed) ||    
                ((_MDInGrdWatRechargeID                = MFVarGetID (MDVarGroundWaterRecharge,      "mm", MFInput, MFFlux,  MFBoundary)) == CMfailed) ||
                ((_MDInGrdWatID                        = MFVarGetID (MDVarGroundWater,              "mm", MFInput, MFState, MFInitial))  == CMfailed) ||
                ((_MDOutClRootZoneID                   = MFVarGetID(MDVarRootZoneCl,              "kg", MFOutput, MFState, MFInitial)) == CMfailed) ||
                ((_MDOutClGrdWatID                     = MFVarGetID(MDVarGrdWatCl,              "kg", MFOutput, MFState, MFInitial)) == CMfailed) ||
                ((_MDOutClSurfaceRunoffPoolID          = MFVarGetID(MDVarSurfROPoolCl,              "kg", MFOutput, MFState, MFInitial)) == CMfailed) ||
                ((_MDOutClFluxID                       = MFVarGetID(MDVarFluxCl,                  "kg/d", MFRoute, MFFlux, MFBoundary)) == CMfailed) ||
                ((_MDOutClStoreID                      = MFVarGetID(MDVarStoreWaterCl,              "kg/d", MFOutput, MFState, MFInitial)) == CMfailed) ||
                ((_MDOutLocalLoadClID                  = MFVarGetID(MDVarLocalLoadCl,              "kg/d", MFOutput, MFFlux, MFBoundary)) == CMfailed) ||
                ((_MDOutPostConc_ClID                  = MFVarGetID(MDVarPostConcCl,               "mg/L", MFOutput, MFState, MFBoundary)) == CMfailed) || 
                ((_MDOutConcClgwPreID                  = MFVarGetID(MDVarConcClgwPre,              "mg/L", MFOutput, MFState, MFInitial)) == CMfailed) || 
                ((_MDOutConcClimmPreID                  = MFVarGetID(MDVarConcClimmPre,            "mg/L", MFOutput, MFState, MFInitial)) == CMfailed) || 
                ((_MDOutCldeicerInputID                 = MFVarGetID(MDVarFluxClimpInput,         "kg/d",MFOutput,MFFlux, MFBoundary)) == CMfailed)  ||
                ((_MDOutCltotalInputID                 = MFVarGetID(MDVarFluxCltotalInput,         "kg/d",MFOutput,MFFlux, MFBoundary)) == CMfailed)  ||
                ((MFModelAddFunction (_MDChloride_PnET) == CMfailed))
                    ) return (CMfailed) ;
            break;
        default: 
            MFOptionMessage (optName, optStr, options); return (CMfailed);
        }
	MFDefLeaving ("Specific Conductance Routing");
    return (_MDOutFlux_SCID);
}
