/******************************************************************************

GHAAS Water Balance/Transport Model V3.0
Global Hydrologic Archive and Analysis System
Copyright 1994-2007, University of New Hampshire

MDNitrogenInputs.c  - Nitrogen Loading

rob.stewart@unh.edu

*******************************************************************************/

#include <cm.h>
#include <MF.h>
#include <MD.h>
#include <math.h>

// Input

static int _MDInLandUseID		= MFUnset;
static int _MDInRunoffID		= MFUnset;
static int _MDInRunoffVolID		= MFUnset;
static int _MDInRiverOrderID		= MFUnset;		//RJS 090508
static int _MDInDINLoadConcID		= MFUnset;		//RJS 102810
static int _MDInRunoffPoolMassRelID     = MFUnset;		//RJS 050511
static int _MDInGrdWatMassRelID		= MFUnset;		//RJS 050511
static int _MDInLawnFractionID          = MFUnset;		//RJS 051111
static int _MDInLoadAdjustID		= MFUnset;		//RJS 112211
static int _MDInTotalPointID		= MFUnset;
static int _MDInLandUseSubID		= MFUnset;
static int _MDInLandUseAgID		= MFUnset;
static int _MDInLocalLoad_DINID         = MFUnset;              //RJS 011414
static int _MDInLocalLoad_Con_DINID     = MFUnset;              //RJS 081814
static int _MDInLocalLoad_Dec_DINID     = MFUnset;              //RJS 081814
static int _MDInLocalLoad_Mix_DINID     = MFUnset;              //RJS 081814
static int _MDInLocalLoad_Lak_DINID     = MFUnset;              //RJS 081914
static int _MDInLocalLoad_Wet_DINID     = MFUnset;
static int _MDInLocalLoad_Imp_DINID     = MFUnset;              // SZ 20150421
static int _MDInLocalLoad_Lawn_DINID    = MFUnset;              // SZ 20150421
static int _MDInLocalLoad_Ag_DINID      = MFUnset;              // SZ 20150421
static int _MDInLandUseConID            = MFUnset;              //RJS 081914
static int _MDInLandUseDecID            = MFUnset;              //RJS 081914
static int _MDInLandUseMixID            = MFUnset;              //RJS 081914
static int _MDInH2OFractionID           = MFUnset;              //RJS 081914
static int _MDInAsymID                  = MFUnset;              //RJS 102314
static int _MDInScaleID                 = MFUnset;              //RJS 102314
static int _MDInXmid_bID                = MFUnset;              //RJS 102314
static int _MDInXmid_mID                = MFUnset;              //RJS 102314
static int _MDInImpFractionID           = MFUnset;
static int _MDInWetlandsID              = MFUnset;
static int _MDInAgMultiplierID		= MFUnset;		//RJS 051815

static int _MDInPropROStormWaterID	= MFUnset; 		//SZ 20160211
static int _MDInPropROSurfaceWaterID	= MFUnset; 		//SZ 20160211
static int _MDInPropROGroundWaterID	= MFUnset; 		//SZ 20160211
// Output

static int _MDOutLocalLoad_DINID                   = MFUnset;		// RJS 090308
static int _MDOutLocalConc_DINID                   = MFUnset;		// RJS 090308
static int _MDOutLocalLoadnew_DINID                = MFUnset;		// RJS 061511	just to CHECK whether the new and old loading are the same
static int _MDOutDINLoadConcID                     = MFUnset;		// KYLE
static int _MDOutLocalLoad_Agr_DINID		   = MFUnset;
static int _MDOutDINSubLoadConcID		   = MFUnset;
static int _MDOutDINAgLoadConcID		   = MFUnset;
static int _MDOutLocalLoad_PnET_DINID              = MFUnset;
static int _MDOutLocalLoad_Lak_DINID               = MFUnset;
static int _MDOutLocalLoad_Wet_DINID               = MFUnset;
static int _MDOutLocalLoad_Con_DINID               = MFUnset;
static int _MDOutLocalLoad_Dec_DINID               = MFUnset;
static int _MDOutLocalLoad_Mix_DINID               = MFUnset;
static int _MDOutLocalLoad_Sub_DINID               = MFUnset;

static int _MDOutLocalConc_Sub_DINID               = MFUnset;
static int _MDOutLocalConc_Agr_DINID               = MFUnset;
static int _MDOutLocalConc_Lak_DINID               = MFUnset;
static int _MDOutLocalConc_Wet_DINID               = MFUnset;
static int _MDOutLocalConc_Con_DINID               = MFUnset;
static int _MDOutLocalConc_Dec_DINID               = MFUnset;
static int _MDOutLocalConc_Mix_DINID               = MFUnset;
static int _MDOutLocalConc_PnET_DINID              = MFUnset;
static int _MDOutTotalRiparianRemoval_DINID	   = MFUnset; // SZ 20160211

static void _MDNitrogenInputsInput (int itemID) {

	//Input
	float runoff;		//mm/day
	float runoffVol; 	//m3/sec
	float LocalConc_DIN;
	float riverOrder;
        float LocalLoad_DIN_in; // kg/d

	//Output
	float LocalLoad_DIN;

//	runoff             = MFVarGetFloat (_MDInRunoffID,         	  itemID, 0.0); // mm / d
//	runoffVol          = MFVarGetFloat (_MDInRunoffVolID,             itemID, 0.0); // m3/sec
//	LocalConc_DIN	   = MFVarGetFloat (_MDInDINLoadConcID,           itemID, 0.0); // mg/L
        LocalLoad_DIN_in   = MFVarGetFloat (_MDInLocalLoad_DINID,         itemID, 0.0); // kg/d
//	riverOrder	   = MFVarGetFloat (_MDInRiverOrderID,            itemID, 0.0);

//	LocalLoad_DIN      = runoffVol * 86400 * LocalConc_DIN / 1000; // kg/day
        LocalLoad_DIN      = LocalLoad_DIN_in;                         // kg/day
        
        
	MFVarSetFloat (_MDOutLocalLoad_DINID,        itemID, LocalLoad_DIN);	  // RJS 090308

}

static void _MDNitrogenInputsPnET (int itemID) {

	//Input
	float runoff                    = 0.0; //mm/day
	float runoffVol                 = 0.0; //m3/sec
	float LocalConc_DIN             = 0.0;
	float riverOrder                = 0.0;
   
        float luLawn                    = 0.0; // proportion lawn
        float luImp                     = 0.0; // proportion impervious
        float luSub                     = 0.0; // PERCENT suburban or developed
        float luAgr                     = 0.0; // PERCENT agriculture
        float luCon                     = 0.0; // proportion suburban or developed
        float luDec                     = 0.0; // proportion agriculture
        float luMix                     = 0.0; // proportion suburban or developed
        float luWet                     = 0.0; // proportion agriculture
        float luLak                     = 0.0; // proportion lake or water
        
        float LocalLoad_Con_DIN_in      = 0.0; // kg/d
        float LocalLoad_Dec_DIN_in      = 0.0; // kg/d
        float LocalLoad_Mix_DIN_in      = 0.0; // kg/d
        float LocalLoad_Lak_DIN_in      = 0.0; // kg/d
        float LocalLoad_Wet_DIN_in      = 0.0; // kg/d
        float LocalLoad_Imp_DIN_in      = 0.0; // kg/d
        float LocalLoad_Lawn_DIN_in     = 0.0; // kg/d
        float LocalLoad_Ag_DIN_in       = 0.0; // kg/d

        float LocalLoad_Con_DIN         = 0.0; // kg/d
        float LocalLoad_Dec_DIN         = 0.0; // kg/d 
        float LocalLoad_Mix_DIN         = 0.0; // kg/d
        float LocalLoad_Sub_DIN         = 0.0; // kg/d
        float LocalLoad_Agr_DIN         = 0.0; // kg/d
        float LocalLoad_Lak_DIN         = 0.0; // kg/d
        float LocalLoad_Wet_DIN         = 0.0; // kg/d
        float LocalLoad_Imp_DIN         = 0.0; // kg/d
        float LocalLoad_Lawn_DIN        = 0.0; // kg/d
        float LocalLoad_Ag_DIN          = 0.0; // kg/d
        
        float LocalLoad_PnET_DIN        = 0.0; // kg/d
        float LocalConc_PnET_DIN        = 0.0; // mg/L
        
        float LocalConc_Con_DIN         = 0.0; // mg/L
        float LocalConc_Dec_DIN         = 0.0; // mg/L
        float LocalConc_Mix_DIN         = 0.0; // mg/L
        float LocalConc_Sub_DIN         = 0.0; // mg/L
        float LocalConc_Agr_DIN         = 0.0; // mg/L
        float LocalConc_Lak_DIN         = 0.0; // mg/L
        float LocalConc_Wet_DIN         = 0.0; // mg/L
        
        float scale                     = 0.0;
        float asym                      = 0.0;
        float xMid                      = 0.0;
        float area_m2                   = 0.0;
	float multiplier		= 0.0;	// multiplier for Agr DIN concentration RJS 051815
        
	//Output
	float LocalLoad_DIN;
        
        if (_MDInLawnFractionID != MFUnset)	  luLawn = MFVarGetFloat (_MDInLawnFractionID,  itemID, 0.0);	
	if (_MDInImpFractionID  != MFUnset)	   luImp = MFVarGetFloat (_MDInImpFractionID,   itemID, 0.0);	
	if (_MDInLandUseAgID   != MFUnset)	   luAgr = MFVarGetFloat (_MDInLandUseAgID,     itemID, 0.0);	// PERCENT
        if (_MDInLandUseConID  != MFUnset)	   luCon = MFVarGetFloat (_MDInLandUseConID,    itemID, 0.0);	
	if (_MDInLandUseDecID  != MFUnset)	   luDec = MFVarGetFloat (_MDInLandUseDecID,    itemID, 0.0);	
	if (_MDInLandUseMixID  != MFUnset)	   luMix = MFVarGetFloat (_MDInLandUseMixID,    itemID, 0.0);	
	if (_MDInH2OFractionID != MFUnset)	   luLak = MFVarGetFloat (_MDInH2OFractionID,   itemID, 0.0);	
	if (_MDInWetlandsID != MFUnset)            luWet = MFVarGetFloat (_MDInWetlandsID,   itemID, 0.0);	

        luSub                  = (luLawn + luImp) * 100;
        area_m2                = MFModelGetArea (itemID);
        LocalLoad_Con_DIN_in   = MFVarGetFloat (_MDInLocalLoad_Con_DINID,         itemID, 0.0) * area_m2 * luCon / 1000; // converting gN/m2/d to kg/d for Con area only 
        LocalLoad_Dec_DIN_in   = MFVarGetFloat (_MDInLocalLoad_Dec_DINID,         itemID, 0.0) * area_m2 * luDec / 1000; // converting gN/m2/d to kg/d for Dec area only
        LocalLoad_Mix_DIN_in   = MFVarGetFloat (_MDInLocalLoad_Mix_DINID,         itemID, 0.0) * area_m2 * luMix / 1000; // converting gN/m2/d to kg/d for Mix area only
        LocalLoad_Lak_DIN_in   = MFVarGetFloat (_MDInLocalLoad_Lak_DINID,         itemID, 0.0) * area_m2 * luLak / 1000; // converting gN/m2/d to kg/d for Lak area only
        LocalLoad_Wet_DIN_in   = MFVarGetFloat (_MDInLocalLoad_Wet_DINID,         itemID, 0.0) * area_m2 * luWet / 1000; // converting gN/m2/d to kg/d for Wet area only
        LocalLoad_Imp_DIN_in    = MFVarGetFloat(_MDInLocalLoad_Imp_DINID,           itemID,0.0) * area_m2 * luImp / 1000; //converting gN/m2/d to kg/d for Imp area only
        LocalLoad_Lawn_DIN_in = MFVarGetFloat(_MDInLocalLoad_Lawn_DINID,        itemID,0.0) * area_m2 * luLawn / 1000; //converting gN/m2/d to kg/d for Lawn area only
        LocalLoad_Ag_DIN_in = MFVarGetFloat(_MDInLocalLoad_Ag_DINID,        itemID,0.0) * area_m2 * ( luAgr /100.) / 1000. ; //converting gN/m2/d to kg/d for Ag area only (Ag in percent)
       
        LocalLoad_Con_DIN_in   = luCon > 0.0 ? LocalLoad_Con_DIN_in : 0.0;
        LocalLoad_Dec_DIN_in   = luDec > 0.0 ? LocalLoad_Dec_DIN_in : 0.0;
        LocalLoad_Mix_DIN_in   = luMix > 0.0 ? LocalLoad_Mix_DIN_in : 0.0;
        LocalLoad_Lak_DIN_in   = luLak > 0.0 ? LocalLoad_Lak_DIN_in : 0.0;
        LocalLoad_Wet_DIN_in   = luWet > 0.0 ? LocalLoad_Wet_DIN_in : 0.0;
        LocalLoad_Imp_DIN_in    = luImp > 0.0 ? LocalLoad_Imp_DIN_in : 0.0;
        LocalLoad_Lawn_DIN_in       = luLawn > 0.0 ? LocalLoad_Lawn_DIN_in : 0.0;
        LocalLoad_Ag_DIN_in         = luAgr > 0.0 ? LocalLoad_Ag_DIN_in : 0.0;
      
        luSub = ((luSub/100) + (luAgr/100) + luCon + luDec + luMix + luLak + luWet) > 1.0 ? (luSub/100) / ((luSub/100) + (luAgr/100) + luCon + luDec + luMix + luLak + luWet) * 100.0 : luSub;           
        luAgr = ((luSub/100) + (luAgr/100) + luCon + luDec + luMix + luLak + luWet) > 1.0 ? (luAgr/100) / ((luSub/100) + (luAgr/100) + luCon + luDec + luMix + luLak + luWet) * 100.0 : luAgr;   
                
        runoff             = MFVarGetFloat (_MDInRunoffID,         itemID, 0.0); 	// mm / d
	runoffVol          = MFVarGetFloat (_MDInRunoffVolID,      itemID, 0.0); 	// m3/sec
	multiplier         = MFVarGetFloat (_MDInAgMultiplierID,   itemID, 0.0); 	// unitless
	asym  		   = MFVarGetFloat (_MDInAsymID,           itemID, 0.0);	// RJS 051815

	/**************  Calculate Riparian Removal ********************
	 *
	 * Assume that runoff from each forest type from FrAMES exits
	 * from the three runoff flow-paths according to the proportion 
	 *  of flow-paths on each day.  
	 *  The Storm Runoff flow-path receives no removal.
	 *  The Surface Runoff Pool release flow-path receives removal with Da=0.511
	 *  The Baseflow Runoff Pool release flow-path receives removal with Da=0.223
	 *  These values are estimated from Mark Green's synthesis.
	 *
	 * ***********************************************************/        

	float f_SurfflowTmp = 0.0; // Fraction of temp surface
	float f_Surfflow = 0.0; // Fraction of surface pool release
	float f_Baseflow = 0.0; // Fraction of stormflow
	float LocalLoad_Con_Storm_DIN = 0.0; // Load from PnET Coniferous through stormflow [kg/d]
	float LocalLoad_Con_Surf_DIN = 0.0; // Load from PnET Coniferous through surface ro pool release flow [kg/d]
	float LocalLoad_Con_Base_DIN = 0.0; // Load from PnET Coniferous through baseflow [kg/d]
	float LocalLoad_Dec_Storm_DIN = 0.0; // Load from PnET Deciduous through stormflow [kg/d]
	float LocalLoad_Dec_Surf_DIN = 0.0; // Load from PnET Deciduous shrough surface ro pool release flow [kg/d]
	float LocalLoad_Dec_Base_DIN = 0.0; // Load from PnET Deciduous through baseflow [kg/d]
	float LocalLoad_Mix_Storm_DIN = 0.0; // Load from PnET Mixed through stormflow [kg/d]
	float LocalLoad_Mix_Surf_DIN = 0.0; // Load from PnET Mixed through surface ro pool release flow [kg/d]
	float LocalLoad_Mix_Base_DIN = 0.0; // Load from PnET Mixed through baseflow [kg/d]
	float LocalLoad_Wet_Surf_DIN = 0.0; // Load from PnET Wetland through surface ro pool release flow [kg/d]
	float LocalLoad_Wet_Base_DIN = 0.0; // Load from PnET Wetland through baseflow [kg/d]


	float Rem_Surf = 0.70; // Damkohler = ? 0.4; //Damkohler = 0.511 for SROpool flow-path
	float Rem_Base = 0.70; // Damkohler = ? 0.2; //Damkohler =0.223 for SROpool flow-path
	float Rem_H2O = 0.00; // ...

	float Removal_Con_Surf = 0.0; // Removal from Coniferous load through Surface RO pool flow-path [kg/d]
	float Removal_Con_Base = 0.0; // Removal from Coniferous load through Baseflow pool flow-path [kg/d]
	float Removal_Dec_Surf = 0.0; // Removal from Deciduous load through Surface RO pool flow-path [kg/d]
	float Removal_Dec_Base = 0.0; // Removal from Deciduous load through Baseflow pool flow-path [kg/d]
	float Removal_Mix_Surf = 0.0; // Removal from Mixed load through Surface RO pool flow-path [kg/d]
	float Removal_Mix_Base = 0.0; // Removal from Mixed load through Baseflow pool flow-path [kg/d]
	float Removal_Wet_Surf = 0.0; // Removal from Wetland load through Surface RO pool flow-path [kg/d]
	float Removal_Wet_Base = 0.0; // Removal from Wetland load through Baseflow pool flow-path [kg/d]
	float Removal_H2O = 0.0; // Removal from Water load through in-lake removal (for small lakes) [kg/d]

	float TotalRiparianRemoval = 0.0; // Total Removal from all Riparian flow-paths
	float SurfaceRiparianRemoval = 0.0; // Riparian removal from Surface runoff pool flow-paths
	float BaseflowRiparianRemoval = 0.0; // Riparian removal from Baseflow pool flow-paths
	float ConiferousRiparianRemoval = 0.0; // Riparian removal from coniferous forests
	float DeciduousRiparianRemoval = 0.0; // Riparian removal from deciduous forests
	float MixedRiparianRemoval = 0.0; // Riparian removal from mixed forests
	float WetRiparianRemoval = 0.0; // Riparian removal from wetalnds
	float H2ORemoval = 0.0; // "Riparian" removal from lakes


    f_SurfflowTmp = MFVarGetFloat( _MDInPropROSurfaceWaterID, itemID,0.0);
	f_Baseflow = MFVarGetFloat( _MDInPropROGroundWaterID, itemID,0.0); 
    f_Surfflow =  f_SurfflowTmp + f_Baseflow > 0.0 ? f_SurfflowTmp / (f_SurfflowTmp + f_Baseflow) : 0.0;
    f_Baseflow =  f_SurfflowTmp + f_Baseflow > 0.0 ? f_Baseflow / (f_SurfflowTmp + f_Baseflow) :  0.0;

	// Calculate proportional fluxes from each flow-path and forest type
	LocalLoad_Con_Surf_DIN = LocalLoad_Con_DIN_in * f_Surfflow; // 
	LocalLoad_Con_Base_DIN = LocalLoad_Con_DIN_in * f_Baseflow; // 
	LocalLoad_Dec_Surf_DIN = LocalLoad_Dec_DIN_in * f_Surfflow; // 
	LocalLoad_Dec_Base_DIN = LocalLoad_Dec_DIN_in * f_Baseflow; //
	LocalLoad_Mix_Surf_DIN = LocalLoad_Mix_DIN_in * f_Surfflow; // 
	LocalLoad_Mix_Base_DIN =  LocalLoad_Mix_DIN_in * f_Baseflow; //
	LocalLoad_Wet_Surf_DIN = LocalLoad_Wet_DIN_in * f_Surfflow; //
	LocalLoad_Wet_Base_DIN = LocalLoad_Wet_DIN_in * f_Baseflow; //
	// Calculate removal
	Removal_Con_Surf = Rem_Surf * LocalLoad_Con_Surf_DIN; //
	Removal_Con_Base = Rem_Base * LocalLoad_Con_Base_DIN; //
	Removal_Dec_Surf = Rem_Surf * LocalLoad_Dec_Surf_DIN; //
	Removal_Dec_Base = Rem_Base * LocalLoad_Dec_Base_DIN; //
	Removal_Mix_Surf = Rem_Surf * LocalLoad_Mix_Surf_DIN; //
	Removal_Mix_Base = Rem_Base * LocalLoad_Mix_Base_DIN; //
	Removal_Wet_Surf = Rem_Surf * LocalLoad_Wet_Surf_DIN; //
	Removal_Wet_Base = Rem_Base * LocalLoad_Wet_Base_DIN; //
	Removal_H2O = LocalLoad_Lak_DIN_in * Rem_H2O; //
	// Summarize removal by type
	ConiferousRiparianRemoval = Removal_Con_Surf + Removal_Con_Base; //
	DeciduousRiparianRemoval = Removal_Dec_Surf + Removal_Dec_Base; //
	MixedRiparianRemoval = Removal_Mix_Surf + Removal_Mix_Base; //
	WetRiparianRemoval = Removal_Wet_Surf + Removal_Wet_Base; //
	H2ORemoval = Removal_H2O; //
	// Summarize removal by flow-path
	SurfaceRiparianRemoval = Removal_Con_Surf + Removal_Dec_Surf + Removal_Mix_Surf + Removal_Wet_Surf; //
	BaseflowRiparianRemoval = Removal_Con_Base + Removal_Dec_Base + Removal_Mix_Base + Removal_Wet_Base; //
	// Summarize all riparian removal
	TotalRiparianRemoval = SurfaceRiparianRemoval + BaseflowRiparianRemoval + H2ORemoval; //

	// Update LocalLoads for the three forest land-covers
	LocalLoad_Con_DIN_in -= ConiferousRiparianRemoval; //
	LocalLoad_Dec_DIN_in -= DeciduousRiparianRemoval; //
	LocalLoad_Mix_DIN_in -= MixedRiparianRemoval; //
	LocalLoad_Wet_DIN_in -= WetRiparianRemoval; //
	LocalLoad_Lak_DIN_in -= H2ORemoval; //


        if (runoff > 0.00001) {    
        LocalConc_Con_DIN  = luCon > 0.0 ? LocalLoad_Con_DIN_in / (runoffVol * luCon * 86400) * 1000 : 0.0;
        LocalConc_Dec_DIN  = luDec > 0.0 ? LocalLoad_Dec_DIN_in / (runoffVol * luDec * 86400) * 1000 : 0.0;
        LocalConc_Mix_DIN  = luMix > 0.0 ? LocalLoad_Mix_DIN_in / (runoffVol * luMix * 86400) * 1000 : 0.0;
        LocalConc_Lak_DIN  = luLak > 0.0 ? LocalLoad_Lak_DIN_in / (runoffVol * luLak * 86400) * 1000 : 0.0;
        LocalConc_Wet_DIN  = luWet > 0.0 ? LocalLoad_Wet_DIN_in / (runoffVol * luWet * 86400) * 1000 : 0.0;
        }
        
      
	scale          = 12.5;

	if (runoff <= 0.00001) {	
		LocalConc_Sub_DIN = 0.0;		
		LocalConc_Agr_DIN = 0.0;	
	}

	else {
            xMid = 51.388 + 8.4511 * log(runoff); 
            LocalConc_Sub_DIN  = (asym / (1 + pow(2.718281828, (xMid - luSub) / scale))) - (asym / (1 + pow(2.718281828, (xMid - 0.0) / scale)));                   // subtracting intercept from loading function
            LocalConc_Agr_DIN  = ((asym * multiplier) / (1 + pow(2.718281828, (xMid - luAgr) / scale))) - ((asym * multiplier) / (1 + pow(2.718281828, (xMid - 0.0) / scale)));   // subtracting intercept from loading function
        }

        /////// New PnET-FrAMES ////////
        
        LocalLoad_Sub_DIN = (luSub / 100) * runoffVol * 86400 * LocalConc_Sub_DIN / 1000;   // kg/d
        LocalLoad_Agr_DIN = (luAgr / 100) * runoffVol * 86400 * LocalConc_Agr_DIN / 1000;   // kg/d
        LocalLoad_Con_DIN = luCon * runoffVol * 86400 * LocalConc_Con_DIN / 1000;           // kg/d
        LocalLoad_Dec_DIN = luDec * runoffVol * 86400 * LocalConc_Dec_DIN / 1000;           // kg/d
        LocalLoad_Mix_DIN = luMix * runoffVol * 86400 * LocalConc_Mix_DIN / 1000;           // kg/d
        LocalLoad_Lak_DIN = luLak * runoffVol * 86400 * LocalConc_Lak_DIN / 1000;           // kg/d
        LocalLoad_Wet_DIN = luWet * runoffVol * 86400 * LocalConc_Wet_DIN / 1000;           // kg/d
     
        LocalLoad_DIN = LocalLoad_Sub_DIN + LocalLoad_Agr_DIN + LocalLoad_Con_DIN + LocalLoad_Dec_DIN + LocalLoad_Mix_DIN + LocalLoad_Lak_DIN + LocalLoad_Wet_DIN;
        LocalConc_DIN = runoffVol > 0.0 ? LocalLoad_DIN / (runoffVol * 86400) * 1000 : 0.0;
        
        LocalLoad_PnET_DIN = LocalLoad_Con_DIN_in + LocalLoad_Dec_DIN_in + LocalLoad_Mix_DIN_in + LocalLoad_Wet_DIN_in + LocalLoad_Lak_DIN_in;
 //     LocalConc_PnET_DIN = runoffVol > 0.00001 ? LocalLoad_PnET_DIN / (runoffVol * 86400) * 1000 : 0.0;
        //if ((itemID % 100) == 0 && (MFDateGetCurrentDay() % 5) == 0) printf("Cell: %d  PnETDINload: %f  PnETDINconc: %f\n",itemID,LocalLoad_PnET_DIN,LocalConc_PnET_DIN);
        /////// Old Option 1 /////////
        
//	LocalLoad_Sub_DIN  = (luSub / 100) * runoffVol * 86400 * LocalConc_Sub_DIN / 1000; 		// kg/day //KAW 2013 05 08 Loading from suburban part of each grid cell
//	LocalLoad_Agr_DIN  = (luAgr / 100) * runoffVol * 86400 * LocalConc_Agr_DIN / 1000; 		// kg/day //KAW 2013 05 08 Loading from the agricultural part of each grid cell
//      LocalLoad_Lak_DIN  = (luLak / 100) * runoffVol * 86400 * LocalConc_Lak_DIN / 1000;              // kg/day
//      LocalLoad_For_DIN  = (luForest / 100) * runoffVol * 86400 * LocalConc_For_DIN / 1000;           // kg/day
	
//      LocalLoad_DIN      = LocalLoad_Agr_DIN + LocalLoad_Sub_DIN + LocalLoad_Lak_DIN + LocalLoad_For_DIN;                // kg/day //KAW 2013 05 08: Sum of all N loadings  Total point is from MMM's Code for point source N inputs
//      LocalConc_DIN      = runoffVol > 0.0 ? LocalLoad_DIN / (runoffVol * 86400) * 1000 : 0.0;

        /////// Old Option 2 /////////
    
//        LocalLoad_Sub_DIN = runoffVol * 86400 * LocalConc_Sub_DIN2 / 1000;                         // kg/day
//        LocalLoad_Agr_DIN = runoffVol * 86400 * LocalConc_Agr_DIN2 / 1000;                         // kg/day
//        LocalLoad_For_DIN = runoffVol * (luForest / 100) * 86400 * LocalConc_For_DIN / 1000;       // kg/day
//        LocalLoad_Lak_DIN = runoffVol * (luLak / 100) * 86400 * LocalConc_Lak_DIN / 1000;          // kg/day     
        
//        LocalLoad_DIN      = LocalLoad_Agr_DIN + LocalLoad_Sub_DIN + LocalLoad_For_DIN + LocalLoad_Lak_DIN;     // kg/day 
//        LocalConc_DIN      = runoffVol > 0.0 ? LocalLoad_DIN / (runoffVol * 86400) * 1000 : 0.0;

        
//     if ((itemID == 4456) | (itemID ==1)|(itemID==336)) printf("ID: %d: area_m2 = %f,  luSub = %f, luAgr = %f, luMix = %f, luLak = %f, luCon = %f, luDec = %f, Load=%f, Conc=%f RipRemov=%f\n", itemID,area_m2, luSub, luAgr, luMix, luLak, luCon, luDec, LocalLoad_DIN, LocalConc_DIN,TotalRiparianRemoval);
   //     if (itemID == 7734) printf("runoffVol = %f, localConc_Lak_DIN = %f, LocalLoad_Lak_DIN = %f, LocalLoad_Lak_DIN_in = %f\n",runoffVol, LocalConc_Lak_DIN, LocalLoad_Lak_DIN, LocalLoad_Lak_DIN_in);
   //     if (itemID == 7734) printf("LocalLoad_Agr_DIN = %f, LocalLoad_Sub_DIN = %f, LocalLoad_Lak_DIN = %f, LocalLoad_For_DIN = %f\n", LocalLoad_Agr_DIN, LocalLoad_Sub_DIN, LocalLoad_Lak_DIN, LocalLoad_For_DIN);
        
	MFVarSetFloat (_MDOutLocalLoad_DINID,         itemID, LocalLoad_DIN);                   // RJS 090308
	MFVarSetFloat (_MDOutLocalConc_DINID,         itemID, LocalConc_DIN);                   // RJS 090308
        MFVarSetFloat (_MDOutLocalLoad_PnET_DINID,    itemID, LocalLoad_PnET_DIN);              // SZ 20150421
//        MFVarSetFloat (_MDOutLocalConc_PnET_DINID,    itemID, LocalConc_PnET_DIN);              // 20150421
        MFVarSetFloat (_MDOutLocalLoad_Sub_DINID,     itemID, LocalLoad_Sub_DIN);	  	// RJS 090308
	MFVarSetFloat (_MDOutLocalLoad_Agr_DINID,     itemID, LocalLoad_Agr_DIN);	  	// RJS 090308
	MFVarSetFloat (_MDOutLocalLoad_Wet_DINID,     itemID, LocalLoad_Wet_DIN);	  	// RJS 090308
	MFVarSetFloat (_MDOutLocalLoad_Lak_DINID,     itemID, LocalLoad_Lak_DIN);	  	// RJS 090308
	MFVarSetFloat (_MDOutLocalLoad_Con_DINID,     itemID, LocalLoad_Con_DIN);	  	// RJS 090308
	MFVarSetFloat (_MDOutLocalLoad_Mix_DINID,     itemID, LocalLoad_Mix_DIN);	  	// RJS 090308
	MFVarSetFloat (_MDOutLocalLoad_Dec_DINID,     itemID, LocalLoad_Dec_DIN);	  	// RJS 090308
        MFVarSetFloat (_MDOutLocalConc_Sub_DINID,     itemID, LocalConc_Sub_DIN);	  	// RJS 090308
	MFVarSetFloat (_MDOutLocalConc_Agr_DINID,     itemID, LocalConc_Agr_DIN);	  	// RJS 090308
	MFVarSetFloat (_MDOutLocalConc_Wet_DINID,     itemID, LocalConc_Wet_DIN);	  	// RJS 090308
	MFVarSetFloat (_MDOutLocalConc_Lak_DINID,     itemID, LocalConc_Lak_DIN);	  	// RJS 090308
	MFVarSetFloat (_MDOutLocalConc_Con_DINID,     itemID, LocalConc_Con_DIN);	  	// RJS 090308
	MFVarSetFloat (_MDOutLocalConc_Mix_DINID,     itemID, LocalConc_Mix_DIN);	  	// RJS 090308
	MFVarSetFloat (_MDOutLocalConc_Dec_DINID,     itemID, LocalConc_Dec_DIN);	  	// RJS 090308
	// Riparian Removal Services
	MFVarSetFloat(_MDOutTotalRiparianRemoval_DINID,itemID,TotalRiparianRemoval); 		//SZ 20160211
	// Other riparian removal services can be defined and output as needed.

}

static void _MDNitrogenInputsPartitioned (int itemID) {

	//Input
	float runoffPoolMassRel = 0.0;	//kg/day	RJS 050511
	float grdWatMassRel		= 0.0;	//kg/day	RJS 050511
	float runoffVol			= 0.0;  //m3/sec	RJS 050511	need this to call runoffVol
	float temp				= 0.0;	//kg/day	RJS 061611

	//Output
	float LocalLoad_DIN;

	runoffVol		   = MFVarGetFloat (_MDInRunoffVolID,         itemID, 0.0);	// m3/s  RJS 050511
	runoffPoolMassRel  = MFVarGetFloat (_MDInRunoffPoolMassRelID, itemID, 0.0); // kg/d  RJS 050511
	grdWatMassRel      = MFVarGetFloat (_MDInGrdWatMassRelID,     itemID, 0.0); // kg/d	 RJS 050511

	LocalLoad_DIN      = runoffPoolMassRel + grdWatMassRel; // kg/day	RJS 050511
	temp               = LocalLoad_DIN;						// kg/day	RJS 061611

	MFVarSetFloat (_MDOutLocalLoad_DINID,    itemID, LocalLoad_DIN);	 // RJS 050511
	MFVarSetFloat (_MDOutLocalLoadnew_DINID, itemID, temp);				 // RJS 061611
}



static void _MDNitrogenInputsCalc (int itemID) {

	// Input

	float luSub = 0.0;
	float luAgr = 0.0;
	float runoff;			// mm/day
	float runoffVol;		// m3/sec

	float scale;
	float asym;
	float xMid;

	float LocalConc_DIN      = 0.0;
	float LocalLoad_DIN      = 0.0;

	float LocalConc_Sub_DIN  = 0.0;		//KAW 2013 05 08 Suburban Loading Concentration
	float LocalLoad_Sub_DIN  = 0.0;		//KAW 2013 05 08 Suburban Loading

	float LocalConc_Agr_DIN   = 0.0;		//KAW 2013 05 08 Agricultural Loading Concentration
	float LocalLoad_Agr_DIN   = 0.0;		//KAW 2013 05 08 Agricultural Loading

	float riverOrder;
	float loadAdjust = 0.0;			// RJS 112211
	float Total_point = 0.0;                // test
	float multiplier = 0.0;			// RJS 051815

// test

	if (_MDInRiverOrderID != MFUnset)  riverOrder = MFVarGetFloat (_MDInRiverOrderID, itemID, 0.0);	// RJS 090508
//	if (_MDInLandUseID    != MFUnset)	luSub = MFVarGetFloat (_MDInLandUseID,    itemID, 0.0);
	if (_MDInLandUseSubID != MFUnset)	luSub = MFVarGetFloat (_MDInLandUseSubID, itemID, 0.0);	//KAW 2013 05 08 Suburban Land Use
	if (_MDInLandUseAgID  != MFUnset)	luAgr = MFVarGetFloat (_MDInLandUseAgID,  itemID, 0.0);	//KAW 2013 05 08 Agricultural Land Use

//	luSub = luSub + luAg > 100 ? 100 : luSub + luAg;                                // commented out 110813
	luSub = luSub + luAgr > 100.0 ? luSub / (luSub + luAgr) * 100.0 : luSub;        // RJS 110813    
        luAgr = luSub + luAgr > 100.0 ? luAgr / (luSub + luAgr) * 100.0 : luAgr;        // RJS 110813

        
	runoff             = MFVarGetFloat (_MDInRunoffID,         itemID, 0.0); 	// mm / d
	runoffVol          = MFVarGetFloat (_MDInRunoffVolID,      itemID, 0.0); 	// m3/sec
	loadAdjust	   = MFVarGetFloat (_MDInLoadAdjustID,     itemID, 1.0);	// RJS 112211, adjusts loads, keep at 1 if nodata
	multiplier	   = MFVarGetFloat (_MDInAgMultiplierID,   itemID, 1.0);	// RJS 051815
	asym		   = MFVarGetFloat (_MDInAsymID,   	   itemID, 1.0);	// RJS 051815
		

	scale          = 12.5;      // 12.2
//	asym           = 1.4;

	if (runoff < 0.00001) {
//		LocalConc_DIN = 0.0;                    // RJS commented out 110813
		LocalConc_Sub_DIN = 0.0;		//KAW 2013 05 08
		LocalConc_Agr_DIN = 0.0;		//KAW 2013 05 08
	}

	else {
     //       xMid = 60.0 + 12.0 * log(runoff);                 / RJS 110713 more accurately represents Wollheim et al. 2008 for runoff
     //         xMid           = 51.388 + 19.459 * log(runoff); //RJS 07-29-08	MMM changed from log10 to log 2013-03-13
     //         xMid           = 60.36 + 7.27 * log(runoff);     //RJS 03-09-15	based on raw data from wollheim
              xMid           = 51.388 + 8.4511 * log(runoff);     //RJS 03-23-15 based on raw data from wollheim and ONLY IS sites
              //		LocalConc_DIN  = asym / (1 + pow(2.718281828, (xMid - luSub) / scale)); // mg/L
            LocalConc_Sub_DIN  = asym / (1 + pow(2.718281828, (xMid - luSub) / scale)); 		// mg/L //KAW 2013 05 08
            LocalConc_Agr_DIN  = ((asym * multiplier) / (1 + pow(2.718281828, (xMid - luAgr) / scale))) - ((asym * multiplier) / (1 + pow(2.718281828, (xMid - 0.0) / scale))); 	// mg/L //KAW 2013 05 08 (to remove double counting of natural inputs) Concentrations in Agriculture headwater stream are 3.5x the concentrations in suburban stream (Price, unpublished data)
	}


//	LocalLoad_Sub_DIN  = runoffVol * 86400 * LocalConc_Sub_DIN / 1000; 		// kg/day //KAW 2013 05 08 Loading from suburban part of each grid cell
//	LocalLoad_Agr_DIN  = runoffVol * 86400 * LocalConc_Agr_DIN / 1000; 		// kg/day //KAW 2013 05 08 Loading from the agricultural part of each grid cell
//	LocalLoad_DIN      = LocalLoad_Agr_DIN + LocalLoad_Sub_DIN;                     // kg/day //KAW 2013 05 08: Sum of all N loadings  Total point is from MMM's Code for point source N inputs
       
 	LocalLoad_Sub_DIN  = (runoffVol * (1 - (luAgr/100))) * 86400 * LocalConc_Sub_DIN / 1000; 	// kg/day //RJS 02-15-2015
	LocalLoad_Agr_DIN  = (runoffVol * (luAgr/100)) * 86400 * LocalConc_Agr_DIN / 1000; 		// kg/day //RJS 02-15-2015 
	LocalLoad_DIN      = LocalLoad_Agr_DIN + LocalLoad_Sub_DIN;                                     // kg/day //RJS 02-15-2015
      
        
	if (loadAdjust > 0.0) LocalLoad_DIN =  LocalLoad_DIN * loadAdjust; 	// RJS 112211

 //       if (itemID == 1357) printf("itemID=%d, luAgr = %f, luSub = %f, load_agr = %f, load_sub = %f, totalLoad = %f\n",itemID,luAgr,luSub,LocalLoad_Agr_DIN,LocalLoad_Sub_DIN,LocalLoad_DIN);
        
        //if(LocalConc_Sub_DIN > 1.4 || LocalConc_Ag_DIN > 4.9) printf("runoff = %f, luSub = %f, luAg = %f, xMid = %f, Sub_conc = %f, Ag_conc = %f\n", runoff, luSub, luAg, xMid, LocalConc_Sub_DIN, LocalConc_Ag_DIN);
        //if (itemID == 809) printf("**** itemID = %d, y=%d, m=%d, d=%d, runoff=%f, luSub=%f, luAg=%f, xMid=%f, Sub_conc=%f, Ag_conc=%f, Subload=%f, Agload=%f, total=%f\n", itemID, MFDateGetCurrentYear(), MFDateGetCurrentMonth(), MFDateGetCurrentDay(), runoff, luSub, luAg, xMid, LocalConc_Sub_DIN, LocalConc_Ag_DIN, LocalLoad_Sub_DIN, LocalLoad_Ag_DIN, LocalLoad_DIN);
	//if (itemID == 31) printf("***** itemID = %d, year = %d, month = %d, day = %d, xMid = %f, luSub = %f, scale = %f, asym = %f\n", itemID, MFDateGetCurrentYear(), MFDateGetCurrentMonth(), MFDateGetCurrentDay(), xMid, luSub, scale, asym);
	//if (itemID == 31) printf("runoffVol = %f, runoff = %f, LocalLoad_DIN = %f, LocalConc_DIN = %f, Total_point= %f, loadAdjust = %f\n", runoffVol, runoff, LocalLoad_DIN, LocalConc_DIN, Total_point, loadAdjust);  //mmm commented out 2013-3-13

 //       printf("loadAdjust = %f\n", loadAdjust);

	MFVarSetFloat (_MDOutLocalLoad_DINID,       itemID, LocalLoad_DIN);	  	// RJS 090308
	MFVarSetFloat (_MDOutLocalConc_Sub_DINID,   itemID, LocalLoad_Sub_DIN);	  	// RJS 090308
	MFVarSetFloat (_MDOutDINSubLoadConcID, 	    itemID, LocalConc_Sub_DIN);  	// KAW 2013/03/15
	MFVarSetFloat (_MDOutLocalLoad_Agr_DINID,   itemID, LocalLoad_Agr_DIN);	  	// RJS 090308
	MFVarSetFloat (_MDOutDINAgLoadConcID, 	    itemID, LocalConc_Agr_DIN);		// KAW 2013/03/15
}


static void _MDNitrogenInputsCalc2 (int itemID) {

	// Input

	float luSub     = 0.0;
	float luAgr     = 0.0;
	float runoff    = 0.0;		// mm/day
	float runoffVol = 0.0;		// m3/sec

	float scale     = 0.0;
	float asym      = 0.0;
	float xMid      = 0.0;
        float Xmid_b    = 0.0;
        float Xmid_m    = 0.0;

	float LocalConc_DIN      = 0.0;
	float LocalLoad_DIN      = 0.0;

	float LocalConc_Sub_DIN  = 0.0;		//KAW 2013 05 08 Suburban Loading Concentration
	float LocalLoad_Sub_DIN  = 0.0;		//KAW 2013 05 08 Suburban Loading

	float LocalConc_Agr_DIN   = 0.0;		//KAW 2013 05 08 Agricultural Loading Concentration
	float LocalLoad_Agr_DIN   = 0.0;		//KAW 2013 05 08 Agricultural Loading

	float riverOrder          = 0.0;
	float loadAdjust          = 0.0;		// RJS 112211
	float Total_point         = 0.0;
	float multiplier	  = 0.0;		// Agr multipler for DIN loading concentration 
// test

	riverOrder = MFVarGetFloat (_MDInRiverOrderID, itemID, 0.0);	// RJS 090508
//	luSub = MFVarGetFloat (_MDInLandUseID,    itemID, 0.0);
	luSub = MFVarGetFloat (_MDInLandUseSubID, itemID, 0.0);	//KAW 2013 05 08 Suburban Land Use
	luAgr = MFVarGetFloat (_MDInLandUseAgID,  itemID, 0.0);	//KAW 2013 05 08 Agricultural Land Use
	asym  = MFVarGetFloat (_MDInAsymID,  itemID, 0.0);	// RJS 102314
        scale = MFVarGetFloat (_MDInScaleID,  itemID, 0.0);	// RJS 102314
        Xmid_b = MFVarGetFloat (_MDInXmid_bID,  itemID, 0.0);	// RJS 102314
        Xmid_m = MFVarGetFloat (_MDInXmid_mID,  itemID, 0.0);	// RJS 102314

	luSub = luSub + luAgr > 100.0 ? luSub / (luSub + luAgr) * 100.0 : luSub;        // RJS 110813    
        luAgr = luSub + luAgr > 100.0 ? luAgr / (luSub + luAgr) * 100.0 : luAgr;        // RJS 110813

        runoff             = MFVarGetFloat (_MDInRunoffID,         itemID, 0.0); 	// mm / d
	runoffVol          = MFVarGetFloat (_MDInRunoffVolID,      itemID, 0.0); 	// m3/sec
	loadAdjust	   = MFVarGetFloat (_MDInLoadAdjustID,     itemID, 1.0);	// RJS 112211, adjusts loads, keep at 1 if nodata
	multiplier	   = MFVarGetFloat (_MDInAgMultiplierID,   itemID, 0.0);	// RJS 051815

	if (runoff < 0.00001) {
		LocalConc_Sub_DIN = 0.0;		//KAW 2013 05 08
		LocalConc_Agr_DIN = 0.0;		//KAW 2013 05 08
	}

	else {
            xMid = Xmid_b + Xmid_m * log(runoff);   // RJS 110713 more accurately represents Wollheim et al. 2008 for runoff

            LocalConc_Sub_DIN  = asym / (1 + pow(2.718281828, (xMid - luSub) / scale)); 		// mg/L //KAW 2013 05 08
            LocalConc_Agr_DIN  = ((asym * multiplier) / (1 + pow(2.718281828, (xMid - luAgr) / scale))) - ((asym * multiplier) / (1 + pow(2.718281828, (xMid - 0.0) / scale))); 	// mg/L //KAW 2013 05 08 (to remove double counting of natural inputs) Concentrations in Agriculture headwater stream are 3.5x the concentrations in suburban stream (Price, unpublished data)
	}

//	LocalLoad_Sub_DIN  = runoffVol * 86400 * LocalConc_Sub_DIN / 1000; 		// kg/day //KAW 2013 05 08 Loading from suburban part of each grid cell
//	LocalLoad_Agr_DIN  = runoffVol * 86400 * LocalConc_Agr_DIN / 1000; 		// kg/day //KAW 2013 05 08 Loading from the agricultural part of each grid cell
//	LocalLoad_DIN      = LocalLoad_Agr_DIN + LocalLoad_Sub_DIN;                     // kg/day //KAW 2013 05 08: Sum of all N loadings  Total point is from MMM's Code for point source N inputs

 	LocalLoad_Sub_DIN  = (runoffVol * (1 - (luAgr/100))) * 86400 * LocalConc_Sub_DIN / 1000; 	// kg/day //RJS 02-15-2015
	LocalLoad_Agr_DIN  = (runoffVol * (luAgr/100)) * 86400 * LocalConc_Agr_DIN / 1000; 		// kg/day //RJS 02-15-2015 
	LocalLoad_DIN      = LocalLoad_Agr_DIN + LocalLoad_Sub_DIN;                                     // kg/day //RJS 02-15-2015
      
        
	if (loadAdjust > 0.0) LocalLoad_DIN =  LocalLoad_DIN * loadAdjust; 	// RJS 112211

        //if (itemID == 662) printf("itemID=%d, luAgr = %f, luSub = %f, load_agr = %f, load_sub = %f, totalLoad = %f\n",itemID,luAgr,luSub,LocalLoad_Agr_DIN,LocalLoad_Sub_DIN,LocalLoad_DIN);
        
        //if(LocalConc_Sub_DIN > 1.4 || LocalConc_Ag_DIN > 4.9) printf("runoff = %f, luSub = %f, luAg = %f, xMid = %f, Sub_conc = %f, Ag_conc = %f\n", runoff, luSub, luAg, xMid, LocalConc_Sub_DIN, LocalConc_Ag_DIN);
        //if (itemID == 809) printf("**** itemID = %d, y=%d, m=%d, d=%d, runoff=%f, luSub=%f, luAg=%f, xMid=%f, Sub_conc=%f, Ag_conc=%f, Subload=%f, Agload=%f, total=%f\n", itemID, MFDateGetCurrentYear(), MFDateGetCurrentMonth(), MFDateGetCurrentDay(), runoff, luSub, luAg, xMid, LocalConc_Sub_DIN, LocalConc_Ag_DIN, LocalLoad_Sub_DIN, LocalLoad_Ag_DIN, LocalLoad_DIN);
	//if (itemID == 662) printf("***** itemID = %d, year = %d, month = %d, day = %d, xMid = %f, luSub = %f, scale = %f, asym = %f\n", itemID, MFDateGetCurrentYear(), MFDateGetCurrentMonth(), MFDateGetCurrentDay(), xMid, luSub, scale, asym);
	//if (itemID == 662) printf("runoffVol = %f, runoff = %f, LocalLoad_DIN = %f, LocalConc_DIN = %f, Total_point= %f, loadAdjust = %f\n", runoffVol, runoff, LocalLoad_DIN, LocalConc_DIN, Total_point, loadAdjust);  //mmm commented out 2013-3-13
	//if (itemID == 662) printf("Xmid_b = %f, Xmid_m = %f, runoff = %f\n", Xmid_b, Xmid_m, runoff);  //mmm commented out 2013-3-13

	MFVarSetFloat (_MDOutLocalLoad_DINID,       itemID, LocalLoad_DIN);	  	// RJS 090308
	MFVarSetFloat (_MDOutLocalConc_Sub_DINID,   itemID, LocalLoad_Sub_DIN);	  	// RJS 090308
	MFVarSetFloat (_MDOutDINSubLoadConcID, 	    itemID, LocalConc_Sub_DIN);  	// KAW 2013/03/15
	MFVarSetFloat (_MDOutLocalLoad_Agr_DINID,   itemID, LocalLoad_Agr_DIN);	  	// RJS 090308
	MFVarSetFloat (_MDOutDINAgLoadConcID, 	    itemID, LocalConc_Agr_DIN);		// KAW 2013/03/15
}


enum {MDcalculate, MDcalculate2, MDinput, MDPnET, MDnone};

int MDNitrogenInputsDef () {

			int  optID = MFUnset;													    //RJS 10-28-10
			const char *optStr, *optName = MDOptDINInputs;								//RJS 10-28-10
			const char *options [] = { MDCalculateStr, MDCalculate2Str, MDInputStr, MDPnETStr, MDNoneStr, (char *) NULL };		//RJS 10-28-10

	MFDefEntering ("Nitrogen Inputs");

	if ((optStr = MFOptionGet (optName)) != (char *) NULL) optID = CMoptLookup (options, optStr, true);  //RJS 02-11-10

	switch (optID) {	//RJS 02-11-10

	case MDnone:
	if (((_MDInLandUseID             = MFVarGetID (MDVarLandUseSpatial,       "-",  MFInput, MFState, MFBoundary)) == CMfailed) ||
	    ((_MDInRunoffID   		 = MFVarGetID (MDVarRunoff,              "mm",  MFInput,  MFFlux, MFBoundary)) == CMfailed) ||
	    ((_MDInRunoffVolID           = MDRunoffVolumeDef ()) == CMfailed) ||
	    ((_MDInLawnFractionID        = MFVarGetID (MDVarLawnFraction,      "-", MFInput, MFState, MFBoundary)) == CMfailed) ||		//RJS 051111
	    ((_MDInRiverOrderID          = MFVarGetID (MDVarRiverOrder,  	  	  "-",  MFInput, MFState, MFBoundary)) == CMfailed) ||	// RJS 090308
	    ((_MDInLoadAdjustID          = MFVarGetID (MDVarLoadAdjust,           "-", MFInput, MFState, MFBoundary)) == CMfailed)  ||  // RJS 112211
	    ((_MDOutLocalLoad_DINID      = MFVarGetID (MDVarLocalLoadDIN,    "kg/day", MFOutput, MFFlux,  MFBoundary)) == CMfailed) ||  // RJS 090308
	    ((_MDOutLocalLoadnew_DINID   = MFVarGetID (MDVarLocalLoadDINnew, "kg/day", MFOutput, MFFlux,  MFBoundary)) == CMfailed) ||	// RJS 061511	created to compare new and old loads ONLY.
	    (MFModelAddFunction (_MDNitrogenInputsPartitioned) == CMfailed)) return (CMfailed);
		break;

	case MDinput:
	if (((_MDInLocalLoad_DINID       = MFVarGetID (MDVarInLocalLoadDIN,       "kg/day", MFInput, MFFlux, MFBoundary)) == CMfailed) || 
//          ((_MDInGrdWatMassRelID       = MFVarGetID (MDVarGroundWaterMassRel,   "kg/day", MFInput, MFFlux, MFBoundary)) == CMfailed) ||	//RJS 050511
//          ((_MDInRunoffPoolMassRelID   = MFVarGetID (MDVarRunoffPoolMassRel,    "kg/day", MFInput, MFFlux, MFBoundary)) == CMfailed) ||	//RJS 050511
//          ((_MDInRunoffVolID           = MDRunoffVolumeDef ()) == CMfailed) ||															//RJS 050511
	    ((_MDOutLocalLoad_DINID      = MFVarGetID (MDVarLocalLoadDIN,         "kg/day", MFOutput, MFFlux, MFBoundary)) == CMfailed) ||	//RJS 050511
	    (MFModelAddFunction (_MDNitrogenInputsInput) == CMfailed)) return (CMfailed);
		break;
                
        case MDPnET:
	if (((_MDInRunoffID   		 = MFVarGetID (MDVarRunoff,              "mm",  MFInput,  MFFlux, MFBoundary)) == CMfailed) ||	//RJS 10-28-10
	    ((_MDInRunoffVolID           = MDRunoffVolumeDef ()) == CMfailed) ||			// RJS 10-28-10
	    ((_MDInLandUseAgID           = MFVarGetID (MDVarLandUseSpatialAg,        "-",  MFInput, MFState, MFBoundary)) == CMfailed) ||
	    ((_MDInLandUseConID          = MFVarGetID (MDVarLandUseSpatialCon,       "-",  MFInput, MFState, MFBoundary)) == CMfailed) ||
	    ((_MDInLandUseDecID          = MFVarGetID (MDVarLandUseSpatialDec,       "-",  MFInput, MFState, MFBoundary)) == CMfailed) ||
	    ((_MDInLandUseMixID          = MFVarGetID (MDVarLandUseSpatialMix,       "-",  MFInput, MFState, MFBoundary)) == CMfailed) ||
            ((_MDInLawnFractionID        = MFVarGetID (MDVarLawnFraction,               "-",     MFInput,  MFState, MFBoundary)) == CMfailed) ||
            ((_MDInImpFractionID         = MFVarGetID (MDVarImpFracSpatial,  "-",    MFInput,  MFState, MFBoundary)) == CMfailed) ||
            ((_MDInWetlandsID            = MFVarGetID (MDVarFracWetlandArea,                             "-",    MFInput,  MFState, MFBoundary))   == CMfailed) ||           // RJS 112513
            ((_MDInH2OFractionID         = MFVarGetID (MDVarH2OFracSpatial,          "-",  MFInput, MFState, MFBoundary)) == CMfailed) ||		//commented out 082812
            ((_MDInLocalLoad_Con_DINID   = MFVarGetID (MDVarInLocalLoadConDIN,    "g/m2",  MFInput, MFFlux, MFBoundary)) == CMfailed) || 
            ((_MDInLocalLoad_Dec_DINID   = MFVarGetID (MDVarInLocalLoadDecDIN,    "g/m2",  MFInput, MFFlux, MFBoundary)) == CMfailed) || 
            ((_MDInLocalLoad_Mix_DINID   = MFVarGetID (MDVarInLocalLoadMixDIN,    "g/m2",  MFInput, MFFlux, MFBoundary)) == CMfailed) || 
            ((_MDInLocalLoad_Lak_DINID   = MFVarGetID (MDVarInLocalLoadLakDIN,    "g/m2",  MFInput, MFFlux, MFBoundary)) == CMfailed) || 
            ((_MDInLocalLoad_Wet_DINID   = MFVarGetID (MDVarInLocalLoadWetDIN,    "g/m2",  MFInput, MFFlux, MFBoundary)) == CMfailed) || 
            ((_MDInLocalLoad_Imp_DINID   = MFVarGetID (MDVarInLocalLoadImpDIN,    "g/m2",  MFInput, MFFlux, MFBoundary)) == CMfailed) || 
            ((_MDInLocalLoad_Lawn_DINID   = MFVarGetID (MDVarInLocalLoadLawnDIN,    "g/m2",  MFInput, MFFlux, MFBoundary)) == CMfailed) || 
            ((_MDInLocalLoad_Ag_DINID   = MFVarGetID (MDVarInLocalLoadAgDIN,    "g/m2",  MFInput, MFFlux, MFBoundary)) == CMfailed) || 
	    ((_MDInAsymID                = MFVarGetID (MDVarAsym,           "-", MFInput, MFState, MFBoundary)) == CMfailed)  ||  // RJS 051815           
	    ((_MDInAgMultiplierID        = MFVarGetID (MDVarAgMultiplier,         "-", MFInput, MFState, MFBoundary)) == CMfailed)  ||  // RJS 051815           
            ((_MDInPropROStormWaterID   = MFVarGetID (MDVarPropROStormWater,    "-",     MFOutput, MFState, MFBoundary)) == CMfailed) ||       // sz 20160211
            ((_MDInPropROSurfaceWaterID = MFVarGetID (MDVarPropROSurfaceWater,  "-",     MFOutput, MFState, MFBoundary)) == CMfailed) ||       // sz 20160211
            ((_MDInPropROGroundWaterID  = MFVarGetID (MDVarPropROGroundWater,   "-",     MFOutput, MFState, MFBoundary)) == CMfailed) ||       // sz 20160211
	    ((_MDOutLocalLoad_DINID      = MFVarGetID (MDVarLocalLoadDIN,         "kg/day", MFOutput, MFFlux, MFBoundary)) == CMfailed) ||	//RJS 050511
            ((_MDOutLocalConc_DINID      = MFVarGetID (MDVarLocalConcDIN,         "mg/L", MFOutput, MFState, MFBoundary)) == CMfailed) ||	//RJS 050511
	    ((_MDOutLocalLoad_PnET_DINID = MFVarGetID (MDVarLocalLoadPnETDIN,         "kg/day", MFOutput, MFFlux, MFBoundary)) == CMfailed) ||	//SZ 20150421
 //         ((_MDOutLocalConc_PnET_DINID = MFVarGetID (MDVarLocalConcPnETDIN,         "mg/L", MFOutput, MFState, MFBoundary)) == CMfailed) ||	//SZ 20150421
            ((_MDOutLocalLoad_Con_DINID  = MFVarGetID (MDVarLocalLoadConDIN,      "kg/day", MFOutput, MFFlux, MFBoundary)) == CMfailed) ||	//RJS 050511
	    ((_MDOutLocalLoad_Agr_DINID  = MFVarGetID (MDVarLocalLoadAgDIN,      "kg/day",  MFOutput, MFFlux, MFBoundary)) == CMfailed) ||	//RJS 050511
            ((_MDOutLocalLoad_Lak_DINID  = MFVarGetID (MDVarLocalLoadLakDIN,      "kg/day", MFOutput, MFFlux, MFBoundary)) == CMfailed) ||	//RJS 050511
            ((_MDOutLocalLoad_Dec_DINID  = MFVarGetID (MDVarLocalLoadDecDIN,      "kg/day", MFOutput, MFFlux, MFBoundary)) == CMfailed) ||	//RJS 050511
            ((_MDOutLocalLoad_Mix_DINID  = MFVarGetID (MDVarLocalLoadMixDIN,      "kg/day", MFOutput, MFFlux, MFBoundary)) == CMfailed) ||	//RJS 050511
            ((_MDOutLocalLoad_Wet_DINID  = MFVarGetID (MDVarLocalLoadWetDIN,      "kg/day", MFOutput, MFFlux, MFBoundary)) == CMfailed) ||	//RJS 050511
            ((_MDOutLocalLoad_Sub_DINID  = MFVarGetID (MDVarLocalLoadSubDIN,      "kg/day", MFOutput, MFFlux, MFBoundary)) == CMfailed) ||	//RJS 050511
	    ((_MDOutLocalConc_Sub_DINID  = MFVarGetID (MDVarLocalConcSubDIN,      "mg/L", MFOutput, MFState, MFBoundary)) == CMfailed) ||	//RJS 050511
	    ((_MDOutLocalConc_Agr_DINID  = MFVarGetID (MDVarLocalConcAgDIN,      "mg/L",  MFOutput, MFState, MFBoundary)) == CMfailed) ||	//RJS 050511
            ((_MDOutLocalConc_Lak_DINID  = MFVarGetID (MDVarLocalConcLakDIN,      "mg/L", MFOutput, MFState, MFBoundary)) == CMfailed) ||	//RJS 050511
            ((_MDOutLocalConc_Dec_DINID  = MFVarGetID (MDVarLocalConcDecDIN,      "mg/L", MFOutput, MFState, MFBoundary)) == CMfailed) ||	//RJS 050511
            ((_MDOutLocalConc_Mix_DINID  = MFVarGetID (MDVarLocalConcMixDIN,      "mg/L", MFOutput, MFState, MFBoundary)) == CMfailed) ||	//RJS 050511
            ((_MDOutLocalConc_Wet_DINID  = MFVarGetID (MDVarLocalConcWetDIN,      "mg/L", MFOutput, MFState, MFBoundary)) == CMfailed) ||	//RJS 050511
            ((_MDOutLocalConc_Con_DINID  = MFVarGetID (MDVarLocalConcConDIN,      "mg/L", MFOutput, MFState, MFBoundary)) == CMfailed) ||	//RJS 050511
            ((_MDOutTotalRiparianRemoval_DINID  = MFVarGetID (MDVarToralRiparianRemDIN,      "kg/d", MFOutput, MFFlux, MFBoundary)) == CMfailed) ||	//SZ 20160211
	    (MFModelAddFunction (_MDNitrogenInputsPnET) == CMfailed)) return (CMfailed);
            break;
                
	case MDcalculate:
	if (((_MDInRunoffID   			 = MFVarGetID (MDVarRunoff,              "mm",  MFInput,  MFFlux, MFBoundary)) == CMfailed) ||	//RJS 10-28-10
	    ((_MDInRunoffVolID           = MDRunoffVolumeDef ()) == CMfailed) ||			// RJS 10-28-10
//	    ((_MDInDINLoadConcID         = MFVarGetID (MDVarDINLoadConc,       "mg/L",  MFInput, MFFlux, MFBoundary))  == CMfailed) ||	//RJS 10-28-10
//	    ((_MDInLandUseID             = MFVarGetID (MDVarLandUseSpatial,       "-",  MFInput, MFState, MFBoundary)) == CMfailed) ||
	    ((_MDInLandUseSubID             = MFVarGetID (MDVarLandUseSpatialSub,       "-",  MFInput, MFState, MFBoundary)) == CMfailed) ||
	    ((_MDInLandUseAgID             	= MFVarGetID (MDVarLandUseSpatialAg,       "-",  MFInput, MFState, MFBoundary)) == CMfailed) ||
//	    ((_MDInTotalPointID			 = MDPointSourceDef ()) == CMfailed) ||															//MMM Added this so that Nitrogeninputs knows it must run MDPointSource in order to calculate localload
	    ((_MDInLoadAdjustID          = MFVarGetID (MDVarLoadAdjust,           "-", MFInput, MFState, MFBoundary)) == CMfailed)  ||  // RJS 112211
	    ((_MDInAsymID                = MFVarGetID (MDVarAsym,           "-", MFInput, MFState, MFBoundary)) == CMfailed)  ||  // RJS 051815           
	    ((_MDInAgMultiplierID        = MFVarGetID (MDVarAgMultiplier,         "-", MFInput, MFState, MFBoundary)) == CMfailed)  ||  // RJS 051815               
	    ((_MDOutLocalConc_Sub_DINID     = MFVarGetID (MDVarLocalLoadSubDIN,    "kg/day", MFOutput, MFFlux, MFBoundary))  == CMfailed) ||	//KAW 2013 05 08
	    ((_MDOutLocalLoad_Agr_DINID      = MFVarGetID (MDVarLocalLoadAgDIN,    "kg/day", MFOutput, MFFlux, MFBoundary))  == CMfailed) ||	//KAW 2013 05 08
	    ((_MDOutDINSubLoadConcID        = MFVarGetID (MDVarDINSubLoadConc,      "mg/L",  MFOutput, MFState, MFBoundary)) == CMfailed) || 	//KAW 2013 05 08
	    ((_MDOutDINAgLoadConcID        	= MFVarGetID (MDVarDINAgLoadConc,      "mg/L",  MFOutput, MFState, MFBoundary)) == CMfailed) || 	//KAW 2013 05 08
	    ((_MDOutLocalLoad_DINID      = MFVarGetID (MDVarLocalLoadDIN,    "kg/day", MFOutput, MFFlux, MFBoundary))  == CMfailed) ||	//RJS 10-28-10
	    ((_MDOutDINLoadConcID        = MFVarGetID (MDVarDINLoadConc,      "mg/L",  MFOutput, MFState, MFBoundary)) == CMfailed) || 	// KYLE
	    ((_MDInRiverOrderID          = MFVarGetID (MDVarRiverOrder,           "-", MFInput,  MFState, MFBoundary)) == CMfailed) ||  //RJS 10-29-10
	    (MFModelAddFunction (_MDNitrogenInputsCalc) == CMfailed)) return (CMfailed);
	   break;
           
        case MDcalculate2:
	if (((_MDInRunoffID   			 = MFVarGetID (MDVarRunoff,              "mm",  MFInput,  MFFlux, MFBoundary)) == CMfailed) ||	//RJS 10-28-10
	    ((_MDInRunoffVolID           = MDRunoffVolumeDef ()) == CMfailed) ||			// RJS 10-28-10
            ((_MDInLandUseSubID             = MFVarGetID (MDVarLandUseSpatialSub,       "-",  MFInput, MFState, MFBoundary)) == CMfailed) ||
	    ((_MDInLandUseAgID             	= MFVarGetID (MDVarLandUseSpatialAg,       "-",  MFInput, MFState, MFBoundary)) == CMfailed) ||
	    ((_MDInLoadAdjustID          = MFVarGetID (MDVarLoadAdjust,           "-", MFInput, MFState, MFBoundary)) == CMfailed)  ||  // RJS 112211
	    ((_MDInAsymID               = MFVarGetID (MDVarAsym,           "-", MFInput, MFState, MFBoundary)) == CMfailed)  ||  // RJS 102314
	    ((_MDInAgMultiplierID        = MFVarGetID (MDVarAgMultiplier,         "-", MFInput, MFState, MFBoundary)) == CMfailed)  ||  // RJS 051815               
  	    ((_MDInScaleID               = MFVarGetID (MDVarScale,           "-", MFInput, MFState, MFBoundary)) == CMfailed)  ||  // RJS 102314
	    ((_MDInXmid_bID               = MFVarGetID (MDVarXmid_b,           "-", MFInput, MFState, MFBoundary)) == CMfailed)  ||  // RJS 102314
	    ((_MDInXmid_mID               = MFVarGetID (MDVarXmid_m,           "-", MFInput, MFState, MFBoundary)) == CMfailed)  ||  // RJS 102314
	    ((_MDOutLocalConc_Sub_DINID     = MFVarGetID (MDVarLocalLoadSubDIN,    "kg/day", MFOutput, MFFlux, MFBoundary))  == CMfailed) ||	//KAW 2013 05 08
	    ((_MDOutLocalLoad_Agr_DINID      = MFVarGetID (MDVarLocalLoadAgDIN,    "kg/day", MFOutput, MFFlux, MFBoundary))  == CMfailed) ||	//KAW 2013 05 08
	    ((_MDOutDINSubLoadConcID        = MFVarGetID (MDVarDINSubLoadConc,      "mg/L",  MFOutput, MFState, MFBoundary)) == CMfailed) || 	//KAW 2013 05 08
	    ((_MDOutDINAgLoadConcID        	= MFVarGetID (MDVarDINAgLoadConc,      "mg/L",  MFOutput, MFState, MFBoundary)) == CMfailed) || 	//KAW 2013 05 08
	    ((_MDOutLocalLoad_DINID      = MFVarGetID (MDVarLocalLoadDIN,    "kg/day", MFOutput, MFFlux, MFBoundary))  == CMfailed) ||	//RJS 10-28-10
	    ((_MDOutDINLoadConcID        = MFVarGetID (MDVarDINLoadConc,      "mg/L",  MFOutput, MFState, MFBoundary)) == CMfailed) || 	// KYLE
	    ((_MDInRiverOrderID          = MFVarGetID (MDVarRiverOrder,           "-", MFInput,  MFState, MFBoundary)) == CMfailed) ||  //RJS 10-29-10
	    (MFModelAddFunction (_MDNitrogenInputsCalc2) == CMfailed)) return (CMfailed);
	   break;
	default: MFOptionMessage (optName, optStr, options); return (CMfailed);
	}

	MFDefLeaving ("Nitrogen Inputs");
	return (_MDOutLocalLoad_DINID);

}


