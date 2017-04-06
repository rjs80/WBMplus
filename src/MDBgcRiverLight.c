/******************************************************************************

GHAAS Water Balance/Transport Model V3.0
Global Hydrologic Archive and Analysis System
Copyright 1994-2011, UNH - CCNY/CUNY

MDWRiverLight.c

wil.wollheim@unh.edu and ken.r.sheehan@gmail.com

Estimate light inputs to water surface and to river bottom

*******************************************************************************/

#include <stdio.h>
#include <string.h>
#include <cm.h>
#include <MF.h>
#include <MD.h>
#include <math.h>

// Model
static int _MDRiverLightID          = MFUnset;

// Input
static int _MDInDischargeID         = MFUnset;
static int _MDInFluxDOCID           = MFUnset; // here to get function call only
static int _MDInRiverbedWidthMeanID = MFUnset;
static int _MDInRiverDepthID        = MFUnset;
static int _MDInRiverWidthID        = MFUnset;
static int _MDInSolarRadID          = MFUnset;
static int _MDInResStorageChangeID  = MFUnset;
static int _MDInResStorageID        = MFUnset;
static int _MDInConcDOCID           = MFUnset;
static int _MDInKoppenID            = MFUnset;
static int _MDInCanopyCoverID       = MFUnset;
static int _MDInTotalPARID          = MFUnset;
static int _MDInSCALERBasinID       = MFUnset;
static int _MDInDrainageAreaID      = MFUnset;
static int _MDInTauID               = MFUnset;
static int _MDInPhiID               = MFUnset;
static int _MDInOrderSwitchID       = MFUnset;
static int _MDInRiverOrderID        = MFUnset;
static int _MDInRiparianOverHangID  = MFUnset;
// Output
static int _MDOutPARBenthicID          = MFUnset;
static int _MDOutPAR2ReachID           = MFUnset;
static int _MDOutPARSurfaceID          = MFUnset;
static int _MDOutPARID                 = MFUnset;
static int _MDOutCanopyCover2ID        = MFUnset;

static void _MDRiverLightCalc (int itemID) {    
  // This function (MDRiverLightCalc) was written by WW and requires DOC to be calculated (this function is untested)
    
    float discharge;
    float channelDepth;
    float channelWidth;
    float channelWidthMean;
    float width_canopy;   // mean width of overhanging trees from each bank
    float canopy_proportion = 0.7; // proportion of total bank length (stream length * 2) with canopy
    float LAI = 0.9;   // proportion of light shaded out
    float DOC;
    int koppen;	 

    //light parameters
    float solarRad;  //MJ/m2/
    float Ecan;
    float rad2par = 0.5;     // proportion of total radiation as PAR
    float reflectance = 0.9; // proportion of PAR transmitted across air-water interface
    float shading;
    float turbidity = 10;    // NTU
    float Kturb = 0.177;     // from Julian et al. 2008
    float Kdoc  = 0.05;      // assumption: assume every mg/l of DOC leads to 1% attenuation (m-1)
    float Kd;                // m-1
    float par2bottom;
    float par2reach;
   	
    discharge        = MFVarGetFloat ( _MDInDischargeID,        itemID, 0.0);
    solarRad         = MFVarGetFloat ( _MDInSolarRadID,         itemID, 0.0); //MJ/m2/d - CHECK UNITS - already accounts for clouds
    channelDepth     = MFVarGetFloat ( _MDInRiverDepthID,       itemID, 0.0);
    channelWidth     = MFVarGetFloat ( _MDInRiverWidthID,       itemID, 0.0);
    channelWidthMean = MFVarGetFloat (_MDInRiverbedWidthMeanID, itemID, 0.0);  //assume mean width = bankfull width
    DOC              = MFVarGetFloat (_MDInConcDOCID,           itemID, 0.0) * 1000; //converts kg/m3 to mg/l
    koppen           = MFVarGetInt   (_MDInKoppenID,            itemID, 0.0); // % TODO I think this is not % but category that should be read as integer FBM

     // KoppenID's : 1=Taiga, 2=semi-arid,3=tundra,4=temperate,5=tropics
     // Canopy width (m), assumption: Tundra=0,Taiga=2,Temperate=5,Semi-arid=0, WetTropics=10, DryTropics=0
     switch (koppen) { // TODO Wil originally read this as float and tested for nan which should never happen particularly, when it is read as integer.
     case 1: width_canopy = 2;  break;
     case 2: width_canopy = 0;  break;
     case 3: width_canopy = 0;  break;
     case 4: width_canopy = 5;  break;
     case 5: width_canopy = 10; break;
     }

	  Ecan = solarRad * rad2par;
	  if (!isnan(discharge) && (channelWidthMean > 0) && (channelWidth > 0)){
		  if ((channelWidthMean / 2) < width_canopy){
			  shading = canopy_proportion * MDMinimum(1, LAI);
		  }
		  else {
			  //proportion with canopy * proportion of half width not shaded
			  // does not include effect of variable width within channel; need to scale down based on actual water surface area
			  shading = canopy_proportion * (width_canopy / (channelWidthMean / 2)) * MDMinimum(1,LAI);
		  }
		  if (discharge > 0){
			  // turbidity = turbidity_const * pow(Q, turbidity_slope);
			  Kd = Kturb * turbidity + Kdoc * DOC; //Figure 5; Julian et al. 2008
			  par2bottom = Ecan * (1 - shading) * reflectance * pow(2.71828, -1 * Kd * channelDepth);
			  if (isnan(par2bottom)){
				  printf("Light: Q %f solarRad %f DOC %f Ecan %f shading %f reflectance %f Kd %f channelWidthMean %f channelDepth %f par2bottom %f \n",
						  discharge, solarRad, DOC, Ecan, shading, reflectance, Kd, channelWidthMean, channelDepth, par2bottom);
			  }
		  }
		  else {
			  par2bottom = Ecan * shading;
		  }
		  par2reach = par2bottom * channelWidth * MFModelGetLength(itemID);   //MJ/reach
		  MFVarSetFloat(_MDOutPARBenthicID, itemID, par2bottom);
		  MFVarSetFloat(_MDOutPAR2ReachID, itemID, par2reach);
	  }
	  else{
		  MFVarSetFloat(_MDOutPARBenthicID, itemID, 0.0);
		  MFVarSetFloat(_MDOutPAR2ReachID, itemID, 0.0);
	  }
}


static void _MDRiverLightInput (int itemID) {    
  // This function (MDRiverLightInput) was written by Ken Sheehan and uses a constant value for DOC)
   
    float PAR                = 0.0;     // mol/m2/day
    float TotalPAR           = 0.0;     // micromol/m2/second
    float Kd                 = 1.6;     // parameter in Julian et al. 2008 (can range from 0.8 to 1.6: 0.8 for CWT, 1.6 for KNZ)
    float refl               = 0.9;     // parameter in Julian et al. 2008 for reflectance, light availability in rivers
//    float c_canopy_shade     = 0.0;   // coniferous canopy cover percentage (KRS calc)
//    float d_canopy_shade     = 0.0;   // deciduous " (KRS calc)
//    float m_canopy_shade     = 0.0;   // mixed "  (KRS calc))
//    float canopy_shade_total = 0.0;   // total canopy cover
    float PAR_surface        = 0.0;     // mol/m2/day reaching 1m above water surface
    float Ecan               = 0.0;     // mol/m2/day total PAR coming in from Sun (includes cloud coverage)
    float PAR_benthic        = 0.0;     // mol/m2/day total PAR reaching benthic surface. For smaller streams where light actually reaches benthic surface
    float canopy_cover       = 0.0;     // proportion of the sky that is covered by canopy based on densiometer
    float depth              = 0.0;     // river depth (m)
    float PAR_surface_MJm2d  = 0.0;     // PAR surface in MJ/m2/d
    float PAR_benthic_MJm2d  = 0.0;     // PAR benthic in MJ/m2/d
    float SCALER_BasinID     = 0.0;     // SCALER Basin ID
    float DA                 = 0.0;     // Drainage Area (km2)  
    float overHang           = 1.0;     // length of canopy overhang on one side of the river
    float width              = 0.0;     // river width
    float tau                = 0.0;
    float phi                = 0.0;
    float Q_out              = 0.0;
    float riverOrder         = 0.0;
    float orderSwitch        = 0.0;
    
    depth                    = MFVarGetFloat (_MDInRiverDepthID,  itemID, 0.0);			// meters
    width                    = MFVarGetFloat (_MDInRiverWidthID,  itemID, 0.0);			// meters
    canopy_cover             = MFVarGetFloat (_MDInCanopyCoverID, itemID, 0.0);			// proportion
    TotalPAR                 = MFVarGetFloat (_MDInTotalPARID,    itemID, 0.0);                 // micromol/m2/second
    SCALER_BasinID           = MFVarGetFloat (_MDInSCALERBasinID, itemID, 0.0);                 // BNZ = 1, CWT = 2, KNZ = 3, LUQ = 4, MER = 5, TLK = 6
    DA                       = MFVarGetFloat (_MDInDrainageAreaID, itemID, 0.0);
    tau                      = MFVarGetFloat (_MDInTauID,   itemID, 0.0);
    phi                      = MFVarGetFloat (_MDInPhiID,   itemID, 0.0);
    Q_out                    = MFVarGetFloat (_MDInDischargeID,       itemID, 0.0);  	// cubic meters/second m3/sec - Discharge leaving this grid cell after routing.
    riverOrder               = MFVarGetFloat (_MDInRiverOrderID,   itemID, 0.0);
    orderSwitch              = MFVarGetFloat (_MDInOrderSwitchID,   itemID, 0.0);

    overHang                 = MFVarGetFloat (_MDInRiparianOverHangID,  itemID, 0.0);
    
//printf("#1 itemID = %d, order = %f, switch = %f\n", itemID, riverOrder, orderSwitch);

    if (riverOrder >= orderSwitch) {
    
    PAR   = TotalPAR / 1000000 * 86400; // converts TotalPAR to mol/m2/day from micromol/m2/second
          
 //    if (MFDateGetCurrentYear() > 2000)  printf("itemID = %d, tau = %f, phi = %f, width = %f\n", itemID, tau, phi, width);
    
        if (SCALER_BasinID == 1.0) {          // BNZ
            canopy_cover = canopy_cover;
        }
    
        else if (SCALER_BasinID == 2.0) {     // CWT
            canopy_cover = (91.5 - 0.309 * DA) / 100;
        }
    
        else if (SCALER_BasinID == 3.0) {     // KNZ
            canopy_cover = canopy_cover;
        }

        else if (SCALER_BasinID == 4.0) {     // LUQ
            canopy_cover = DA >= 1.0 ? (79.5 - 13.1 * log10(DA)) / 100 : (91.9 - 12.4 * DA) / 100;  // changed from 86.8 to 79.5
            canopy_cover = canopy_cover >= 1.0 ? 1.0 : canopy_cover;
        }
    
        else if (SCALER_BasinID == 5.0) {     // MER
          canopy_cover = MDMinimum((2.0 * overHang) / width, 1.0);
        }
    
        else if (SCALER_BasinID == 6.0) {     // TLK
            canopy_cover = canopy_cover;
        }
  
//      The 5 lines below calculate canopy cover using a very site-specific function derived in Coweeta.  Comment out these lines if (1) reading in canopy cover or (2) using a more general calculation of canopy cover)
//	c_canopy_shade		= 0.1 * (((2 * Length)/(Length * Width2) * 0.02 + (1 - (2 * Length)/(Length * Width2)))); /* 0.2 is a place holder for proportion of canopy that is coniferous in the reach*/
//	d_canopy_shade		= (month > 4) && (month < 10) ? 0.3: 0.7 * (((2 * Length)/(Length * Width2) * 0.02  + (1 - (2 * Length)/(Length * Width)))); /* 0.6 is a place holder for proportion of canopy that is deciduous in the reach*/
//	m_canopy_shade		= 0.2 * (((2 * Length)/(Length * Width2) * 0.02 + (1 - (2 * Length)/(Length * Width2))));
//	canopy_shade_total	= Width2 < 0.2 ? 0.98 : c_canopy_shade + d_canopy_shade + m_canopy_shade;
//      PAR_surface		= PAR * (1 - canopy_shade_total); // proportion of PAR reaching 1 meter above river surface for Coweeta
	
        PAR_surface		= PAR * (1 - 0.85 * canopy_cover); // Using Canopy Cover input grid.  // Canopy absorbs 85% of PAR: SZ 20170406 Grahm and Casper 2002
	Ecan			= PAR;
	PAR_benthic		= (Ecan * (1 - canopy_cover) * refl * exp(-Kd * depth)); //<-- Coweeta Function for Benthic PAR
//	PAR_benthic		= (Ecan * (1 - canopy_cover) * refl * exp(-Kd * depth)) / 24 / 0.0036; //<-- Coweeta Function for Benthic PAR

	
//       printf("#2 ID = %d, %d-%d-%d, depth = %f, Kd = %f, canopy_cover = %f, DA = %f, PAR = %f, PAR_surface = %f, PAR_benthic = %f\n", itemID, MFDateGetCurrentYear(), MFDateGetCurrentMonth(), MFDateGetCurrentDay(), depth, -Kd, canopy_cover, DA, PAR, PAR_surface, PAR_benthic);
    }
//       printf("#3 ID = %d, %d-%d-%d, depth = %f, Kd = %f, canopy_cover = %f, DA = %f, PAR = %f, PAR_surface = %f, PAR_benthic = %f\n", itemID, MFDateGetCurrentYear(), MFDateGetCurrentMonth(), MFDateGetCurrentDay(), depth, -Kd, canopy_cover, DA, PAR, PAR_surface, PAR_benthic);

        PAR_surface_MJm2d       = PAR_surface * 1000000 / 2.0513;   // MJ/m2/d
        PAR_benthic_MJm2d       = PAR_benthic * 1000000 / 2.0513;   // MJ/m2/d
        PAR                     = PAR * 1000000 / 2.0513;           // MJ/m2/d
        
  //      if (itemID == 17738) printf("canopy_cover = %f, PAR_benthic = %f\n", canopy_cover, PAR_benthic);

        
//   	MFVarSetFloat (_MDOutCanopyShadeTotID, 		  itemID, 			  canopy_shade_total);     
    	MFVarSetFloat (_MDOutPARBenthicID, 		  itemID, 			  PAR_benthic_MJm2d);
	MFVarSetFloat (_MDOutPARSurfaceID, 		  itemID, 			  PAR_surface_MJm2d);
 	MFVarSetFloat (_MDOutPARID,                       itemID, 			  PAR);
        MFVarSetFloat (_MDOutCanopyCover2ID,              itemID,                         canopy_cover);

}






enum { MDinput, MDcalc };

int MDBgcRiverLightDef () {
    int  optID = MFUnset;													    
    const char *optStr, *optName = MDOptRiverLight;
    const char *options [] = { MDInputStr, MDCalculateStr, (char *) NULL };

    if (_MDRiverLightID != MFUnset) return (_MDRiverLightID);

	MFDefEntering ("Calculate river light");
	
        if ((optStr = MFOptionGet (optName)) != (char *) NULL) optID = CMoptLookup (options, optStr, true);
	
  	switch (optID) {

        case MDcalc:
     if (((_MDInFluxDOCID           = MDBgcDOCRoutingDef   ()) == CMfailed) ||
         ((_MDInSolarRadID          = MFVarGetID (MDVarSolarRadiation,     "MJ/m2/d", MFInput,  MFFlux,  MFBoundary)) == CMfailed) ||
         ((_MDInDischargeID         = MFVarGetID (MDVarDischarge,          "m3/s",    MFInput,  MFState, MFBoundary)) == CMfailed) ||
         ((_MDInConcDOCID           = MFVarGetID (MDVarDOCConcentration,   "kg/m3",   MFInput,  MFState, MFBoundary)) == CMfailed) ||
         ((_MDInRiverWidthID        = MFVarGetID (MDVarRiverWidth,         "m",       MFInput,  MFState, MFBoundary)) == CMfailed) ||
         ((_MDInRiverbedWidthMeanID = MFVarGetID (MDVarRiverbedWidthMean,  "m",       MFInput,  MFState, MFBoundary)) == CMfailed) ||
 	 ((_MDInRiverDepthID        = MFVarGetID (MDVarRiverDepth,         "m",       MFInput,  MFState, MFBoundary)) == CMfailed) ||
 	 ((_MDInKoppenID            = MFVarGetID (MDVarKoppen,             MFNoUnit,  MFInput,  MFState, MFBoundary)) == CMfailed) ||
         ((_MDOutPARBenthicID       = MFVarGetID (MDVarPARBenthic,         "MJ/m2/d", MFOutput, MFState, MFBoundary))  == CMfailed) ||
         ((_MDOutPAR2ReachID        = MFVarGetID (MDVarPAR2Reach,          "MJ/d",    MFOutput, MFState, MFBoundary))  == CMfailed) ||

         (MFModelAddFunction (_MDRiverLightCalc) == CMfailed)) return (CMfailed);
      	 break;

          case MDinput:
     if (((_MDInCanopyCoverID       = MFVarGetID (MDVarCanopyCover,        "-",   MFInput, MFState, MFBoundary)) == CMfailed) ||
      	 ((_MDInRiverDepthID        = MFVarGetID (MDVarRiverDepth,         "m",   MFInput, MFState, MFBoundary)) == CMfailed) ||
      	 ((_MDInRiverWidthID        = MFVarGetID (MDVarRiverWidth,         "m",   MFInput, MFState, MFBoundary)) == CMfailed) ||
      	 ((_MDInTotalPARID          = MFVarGetID (MDVarTotalPAR,   "um/m2/sec",   MFInput, MFState, MFBoundary)) == CMfailed) ||
      	 ((_MDInSCALERBasinID       = MFVarGetID (MDVarSCALERBasinID,   "-",      MFInput, MFState, MFBoundary)) == CMfailed) ||
       	 ((_MDInDrainageAreaID      = MFVarGetID (MDVarDrainageArea,  "km2",      MFInput, MFState, MFBoundary)) == CMfailed) ||
         ((_MDInDischargeID         = MFVarGetID (MDVarDischarge,    "m3/s",      MFInput, MFState, MFBoundary)) == CMfailed) ||
         ((_MDInTauID               = MFVarGetID (MDVarTau,             "-",      MFInput, MFState, MFBoundary)) == CMfailed) ||
         ((_MDInPhiID               = MFVarGetID (MDVarPhi,             "-",      MFInput, MFState, MFBoundary)) == CMfailed) ||                             
         ((_MDInRiverOrderID        = MFVarGetID (MDVarRiverOrder,      "-",      MFInput, MFState, MFBoundary)) == CMfailed) ||
         ((_MDInOrderSwitchID       = MFVarGetID (MDVarOrderSwitch,     "-",      MFInput, MFState, MFBoundary)) == CMfailed) ||    
         ((_MDInRiparianOverHangID  = MFVarGetID (MDVarRiparianOverHang,"-",      MFInput, MFState, MFBoundary)) == CMfailed) ||
         ((_MDOutPARBenthicID       = MFVarGetID (MDVarPARBenthic,   "MJ/m2/d",  MFOutput, MFState, MFBoundary)) == CMfailed) ||
         ((_MDOutPARSurfaceID       = MFVarGetID (MDVarPARSurface,   "MJ/m2/d",  MFOutput, MFState, MFBoundary)) == CMfailed) ||
         ((_MDOutPARID              = MFVarGetID (MDVarPAR,          "MJ/m2/d",  MFOutput, MFState, MFBoundary)) == CMfailed) ||
         ((_MDOutCanopyCover2ID     = MFVarGetID (MDVarCanopyCover2,       "-",  MFOutput, MFState, MFBoundary)) == CMfailed) ||

         (MFModelAddFunction (_MDRiverLightInput) == CMfailed)) return (CMfailed);
      	 break;
         
        }
       
	   MFDefLeaving ("Calculate river light");
	   return (_MDOutPARBenthicID);
       
}
