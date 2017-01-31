/******************************************************************************

GHAAS Water Balance/Transport Model V3.0
Global Hydrologic Archive and Analysis System
Copyright 1994-2007, University of New Hampshire

MDPointSource.c

rob.stewart@unh.edu
m.m.mineau@unh.edu

*******************************************************************************/


#include <stdio.h>
#include <cm.h>
#include <MF.h>
#include <MD.h>

// Inputs

static int _MDInWWTP_InKgPerDayID    = MFUnset;
static int _MDInWWTP_PopServedID     = MFUnset;
static int _MDInWWTP_RemovalID	     = MFUnset;
static int _MDInWWTP_TreatmentID     = MFUnset;
static int _MDInWWTP_FracDINID       = MFUnset;
static int _MDInWWTP_DeteriorationID = MFUnset;
static int _MDInWWTP_KgTNperPersonID = MFUnset;
static int _MDInWWTP_MaxRemovalID    = MFUnset;

// Outputs

static int _MDOutWWTP_OutKgPerDayID   = MFUnset;
static int _MDOutWWTP_InKgPerDay2ID    = MFUnset;

static void _MDPointSourceRychtecka (int itemID) {

	float WWTP_InKgPerDay   = 0.0;	// Total daily DIN input 
	float WWTP_PopServed	= 0.0;	// Total population served
	float WWTP_Treatment	= 0.0;	// Type of treatment
	float WWTP_Removal	= 0.0;	// Proportion removal at WWTP
	float WWTP_OutKgPerDay	= 0.0;	// Total Kg effluent at WWTP
        float WWTP_FracDIN      = 0.0;  // Fraction of total effluent N that is DIN
        float deterioration     = 0.0;  // deterioration of WWTP treatment
        float WWTP_Removal2     = 0.0;  // proportion removal at WWTP

	int day, month, year;
        
        WWTP_FracDIN               = MFVarGetFloat (_MDInWWTP_FracDINID,       itemID, 0.0);
	WWTP_InKgPerDay    	   = MFVarGetFloat (_MDInWWTP_InKgPerDayID,    itemID, 0.0) * 365 * WWTP_FracDIN;   // change for Monte Carlo
	WWTP_PopServed		   = MFVarGetFloat (_MDInWWTP_PopServedID,     itemID, 0.0);
	WWTP_Treatment		   = MFVarGetFloat (_MDInWWTP_TreatmentID,     itemID, 0.0);
        deterioration              = MFVarGetFloat (_MDInWWTP_DeteriorationID, itemID, 0.0);
        
	WWTP_Removal = WWTP_Treatment == 1.0 ? 0.311765 : WWTP_Removal;		// Primary treatment
	WWTP_Removal = WWTP_Treatment == 2.0 ? 0.760735 : WWTP_Removal;		// Secondary treatment
	WWTP_Removal = WWTP_Treatment == 3.0 ? 0.808088 : WWTP_Removal;		// Advanced treatment I
	WWTP_Removal = WWTP_Treatment == 4.0 ? 0.845588 : WWTP_Removal;		// Advanced treatment II
	WWTP_Removal = WWTP_Treatment == 5.0 ? 0.900000 : WWTP_Removal;		// Community Amenities

        WWTP_Removal2 = WWTP_Removal * (1.0 - deterioration) > 0.9999 ? 0.9999 : WWTP_Removal * (1.0 - deterioration);    // This is for individual runs of deterioration (NOT MONTECARLO)
 //         WWTP_Removal2 = WWTP_Removal + deterioration > 0.9999 ? 0.9999 : WWTP_Removal + deterioration;                    // This is for MonteCarlo that uses standard deviations around mean

        
	WWTP_OutKgPerDay = WWTP_InKgPerDay * (1.0 - WWTP_Removal2) * 1.0;        // 50% of N is DIN (added 110413)

 //       if (WWTP_InKgPerDay > 0) {
 //       printf("deterioration = %f, WWTPout = %f, WWTPin = %f, WWTPrem = %f, WWTPrem2 = %f, WWTPfrac = %f\n", deterioration, WWTP_OutKgPerDay, WWTP_InKgPerDay, WWTP_Removal, WWTP_Removal2, WWTP_FracDIN);         
 //       }
        
	day   = MFDateGetCurrentDay ();
	month = MFDateGetCurrentMonth ();
	year  = MFDateGetCurrentYear ();

	if ((WWTP_InKgPerDay > 0.0) && (WWTP_Treatment == 0.0)) printf("itemID = %f, %d-%d-%d, WWTP_InKgPerDay = %f, WWTP_PopServed = %f, WWTP_Treatment = %f, WWTP_Removal = %f\n", itemID, year, month, day, WWTP_InKgPerDay, WWTP_PopServed, WWTP_Treatment, WWTP_Removal);

	MFVarSetFloat (_MDOutWWTP_OutKgPerDayID,   itemID, WWTP_OutKgPerDay);
	MFVarSetFloat (_MDOutWWTP_InKgPerDay2ID,   itemID, WWTP_InKgPerDay);

}

static void _MDPointSourceVanDrecht (int itemID) {

	int   WWTP_PopServed	 = 0.0;	// Total population served
	float WWTP_Treatment	 = 0.0;	// Type of treatment
	float WWTP_Removal	 = 0.0;	// Proportion removal at WWTP
	float WWTP_InKgPerDay	 = 0.0;	// Total Kg influent at WWTP
	float WWTP_OutKgPerDay	 = 0.0;	// Total Kg effluent at WWTP
        float WWTP_FracDIN       = 0.0;  // Fraction of total effluent N that is DIN
        float deterioration      = 0.0;  // deterioration of WWTP treatment
        float WWTP_Removal2      = 0.0;  // proportion removal at WWTP
        float WWTP_KgTNperPerson = 0.0;  // Kg TN per person per day
        float WWTP_MaxRemoval    = 0.0;  // highest proportion removal of WWTP in model

	int day, month, year;
        
	WWTP_PopServed		   = MFVarGetFloat (_MDInWWTP_PopServedID,     itemID, 0.0);
	WWTP_Treatment		   = MFVarGetFloat (_MDInWWTP_TreatmentID,     itemID, 0.0);
	WWTP_KgTNperPerson	   = MFVarGetFloat (_MDInWWTP_KgTNperPersonID, itemID, 0.0);
 	WWTP_MaxRemoval 	   = MFVarGetFloat (_MDInWWTP_MaxRemovalID,    itemID, 0.0);
        deterioration              = MFVarGetFloat (_MDInWWTP_DeteriorationID, itemID, 0.0);
        
	WWTP_Removal = WWTP_Treatment == 1.0 ? 0.00 : WWTP_Removal;		// Primary treatment
	WWTP_Removal = WWTP_Treatment == 2.0 ? 0.35 : WWTP_Removal;		// Secondary treatment
	WWTP_Removal = WWTP_Treatment == 3.0 ? 0.80 : WWTP_Removal;		// Advanced treatment I
	WWTP_Removal = WWTP_Treatment == 4.0 ? 0.80 : WWTP_Removal;		// Advanced treatment II
	WWTP_Removal = WWTP_Treatment == 5.0 ? 0.90 : WWTP_Removal;		// Community Amenities

        WWTP_PopServed = abs(WWTP_PopServed);
        WWTP_Removal2 = WWTP_Removal * (1.0 - deterioration) > 0.9999 ? 0.9999 : WWTP_Removal * (1.0 - deterioration);    // This is for individual runs of deterioration (NOT MONTECARLO)
 //         WWTP_Removal2 = WWTP_Removal + deterioration > 0.9999 ? 0.9999 : WWTP_Removal + deterioration;                    // This is for MonteCarlo that uses standard deviations around mean

        WWTP_InKgPerDay = (WWTP_KgTNperPerson * WWTP_PopServed) * (0.485 + ((WWTP_Removal2 / WWTP_MaxRemoval) * 0.255));    // Dumont et al. 2005 for calculation of DIN from TN       
	WWTP_OutKgPerDay = WWTP_InKgPerDay * (1.0 - WWTP_Removal2);        

//        if (WWTP_InKgPerDay > 0) {
//        printf("deterioration = %f, Pop = %f, WWTPout = %f, WWTPin = %f, WWTPrem = %f, WWTPrem2 = %f, WWTPfrac = %f\n", deterioration, WWTP_PopServed, WWTP_OutKgPerDay, WWTP_InKgPerDay, WWTP_Removal, WWTP_Removal2);         
//        }
        
	day   = MFDateGetCurrentDay ();
	month = MFDateGetCurrentMonth ();
	year  = MFDateGetCurrentYear ();


	MFVarSetFloat (_MDOutWWTP_OutKgPerDayID,   itemID, WWTP_OutKgPerDay);
	MFVarSetFloat (_MDOutWWTP_InKgPerDay2ID,   itemID, WWTP_InKgPerDay);

}


enum {MDRychtecka, MDVanDrecht, MDnone};

int MDPointSourceDef() {

    	int  optID = MFUnset;													    
	const char *optStr, *optName = MDOptWWTP;								
	const char *options [] = { MDRychteckaStr, MDVanDrechtStr, MDNoneStr, (char *) NULL };		

	MFDefEntering ("WWTP Point Source");

	if ((optStr = MFOptionGet (optName)) != (char *) NULL) optID = CMoptLookup (options, optStr, true);  

	switch (optID) {

     	case MDRychtecka:
	 if (((_MDInWWTP_InKgPerDayID	    = MFVarGetID (MDVarWWTPInKgPerDay,	      "kg/d",  MFInput,   MFFlux,  MFBoundary)) == CMfailed) ||
	     ((_MDInWWTP_PopServedID	    = MFVarGetID (MDVarWWTPPopServed,         "-",     MFInput,   MFFlux,  MFBoundary)) == CMfailed) ||
             ((_MDInWWTP_TreatmentID	    = MFVarGetID (MDVarWWTPTreatment,         "-",     MFInput,   MFState, MFBoundary)) == CMfailed) || 
             ((_MDInWWTP_FracDINID          = MFVarGetID (MDVarWWTPFracDIN,           "-",     MFInput,   MFState, MFBoundary)) == CMfailed) ||
             ((_MDOutWWTP_OutKgPerDayID	    = MFVarGetID (MDVarWWTPOutKgPerDay,       "kg/d",  MFOutput,  MFFlux,  MFBoundary)) == CMfailed) ||
             ((_MDOutWWTP_InKgPerDay2ID	    = MFVarGetID (MDVarWWTPInKgPerDay2,       "kg/d",  MFOutput,  MFFlux,  MFBoundary)) == CMfailed) ||
             ((_MDInWWTP_DeteriorationID    = MFVarGetID (MDVarWWTPDeterioration,      "-",     MFInput,   MFState, MFBoundary)) == CMfailed) ||
             (MFModelAddFunction(_MDPointSourceRychtecka) == CMfailed)) return (CMfailed);
	break;
        
  	case MDVanDrecht:
	 if (((_MDInWWTP_PopServedID	    = MFVarGetID (MDVarWWTPPopServed,         "-",     MFInput,   MFState,  MFBoundary)) == CMfailed) ||
             ((_MDInWWTP_TreatmentID	    = MFVarGetID (MDVarWWTPTreatment,         "-",     MFInput,   MFState, MFBoundary)) == CMfailed) || 
             ((_MDInWWTP_KgTNperPersonID    = MFVarGetID (MDVarWWTPKgTNperPerson,     "-",     MFInput,   MFState, MFBoundary)) == CMfailed) || 
             ((_MDInWWTP_MaxRemovalID	    = MFVarGetID (MDVarWWTPMaxRemoval,        "-",     MFInput,   MFState, MFBoundary)) == CMfailed) || 
             ((_MDOutWWTP_OutKgPerDayID	    = MFVarGetID (MDVarWWTPOutKgPerDay,       "kg/d",  MFOutput,  MFFlux,  MFBoundary)) == CMfailed) ||
             ((_MDOutWWTP_InKgPerDay2ID	    = MFVarGetID (MDVarWWTPInKgPerDay2,       "kg/d",  MFOutput,  MFFlux,  MFBoundary)) == CMfailed) ||
             ((_MDInWWTP_DeteriorationID    = MFVarGetID (MDVarWWTPDeterioration,      "-",     MFInput,   MFState, MFBoundary)) == CMfailed) ||
             (MFModelAddFunction(_MDPointSourceVanDrecht) == CMfailed)) return (CMfailed);
	break;
        default: MFOptionMessage (optName, optStr, options); return (CMfailed);
        }
        
	MFDefLeaving ("WWTP Point Source");
	return (_MDOutWWTP_OutKgPerDayID);

}

