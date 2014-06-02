/******************************************************************************
GHAAS Water Balance/Transport Model V3.0
Global Hydrologic Archive and Analysis System
Copyright 1994-2011, UNH - CCNY/CUNY

MDReservoirs.c

nehsani00@ccny.cuny.edu
dominik.wisser@unh.edu

Updated with a Neural Network function for
Reservoir Operation.
2013 NE
*******************************************************************************/
#include <stdio.h>
#include <cm.h>
#include <MF.h>
#include <MD.h>
#include <math.h>

// Input
static int _MDInDischargeID   = MFUnset;
static int _MDInDischMeanID   = MFUnset;
static int _MDInResCapacityID = MFUnset;

// Output
static int _MDOutResStorageID      = MFUnset;
static int _MDOutPreResStorageID   = MFUnset;
static int _MDOutResStorageChgID   = MFUnset;
static int _MDOutResReleaseID      = MFUnset;
static int _MDOutResRelease_t_1_ID = MFUnset;
static int _MDOutResRelease_t_2_ID = MFUnset;
static int _MDOutResRelease_t_3_ID = MFUnset;
static int _MDOutDisch_t_1_ID      = MFUnset;
static int _MDOutDisch_t_2_ID      = MFUnset;
static int _MDOutDisch_t_3_ID      = MFUnset;


static float ANNOUTPUT(float I1[3][1], float I2[2][1], float I3) {

    float FirstLayerBias[6][1] = {
        {2.5050069210611214},
        {0.56324962846021676},
        {2.9934573505209836},
        {-0.32326573637600164},
        {4.6033876322491265},
        {1.2269056039304829},

    };

    float FirstLayerWeight [6][3] = {
        {-2.8105005056468633,-1.0077666925314335,4.1021996487606796},
        {-1.1210513651438023,1.185386728843735,5.4851321695520836},
        {-2.8216586095746585,1.8748037193174842,-3.7739108983940559},
        {-2.3684374806144071,3.2075926659486464,3.3303063067244341},
        {0.093124111828116238,-1.35181399260482,-4.9053598052877607},
        {4.4406543872673385,2.7075640683000683,0.4366000483306508},
    };

    float SecondLayerWeight[4][2] = {
        {1.5389008605237542,5.3844216437162515},
        {0.79857647496162087,-9.4352753275067478},
        {-1.7449580935377649,-5.3580107251214919},
        {4.3899329977858814,9.4021829034453575},
    };

    float SecondLayerBias [4][1] = {
        {-6.2614696535343883},
        {-0.019643179146663855},
        {2.2987927416795513},
        {0.19030960925498394},
    };

    float ThirdLayerWeight[2][1] = {
        {-5.848781527524336},
        {-5.6209710811412403},
    };

    float ThirdLayerBias [2][1] = {
        {4.994739278766394},
        {-0.68349310166014743},
    };

    float FourthLayerWeight1[6][6] = {
        {0.44877903382792333, -0.18254257814405594, -0.5473500863373173, -0.5169833821838733, -0.51032491282695858, -0.16942474507699373},
        {-0.18790036759284953, -0.26304777707464544,-0.33261319704729081,1.1576434730697813,-1.2867431339504793,-0.077830787536354085},
        {0.3364691773905607,-0.23529779106560525,0.23479515988057528,-0.52679632405952947,-0.32193897061867188,0.99983912186290358},
        {-0.1109396666084288,0.12589761050901063,0.39768166255874171,-0.95504719167885643,-0.48378034629670724,-0.1460218192386665},
        {-0.56781929331280412,-0.1340982754278453,-0.27227312962879308,0.65829463376767638,-0.60128020573007512,-0.20450478250200113},
        {-0.33539837834969327,2.3587877179318752,-0.18965700642344413,0.53897717785022026,-0.010261133397695224,2.0066601309892187},
      };

    float FourthLayerWeight2[6][4] = {
        {0.85322762061547153,0.38151487077552471,0.013375224604868622,0.47030099532818342},
        {0.12484465422561394,0.032015171088060207,-0.9782598080355368,1.1911283432788269},
        {0.86672236918271872,-7.537237294147964,-0.66984774984780726,8.3340732800770603},
        {0.018358440088093841,-1.8631112090166739,-0.41297818840840744,0.98408622321802774},
        {0.58962219530808313,-0.6998245214671831,-0.27655082884087873,0.59878151983249128},
        {0.335063480868177,0.51247635784616929,0.11022166259630234,-2.2748227165955215},
    };

    float FourthLayerWeight3[6][2] = {
        {-0.023965705562068458,-0.25539237922783004},
        {0.53012265293183103,0.036057538999157722},
        {0.01683904990344131,-0.16364599191398779},
        {0.079757626110007598,-0.26961832779387657},
        {0.12382350890396553,0.17415719333320617},
        {-0.19597829571397554,0.58209152411680332},

    };

    float FourthLayerBias[6][1] = {
        {-1.8339298975128442},
        {-1.5816808464939591},
        {-0.77466247684221412},
        {-0.11298923019448159},
        {-0.95466593701455416},
        {-1.7755738768746663},
};

    float FifthLayerWeight [1][6] = {
        {1.6857476163957275,3.395168187570675,6.6872284910476232,3.6176310135802896,1.3334351611372237,3.8604197851257109},
    };

    float FifthLayerBias = -2.6538309721236089;
    // LAYER ONE
    float FirstLayerOut [6][1];
    int i;
    float n;

    for (i = 0; i < 6; i++) {
        n = I1[0][0] * FirstLayerWeight[i][0] + I1[1][0] * FirstLayerWeight[i][1] + I1[2][0] * FirstLayerWeight[i][2] + FirstLayerBias [i][0];
        FirstLayerOut[i][0] = 2 / (1 + exp(-2 * n)) - 1;
    }
    // LAYER TWO
    float SecondLayerOut[4][1];

    for (i = 0; i < 4; i++) {
        n = I2[0][0] * SecondLayerWeight[i][0] + I2[1][0] * SecondLayerWeight[i][1] + SecondLayerBias [i][0];
        SecondLayerOut[i][0] = 2 / (1 + exp(-2 * n)) - 1;
    }
    // LAYER THREE
    float ThirdLayerOut[2][1];

    for (i = 0; i < 2; i++) {
        n = I3 * ThirdLayerWeight[i][0] + ThirdLayerBias [i][0];
        ThirdLayerOut[i][0] = 2 / (1 + exp(-2 * n)) - 1;
    }
    // LAYER FOUR //
    float FourthLayerOut[6][1];

    for (i = 0; i < 6; i++) {
        n = FirstLayerOut[0][0] * FourthLayerWeight1[i][0] + FirstLayerOut[1][0] * FourthLayerWeight1[i][1] + FirstLayerOut[2][0] * FourthLayerWeight1[i][2] + FirstLayerOut[3][0] * FourthLayerWeight1[i][3] + FirstLayerOut[4][0] * FourthLayerWeight1[i][4] + FirstLayerOut[5][0] * FourthLayerWeight1[i][5]
                + SecondLayerOut[0][0] * FourthLayerWeight2[i][0] + SecondLayerOut[1][0] * FourthLayerWeight2[i][1] + SecondLayerOut[2][0] * FourthLayerWeight2[i][2] + SecondLayerOut[3][0] * FourthLayerWeight2[i][3]
                + ThirdLayerOut[0][0] * FourthLayerWeight3[i][0] + ThirdLayerOut[1][0] * FourthLayerWeight3[i][1]
                + FourthLayerBias [i][0];
        FourthLayerOut[i][0] = 2 / (1 + exp(-2 * n)) - 1;
    }
    // LAYER FIVE
    float FifthLayerOut;

    n = FourthLayerOut[0][0] * FifthLayerWeight[0][0] + FourthLayerOut[1][0] * FifthLayerWeight[0][1] + FourthLayerOut[2][0] * FifthLayerWeight[0][2] + FourthLayerOut[3][0] * FifthLayerWeight[0][3] + FourthLayerOut[4][0] * FifthLayerWeight[0][4] + FourthLayerOut[5][0] * FifthLayerWeight[0][5]
            + FifthLayerBias;
    FifthLayerOut = 1 / (1 + exp(-n));

    float ANNOUTPUT;
    ANNOUTPUT = FifthLayerOut;

    return (ANNOUTPUT);
}
 
static void _MDReservoirNeuralNet(int itemID) {

    float discharge; // Current discharge [m3/s]
    float resCapacity; // Reservoir capacity [m3]
    float minresStorage;
    float discharge_t_1;
    float discharge_t_2;
    float discharge_t_3;

    float I1[3][1]; // Input to ANN (ANNOUTPUT.c)
    float I2[2][1]; // Input to ANN (ANNOUTPUT.c)
    float I3; // Input to ANN (ANNOUTPUT.c)
    float ANN;
    float resStorage; // Reservoir storage [m3]
    float resStorageChg; // Reservoir storage change [m3/dt]
    float resRelease; // Reservoir release [m3/s]
    float res_release_t_1; // STANDARD Reservoir release (resRelease) [m3/s]
    float res_release_t_2;
    float res_release_t_3;
    float SIMOUT;
    float prevResStorage;
    float SD_t_3;
    float SD_t_2;
    float SD_t_1;
    float SR_t_3;
    float SR_t_2;
    int y = MFDateGetCurrentYear();

    discharge = MFVarGetFloat(_MDInDischargeID, itemID, 0.0);
   
    if (((resCapacity = MFVarGetFloat(_MDInResCapacityID, itemID, 0.0)) <= 0.0) || y<=1900 ){
                       MFVarSetFloat(_MDOutResStorageID, itemID, 0.0);
                       MFVarSetFloat(_MDOutResStorageChgID, itemID, 0.0);
                       MFVarSetFloat(_MDOutResReleaseID, itemID, discharge);
        return;
    }
        resCapacity = 0.90*MFVarGetFloat(_MDInResCapacityID, itemID, 0.0); //Assuming 30% dead storage; Should add 0.25*MaxCap to final result while presenting
        prevResStorage = MFVarGetFloat(_MDOutPreResStorageID, itemID, 0.03*resCapacity);
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
* !!!!!! 0.3*resCapacity should be added to Reservoir storage result file !!!!!!!!
* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/


        discharge_t_1 = discharge;
        discharge_t_2 = MFVarGetFloat(_MDOutDisch_t_2_ID, itemID, 0); // Last Month
        discharge_t_3 = MFVarGetFloat(_MDOutDisch_t_3_ID, itemID, 0); // Two Month Ago
        res_release_t_2 = MFVarGetFloat(_MDOutResRelease_t_2_ID, itemID, 0); // Last Month
        res_release_t_3 = MFVarGetFloat(_MDOutResRelease_t_3_ID, itemID, 0); // Two Month Ago

        SD_t_1 = discharge_t_1 *24*3600/resCapacity;
        SD_t_2 = discharge_t_2 *24*3600/resCapacity;
        SD_t_3 = discharge_t_3 *24*3600/resCapacity;
        
        SR_t_3 = res_release_t_3 *24*3600/resCapacity;
        SR_t_2 = res_release_t_2 *24*3600/resCapacity;
        
        I1[0][0] = SD_t_3;
        I1[1][0] = SD_t_2;
        I1[2][0] = SD_t_1;

        I2[0][0] = SR_t_3;
        I2[1][0] = SR_t_2;
        

        I3 = prevResStorage/resCapacity;

 
        ANN = ANNOUTPUT (I1, I2, I3)*resCapacity/(24*3600) ;
        


        

        resStorageChg = (discharge - ANN)*3600 * 24;
        minresStorage = 0;
        

        if ((prevResStorage + resStorageChg <= resCapacity) && (prevResStorage + resStorageChg >= minresStorage)) {
            SIMOUT = ANN;
             
            if (SIMOUT < 0) {
                printf("Error: Negative release (1)! \n");
                printf("%f %f %f %f %f \n", SIMOUT, resCapacity, resStorageChg, minresStorage, discharge);
            }
            resStorage = prevResStorage + resStorageChg;
        } else {
            if ((prevResStorage + resStorageChg) > resCapacity) {
                SIMOUT = ((discharge * 3600 * 24)-(resCapacity - prevResStorage)) / (3600 * 24);
                
                if (SIMOUT < 0) {
                    printf("Error: Negative release (2)! \n");
                        printf("%f %f %f %f %f \n", SIMOUT, resCapacity, resStorageChg, minresStorage, discharge); }
                resStorage = resCapacity;
            } else {
                SIMOUT = (prevResStorage - minresStorage + (discharge * 3600 * 24)) / (3600 * 24);
               
                if (SIMOUT < 0) {
                    printf("Error: Negative release (3)! \n");
                        printf("%f %f %f %f %f \n", SIMOUT, resCapacity, resStorageChg, minresStorage, discharge); }
                resStorage = minresStorage;
            }
        }
                // printf("%f %f %f %f %f %f %f\n",SIMOUT, ANN, discharge, resCapacity, resStorageChg, minresStorage, resStorage);
     
///////////////////////////////////////////////////////////////////
/////////////////// Maximum Discharge ///////////////////////
/* if (SIMOUT > (resCapacity/(3.21*24*3600))){
resRelease=(resCapacity/(3.21*24*3600));
resStorage = prevResStorage + (discharge * 3600 * 24)-(resRelease * 3600 * 24);
if (resStorage > resCapacity) {
resRelease = ((discharge * 3600 * 24)-(resCapacity - prevResStorage)) / (3600 * 24);
resStorage = resCapacity;
}
} else {
resRelease = SIMOUT;
}
*/
////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////

                      resRelease = SIMOUT;

        res_release_t_1 = resRelease;

        MFVarSetFloat(_MDOutResReleaseID, itemID, resRelease);
        MFVarSetFloat(_MDOutDisch_t_3_ID, itemID, discharge_t_2); //You Should Set t-2 before t-1
        MFVarSetFloat(_MDOutDisch_t_2_ID, itemID, discharge_t_1); //You Should Set t-2 before t-1
        MFVarSetFloat(_MDOutResRelease_t_3_ID, itemID, res_release_t_2); //You Should Set t-2 before t-1
        MFVarSetFloat(_MDOutResRelease_t_2_ID, itemID, res_release_t_1); //You Should Set t-2 before t-1
        MFVarSetFloat(_MDOutResStorageChgID, itemID, resStorageChg);
        MFVarSetFloat(_MDOutResStorageID, itemID, resStorage);
        MFVarSetFloat(_MDOutPreResStorageID, itemID, resStorage);
}

static void _MDReservoirDW(int itemID) {

    // Input
    float discharge; // Current discharge [m3/s]
    float meanDischarge; // Long-term mean annual discharge [m3/s]
    float resCapacity; // Reservoir capacity [km3]

    // Output
    float resStorage; // Reservoir storage [km3]
    float resStorageChg; // Reservoir storage change [km3/dt]
    float resRelease; // Reservoir release [m3/s]

    // local
    float prevResStorage; // Reservoir storage from the previous time step [km3]
    float dt; // Time step length [s]
    float balance; // water balance [m3/s]

    // Parameters
    float drySeasonPct = 0.60; // RJS 071511
    float wetSeasonPct = 0.16; // RJS 071511
    float year = 0; // RJS 082311

    discharge = MFVarGetFloat(_MDInDischargeID, itemID, 0.0);
    //printf("discharge= %f \n", discharge);
    meanDischarge = MFVarGetFloat(_MDInDischMeanID, itemID, discharge);
    year = MFDateGetCurrentYear();

    if ((resCapacity = MFVarGetFloat(_MDInResCapacityID, itemID, 0.0)) <= 0.0) {
        MFVarSetFloat(_MDOutResStorageID, itemID, 0.0);
        MFVarSetFloat(_MDOutResStorageChgID, itemID, 0.0);
        MFVarSetFloat(_MDOutResReleaseID, itemID, discharge);
        printf("discharge NO Res= %f \n" , discharge);
        return;
    }
    printf("discharge WT Res= %f \n" , discharge);
    dt = MFModelGet_dt();
    prevResStorage = MFVarGetFloat(_MDOutResStorageID, itemID, 0.0);

    resRelease = discharge > meanDischarge ? wetSeasonPct * discharge : drySeasonPct * discharge + (meanDischarge - discharge);

    resStorage = prevResStorage + (discharge - resRelease) * 86400.0 / 1e9;

    if (resStorage > resCapacity) {
        resRelease = discharge * dt / 1e9 + prevResStorage - resCapacity;
        resRelease = resRelease * 1e9 / dt;
        resStorage = resCapacity;
    } else if (resStorage < 0.0) {
        resRelease = prevResStorage + discharge * dt / 1e9;
        resRelease = resRelease * 1e9 / dt;
        resStorage = 0;
    }

    resStorageChg = resStorage - prevResStorage;

    balance = discharge - resRelease - (resStorageChg / 86400 * 1e9); // water balance

    MFVarSetFloat(_MDOutResStorageID, itemID, resStorage);
    MFVarSetFloat(_MDOutResStorageChgID, itemID, resStorageChg);
    MFVarSetFloat(_MDOutResReleaseID, itemID, resRelease);
}

enum {
    MDnone, MDcalculate, MDneuralnet
};

int MDReservoirDef() {
    int optID = MFUnset;
    const char *optStr, *optName = MDOptReservoirs;
    const char *options [] = {MDNoneStr, MDCalculateStr, "neuralnet", (char *) NULL};

    if ((optStr = MFOptionGet(optName)) != (char *) NULL) optID = CMoptLookup(options, optStr, true);

    if ((optID == MDnone) || (_MDOutResReleaseID != MFUnset)) return (_MDOutResReleaseID);

    MFDefEntering("Reservoirs");
    switch (optID) {
        case MDcalculate:
            if ( ((_MDInDischMeanID = MDDischMeanDef()) == CMfailed) ||
                    ((_MDInDischargeID = MDDischLevel2Def()) == CMfailed) ||
                    ((_MDInResCapacityID = MFVarGetID(MDVarReservoirCapacity, "km3", MFInput, MFState, MFBoundary)) == CMfailed) ||
                    ((_MDOutResStorageID = MFVarGetID(MDVarReservoirStorage, "km3", MFOutput, MFState, MFInitial)) == CMfailed) ||
                    ((_MDOutResStorageChgID = MFVarGetID(MDVarReservoirStorageChange, "km3", MFOutput, MFState, MFBoundary)) == CMfailed) ||
                    ((_MDOutResReleaseID = MFVarGetID(MDVarReservoirRelease, "m3/s", MFOutput, MFFlux, MFBoundary)) == CMfailed) ||
                     ( MDHydroPowerDef() == CMfailed) ||
                    (MFModelAddFunction(_MDReservoirDW) == CMfailed)) return (CMfailed);
            break;
        case MDneuralnet:

            if ( ((_MDInDischargeID = MDDischLevel2Def()) == CMfailed) ||
                    ((_MDInResCapacityID = MFVarGetID(MDVarReservoirCapacity, "m3", MFInput, MFState, MFBoundary)) == CMfailed) ||
                    ((_MDOutDisch_t_1_ID = MFVarGetID(MDVarDisch_t_1_, "m3/s", MFOutput, MFState, MFInitial)) == CMfailed) ||
                    ((_MDOutDisch_t_2_ID = MFVarGetID(MDVarDisch_t_2_, "m3/s", MFOutput, MFState, MFInitial)) == CMfailed) ||
                    ((_MDOutDisch_t_3_ID = MFVarGetID(MDVarDisch_t_3_, "m3/s", MFOutput, MFState, MFInitial)) == CMfailed) ||
                    ((_MDOutResStorageID = MFVarGetID(MDVarReservoirStorage, "m3" , MFOutput, MFState, MFInitial)) == CMfailed) ||
                    ((_MDOutPreResStorageID = MFVarGetID(MDVarPreResStorage, "m3" , MFOutput, MFState, MFInitial)) == CMfailed) ||
                    ((_MDOutResStorageChgID = MFVarGetID(MDVarReservoirStorageChange, "m3" , MFOutput, MFState, MFInitial)) == CMfailed) ||
                    ((_MDOutResReleaseID = MFVarGetID(MDVarReservoirRelease, "m3/s", MFOutput, MFState, MFInitial)) == CMfailed) ||
                    ((_MDOutResRelease_t_1_ID = MFVarGetID(MDVarResRelease_t_1_, "m3/s", MFOutput, MFState, MFInitial)) == CMfailed) ||
                    ((_MDOutResRelease_t_2_ID = MFVarGetID(MDVarResRelease_t_2_, "m3/s", MFOutput, MFState, MFInitial)) == CMfailed) ||
                    ((_MDOutResRelease_t_3_ID = MFVarGetID(MDVarResRelease_t_3_, "m3/s", MFOutput, MFState, MFInitial)) == CMfailed) ||
                    ( MDHydroPowerDef() == CMfailed) ||
                    (MFModelAddFunction(_MDReservoirNeuralNet) == CMfailed)) return (CMfailed);
              break;
        default: MFOptionMessage(optName, optStr, options);
            return (CMfailed);
    }
    MFDefLeaving("Reservoirs");
    return (_MDOutResReleaseID);
}

