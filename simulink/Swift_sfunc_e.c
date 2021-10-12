/* This file creates simulink interface for qpSWIFT solver
Author - Abhishek Pandala
*/

#include "qpSWIFT.h"

#define S_FUNCTION_NAME Swift_sfunc_e /**< Name of the S function. */
#define S_FUNCTION_LEVEL 2            /**< S function level. */

#define SAMPLINGTIME -1 /**< Sampling time. */
#include "simstruc.h"

#define NV 12 /*Output Needed */

static void mdlInitializeSizes(SimStruct *S) /* Init sizes array */
{

    /* Specify the number of continuous and discrete states */
    ssSetNumContStates(S, 0);
    ssSetNumDiscStates(S, 0);

    /* Specify the number of parameters */
    ssSetNumSFcnParams(S, 0);
    if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S))
        return;

    /* Specify the number of input ports */
    if (!ssSetNumInputPorts(S, 7))
        return;

    /* Specify the number of output ports */
    if (!ssSetNumOutputPorts(S, 5))
        return;

    /* Specify dimension information for the input ports */

    ssSetInputPortVectorDimension(S, 0, DYNAMICALLY_SIZED); /* P - Matrix */
    ssSetInputPortVectorDimension(S, 1, DYNAMICALLY_SIZED); /* c - Vector */
    ssSetInputPortVectorDimension(S, 2, DYNAMICALLY_SIZED); /* A - Matrix */
    ssSetInputPortVectorDimension(S, 3, DYNAMICALLY_SIZED); /* b - Vector */
    ssSetInputPortVectorDimension(S, 4, DYNAMICALLY_SIZED); /* G - Matrix */
    ssSetInputPortVectorDimension(S, 5, DYNAMICALLY_SIZED); /* h - Vector */
    ssSetInputPortVectorDimension(S, 6, DYNAMICALLY_SIZED); /* sigma-d */

    /* Specify dimension information for the output ports */

    ssSetOutputPortVectorDimension(S, 0, NV); /* x */
    ssSetOutputPortVectorDimension(S, 1, 1);  /* exitflag */
    ssSetOutputPortVectorDimension(S, 2, 1);  /* iter */
    ssSetOutputPortVectorDimension(S, 3, 1);  /* tsetup */
    ssSetOutputPortVectorDimension(S, 4, 1);  /* tsolve */

    /* Specify the direct feedthrough status */

    ssSetInputPortDirectFeedThrough(S, 0, 1);
    ssSetInputPortDirectFeedThrough(S, 1, 1);
    ssSetInputPortDirectFeedThrough(S, 2, 1);
    ssSetInputPortDirectFeedThrough(S, 3, 1);
    ssSetInputPortDirectFeedThrough(S, 4, 1);
    ssSetInputPortDirectFeedThrough(S, 5, 1);
    ssSetInputPortDirectFeedThrough(S, 6, 1);

    /* One sample time */

    ssSetNumSampleTimes(S, 1);
}

#if defined(MATLAB_MEX_FILE)

#define MDL_SET_INPUT_PORT_DIMENSION_INFO
#define MDL_SET_OUTPUT_PORT_DIMENSION_INFO

static void mdlSetInputPortDimensionInfo(SimStruct *S, int_T port, const DimsInfo_T *dimsInfo)
{
    if (!ssSetInputPortDimensionInfo(S, port, dimsInfo))
        return;
}

static void mdlSetOutputPortDimensionInfo(SimStruct *S, int_T port, const DimsInfo_T *dimsInfo)
{
    if (!ssSetOutputPortDimensionInfo(S, port, dimsInfo))
        return;
}

#endif

static void mdlInitializeSampleTimes(SimStruct *S)
{
    ssSetSampleTime(S, 0, SAMPLINGTIME);
    ssSetOffsetTime(S, 0, 0.0);
}

static void mdlOutputs(SimStruct *S, int_T tid)
{

    real *c, *h, *b;
    idxint *Pir, *Pjc, *Gir, *Gjc, *Air, *Ajc;
    real *Ppr, *Gpr, *Apr;
    QP *myQP;
    idxint m, n, p;
    real sigma_d;
    idxint i, l;
    idxint a1, b1, c1, d1, e1, f1, g1;

    real_T *out_x, *out_exitflag, *out_iter, *out_tsetup, *out_tsolve;

    out_x = ssGetOutputPortRealSignal(S, 0);        /* x */
    out_exitflag = ssGetOutputPortRealSignal(S, 1); /* exitflag */
    out_iter = ssGetOutputPortRealSignal(S, 2);     /* iter */
    out_tsetup = ssGetOutputPortRealSignal(S, 3);   /* tsetup */
    out_tsolve = ssGetOutputPortRealSignal(S, 4);   /* tsolve */

    InputRealPtrsType in_P, in_c, in_A, in_b, in_G, in_h, in_sigma_d;

    in_P = ssGetInputPortRealSignalPtrs(S, 0);
    in_c = ssGetInputPortRealSignalPtrs(S, 1);
    in_A = ssGetInputPortRealSignalPtrs(S, 2);
    in_b = ssGetInputPortRealSignalPtrs(S, 3);
    in_G = ssGetInputPortRealSignalPtrs(S, 4);
    in_h = ssGetInputPortRealSignalPtrs(S, 5);
    in_sigma_d = ssGetInputPortRealSignalPtrs(S, 6);

    idxint Permut[504] = {
        224,
        223,
        222,
        221,
        220,
        219,
        218,
        217,
        216,
        240,
        241,
        242,
        234,
        235,
        236,
        237,
        238,
        239,
        59,
        58,
        57,
        260,
        259,
        258,
        89,
        88,
        87,
        80,
        79,
        78,
        278,
        277,
        276,
        83,
        82,
        81,
        119,
        118,
        117,
        110,
        109,
        108,
        296,
        295,
        294,
        113,
        112,
        111,
        448,
        449,
        434,
        435,
        121,
        432,
        433,
        120,
        122,
        450,
        451,
        438,
        439,
        124,
        436,
        437,
        123,
        125,
        452,
        453,
        442,
        443,
        127,
        440,
        441,
        126,
        128,
        454,
        455,
        446,
        447,
        130,
        444,
        445,
        129,
        131,
        106,
        283,
        197,
        196,
        193,
        319,
        166,
        136,
        301,
        163,
        135,
        300,
        195,
        192,
        318,
        165,
        496,
        497,
        482,
        483,
        181,
        480,
        481,
        180,
        182,
        498,
        499,
        486,
        487,
        184,
        484,
        485,
        183,
        185,
        500,
        501,
        490,
        491,
        187,
        488,
        489,
        186,
        188,
        502,
        503,
        494,
        495,
        190,
        492,
        493,
        189,
        191,
        321,
        162,
        322,
        472,
        473,
        458,
        459,
        151,
        456,
        457,
        150,
        152,
        474,
        475,
        462,
        463,
        154,
        460,
        461,
        153,
        155,
        476,
        477,
        466,
        467,
        157,
        464,
        465,
        156,
        158,
        478,
        479,
        470,
        471,
        160,
        468,
        469,
        159,
        161,
        194,
        320,
        167,
        137,
        302,
        164,
        305,
        303,
        304,
        323,
        133,
        286,
        105,
        282,
        75,
        264,
        45,
        246,
        72,
        102,
        107,
        284,
        424,
        425,
        410,
        411,
        91,
        408,
        409,
        90,
        92,
        426,
        427,
        414,
        415,
        94,
        412,
        413,
        93,
        95,
        428,
        429,
        418,
        419,
        97,
        416,
        417,
        96,
        98,
        430,
        431,
        422,
        423,
        100,
        420,
        421,
        99,
        101,
        77,
        266,
        47,
        248,
        74,
        269,
        104,
        267,
        132,
        285,
        287,
        134,
        376,
        377,
        362,
        363,
        31,
        360,
        361,
        30,
        32,
        378,
        379,
        366,
        367,
        34,
        364,
        365,
        33,
        35,
        380,
        381,
        370,
        371,
        37,
        368,
        369,
        36,
        38,
        382,
        383,
        374,
        375,
        40,
        372,
        373,
        39,
        41,
        211,
        16,
        229,
        13,
        210,
        15,
        228,
        12,
        352,
        353,
        338,
        339,
        1,
        336,
        337,
        0,
        2,
        354,
        355,
        342,
        343,
        4,
        340,
        341,
        3,
        5,
        356,
        357,
        346,
        347,
        7,
        344,
        345,
        6,
        8,
        358,
        359,
        350,
        351,
        10,
        348,
        349,
        9,
        11,
        227,
        226,
        225,
        212,
        17,
        230,
        14,
        215,
        213,
        214,
        400,
        401,
        386,
        387,
        61,
        384,
        385,
        60,
        62,
        402,
        403,
        390,
        391,
        64,
        388,
        389,
        63,
        65,
        404,
        405,
        394,
        395,
        67,
        392,
        393,
        66,
        68,
        406,
        407,
        398,
        399,
        70,
        396,
        397,
        69,
        71,
        76,
        265,
        46,
        247,
        73,
        250,
        43,
        232,
        42,
        44,
        103,
        231,
        233,
        249,
        268,
        251,
        149,
        148,
        179,
        207,
        208,
        209,
        206,
        205,
        204,
        203,
        202,
        201,
        200,
        199,
        198,
        178,
        177,
        170,
        330,
        331,
        332,
        324,
        325,
        326,
        327,
        328,
        329,
        169,
        168,
        173,
        172,
        171,
        176,
        175,
        174,
        147,
        140,
        139,
        138,
        312,
        313,
        314,
        306,
        307,
        308,
        309,
        310,
        311,
        143,
        142,
        141,
        315,
        316,
        317,
        334,
        335,
        333,
        288,
        289,
        291,
        292,
        293,
        298,
        299,
        290,
        144,
        145,
        297,
        146,
        270,
        271,
        273,
        274,
        275,
        280,
        281,
        272,
        114,
        115,
        116,
        279,
        252,
        253,
        254,
        255,
        256,
        257,
        263,
        262,
        84,
        85,
        86,
        261,
        29,
        18,
        19,
        20,
        21,
        22,
        23,
        24,
        25,
        26,
        28,
        48,
        49,
        50,
        27,
        51,
        52,
        53,
        54,
        55,
        56,
        245,
        244,
        243};
    a1 = ssGetInputPortWidth(S, 0); // P
    b1 = ssGetInputPortWidth(S, 1); // c
    c1 = ssGetInputPortWidth(S, 2); // A
    d1 = ssGetInputPortWidth(S, 3); // b
    e1 = ssGetInputPortWidth(S, 4); // G
    f1 = ssGetInputPortWidth(S, 5); // h
    g1 = ssGetInputPortWidth(S, 6); // sigma_d

    sigma_d = (real)(*in_sigma_d)[0];
    n = (idxint)b1;
    p = (idxint)d1;
    m = (idxint)f1;

    c = (real *)malloc(n * sizeof(real));
    b = (real *)malloc(p * sizeof(real));
    h = (real *)malloc(m * sizeof(real));

    for (i = 0; i < n; i++)
        c[i] = (*in_c)[i];

    for (i = 0; i < p; i++)
        b[i] = (*in_b)[i];

    for (i = 0; i < m; i++)
        h[i] = (*in_h)[i];

    float sparsity = 0.1;
    int nzemax = (int)(sparsity * m * n);
    int jc_index = 0;

    int nze_r = 0;
    Gjc = (idxint *)malloc((n + 1) * sizeof(idxint));
    Gir = (idxint *)malloc(nzemax * sizeof(idxint));
    Gpr = (real *)malloc(nzemax * sizeof(real));
    Gjc[0] = 0;

    for (i = 0; i < e1; i++)
    {
        if (((*in_G)[i]) != 0)
        {
            nze_r++;
            jc_index++;
            if (nze_r > nzemax)
            {
                sparsity += 0.1;
                nzemax = (int)(sparsity * m * n);
                Gpr = (real *)realloc(Gpr, nzemax * sizeof(real));
                Gir = (idxint *)realloc(Gir, nzemax * sizeof(idxint));
            }
            Gpr[nze_r - 1] = (*in_G)[i];
            Gir[nze_r - 1] = i % m;
        }

        if ((i % m))
        {
            l = i / m;
            Gjc[l + 1] = jc_index;
        }
    }

    sparsity = 0.1;
    nzemax = (idxint)(sparsity * n * n);
    jc_index = 0;

    nze_r = 0;
    Pjc = (idxint *)malloc((n + 1) * sizeof(idxint));
    Pir = (idxint *)malloc(nzemax * sizeof(idxint));
    Ppr = (real *)malloc(nzemax * sizeof(real));
    Pjc[0] = 0;

    for (i = 0; i < a1; i++)
    {
        if (((*in_P)[i]) != 0)
        {
            nze_r++;
            jc_index++;
            if (nze_r > nzemax)
            {
                sparsity += 0.1;
                nzemax = (int)(sparsity * n * n);
                Ppr = (real *)realloc(Ppr, nzemax * sizeof(real));
                Pir = (idxint *)realloc(Pir, nzemax * sizeof(idxint));
            }
            Ppr[nze_r - 1] = (*in_P)[i];
            Pir[nze_r - 1] = i % n;
        }

        if ((i % n))
        {
            l = i / n;
            Pjc[l + 1] = jc_index;
        }
    }

    sparsity = 0.1;
    nzemax = (idxint)(sparsity * p * n);
    jc_index = 0;

    nze_r = 0;
    Ajc = (idxint *)malloc((n + 1) * sizeof(idxint));
    Air = (idxint *)malloc(nzemax * sizeof(idxint));
    Apr = (real *)malloc(nzemax * sizeof(real));
    Ajc[0] = 0;

    for (i = 0; i < c1; i++)
    {
        if (((*in_A)[i]) != 0)
        {
            nze_r++;
            jc_index++;
            if (nze_r > nzemax)
            {
                sparsity += 0.1;
                nzemax = (int)(sparsity * p * n);
                Apr = (real *)realloc(Apr, nzemax * sizeof(real));
                Air = (idxint *)realloc(Air, nzemax * sizeof(idxint));
            }
            Apr[nze_r - 1] = (*in_A)[i];
            Air[nze_r - 1] = i % p;
        }

        if ((i % p))
        {
            l = i / p;
            Ajc[l + 1] = jc_index;
        }
    }

    myQP = QP_SETUP(n, m, p, Pjc, Pir, Ppr, Ajc, Air, Apr, Gjc, Gir, Gpr, c, h, b, sigma_d, Permut);
    if (myQP != NULL)
    {

        out_exitflag[0] = (real_T)QP_SOLVE(myQP);
        out_iter[0] = (real_T)(myQP->stats->IterationCount);
        for (i = 0; i < NV; i++)
        {
            out_x[i] = (real_T)(myQP->x[i]);
        }
        out_tsetup[0] = (real_T)myQP->stats->tsetup;
        out_tsolve[0] = (real_T)myQP->stats->tsolve;
        QP_CLEANUP(myQP);
    }
    else
    {
        for (i = 0; i < NV; i++)
        {
            out_x[i] = (real_T)(b[i]);
        }
        out_exitflag[0] = (real_T)(n);
        out_iter[0] = (real_T)(m);
        out_tsetup[0] = (real_T)(p);
        out_tsolve[0] = (real_T)(sigma_d);
        QP_CLEANUP(myQP);
    }

    free(Pjc);
    free(Pir);
    free(Ppr);
    free(Gpr);
    free(Gir);
    free(Gjc);
    free(Air);
    free(Ajc);
    free(Apr);
    free(c);
    free(h);
    free(b);
}

static void mdlTerminate(SimStruct *S)
{
}

#ifdef MATLAB_MEX_FILE
#include "simulink.c"
#else
#include "cg_sfun.h"
#endif
