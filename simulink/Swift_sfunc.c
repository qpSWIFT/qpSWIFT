/* This file creates simulink interface for qpSWIFT solver
Author - Abhishek Pandala
*/

#include "qpSWIFT.h"

#define S_FUNCTION_NAME Swift_sfunc /**< Name of the S function. */
#define S_FUNCTION_LEVEL 2          /**< S function level. */

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
    if (!ssSetNumInputPorts(S, 5))
        return;

    /* Specify the number of output ports */
    if (!ssSetNumOutputPorts(S, 5))
        return;

    /* Specify dimension information for the input ports */

    ssSetInputPortVectorDimension(S, 0, DYNAMICALLY_SIZED); /* P - Matrix */
    ssSetInputPortVectorDimension(S, 1, DYNAMICALLY_SIZED); /* c - Vector */
    ssSetInputPortVectorDimension(S, 2, DYNAMICALLY_SIZED); /* G - Matrix */
    ssSetInputPortVectorDimension(S, 3, DYNAMICALLY_SIZED); /* h - Vector */
    ssSetInputPortVectorDimension(S, 4, DYNAMICALLY_SIZED); /* sigma-d */

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

    real *c, *h;
    idxint *Pir, *Pjc, *Gir, *Gjc;
    real *Ppr, *Gpr;
    QP *myQP;
    idxint m, n;
    real sigma_d;
    idxint i, l;
    idxint a1, b1, c1, d1, e1;

    real_T *out_x, *out_exitflag, *out_iter, *out_tsetup, *out_tsolve;

    out_x = ssGetOutputPortRealSignal(S, 0);        /* x */
    out_exitflag = ssGetOutputPortRealSignal(S, 1); /* exitflag */
    out_iter = ssGetOutputPortRealSignal(S, 2);     /* iter */
    out_tsetup = ssGetOutputPortRealSignal(S, 3);   /* tsetup */
    out_tsolve = ssGetOutputPortRealSignal(S, 4);   /* tsolve */

    InputRealPtrsType in_P, in_c, in_G, in_h, in_sigma_d;

    in_P = ssGetInputPortRealSignalPtrs(S, 0);
    in_c = ssGetInputPortRealSignalPtrs(S, 1);
    in_G = ssGetInputPortRealSignalPtrs(S, 2);
    in_h = ssGetInputPortRealSignalPtrs(S, 3);
    in_sigma_d = ssGetInputPortRealSignalPtrs(S, 4);

    idxint Permut[504] = {
        78, 81, 84, 80, 79, 86, 85, 83, 82, 87, 88, 89, 117, 108, 111, 114, 110, 109, 116, 115, 113, 112, 118, 119, 147, 138, 141, 144, 140, 139, 146, 145, 143, 142, 148, 149, 195, 192, 318, 165, 300, 135, 162, 302, 137, 197, 194, 320, 167, 207, 208, 209, 481, 480, 180, 497, 496, 483, 482, 181, 485, 484, 183, 499, 498, 487, 486, 184, 489, 488, 186, 501, 500, 491, 490, 187, 493, 492, 189, 503, 502, 495, 494, 190, 191, 188, 185, 182, 323, 164, 321, 457, 456, 150, 473, 472, 459, 458, 151, 461, 460, 153, 475, 474, 463, 462, 154, 465, 464, 156, 477, 476, 467, 466, 157, 469, 468, 159, 479, 478, 471, 470, 160, 161, 158, 155, 152, 196, 193, 319, 166, 301, 136, 163, 304, 305, 322, 303, 198, 179, 201, 204, 200, 199, 168, 171, 174, 170, 169, 206, 205, 203, 202, 177, 178, 176, 175, 173, 172, 330, 332, 331, 324, 325, 326, 327, 329, 328, 333, 334, 335, 307, 308, 306, 309, 310, 311, 312, 313, 314, 315, 316, 317, 288, 290, 289, 291, 292, 293, 295, 296, 294, 433, 432, 120, 449, 448, 435, 434, 121, 437, 436, 123, 451, 450, 439, 438, 124, 441, 440, 126, 453, 452, 443, 442, 127, 445, 444, 129, 455, 454, 447, 446, 130, 131, 128, 125, 122, 266, 77, 284, 107, 104, 409, 408, 90, 425, 424, 411, 410, 91, 413, 412, 93, 427, 426, 415, 414, 94, 417, 416, 96, 429, 428, 419, 418, 97, 421, 420, 99, 431, 430, 423, 422, 100, 101, 98, 95, 92, 264, 75, 282, 105, 102, 267, 269, 385, 384, 60, 401, 400, 387, 386, 61, 389, 388, 63, 403, 402, 391, 390, 64, 393, 392, 66, 405, 404, 395, 394, 67, 397, 396, 69, 407, 406, 399, 398, 70, 71, 68, 65, 62, 247, 46, 265, 76, 283, 106, 103, 73, 250, 268, 211, 229, 16, 13, 337, 336, 0, 353, 352, 339, 338, 1, 341, 340, 3, 355, 354, 343, 342, 4, 345, 344, 6, 357, 356, 347, 346, 7, 349, 348, 9, 359, 358, 351, 350, 10, 11, 8, 5, 2, 225, 227, 226, 214, 248, 47, 212, 230, 17, 14, 44, 210, 228, 15, 12, 246, 45, 361, 360, 30, 377, 376, 363, 362, 31, 365, 364, 33, 379, 378, 367, 366, 34, 369, 368, 36, 381, 380, 371, 370, 37, 373, 372, 39, 383, 382, 375, 374, 40, 41, 38, 35, 32, 231, 42, 233, 213, 215, 48, 57, 58, 59, 220, 221, 223, 224, 222, 217, 218, 219, 216, 18, 21, 20, 24, 19, 23, 25, 26, 22, 242, 241, 240, 239, 238, 237, 236, 235, 50, 51, 52, 53, 54, 55, 56, 234, 49, 27, 28, 243, 244, 245, 29, 43, 72, 232, 249, 251, 74, 132, 133, 134, 252, 253, 254, 255, 256, 257, 258, 259, 260, 261, 262, 263, 270, 271, 272, 273, 274, 275, 276, 277, 278, 279, 280, 281, 285, 286, 287, 298, 299, 297};

    a1 = ssGetInputPortWidth(S, 0); // P
    b1 = ssGetInputPortWidth(S, 1); // c
    c1 = ssGetInputPortWidth(S, 2); // G
    d1 = ssGetInputPortWidth(S, 3); // h
    e1 = ssGetInputPortWidth(S, 4); // sigma_d

    sigma_d = (real)(*in_sigma_d)[0];
    n = (idxint)b1;
    m = (idxint)d1;

    c = (real *)malloc(n * sizeof(real));
    h = (real *)malloc(m * sizeof(real));

    for (i = 0; i < n; i++)
        c[i] = (*in_c)[i];

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

    for (i = 0; i < c1; i++)
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

    myQP = QP_SETUP(n, m, 0, Pjc, Pir, Ppr, NULL, NULL, NULL, Gjc, Gir, Gpr, c, h, NULL, sigma_d, Permut);
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
            out_x[i] = (real_T)(0);
        }
        out_exitflag[0] = (real_T)(1);
        out_iter[0] = (real_T)(2);
        out_tsetup[0] = (real_T)(3);
        out_tsolve[0] = (real_T)(4);
        QP_CLEANUP(myQP);
    }

    free(Pjc);
    free(Pir);
    free(Ppr);
    free(Gpr);
    free(Gir);
    free(Gjc);
    free(c);
    free(h);
}

static void mdlTerminate(SimStruct *S)
{
}

#ifdef MATLAB_MEX_FILE
#include "simulink.c"
#else
#include "cg_sfun.h"
#endif
