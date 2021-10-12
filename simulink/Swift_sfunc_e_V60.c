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

    idxint Permut[2] = {1, 2};
    idxint Permut1[2] = {1, 2};
    idxint Permut2[2] = {1, 2};
    idxint Permut3[2] = {1, 2};

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

    if (sigma_d == 0)
    {
        myQP = QP_SETUP(n, m, p, Pjc, Pir, Ppr, Ajc, Air, Apr, Gjc, Gir, Gpr, c, h, b, 0, Permut);
    }
    else if (sigma_d == 1)
    {
        myQP = QP_SETUP(n, m, p, Pjc, Pir, Ppr, Ajc, Air, Apr, Gjc, Gir, Gpr, c, h, b, 0, Permut1);
    }
    else if (sigma_d == 2)
    {
        myQP = QP_SETUP(n, m, p, Pjc, Pir, Ppr, Ajc, Air, Apr, Gjc, Gir, Gpr, c, h, b, 0, Permut2);
    }
    else if (sigma_d == 3)
    {
        myQP = QP_SETUP(n, m, p, Pjc, Pir, Ppr, Ajc, Air, Apr, Gjc, Gir, Gpr, c, h, b, 0, Permut3);
    }
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
