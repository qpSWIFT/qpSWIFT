#ifndef __TEST_DATA_STRUCTURES__HPP
#define __TEST_DATA_STRUCTURES__HPP
//#include "Prime.h"
#define qp_real double
#define qp_int long


typedef struct data
{
    qp_int m;
    qp_int n;
    qp_int p;

    qp_int P_nnz;
    qp_int A_nnz;
    qp_int G_nnz;

    /* Data Pointers */

    qp_real *Ppr;
    qp_int *Pir;
    qp_int *Pjc;

    qp_real *Apr;
    qp_int *Air;
    qp_int *Ajc;

    qp_real *Gpr;
    qp_int *Gir;
    qp_int *Gjc;

    qp_real *s;
    qp_real *delta_s;
    qp_real *z;
    qp_real *delta_z;

    qp_int *Permut;

} standalone_data;



#endif