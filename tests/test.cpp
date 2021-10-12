#include "test.hpp"
#define GTEST_COUT std::cerr << "[          ] [ INFO ]"

standAloneTests ::standAloneTests(/* args */)
{
    parser psr;

    if (!psr.readfile("some.txt"))
    {
        printf("Data cannot be read\n");
        // return 0;
    }
    else
    {
        GTEST_COUT << "Data read\n";
        data = psr.getdateptr();

        SparseMatrixSetup(data->n, data->n, data->P_nnz, data->Pjc, data->Pir, data->Ppr, &P);
        SparseMatrixSetup(data->p, data->n, data->A_nnz, data->Ajc, data->Air, data->Apr, &A);
        SparseMatrixSetup(data->m, data->n, data->G_nnz, data->Gjc, data->Gir, data->Gpr, &G);
    }
}

standAloneTests ::~standAloneTests()
{
}
// void standAloneTests::SetUp(standalone_data *new_data)
// {
//     data = new_data;
//     SparseMatrixSetup(data->n, data->n, data->P_nnz, data->Pjc, data->Pir, data->Ppr, &P);
//     SparseMatrixSetup(data->p, data->n, data->A_nnz, data->Ajc, data->Air, data->Apr, &A);
//     SparseMatrixSetup(data->m, data->n, data->G_nnz, data->Gjc, data->Gir, data->Gpr, &G);
// }

// TEST_F(standAloneTests, func_obj_value)
// {

//     qp_real c[2] = {1.0, 1.0};

//     qp_real temp[2];

//     qp_real x[2] = {1.0, 1.0};

//     qp_real val = obj_value(&P, c, x, temp);
//     x[0] = -1.0;
//     x[1] = -1.0;

//     qp_real val1 = obj_value(&P, c, x, temp);

//     ASSERT_EQ(3.0, val);
//     ASSERT_EQ(-1.0, val1);
// }

TEST_F(standAloneTests, func_SparseMatrixTranspose)
{
    /* Square Matrix Test */
    // qp_int *Ptir, *Ptjc;
    // qp_real *Ptpr;
    // Ptir = (qp_int *)MALLOC(P.nnz * sizeof(qp_int));
    // Ptjc = (qp_int *)MALLOC((P.m + 1) * sizeof(qp_int));
    // Ptpr = (qp_real *)MALLOC(P.nnz * sizeof(qp_real));
    // SparseMatrixSetup(P.n, P.m, P.nnz, Ptjc, Ptir, Ptpr, &Pt);

    // SparseMatrixTranspose(&P, &Pt);
    // ASSERT_TRUE(validate_smat(&Pt, diag_msg)) << diag_msg;
    // ASSERT_TRUE(validate_smat_smatT(&P, &Pt, diag_msg)) << diag_msg;
    // FREE(Ptir);
    // FREE(Ptjc);
    // FREE(Ptpr);
    /* Square Matrix Test */

    /* Rectangular Matrix Test*/
    // qp_int *Atir, *Atjc;
    // qp_real *Atpr;
    // Atir = (qp_int *)MALLOC(A.nnz * sizeof(qp_int));
    // Atjc = (qp_int *)MALLOC(A.m + 1 * sizeof(qp_int));
    // Atpr = (qp_real *)MALLOC(A.nnz * sizeof(qp_real));
    // SparseMatrixSetup(A.n, A.m, A.nnz, Atjc, Atir, Atpr, &At);

    // SparseMatrixTranspose(&A, &At);
    // ASSERT_TRUE(validate_smat(&At, diag_msg)) << diag_msg;
    // //ASSERT_TRUE(validate_smat_smatT(&A, &At, diag_msg)) << diag_msg;
    // FREE(Atjc);
    // FREE(Atpr);
    // FREE(Atir);

    /* Rectangular Matrix Test*/
}

TEST_F(standAloneTests, func_innerproduct)
{
    testing::internal::CaptureStdout();
    std::cout << "Data read\n";
    std::string output = testing::internal::GetCapturedStdout();
    qp_real x[3] = {1.0, 2.0, 3.0};
    qp_real y[3] = {-1.0, -2.0, -3.0};
    ASSERT_EQ(14, innerproduct(x, x, 3));
    ASSERT_EQ(-14, innerproduct(x, y, 3));
}

TEST_F(standAloneTests, func_checksign)
{
    double x[2] = {1, 2};

    ASSERT_EQ(1, checksign(x, x, -2, 2));
    ASSERT_EQ(0, checksign(x, x, 1, 2));
}

// TEST_F(standAloneTests, func_SparseMatrixMultiply)
// {

//     SparseMatrixSetup(2, 2, 4, A_jc, A_ir, A_pr, &A);

//     double x[2], y[2];

//     x[0] = 1;
//     x[1] = 1;

//     SparseMatrixMultiply(&A, x, y, 1);

//     ASSERT_EQ(-3, y[0]);
//     ASSERT_EQ(-7, y[1]);

//     y[0] = -y[0];
//     y[1] = -y[1];

//     SparseMatrixMultiply(&A, x, y, 0);

//     ASSERT_EQ(0, y[0]);
//     ASSERT_EQ(0, y[1]);
// }

// TEST_F(standAloneTests, func_SparseMatrixTransMultiply)
// {
//     SparseMatrixSetup(2, 2, 4, A_jc, A_ir, A_pr, &A);

//     double x[2], y[2];

//     x[0] = 1;
//     x[1] = 1;

//     SparseMatrixTransMultiply(&A, x, y, 1);

//     ASSERT_EQ(-4, y[0]);
//     ASSERT_EQ(-6, y[1]);

//     y[0] = -y[0];
//     y[1] = -y[1];

//     SparseMatrixTransMultiply(&A, x, y, 0);

//     ASSERT_EQ(0, y[0]);
//     ASSERT_EQ(0, y[1]);
// }

// TEST_F(standAloneTests, func_densetosparse)
// {

//     /* Sparse Matrices; square and rectangular */
//     smat sq, rec;

//     /* Square Matrix Pointer */
//     qp_int sq_m = 3, sq_n = 3;
//     qp_real sq_pr[9] = {1.0, 0.0, 2.0, 3.0, 0.0, 1.1, 0.1, 0.0, 3.0};

//     /* Rectangular Matrix Pointer */
//     qp_int rec_m = 3, rec_n = 1;
//     qp_real rec_pr[3] = {0.0, 1.0, 0.0};

//     /*  Perform computations */
//     densetosparse(sq_m, sq_n, sq_pr, &sq);
//     densetosparse(rec_m, rec_n, rec_pr, &rec);

//     /* Validate the Matrix */
//     ASSERT_TRUE(validate_smat(&sq, diag_msg)) << diag_msg;
//     ASSERT_TRUE(validate_smat(&rec, diag_msg)) << diag_msg;
//     //std::cout << msg << std::endl;

//     /* Data Integrity Checks */
//     ASSERT_TRUE(compare_dmat_C_smat(&sq, sq_pr, sq_m, sq_n));
//     ASSERT_TRUE(compare_dmat_C_smat(&rec, rec_pr, rec_m, rec_n));
// }

int main(int argc, char **argv)
{

    // data = psr.getdateptr();

    testing::InitGoogleTest(&argc, argv);
    // testing::InitGoogleTest(&argc, argv);
    // obj.SetUp(psr.getdateptr());

    //::testing::UnitTest::GetInstance()->Setup(psr.getdateptr());
    // testing::UnitTest::Setup(psr.getdateptr());
    // standAloneTests::init(psr.getdateptr());

    //::testing::UnitTest::GetInstance()->Setup(psr.getdateptr());

    return ::testing::UnitTest::GetInstance()->Run();
    // return RUN_ALL_TESTS();
}