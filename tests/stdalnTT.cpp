#include "test.hpp"
#define GTEST_COUT std::cerr << "[          ] [ INFO ]"

TEST(standAloneTests, func_updatekktmatrix)
{
    smat P;
    SparseMatrixSetup(P_m, P_n, P_nnz, P_jc, P_ir, P_pr, &P);
    char diag_msg[50];
    ASSERT_TRUE(check_matrix_symmetric(&P, diag_msg)) << diag_msg;
}
TEST(standAloneTests, func_updatekktmatrix_b)
{
}

TEST(standAloneTests, func_computeresiduals)
{

    QP *myQP;

    myQP = QP_SETUP(P_n, G_m, A_m, Pjc, Pir, Ppr, Ajc, Air, Apr, Gjc, Gir, Gpr, c, h, b, 0.0, NULL);

    computeresiduals(myQP);

    ASSERT_EQ(myQP->stats->n_rx, );
    ASSERT_EQ(myQP->stats->n_rx, );
    ASSERT_EQ(myQP->stats->n_rx, );
}

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return ::testing::UnitTest::GetInstance()->Run();
    return RUN_ALL_TESTS();
}