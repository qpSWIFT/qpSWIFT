#include <iostream>
#include "test_data_structures.hpp"
#include <fstream>
#include <sstream>
#include <cstring>

class parser
{
private:
    standalone_data data;
    std::ifstream fin;
    char msg[50];

public:
    bool File_read_status;
    parser();
    ~parser();

    bool readfile(const char *fileName);

    bool readvectorInt(qp_int *pr, qp_int len, char *delimiter);
    bool readvectorFloat(qp_real *pr, qp_int len, char *delimiter);

    standalone_data *getdateptr() { return &data; };

    char *getdiag_msg() { return msg; };

    bool readinitData();
};