
#include "parser.hpp"

parser::parser()
{
}

bool parser::readinitData()
{
    bool status = false;
    std::string str;
    std::getline(fin, str);
    std::stringstream buff(str);

    //buff >> data.n >> data.m >> data.p >> data.P_nnz >> data.A_nnz >> data.G_nnz;

    buff >> data.n;
    buff.ignore();
    buff >> data.m;
    buff.ignore();
    buff >> data.p;
    buff.ignore();
    buff >> data.P_nnz;
    buff.ignore();
    buff >> data.A_nnz;
    buff.ignore();
    buff >> data.G_nnz;

    status = true;
    return status;
}

bool parser::readvectorInt(qp_int *pr, qp_int len, char *delimiter)
{
    bool status = false;
    std::string str;

    qp_int idx = 0;

    std::getline(fin, str);
    std::stringstream buff(str);

    qp_int temp;
    while (buff >> temp && idx < len)

    {
        pr[idx++] = temp;
        if (buff.peek() == delimiter[0])
        {
            buff.ignore();
        }
    }

    status = true;
    return status;
}

bool parser::readvectorFloat(qp_real *pr, qp_int len, char *delimiter)
{
    bool status = false;
    std::string str;

    qp_int idx = 0;

    std::getline(fin, str);
    std::stringstream buff(str);

    qp_real temp;
    while (buff >> temp && idx < len)
    {
        pr[idx++] = temp;
        if (buff.peek() == delimiter[0])
        {
            buff.ignore();
        }
    }

    status = true;
    return status;
}

bool parser::readfile(const char *fileName)
{

    fin.open(fileName);
    if (fin.fail())
    {
        printf("Cannot open the specified file\n");
        File_read_status = false;
        return File_read_status;
    }
    else
    {
        File_read_status = true;
        char delim[] = {','};
        if (!readinitData())
        {

            File_read_status = false;
            return File_read_status;
        }
        else
        {

            data.Ppr = (qp_real *)malloc(data.P_nnz * sizeof(qp_real));
            data.Pir = (qp_int *)malloc(data.P_nnz * sizeof(qp_int));
            data.Pjc = (qp_int *)malloc((data.n + 1) * sizeof(qp_int));

            data.Apr = (qp_real *)malloc(data.A_nnz * sizeof(qp_real));
            data.Air = (qp_int *)malloc(data.A_nnz * sizeof(qp_int));
            data.Ajc = (qp_int *)malloc((data.n + 1) * sizeof(qp_int));

            data.Gpr = (qp_real *)malloc(data.G_nnz * sizeof(qp_real));
            data.Gir = (qp_int *)malloc(data.G_nnz * sizeof(qp_int));
            data.Gjc = (qp_int *)malloc((data.n + 1) * sizeof(qp_int));

            data.s = (qp_real *)malloc(data.m * sizeof(qp_real));
            data.z = (qp_real *)malloc(data.m * sizeof(qp_real));
            data.delta_s = (qp_real *)malloc(data.m * sizeof(qp_real));
            data.delta_z = (qp_real *)malloc(data.m * sizeof(qp_real));

            if (!readvectorFloat(data.Ppr, data.P_nnz, delim))
            {
                File_read_status = false;
                strcpy(msg, "Ppr data read failed\n");
                return File_read_status;
            }

            if (!readvectorInt(data.Pir, data.P_nnz, delim))
            {
                File_read_status = false;
                strcpy(msg, "Pir data read failed\n");
                return File_read_status;
            }

            if (!readvectorInt(data.Pjc, data.n + 1, delim))
            {
                File_read_status = false;
                strcpy(msg, "Pjc data read failed\n");
                return File_read_status;
            }

            if (!readvectorFloat(data.Apr, data.A_nnz, delim))
            {
                File_read_status = false;
                strcpy(msg, "Apr data read failed\n");
                return File_read_status;
            }

            if (!readvectorInt(data.Air, data.A_nnz, delim))
            {
                File_read_status = false;
                strcpy(msg, "Air data read failed\n");
                return File_read_status;
            }

            if (!readvectorInt(data.Ajc, data.n + 1, delim))
            {
                File_read_status = false;
                strcpy(msg, "Ajc data read failed\n");
                return File_read_status;
            }

            if (!readvectorFloat(data.Gpr, data.G_nnz, delim))
            {
                File_read_status = false;
                strcpy(msg, "Gpr data read failed\n");
                return File_read_status;
            }

            if (!readvectorInt(data.Gir, data.G_nnz, delim))
            {
                File_read_status = false;
                strcpy(msg, "Gir data read failed\n");
                return File_read_status;
            }

            if (!readvectorInt(data.Gjc, data.n + 1, delim))
            {
                File_read_status = false;
                strcpy(msg, "Gjc data read failed\n");
                return File_read_status;
            }

            if (!readvectorFloat(data.s, data.m, delim))
            {
                File_read_status = false;
                strcpy(msg, "s data read failed\n");
                return File_read_status;
            }

            if (!readvectorFloat(data.delta_s, data.m, delim))
            {
                File_read_status = false;
                strcpy(msg, "delta_s data read failed\n");
                return File_read_status;
            }

            if (!readvectorFloat(data.z, data.m, delim))
            {
                File_read_status = false;
                strcpy(msg, "z data read failed\n");
                return File_read_status;
            }
            if (!readvectorFloat(data.delta_z, data.m, delim))
            {
                File_read_status = false;
                strcpy(msg, "delta_z data read failed\n");
                return File_read_status;
            }
        }
    }

    return File_read_status;
}

parser::~parser()
{
    if (File_read_status)
    {
        if (data.Pir)
            free(data.Pir);

        if (data.Pjc)
            free(data.Pjc);

        if (data.Ppr)
            free(data.Ppr);

        if (data.Air)
            free(data.Air);

        if (data.Ajc)
            free(data.Ajc);

        if (data.Apr)
            free(data.Apr);

        if (data.Gir)
            free(data.Gir);

        if (data.Gjc)
            free(data.Gjc);

        if (data.Gpr)
            free(data.Gpr);
    }
}