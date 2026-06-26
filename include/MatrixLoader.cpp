#include "MatrixLoader.h"
#include <fstream>
#include <iostream>
#include <cstdlib>

static void checkFile(std::ifstream& file, const std::string& name)
{
    if (!file.is_open())
    {
        std::cerr << "Cannot open file: " << name << std::endl;
        std::exit(EXIT_FAILURE);
    }
}

template<typename T>
static void readChecked(std::ifstream& in, T& value, const std::string& msg)
{
    if (!(in >> value))
    {
        std::cerr << "ERROR reading: " << msg << std::endl;
        std::exit(EXIT_FAILURE);
    }
}

CSRMatrix MatrixLoader::load(const std::string& prefix)
{
    CSRMatrix A;

    // ----------------------------
    // INFO
    // ----------------------------
    std::ifstream info(prefix + "_info.txt");
    checkFile(info, prefix + "_info.txt");

    readChecked(info, A.n, "n");
    readChecked(info, A.nnz, "nnz");

    // ----------------------------
    // VALUES
    // ----------------------------
    A.values.resize(A.nnz);

    std::ifstream values(prefix + "_values.txt");
    checkFile(values, prefix + "_values.txt");

    for (int i = 0; i < A.nnz; i++)
        readChecked(values, A.values[i], "values[" + std::to_string(i) + "]");

    // ----------------------------
    // COLUMNS
    // ----------------------------
    A.columns.resize(A.nnz);

    std::ifstream columns(prefix + "_columns.txt");
    checkFile(columns, prefix + "_columns.txt");

    for (int i = 0; i < A.nnz; i++)
    {
        readChecked(columns, A.columns[i], "columns[" + std::to_string(i) + "]");

        if (A.columns[i] < 1 || A.columns[i] > A.n)
        {
            std::cerr << "ERROR: column index out of range at " << i
                      << " value=" << A.columns[i] << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }

    // ----------------------------
    // ROW POINTERS
    // ----------------------------
    A.rowPointers.resize(A.n + 1);

    std::ifstream rowptr(prefix + "_rowptr.txt");
    checkFile(rowptr, prefix + "_rowptr.txt");

    for (int i = 0; i <= A.n; i++)
        readChecked(rowptr, A.rowPointers[i], "rowPtr[" + std::to_string(i) + "]");

    // CSR sanity checks
    for (int i = 0; i < A.n; i++)
    {
        if (A.rowPointers[i] < 1 || A.rowPointers[i] > A.nnz + 1)
        {
            std::cerr << "ERROR: rowPtr out of bounds at " << i << std::endl;
            std::exit(EXIT_FAILURE);
        }

        if (A.rowPointers[i] > A.rowPointers[i + 1])
        {
            std::cerr << "ERROR: rowPtr not monotonic at " << i << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }

    if (A.rowPointers[A.n] != A.nnz + 1)
    {
        std::cerr << "ERROR: rowPtr[n] != nnz+1" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    // ----------------------------
    // RHS
    // ----------------------------
    A.rhs.resize(A.n);

    std::ifstream rhs(prefix + "_rhs.txt");
    checkFile(rhs, prefix + "_rhs.txt");

    for (int i = 0; i < A.n; i++)
        readChecked(rhs, A.rhs[i], "rhs[" + std::to_string(i) + "]");

    return A;
}