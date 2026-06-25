#include "MatrixLoader.h"
#include <fstream>
#include <iostream>

CSRMatrix MatrixLoader::load(const std::string &filename)
{
    CSRMatrix A;

    std::ifstream in(filename);
    if (!in)
    {
        std::cerr << "ERROR: cannot open file " << filename << "\n";
        exit(1);
    }

    in >> A.n >> A.nnz;

    std::string tag;

    // ---------------- values ----------------
    in >> tag; // values
    A.values.resize(A.nnz);

    for (int i = 0; i < A.nnz; i++)
        in >> A.values[i];

    // ---------------- columns ----------------
    in >> tag; // columns
    A.columns.resize(A.nnz);

    for (int i = 0; i < A.nnz; i++)
        in >> A.columns[i];

    // ---------------- rowPointers ----------------
    in >> tag; // rowPointers
    A.rowPointers.resize(A.n + 1);

    for (int i = 0; i < A.n + 1; i++)
        in >> A.rowPointers[i];

    // ---------------- rhs ----------------
    in >> tag; // rhs
    A.rhs.resize(A.n);

    for (int i = 0; i < A.n; i++)
        in >> A.rhs[i];

    return A;
}