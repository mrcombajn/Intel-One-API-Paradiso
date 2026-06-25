#pragma once

#include <vector>
#include <string>
#include "mkl.h"

struct CSRMatrix
{
    int n;
    int nnz;

    std::vector<double> values;
    std::vector<MKL_INT> columns;
    std::vector<MKL_INT> rowPointers;

    std::vector<double> rhs;
};

class MatrixLoader
{
public:
    static CSRMatrix load(const std::string &filename);
};