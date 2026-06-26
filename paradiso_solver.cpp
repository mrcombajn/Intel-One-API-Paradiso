#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <chrono>
#include "mkl.h"
#include "include/MatrixLoader.h"

void printStats(MKL_INT *iparm)
{
    std::cout << "\n--- PARDISO STATS ---\n";
    std::cout << "Peak memory (internal): " << iparm[14] << "\n";
    std::cout << "Peak memory (external): " << iparm[15] << "\n";
    std::cout << "Number of nonzeros in factors: " << iparm[17] << "\n";
    std::cout << "Peak memory symbolic: " << iparm[18] << "\n";
}

int main(int argc, char **argv)
{

    std::string filename = argv[1];
    CSRMatrix A =  MatrixLoader::load(filename);

    std::cout << "n = " << A.n << std::endl;
    std::cout << "nnz = " << A.nnz << std::endl;
    std::cout << "rowPointers[n] = " << A.rowPointers[A.n] << std::endl;

    int emptyRows = 0;
    int missingDiagonal = 0;

    for (int i = 0; i < A.n; i++)
    {
        if (A.rowPointers[i] == A.rowPointers[i + 1])
            emptyRows++;

        bool diag = false;
        for (int k = A.rowPointers[i] - 1; k < A.rowPointers[i + 1] - 1; k++)
        {
            if (A.columns[k] == i + 1)
            {
                diag = true;
                break;
            }
        }

        if (!diag)
            missingDiagonal++;
    }

    std::cout << "Empty rows      = " << emptyRows << '\n';
    std::cout << "Missing diagonal= " << missingDiagonal << '\n';


    int n = A.n;
    int nrhs = 1;

    std::cout << A.rowPointers[0] << " ... " << A.rowPointers[A.n] << "\n";
    std::vector<double> x(n, 0.0);

    void *pt[64] = { nullptr };
    MKL_INT iparm[64];

    for (int & i : iparm)
        i = 0;

    iparm[0] = 1;
    iparm[1] = 2;

    MKL_INT maxfct = 1;
    MKL_INT mnum = 1;
    MKL_INT mtype = 2; // SPD real matrix
    MKL_INT msglvl = 2;
    MKL_INT error = 0;

    MKL_INT phase;

    phase = 11;

    auto t1 = std::chrono::high_resolution_clock::now();

    pardiso(pt, &maxfct, &mnum, &mtype, &phase,
            &n,
            A.values.data(),
            A.rowPointers.data(),
            A.columns.data(),
            nullptr,
            &nrhs,
            iparm,
            &msglvl,
            nullptr,
            nullptr,
            &error);

    auto t2 = std::chrono::high_resolution_clock::now();

    if (error)
    {
        std::cerr << "Error in analysis: " << error << "\n";
        return -1;
    }

    phase = 22;

    auto t3 = std::chrono::high_resolution_clock::now();

    pardiso(pt, &maxfct, &mnum, &mtype, &phase,
            &n,
            A.values.data(),
            A.rowPointers.data(),
            A.columns.data(),
            NULL,
            &nrhs,
            iparm,
            &msglvl,
            A.rhs.data(),
            x.data(),
            &error);

    auto t4 = std::chrono::high_resolution_clock::now();

    if (error)
    {
        std::cerr << "Error in factorization: " << error << "\n";
        return -1;
    }

    phase = 33;

    auto t5 = std::chrono::high_resolution_clock::now();

    pardiso(pt, &maxfct, &mnum, &mtype, &phase,
            &n,
            A.values.data(),
            A.rowPointers.data(),
            A.columns.data(),
            NULL,
            &nrhs,
            iparm,
            &msglvl,
            A.rhs.data(),
            x.data(),
            &error);

    auto t6 = std::chrono::high_resolution_clock::now();

    if (error)
    {
        std::cerr << "Error in solve: " << error << "\n";
        return -1;
    }

    phase = -1;

    pardiso(pt, &maxfct, &mnum, &mtype, &phase,
            &n,
            NULL,
            A.rowPointers.data(),
            A.columns.data(),
            NULL,
            &nrhs,
            iparm,
            &msglvl,
            NULL,
            NULL,
            &error);

    auto t7 = std::chrono::high_resolution_clock::now();

    std::cout << "\nFILE: " << filename << "\n";
    std::cout << "N: " << n << "  NNZ: " << A.nnz << "\n";

    std::cout << "Analysis time   : "
              << std::chrono::duration<double>(t2 - t1).count()
              << " s\n";

    std::cout << "Factorization   : "
              << std::chrono::duration<double>(t4 - t3).count()
              << " s\n";

    std::cout << "Solve           : "
              << std::chrono::duration<double>(t6 - t5).count()
              << " s\n";

    std::cout << "Total PARDISO   : "
              << std::chrono::duration<double>(t6 - t1).count()
              << " s\n";

    printStats(iparm);

    std::cout << "\nSample x[0] = " << x[0] << "\n";

    return 0;
}