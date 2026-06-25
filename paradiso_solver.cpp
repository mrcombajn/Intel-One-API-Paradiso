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
    if (argc < 2)
    {
        std::cout << "Usage: ./solver matrix_file.txt\n";
        return 0;
    }

    std::string filename = argv[1];
    CSRMatrix A =  MatrixLoader::load(filename);

    std::cout << "min col = "
          << *std::min_element(A.columns.begin(), A.columns.end())
          << "\n";

    std::cout << "max col = "
              << *std::max_element(A.columns.begin(), A.columns.end())
              << "\n";

    int n = A.n;
    int nrhs = 1;

    std::cout << A.rowPointers[0] << " ... " << A.rowPointers[A.n] << "\n";
    std::vector<double> x(n, 0.0);

    void *pt[64] = {0};
    MKL_INT iparm[64];

    for (int i = 0; i < 64; i++)
        iparm[i] = 0;

    // --------------------------
    // KONFIGURACJA PARDISO
    // --------------------------
    iparm[0] = 0;   // no solver default
    iparm[1] = 2;   // METIS ordering
    iparm[7] = 2;   // iterative refinement steps
    iparm[9] = 13;  // pivot perturbation
    iparm[10] = 1;  // scaling
    iparm[17] = -1; // output nnz in LU
    iparm[18] = -1; // memory estimate
    iparm[26] = 1;
    iparm[27] = 1;

    MKL_INT maxfct = 1;
    MKL_INT mnum = 1;
    MKL_INT mtype = 2; // SPD real matrix
    MKL_INT msglvl = 1;
    MKL_INT error = 0;

    MKL_INT phase;

    // ==========================
    // 1. ANALIZA SYMBOLICZNA
    // ==========================
    phase = 11;

    auto t1 = std::chrono::high_resolution_clock::now();

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

    auto t2 = std::chrono::high_resolution_clock::now();

    if (error)
    {
        std::cerr << "Error in analysis: " << error << "\n";
        return -1;
    }

    // ==========================
    // 2. FAKTORYZACJA
    // ==========================
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

    // ==========================
    // 3. ROZWIĄZANIE UKŁADU
    // ==========================
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

    // ==========================
    // 4. ZWOLNIENIE PAMIĘCI
    // ==========================
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

    // ==========================
    // WYNIKI
    // ==========================
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

    // opcjonalna kontrola wyniku
    std::cout << "\nSample x[0] = " << x[0] << "\n";

    return 0;
}