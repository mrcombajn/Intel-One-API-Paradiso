#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <unordered_map>
#include <set>

int main()
{
    int n = 1000;
    double density = 0.005; // 0.5%

    std::mt19937 gen(123);
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    // --- temporary dense matrix for SPD construction ---
    std::vector<std::vector<double>> A(n, std::vector<double>(n, 0.0));

    // -------------------------
    // STEP 1: build random B
    // -------------------------
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (dist(gen) < density)
            {
                double v = dist(gen);
                A[i][j] += v;
            }
        }
    }

    // -------------------------
    // STEP 2: A = B^T * B
    // -------------------------
    std::vector<std::vector<double>> SPD(n, std::vector<double>(n, 0.0));

    for (int i = 0; i < n; i++)
        for (int k = 0; k < n; k++)
            if (A[i][k] != 0)
                for (int j = 0; j < n; j++)
                    SPD[i][j] += A[k][i] * A[k][j];

    // -------------------------
    // STEP 3: add alpha*I (stabilization)
    // -------------------------
    double alpha = 1e-3;
    for (int i = 0; i < n; i++)
        SPD[i][i] += alpha;

    // -------------------------
    // STEP 4: convert to CSR (upper triangle only)
    // -------------------------
    std::vector<double> values;
    std::vector<int> columns;
    std::vector<int> rowPtr(n + 1, 1);

    int nnz = 0;

    for (int i = 0; i < n; i++)
    {
        rowPtr[i] = nnz + 1;

        for (int j = 0; j < n; j++)
        {
            if (SPD[i][j] != 0.0)
            {
                values.push_back(SPD[i][j]);
                columns.push_back(j + 1); // PARDISO = 1-based
                nnz++;
            }
        }
    }

    rowPtr[n] = nnz + 1;

    // -------------------------
    // RHS
    // -------------------------
    std::vector<double> rhs(n, 1.0);

    // -------------------------
    // SAVE
    // -------------------------
    std::ofstream out("matrix_small.txt");

    out << n << " " << nnz << "\n\n";

    out << "values\n";
    for (auto v : values) out << v << "\n";

    out << "columns\n";
    for (auto c : columns) out << c << "\n";

    out << "rowPointers\n";
    for (auto r : rowPtr) out << r << "\n";

    out << "rhs\n";
    for (auto r : rhs) out << r << "\n";

    std::cout << "SPD matrix generated. nnz = " << nnz << "\n";

    return 0;
}