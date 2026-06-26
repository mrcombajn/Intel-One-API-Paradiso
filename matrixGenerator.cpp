#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <filesystem>
#include <cmath>
#include <algorithm>

namespace fs = std::filesystem;

int main(int argc, char* argv[])
{
    if (argc < 3)
    {
        std::cout << "Usage:\n";
        std::cout << argv[0] << " matrix_size bandwidth\n";
        return 0;
    }

    const int n = std::stoi(argv[1]);
    const int bandwidth = std::stoi(argv[2]);

    std::mt19937 rng(12345);
    std::uniform_real_distribution<double> valueDist(-1.0, 1.0);
    std::uniform_real_distribution<double> prob(0.0, 1.0);

    fs::create_directories("generated");

    std::vector<std::vector<std::pair<int,double>>> L(n);

    for (int i = 0; i < n; i++)
    {
        int begin = std::max(0, i - bandwidth / 2);

        for (int j = begin; j <= i; j++)
        {
            if (i == j)
            {
                // diagonal strictly positive -> SPD guarantee
                L[i].push_back({j, 1.0 + std::abs(valueDist(rng))});
            }
            else if (prob(rng) < 0.3)
            {
                L[i].push_back({j, valueDist(rng)});
            }
        }
    }

    std::vector<std::vector<double>> A(n, std::vector<double>(n, 0.0));

    for (int i = 0; i < n; i++)
    {
        for (auto [k1, val1] : L[i])
        {
            for (int j = 0; j <= i; j++)
            {
                for (auto [k2, val2] : L[j])
                {
                    if (k1 == k2)
                    {
                        A[i][j] += val1 * val2;
                    }
                }
            }
        }
    }

    std::vector<double> values;
    std::vector<int> columns;
    std::vector<int> rowPtr(n + 1);

    int nnz = 0;
    rowPtr[0] = 1;

    const double EPS = 1e-12;

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j <= i; j++)
        {
            if (std::abs(A[i][j]) > EPS)
            {
                values.push_back(A[i][j]);
                columns.push_back(j + 1); // 1-based for PARDISO
                nnz++;
            }
        }

        rowPtr[i + 1] = nnz + 1;
    }

    std::vector<double> rhs(n, 1.0);

    std::string prefix = "generated/matrix_" + std::to_string(n);

    std::ofstream info(prefix + "_info.txt");
    info << n << "\n" << nnz << "\n";

    std::ofstream val(prefix + "_values.txt");
    for (auto v : values) val << v << "\n";

    std::ofstream col(prefix + "_columns.txt");
    for (auto c : columns) col << c << "\n";

    std::ofstream row(prefix + "_rowptr.txt");
    for (auto r : rowPtr) row << r << "\n";

    std::ofstream rhsFile(prefix + "_rhs.txt");
    for (auto b : rhs) rhsFile << b << "\n";

    std::cout << "SPD matrix generated WITHOUT dense storage\n";
    std::cout << "Size : " << n << "\n";
    std::cout << "NNZ  : " << nnz << "\n";

    return 0;
}