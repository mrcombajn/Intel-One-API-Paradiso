#include <iostream>
#include "mkl.h"

int main() {

    MKL_INT n = 3;
    MKL_INT nrhs = 1;

    double a[] = {4,-1,-1,4,-1,-1,4};

    MKL_INT ia[] = {1,3,6,8};
    MKL_INT ja[] = {1,2,1,2,3,2,3};

    double b[] = {15,10,10};
    double x[] = {0,0,0};

    void *pt[64] = {0};
    MKL_INT iparm[64] = {0};

    iparm[0] = 1;
    iparm[1] = 2;

    MKL_INT maxfct = 1;
    MKL_INT mnum = 1;
    MKL_INT mtype = 2;
    MKL_INT phase = 13;
    MKL_INT error = 0;
    MKL_INT msglvl = 1;

    pardiso(
        pt,
        &maxfct,
        &mnum,
        &mtype,
        &phase,
        &n,
        a,
        ia,
        ja,
        NULL,
        &nrhs,
        iparm,
        &msglvl,
        b,
        x,
        &error
    );

    std::cout << "error = " << error << std::endl;

    if(error == 0) {
        for(int i=0;i<n;i++)
            std::cout << x[i] << std::endl;
    }

    phase = -1;

    pardiso(
        pt,
        &maxfct,
        &mnum,
        &mtype,
        &phase,
        &n,
        NULL,
        ia,
        ja,
        NULL,
        &nrhs,
        iparm,
        &msglvl,
        NULL,
        NULL,
        &error
    );

    return 0;
}