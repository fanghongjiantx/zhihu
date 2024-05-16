/*
* @Description: 计算力学(七)单元刚度矩阵中代码
* @Author: FangHongjianTX
* @Date: 2024-04-30
* @Email: shijunc@foxmail.com
*
* Copyright (c) 2024 fanghongjiantx
* SPDX-License-Identifier: BSD-3-Clause
*/
#include <iostream>
using namespace std;

/// 二维应力应变矩阵D
/// \param Iopt =0为平面应力情况，=1为平面应变情况，=2为轴对称情况
void eld(double E, double u, double(*D)[4], int Iopt) {
    switch (Iopt) {
        case 0:
            D[0][0] = D[1][1] = E / (1 - u * u);
            D[2][2] = E / (2 + 2 * u);
            D[0][1] = D[1][0] = u * E / (1 - u * u);
            break;
        case 1:
            u /= (1 - u);
            E /= (1 - u * u);
            D[0][0] = D[1][1] = E / (1 - u * u);
            D[2][2] = E / (2 + 2 * u);
            D[0][1] = D[1][0] = u * E / (1 - u * u);
        case 2:
            D[0][1] = D[1][0] = D[0][3] = D[3][0] = D[1][3] = D[3][1] = u * E / (1 - u * u);
            break;
        default:
            cout << "Please enter the correct Iopt" << endl;
            exit(0);
    }
}

/// 三维应力应变关系矩阵D
void eld_cube(double E, double u, double(*D)[6]) {
    for (int i = 0; i < 6; i++) { for (int j = 0; j < 6; j++) { D[i][j] = 0; }}
    double coff = E / ((1 + u) * (1 - 2 * u));
    D[0][0] = D[1][1] = D[2][2] = (1 - u) * coff;
    D[0][1] = D[0][2] = D[1][0] = D[1][2] = D[2][0] = D[2][1] = u * coff;
    D[3][3] = D[4][4] = D[5][5] = (1 - 2 * u) * coff / 2;
}

/// 二维刚度矩阵
/// \param Thick 单元厚度，仅用于平面应力情况
/// \param xo 插值节点总体x坐标，长度为8的数组指针
/// \param yo 插值节点总体y坐标，长度为8的数组指针
/// \param ek 按行排列的单元刚度矩阵，长度为256的数组指针
/// \param Ele 单元节点编号，长度为8的数组的指针
/// \param Iopt =0为平面应力情况，=1为平面应变情况，=2为轴对称情况
/// \param nGauss =3时采用3*3的高斯积分，否测采用2*2的高斯积分
void CacuEk(double E, double u, double Thick, double *xo, double *yo, double *ek, int *Ele, int Iopt, int nGauss) {
    double Nx[8] = {0}, Ny[8] = {0}, N[8] = {0}, B[4][16] = {0}, D[4][4] = {0}, AB[3] = {0}, H[3] = {
            0}, det, coeff, x;
    int ij, i, j, k, l, m, n;
    for (i = 0; i < 256; i++) { ek[i] = 0; }
    for (i = 0; i < 4; i++) { for (j = 0; j < 16; j++) { B[i][j] = 0; }}
    for (i = 0; i < 4; i++) { for (j = 0; j < 4; j++) { D[i][j] = 0; }}
    eld(E, u, D, Iopt);
    if (nGauss == 3) {
        AB[0] = -sqrt(0.6);
        AB[1] = 0;
        AB[2] = -AB[0];
        H[0] = 5.0 / 9.0;
        H[1] = 8.0 / 9.0;
        H[2] = H[0];
        ij = 1;
    } else {
        AB[0] = -1 / sqrt(3);
        AB[1] = 0;
        AB[2] = -AB[0];
        H[0] = 1;
        H[1] = 0;
        H[2] = H[0];
        ij = 2;
    }

    for (i = 0; i < 3; i += ij) {
        for (j = 0; j < 3; j += ij) {
            det = Inv_jaco(AB[i], AB[j], xo, yo, Nx, Ny, Ele);
            x = 0;
            coeff = det * H[j] * H[i];
            if (Iopt == 2) {
                FN(AB[i], AB[j], N, Ele);
                for (k = 0; k < 8; k++) { x += xo[k] * N[k]; }
            }
            if (Iopt == 0) { coeff *= Thick; }
            else if (Iopt == 1) {}
            else if (Iopt == 2) { coeff *= 2 * PI * x; }
            for (k = 0; k < 8; k++) {
                B[0][2 * k] = Nx[k];
                B[2][2 * k] = Ny[k];
                B[1][2 * k + 1] = Ny[k];
                B[2][2 * k + 1] = Nx[k];
                if (Iopt == 2) { B[3][2 * k] = N[k] / x; }
            }
            for (k = 0; k < 16; k++) {
                for (l = 0; l <= k; l++) {
                    for (m = 0; m < 4; m++) {
                        for (n = 0; n < 4; n++) {
                            ek[16 * k + l] += B[m][k] * D[m][n] * B[n][l] * coeff;
                        }
                    }
                }
            }
        }
    }
    for (i = 0; i < 16; i++) {
        for (j = i + 1; j < 16; j++) {
            ek[16 * i + j] = ek[16 * j + i];
        }
    }
}

/// 三角单元刚度矩阵
void CacuEk(double E, double u, double Thick, double *xo, double *yo, double *ek, int *Ele, int Iopt) {
    double(*D)[4] = new double[4][4];
    int i, j, k, l, m, n;
    for (i = 0; i < 4; i++) { for (j = 0; j < 4; j++) { D[i][j] = 0; }}
    double Nx[6] = {0}, Ny[6] = {0}, N[6] = {0}, B[4][12] = {0}, det, coeff, x;
    eld(E, u, D, Iopt);
    double integral[4][4] = {{1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, -27.0 / 48.0},
                             {0.6,       0.2,       0.2,       25.0 / 48.0},
                             {0.2,       0.6,       0.2,       25.0 / 48.0},
                             {0.2,       0.2,       0.6,       25.0 / 48.0}};

    for (i = 0; i < 4; i++) {
        x = 0;
        if (Iopt == 2) {
            FN(integral[i][0], integral[i][1], integral[i][2], N, Ele);
            for (k = 0; k < 6; k++) { x += xo[k] * N[k]; }
        }
        det = Inv_jaco(integral[i][0], integral[i][1], integral[i][2], xo, yo, Nx, Ny, Ele);
        coeff = det * integral[i][3] / 2.0;
        if (Iopt == 0) { coeff *= Thick; }
        else if (Iopt == 1) {}
        else if (Iopt == 2) { coeff *= 2 * PI * x; }
        for (k = 0; k < 6; k++) {
            B[0][2 * k] = Nx[k];
            B[2][2 * k] = Ny[k];
            B[1][2 * k + 1] = Ny[k];
            B[2][2 * k + 1] = Nx[k];
            if (Iopt == 2) { B[3][2 * k] = N[k] / x; }
        }
        for (k = 0; k < 12; k++) {
            for (l = 0; l <= k; l++) {
                for (m = 0; m < 4; m++) {
                    for (n = 0; n < 4; n++) {
                        ek[12 * k + l] += B[m][k] * D[m][n] * B[n][l] * coeff;
                    }
                }
            }
        }
    }
    for (i = 0; i < 12; i++) {
        for (j = i + 1; j < 12; j++) {
            ek[12 * i + j] = ek[12 * j + i];
        }
    }
    delete[]D;
}

/// 三维单元刚度矩阵
void CacuEk_cube(double E, double u, double *xo, double *yo, double *zo, double *ek, int *Ele) {
    double Nx[20] = {0}, Ny[20] = {0}, Nz[20] = {0}, B[6][60] = {0}, D[6][6] = {0}, det, a1 = 320.0 / 361.0, a2 = sqrt(19.0 / 30.0), b1 = 121.0 / 361.0, b2 = sqrt(19.0 / 33.0);
    int i, j, k, l, m, n;
    double AB[14][4] = {{a1, -a2, 0,   0},
                        {a1, a2,  0,   0},
                        {a1, 0,   -a2, 0},
                        {a1, 0,   a2,  0},
                        {a1, 0,   0,   -a2},
                        {a1, 0,   0,   a2},
                        {b1, -b2, -b2, -b2},
                        {b1, -b2, -b2, b2},
                        {b1, -b2, b2,  -b2},
                        {b1, -b2, b2,  b2},
                        {b1, b2,  -b2, -b2},
                        {b1, b2,  -b2, b2},
                        {b1, b2,  b2,  -b2},
                        {b1, b2,  b2,  b2}};
    for (i = 0; i < 6; i++) { for (j = 0; j < 6; j++) { D[i][j] = 0; }}
    for (i = 0; i < 6; i++) { for (j = 0; j < 60; j++) { B[i][j] = 0; }}
    eld_cube(E, u, D);
    for (i = 0; i < 14; i++) {
        det = Inv_jaco_cube(AB[i][1], AB[i][2], AB[i][3], xo, yo, zo, Nx, Ny, Nz, Ele) * AB[i][0];
        for (k = 0; k < 20; k++) {
            B[0][3 * k] = Nx[k];
            B[3][3 * k] = Ny[k];
            B[5][3 * k] = Nz[k];
            B[1][3 * k + 1] = Ny[k];
            B[3][3 * k + 1] = Nx[k];
            B[4][3 * k + 1] = Nz[k];
            B[2][3 * k + 2] = Nz[k];
            B[4][3 * k + 2] = Ny[k];
            B[5][3 * k + 2] = Nx[k];
        }
        for (k = 0; k < 60; k++) {
            for (l = 0; l <= k; l++) {
                for (m = 0; m < 6; m++) {
                    for (n = 0; n < 6; n++) {
                        ek[60 * k + l] += B[m][k] * D[m][n] * B[n][l] * det;
                    }
                }
            }
        }
    }

    for (i = 0; i < 60; i++) {
        for (j = i + 1; j < 60; j++) {
            if (ek[60 * j + i] < 1e-14 && ek[60 * j + i] > -1e-14) { ek[60 * j + i] = 0; }
            ek[60 * i + j] = ek[60 * j + i];
        }
    }
}

