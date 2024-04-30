/*
* @Description: 
* @Author: FangHongjianTX
* @Date: 2024-04-30
* @Email: shijunc@foxmail.com
*
* Copyright (c) 2024 fanghongjiantx
* SPDX-License-Identifier: BSD-3-Clause
*/
#include <iostream>
using namespace std;
/// 二维插值基函数的偏导数
void FNA(double A, double B, double *N, const int *Ele) {
    if (Ele[4] > 0) N[4] = -A * (1 - B); else N[4] = 0;
    if (Ele[5] > 0) N[5] = (1 - B * B) / 2; else N[5] = 0;
    if (Ele[6] > 0) N[6] = -A * (1 + B); else N[6] = 0;
    if (Ele[7] > 0) N[7] = -(1 - B * B) / 2; else N[7] = 0;
    N[0] = -(1 - B) / 4 - (N[4] + N[7]) / 2;
    N[1] = (1 - B) / 4 - (N[4] + N[5]) / 2;
    N[2] = (1 + B) / 4 - (N[5] + N[6]) / 2;
    N[3] = -(1 + B) / 4 - (N[6] + N[7]) / 2;
}

void FNB(double A, double B, double *N, const int *Ele) {
    if (Ele[4] > 0) N[4] = -(1 - A * A) / 2; else N[4] = 0;
    if (Ele[5] > 0) N[5] = -(1 + A) * B; else N[5] = 0;
    if (Ele[6] > 0) N[6] = (1 - A * A) / 2; else N[6] = 0;
    if (Ele[7] > 0) N[7] = -(1 - A) * B; else N[7] = 0;
    N[0] = -(1 - A) / 4 - (N[4] + N[7]) / 2;
    N[1] = -(1 + A) / 4 - (N[4] + N[5]) / 2;
    N[2] = +(1 + A) / 4 - (N[5] + N[6]) / 2;
    N[3] = +(1 - A) / 4 - (N[6] + N[7]) / 2;
}


/// 二维三角单元的偏导数
void FNA(double A, double B, double C, double *N, const int *Ele) {
    if (Ele[3] > 0) N[3] = 4 * B; else N[3] = 0;
    N[4] = 0;
    if (Ele[5] > 0) N[5] = 4 * C; else N[5] = 0;
    N[0] = 1 - N[5] / 2 - N[3] / 2;
    N[1] = -N[3] / 2 - N[4] / 2;
    N[2] = -N[5] / 2 - N[4] / 2;
}

void FNB(double A, double B, double C, double *N, const int *Ele) {
    if (Ele[3] > 0) N[3] = 4 * A; else N[3] = 0;
    if (Ele[4] > 0) N[4] = 4 * C; else N[4] = 0;
    N[5] = 0;
    N[0] = -N[5] / 2 - N[3] / 2;
    N[1] = 1 - N[3] / 2 - N[4] / 2;
    N[2] = -N[4] / 2 - N[5] / 2;
}

void FNC(double A, double B, double C, double *N, const int *Ele) {
    N[3] = 0;
    if (Ele[4] > 0) N[4] = 4 * B; else N[4] = 0;
    if (Ele[5] > 0) N[5] = 4 * A; else N[5] = 0;
    N[0] = -N[5] / 2 - N[3] / 2;
    N[1] = -N[3] / 2 - N[4] / 2;
    N[2] = 1 - N[4] / 2 - N[5] / 2;
}


/// 三维插值基函数的偏导数
void FNA_cube(double A, double B, double C, double *N, const int *Ele) {
    if (Ele[8]) N[8] = +(1 - B * B) * (1 - C) / 4; else N[8] = 0;
    if (Ele[9]) N[9] = -A * (1 + B) * (1 - C) / 2; else N[9] = 0;
    if (Ele[10]) N[10] = -(1 - B * B) * (1 - C) / 4; else N[10] = 0;
    if (Ele[11]) N[11] = -A * (1 - B) * (1 - C) / 2; else N[11] = 0;
    if (Ele[12]) N[12] = +(1 - B * B) * (1 + C) / 4; else N[12] = 0;
    if (Ele[13]) N[13] = -A * (1 + B) * (1 + C) / 2; else N[13] = 0;
    if (Ele[14]) N[14] = -(1 - B * B) * (1 + C) / 4; else N[14] = 0;
    if (Ele[15]) N[15] = -A * (1 - B) * (1 + C) / 2; else N[15] = 0;
    if (Ele[16]) N[16] = +(1 - B) * (1 - C * C) / 4; else N[16] = 0;
    if (Ele[17]) N[17] = +(1 + B) * (1 - C * C) / 4; else N[17] = 0;
    if (Ele[18]) N[18] = -(1 + B) * (1 - C * C) / 4; else N[18] = 0;
    if (Ele[19]) N[19] = -(1 - B) * (1 - C * C) / 4; else N[19] = 0;
    N[0] = +(1 - B) * (1 - C) / 8 - (N[8] + N[11] + N[16]) / 2;
    N[1] = +(1 + B) * (1 - C) / 8 - (N[8] + N[9] + N[17]) / 2;
    N[2] = -(1 + B) * (1 - C) / 8 - (N[9] + N[10] + N[18]) / 2;
    N[3] = -(1 - B) * (1 - C) / 8 - (N[10] + N[11] + N[19]) / 2;
    N[4] = +(1 - B) * (1 + C) / 8 - (N[12] + N[15] + N[16]) / 2;
    N[5] = +(1 + B) * (1 + C) / 8 - (N[12] + N[13] + N[17]) / 2;
    N[6] = -(1 + B) * (1 + C) / 8 - (N[13] + N[14] + N[18]) / 2;
    N[7] = -(1 - B) * (1 + C) / 8 - (N[14] + N[15] + N[19]) / 2;
}

void FNB_cube(double A, double B, double C, double *N, const int *Ele) {
    if (Ele[8]) N[8] = -(1 + A) * B * (1 - C) / 2; else N[8] = 0;
    if (Ele[9]) N[9] = (1 - A * A) * (1 - C) / 4; else N[9] = 0;
    if (Ele[10]) N[10] = -(1 - A) * B * (1 - C) / 2; else N[10] = 0;
    if (Ele[11]) N[11] = -(1 - A * A) * (1 - C) / 4; else N[11] = 0;
    if (Ele[12]) N[12] = -(1 + A) * B * (1 + C) / 2; else N[12] = 0;
    if (Ele[13]) N[13] = (1 - A * A) * (1 + C) / 4; else N[13] = 0;
    if (Ele[14]) N[14] = -(1 - A) * B * (1 + C) / 2; else N[14] = 0;
    if (Ele[15]) N[15] = -(1 - A * A) * (1 + C) / 4; else N[15] = 0;
    if (Ele[16]) N[16] = -(1 + A) * (1 - C * C) / 4; else N[16] = 0;
    if (Ele[17]) N[17] = +(1 + A) * (1 - C * C) / 4; else N[17] = 0;
    if (Ele[18]) N[18] = +(1 - A) * (1 - C * C) / 4; else N[18] = 0;
    if (Ele[19]) N[19] = -(1 - A) * (1 - C * C) / 4; else N[19] = 0;
    N[0] = -(1 + A) * (1 - C) / 8 - (N[8] + N[11] + N[16]) / 2;
    N[1] = +(1 + A) * (1 - C) / 8 - (N[8] + N[9] + N[17]) / 2;
    N[2] = +(1 - A) * (1 - C) / 8 - (N[9] + N[10] + N[18]) / 2;
    N[3] = -(1 - A) * (1 - C) / 8 - (N[10] + N[11] + N[19]) / 2;
    N[4] = -(1 + A) * (1 + C) / 8 - (N[12] + N[15] + N[16]) / 2;
    N[5] = +(1 + A) * (1 + C) / 8 - (N[12] + N[13] + N[17]) / 2;
    N[6] = +(1 - A) * (1 + C) / 8 - (N[13] + N[14] + N[18]) / 2;
    N[7] = -(1 - A) * (1 + C) / 8 - (N[14] + N[15] + N[19]) / 2;
}

void FNC_cube(double A, double B, double C, double *N, const int *Ele) {
    if (Ele[8]) N[8] = -(1 + A) * (1 - B * B) / 4; else N[8] = 0;
    if (Ele[9]) N[9] = -(1 - A * A) * (1 + B) / 4; else N[9] = 0;
    if (Ele[10]) N[10] = -(1 - A) * (1 - B * B) / 4; else N[10] = 0;
    if (Ele[11]) N[11] = -(1 - A * A) * (1 - B) / 4; else N[11] = 0;
    if (Ele[12]) N[12] = +(1 + A) * (1 - B * B) / 4; else N[12] = 0;
    if (Ele[13]) N[13] = +(1 - A * A) * (1 + B) / 4; else N[13] = 0;
    if (Ele[14]) N[14] = +(1 - A) * (1 - B * B) / 4; else N[14] = 0;
    if (Ele[15]) N[15] = +(1 - A * A) * (1 - B) / 4; else N[15] = 0;
    if (Ele[16]) N[16] = -(1 + A) * (1 - B) * C / 2; else N[16] = 0;
    if (Ele[17]) N[17] = -(1 + A) * (1 + B) * C / 2; else N[17] = 0;
    if (Ele[18]) N[18] = -(1 - A) * (1 + B) * C / 2; else N[18] = 0;
    if (Ele[19]) N[19] = -(1 - A) * (1 - B) * C / 2; else N[19] = 0;
    N[0] = -(1 + A) * (1 - B) / 8 - (N[8] + N[11] + N[16]) / 2;
    N[1] = -(1 + A) * (1 + B) / 8 - (N[8] + N[9] + N[17]) / 2;
    N[2] = -(1 - A) * (1 + B) / 8 - (N[9] + N[10] + N[18]) / 2;
    N[3] = -(1 - A) * (1 - B) / 8 - (N[10] + N[11] + N[19]) / 2;
    N[4] = +(1 + A) * (1 - B) / 8 - (N[12] + N[15] + N[16]) / 2;
    N[5] = +(1 + A) * (1 + B) / 8 - (N[12] + N[13] + N[17]) / 2;
    N[6] = +(1 - A) * (1 + B) / 8 - (N[13] + N[14] + N[18]) / 2;
    N[7] = +(1 - A) * (1 - B) / 8 - (N[14] + N[15] + N[19]) / 2;
}

/// 雅各比矩阵的行列式
/// \param A 局部坐标
/// \param B 局部坐标
/// \param xo 插值节点总体x坐标，长度为8的数组指针，下同
/// \param yo 插值节点总体y坐标
/// \param Nx 插值基函数对x偏导数
/// \param Ny 插值基函数对y偏导数
/// \param Ele 存放插值区域的节点构成，某元素为0，则意味该节点不存在
/// \return 雅各比矩阵行列式的值
double Inv_jaco(double A, double B, const double *xo, const double *yo, double *Nx, double *Ny, const int *Ele) {
    double jaco11 = 0, jaco12 = 0, jaco21 = 0, jaco22 = 0, det, NA[8] = {0}, NB[8] = {0};
    int i;
    FNA(A, B, NA, Ele);
    FNB(A, B, NB, Ele);
    for (i = 0; i < 8; i++) {
        if (Ele[i] < 0) { continue; }
        jaco11 += xo[i] * NA[i];
        jaco12 += yo[i] * NA[i];
        jaco21 += xo[i] * NB[i];
        jaco22 += yo[i] * NB[i];
    }
    det = jaco11 * jaco22 - jaco21 * jaco12;
    if (det < 0) {
        cout << "雅各比矩阵为负" << endl;
        exit(0);
    }
    for (i = 0; i < 8; i++) {
        if (Ele[i] < 0) {
            Nx[i] = Ny[i] = 0;
            continue;
        }
        Nx[i] = (jaco22 * NA[i] - jaco12 * NB[i]) / det;
        Ny[i] = (-jaco21 * NA[i] + jaco11 * NB[i]) / det;
    }
    return det;
}

double Inv_jaco(double A, double B, double C, const double *xo, const double *yo, double *Nx, double *Ny, int *Ele) {
    double jaco11 = 0, jaco12 = 0, jaco21 = 0, jaco22 = 0, det, NA[6] = {0.0}, NB[6] = {0.0}, NC[6] = {0.0};
    int i;
    FNA(A, B, C, NA, Ele);
    FNB(A, B, C, NB, Ele);
    FNC(A, B, C, NC, Ele);
    for (i = 0; i < 6; i++)
        if (Ele[i] > 0) {
            jaco11 += xo[i] * (NA[i] - NC[i]);
            jaco12 += yo[i] * (NA[i] - NC[i]);
            jaco21 += xo[i] * (NB[i] - NC[i]);
            jaco22 += yo[i] * (NB[i] - NC[i]);
        }
    det = jaco11 * jaco22 - jaco21 * jaco12;
    if (det < 0) {

        cout << "雅各比矩阵为负" << det << endl;

        exit(0);
    }
    for (i = 0; i < 6; i++)
        if (Ele[i] > 0) {
            Nx[i] = (jaco22 * (NA[i] - NC[i]) - jaco12 * (NB[i] - NC[i])) / det;
            Ny[i] = (-jaco21 * (NA[i] - NC[i]) + jaco11 * (NB[i] - NC[i])) / det;
        }
    return det;
}


double Inv_jaco_cube(double A, double B, double C, const double *xo, const double *yo, const double *zo, double *Nx, double *Ny, double *Nz, int *Ele) {
    double jaco11 = 0, jaco12 = 0, jaco13 = 0, jaco21 = 0, jaco22 = 0, jaco23 = 0, jaco31 = 0, jaco32 = 0, jaco33 = 0, det, NA[20] = {
            0}, NB[20] = {0}, NC[20] = {0};
    int i;
    FNA_cube(A, B, C, NA, Ele);
    FNB_cube(A, B, C, NB, Ele);
    FNC_cube(A, B, C, NC, Ele);
    for (int i = 0; i < 20; i++) {
        Nx[i] = Ny[i] = Nz[i] = 0;
    }
    for (i = 0; i < 20; i++) {
        if (!Ele[i])continue;
        jaco11 += xo[i] * NA[i];
        jaco12 += yo[i] * NA[i];
        jaco13 += zo[i] * NA[i];
        jaco21 += xo[i] * NB[i];
        jaco22 += yo[i] * NB[i];
        jaco23 += zo[i] * NB[i];
        jaco31 += xo[i] * NC[i];
        jaco32 += yo[i] * NC[i];
        jaco33 += zo[i] * NC[i];
    }
    det = jaco11 * jaco22 * jaco33 + jaco12 * jaco23 * jaco31 + jaco13 * jaco21 * jaco32 - jaco13 * jaco22 * jaco31 -
          jaco12 * jaco21 * jaco33 - jaco11 * jaco23 * jaco32;
    if (det < 0) {
        cout << "雅各比矩阵为负" << endl;
        exit(0);
    }
    for (i = 0; i < 20; i++) {
        if (!Ele[i]) {
            Nx[i] = 0;
            Ny[i] = 0;
            Nz[i] = 0;
            continue;
        }
        Nx[i] = ((jaco22 * jaco33 - jaco23 * jaco32) * NA[i] - (jaco12 * jaco33 - jaco13 * jaco32) * NB[i] +
                    (jaco12 * jaco23 - jaco13 * jaco22) * NC[i]) / det;
        Ny[i] = (-(jaco12 * jaco33 - jaco13 * jaco31) * NA[i] + (jaco11 * jaco33 - jaco13 * jaco31) * NB[i] -
                    (jaco11 * jaco23 - jaco13 * jaco21) * NC[i]) / det;
        Nz[i] = ((jaco21 * jaco32 - jaco22 * jaco31) * NA[i] - (jaco11 * jaco32 - jaco12 * jaco31) * NB[i] +
                    (jaco11 * jaco22 - jaco12 * jaco21) * NC[i]) / det;
    }
    return det;
}
