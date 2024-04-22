/*
* @Description: 计算力学(五)有限元简介和插值函数中代码
* @Author: FangHongjianTX
* @Date: 2024-04-22
* @Email: shijunc@foxmail.com
*
* Copyright (c) 2024 fanghongjiantx
* SPDX-License-Identifier: BSD-3-Clause
*/

void FN(double A, double B, double *N, const int *Ele) {
    if (Ele[4] > 0) N[4] = (1 - A * A) * (1 - B) / 2; else N[4] = 0;
    if (Ele[5] > 0) N[5] = (1 + A) * (1 - B * B) / 2; else N[5] = 0;
    if (Ele[6] > 0) N[6] = (1 - A * A) * (1 + B) / 2; else N[6] = 0;
    if (Ele[7] > 0) N[7] = (1 - A) * (1 - B * B) / 2; else N[7] = 0;
    N[0] = (1 - A) * (1 - B) / 4 - (N[4] + N[7]) / 2;
    N[1] = (1 + A) * (1 - B) / 4 - (N[4] + N[5]) / 2;
    N[2] = (1 + A) * (1 + B) / 4 - (N[5] + N[6]) / 2;
    N[3] = (1 - A) * (1 + B) / 4 - (N[6] + N[7]) / 2;
}

void FN(double A, double B, double C, double *N, const int *Ele) {
    if (Ele[3] > 0) N[3] = 4 * A * B; else N[3] = 0;
    if (Ele[4] > 0) N[4] = 4 * B * C; else N[4] = 0;
    if (Ele[5] > 0) N[5] = 4 * C * A; else N[5] = 0;
    N[0] = A - N[5] / 2 - N[3] / 2;
    N[1] = B - N[3] / 2 - N[4] / 2;
    N[2] = C - N[4] / 2 - N[5] / 2;
}

void FN_cube(double A, double B, double C, double *N, const int *Ele) {
    if (Ele[8]) N[8] = (1 + A) * (1 - B * B) * (1 - C) / 4; else N[8] = 0;
    if (Ele[9]) N[9] = (1 - A * A) * (1 + B) * (1 - C) / 4; else N[9] = 0;
    if (Ele[10]) N[10] = (1 - A) * (1 - B * B) * (1 - C) / 4; else N[10] = 0;
    if (Ele[11]) N[11] = (1 - A * A) * (1 - B) * (1 - C) / 4; else N[11] = 0;
    if (Ele[12]) N[12] = (1 + A) * (1 - B * B) * (1 + C) / 4; else N[12] = 0;
    if (Ele[13]) N[13] = (1 - A * A) * (1 + B) * (1 + C) / 4; else N[13] = 0;
    if (Ele[14]) N[14] = (1 - A) * (1 - B * B) * (1 + C) / 4; else N[14] = 0;
    if (Ele[15]) N[15] = (1 - A * A) * (1 - B) * (1 + C) / 4; else N[15] = 0;
    if (Ele[16]) N[16] = (1 + A) * (1 - B) * (1 - C * C) / 4; else N[16] = 0;
    if (Ele[17]) N[17] = (1 + A) * (1 + B) * (1 - C * C) / 4; else N[17] = 0;
    if (Ele[18]) N[18] = (1 - A) * (1 + B) * (1 - C * C) / 4; else N[18] = 0;
    if (Ele[19]) N[19] = (1 - A) * (1 - B) * (1 - C * C) / 4; else N[19] = 0;
    N[0] = (1 + A) * (1 - B) * (1 - C) / 8 - (N[8] + N[11] + N[16]) / 2;
    N[1] = (1 + A) * (1 + B) * (1 - C) / 8 - (N[8] + N[9] + N[17]) / 2;
    N[2] = (1 - A) * (1 + B) * (1 - C) / 8 - (N[9] + N[10] + N[18]) / 2;
    N[3] = (1 - A) * (1 - B) * (1 - C) / 8 - (N[10] + N[11] + N[19]) / 2;
    N[4] = (1 + A) * (1 - B) * (1 + C) / 8 - (N[12] + N[15] + N[16]) / 2;
    N[5] = (1 + A) * (1 + B) * (1 + C) / 8 - (N[12] + N[13] + N[17]) / 2;
    N[6] = (1 - A) * (1 + B) * (1 + C) / 8 - (N[13] + N[14] + N[18]) / 2;
    N[7] = (1 - A) * (1 - B) * (1 + C) / 8 - (N[14] + N[15] + N[19]) / 2;
}
