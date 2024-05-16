/*
* @Description: 计算力学(九)方程组求解和后处理
* @Author: FangHongjianTX
* @Date: 2024-05-16
* @Email: shijunc@foxmail.com
*
* Copyright (c) 2024 fanghongjiantx
* SPDX-License-Identifier: BSD-3-Clause
*/
#include <iostream>
using namespace std;

/// 对刚度矩阵制定位移约束
/// \param k 制定位移约束对应的总体自由度编号
/// \param a 制定的位移约束值
/// \param Stiff 一维存储的刚度矩阵，长度为jd[NF]+1，Stiff[0]未用
/// \param Load 存储节点等效载荷的数组,长度为NF+1
/// \param jd 存储刚度矩阵主对角元在以为存储中的地址的数组
/// \param NF 系统的自由度个数
void Fix(int k, double a, double *Stiff, double *Load, int *jd, int NF) {
    int i;
    for (i = 1; i <= NF; i++) { Load[i] -= GetStiff(i, k, Stiff, jd) * a; }
    Load[k] = a;
    for (i = 0; i < NF; i++) {
        SetStiff(i, k, 0, Stiff, jd);
        SetStiff(k, i, 0, Stiff, jd);
    }
    SetStiff(k, k, 1, Stiff, jd);
}

/// 平方根
/// \param jd 储刚度矩阵主对角元在一维存储中的地址，长度为NF(自由度总数)
/// \param zk 一维存储的刚度矩阵,zk[0]未用
/// \param F 存储节点等效载荷,长度NF+1
/// \param NF 系统的自由度总数
void LLT(const int *jd, double *zk, double *F, int NF) {
    int I, L, MI, LI, J, M, LJ, J2, MJ, IJ, J1, K, II, I1;
    double z;
    for (I = 1; I <= NF; I++) {
        L = jd[I];
        MI = jd[I - 1] + I - L + 1;
        LI = L - I;
        for (J = MI; J <= I; J++) {
            M = jd[J];
            LJ = LI + J;
            J2 = J - 1;
            MJ = jd[J2] + J + 1 - M;
            IJ = MJ;
            if (MJ < MI) IJ = MI;
            z = 0;
            J1 = M - J;
            if (IJ <= J2) for (K = IJ; K <= J2; K++) z += zk[K + LI] * zk[K + J1];
            z = zk[LJ] - z;
            if (I == J) {
                if (z < 1.0e-16) cout << "Matrix is no positive determined\n" << "\t" << z << "\t" << I;
                zk[LJ] = sqrt(z);
            } else {
                zk[LJ] = z / zk[M];
            }
        }
    }

    for (I = 1; I <= NF; I++) {
        L = jd[I];
        MI = jd[I - 1] + I - L + 1;
        LI = L - I;
        J1 = I - 1;
        z = 0;
        if (MI <= J1) for (K = MI; K <= J1; K++) z += zk[LI + K] * F[K];
        F[I] = (F[I] - z) / zk[L];
    }
    for (II = 1; II <= NF; II++) {
        I = NF + 1 - II;
        L = jd[I];
        LI = L - I;
        MI = jd[I - 1] - LI + 1;
        z = F[I] / zk[L];
        F[I] = z;
        I1 = I - 1;
        if (MI <= I1) for (J = MI; J <= I1; J++) F[J] = F[J] - zk[LI + J] * z;
    }
}


/// 计算节点应力
/// \param E 杨氏模量
/// \param mu 泊松比
/// \param xo 节点x坐标,长度为8的数组指针
/// \param yo 节点y坐标,长度为8的数组指针
/// \param uo 单元x位移,长度为8的数组指针
/// \param vo 单元y位移,长度为8的数组指针
/// \param SS 节点应力
/// \param Ele 节点编号
/// \param Iopt =0平面应力,=1平面应变,=2轴对称
void CacuSS(double E, double mu, double *xo, double *yo, const double *uo, const double *vo, double(*SS)[5], int *Ele, int Iopt) {
    int i, j, k;
    double Nx[8] = {0}, Ny[8] = {0}, strain[4] = {0}, sum_u = 0, sum_x = 0, A[8] = {-1.0, 1.0, 1.0, -1.0, 0.0, 1.0, 0.0, -1.0},
            B[8] = {-1.0, -1.0, 1.0, 1.0, -1.0, 0.0, 1.0, 0.0}, D[4][4] = {};
    for (i = 0; i < 4; i++) { for (j = 0; j < 4; j++) { D[i][j] = 0; }}
    eld(E, mu, D, Iopt);
    for (i = 0; i < 8; i++)
        if (Ele[i] > 0) {
            SS[Ele[i]][0] += 1;
            Inv_jaco(A[i], B[i], xo, yo, Nx, Ny, Ele);
            for (int i = 0; i < 4; i++) { strain[i] = 0; }
            for (j = 0; j < 8; j++) {
                strain[0] += Nx[j] * uo[j];
                strain[1] += Ny[j] * vo[j];
                strain[2] += Nx[j] * vo[j];
                strain[2] += Ny[j] * uo[j];
                if (Iopt == 2) {
                    if (xo[j] < 1e-12 && xo[j] > -1e-12) { strain[3] = strain[0]; }
                    else { strain[3] += uo[j] / xo[j]; }
                }
            }

            if (Iopt == 2) { strain[3] = sum_u / sum_x; }
            for (k = 0; k < 4; k++) {
                for (j = 0; j < 4; j++) {
                    SS[Ele[i]][k + 1] += D[k][j] * strain[j];
                }
            }
        }
}

/// 计算节点应力
/// \param E 杨氏模量
/// \param mu 泊松比
/// \param xo 节点x坐标,长度为6的数组指针
/// \param yo 节点y坐标,长度为6的数组指针
/// \param uo 单元x位移,长度为6的数组指针
/// \param vo 单元y位移,长度为6的数组指针
/// \param SS 节点应力分量
/// \param Ele 单元节点编号
/// \param Iopt =0平面应力,=1平面应变,=2轴对称
void CacuSS_tri(double E, double mu, double *xo, double *yo, const double *uo, const double *vo, double(*SS)[5], int *Ele, int Iopt) {
    int i, j, k;
    double Nx[6] = {0}, Ny[6] = {0}, strain[4] = {0}, sum_u = 0, sum_x = 0, A[6] = {1, 0, 0, 0.5, 0,
                                                                                    0.5}, B[6] = {0, 1, 0,
                                                                                                  0.5, 0.5,
                                                                                                  0}, C[6] = {
            0, 0, 1, 0, 0.5, 0.5}, D[4][4] = {};
    for (i = 0; i < 4; i++) { for (j = 0; j < 4; j++) { D[i][j] = 0; }}
    eld(E, mu, D, Iopt);
    for (i = 0; i < 6; i++) {
        if (Ele[i] > 0) {
            SS[Ele[i]][0] += 1;
            Inv_jaco(A[i], B[i], C[i], xo, yo, Nx, Ny, Ele);
            for (int i = 0; i < 4; i++) { strain[i] = 0; }
            for (j = 0; j < 6; j++) {
                strain[0] += Nx[j] * uo[j];
                strain[1] += Ny[j] * vo[j];
                strain[2] += Nx[j] * vo[j];
                strain[2] += Ny[j] * uo[j];
                if (Iopt == 2) {
                    if (xo[j] < 1e-12 && xo[j] > -1e-12) { strain[3] = strain[0]; }
                    else { strain[3] += uo[j] / xo[j]; }
                }
            }
            if (Iopt == 2) { strain[3] = sum_u / sum_x; }
            for (k = 0; k < 4; k++) {
                for (j = 0; j < 4; j++) {
                    SS[Ele[i]][k + 1] += D[k][j] * strain[j];
                }
            }
        }
    }
}

/// 计算节点应力
/// \param E 杨氏模量
/// \param mu 泊松比
/// \param xo 节点x坐标,长度为6的数组指针
/// \param yo 节点y坐标,长度为6的数组指针
/// \param zo 节点z坐标,长度为6的数组指针
/// \param uo 单元x位移,长度为6的数组指针
/// \param vo 单元y位移,长度为6的数组指针
/// \param wo 单元z位移,长度为6的数组指针
/// \param SS 节点应力分量
/// \param Ele 单元节点编号
void CacuSS_cube(double E, double mu, double *xo, double *yo, double *zo, const double *uo, const double *vo, const double *wo, double(*SS)[7], int *Ele) {
    int i, j, k;
    double Nx[20] = {0}, Ny[20] = {0}, Nz[20] = {0}, strain[6] = {0}, D[6][6] = {},
            A[20] = {1, 1, -1, -1, 1, 1, -1, -1, 1, 0, -1, 0, 1, 0, -1, 0, 1, 1, -1, -1},
            B[20] = {-1, 1, 1, -1, -1, 1, 1, -1, 0, 1, 0, -1, 0, 1, 0, -1, -1, 1, 1, -1},
            C[20] = {-1, -1, -1, -1, 1, 1, 1, 1, -1, -1, -1, -1, 1, 1, 1, 1, 0, 0, 0, 0};
    for (i = 0; i < 6; i++) { for (j = 0; j < 6; j++) { D[i][j] = 0; }}
    eld_cube(E, mu, D);
    for (i = 0; i < 20; i++) {
        if (Ele[i] < 0) { continue; }
        SS[Ele[i]][0] += 1;
        Inv_jaco_cube(A[i], B[i], C[i], xo, yo, zo, Nx, Ny, Nz, Ele);
        for (j = 0; j < 6; j++) { strain[j] = 0; }
        for (j = 0; j < 6; j++) {
            strain[0] += Nx[j] * uo[j];
            strain[1] += Ny[j] * vo[j];
            strain[2] += Nz[j] * wo[j];
            strain[3] += Nx[j] * vo[j];
            strain[3] += Ny[j] * uo[j];
            strain[4] += Ny[j] * wo[j];
            strain[4] += Nz[j] * vo[j];
            strain[5] += Nz[j] * uo[j];
            strain[5] += Nx[j] * wo[j];
        }
        for (j = 0; j < 4; j++) {
            for (k = 0; k < 6; k++) {
                SS[Ele[i]][k + 1] += D[k][j] * strain[j];
            }
        }
    }
}

/// 计算体积力节点等效载荷
/// \param Density 单元密度
/// \param Thick 单元厚度,仅用于子平面应力情况
/// \param ax 在x方向的惯性加速度,Iopt=2时,则为转动角速度
/// \param ay 在y方向的惯性加速度
/// \param xo 单元x坐标,长度为8的数组指针
/// \param yo 单元y坐标,长度为8的数组指针
/// \param Load 长度为NF+1的数组,存放节点等效载荷
/// \param Ele 存放单元节点编号,长度为8的数组指针
/// \param Iopt =0为平面应力情况，=1为平面应变情况，=2为轴对称情况
/// \param nGauss =3时采用3*3的高斯积分，否测采用2*2的高斯积分
void CacuBF(double Density, double Thick, double ax, double ay, double *xo, double *yo, double *Load, int *Ele, int Iopt, int nGauss) {
    int ij, i, j;
    double AB[3] = {0}, H[3] = {0}, N[8] = {0}, Nx[8] = {0}, Ny[8] = {0}, coff, x, fx;
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
            FN(AB[i], AB[j], N, Ele);
            coff = Density * Inv_jaco(AB[i], AB[j], xo, yo, Nx, Ny, Ele) * H[i] * H[j];
            if (Iopt == 0) {
                coff *= Thick;
                fx = ax;
            } else if (Iopt == 1) { fx = ax; }
            else if (Iopt == 2) {
                x = 0;
                for (int i = 0; i < 8; i++) { x += N[i] * xo[i]; }
                coff *= 2 * PI * x;
                fx = ax * ax * x;
            } else {
                cout << "Please enter the correct Iopt" << endl;
                exit(0);
            }
            for (int i = 0; i < 8; i++) {
                if (Ele[i] < 0) { continue; }
                Load[2 * Ele[i] - 1] += coff * fx * N[i];
                Load[2 * Ele[i]] += coff * ay * N[i];
            }
        }
    }
}

/// 计算体积力节点等效载荷
/// \param Density 单元密度
/// \param Thick 单元厚度,仅用于子平面应力情况
/// \param ax 在x方向的惯性加速度,Iopt=2时,则为转动角速度
/// \param ay 在y方向的惯性加速度
/// \param xo 单元x坐标,长度为6的数组指针
/// \param yo 单元y坐标,长度为6的数组指针
/// \param Load 长度为NF+1的数组指针,存放节点等效载荷
/// \param Ele 存放单元节点编号,长度为6的数组指针
/// \param Iopt =0为平面应力情况，=1为平面应变情况，=2为轴对称情况
void CacuBF(double Density, double Thick, double ax, double ay, double *xo, double *yo, double *Load, int *Ele, int Iopt) {
    double N[8] = {0}, Nx[8] = {0}, Ny[8] = {0}, coff, x, fx;
    double integral[4][4] = {{1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, -27.0 / 48.0},
                             {0.6,       0.2,       0.2,       25.0 / 48.0},
                             {0.2,       0.6,       0.2,       25.0 / 48.0},
                             {0.2,       0.2,       0.6,       25.0 / 48.0}};
    for (int i = 0; i < 4; i++) {
        FN(integral[i][0], integral[i][1], integral[i][2], N, Ele);
        if (Iopt == 2) {
            x = 0;
            for (int i = 0; i < 6; i++) { x += xo[i] * N[i]; }
        }
        coff = Density * Inv_jaco(integral[i][0], integral[i][1], integral[i][2], xo, yo, Nx, Ny, Ele) *
               integral[i][3] / 2.0;
        if (Iopt == 0) {
            coff *= Thick;
            fx = ax;
        } else if (Iopt == 1) { fx = ax; }
        else if (Iopt == 2) {
            coff *= 2 * PI * x;
            fx = ax * ax * x;
        } else {
            cout << "Please enter the correct Iopt" << endl;
            exit(0);
        }
        for (int j = 0; j < 6; j++) {
            if (Ele[j] < 0) { continue; }
            Load[2 * Ele[j] - 1] += coff * fx * N[j];
            Load[2 * Ele[j]] += coff * ay * N[j];
        }
    }
}

/// 计算体积力节点等效载荷
/// \param Density 单元密度
/// \param ax 在x方向的惯性加速度
/// \param ay 在y方向的惯性加速度
/// \param az 在z方向的惯性加速度
/// \param xo 节点x坐标,长度为20的数组指针
/// \param yo 节点y坐标,长度为20的数组指针
/// \param zo 节点y坐标,长度为20的数组指针
/// \param Load 长度为NF+1的数组指针,存放节点等效载荷
/// \param Ele 存放单元节点编号,长度为20的数组指针
void CacuBF_cube(double Density, double ax, double ay, double az, double *xo, double *yo, double *zo, double *Load, int *Ele) {
    int ij, i, j, k;
    double AB[3] = {0}, H[3] = {0}, N[20] = {0}, Nx[20] = {0}, Ny[20] = {0}, Nz[20] = {0}, coff;
    AB[0] = -sqrt(0.6);
    AB[1] = 0;
    AB[2] = -AB[0];
    H[0] = 5.0 / 9.0;
    H[1] = 8.0 / 9.0;
    H[2] = H[0];
    ij = 1;
    for (i = 0; i < 3; i += ij) {
        for (j = 0; j < 3; j += ij) {
            for (k = 0; k < 3; k += ij) {
                FN_cube(AB[i], AB[j], AB[k], N, Ele);
                coff = Density * Inv_jaco_cube(AB[i], AB[j], AB[k], xo, yo, zo, Nx, Ny, Nz, Ele) * H[i] *
                       H[j] * H[k];
                for (int i = 0; i < 20; i++) {
                    if (Ele[i] < 0) { continue; }
                    Load[3 * Ele[i] - 2] += coff * ax * N[i];
                    Load[3 * Ele[i] - 1] += coff * ay * N[i];
                    Load[3 * Ele[i]] += coff * az * N[i];
                }
            }
        }
    }
}

/// 计算压力等效节点载荷
/// \param Thick 单元厚度，仅用于平面应力情况
/// \param po 边界节点压力值，长度为3的数组指针
/// \param xo 边界节点x坐标，长度为3的数组指针
/// \param yo 边界节点y坐标，长度为3的数组指针
/// \param Load 存储节点等效载荷的数组,长度为NF+1
/// \param Side 单元边界节点编号，长度为3的数组指针
/// \param Iopt =0为平面应力情况，=1为平面应变情况，=2为轴对称情况
/// \param nGauss =3时采用3*3的高斯积分，否测采用2*2的高斯积分
void CacuSF(double Thick, const double *po, const double *xo, const double *yo, double *Load, const int *Side, int Iopt, int nGauss) {
    double AB[3] = {0}, H[3] = {0}, qx, qy, coff, x, N[3] = {0}, NA[3] = {0};
    int ij, i, i;
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
    if (Iopt == 0) { coff = Thick; }
    else if (Iopt == 1) { coff = 1; }
    else if (Iopt == 2) { coff = 2 * PI; }
    else {
        cout << "Please enter the correct Iopt" << endl;
        exit(0);
    }
    for (i = 0; i < 3; i += ij) {
        if (Side[2]) {
            N[2] = 1 - AB[i] * AB[i];
            NA[2] = -2 * AB[i];
        } else {
            N[2] = 0;
            NA[2] = 0;
        }
        N[0] = (1 - AB[i] - N[2]) / 2;
        NA[0] = (-1 - NA[2]) / 2;
        N[1] = (1 + AB[i] - N[2]) / 2;
        NA[1] = (1 - NA[2]) / 2;
        if (Iopt == 2) { x = N[0] * xo[0] + N[1] * xo[1] + N[2] * xo[2]; }
        else { x = 1; }
        for (i = 0; i < 3; i++) {
            qx = -(po[0] * N[0] + po[1] * N[1] + po[2] * N[2]) * (NA[0] * yo[0] + NA[1] * yo[1] + NA[2] * yo[2]);
            qy = (po[0] * N[0] + po[1] * N[1] + po[2] * N[2]) * (NA[0] * xo[0] + NA[1] * xo[1] + NA[2] * xo[2]);
            Load[2 * Side[i] - 1] += coff * qx * x * N[i] * H[i];
            Load[2 * Side[i]] += coff * qy * x * N[i] * H[i];
        }
    }
}

