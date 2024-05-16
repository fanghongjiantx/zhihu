/*
* @Description: 计算力学(八)总体刚度矩阵和边界条件中代码
* @Author: FangHongjianTX
* @Date: 2024-05-16
* @Email: shijunc@foxmail.com
*
* Copyright (c) 2024 fanghongjiantx
* SPDX-License-Identifier: BSD-3-Clause
*/
#include <iostream>
using namespace std;

/// 刚度矩阵半带宽计算
/// \param jd 长度为ND*NN+1一维整型数组指针，jd[i]为总刚度矩阵中第i行主对角元在一维存储中的地址
/// \param Element NE*8的存储单元的节点编号的二维数组指针，Element[i][j]为第i个单元的第j个节点编号
/// \param NN 节点总数
/// \param ND 每个节点自由度数，所有节点相同
/// \param NE 单元总数
void skdd(int *jd, int(*Element)[8], int NN, int ND, int NE) {
    int i, j, d, *min, *band, min_node;
    min = new int[NN + 1]{0};
    for (i = 1; i <= NN; i++) { min[i] = NN; }
    for (i = 0; i < NE; i++) {
        min_node = Element[i][0];
        for (j = 0; j < 8; j++) {
            if (Element[i][j] > 0 && Element[i][j] < min_node) { min_node = Element[i][j]; }
        }
        for (j = 0; j < 8; j++) {
            if (Element[i][j] > 0 && min_node < min[Element[i][j]]) { min[Element[i][j]] = min_node; }
        }
    }
    band = new int[ND * NN + 1];
    for (i = 1; i <= NN; i++) {
        d = i - min[i];
        for (j = 1; j <= ND; j++) { band[ND * (i - 1) + j] = ND * d + j; }
    }
    jd[1] = 1;
    for (i = 2; i <= ND * NN; i++) { jd[i] = jd[i - 1] + band[i]; }
    jd[0] = 0;
    delete[]min;
    delete[]band;
}

/// 刚度矩阵半带宽计算
/// \param jd 长度为ND*NN+1一维整型数组指针，jd[i]为总刚度矩阵中第i行主对角元在一维存储中的地址
/// \param Element NE*20的存储单元的节点编号的二维数组指针，Element[i][j]为第i个单元的第j个节点编号
/// \param NN 节点总数
/// \param ND 每个节点自由度个数，对所有节点相同
/// \param NE 单元总数
void skdd_cube(int *jd, int(*Element)[20], int NN, int ND, int NE) {
    int i, j, d, *min, *band, min_node;
    min = new int[NN + 1];

    for (i = 1; i <= NN; i++) {
        min[i] = NN;
    }

    for (i = 0; i < NE; i++) {
        min_node = Element[i][0];
        //寻找第i个单元的最小节点号
        for (j = 0; j < 20; j++) {
            if (Element[i][j] > 0 && Element[i][j] < min_node) {
                min_node = Element[i][j];
            }
        }
        //更新第i个单元所有节点的相连的最小节点号
        for (j = 0; j < 20; j++) {
            if (Element[i][j] > 0 && min_node < min[Element[i][j]]) {
                min[Element[i][j]] = min_node;
            }
        }
    }

    band = new int[ND * NN + 1];
    //完成所有半带宽的计算
    for (i = 1; i <= NN; i++) {
        d = i - min[i];
        for (j = 1; j <= ND; j++) {
            band[ND * (i - 1) + j] = ND * d + j;
        }
    }

    jd[1] = 1;

    for (i = 2; i <= ND * NN; i++) {
        jd[i] = jd[i - 1] + band[i];
    }

    jd[0] = 0;
    delete[]min;
    delete[]band;
}

/// 刚度矩阵的调用函数
/// \param i 刚度矩阵元素的行号
/// \param j 刚度矩阵元素的列号
/// \param Value 刚度矩阵元素的值
/// \param Stiff 一维存储的刚度矩阵，长度为jd[NF]+1，Stiff[0]未用
/// \param jd 存储刚度矩阵主对角元在一维存储中的地址的数组，长度为NF(自由度总数)
void SetStiff(int i, int j, double Value, double *Stiff, const int *jd) {
    if (i < j) { return; }
    if (j < (i - jd[i] + jd[i - 1] + 1)) { return; }
    Stiff[jd[i] + j - i] = Value;
}

/// 刚度矩阵的调用函数
/// \param i 刚度矩阵元素的行号
/// \param j 刚度矩阵元素的列号
/// \param Stiff 一维存储的刚度矩阵，长度为jd[NF]+1，Stiff[0]未用
/// \param jd 存储刚度矩阵主对角元在一维存储中的地址的数组，长度为NF(自由度总数)
/// \return 刚度矩阵元素
double GetStiff(int i, int j, double *Stiff, const int *jd) {
    int ii, jj;
    if (i < j) {
        ii = j;
        jj = i;
    } else {
        ii = i;
        jj = j;
    }
    if (jj < ii - jd[ii] + jd[ii - 1] + 1) { return 0; }
    else { return Stiff[jd[ii] + jj - ii]; }
}

/// 单元矩阵的程序
/// \param ek 存放一维按行存储的单元刚度矩阵的数组，长度为EFN*EFN
/// \param je 存放单元自由度总体编号的数组，长度为EFN，若je[i]=0则表明单元中不存在第i个自由度
/// \param EFN 标准单元自由度的个数
/// \param Stiff 一维存储的刚度矩阵，长度为jd[NF]+1,Stiff[0]未用
/// \param jd 存储刚度矩阵主对角元在一维存储中的地址的数组，长度为自由度总数
void Assemble(const double *ek, int *je, int EFN, double *Stiff, int *jd) {
    int ii, ij;
    double value;
    for (ii = 0; ii < EFN; ii++) {
        if (je[ii] <= 0) { continue; }
        for (ij = 0; ij < EFN; ij++) {
            if (je[ij] <= 0) { continue; }
            value = GetStiff(je[ii], je[ij], Stiff, jd) + ek[ii * EFN + ij];
            SetStiff(je[ii], je[ij], value, Stiff, jd);
        }
    }
}



