/*
###################################################################################
#
# BCMTools
#
# Copyright (c) 2011-2014 Institute of Industrial Science, The University of Tokyo.
# All rights reserved.
#
# Copyright (c) 2012-2016 Advanced Institute for Computational Science (AICS), RIKEN.
# All rights reserved.
#
# Copyright (c) 2017 Research Institute for Information Technology (RIIT), Kyushu University.
# All rights reserved.
#
###################################################################################
*/

#include "BlockScalar3D.h"

#include "gv.h"

template <>
void BlockScalar3D<real>::ImposeBlockBoundaryCondition_X_M_POISEUILLE_U() {
	::Vec3r c = g_pconf->boundaryValuePoiseuilleCenter;
	real center[3] = {c.x, c.y, c.z};
	bc_x3_poiseuille_u_(this->blockData, &(this->blockBoundaryValue[X_M]), this->size, (int*)&(this->vc), center, this->origin, this->blockSize, this->cellSize);
}
template <>
void BlockScalar3D<real>::ImposeBlockBoundaryCondition_Aw_POISEUILLE_U(real* Ap, real* Aw, real* Ae, real* b) {
	::Vec3r c = g_pconf->boundaryValuePoiseuilleCenter;
	real center[3] = {c.x, c.y, c.z};
	bc_aw_poiseuille_u_(Ap, Aw, b, &(this->blockBoundaryValue[X_M]), this->size, (int*)&(this->vc), center, this->origin, this->blockSize, this->cellSize);
}
template <>
void BlockScalar3D<real>::ImposeBlockBoundaryCondition_X_M_POISEUILLE_P() {
	::Vec3r c = g_pconf->boundaryValuePoiseuilleCenter;
	real center[3] = {c.x, c.y, c.z};
	bc_x3_poiseuille_p_(this->blockData, &(this->blockBoundaryValue[X_M]), this->size, (int*)&(this->vc), center, this->origin, this->blockSize, this->cellSize);
}
template <>
void BlockScalar3D<real>::ImposeBlockBoundaryCondition_Aw_POISEUILLE_P(real* Ap, real* Aw, real* Ae, real* b) {
	::Vec3r c = g_pconf->boundaryValuePoiseuilleCenter;
	real center[3] = {c.x, c.y, c.z};
	bc_aw_poiseuille_p_(Ap, Aw, b, &(this->blockBoundaryValue[X_M]), this->size, (int*)&(this->vc), center, this->origin, this->blockSize, this->cellSize);
}

template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_Dummy() {
}
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_X_M_D() {
}
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_X_P_D() {
}
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_Y_M_D() {
}
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_Y_P_D() {
}
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_Z_M_D() {
}
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_Z_P_D() {
}
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_X_M_N() {
}
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_X_P_N() {
}
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_Y_M_N() {
}
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_Y_P_N() {
}
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_Z_M_N() {
}
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_Z_P_N() {
}
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_X_M_P() {
}
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_X_P_P() {
}
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_Y_M_P() {
}
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_Y_P_P() {
}
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_Z_M_P() {
}
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_Z_P_P() {
}

template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_Dummy(real* Ap, real* Aw, real* Ae, real* b) {
}
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_Aw_D(real* Ap, real* Aw, real* Ae, real* b) {
}
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_Ae_D(real* Ap, real* Aw, real* Ae, real* b) {
}
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_As_D(real* Ap, real* As, real* An, real* b) {
}
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_An_D(real* Ap, real* As, real* An, real* b) {
}
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_Ab_D(real* Ap, real* Ab, real* At, real* b) {
}
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_At_D(real* Ap, real* Ab, real* At, real* b) {
}

template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_Aw_N(real* Ap, real* Aw, real* Ae, real* b) {
}
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_Ae_N(real* Ap, real* Aw, real* Ae, real* b) {
}
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_As_N(real* Ap, real* As, real* An, real* b) {
}
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_An_N(real* Ap, real* As, real* An, real* b) {
}
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_Ab_N(real* Ap, real* Ab, real* At, real* b) {
}
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_At_N(real* Ap, real* Ab, real* At, real* b) {
}

template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_Aw_P(real* Ap, real* Aw, real* Ae, real* b) {
}
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_Ae_P(real* Ap, real* Aw, real* Ae, real* b) {
}
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_As_P(real* Ap, real* As, real* An, real* b) {
}
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_An_P(real* Ap, real* As, real* An, real* b) {
}
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_Ab_P(real* Ap, real* Ab, real* At, real* b) {
}
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_At_P(real* Ap, real* Ab, real* At, real* b) {
}

template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_X_M_POISEUILLE_U() {
	Exit(0);
}
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_Aw_POISEUILLE_U(real* Ap, real* Aw, real* Ae, real* b) {
	Exit(0);
}

template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_X_M_POISEUILLE_P() {
	Exit(0);
}
template <>
void BlockScalar3D<int>::ImposeBlockBoundaryCondition_Aw_POISEUILLE_P(real* Ap, real* Aw, real* Ae, real* b) {
	Exit(0);
}
