#include "BlockScalar3D.h"

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

