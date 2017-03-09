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

#ifndef __BCUT_H__
#define __BCUT_H__

#include "real.h"

extern "C" {
	void bcut_calc_c_f_(
				real* fc,
				real* f,
				real* vw, real* ve, real* vs, real* vn, real* vb, real* vt,
				real* c0, real* c1, real* c2, real* c3, real* c4, real* c5,
				int* cid0, int* cid1, int* cid2, int* cid3, int* cid4, int* cid5,
				int* pid,
				real* dx, real* dt,
				real* fg,
				int *sz, int *g);
	void bcut_calc_c_f_c2_(
				real* fc,
				real* f,
				real* vw, real* ve, real* vs, real* vn, real* vb, real* vt,
				real* c0, real* c1, real* c2, real* c3, real* c4, real* c5,
				int* cid0, int* cid1, int* cid2, int* cid3, int* cid4, int* cid5,
				int* pid,
				real* dx, real* dt,
				real* fg,
				int *sz, int *g);
	void bcut_calc_c_f_e3_(
				real* fc,
				real* f,
				real* vw, real* ve, real* vs, real* vn, real* vb, real* vt,
				real* c0, real* c1, real* c2, real* c3, real* c4, real* c5,
				int* cid0, int* cid1, int* cid2, int* cid3, int* cid4, int* cid5,
				int* pid,
				real* dx, real* dt,
				real* fg,
				int *sz, int *g);
	void bcut_calc_c_f_w3_(
				real* fc,
				real* f,
				real* vw, real* ve, real* vs, real* vn, real* vb, real* vt,
				real* c0, real* c1, real* c2, real* c3, real* c4, real* c5,
				int* cid0, int* cid1, int* cid2, int* cid3, int* cid4, int* cid5,
				int* pid,
				real* dx, real* dt,
				real* fg,
				int *sz, int *g);
	void bcut_calc_c_f_u1_(
				real* fc,
				real* f,
				real* vw, real* ve, real* vs, real* vn, real* vb, real* vt,
				real* c0, real* c1, real* c2, real* c3, real* c4, real* c5,
				int* cid0, int* cid1, int* cid2, int* cid3, int* cid4, int* cid5,
				int* pid,
				real* dx, real* dt,
				real* fg,
				int *sz, int *g);
	void bcut_calc_c_f_blend_(
				real* fc,
				real* f,
				real* vw, real* ve, real* vs, real* vn, real* vb, real* vt,
				real* c0, real* c1, real* c2, real* c3, real* c4, real* c5,
				int* cid0, int* cid1, int* cid2, int* cid3, int* cid4, int* cid5,
				int* pid,
				real* dx, real* dt,
				real* fg,
				real* alpha,
				int *sz, int *g);

	void bcut_calc_d_u_(
				real* uxd0_, real* uyd0_, real* uzd0_,
				real* ux0_, real* uy0_, real* uz0_,
				real* c0, real* c1, real* c2, real* c3, real* c4, real* c5,
				int* cid0, int* cid1, int* cid2, int* cid3, int* cid4, int* cid5,
				int* pid,
				real* rho,
				real* mu,
				real* dx, real* dt,
				real* Us,
				int *sz, int *g);

	void bcut_calc_d_t_(
				real* td0_,
				real* t0_,
				real* c0, real* c1, real* c2, real* c3, real* c4, real* c5,
				int* cid0, int* cid1, int* cid2, int* cid3, int* cid4, int* cid5,
				int* pid,
				real* rhof, real* rhos,
				real* cpf, real* cps,
				real* kf, real* ks,
				real* dx, real* dt,
				real* Tc,
				int *sz, int *g);

	void bcut_calc_ab_u_(
				real* Ap, real* Aw, real* Ae, real* As, real* An, real* Ab, real* At, real* b,
				real* ux0_, real* uy0_, real* uz0_,
				real* uc0_, real* ucp_,
				real* ud0_,
				real* p0_,
				real* c0, real* c1, real* c2, real* c3, real* c4, real* c5,
				int* cid0, int* cid1, int* cid2, int* cid3, int* cid4, int* cid5,
				int* pid,
				int* axis,
				real* rhof,
				real* mu,
				real* dx, real* dt,
				real* Us,
				real* gx, real* gy, real* gz,
				int *sz, int *g);

	void bcut_calc_ab_u_e1_(
				real* Ap, real* Aw, real* Ae, real* As, real* An, real* Ab, real* At, real* b,
				real* ux0_, real* uy0_, real* uz0_,
				real* uc0_, real* ucp_,
				real* ud0_,
				real* p0_,
				real* c0, real* c1, real* c2, real* c3, real* c4, real* c5,
				int* cid0, int* cid1, int* cid2, int* cid3, int* cid4, int* cid5,
				int* pid,
				int* axis,
				real* rhof,
				real* mu,
				real* dx, real* dt,
				real* Us,
				int *sz, int *g);

	void bcut_remove_p_(
				real* ux, real* uy, real* uz,
				real* p0_,
				real* c0, real* c1, real* c2, real* c3, real* c4, real* c5,
				int* cid0, int* cid1, int* cid2, int* cid3, int* cid4, int* cid5,
				int* pid,
				real* rho,
				real* dx, real* dt,
				int *sz, int *g);

	void bcut_calc_ab_p_(
				real* Ap, real* Aw, real* Ae, real* As, real* An, real* Ab, real* At, real* b,
				real* vw, real* ve, real* vs, real* vn, real* vb, real* vt,
				real* p0_,
				real* ux, real* uy, real* uz,
				real* c0, real* c1, real* c2, real* c3, real* c4, real* c5,
				int* cid0, int* cid1, int* cid2, int* cid3, int* cid4, int* cid5,
				int* pid,
				real* rho,
				real* dx, real* dt,
				int *sz, int *g);

	void bcut_corr_u_(
				real* ux0_, real* uy0_, real* uz0_,
				real* vw, real* ve, real* vs, real* vn, real* vb, real* vt,
				real* lapp,
				real* p0_,
				real* c0, real* c1, real* c2, real* c3, real* c4, real* c5,
				int* cid0, int* cid1, int* cid2, int* cid3, int* cid4, int* cid5,
				int* pid,
				real* rho,
				real* dx, real* dt,
				int *sz, int *g);

	void bcut_calc_ab_t_1st_(
				real* Ap, real* Aw, real* Ae, real* As, real* An, real* Ab, real* At, real* b,
				real* t0_,
				real* tc0_,
				real* c0, real* c1, real* c2, real* c3, real* c4, real* c5,
				int* cid0, int* cid1, int* cid2, int* cid3, int* cid4, int* cid5,
				int* pid,
				real* rhof, real* rhos,
				real* cpf, real* cps,
				real* kf, real* ks,
				real* dx, real* dt,
				real* Tc,
				int *sz, int *g);

	void bcut_calc_ab_t_(
				real* Ap, real* Aw, real* Ae, real* As, real* An, real* Ab, real* At, real* b,
				real* t0_,
				real* tc0_, real* tcp_,
				real* td0_,
				real* c0, real* c1, real* c2, real* c3, real* c4, real* c5,
				int* cid0, int* cid1, int* cid2, int* cid3, int* cid4, int* cid5,
				int* pid,
				real* rhof, real* rhos,
				real* cpf, real* cps,
				real* kf, real* ks,
				real* dx, real* dt,
				real* Tc,
				int *sz, int *g);

	void bcut_calc_f_p_(
				real *fspx,
				real *fspy,
				real *fspz,
				real *fsp,
				real *p,
				real* c0, real* c1, real* c2, real* c3, real* c4, real* c5,
				int* cid0, int* cid1, int* cid2, int* cid3, int* cid4, int* cid5,
				int* pid,
				real* dx, real* dt,
				int *sz, int *g);

	void bcut_calc_f_v_(
				real *fsvx,
				real *fsvy,
				real *fsvz,
				real *fsv,
				real *ux,
				real *uy,
				real *uz,
				real* c0, real* c1, real* c2, real* c3, real* c4, real* c5,
				int* cid0, int* cid1, int* cid2, int* cid3, int* cid4, int* cid5,
				int* pid,
				real* rho,
				real* mu,
				real* dx, real* dt,
				real* Us,
				int *sz, int *g);

	void bcut_set_fluidseed_(
				int* pid,
				real* xs, real* ys, real *zs,
				real* dx,
				real* org,
				int *sz, int *g);
}

#endif
