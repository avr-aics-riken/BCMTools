!
! BCMTools
!
! Copyright (C) 2011-2014 Institute of Industrial Science, The University of Tokyo.
! All rights reserved.
!
! Copyright (c) 2012-2015 Advanced Institute for Computational Science, RIKEN.
! All rights reserved.
!

function bcut_getminmod(a, b)
	implicit none
	real					:: bcut_getminmod
	real					:: a, b
!	bcut_getminmod = 0.5*(sign(a, a) + sign(b, b))*min(abs(a), abs(b))
	bcut_getminmod = 0.0d0
	if( abs(a) .lt. abs(b) ) then
		bcut_getminmod = a
	else if( abs(a) .gt. abs(b) ) then
		bcut_getminmod = b
	endif
	return
end function bcut_getminmod

function bcut_getupwind(u, n, p)
	implicit none
	real					:: bcut_getupwind
	real					:: u
	real					:: n, p
	bcut_getupwind = 0.0d0
	if( u .gt. 0 ) then
		bcut_getupwind = n
	else if( u .lt. 0 ) then
		bcut_getupwind = p
	endif
	return
end function bcut_getupwind

function bcut_getweno3(v1, v2, v3)
	implicit none
	real					:: bcut_getweno3
	real					:: v1, v2, v3
	real					:: r
	real					:: w
	real					:: eps
	eps = 0.001
	r = (eps + (v2-v1)**2)/(eps + (v3-v2)**2)
	w = 1.0/(1.0 + 2*r*r)
	bcut_getweno3 = 0.5*(v2 + v3) - 0.5*w*(v1 - 2.0*v2 + v3)
	return
end function bcut_getweno3

subroutine bcut_calc_c_f_c2( &
								fc, &
								f, &
								vw, ve, vs, vn, vb, vt, &
								c0, c1, c2, c3, c4, c5, &
								cid0, cid1, cid2, cid3, cid4, cid5, &
								pid, &
								dx, dt, &
								fi, &
								sz, g)
  implicit none
  integer                 :: i, j, k
  integer                 :: ix, jx, kx
  integer                 :: g
  integer, dimension(3)   :: sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: fc
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: f
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: vw, ve, vs, vn, vb, vt
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: c0, c1, c2, c3, c4, c5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: cid0, cid1, cid2, cid3, cid4, cid5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: pid
	integer									:: cidp
	integer									:: cidp0, cidp1, cidp2, cidp3, cidp4, cidp5
	integer									:: cidw, cide, cids, cidn, cidb, cidt
	integer									:: cidw0, cide1, cids2, cidn3, cidb4, cidt5
	integer									:: pidp, pidw, pide, pids, pidn, pidb, pidt
  real                    :: dx, dt
	real										:: fi
	real										:: d0, d1, d2, d3, d4, d5
	real										:: m0, m1, m2, m3, m4, m5
	real										:: l0, l1, l2, l3, l4, l5
	real										:: vx0, vx1, vy2, vy3, vz4, vz5
	real										:: fp, fw, fe, fs, fn, fb, ft
	real										:: fww, fee, fss, fnn, fbb, ftt
	real										:: f0, f1, f2, f3, f4, f5
	real										:: dfx_p, dfx_c, dfx_n
	real										:: dfy_p, dfy_c, dfy_n
	real										:: dfz_p, dfz_c, dfz_n
	real										:: bcut_getupwind
	real										:: bcut_getminmod
	real										:: q0, q1, q2, q3, q4, q5
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) &
!$omp					 private(cidp) &
!$omp					 private(cidp0, cidp1, cidp2, cidp3, cidp4, cidp5) &
!$omp					 private(cidw, cide, cids, cidn, cidb, cidt) &
!$omp					 private(cidw0, cide1, cids2, cidn3, cidb4, cidt5) &
!$omp					 private(pidp, pidw, pide, pids, pidn, pidb, pidt) &
!$omp					 private(d0, d1, d2, d3, d4, d5) &
!$omp					 private(m0, m1, m2, m3, m4, m5) &
!$omp					 private(l0, l1, l2, l3, l4, l5) &
!$omp					 private(vx0, vx1, vy2, vy3, vz4, vz5) &
!$omp					 private(fp, fw, fe, fs, fn, fb, ft) &
!$omp					 private(fww, fee, fss, fnn, fbb, ftt) &
!$omp					 private(f0, f1, f2, f3, f4, f5) &
!$omp					 private(dfx_p, dfx_c, dfx_n) &
!$omp					 private(dfy_p, dfy_c, dfy_n) &
!$omp					 private(dfz_p, dfz_c, dfz_n) &
!$omp					 private(q0, q1, q2, q3, q4, q5)
!$omp do schedule(static, 1)
#else
#endif
  do k=1, kx
  do j=1, jx
!ocl nouxsimd
  do i=1, ix
		vx0 = vw(i, j, k)
		vx1 = ve(i, j, k)
		vy2 = vs(i, j, k)
		vy3 = vn(i, j, k)
		vz4 = vb(i, j, k)
		vz5 = vt(i, j, k)

		fp = f(i, j, k)
		fw = f(i-1, j, k)
		fe = f(i+1, j, k)
		fs = f(i, j-1, k)
		fn = f(i, j+1, k)
		fb = f(i, j, k-1)
		ft = f(i, j, k+1)

		f0 = 0.5d0*(fp + fw)
		f1 = 0.5d0*(fp + fe)
		f2 = 0.5d0*(fp + fs)
		f3 = 0.5d0*(fp + fn)
		f4 = 0.5d0*(fp + fb)
		f5 = 0.5d0*(fp + ft)

		q0 = vx0*f0
		q1 = vx1*f1
		q2 = vy2*f2
		q3 = vy3*f3
		q4 = vz4*f4
		q5 = vz5*f5

		m0 = 0.0d0
		m1 = 0.0d0
		m2 = 0.0d0
		m3 = 0.0d0
		m4 = 0.0d0
		m5 = 0.0d0

		d0 = c0(i, j, k)
		d1 = c1(i, j, k)
		d2 = c2(i, j, k)
		d3 = c3(i, j, k)
		d4 = c4(i, j, k)
		d5 = c5(i, j, k)

		cidp0 = cid0(i, j, k)
		cidp1 = cid1(i, j, k)
		cidp2 = cid2(i, j, k)
		cidp3 = cid3(i, j, k)
		cidp4 = cid4(i, j, k)
		cidp5 = cid5(i, j, k)

		cidw0 = cid0(i-1, j, k)
		cide1 = cid1(i+1, j, k)
		cids2 = cid2(i, j-1, k)
		cidn3 = cid3(i, j+1, k)
		cidb4 = cid4(i, j, k-1)
		cidt5 = cid5(i, j, k+1)

		pidp = pid(i, j, k)

		if( cidp0 /= 0 ) then
			fw = fe
			fw = (1.0d0 - 1.0d0/d0)*fp + (1.0d0/d0)*fi
			fw = (1.0d0 - 1.0d0/d0)*fp 
			q0 = (d0 - 0.5d0)/(d0 + 0.5d0)*q1
			m0 = 1.0d0
		endif
		if( cidp1 /= 0 ) then
			fe = fw
			fe = (1.0d0 - 1.0d0/d1)*fp + (1.0d0/d1)*fi
			fe = (1.0d0 - 1.0d0/d1)*fp 
			q1 = (d1 - 0.5d0)/(d1 + 0.5d0)*q0
			m1 = 1.0d0
		endif
		if( cidp2 /= 0 ) then
			fs = fn
			fs = (1.0d0 - 1.0d0/d2)*fp + (1.0d0/d2)*fi
			fs = (1.0d0 - 1.0d0/d2)*fp
			q2 = (d2 - 0.5d0)/(d2 + 0.5d0)*q3
			m2 = 1.0d0
		endif
		if( cidp3 /= 0 ) then
			fn = fs
			fn = (1.0d0 - 1.0d0/d3)*fp + (1.0d0/d3)*fi
			fn = (1.0d0 - 1.0d0/d3)*fp
			q3 = (d3 - 0.5d0)/(d3 + 0.5d0)*q2
			m3 = 1.0d0
		endif
		if( cidp4 /= 0 ) then
			fb = ft
			fb = (1.0d0 - 1.0d0/d4)*fp + (1.0d0/d4)*fi
			fb = (1.0d0 - 1.0d0/d4)*fp
			q4 = (d4 - 0.5d0)/(d4 + 0.5d0)*q5
			m4 = 1.0d0
		endif
		if( cidp5 /= 0 ) then
			ft = fb
			ft = (1.0d0 - 1.0d0/d5)*fp + (1.0d0/d5)*fi
			ft = (1.0d0 - 1.0d0/d5)*fp
			q5 = (d5 - 0.5d0)/(d5 + 0.5d0)*q4
			m5 = 1.0d0
		endif

!		f0 = 0.5d0*(fp + fw)
!		f1 = 0.5d0*(fp + fe)
!		f2 = 0.5d0*(fp + fs)
!		f3 = 0.5d0*(fp + fn)
!		f4 = 0.5d0*(fp + fb)
!		f5 = 0.5d0*(fp + ft)
!		fc(i, j, k) = (f1*vx1 - f0*vx0)/dx &
!								+ (f3*vy3 - f2*vy2)/dx &
!								+ (f5*vz5 - f4*vz4)/dx

		fc(i, j, k) = (q1 - q0)/dx &
								+ (q3 - q2)/dx &
								+ (q5 - q4)/dx

!		fc(i, j, k) = (f1*vx1 - f0*vx0)/dx 
!		fc(i, j, k) = 0.0

		if( pidp /= 1 ) then
			fc(i, j, k) = 0.0d0
		endif

  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bcut_calc_c_f_c2

subroutine bcut_calc_c_f_e3( &
								fc, &
								f, &
								vw, ve, vs, vn, vb, vt, &
								c0, c1, c2, c3, c4, c5, &
								cid0, cid1, cid2, cid3, cid4, cid5, &
								pid, &
								dx, dt, &
								fi, &
								sz, g)
  implicit none
  integer                 :: i, j, k
  integer                 :: ix, jx, kx
  integer                 :: g
  integer, dimension(3)   :: sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: fc
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: f
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: vw, ve, vs, vn, vb, vt
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: c0, c1, c2, c3, c4, c5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: cid0, cid1, cid2, cid3, cid4, cid5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: pid
	integer									:: cidp
	integer									:: cidp0, cidp1, cidp2, cidp3, cidp4, cidp5
	integer									:: cidw, cide, cids, cidn, cidb, cidt
	integer									:: cidw0, cide1, cids2, cidn3, cidb4, cidt5
	integer									:: pidp, pidw, pide, pids, pidn, pidb, pidt
  real                    :: dx, dt
	real										:: fi
	real										:: d0, d1, d2, d3, d4, d5
	real										:: m0, m1, m2, m3, m4, m5
	real										:: l0, l1, l2, l3, l4, l5
	real										:: vx0, vx1, vy2, vy3, vz4, vz5
	real										:: fp, fw, fe, fs, fn, fb, ft
	real										:: fww, fee, fss, fnn, fbb, ftt
	real										:: f0, f1, f2, f3, f4, f5
	real										:: dfx_p, dfx_c, dfx_n
	real										:: dfy_p, dfy_c, dfy_n
	real										:: dfz_p, dfz_c, dfz_n
	real										:: bcut_getupwind
	real										:: bcut_getminmod
	real										:: q0, q1, q2, q3, q4, q5
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) &
!$omp					 private(cidp) &
!$omp					 private(cidp0, cidp1, cidp2, cidp3, cidp4, cidp5) &
!$omp					 private(cidw, cide, cids, cidn, cidb, cidt) &
!$omp					 private(cidw0, cide1, cids2, cidn3, cidb4, cidt5) &
!$omp					 private(pidp, pidw, pide, pids, pidn, pidb, pidt) &
!$omp					 private(d0, d1, d2, d3, d4, d5) &
!$omp					 private(m0, m1, m2, m3, m4, m5) &
!$omp					 private(l0, l1, l2, l3, l4, l5) &
!$omp					 private(vx0, vx1, vy2, vy3, vz4, vz5) &
!$omp					 private(fp, fw, fe, fs, fn, fb, ft) &
!$omp					 private(fww, fee, fss, fnn, fbb, ftt) &
!$omp					 private(f0, f1, f2, f3, f4, f5) &
!$omp					 private(dfx_p, dfx_c, dfx_n) &
!$omp					 private(dfy_p, dfy_c, dfy_n) &
!$omp					 private(dfz_p, dfz_c, dfz_n) &
!$omp					 private(q0, q1, q2, q3, q4, q5)
!$omp do schedule(static, 1)
#else
#endif
  do k=1, kx
  do j=1, jx
!ocl nouxsimd
  do i=1, ix
		vx0 = vw(i, j, k)
		vx1 = ve(i, j, k)
		vy2 = vs(i, j, k)
		vy3 = vn(i, j, k)
		vz4 = vb(i, j, k)
		vz5 = vt(i, j, k)

		fp = f(i, j, k)
		fw = f(i-1, j, k)
		fe = f(i+1, j, k)
		fs = f(i, j-1, k)
		fn = f(i, j+1, k)
		fb = f(i, j, k-1)
		ft = f(i, j, k+1)
		fww = f(i-2, j, k)
		fee = f(i+2, j, k)
		fss = f(i, j-2, k)
		fnn = f(i, j+2, k)
		fbb = f(i, j, k-2)
		ftt = f(i, j, k+2)

		m0 = 0.0d0
		m1 = 0.0d0
		m2 = 0.0d0
		m3 = 0.0d0
		m4 = 0.0d0
		m5 = 0.0d0

		d0 = c0(i, j, k)
		d1 = c1(i, j, k)
		d2 = c2(i, j, k)
		d3 = c3(i, j, k)
		d4 = c4(i, j, k)
		d5 = c5(i, j, k)

		cidp0 = cid0(i, j, k)
		cidp1 = cid1(i, j, k)
		cidp2 = cid2(i, j, k)
		cidp3 = cid3(i, j, k)
		cidp4 = cid4(i, j, k)
		cidp5 = cid5(i, j, k)

		cidw0 = cid0(i-1, j, k)
		cide1 = cid1(i+1, j, k)
		cids2 = cid2(i, j-1, k)
		cidn3 = cid3(i, j+1, k)
		cidb4 = cid4(i, j, k-1)
		cidt5 = cid5(i, j, k+1)

		pidp = pid(i, j, k)

		if( cidp0 /= 0 ) then
			fw  = (1.0d0 - 1.0d0/d0)*fp + (1.0d0/d0)*fi
			fww = (1.0d0 - 2.0d0/d0)*fp + (2.0d0/d0)*fi
			m0 = 1.0d0
		endif
		if( cidp1 /= 0 ) then
			fe  = (1.0d0 - 1.0d0/d1)*fp + (1.0d0/d1)*fi
			fee = (1.0d0 - 2.0d0/d1)*fp + (2.0d0/d1)*fi
			m1 = 1.0d0
		endif
		if( cidp2 /= 0 ) then
			fs  = (1.0d0 - 1.0d0/d2)*fp + (1.0d0/d2)*fi
			fss = (1.0d0 - 2.0d0/d2)*fp + (2.0d0/d2)*fi
			m2 = 1.0d0
		endif
		if( cidp3 /= 0 ) then
			fn  = (1.0d0 - 1.0d0/d3)*fp + (1.0d0/d3)*fi
			fnn = (1.0d0 - 2.0d0/d3)*fp + (2.0d0/d3)*fi
			m3 = 1.0d0
		endif
		if( cidp4 /= 0 ) then
			fb  = (1.0d0 - 1.0d0/d4)*fp + (1.0d0/d4)*fi
			fbb = (1.0d0 - 2.0d0/d4)*fp + (2.0d0/d4)*fi
			m4 = 1.0d0
		endif
		if( cidp5 /= 0 ) then
			ft  = (1.0d0 - 1.0d0/d5)*fp + (1.0d0/d5)*fi
			ftt = (1.0d0 - 2.0d0/d5)*fp + (2.0d0/d5)*fi
			m5 = 1.0d0
		endif

		if( cidw0 /= 0 ) then
			fww = fi
		endif
		if( cide1 /= 0 ) then
			fee = fi
		endif
		if( cids2 /= 0 ) then
			fss = fi
		endif
		if( cidn3 /= 0 ) then
			fnn = fi
		endif
		if( cidb4 /= 0 ) then
			fbb = fi
		endif
		if( cidt5 /= 0 ) then
			ftt = fi
		endif

		dfx_p = bcut_getminmod(fee-fe, fe-fp )
		dfx_c = bcut_getminmod(fe -fp, fp-fw )
		dfx_n = bcut_getminmod(fp -fw, fw-fww)
		dfy_p = bcut_getminmod(fnn-fn, fn-fp )
		dfy_c = bcut_getminmod(fn -fp, fp-fs )
		dfy_n = bcut_getminmod(fp -fs, fs-fss)
		dfz_p = bcut_getminmod(ftt-ft, ft-fp )
		dfz_c = bcut_getminmod(ft -fp, fp-fb )
		dfz_n = bcut_getminmod(fp -fb, fb-fbb)
		f0 = bcut_getupwind(vx0, fw + 0.5d0*dfx_n, fp - 0.5d0*dfx_c)
		f1 = bcut_getupwind(vx1, fp + 0.5d0*dfx_c, fe - 0.5d0*dfx_p)
		f2 = bcut_getupwind(vy2, fs + 0.5d0*dfy_n, fp - 0.5d0*dfy_c)
		f3 = bcut_getupwind(vy3, fp + 0.5d0*dfy_c, fn - 0.5d0*dfy_p)
		f4 = bcut_getupwind(vz4, fb + 0.5d0*dfz_n, fp - 0.5d0*dfz_c)
		f5 = bcut_getupwind(vz5, fp + 0.5d0*dfz_c, ft - 0.5d0*dfz_p)

!		f0 = bcut_getupwind(vx0, fw, fp)
!		f1 = bcut_getupwind(vx1, fp, fe)
!		f2 = bcut_getupwind(vy2, fs, fp)
!		f3 = bcut_getupwind(vy3, fp, fn)
!		f4 = bcut_getupwind(vz4, fb, fp)
!		f5 = bcut_getupwind(vz5, fp, ft)

!		f0 = 0.5d0*(fp + fw)
!		f1 = 0.5d0*(fp + fe)
!		f2 = 0.5d0*(fp + fs)
!		f3 = 0.5d0*(fp + fn)
!		f4 = 0.5d0*(fp + fb)
!		f5 = 0.5d0*(fp + ft)

!		f0 = 0.5d0*(fp + fw)*0.95 + bcut_getupwind(vx0, fw, fp)*0.05
!		f1 = 0.5d0*(fp + fe)*0.95 + bcut_getupwind(vx1, fp, fe)*0.05
!		f2 = 0.5d0*(fp + fs)*0.95 + bcut_getupwind(vy2, fs, fp)*0.05
!		f3 = 0.5d0*(fp + fn)*0.95 + bcut_getupwind(vy3, fp, fn)*0.05
!		f4 = 0.5d0*(fp + fb)*0.95 + bcut_getupwind(vz4, fb, fp)*0.05
!		f5 = 0.5d0*(fp + ft)*0.95 + bcut_getupwind(vz5, fp, ft)*0.05

		fc(i, j, k) = (f1*vx1 - f0*vx0)/dx &
								+ (f3*vy3 - f2*vy2)/dx &
								+ (f5*vz5 - f4*vz4)/dx


		q0 = vx0*f0
		q1 = vx1*f1
		q2 = vy2*f2
		q3 = vy3*f3
		q4 = vz4*f4
		q5 = vz5*f5

		if( cidp0 /= 0 ) then
			q0 = (d0 - 0.5d0)/(d0 + 0.5d0)*q1
			m0 = 1.0d0
		endif
		if( cidp1 /= 0 ) then
			q1 = (d1 - 0.5d0)/(d1 + 0.5d0)*q0
			m1 = 1.0d0
		endif
		if( cidp2 /= 0 ) then
			q2 = (d2 - 0.5d0)/(d2 + 0.5d0)*q3
			m2 = 1.0d0
		endif
		if( cidp3 /= 0 ) then
			q3 = (d3 - 0.5d0)/(d3 + 0.5d0)*q2
			m3 = 1.0d0
		endif
		if( cidp4 /= 0 ) then
			q4 = (d4 - 0.5d0)/(d4 + 0.5d0)*q5
			m4 = 1.0d0
		endif
		if( cidp5 /= 0 ) then
			q5 = (d5 - 0.5d0)/(d5 + 0.5d0)*q4
			m5 = 1.0d0
		endif

		fc(i, j, k) = (q1 - q0)/dx &
								+ (q3 - q2)/dx &
								+ (q5 - q4)/dx

!		fc(i, j, k) = (f1*vx1 - f0*vx0)/dx 
!		fc(i, j, k) = 0.0

		if( pidp /= 1 ) then
			fc(i, j, k) = 0.0d0
		endif

  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bcut_calc_c_f_e3

subroutine bcut_calc_c_f_w3( &
								fc, &
								f, &
								vw, ve, vs, vn, vb, vt, &
								c0, c1, c2, c3, c4, c5, &
								cid0, cid1, cid2, cid3, cid4, cid5, &
								pid, &
								dx, dt, &
								fi, &
								sz, g)
  implicit none
  integer                 :: i, j, k
  integer                 :: ix, jx, kx
  integer                 :: g
  integer, dimension(3)   :: sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: fc
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: f
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: vw, ve, vs, vn, vb, vt
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: c0, c1, c2, c3, c4, c5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: cid0, cid1, cid2, cid3, cid4, cid5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: pid
	integer									:: cidp
	integer									:: cidp0, cidp1, cidp2, cidp3, cidp4, cidp5
	integer									:: cidw, cide, cids, cidn, cidb, cidt
	integer									:: cidw0, cide1, cids2, cidn3, cidb4, cidt5
	integer									:: pidp, pidw, pide, pids, pidn, pidb, pidt
  real                    :: dx, dt
	real										:: fi
	real										:: d0, d1, d2, d3, d4, d5
	real										:: m0, m1, m2, m3, m4, m5
	real										:: l0, l1, l2, l3, l4, l5
	real										:: vx0, vx1, vy2, vy3, vz4, vz5
	real										:: fp, fw, fe, fs, fn, fb, ft
	real										:: fww, fee, fss, fnn, fbb, ftt
	real										:: f0, f1, f2, f3, f4, f5
	real										:: dfx_p, dfx_c, dfx_n
	real										:: dfy_p, dfy_c, dfy_n
	real										:: dfz_p, dfz_c, dfz_n
	real										:: bcut_getupwind
	real										:: bcut_getminmod
	real										:: q0, q1, q2, q3, q4, q5
	real										:: vx, vy, vz
	real										:: bcut_getweno3
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) &
!$omp					 private(cidp) &
!$omp					 private(cidp0, cidp1, cidp2, cidp3, cidp4, cidp5) &
!$omp					 private(cidw, cide, cids, cidn, cidb, cidt) &
!$omp					 private(cidw0, cide1, cids2, cidn3, cidb4, cidt5) &
!$omp					 private(pidp, pidw, pide, pids, pidn, pidb, pidt) &
!$omp					 private(d0, d1, d2, d3, d4, d5) &
!$omp					 private(m0, m1, m2, m3, m4, m5) &
!$omp					 private(l0, l1, l2, l3, l4, l5) &
!$omp					 private(vx0, vx1, vy2, vy3, vz4, vz5) &
!$omp					 private(fp, fw, fe, fs, fn, fb, ft) &
!$omp					 private(fww, fee, fss, fnn, fbb, ftt) &
!$omp					 private(f0, f1, f2, f3, f4, f5) &
!$omp					 private(dfx_p, dfx_c, dfx_n) &
!$omp					 private(dfy_p, dfy_c, dfy_n) &
!$omp					 private(dfz_p, dfz_c, dfz_n) &
!$omp					 private(q0, q1, q2, q3, q4, q5) &
!$omp					 private(vx, vy, vz)
!$omp do schedule(static, 1)
#else
#endif
  do k=1, kx
  do j=1, jx
!ocl nouxsimd
  do i=1, ix
		vx0 = vw(i, j, k)
		vx1 = ve(i, j, k)
		vy2 = vs(i, j, k)
		vy3 = vn(i, j, k)
		vz4 = vb(i, j, k)
		vz5 = vt(i, j, k)

		fp = f(i, j, k)
		fw = f(i-1, j, k)
		fe = f(i+1, j, k)
		fs = f(i, j-1, k)
		fn = f(i, j+1, k)
		fb = f(i, j, k-1)
		ft = f(i, j, k+1)
		fww = f(i-2, j, k)
		fee = f(i+2, j, k)
		fss = f(i, j-2, k)
		fnn = f(i, j+2, k)
		fbb = f(i, j, k-2)
		ftt = f(i, j, k+2)

		m0 = 0.0d0
		m1 = 0.0d0
		m2 = 0.0d0
		m3 = 0.0d0
		m4 = 0.0d0
		m5 = 0.0d0

		d0 = c0(i, j, k)
		d1 = c1(i, j, k)
		d2 = c2(i, j, k)
		d3 = c3(i, j, k)
		d4 = c4(i, j, k)
		d5 = c5(i, j, k)

		cidp0 = cid0(i, j, k)
		cidp1 = cid1(i, j, k)
		cidp2 = cid2(i, j, k)
		cidp3 = cid3(i, j, k)
		cidp4 = cid4(i, j, k)
		cidp5 = cid5(i, j, k)

		cidw0 = cid0(i-1, j, k)
		cide1 = cid1(i+1, j, k)
		cids2 = cid2(i, j-1, k)
		cidn3 = cid3(i, j+1, k)
		cidb4 = cid4(i, j, k-1)
		cidt5 = cid5(i, j, k+1)

		pidp = pid(i, j, k)

		if( cidp0 /= 0 ) then
			fw  = (1.0d0 - 1.0d0/d0)*fp + (1.0d0/d0)*fi
			fww = (1.0d0 - 2.0d0/d0)*fp + (2.0d0/d0)*fi
			m0 = 1.0d0
		endif
		if( cidp1 /= 0 ) then
			fe  = (1.0d0 - 1.0d0/d1)*fp + (1.0d0/d1)*fi
			fee = (1.0d0 - 2.0d0/d1)*fp + (2.0d0/d1)*fi
			m1 = 1.0d0
		endif
		if( cidp2 /= 0 ) then
			fs  = (1.0d0 - 1.0d0/d2)*fp + (1.0d0/d2)*fi
			fss = (1.0d0 - 2.0d0/d2)*fp + (2.0d0/d2)*fi
			m2 = 1.0d0
		endif
		if( cidp3 /= 0 ) then
			fn  = (1.0d0 - 1.0d0/d3)*fp + (1.0d0/d3)*fi
			fnn = (1.0d0 - 2.0d0/d3)*fp + (2.0d0/d3)*fi
			m3 = 1.0d0
		endif
		if( cidp4 /= 0 ) then
			fb  = (1.0d0 - 1.0d0/d4)*fp + (1.0d0/d4)*fi
			fbb = (1.0d0 - 2.0d0/d4)*fp + (2.0d0/d4)*fi
			m4 = 1.0d0
		endif
		if( cidp5 /= 0 ) then
			ft  = (1.0d0 - 1.0d0/d5)*fp + (1.0d0/d5)*fi
			ftt = (1.0d0 - 2.0d0/d5)*fp + (2.0d0/d5)*fi
			m5 = 1.0d0
		endif

		if( cidw0 /= 0 ) then
			fww = fi
		endif
		if( cide1 /= 0 ) then
			fee = fi
		endif
		if( cids2 /= 0 ) then
			fss = fi
		endif
		if( cidn3 /= 0 ) then
			fnn = fi
		endif
		if( cidb4 /= 0 ) then
			fbb = fi
		endif
		if( cidt5 /= 0 ) then
			ftt = fi
		endif

		dfx_p = bcut_getweno3(fee - fe, fe - fp, fp - fw)/dx
		dfx_n = bcut_getweno3(fw - fww, fp - fw, fe - fp)/dx
		dfy_p = bcut_getweno3(fnn - fn, fn - fp, fp - fs)/dx
		dfy_n = bcut_getweno3(fs - fss, fp - fs, fn - fp)/dx
		dfz_p = bcut_getweno3(ftt - ft, ft - fp, fp - fb)/dx
		dfz_n = bcut_getweno3(fb - fbb, fp - fb, ft - fp)/dx

		vx = 0.5*(vx0 + vx1)
		vy = 0.5*(vy2 + vy3)
		vz = 0.5*(vz4 + vz5)

		fc(i, j, k) = (max(vx, 0.0)*dfx_n + min(vx, 0.0)*dfx_p) &
								+ (max(vy, 0.0)*dfy_n + min(vy, 0.0)*dfy_p) &
								+ (max(vz, 0.0)*dfz_n + min(vz, 0.0)*dfz_p) 

!		dfx_p = bcut_getminmod(fee-fe, fe-fp )
!		dfx_c = bcut_getminmod(fe -fp, fp-fw )
!		dfx_n = bcut_getminmod(fp -fw, fw-fww)
!		dfy_p = bcut_getminmod(fnn-fn, fn-fp )
!		dfy_c = bcut_getminmod(fn -fp, fp-fs )
!		dfy_n = bcut_getminmod(fp -fs, fs-fss)
!		dfz_p = bcut_getminmod(ftt-ft, ft-fp )
!		dfz_c = bcut_getminmod(ft -fp, fp-fb )
!		dfz_n = bcut_getminmod(fp -fb, fb-fbb)
!		f0 = bcut_getupwind(vx0, fw + 0.5d0*dfx_n, fp - 0.5d0*dfx_c)
!		f1 = bcut_getupwind(vx1, fp + 0.5d0*dfx_c, fe - 0.5d0*dfx_p)
!		f2 = bcut_getupwind(vy2, fs + 0.5d0*dfy_n, fp - 0.5d0*dfy_c)
!		f3 = bcut_getupwind(vy3, fp + 0.5d0*dfy_c, fn - 0.5d0*dfy_p)
!		f4 = bcut_getupwind(vz4, fb + 0.5d0*dfz_n, fp - 0.5d0*dfz_c)
!		f5 = bcut_getupwind(vz5, fp + 0.5d0*dfz_c, ft - 0.5d0*dfz_p)

!		f0 = bcut_getupwind(vx0, fw, fp)
!		f1 = bcut_getupwind(vx1, fp, fe)
!		f2 = bcut_getupwind(vy2, fs, fp)
!		f3 = bcut_getupwind(vy3, fp, fn)
!		f4 = bcut_getupwind(vz4, fb, fp)
!		f5 = bcut_getupwind(vz5, fp, ft)

!		f0 = 0.5d0*(fp + fw)
!		f1 = 0.5d0*(fp + fe)
!		f2 = 0.5d0*(fp + fs)
!		f3 = 0.5d0*(fp + fn)
!		f4 = 0.5d0*(fp + fb)
!		f5 = 0.5d0*(fp + ft)

!		f0 = 0.5d0*(fp + fw)*0.95 + bcut_getupwind(vx0, fw, fp)*0.05
!		f1 = 0.5d0*(fp + fe)*0.95 + bcut_getupwind(vx1, fp, fe)*0.05
!		f2 = 0.5d0*(fp + fs)*0.95 + bcut_getupwind(vy2, fs, fp)*0.05
!		f3 = 0.5d0*(fp + fn)*0.95 + bcut_getupwind(vy3, fp, fn)*0.05
!		f4 = 0.5d0*(fp + fb)*0.95 + bcut_getupwind(vz4, fb, fp)*0.05
!		f5 = 0.5d0*(fp + ft)*0.95 + bcut_getupwind(vz5, fp, ft)*0.05

!		fc(i, j, k) = (f1*vx1 - f0*vx0)/dx &
!								+ (f3*vy3 - f2*vy2)/dx &
!								+ (f5*vz5 - f4*vz4)/dx
!
!
!		q0 = vx0*f0
!		q1 = vx1*f1
!		q2 = vy2*f2
!		q3 = vy3*f3
!		q4 = vz4*f4
!		q5 = vz5*f5
!
!		if( cidp0 /= 0 ) then
!			q0 = (d0 - 0.5d0)/(d0 + 0.5d0)*q1
!			m0 = 1.0d0
!		endif
!		if( cidp1 /= 0 ) then
!			q1 = (d1 - 0.5d0)/(d1 + 0.5d0)*q0
!			m1 = 1.0d0
!		endif
!		if( cidp2 /= 0 ) then
!			q2 = (d2 - 0.5d0)/(d2 + 0.5d0)*q3
!			m2 = 1.0d0
!		endif
!		if( cidp3 /= 0 ) then
!			q3 = (d3 - 0.5d0)/(d3 + 0.5d0)*q2
!			m3 = 1.0d0
!		endif
!		if( cidp4 /= 0 ) then
!			q4 = (d4 - 0.5d0)/(d4 + 0.5d0)*q5
!			m4 = 1.0d0
!		endif
!		if( cidp5 /= 0 ) then
!			q5 = (d5 - 0.5d0)/(d5 + 0.5d0)*q4
!			m5 = 1.0d0
!		endif
!
!		fc(i, j, k) = (q1 - q0)/dx &
!								+ (q3 - q2)/dx &
!								+ (q5 - q4)/dx

!		fc(i, j, k) = (f1*vx1 - f0*vx0)/dx 
!		fc(i, j, k) = 0.0

		if( pidp /= 1 ) then
			fc(i, j, k) = 0.0d0
		endif

  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bcut_calc_c_f_w3

subroutine bcut_calc_c_f_u1( &
								fc, &
								f, &
								vw, ve, vs, vn, vb, vt, &
								c0, c1, c2, c3, c4, c5, &
								cid0, cid1, cid2, cid3, cid4, cid5, &
								pid, &
								dx, dt, &
								fi, &
								sz, g)
  implicit none
  integer                 :: i, j, k
  integer                 :: ix, jx, kx
  integer                 :: g
  integer, dimension(3)   :: sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: fc
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: f
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: vw, ve, vs, vn, vb, vt
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: c0, c1, c2, c3, c4, c5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: cid0, cid1, cid2, cid3, cid4, cid5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: pid
	integer									:: cidp
	integer									:: cidp0, cidp1, cidp2, cidp3, cidp4, cidp5
	integer									:: cidw, cide, cids, cidn, cidb, cidt
	integer									:: cidw0, cide1, cids2, cidn3, cidb4, cidt5
	integer									:: pidp, pidw, pide, pids, pidn, pidb, pidt
  real                    :: dx, dt
	real										:: fi
	real										:: d0, d1, d2, d3, d4, d5
	real										:: m0, m1, m2, m3, m4, m5
	real										:: l0, l1, l2, l3, l4, l5
	real										:: vx0, vx1, vy2, vy3, vz4, vz5
	real										:: fp, fw, fe, fs, fn, fb, ft
	real										:: fww, fee, fss, fnn, fbb, ftt
	real										:: f0, f1, f2, f3, f4, f5
	real										:: dfx_p, dfx_c, dfx_n
	real										:: dfy_p, dfy_c, dfy_n
	real										:: dfz_p, dfz_c, dfz_n
	real										:: bcut_getupwind
	real										:: bcut_getminmod
	real										:: q0, q1, q2, q3, q4, q5
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) &
!$omp					 private(cidp) &
!$omp					 private(cidp0, cidp1, cidp2, cidp3, cidp4, cidp5) &
!$omp					 private(cidw, cide, cids, cidn, cidb, cidt) &
!$omp					 private(cidw0, cide1, cids2, cidn3, cidb4, cidt5) &
!$omp					 private(pidp, pidw, pide, pids, pidn, pidb, pidt) &
!$omp					 private(d0, d1, d2, d3, d4, d5) &
!$omp					 private(m0, m1, m2, m3, m4, m5) &
!$omp					 private(l0, l1, l2, l3, l4, l5) &
!$omp					 private(vx0, vx1, vy2, vy3, vz4, vz5) &
!$omp					 private(fp, fw, fe, fs, fn, fb, ft) &
!$omp					 private(fww, fee, fss, fnn, fbb, ftt) &
!$omp					 private(f0, f1, f2, f3, f4, f5) &
!$omp					 private(dfx_p, dfx_c, dfx_n) &
!$omp					 private(dfy_p, dfy_c, dfy_n) &
!$omp					 private(dfz_p, dfz_c, dfz_n) &
!$omp					 private(q0, q1, q2, q3, q4, q5)
!$omp do schedule(static, 1)
#else
#endif
  do k=1, kx
  do j=1, jx
!ocl nouxsimd
  do i=1, ix
		vx0 = vw(i, j, k)
		vx1 = ve(i, j, k)
		vy2 = vs(i, j, k)
		vy3 = vn(i, j, k)
		vz4 = vb(i, j, k)
		vz5 = vt(i, j, k)

		fp = f(i, j, k)
		fw = f(i-1, j, k)
		fe = f(i+1, j, k)
		fs = f(i, j-1, k)
		fn = f(i, j+1, k)
		fb = f(i, j, k-1)
		ft = f(i, j, k+1)

		m0 = 0.0d0
		m1 = 0.0d0
		m2 = 0.0d0
		m3 = 0.0d0
		m4 = 0.0d0
		m5 = 0.0d0

		d0 = c0(i, j, k)
		d1 = c1(i, j, k)
		d2 = c2(i, j, k)
		d3 = c3(i, j, k)
		d4 = c4(i, j, k)
		d5 = c5(i, j, k)

		cidp0 = cid0(i, j, k)
		cidp1 = cid1(i, j, k)
		cidp2 = cid2(i, j, k)
		cidp3 = cid3(i, j, k)
		cidp4 = cid4(i, j, k)
		cidp5 = cid5(i, j, k)

		cidw0 = cid0(i-1, j, k)
		cide1 = cid1(i+1, j, k)
		cids2 = cid2(i, j-1, k)
		cidn3 = cid3(i, j+1, k)
		cidb4 = cid4(i, j, k-1)
		cidt5 = cid5(i, j, k+1)

		pidp = pid(i, j, k)

		if( cidp0 /= 0 ) then
			fw  = fe
			m0 = 1.0d0
		endif
		if( cidp1 /= 0 ) then
			fe  = fw
			m1 = 1.0d0
		endif
		if( cidp2 /= 0 ) then
			fs  = fn
			m2 = 1.0d0
		endif
		if( cidp3 /= 0 ) then
			fn  = fs
			m3 = 1.0d0
		endif
		if( cidp4 /= 0 ) then
			fb  = ft
			m4 = 1.0d0
		endif
		if( cidp5 /= 0 ) then
			ft  = fb
			m5 = 1.0d0
		endif

		f0 = bcut_getupwind(vx0, fw, fp)
		f1 = bcut_getupwind(vx1, fp, fe)
		f2 = bcut_getupwind(vy2, fs, fp)
		f3 = bcut_getupwind(vy3, fp, fn)
		f4 = bcut_getupwind(vz4, fb, fp)
		f5 = bcut_getupwind(vz5, fp, ft)

		q0 = vx0*f0
		q1 = vx1*f1
		q2 = vy2*f2
		q3 = vy3*f3
		q4 = vz4*f4
		q5 = vz5*f5

		if( cidp0 /= 0 ) then
			q0 = (d0 - 0.5d0)/(d0 + 0.5d0)*q1
			m0 = 1.0d0
		endif
		if( cidp1 /= 0 ) then
			q1 = (d1 - 0.5d0)/(d1 + 0.5d0)*q0
			m1 = 1.0d0
		endif
		if( cidp2 /= 0 ) then
			q2 = (d2 - 0.5d0)/(d2 + 0.5d0)*q3
			m2 = 1.0d0
		endif
		if( cidp3 /= 0 ) then
			q3 = (d3 - 0.5d0)/(d3 + 0.5d0)*q2
			m3 = 1.0d0
		endif
		if( cidp4 /= 0 ) then
			q4 = (d4 - 0.5d0)/(d4 + 0.5d0)*q5
			m4 = 1.0d0
		endif
		if( cidp5 /= 0 ) then
			q5 = (d5 - 0.5d0)/(d5 + 0.5d0)*q4
			m5 = 1.0d0
		endif

		fc(i, j, k) = (q1 - q0)/dx &
								+ (q3 - q2)/dx &
								+ (q5 - q4)/dx

!		fc(i, j, k) = (f1*vx1 - f0*vx0)/dx 
!		fc(i, j, k) = 0.0

		if( pidp /= 1 ) then
			fc(i, j, k) = 0.0d0
		endif

  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bcut_calc_c_f_u1

subroutine bcut_calc_c_f_blend( &
								fc, &
								f, &
								vw, ve, vs, vn, vb, vt, &
								c0, c1, c2, c3, c4, c5, &
								cid0, cid1, cid2, cid3, cid4, cid5, &
								pid, &
								dx, dt, &
								fi, &
								alpha, &
								sz, g)
  implicit none
  integer                 :: i, j, k
  integer                 :: ix, jx, kx
  integer                 :: g
  integer, dimension(3)   :: sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: fc
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: f
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: vw, ve, vs, vn, vb, vt
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: c0, c1, c2, c3, c4, c5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: cid0, cid1, cid2, cid3, cid4, cid5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: pid
	integer									:: cidp
	integer									:: cidp0, cidp1, cidp2, cidp3, cidp4, cidp5
	integer									:: cidw, cide, cids, cidn, cidb, cidt
	integer									:: cidw0, cide1, cids2, cidn3, cidb4, cidt5
	integer									:: pidp, pidw, pide, pids, pidn, pidb, pidt
  real                    :: dx, dt
	real										:: fi
	real										:: alpha
	real										:: d0, d1, d2, d3, d4, d5
	real										:: m0, m1, m2, m3, m4, m5
	real										:: l0, l1, l2, l3, l4, l5
	real										:: vx0, vx1, vy2, vy3, vz4, vz5
	real										:: fp, fw, fe, fs, fn, fb, ft
	real										:: fww, fee, fss, fnn, fbb, ftt
	real										:: f0, f1, f2, f3, f4, f5
	real										:: dfx_p, dfx_c, dfx_n
	real										:: dfy_p, dfy_c, dfy_n
	real										:: dfz_p, dfz_c, dfz_n
	real										:: bcut_getupwind
	real										:: bcut_getminmod
	real										:: q0, q1, q2, q3, q4, q5
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) &
!$omp					 private(cidp) &
!$omp					 private(cidp0, cidp1, cidp2, cidp3, cidp4, cidp5) &
!$omp					 private(cidw, cide, cids, cidn, cidb, cidt) &
!$omp					 private(cidw0, cide1, cids2, cidn3, cidb4, cidt5) &
!$omp					 private(pidp, pidw, pide, pids, pidn, pidb, pidt) &
!$omp					 private(d0, d1, d2, d3, d4, d5) &
!$omp					 private(m0, m1, m2, m3, m4, m5) &
!$omp					 private(l0, l1, l2, l3, l4, l5) &
!$omp					 private(vx0, vx1, vy2, vy3, vz4, vz5) &
!$omp					 private(fp, fw, fe, fs, fn, fb, ft) &
!$omp					 private(fww, fee, fss, fnn, fbb, ftt) &
!$omp					 private(f0, f1, f2, f3, f4, f5) &
!$omp					 private(dfx_p, dfx_c, dfx_n) &
!$omp					 private(dfy_p, dfy_c, dfy_n) &
!$omp					 private(dfz_p, dfz_c, dfz_n) &
!$omp					 private(q0, q1, q2, q3, q4, q5)
!$omp do schedule(static, 1)
#else
#endif
  do k=1, kx
  do j=1, jx
!ocl nouxsimd
  do i=1, ix
		vx0 = vw(i, j, k)
		vx1 = ve(i, j, k)
		vy2 = vs(i, j, k)
		vy3 = vn(i, j, k)
		vz4 = vb(i, j, k)
		vz5 = vt(i, j, k)

		fp = f(i, j, k)
		fw = f(i-1, j, k)
		fe = f(i+1, j, k)
		fs = f(i, j-1, k)
		fn = f(i, j+1, k)
		fb = f(i, j, k-1)
		ft = f(i, j, k+1)

		m0 = 0.0d0
		m1 = 0.0d0
		m2 = 0.0d0
		m3 = 0.0d0
		m4 = 0.0d0
		m5 = 0.0d0

		d0 = c0(i, j, k)
		d1 = c1(i, j, k)
		d2 = c2(i, j, k)
		d3 = c3(i, j, k)
		d4 = c4(i, j, k)
		d5 = c5(i, j, k)

		cidp0 = cid0(i, j, k)
		cidp1 = cid1(i, j, k)
		cidp2 = cid2(i, j, k)
		cidp3 = cid3(i, j, k)
		cidp4 = cid4(i, j, k)
		cidp5 = cid5(i, j, k)

		cidw0 = cid0(i-1, j, k)
		cide1 = cid1(i+1, j, k)
		cids2 = cid2(i, j-1, k)
		cidn3 = cid3(i, j+1, k)
		cidb4 = cid4(i, j, k-1)
		cidt5 = cid5(i, j, k+1)

		pidp = pid(i, j, k)

		if( cidp0 /= 0 ) then
			fw  = fe
			m0 = 1.0d0
		endif
		if( cidp1 /= 0 ) then
			fe  = fw
			m1 = 1.0d0
		endif
		if( cidp2 /= 0 ) then
			fs  = fn
			m2 = 1.0d0
		endif
		if( cidp3 /= 0 ) then
			fn  = fs
			m3 = 1.0d0
		endif
		if( cidp4 /= 0 ) then
			fb  = ft
			m4 = 1.0d0
		endif
		if( cidp5 /= 0 ) then
			ft  = fb
			m5 = 1.0d0
		endif

		f0 = 0.5d0*(fp + fw)*alpha + bcut_getupwind(vx0, fw, fp)*(1.0d0 - alpha)
		f1 = 0.5d0*(fp + fe)*alpha + bcut_getupwind(vx1, fp, fe)*(1.0d0 - alpha)
		f2 = 0.5d0*(fp + fs)*alpha + bcut_getupwind(vy2, fs, fp)*(1.0d0 - alpha)
		f3 = 0.5d0*(fp + fn)*alpha + bcut_getupwind(vy3, fp, fn)*(1.0d0 - alpha)
		f4 = 0.5d0*(fp + fb)*alpha + bcut_getupwind(vz4, fb, fp)*(1.0d0 - alpha)
		f5 = 0.5d0*(fp + ft)*alpha + bcut_getupwind(vz5, fp, ft)*(1.0d0 - alpha)

		q0 = vx0*f0
		q1 = vx1*f1
		q2 = vy2*f2
		q3 = vy3*f3
		q4 = vz4*f4
		q5 = vz5*f5

		if( cidp0 /= 0 ) then
			q0 = (d0 - 0.5d0)/(d0 + 0.5d0)*q1
			m0 = 1.0d0
		endif
		if( cidp1 /= 0 ) then
			q1 = (d1 - 0.5d0)/(d1 + 0.5d0)*q0
			m1 = 1.0d0
		endif
		if( cidp2 /= 0 ) then
			q2 = (d2 - 0.5d0)/(d2 + 0.5d0)*q3
			m2 = 1.0d0
		endif
		if( cidp3 /= 0 ) then
			q3 = (d3 - 0.5d0)/(d3 + 0.5d0)*q2
			m3 = 1.0d0
		endif
		if( cidp4 /= 0 ) then
			q4 = (d4 - 0.5d0)/(d4 + 0.5d0)*q5
			m4 = 1.0d0
		endif
		if( cidp5 /= 0 ) then
			q5 = (d5 - 0.5d0)/(d5 + 0.5d0)*q4
			m5 = 1.0d0
		endif

		fc(i, j, k) = (q1 - q0)/dx &
								+ (q3 - q2)/dx &
								+ (q5 - q4)/dx

!		fc(i, j, k) = (f1*vx1 - f0*vx0)/dx 
!		fc(i, j, k) = 0.0

		if( pidp /= 1 ) then
			fc(i, j, k) = 0.0d0
		endif

  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bcut_calc_c_f_blend


subroutine bcut_calc_d_u( &
								uxd0_, uyd0_, uzd0_, &
								ux0_, uy0_, uz0_, &
								c0, c1, c2, c3, c4, c5, &
								cid0, cid1, cid2, cid3, cid4, cid5, &
								pid, &
								rho, &
								mu, &
								dx, dt, &
								Us, &
								sz, g)
  implicit none
  integer                 :: i, j, k
  integer                 :: ix, jx, kx
  integer                 :: g
  integer, dimension(3)   :: sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: uxd0_, uyd0_, uzd0_
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: ux0_, uy0_, uz0_
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: c0, c1, c2, c3, c4, c5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: cid0, cid1, cid2, cid3, cid4, cid5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: pid
	integer									:: cidp
	integer									:: cidp0, cidp1, cidp2, cidp3, cidp4, cidp5
	integer									:: pidp, pidw, pide, pids, pidn, pidb, pidt
	real										:: rho
	real										:: mu
  real                    :: dx, dt
	real										:: Us
	real										:: d0, d1, d2, d3, d4, d5
	real										:: mu0, mu1, mu2, mu3, mu4, mu5
	real										:: m0, m1, m2, m3, m4, m5
	real										:: l0, l1, l2, l3, l4, l5
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) &
!$omp					 private(cidp) &
!$omp					 private(cidp0, cidp1, cidp2, cidp3, cidp4, cidp5) &
!$omp					 private(pidp, pidw, pide, pids, pidn, pidb, pidt) &
!$omp					 private(d0, d1, d2, d3, d4, d5) &
!$omp					 private(mu0, mu1, mu2, mu3, mu4, mu5) &
!$omp					 private(m0, m1, m2, m3, m4, m5) &
!$omp					 private(l0, l1, l2, l3, l4, l5) 
!$omp do schedule(static, 1)
#else
#endif
  do k=1, kx
  do j=1, jx
!ocl nouxsimd
  do i=1, ix
		mu0 = mu
		mu1 = mu
		mu2 = mu
		mu3 = mu
		mu4 = mu
		mu5 = mu

		m0 = 0.0d0
		m1 = 0.0d0
		m2 = 0.0d0
		m3 = 0.0d0
		m4 = 0.0d0
		m5 = 0.0d0

		d0 = c0(i, j, k)
		d1 = c1(i, j, k)
		d2 = c2(i, j, k)
		d3 = c3(i, j, k)
		d4 = c4(i, j, k)
		d5 = c5(i, j, k)

		cidp0 = cid0(i, j, k)
		cidp1 = cid1(i, j, k)
		cidp2 = cid2(i, j, k)
		cidp3 = cid3(i, j, k)
		cidp4 = cid4(i, j, k)
		cidp5 = cid5(i, j, k)

		pidp = pid(i, j, k)

		if( cidp0 /= 0 ) then
			mu0 = mu/d0*2.0d0/(d0 + d1)
			mu1 = mu/d1*2.0d0/(d0 + d1)
			m0 = 1.0d0
		endif

		if( cidp1 /= 0 ) then
			mu0 = mu/d0*2.0d0/(d0 + d1)
			mu1 = mu/d1*2.0d0/(d0 + d1)
			m1 = 1.0d0
		endif

		if( cidp2 /= 0 ) then
			mu2 = mu/d2*2.0d0/(d2 + d3)
			mu3 = mu/d3*2.0d0/(d2 + d3)
			m2 = 1.0d0
		endif

		if( cidp3 /= 0 ) then
			mu2 = mu/d2*2.0d0/(d2 + d3)
			mu3 = mu/d3*2.0d0/(d2 + d3)
			m3 = 1.0d0
		endif

		if( cidp4 /= 0 ) then
			mu4 = mu/d4*2.0d0/(d4 + d5)
			mu5 = mu/d5*2.0d0/(d4 + d5)
			m4 = 1.0d0
		endif

		if( cidp5 /= 0 ) then
			mu4 = mu/d4*2.0d0/(d4 + d5)
			mu5 = mu/d5*2.0d0/(d4 + d5)
			m5 = 1.0d0
		endif

		l0 = mu0/(rho)/(dx*dx)
		l1 = mu1/(rho)/(dx*dx)
		l2 = mu2/(rho)/(dx*dx)
		l3 = mu3/(rho)/(dx*dx)
		l4 = mu4/(rho)/(dx*dx)
		l5 = mu5/(rho)/(dx*dx)

		uxd0_(i, j, k) = ( &
												l1*(ux0_(i+1, j, k) - ux0_(i, j, k)) &
											- l0*(ux0_(i, j, k) - ux0_(i-1, j, k)) &
											+ l3*(ux0_(i, j+1, k) - ux0_(i, j, k)) &
											- l2*(ux0_(i, j, k) - ux0_(i, j-1, k)) &
											+ l5*(ux0_(i, j, k+1) - ux0_(i, j, k)) &
											- l4*(ux0_(i, j, k) - ux0_(i, j, k-1)) &
											)
		uyd0_(i, j, k) = ( &
												l1*(uy0_(i+1, j, k) - uy0_(i, j, k)) &
											- l0*(uy0_(i, j, k) - uy0_(i-1, j, k)) &
											+ l3*(uy0_(i, j+1, k) - uy0_(i, j, k)) &
											- l2*(uy0_(i, j, k) - uy0_(i, j-1, k)) &
											+ l5*(uy0_(i, j, k+1) - uy0_(i, j, k)) &
											- l4*(uy0_(i, j, k) - uy0_(i, j, k-1)) &
											) 
		uzd0_(i, j, k) = ( &
												l1*(uz0_(i+1, j, k) - uz0_(i, j, k)) &
											- l0*(uz0_(i, j, k) - uz0_(i-1, j, k)) &
											+ l3*(uz0_(i, j+1, k) - uz0_(i, j, k)) &
											- l2*(uz0_(i, j, k) - uz0_(i, j-1, k)) &
											+ l5*(uz0_(i, j, k+1) - uz0_(i, j, k)) &
											- l4*(uz0_(i, j, k) - uz0_(i, j, k-1)) &
											) 

		if( pidp /= 1 ) then
			uxd0_(i, j, k) = 0.0d0
			uyd0_(i, j, k) = 0.0d0
			uzd0_(i, j, k) = 0.0d0
		endif

  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bcut_calc_d_u

subroutine bcut_calc_d_t( &
								td0_, &
								t0_, &
								c0, c1, c2, c3, c4, c5, &
								cid0, cid1, cid2, cid3, cid4, cid5, &
								pid, &
								rhof, rhos, &
								cpf, cps, &
								kf, ks, &
								dx, dt, &
								Tc, &
								sz, g)
  implicit none
  integer                 :: i, j, k
  integer                 :: ix, jx, kx
  integer                 :: g
  integer, dimension(3)   :: sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: td0_
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: t0_
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: c0, c1, c2, c3, c4, c5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: cid0, cid1, cid2, cid3, cid4, cid5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: pid
	integer									:: cidp
	integer									:: cidp0, cidp1, cidp2, cidp3, cidp4, cidp5
	integer									:: pidp, pidw, pide, pids, pidn, pidb, pidt
	real										:: rhof, rhos
	real										:: cpf, cps
	real										:: kf, ks
  real                    :: dx, dt
	real										:: Tc
	real										:: d0, d1, d2, d3, d4, d5
	real										:: k0, k1, k2, k3, k4, k5
	real										:: m0, m1, m2, m3, m4, m5
	real										:: l0, l1, l2, l3, l4, l5
	real										:: tp, tw, te, ts, tn, tb, tt
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) &
!$omp					 private(cidp) &
!$omp					 private(cidp0, cidp1, cidp2, cidp3, cidp4, cidp5) &
!$omp					 private(pidp, pidw, pide, pids, pidn, pidb, pidt) &
!$omp					 private(d0, d1, d2, d3, d4, d5) &
!$omp					 private(k0, k1, k2, k3, k4, k5) &
!$omp					 private(m0, m1, m2, m3, m4, m5) &
!$omp					 private(l0, l1, l2, l3, l4, l5) &
!$omp					 private(tp, tw, te, ts, tn, tb, tt)
!$omp do schedule(static, 1)
#else
#endif
  do k=1, kx
  do j=1, jx
!ocl nouxsimd
  do i=1, ix
		k0 = kf
		k1 = kf
		k2 = kf
		k3 = kf
		k4 = kf
		k5 = kf

		m0 = 0.0d0
		m1 = 0.0d0
		m2 = 0.0d0
		m3 = 0.0d0
		m4 = 0.0d0
		m5 = 0.0d0

		tp = t0_(i, j, k)
		tw = t0_(i-1, j, k)
		te = t0_(i+1, j, k)
		ts = t0_(i, j-1, k)
		tn = t0_(i, j+1, k)
		tb = t0_(i, j, k-1)
		tt = t0_(i, j, k+1)

		d0 = c0(i, j, k)
		d1 = c1(i, j, k)
		d2 = c2(i, j, k)
		d3 = c3(i, j, k)
		d4 = c4(i, j, k)
		d5 = c5(i, j, k)

		cidp0 = cid0(i, j, k)
		cidp1 = cid1(i, j, k)
		cidp2 = cid2(i, j, k)
		cidp3 = cid3(i, j, k)
		cidp4 = cid4(i, j, k)
		cidp5 = cid5(i, j, k)

		pidp = pid(i, j, k)

		if( cidp0 /= 0 ) then
			k0 = kf/d0*2.0/(d0 + d1)
			k1 = kf/d1*2.0/(d0 + d1)
			tw = Tc
			m0 = 1.0
		endif

		if( cidp1 /= 0 ) then
			k0 = kf/d0*2.0/(d0 + d1)
			k1 = kf/d1*2.0/(d0 + d1)
			te = Tc
			m1 = 1.0
		endif

		if( cidp2 /= 0 ) then
			k2 = kf/d2*2.0/(d2 + d3)
			k3 = kf/d3*2.0/(d2 + d3)
			ts = Tc
			m2 = 1.0
		endif

		if( cidp3 /= 0 ) then
			k2 = kf/d2*2.0/(d2 + d3)
			k3 = kf/d3*2.0/(d2 + d3)
			tn = Tc
			m3 = 1.0
		endif

		if( cidp4 /= 0 ) then
			k4 = kf/d4*2.0/(d4 + d5)
			k5 = kf/d5*2.0/(d4 + d5)
			tb = Tc
			m4 = 1.0
		endif

		if( cidp5 /= 0 ) then
			k4 = kf/d4*2.0/(d4 + d5)
			k5 = kf/d5*2.0/(d4 + d5)
			tt = Tc
			m5 = 1.0
		endif

		l0 = k0/(rhof*cpf)*dt/(dx*dx)
		l1 = k1/(rhof*cpf)*dt/(dx*dx)
		l2 = k2/(rhof*cpf)*dt/(dx*dx)
		l3 = k3/(rhof*cpf)*dt/(dx*dx)
		l4 = k4/(rhof*cpf)*dt/(dx*dx)
		l5 = k5/(rhof*cpf)*dt/(dx*dx)

		td0_(i, j, k) = ( &
												l1*(te - tp) &
											- l0*(tp - tw) &
											+ l3*(tn - tp) &
											- l2*(tp - ts) &
											+ l5*(tt - tp) &
											- l4*(tp - tb) &
										)

		if( pidp /= 1 ) then
			td0_(i, j, k) = 0.0d0
		endif

  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bcut_calc_d_t

subroutine bcut_calc_ab_u_1st( &
								Ap, Aw, Ae, As, An, Ab, At, b, &
								u0_, &
								uc0_, &
								c0, c1, c2, c3, c4, c5, &
								cid, &
								pid, &
								rho, &
								mu, &
								dx, dt, &
								Us, &
								sz, g)
  implicit none
  integer                 :: i, j, k
  integer                 :: ix, jx, kx
  integer                 :: g
  integer, dimension(3)   :: sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: Ap, Aw, Ae, As, An, Ab, At, b
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: u0_
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: uc0_
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: c0, c1, c2, c3, c4, c5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: cid
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: pid
	integer									:: cidp
	integer									:: cidp0, cidp1, cidp2, cidp3, cidp4, cidp5
	integer									:: pidp, pidw, pide, pids, pidn, pidb, pidt
	real										:: rho
	real										:: mu
  real                    :: dx, dt
	real										:: Us
	real										:: d0, d1, d2, d3, d4, d5
	real										:: mu0, mu1, mu2, mu3, mu4, mu5
	real										:: m0, m1, m2, m3, m4, m5
	real										:: l0, l1, l2, l3, l4, l5
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) &
!$omp					 private(cidp) &
!$omp					 private(cidp0, cidp1, cidp2, cidp3, cidp4, cidp5) &
!$omp					 private(pidp, pidw, pide, pids, pidn, pidb, pidt) &
!$omp					 private(d0, d1, d2, d3, d4, d5) &
!$omp					 private(mu0, mu1, mu2, mu3, mu4, mu5) &
!$omp					 private(m0, m1, m2, m3, m4, m5) &
!$omp					 private(l0, l1, l2, l3, l4, l5) 
!$omp do schedule(static, 1)
#else
#endif
  do k=1, kx
  do j=1, jx
!ocl nouxsimd
  do i=1, ix
		mu0 = mu
		mu1 = mu
		mu2 = mu
		mu3 = mu
		mu4 = mu
		mu5 = mu

		m0 = 0.0d0
		m1 = 0.0d0
		m2 = 0.0d0
		m3 = 0.0d0
		m4 = 0.0d0
		m5 = 0.0d0

		d0 = c0(i, j, k)
		d1 = c1(i, j, k)
		d2 = c2(i, j, k)
		d3 = c3(i, j, k)
		d4 = c4(i, j, k)
		d5 = c5(i, j, k)

		cidp = cid(i, j, k)
		cidp0 = ibits(cidp,  0, 5)
		cidp1 = ibits(cidp,  5, 5)
		cidp2 = ibits(cidp, 10, 5)
		cidp3 = ibits(cidp, 15, 5)
		cidp4 = ibits(cidp, 20, 5)
		cidp5 = ibits(cidp, 25, 5)

		pidp = pid(i, j, k)

		if( cidp0 /= 0 ) then
			mu0 = mu/d0*2.0d0/(d0 + d1)
			mu1 = mu/d1*2.0d0/(d0 + d1)
			m0 = 1.0d0
		endif

		if( cidp1 /= 0 ) then
			mu0 = mu/d0*2.0d0/(d0 + d1)
			mu1 = mu/d1*2.0d0/(d0 + d1)
			m1 = 1.0d0
		endif

		if( cidp2 /= 0 ) then
			mu2 = mu/d2*2.0d0/(d2 + d3)
			mu3 = mu/d3*2.0d0/(d2 + d3)
			m2 = 1.0d0
		endif

		if( cidp3 /= 0 ) then
			mu2 = mu/d2*2.0d0/(d2 + d3)
			mu3 = mu/d3*2.0d0/(d2 + d3)
			m3 = 1.0d0
		endif

		if( cidp4 /= 0 ) then
			mu4 = mu/d4*2.0d0/(d4 + d5)
			mu5 = mu/d5*2.0d0/(d4 + d5)
			m4 = 1.0d0
		endif

		if( cidp5 /= 0 ) then
			mu4 = mu/d4*2.0d0/(d4 + d5)
			mu5 = mu/d5*2.0d0/(d4 + d5)
			m5 = 1.0d0
		endif

		l0 = mu0/(rho)*dt/(dx*dx)
		l1 = mu1/(rho)*dt/(dx*dx)
		l2 = mu2/(rho)*dt/(dx*dx)
		l3 = mu3/(rho)*dt/(dx*dx)
		l4 = mu4/(rho)*dt/(dx*dx)
		l5 = mu5/(rho)*dt/(dx*dx)

		Ap(i, j, k) = 1.0d0 + (l0 + l1) &
											  + (l2 + l3) &
											  + (l4 + l5)
		Aw(i, j, k) = -l0*(1.0d0 - m0)
		Ae(i, j, k) = -l1*(1.0d0 - m1)
		As(i, j, k) = -l2*(1.0d0 - m2)
		An(i, j, k) = -l3*(1.0d0 - m3)
		Ab(i, j, k) = -l4*(1.0d0 - m4)
		At(i, j, k) = -l5*(1.0d0 - m5)
		 b(i, j, k) = u0_(i, j, k) &
									+ l0*Us*m0 &
									+ l1*Us*m1 &
									+ l2*Us*m2 &
									+ l3*Us*m3 &
									+ l4*Us*m4 &
									+ l5*Us*m5 &
									- uc0_(i, j, k)*dt 

		if( pidp /= 1 ) then
			Ap(i, j, k) = 1.0d0
			Aw(i, j, k) = 0.0d0
			Ae(i, j, k) = 0.0d0
			As(i, j, k) = 0.0d0
			An(i, j, k) = 0.0d0
			Ab(i, j, k) = 0.0d0
			At(i, j, k) = 0.0d0
			b (i, j, k) = Us
		endif

  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bcut_calc_ab_u_1st

subroutine bcut_calc_ab_u_2nd( &
								Ap, Aw, Ae, As, An, Ab, At, b, &
								u0_, &
								uc0_, ucp_, &
								ud0_, &
								p0_, &
								c0, c1, c2, c3, c4, c5, &
								cid, &
								pid, &
								dir, &
								rhof, &
								mu, &
								dx, dt, &
								Us, &
								sz, g)
  implicit none
  integer                 :: i, j, k
  integer                 :: ix, jx, kx
  integer                 :: g
  integer, dimension(3)   :: sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: Ap, Aw, Ae, As, An, Ab, At, b
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: u0_
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: uc0_, ucp_
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: ud0_
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: p0_
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: c0, c1, c2, c3, c4, c5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: cid
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: pid
	integer									:: cidp
	integer									:: cidp0, cidp1, cidp2, cidp3, cidp4, cidp5
	integer									:: pidp, pidw, pide, pids, pidn, pidb, pidt
	integer									:: dir
	real										:: rhof
	real										:: mu
  real                    :: dx, dt
	real										:: Us
	real										:: d0, d1, d2, d3, d4, d5
	real										:: mu0, mu1, mu2, mu3, mu4, mu5
	real										:: m0, m1, m2, m3, m4, m5
	real										:: l0, l1, l2, l3, l4, l5
	real										:: r0, r1, r2, r3, r4, r5
	real										:: dpx0, dpx1, dpy2, dpy3, dpz4, dpz5
	real										:: rdpx0, rdpx1, rdpy2, rdpy3, rdpz4, rdpz5
	real										:: rdpx, rdpy, rdpz, rdp
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) &
!$omp					 private(cidp) &
!$omp					 private(cidp0, cidp1, cidp2, cidp3, cidp4, cidp5) &
!$omp					 private(pidp, pidw, pide, pids, pidn, pidb, pidt) &
!$omp					 private(d0, d1, d2, d3, d4, d5) &
!$omp					 private(mu0, mu1, mu2, mu3, mu4, mu5) &
!$omp					 private(m0, m1, m2, m3, m4, m5) &
!$omp					 private(l0, l1, l2, l3, l4, l5) &
!$omp					 private(r0, r1, r2, r3, r4, r5) &
!$omp					 private(dpx0, dpx1, dpy2, dpy3, dpz4, dpz5) &
!$omp					 private(rdpx0, rdpx1, rdpy2, rdpy3, rdpz4, rdpz5) &
!$omp					 private(rdpx, rdpy, rdpz, rdp)
!$omp do schedule(static, 1)
#else
#endif
  do k=1, kx
  do j=1, jx
!ocl nouxsimd
  do i=1, ix
		mu0 = mu
		mu1 = mu
		mu2 = mu
		mu3 = mu
		mu4 = mu
		mu5 = mu

		r0 = 1.0d0/rhof
		r1 = 1.0d0/rhof
		r2 = 1.0d0/rhof
		r3 = 1.0d0/rhof
		r4 = 1.0d0/rhof
		r5 = 1.0d0/rhof

		dpx0 = (p0_(i, j, k) - p0_(i-1, j, k))/dx
		dpx1 = (p0_(i+1, j, k) - p0_(i, j, k))/dx
		dpy2 = (p0_(i, j, k) - p0_(i, j-1, k))/dx
		dpy3 = (p0_(i, j+1, k) - p0_(i, j, k))/dx
		dpz4 = (p0_(i, j, k) - p0_(i, j, k-1))/dx
		dpz5 = (p0_(i, j, k+1) - p0_(i, j, k))/dx

		rdpx0 = r0*dpx0
		rdpx1 = r1*dpx1
		rdpy2 = r2*dpy2
		rdpy3 = r3*dpy3
		rdpz4 = r4*dpz4
		rdpz5 = r5*dpz5

		m0 = 0.0d0
		m1 = 0.0d0
		m2 = 0.0d0
		m3 = 0.0d0
		m4 = 0.0d0
		m5 = 0.0d0

		d0 = c0(i, j, k)
		d1 = c1(i, j, k)
		d2 = c2(i, j, k)
		d3 = c3(i, j, k)
		d4 = c4(i, j, k)
		d5 = c5(i, j, k)

		cidp = cid(i, j, k)
		cidp0 = ibits(cidp,  0, 5)
		cidp1 = ibits(cidp,  5, 5)
		cidp2 = ibits(cidp, 10, 5)
		cidp3 = ibits(cidp, 15, 5)
		cidp4 = ibits(cidp, 20, 5)
		cidp5 = ibits(cidp, 25, 5)

		pidp = pid(i, j, k)

		if( cidp0 /= 0 ) then
			mu0 = mu/d0*2.0d0/(d0 + d1)
			mu1 = mu/d1*2.0d0/(d0 + d1)
			rdpx0 = rdpx1*(d0 - 0.5d0)/(d0 + 0.5d0)
			dpx0 = dpx1*(d0 - 0.5d0)/(d0 + 0.5d0)
			m0 = 1.0d0
		endif

		if( cidp1 /= 0 ) then
			mu0 = mu/d0*2.0d0/(d0 + d1)
			mu1 = mu/d1*2.0d0/(d0 + d1)
			rdpx1 = rdpx0*(d1 - 0.5d0)/(d1 + 0.5d0)
			dpx1 = dpx0*(d1 - 0.5d0)/(d1 + 0.5d0)
			m1 = 1.0d0
		endif

		if( cidp2 /= 0 ) then
			mu2 = mu/d2*2.0d0/(d2 + d3)
			mu3 = mu/d3*2.0d0/(d2 + d3)
			rdpy2 = rdpy3*(d2 - 0.5d0)/(d2 + 0.5d0)
			dpy2 = dpy3*(d2 - 0.5d0)/(d2 + 0.5d0)
			m2 = 1.0d0
		endif

		if( cidp3 /= 0 ) then
			mu2 = mu/d2*2.0d0/(d2 + d3)
			mu3 = mu/d3*2.0d0/(d2 + d3)
			rdpy3 = rdpy2*(d3 - 0.5d0)/(d3 + 0.5d0)
			dpy3 = dpy2*(d3 - 0.5d0)/(d3 + 0.5d0)
			m3 = 1.0d0
		endif

		if( cidp4 /= 0 ) then
			mu4 = mu/d4*2.0d0/(d4 + d5)
			mu5 = mu/d5*2.0d0/(d4 + d5)
			rdpz4 = rdpz5*(d4 - 0.5d0)/(d4 + 0.5d0)
			dpz4 = dpz5*(d4 - 0.5d0)/(d4 + 0.5d0)
			m4 = 1.0d0
		endif

		if( cidp5 /= 0 ) then
			mu4 = mu/d4*2.0d0/(d4 + d5)
			mu5 = mu/d5*2.0d0/(d4 + d5)
			rdpz5 = rdpz4*(d5 - 0.5d0)/(d5 + 0.5d0)
			dpz5 = dpz4*(d5 - 0.5d0)/(d5 + 0.5d0)
			m5 = 1.0d0
		endif

		if( cidp0 /= 0 .and. cidp1 /= 0 ) then
			rdpx0 = 0.0
			rdpx1 = 0.0
			dpx0 = 0.0
			dpx1 = 0.0
			m0 = 1.0d0
			m1 = 1.0d0
		endif

		if( cidp2 /= 0 .and. cidp3 /= 0 ) then
			rdpy2 = 0.0
			rdpy3 = 0.0
			dpy2 = 0.0
			dpy3 = 0.0
			m2 = 1.0d0
			m3 = 1.0d0
		endif

		if( cidp4 /= 0 .and. cidp5 /= 0 ) then
			rdpz4 = 0.0
			rdpz5 = 0.0
			dpz4 = 0.0
			dpz5 = 0.0
			m4 = 1.0d0
			m5 = 1.0d0
		endif

		l0 = mu0/(rhof)*dt/(dx*dx)
		l1 = mu1/(rhof)*dt/(dx*dx)
		l2 = mu2/(rhof)*dt/(dx*dx)
		l3 = mu3/(rhof)*dt/(dx*dx)
		l4 = mu4/(rhof)*dt/(dx*dx)
		l5 = mu5/(rhof)*dt/(dx*dx)

		rdpx = 0.5d0*(rdpx0 + rdpx1)
		rdpy = 0.5d0*(rdpy2 + rdpy3)
		rdpz = 0.5d0*(rdpz4 + rdpz5)

		rdp = 0.0
		if( dir==0 ) then
			rdp = rdpx
		else if( dir==1 ) then
			rdp = rdpy
		else if( dir==2 ) then
			rdp = rdpz
		end if

		Ap(i, j, k) = 1.0d0 + 0.5d0*(l0 + l1) &
											  + 0.5d0*(l2 + l3) &
											  + 0.5d0*(l4 + l5)
		Aw(i, j, k) = -0.5d0*l0*(1.0d0 - m0)
		Ae(i, j, k) = -0.5d0*l1*(1.0d0 - m1)
		As(i, j, k) = -0.5d0*l2*(1.0d0 - m2)
		An(i, j, k) = -0.5d0*l3*(1.0d0 - m3)
		Ab(i, j, k) = -0.5d0*l4*(1.0d0 - m4)
		At(i, j, k) = -0.5d0*l5*(1.0d0 - m5)
		 b(i, j, k) = u0_(i, j, k) &
									+ 0.5d0*l0*Us*m0 &
									+ 0.5d0*l1*Us*m1 &
									+ 0.5d0*l2*Us*m2 &
									+ 0.5d0*l3*Us*m3 &
									+ 0.5d0*l4*Us*m4 &
									+ 0.5d0*l5*Us*m5 &
									- (1.5d0*uc0_(i, j, k) - 0.5d0*ucp_(i, j, k) )*dt &
									- rdp*dt &
									+ 0.5d0*ud0_(i, j, k)*dt

		if( pidp /= 1 ) then
			Ap(i, j, k) = 1.0d0
			Aw(i, j, k) = 0.0d0
			Ae(i, j, k) = 0.0d0
			As(i, j, k) = 0.0d0
			An(i, j, k) = 0.0d0
			Ab(i, j, k) = 0.0d0
			At(i, j, k) = 0.0d0
			b (i, j, k) = Us
		endif

  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bcut_calc_ab_u_2nd

subroutine bcut_calc_ab_u( &
								Ap, Aw, Ae, As, An, Ab, At, b, &
								ux0_, uy0_, uz0_, &
								uc0_, ucp_, &
								ud0_, &
								p0_, &
								c0, c1, c2, c3, c4, c5, &
								cid0, cid1, cid2, cid3, cid4, cid5, &
								pid, &
								dir, &
								rhof, &
								mu, &
								dx, dt, &
								Us, &
								gx, gy, gz, &
								sz, g)
  implicit none
  integer                 :: i, j, k
  integer                 :: ix, jx, kx
  integer                 :: g
  integer, dimension(3)   :: sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: Ap, Aw, Ae, As, An, Ab, At, b
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: ux0_, uy0_, uz0_
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: uc0_, ucp_
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: ud0_
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: p0_
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: c0, c1, c2, c3, c4, c5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: cid0, cid1, cid2, cid3, cid4, cid5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: pid
	integer									:: cidp
	integer									:: cidp0, cidp1, cidp2, cidp3, cidp4, cidp5
	integer									:: pidp, pidw, pide, pids, pidn, pidb, pidt
	integer									:: dir
	real										:: rhof
	real										:: mu
  real                    :: dx, dt
	real										:: Us
	real										:: gx, gy, gz
	real										:: d0, d1, d2, d3, d4, d5
	real										:: mu0, mu1, mu2, mu3, mu4, mu5
	real										:: m0, m1, m2, m3, m4, m5
	real										:: l0, l1, l2, l3, l4, l5
	real										:: r0, r1, r2, r3, r4, r5
	real										:: dpx0, dpx1, dpy2, dpy3, dpz4, dpz5
	real										:: rdpx0, rdpx1, rdpy2, rdpy3, rdpz4, rdpz5
	real										:: rdpx, rdpy, rdpz, rdp
	real										:: u0, ux0, uy0, uz0
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) &
!$omp					 private(cidp) &
!$omp					 private(cidp0, cidp1, cidp2, cidp3, cidp4, cidp5) &
!$omp					 private(pidp, pidw, pide, pids, pidn, pidb, pidt) &
!$omp					 private(d0, d1, d2, d3, d4, d5) &
!$omp					 private(mu0, mu1, mu2, mu3, mu4, mu5) &
!$omp					 private(m0, m1, m2, m3, m4, m5) &
!$omp					 private(l0, l1, l2, l3, l4, l5) &
!$omp					 private(r0, r1, r2, r3, r4, r5) &
!$omp					 private(dpx0, dpx1, dpy2, dpy3, dpz4, dpz5) &
!$omp					 private(rdpx0, rdpx1, rdpy2, rdpy3, rdpz4, rdpz5) &
!$omp					 private(rdpx, rdpy, rdpz, rdp) &
!$omp					 private(u0, ux0, uy0, uz0) 
!$omp do schedule(static, 1)
#else
#endif
  do k=1, kx
  do j=1, jx
!ocl nouxsimd
  do i=1, ix
		dpx0 = (p0_(i, j, k) - p0_(i-1, j, k))/dx
		dpx1 = (p0_(i+1, j, k) - p0_(i, j, k))/dx
		dpy2 = (p0_(i, j, k) - p0_(i, j-1, k))/dx
		dpy3 = (p0_(i, j+1, k) - p0_(i, j, k))/dx
		dpz4 = (p0_(i, j, k) - p0_(i, j, k-1))/dx
		dpz5 = (p0_(i, j, k+1) - p0_(i, j, k))/dx

		d0 = c0(i, j, k)
		d1 = c1(i, j, k)
		d2 = c2(i, j, k)
		d3 = c3(i, j, k)
		d4 = c4(i, j, k)
		d5 = c5(i, j, k)

		cidp0 = cid0(i, j, k)
		cidp1 = cid1(i, j, k)
		cidp2 = cid2(i, j, k)
		cidp3 = cid3(i, j, k)
		cidp4 = cid4(i, j, k)
		cidp5 = cid5(i, j, k)

		pidp = pid(i, j, k)

    ux0 = ux0_(i, j, k)
    uy0 = uy0_(i, j, k)
    uz0 = uz0_(i, j, k)

		mu0 = mu
		mu1 = mu
		mu2 = mu
		mu3 = mu
		mu4 = mu
		mu5 = mu

		r0 = 1.0d0/rhof
		r1 = 1.0d0/rhof
		r2 = 1.0d0/rhof
		r3 = 1.0d0/rhof
		r4 = 1.0d0/rhof
		r5 = 1.0d0/rhof

		rdpx0 = r0*dpx0 + gx
		rdpx1 = r1*dpx1 + gx
		rdpy2 = r2*dpy2 + gy
		rdpy3 = r3*dpy3 + gy
		rdpz4 = r4*dpz4 + gz
		rdpz5 = r5*dpz5 + gz

		m0 = 0.0d0
		m1 = 0.0d0
		m2 = 0.0d0
		m3 = 0.0d0
		m4 = 0.0d0
		m5 = 0.0d0

		if( cidp0 /= 0 ) then
			mu0 = mu/d0*2.0d0/(d0 + d1)
			mu1 = mu/d1*2.0d0/(d0 + d1)
			rdpx0 = rdpx1*(d0 - 0.5d0)/(d0 + 0.5d0)
			dpx0 = dpx1*(d0 - 0.5d0)/(d0 + 0.5d0)
			m0 = 1.0d0
		endif

		if( cidp1 /= 0 ) then
			mu0 = mu/d0*2.0d0/(d0 + d1)
			mu1 = mu/d1*2.0d0/(d0 + d1)
			rdpx1 = rdpx0*(d1 - 0.5d0)/(d1 + 0.5d0)
			dpx1 = dpx0*(d1 - 0.5d0)/(d1 + 0.5d0)
			m1 = 1.0d0
		endif

		if( cidp2 /= 0 ) then
			mu2 = mu/d2*2.0d0/(d2 + d3)
			mu3 = mu/d3*2.0d0/(d2 + d3)
			rdpy2 = rdpy3*(d2 - 0.5d0)/(d2 + 0.5d0)
			dpy2 = dpy3*(d2 - 0.5d0)/(d2 + 0.5d0)
			m2 = 1.0d0
		endif

		if( cidp3 /= 0 ) then
			mu2 = mu/d2*2.0d0/(d2 + d3)
			mu3 = mu/d3*2.0d0/(d2 + d3)
			rdpy3 = rdpy2*(d3 - 0.5d0)/(d3 + 0.5d0)
			dpy3 = dpy2*(d3 - 0.5d0)/(d3 + 0.5d0)
			m3 = 1.0d0
		endif

		if( cidp4 /= 0 ) then
			mu4 = mu/d4*2.0d0/(d4 + d5)
			mu5 = mu/d5*2.0d0/(d4 + d5)
			rdpz4 = rdpz5*(d4 - 0.5d0)/(d4 + 0.5d0)
			dpz4 = dpz5*(d4 - 0.5d0)/(d4 + 0.5d0)
			m4 = 1.0d0
		endif

		if( cidp5 /= 0 ) then
			mu4 = mu/d4*2.0d0/(d4 + d5)
			mu5 = mu/d5*2.0d0/(d4 + d5)
			rdpz5 = rdpz4*(d5 - 0.5d0)/(d5 + 0.5d0)
			dpz5 = dpz4*(d5 - 0.5d0)/(d5 + 0.5d0)
			m5 = 1.0d0
		endif

		if( cidp0 /= 0 .and. cidp1 /= 0 ) then
			rdpx0 = 0.0
			rdpx1 = 0.0
			dpx0 = 0.0
			dpx1 = 0.0
			m0 = 1.0d0
			m1 = 1.0d0
		endif

		if( cidp2 /= 0 .and. cidp3 /= 0 ) then
			rdpy2 = 0.0
			rdpy3 = 0.0
			dpy2 = 0.0
			dpy3 = 0.0
			m2 = 1.0d0
			m3 = 1.0d0
		endif

		if( cidp4 /= 0 .and. cidp5 /= 0 ) then
			rdpz4 = 0.0
			rdpz5 = 0.0
			dpz4 = 0.0
			dpz5 = 0.0
			m4 = 1.0d0
			m5 = 1.0d0
		endif

		l0 = mu0/(rhof)*dt/(dx*dx)
		l1 = mu1/(rhof)*dt/(dx*dx)
		l2 = mu2/(rhof)*dt/(dx*dx)
		l3 = mu3/(rhof)*dt/(dx*dx)
		l4 = mu4/(rhof)*dt/(dx*dx)
		l5 = mu5/(rhof)*dt/(dx*dx)

		rdpx = 0.5d0*(rdpx0 + rdpx1)
		rdpy = 0.5d0*(rdpy2 + rdpy3)
		rdpz = 0.5d0*(rdpz4 + rdpz5)

		rdp = 0.0
    u0 = 0.0
		if( dir==0 ) then
			rdp = rdpx 
			u0 = ux0
		else if( dir==1 ) then
			rdp = rdpy
			u0 = uy0
		else if( dir==2 ) then
			rdp = rdpz
			u0 = uz0
		end if

		Ap(i, j, k) = 1.0d0 + 0.5d0*(l0 + l1) &
											  + 0.5d0*(l2 + l3) &
											  + 0.5d0*(l4 + l5)
		Aw(i, j, k) = -0.5d0*l0*(1.0d0 - m0)
		Ae(i, j, k) = -0.5d0*l1*(1.0d0 - m1)
		As(i, j, k) = -0.5d0*l2*(1.0d0 - m2)
		An(i, j, k) = -0.5d0*l3*(1.0d0 - m3)
		Ab(i, j, k) = -0.5d0*l4*(1.0d0 - m4)
		At(i, j, k) = -0.5d0*l5*(1.0d0 - m5)
		 b(i, j, k) = u0 &
									+ 0.5d0*l0*Us*m0 &
									+ 0.5d0*l1*Us*m1 &
									+ 0.5d0*l2*Us*m2 &
									+ 0.5d0*l3*Us*m3 &
									+ 0.5d0*l4*Us*m4 &
									+ 0.5d0*l5*Us*m5 &
									- (1.5d0*uc0_(i, j, k) - 0.5d0*ucp_(i, j, k) )*dt &
									- rdp*dt &
									+ 0.5d0*ud0_(i, j, k)*dt

		if( pidp /= 1 ) then
			Ap(i, j, k) = 1.0d0
			Aw(i, j, k) = 0.0d0
			Ae(i, j, k) = 0.0d0
			As(i, j, k) = 0.0d0
			An(i, j, k) = 0.0d0
			Ab(i, j, k) = 0.0d0
			At(i, j, k) = 0.0d0
			b (i, j, k) = Us
		endif

  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bcut_calc_ab_u

subroutine bcut_calc_ab_u_e1( &
								Ap, Aw, Ae, As, An, Ab, At, b, &
								ux0_, uy0_, uz0_, &
								uc0_, ucp_, &
								ud0_, &
								p0_, &
								c0, c1, c2, c3, c4, c5, &
								cid0, cid1, cid2, cid3, cid4, cid5, &
								pid, &
								dir, &
								rhof, &
								mu, &
								dx, dt, &
								Us, &
								sz, g)
  implicit none
  integer                 :: i, j, k
  integer                 :: ix, jx, kx
  integer                 :: g
  integer, dimension(3)   :: sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: Ap, Aw, Ae, As, An, Ab, At, b
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: ux0_, uy0_, uz0_
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: uc0_, ucp_
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: ud0_
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: p0_
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: c0, c1, c2, c3, c4, c5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: cid0, cid1, cid2, cid3, cid4, cid5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: pid
	integer									:: cidp
	integer									:: cidp0, cidp1, cidp2, cidp3, cidp4, cidp5
	integer									:: pidp, pidw, pide, pids, pidn, pidb, pidt
	integer									:: dir
	real										:: rhof
	real										:: mu
  real                    :: dx, dt
	real										:: Us
	real										:: d0, d1, d2, d3, d4, d5
	real										:: mu0, mu1, mu2, mu3, mu4, mu5
	real										:: m0, m1, m2, m3, m4, m5
	real										:: l0, l1, l2, l3, l4, l5
	real										:: r0, r1, r2, r3, r4, r5
	real										:: dpx0, dpx1, dpy2, dpy3, dpz4, dpz5
	real										:: rdpx0, rdpx1, rdpy2, rdpy3, rdpz4, rdpz5
	real										:: rdpx, rdpy, rdpz, rdp
	real										:: u0, ux0, uy0, uz0
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) &
!$omp					 private(cidp) &
!$omp					 private(cidp0, cidp1, cidp2, cidp3, cidp4, cidp5) &
!$omp					 private(pidp, pidw, pide, pids, pidn, pidb, pidt) &
!$omp					 private(d0, d1, d2, d3, d4, d5) &
!$omp					 private(mu0, mu1, mu2, mu3, mu4, mu5) &
!$omp					 private(m0, m1, m2, m3, m4, m5) &
!$omp					 private(l0, l1, l2, l3, l4, l5) &
!$omp					 private(r0, r1, r2, r3, r4, r5) &
!$omp					 private(dpx0, dpx1, dpy2, dpy3, dpz4, dpz5) &
!$omp					 private(rdpx0, rdpx1, rdpy2, rdpy3, rdpz4, rdpz5) &
!$omp					 private(rdpx, rdpy, rdpz, rdp) &
!$omp					 private(u0, ux0, uy0, uz0) 
!$omp do schedule(static, 1)
#else
#endif
  do k=1, kx
  do j=1, jx
!ocl nouxsimd
  do i=1, ix
		dpx0 = (p0_(i, j, k) - p0_(i-1, j, k))/dx
		dpx1 = (p0_(i+1, j, k) - p0_(i, j, k))/dx
		dpy2 = (p0_(i, j, k) - p0_(i, j-1, k))/dx
		dpy3 = (p0_(i, j+1, k) - p0_(i, j, k))/dx
		dpz4 = (p0_(i, j, k) - p0_(i, j, k-1))/dx
		dpz5 = (p0_(i, j, k+1) - p0_(i, j, k))/dx

		d0 = c0(i, j, k)
		d1 = c1(i, j, k)
		d2 = c2(i, j, k)
		d3 = c3(i, j, k)
		d4 = c4(i, j, k)
		d5 = c5(i, j, k)

		cidp0 = cid0(i, j, k)
		cidp1 = cid1(i, j, k)
		cidp2 = cid2(i, j, k)
		cidp3 = cid3(i, j, k)
		cidp4 = cid4(i, j, k)
		cidp5 = cid5(i, j, k)

		pidp = pid(i, j, k)

    ux0 = ux0_(i, j, k)
    uy0 = uy0_(i, j, k)
    uz0 = uz0_(i, j, k)

		mu0 = mu
		mu1 = mu
		mu2 = mu
		mu3 = mu
		mu4 = mu
		mu5 = mu

		r0 = 1.0d0/rhof
		r1 = 1.0d0/rhof
		r2 = 1.0d0/rhof
		r3 = 1.0d0/rhof
		r4 = 1.0d0/rhof
		r5 = 1.0d0/rhof

		rdpx0 = r0*dpx0
		rdpx1 = r1*dpx1
		rdpy2 = r2*dpy2
		rdpy3 = r3*dpy3
		rdpz4 = r4*dpz4
		rdpz5 = r5*dpz5

		m0 = 0.0d0
		m1 = 0.0d0
		m2 = 0.0d0
		m3 = 0.0d0
		m4 = 0.0d0
		m5 = 0.0d0

		if( cidp0 /= 0 ) then
			mu0 = mu/d0*2.0d0/(d0 + d1)
			mu1 = mu/d1*2.0d0/(d0 + d1)
			rdpx0 = rdpx1*(d0 - 0.5d0)/(d0 + 0.5d0)
			dpx0 = dpx1*(d0 - 0.5d0)/(d0 + 0.5d0)
			m0 = 1.0d0
		endif

		if( cidp1 /= 0 ) then
			mu0 = mu/d0*2.0d0/(d0 + d1)
			mu1 = mu/d1*2.0d0/(d0 + d1)
			rdpx1 = rdpx0*(d1 - 0.5d0)/(d1 + 0.5d0)
			dpx1 = dpx0*(d1 - 0.5d0)/(d1 + 0.5d0)
			m1 = 1.0d0
		endif

		if( cidp2 /= 0 ) then
			mu2 = mu/d2*2.0d0/(d2 + d3)
			mu3 = mu/d3*2.0d0/(d2 + d3)
			rdpy2 = rdpy3*(d2 - 0.5d0)/(d2 + 0.5d0)
			dpy2 = dpy3*(d2 - 0.5d0)/(d2 + 0.5d0)
			m2 = 1.0d0
		endif

		if( cidp3 /= 0 ) then
			mu2 = mu/d2*2.0d0/(d2 + d3)
			mu3 = mu/d3*2.0d0/(d2 + d3)
			rdpy3 = rdpy2*(d3 - 0.5d0)/(d3 + 0.5d0)
			dpy3 = dpy2*(d3 - 0.5d0)/(d3 + 0.5d0)
			m3 = 1.0d0
		endif

		if( cidp4 /= 0 ) then
			mu4 = mu/d4*2.0d0/(d4 + d5)
			mu5 = mu/d5*2.0d0/(d4 + d5)
			rdpz4 = rdpz5*(d4 - 0.5d0)/(d4 + 0.5d0)
			dpz4 = dpz5*(d4 - 0.5d0)/(d4 + 0.5d0)
			m4 = 1.0d0
		endif

		if( cidp5 /= 0 ) then
			mu4 = mu/d4*2.0d0/(d4 + d5)
			mu5 = mu/d5*2.0d0/(d4 + d5)
			rdpz5 = rdpz4*(d5 - 0.5d0)/(d5 + 0.5d0)
			dpz5 = dpz4*(d5 - 0.5d0)/(d5 + 0.5d0)
			m5 = 1.0d0
		endif

		if( cidp0 /= 0 .and. cidp1 /= 0 ) then
			rdpx0 = 0.0
			rdpx1 = 0.0
			dpx0 = 0.0
			dpx1 = 0.0
			m0 = 1.0d0
			m1 = 1.0d0
		endif

		if( cidp2 /= 0 .and. cidp3 /= 0 ) then
			rdpy2 = 0.0
			rdpy3 = 0.0
			dpy2 = 0.0
			dpy3 = 0.0
			m2 = 1.0d0
			m3 = 1.0d0
		endif

		if( cidp4 /= 0 .and. cidp5 /= 0 ) then
			rdpz4 = 0.0
			rdpz5 = 0.0
			dpz4 = 0.0
			dpz5 = 0.0
			m4 = 1.0d0
			m5 = 1.0d0
		endif

		l0 = mu0/(rhof)*dt/(dx*dx)
		l1 = mu1/(rhof)*dt/(dx*dx)
		l2 = mu2/(rhof)*dt/(dx*dx)
		l3 = mu3/(rhof)*dt/(dx*dx)
		l4 = mu4/(rhof)*dt/(dx*dx)
		l5 = mu5/(rhof)*dt/(dx*dx)

		rdpx = 0.5d0*(rdpx0 + rdpx1)
		rdpy = 0.5d0*(rdpy2 + rdpy3)
		rdpz = 0.5d0*(rdpz4 + rdpz5)

		rdp = 0.0
    u0 = 0.0
		if( dir==0 ) then
			rdp = rdpx 
			u0 = ux0
		else if( dir==1 ) then
			rdp = rdpy
			u0 = uy0
		else if( dir==2 ) then
			rdp = rdpz
			u0 = uz0
		end if

		Ap(i, j, k) = 1.0d0 + (l0 + l1) &
											  + (l2 + l3) &
											  + (l4 + l5)
		Aw(i, j, k) = -l0*(1.0d0 - m0)
		Ae(i, j, k) = -l1*(1.0d0 - m1)
		As(i, j, k) = -l2*(1.0d0 - m2)
		An(i, j, k) = -l3*(1.0d0 - m3)
		Ab(i, j, k) = -l4*(1.0d0 - m4)
		At(i, j, k) = -l5*(1.0d0 - m5)
		 b(i, j, k) = u0 &
									+ l0*Us*m0 &
									+ l1*Us*m1 &
									+ l2*Us*m2 &
									+ l3*Us*m3 &
									+ l4*Us*m4 &
									+ l5*Us*m5 &
									- uc0_(i, j, k)*dt &
									- rdp*dt 

		if( pidp /= 1 ) then
			Ap(i, j, k) = 1.0d0
			Aw(i, j, k) = 0.0d0
			Ae(i, j, k) = 0.0d0
			As(i, j, k) = 0.0d0
			An(i, j, k) = 0.0d0
			Ab(i, j, k) = 0.0d0
			At(i, j, k) = 0.0d0
			b (i, j, k) = Us
		endif

  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bcut_calc_ab_u_e1

subroutine bcut_remove_p( &
								ux, uy, uz, &
								p0_, &
								c0, c1, c2, c3, c4, c5, &
								cid0, cid1, cid2, cid3, cid4, cid5, &
								pid, &
								rhof, &
								dx, dt, &
								sz, g)
  implicit none
  integer                 :: i, j, k
  integer                 :: ix, jx, kx
  integer                 :: g
  integer, dimension(3)   :: sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: ux, uy, uz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: p0_
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: c0, c1, c2, c3, c4, c5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: cid0, cid1, cid2, cid3, cid4, cid5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: pid
	integer									:: cidp
	integer									:: cidp0, cidp1, cidp2, cidp3, cidp4, cidp5
	integer									:: pidp, pidw, pide, pids, pidn, pidb, pidt
	real										:: rhof
  real                    :: dx, dt
	real										:: d0, d1, d2, d3, d4, d5
	real										:: r0, r1, r2, r3, r4, r5
	real										:: m0, m1, m2, m3, m4, m5
	real										:: l0, l1, l2, l3, l4, l5
	real										:: dpx0, dpx1, dpy2, dpy3, dpz4, dpz5
	real										:: rdpx0, rdpx1, rdpy2, rdpy3, rdpz4, rdpz5
	real										:: rdpx, rdpy, rdpz, rdp
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) &
!$omp					 private(cidp) &
!$omp					 private(cidp0, cidp1, cidp2, cidp3, cidp4, cidp5) &
!$omp					 private(pidp, pidw, pide, pids, pidn, pidb, pidt) &
!$omp					 private(d0, d1, d2, d3, d4, d5) &
!$omp					 private(r0, r1, r2, r3, r4, r5) &
!$omp					 private(m0, m1, m2, m3, m4, m5) &
!$omp					 private(l0, l1, l2, l3, l4, l5) &
!$omp					 private(dpx0, dpx1, dpy2, dpy3, dpz4, dpz5) &
!$omp					 private(rdpx0, rdpx1, rdpy2, rdpy3, rdpz4, rdpz5) &
!$omp					 private(rdpx, rdpy, rdpz, rdp)
!$omp do schedule(static, 1)
#else
#endif
  do k=1, kx
  do j=1, jx
!ocl nouxsimd
  do i=1, ix
		r0 = 1.0d0/rhof
		r1 = 1.0d0/rhof
		r2 = 1.0d0/rhof
		r3 = 1.0d0/rhof
		r4 = 1.0d0/rhof
		r5 = 1.0d0/rhof

		dpx0 = (p0_(i, j, k) - p0_(i-1, j, k))/dx
		dpx1 = (p0_(i+1, j, k) - p0_(i, j, k))/dx
		dpy2 = (p0_(i, j, k) - p0_(i, j-1, k))/dx
		dpy3 = (p0_(i, j+1, k) - p0_(i, j, k))/dx
		dpz4 = (p0_(i, j, k) - p0_(i, j, k-1))/dx
		dpz5 = (p0_(i, j, k+1) - p0_(i, j, k))/dx

		rdpx0 = r0*dpx0
		rdpx1 = r1*dpx1
		rdpy2 = r2*dpy2
		rdpy3 = r3*dpy3
		rdpz4 = r4*dpz4
		rdpz5 = r5*dpz5

		m0 = 0.0d0
		m1 = 0.0d0
		m2 = 0.0d0
		m3 = 0.0d0
		m4 = 0.0d0
		m5 = 0.0d0

		d0 = c0(i, j, k)
		d1 = c1(i, j, k)
		d2 = c2(i, j, k)
		d3 = c3(i, j, k)
		d4 = c4(i, j, k)
		d5 = c5(i, j, k)

		cidp0 = cid0(i, j, k)
		cidp1 = cid1(i, j, k)
		cidp2 = cid2(i, j, k)
		cidp3 = cid3(i, j, k)
		cidp4 = cid4(i, j, k)
		cidp5 = cid5(i, j, k)

		pidp = pid(i, j, k)

		if( cidp0 /= 0 ) then
			rdpx0 = rdpx1*(d0 - 0.5d0)/(d0 + 0.5d0)
			dpx0 = dpx1*(d0 - 0.5d0)/(d0 + 0.5d0)
			m0 = 1.0d0
		endif

		if( cidp1 /= 0 ) then
			rdpx1 = rdpx0*(d1 - 0.5d0)/(d1 + 0.5d0)
			dpx1 = dpx0*(d1 - 0.5d0)/(d1 + 0.5d0)
			m1 = 1.0d0
		endif

		if( cidp2 /= 0 ) then
			rdpy2 = rdpy3*(d2 - 0.5d0)/(d2 + 0.5d0)
			dpy2 = dpy3*(d2 - 0.5d0)/(d2 + 0.5d0)
			m2 = 1.0d0
		endif

		if( cidp3 /= 0 ) then
			rdpy3 = rdpy2*(d3 - 0.5d0)/(d3 + 0.5d0)
			dpy3 = dpy2*(d3 - 0.5d0)/(d3 + 0.5d0)
			m3 = 1.0d0
		endif

		if( cidp4 /= 0 ) then
			rdpz4 = rdpz5*(d4 - 0.5d0)/(d4 + 0.5d0)
			dpz4 = dpz5*(d4 - 0.5d0)/(d4 + 0.5d0)
			m4 = 1.0d0
		endif

		if( cidp5 /= 0 ) then
			rdpz5 = rdpz4*(d5 - 0.5d0)/(d5 + 0.5d0)
			dpz5 = dpz4*(d5 - 0.5d0)/(d5 + 0.5d0)
			m5 = 1.0d0
		endif

		if( cidp0 /= 0 .and. cidp1 /= 0 ) then
			rdpx0 = 0.0
			rdpx1 = 0.0
			dpx0 = 0.0
			dpx1 = 0.0
			m0 = 1.0d0
			m1 = 1.0d0
		endif

		if( cidp2 /= 0 .and. cidp3 /= 0 ) then
			rdpy2 = 0.0
			rdpy3 = 0.0
			dpy2 = 0.0
			dpy3 = 0.0
			m2 = 1.0d0
			m3 = 1.0d0
		endif

		if( cidp4 /= 0 .and. cidp5 /= 0 ) then
			rdpz4 = 0.0
			rdpz5 = 0.0
			dpz4 = 0.0
			dpz5 = 0.0
			m4 = 1.0d0
			m5 = 1.0d0
		endif

		rdpx = 0.5d0*(rdpx0 + rdpx1)
		rdpy = 0.5d0*(rdpy2 + rdpy3)
		rdpz = 0.5d0*(rdpz4 + rdpz5)

		ux(i, j, k) = ux(i, j, k) + rdpx*dt
		uy(i, j, k) = uy(i, j, k) + rdpy*dt
		uz(i, j, k) = uz(i, j, k) + rdpz*dt

		if( pidp /= 1 ) then
			ux(i, j, k) = 0.0
			uy(i, j, k) = 0.0
			uz(i, j, k) = 0.0
		endif

  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bcut_remove_p

subroutine bcut_calc_ab_p( &
								Ap, Aw, Ae, As, An, Ab, At, b, &
								vw, ve, vs, vn, vb, vt, &
								p0_, &
								ux, uy, uz, &
								c0, c1, c2, c3, c4, c5, &
								cid0, cid1, cid2, cid3, cid4, cid5, &
								pid, &
								rhof, &
								dx, dt, &
								sz, g)
  implicit none
  integer                 :: i, j, k
  integer                 :: ix, jx, kx
  integer                 :: g
  integer, dimension(3)   :: sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: Ap, Aw, Ae, As, An, Ab, At, b
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: vw, ve, vs, vn, vb, vt
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: p0_
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: ux, uy, uz
	real																										:: vx0, vx1, vy2, vy3, vz4, vz5
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: c0, c1, c2, c3, c4, c5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: cid0, cid1, cid2, cid3, cid4, cid5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: pid
	integer									:: cidp
	integer									:: cidp0, cidp1, cidp2, cidp3, cidp4, cidp5
	integer									:: pidp, pidw, pide, pids, pidn, pidb, pidt
	real										:: rhof
  real                    :: dx, dt
	real										:: d0, d1, d2, d3, d4, d5
	real										:: r0, r1, r2, r3, r4, r5
	real										:: m0, m1, m2, m3, m4, m5
	real										:: l0, l1, l2, l3, l4, l5
	real										:: divv
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) &
!$omp					 private(vx0, vx1, vy2, vy3, vz4, vz5) &
!$omp					 private(cidp) &
!$omp					 private(cidp0, cidp1, cidp2, cidp3, cidp4, cidp5) &
!$omp					 private(pidp, pidw, pide, pids, pidn, pidb, pidt) &
!$omp					 private(d0, d1, d2, d3, d4, d5) &
!$omp					 private(r0, r1, r2, r3, r4, r5) &
!$omp					 private(m0, m1, m2, m3, m4, m5) &
!$omp					 private(l0, l1, l2, l3, l4, l5) &
!$omp					 private(divv)
!$omp do schedule(static, 1)
#else
#endif
  do k=1, kx
  do j=1, jx
!ocl nouxsimd
  do i=1, ix
		r0 = 1.0d0/rhof
		r1 = 1.0d0/rhof
		r2 = 1.0d0/rhof
		r3 = 1.0d0/rhof
		r4 = 1.0d0/rhof
		r5 = 1.0d0/rhof

		vx0 = 0.5d0*(ux(i, j, k) + ux(i-1, j, k))
		vx1 = 0.5d0*(ux(i, j, k) + ux(i+1, j, k))
		vy2 = 0.5d0*(uy(i, j, k) + uy(i, j-1, k))
		vy3 = 0.5d0*(uy(i, j, k) + uy(i, j+1, k))
		vz4 = 0.5d0*(uz(i, j, k) + uz(i, j, k-1))
		vz5 = 0.5d0*(uz(i, j, k) + uz(i, j, k+1))

		m0 = 0.0d0
		m1 = 0.0d0
		m2 = 0.0d0
		m3 = 0.0d0
		m4 = 0.0d0
		m5 = 0.0d0

		d0 = c0(i, j, k)
		d1 = c1(i, j, k)
		d2 = c2(i, j, k)
		d3 = c3(i, j, k)
		d4 = c4(i, j, k)
		d5 = c5(i, j, k)

		cidp0 = cid0(i, j, k)
		cidp1 = cid1(i, j, k)
		cidp2 = cid2(i, j, k)
		cidp3 = cid3(i, j, k)
		cidp4 = cid4(i, j, k)
		cidp5 = cid5(i, j, k)

		pidp = pid(i, j, k)

		if( cidp0 /= 0 ) then
			r0 = 0.0d0
			vx0 = 0.0d0
			vx0 = vx1*(d0 - 0.5d0)/(d0 + 0.5d0)
			r1 = (1.0d0/rhof)/(d0 + 0.5d0)
			m0 = 1.0d0
		endif

		if( cidp1 /= 0 ) then
			r1 = 0.0d0
			vx1 = 0.0d0
			vx1 = vx0*(d1 - 0.5d0)/(d1 + 0.5d0)
			r0 = (1.0d0/rhof)/(d1 + 0.5d0)
			m1 = 1.0d0
		endif

		if( cidp2 /= 0 ) then
			r2 = 0.0d0
			vy2 = 0.0d0
			vy2 = vy3*(d2 - 0.5d0)/(d2 + 0.5d0)
			r3 = (1.0d0/rhof)/(d2 + 0.5d0)
			m2 = 1.0d0
		endif

		if( cidp3 /= 0 ) then
			r3 = 0.0d0
			vy3 = 0.0d0
			vy3 = vy2*(d3 - 0.5d0)/(d3 + 0.5d0)
			r2 = (1.0d0/rhof)/(d3 + 0.5d0)
			m3 = 1.0d0
		endif

		if( cidp4 /= 0 ) then
			r4 = 0.0d0
			vz4 = 0.0d0
			vz4 = vz5*(d4 - 0.5d0)/(d4 + 0.5d0)
			r5 = (1.0d0/rhof)/(d4 + 0.5d0)
			m4 = 1.0d0
		endif

		if( cidp5 /= 0 ) then
			r5 = 0.0d0
			vz5 = 0.0d0
			vz5 = vz4*(d5 - 0.5d0)/(d5 + 0.5d0)
			r4 = (1.0d0/rhof)/(d5 + 0.5d0)
			m5 = 1.0d0
		endif

		if( cidp0 /= 0 .and. cidp1 /= 0 ) then
			r0 = 0.0d0
			r1 = 0.0d0
			vx0 = 0.5d0*(2.0d0 - 1.0d0/d0)*ux(i, j, k)
			vx1 = 0.5d0*(2.0d0 - 1.0d0/d1)*ux(i, j, k)
			vx0 = 0.0d0
			vx1 = 0.0d0
			m0 = 1.0d0
			m1 = 1.0d0
		endif

		if( cidp2 /= 0 .and. cidp3 /= 0 ) then
			r2 = 0.0d0
			r3 = 0.0d0
			vy2 = 0.5d0*(2.0d0 - 1.0d0/d2)*uy(i, j, k)
			vy3 = 0.5d0*(2.0d0 - 1.0d0/d3)*uy(i, j, k)
			vy2 = 0.0d0
			vy3 = 0.0d0
			m2 = 1.0d0
			m3 = 1.0d0
		endif

		if( cidp4 /= 0 .and. cidp5 /= 0 ) then
			r4 = 0.0d0
			r5 = 0.0d0
			vz4 = 0.5d0*(2.0d0 - 1.0d0/d4)*uz(i, j, k)
			vz5 = 0.5d0*(2.0d0 - 1.0d0/d5)*uz(i, j, k)
			vz4 = 0.0d0
			vz5 = 0.0d0
			m4 = 1.0d0
			m5 = 1.0d0
		endif

		l0 = r0/(dx*dx)
		l1 = r1/(dx*dx)
		l2 = r2/(dx*dx)
		l3 = r3/(dx*dx)
		l4 = r4/(dx*dx)
		l5 = r5/(dx*dx)

		divv = (vx1 - vx0 + vy3 - vy2 + vz5 - vz4)/dx

		Ap(i, j, k) = - (l0 + l1) &
									- (l2 + l3) &
									- (l4 + l5)
		Aw(i, j, k) = l0
		Ae(i, j, k) = l1
		As(i, j, k) = l2
		An(i, j, k) = l3
		Ab(i, j, k) = l4
		At(i, j, k) = l5
		 b(i, j, k) = divv/dt
		vw(i, j, k) = vx0
		ve(i, j, k) = vx1
		vs(i, j, k) = vy2
		vn(i, j, k) = vy3
		vb(i, j, k) = vz4
		vt(i, j, k) = vz5

		if( pidp /= 1 ) then
			Ap(i, j, k) =-6.0d0
			Aw(i, j, k) = 1.0d0
			Ae(i, j, k) = 1.0d0
			As(i, j, k) = 1.0d0
			An(i, j, k) = 1.0d0
			Ab(i, j, k) = 1.0d0
			At(i, j, k) = 1.0d0
			b (i, j, k) = 0.0d0
			vw(i, j, k) = 0.0d0
			ve(i, j, k) = 0.0d0
			vs(i, j, k) = 0.0d0
			vn(i, j, k) = 0.0d0
			vb(i, j, k) = 0.0d0
			vt(i, j, k) = 0.0d0
		endif

  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bcut_calc_ab_p

subroutine bcut_corr_u( &
								ux0_, uy0_, uz0_, &
								vw, ve, vs, vn, vb, vt, &
								lapp, &
								p0_, &
								c0, c1, c2, c3, c4, c5, &
								cid0, cid1, cid2, cid3, cid4, cid5, &
								pid, &
								rhof, &
								dx, dt, &
								sz, g)
  implicit none
  integer                 :: i, j, k
  integer                 :: ix, jx, kx
  integer                 :: g
  integer, dimension(3)   :: sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: ux0_, uy0_, uz0_
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: vw, ve, vs, vn, vb, vt
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: lapp
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: p0_
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: c0, c1, c2, c3, c4, c5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: cid0, cid1, cid2, cid3, cid4, cid5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: pid
	integer									:: cidp
	integer									:: cidp0, cidp1, cidp2, cidp3, cidp4, cidp5
	integer									:: pidp, pidw, pide, pids, pidn, pidb, pidt
	real										:: rhof
  real                    :: dx, dt
	real										:: d0, d1, d2, d3, d4, d5
	real										:: r0, r1, r2, r3, r4, r5
	real										:: m0, m1, m2, m3, m4, m5
	real										:: l0, l1, l2, l3, l4, l5
	real										:: dpx0, dpx1, dpy2, dpy3, dpz4, dpz5
	real										:: rdpx0, rdpx1, rdpy2, rdpy3, rdpz4, rdpz5
	real										:: rdpx, rdpy, rdpz
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) &
!$omp					 private(cidp) &
!$omp					 private(cidp0, cidp1, cidp2, cidp3, cidp4, cidp5) &
!$omp					 private(pidp, pidw, pide, pids, pidn, pidb, pidt) &
!$omp					 private(d0, d1, d2, d3, d4, d5) &
!$omp					 private(r0, r1, r2, r3, r4, r5) &
!$omp					 private(m0, m1, m2, m3, m4, m5) &
!$omp					 private(l0, l1, l2, l3, l4, l5) &
!$omp					 private(dpx0, dpx1, dpy2, dpy3, dpz4, dpz5) &
!$omp					 private(rdpx0, rdpx1, rdpy2, rdpy3, rdpz4, rdpz5) &
!$omp					 private(rdpx, rdpy, rdpz)
!$omp do schedule(static, 1)
#else
#endif
  do k=1, kx
  do j=1, jx
!ocl nouxsimd
  do i=1, ix
		r0 = 1.0d0/rhof
		r1 = 1.0d0/rhof
		r2 = 1.0d0/rhof
		r3 = 1.0d0/rhof
		r4 = 1.0d0/rhof
		r5 = 1.0d0/rhof

		dpx0 = (p0_(i, j, k) - p0_(i-1, j, k))/dx
		dpx1 = (p0_(i+1, j, k) - p0_(i, j, k))/dx
		dpy2 = (p0_(i, j, k) - p0_(i, j-1, k))/dx
		dpy3 = (p0_(i, j+1, k) - p0_(i, j, k))/dx
		dpz4 = (p0_(i, j, k) - p0_(i, j, k-1))/dx
		dpz5 = (p0_(i, j, k+1) - p0_(i, j, k))/dx

		rdpx0 = r0*dpx0
		rdpx1 = r1*dpx1
		rdpy2 = r2*dpy2
		rdpy3 = r3*dpy3
		rdpz4 = r4*dpz4
		rdpz5 = r5*dpz5

		m0 = 0.0d0
		m1 = 0.0d0
		m2 = 0.0d0
		m3 = 0.0d0
		m4 = 0.0d0
		m5 = 0.0d0

		d0 = c0(i, j, k)
		d1 = c1(i, j, k)
		d2 = c2(i, j, k)
		d3 = c3(i, j, k)
		d4 = c4(i, j, k)
		d5 = c5(i, j, k)

		cidp0 = cid0(i, j, k)
		cidp1 = cid1(i, j, k)
		cidp2 = cid2(i, j, k)
		cidp3 = cid3(i, j, k)
		cidp4 = cid4(i, j, k)
		cidp5 = cid5(i, j, k)

		pidp = pid(i, j, k)

		if( cidp0 /= 0 ) then
			rdpx0 = rdpx1
			rdpx0 = 0.0d0
			rdpx0 = rdpx1*(d0 - 0.5d0)/(d0 + 0.5d0)
			dpx0 = 0.0d0
			dpx0 = dpx1*(d0 - 0.5d0)/(d0 + 0.5d0)
			m0 = 1.0d0
		endif

		if( cidp1 /= 0 ) then
			rdpx1 = rdpx0
			rdpx1 = 0.0d0
			rdpx1 = rdpx0*(d1 - 0.5d0)/(d1 + 0.5d0)
			dpx1 = 0.0d0
			dpx1 = dpx0*(d1 - 0.5d0)/(d1 + 0.5d0)
			m1 = 1.0d0
		endif

		if( cidp2 /= 0 ) then
			rdpy2 = rdpy3
			rdpy2 = 0.0d0
			rdpy2 = rdpy3*(d2 - 0.5d0)/(d2 + 0.5d0)
			dpy2 = 0.0d0
			dpy2 = dpy3*(d2 - 0.5d0)/(d2 + 0.5d0)
			m2 = 1.0d0
		endif

		if( cidp3 /= 0 ) then
			rdpy3 = rdpy2
			rdpy3 = 0.0d0
			rdpy3 = rdpy2*(d3 - 0.5d0)/(d3 + 0.5d0)
			dpy3 = 0.0d0
			dpy3 = dpy2*(d3 - 0.5d0)/(d3 + 0.5d0)
			m3 = 1.0d0
		endif

		if( cidp4 /= 0 ) then
			rdpz4 = rdpz5
			rdpz4 = 0.0d0
			rdpz4 = rdpz5*(d4 - 0.5d0)/(d4 + 0.5d0)
			dpz4 = 0.0d0
			dpz4 = dpz5*(d4 - 0.5d0)/(d4 + 0.5d0)
			m4 = 1.0d0
		endif

		if( cidp5 /= 0 ) then
			rdpz5 = rdpz4
			rdpz5 = 0.0d0
			rdpz5 = rdpz4*(d5 - 0.5d0)/(d5 + 0.5d0)
			dpz5 = 0.0d0
			dpz5 = dpz4*(d5 - 0.5d0)/(d5 + 0.5d0)
			m5 = 1.0d0
		endif

		if( cidp0 /= 0 .and. cidp1 /= 0 ) then
			rdpx0 = 0.0
			rdpx1 = 0.0
			dpx0 = 0.0
			dpx1 = 0.0
			m0 = 1.0d0
			m1 = 1.0d0
		endif

		if( cidp2 /= 0 .and. cidp3 /= 0 ) then
			rdpy2 = 0.0
			rdpy3 = 0.0
			dpy2 = 0.0
			dpy3 = 0.0
			m2 = 1.0d0
			m3 = 1.0d0
		endif

		if( cidp4 /= 0 .and. cidp5 /= 0 ) then
			rdpz4 = 0.0
			rdpz5 = 0.0
			dpz4 = 0.0
			dpz5 = 0.0
			m4 = 1.0d0
			m5 = 1.0d0
		endif

		rdpx = 0.5d0*(rdpx0 + rdpx1)
		rdpy = 0.5d0*(rdpy2 + rdpy3)
		rdpz = 0.5d0*(rdpz4 + rdpz5)

		ux0_(i, j, k) = ux0_(i, j, k) - rdpx*dt
		uy0_(i, j, k) = uy0_(i, j, k) - rdpy*dt
		uz0_(i, j, k) = uz0_(i, j, k) - rdpz*dt

		vw(i, j, k) = vw(i, j, k) - rdpx0*dt
		ve(i, j, k) = ve(i, j, k) - rdpx1*dt
		vs(i, j, k) = vs(i, j, k) - rdpy2*dt
		vn(i, j, k) = vn(i, j, k) - rdpy3*dt
		vb(i, j, k) = vb(i, j, k) - rdpz4*dt
		vt(i, j, k) = vt(i, j, k) - rdpz5*dt

		lapp(i, j, k) = (rdpx1 - rdpx0 + rdpy3 - rdpy2 + rdpz5 - rdpz4)/dx*rhof
		lapp(i, j, k) = (dpx1 - dpx0 + dpy3 - dpy2 + dpz5 - dpz4)/dx
		lapp(i, j, k) = ( p0_(i+1, j, k) + p0_(i-1, j, k) &
										+ p0_(i, j+1, k) + p0_(i, j-1, k) &
										+ p0_(i, j, k+1) + p0_(i, j, k-1) &
										- 6.0*p0_(i, j, k) )/(dx*dx)

		if( pidp /= 1 ) then
			ux0_(i, j, k) = 0.0d0
			uy0_(i, j, k) = 0.0d0
			uz0_(i, j, k) = 0.0d0

			vw(i, j, k) = 0.0d0
			ve(i, j, k) = 0.0d0
			vs(i, j, k) = 0.0d0
			vn(i, j, k) = 0.0d0
			vb(i, j, k) = 0.0d0
			vt(i, j, k) = 0.0d0
			lapp(i, j, k) = 0.0d0
		endif

		if( cidp0 /= 0 .or. &
				cidp1 /= 0 .or. &
				cidp2 /= 0 .or. &
				cidp3 /= 0 .or. &
				cidp4 /= 0 .or. &
				cidp5 /= 0 ) then
			lapp(i, j, k) = 0.0d0
		endif

  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bcut_corr_u

subroutine bcut_calc_ab_t_1st( &
								Ap, Aw, Ae, As, An, Ab, At, b, &
								t0_, &
								tc0_, &
								c0, c1, c2, c3, c4, c5, &
								cid0, cid1, cid2, cid3, cid4, cid5, &
								pid, &
								rhof, rhos, &
								cpf, cps, &
								kf, ks, &
								dx, dt, &
								Ts, &
								sz, g)
  implicit none
  integer                 :: i, j, k
  integer                 :: ix, jx, kx
  integer                 :: g
  integer, dimension(3)   :: sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: Ap, Aw, Ae, As, An, Ab, At, b
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: t0_
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: tc0_
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: c0, c1, c2, c3, c4, c5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: cid0, cid1, cid2, cid3, cid4, cid5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: pid
	integer									:: cidp
	integer									:: cidp0, cidp1, cidp2, cidp3, cidp4, cidp5
	integer									:: pidp, pidw, pide, pids, pidn, pidb, pidt
	real										:: rhof, rhos
	real										:: cpf, cps
	real										:: kf, ks
  real                    :: dx, dt
	real										:: Ts
	real										:: d0, d1, d2, d3, d4, d5
	real										:: k0, k1, k2, k3, k4, k5
	real										:: m0, m1, m2, m3, m4, m5
	real										:: l0, l1, l2, l3, l4, l5
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) &
!$omp					 private(cidp) &
!$omp					 private(cidp0, cidp1, cidp2, cidp3, cidp4, cidp5) &
!$omp					 private(pidp, pidw, pide, pids, pidn, pidb, pidt) &
!$omp					 private(d0, d1, d2, d3, d4, d5) &
!$omp					 private(k0, k1, k2, k3, k4, k5) &
!$omp					 private(m0, m1, m2, m3, m4, m5) &
!$omp					 private(l0, l1, l2, l3, l4, l5) 
!$omp do schedule(static, 1)
#else
#endif
  do k=1, kx
  do j=1, jx
!ocl nouxsimd
  do i=1, ix
		k0 = kf
		k1 = kf
		k2 = kf
		k3 = kf
		k4 = kf
		k5 = kf

		m0 = 0.0
		m1 = 0.0
		m2 = 0.0
		m3 = 0.0
		m4 = 0.0
		m5 = 0.0

		d0 = c0(i, j, k)
		d1 = c1(i, j, k)
		d2 = c2(i, j, k)
		d3 = c3(i, j, k)
		d4 = c4(i, j, k)
		d5 = c5(i, j, k)

		cidp0 = cid0(i, j, k)
		cidp1 = cid1(i, j, k)
		cidp2 = cid2(i, j, k)
		cidp3 = cid3(i, j, k)
		cidp4 = cid4(i, j, k)
		cidp5 = cid5(i, j, k)

		pidp = pid(i, j, k)

		if( cidp0 /= 0 ) then
			k0 = kf/d0*2.0/(d0 + d1)
			k1 = kf/d1*2.0/(d0 + d1)
			m0 = 1.0
		endif

		if( cidp1 /= 0 ) then
			k0 = kf/d0*2.0/(d0 + d1)
			k1 = kf/d1*2.0/(d0 + d1)
			m1 = 1.0
		endif

		if( cidp2 /= 0 ) then
			k2 = kf/d2*2.0/(d2 + d3)
			k3 = kf/d3*2.0/(d2 + d3)
			m2 = 1.0
		endif

		if( cidp3 /= 0 ) then
			k2 = kf/d2*2.0/(d2 + d3)
			k3 = kf/d3*2.0/(d2 + d3)
			m3 = 1.0
		endif

		if( cidp4 /= 0 ) then
			k4 = kf/d4*2.0/(d4 + d5)
			k5 = kf/d5*2.0/(d4 + d5)
			m4 = 1.0
		endif

		if( cidp5 /= 0 ) then
			k4 = kf/d4*2.0/(d4 + d5)
			k5 = kf/d5*2.0/(d4 + d5)
			m5 = 1.0
		endif

		l0 = k0/(rhof*cpf)*dt/(dx*dx)
		l1 = k1/(rhof*cpf)*dt/(dx*dx)
		l2 = k2/(rhof*cpf)*dt/(dx*dx)
		l3 = k3/(rhof*cpf)*dt/(dx*dx)
		l4 = k4/(rhof*cpf)*dt/(dx*dx)
		l5 = k5/(rhof*cpf)*dt/(dx*dx)

		Ap(i, j, k) = 1.0d0 + (l0 + l1) &
												+ (l2 + l3) &
												+ (l4 + l5)
		Aw(i, j, k) = -l0*(1.0d0 - m0)
		Ae(i, j, k) = -l1*(1.0d0 - m1)
		As(i, j, k) = -l2*(1.0d0 - m2)
		An(i, j, k) = -l3*(1.0d0 - m3)
		Ab(i, j, k) = -l4*(1.0d0 - m4)
		At(i, j, k) = -l5*(1.0d0 - m5)
		 b(i, j, k) = t0_(i, j, k) &
									+ l0*Ts*m0 &
									+ l1*Ts*m1 &
									+ l2*Ts*m2 &
									+ l3*Ts*m3 &
									+ l4*Ts*m4 &
									+ l5*Ts*m5 &
									- tc0_(i, j, k)*dt

		if( pidp /= 1 ) then
			Ap(i, j, k) = 1.0
			Aw(i, j, k) = 0.0
			Ae(i, j, k) = 0.0
			As(i, j, k) = 0.0
			An(i, j, k) = 0.0
			Ab(i, j, k) = 0.0
			At(i, j, k) = 0.0
			b (i, j, k) = Ts
		endif

  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bcut_calc_ab_t_1st

subroutine bcut_calc_ab_t( &
								Ap, Aw, Ae, As, An, Ab, At, b, &
								t0_, &
								tc0_, tcp_, &
								td0_, &
								c0, c1, c2, c3, c4, c5, &
								cid0, cid1, cid2, cid3, cid4, cid5, &
								pid, &
								rhof, rhos, &
								cpf, cps, &
								kf, ks, &
								dx, dt, &
								Tc, &
								sz, g)
  implicit none
  integer                 :: i, j, k
  integer                 :: ix, jx, kx
  integer                 :: g
  integer, dimension(3)   :: sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: Ap, Aw, Ae, As, An, Ab, At, b
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: t0_
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: tc0_, tcp_
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: td0_
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: c0, c1, c2, c3, c4, c5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: cid0, cid1, cid2, cid3, cid4, cid5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: pid
	integer									:: cidp
	integer									:: cidp0, cidp1, cidp2, cidp3, cidp4, cidp5
	integer									:: pidp, pidw, pide, pids, pidn, pidb, pidt
	real										:: rhof, rhos
	real										:: cpf, cps
	real										:: kf, ks
  real                    :: dx, dt
	real										:: Tc
	real										:: d0, d1, d2, d3, d4, d5
	real										:: k0, k1, k2, k3, k4, k5
	real										:: m0, m1, m2, m3, m4, m5
	real										:: l0, l1, l2, l3, l4, l5
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) &
!$omp					 private(cidp) &
!$omp					 private(cidp0, cidp1, cidp2, cidp3, cidp4, cidp5) &
!$omp					 private(pidp, pidw, pide, pids, pidn, pidb, pidt) &
!$omp					 private(d0, d1, d2, d3, d4, d5) &
!$omp					 private(k0, k1, k2, k3, k4, k5) &
!$omp					 private(m0, m1, m2, m3, m4, m5) &
!$omp					 private(l0, l1, l2, l3, l4, l5) 
!$omp do schedule(static, 1)
#else
#endif
  do k=1, kx
  do j=1, jx
!ocl nouxsimd
  do i=1, ix
		k0 = kf
		k1 = kf
		k2 = kf
		k3 = kf
		k4 = kf
		k5 = kf

		m0 = 0.0
		m1 = 0.0
		m2 = 0.0
		m3 = 0.0
		m4 = 0.0
		m5 = 0.0

		d0 = c0(i, j, k)
		d1 = c1(i, j, k)
		d2 = c2(i, j, k)
		d3 = c3(i, j, k)
		d4 = c4(i, j, k)
		d5 = c5(i, j, k)

		cidp0 = cid0(i, j, k)
		cidp1 = cid1(i, j, k)
		cidp2 = cid2(i, j, k)
		cidp3 = cid3(i, j, k)
		cidp4 = cid4(i, j, k)
		cidp5 = cid5(i, j, k)

		pidp = pid(i, j, k)

		if( cidp0 /= 0 ) then
			k0 = kf/d0*2.0/(d0 + d1)
			k1 = kf/d1*2.0/(d0 + d1)
			m0 = 1.0
		endif

		if( cidp1 /= 0 ) then
			k0 = kf/d0*2.0/(d0 + d1)
			k1 = kf/d1*2.0/(d0 + d1)
			m1 = 1.0
		endif

		if( cidp2 /= 0 ) then
			k2 = kf/d2*2.0/(d2 + d3)
			k3 = kf/d3*2.0/(d2 + d3)
			m2 = 1.0
		endif

		if( cidp3 /= 0 ) then
			k2 = kf/d2*2.0/(d2 + d3)
			k3 = kf/d3*2.0/(d2 + d3)
			m3 = 1.0
		endif

		if( cidp4 /= 0 ) then
			k4 = kf/d4*2.0/(d4 + d5)
			k5 = kf/d5*2.0/(d4 + d5)
			m4 = 1.0
		endif

		if( cidp5 /= 0 ) then
			k4 = kf/d4*2.0/(d4 + d5)
			k5 = kf/d5*2.0/(d4 + d5)
			m5 = 1.0
		endif

		l0 = k0/(rhof*cpf)*dt/(dx*dx)
		l1 = k1/(rhof*cpf)*dt/(dx*dx)
		l2 = k2/(rhof*cpf)*dt/(dx*dx)
		l3 = k3/(rhof*cpf)*dt/(dx*dx)
		l4 = k4/(rhof*cpf)*dt/(dx*dx)
		l5 = k5/(rhof*cpf)*dt/(dx*dx)

!		Ap(i, j, k) = 1.0d0 + (l0 + l1) &
!												+ (l2 + l3) &
!												+ (l4 + l5)
!		Aw(i, j, k) = -l0*(1.0d0 - m0)
!		Ae(i, j, k) = -l1*(1.0d0 - m1)
!		As(i, j, k) = -l2*(1.0d0 - m2)
!		An(i, j, k) = -l3*(1.0d0 - m3)
!		Ab(i, j, k) = -l4*(1.0d0 - m4)
!		At(i, j, k) = -l5*(1.0d0 - m5)
!		 b(i, j, k) = t0_(i, j, k) &
!									+ l0*Tc*m0 &
!									+ l1*Tc*m1 &
!									+ l2*Tc*m2 &
!									+ l3*Tc*m3 &
!									+ l4*Tc*m4 &
!									+ l5*Tc*m5 &
!									- (1.5d0*tc0_(i, j, k) - 0.5d0*tcp_(i, j, k))*dt 

		Ap(i, j, k) = 1.0d0 + 0.5d0*(l0 + l1) &
												+ 0.5d0*(l2 + l3) &
												+ 0.5d0*(l4 + l5)
		Aw(i, j, k) = -0.5d0*l0*(1.0d0 - m0)
		Ae(i, j, k) = -0.5d0*l1*(1.0d0 - m1)
		As(i, j, k) = -0.5d0*l2*(1.0d0 - m2)
		An(i, j, k) = -0.5d0*l3*(1.0d0 - m3)
		Ab(i, j, k) = -0.5d0*l4*(1.0d0 - m4)
		At(i, j, k) = -0.5d0*l5*(1.0d0 - m5)
		 b(i, j, k) = t0_(i, j, k) &
									+ 0.5d0*l0*Tc*m0 &
									+ 0.5d0*l1*Tc*m1 &
									+ 0.5d0*l2*Tc*m2 &
									+ 0.5d0*l3*Tc*m3 &
									+ 0.5d0*l4*Tc*m4 &
									+ 0.5d0*l5*Tc*m5 &
									- (1.5d0*tc0_(i, j, k) - 0.5d0*tcp_(i, j, k))*dt &
									+ 0.5d0*td0_(i, j, k)*dt

		if( pidp /= 1 ) then
			Ap(i, j, k) = 1.0
			Aw(i, j, k) = 0.0
			Ae(i, j, k) = 0.0
			As(i, j, k) = 0.0
			An(i, j, k) = 0.0
			Ab(i, j, k) = 0.0
			At(i, j, k) = 0.0
			b (i, j, k) = Tc
		endif

  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bcut_calc_ab_t

subroutine bcut_calc_f_p( &
								fspx, &
								fspy, &
								fspz, &
								fsp, &
								p, &
								c0, c1, c2, c3, c4, c5, &
								cid0, cid1, cid2, cid3, cid4, cid5, &
								pid, &
								dx, dt, &
								sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: fspx, fspy, fspz
	real, dimension(1:3)		:: fsp
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: p
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: c0, c1, c2, c3, c4, c5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: cid0, cid1, cid2, cid3, cid4, cid5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: pid
	integer									:: cidp
	integer									:: cidp0, cidp1, cidp2, cidp3, cidp4, cidp5
	integer									:: pidp, pidw, pide, pids, pidn, pidb, pidt
  real                    :: dx, dt
	real										:: fsx, fsy, fsz
	real										:: d0, d1, d2, d3, d4, d5
	real										:: m0, m1, m2, m3, m4, m5
	real										:: p0, p1, p2, p3, p4, p5
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	fsx = 0.0
	fsy = 0.0
	fsz = 0.0
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) &
!$omp					 private(cidp) &
!$omp					 private(cidp0, cidp1, cidp2, cidp3, cidp4, cidp5) &
!$omp					 private(pidp, pidw, pide, pids, pidn, pidb, pidt) &
!$omp					 private(d0, d1, d2, d3, d4, d5) &
!$omp					 private(m0, m1, m2, m3, m4, m5) &
!$omp					 private(p0, p1, p2, p3, p4, p5) 
!$omp do schedule(static, 1), &
!$omp		 reduction(+:fsx, fsy, fsz)
#else
#endif
	do k=1, kx
	do j=1, jx
!ocl nouxsimd
	do i=1, ix
		p0 = 0.5*(p(i, j, k) + p(i-1, j, k))
		p1 = 0.5*(p(i, j, k) + p(i+1, j, k))
		p2 = 0.5*(p(i, j, k) + p(i, j-1, k))
		p3 = 0.5*(p(i, j, k) + p(i, j+1, k))
		p4 = 0.5*(p(i, j, k) + p(i, j, k-1))
		p5 = 0.5*(p(i, j, k) + p(i, j, k+1))

		m0 = 0.0d0
		m1 = 0.0d0
		m2 = 0.0d0
		m3 = 0.0d0
		m4 = 0.0d0
		m5 = 0.0d0

		d0 = c0(i, j, k)
		d1 = c1(i, j, k)
		d2 = c2(i, j, k)
		d3 = c3(i, j, k)
		d4 = c4(i, j, k)
		d5 = c5(i, j, k)

		cidp0 = cid0(i, j, k)
		cidp1 = cid1(i, j, k)
		cidp2 = cid2(i, j, k)
		cidp3 = cid3(i, j, k)
		cidp4 = cid4(i, j, k)
		cidp5 = cid5(i, j, k)

		pidp = pid(i, j, k)

		if( cidp0 /= 0 ) then
			p0 = p(i, j, k)
			p0 = p(i, j, k) + (p(i, j, k) - p(i+1, j, k))*d0*d0/(d0 + 0.5d0)
			p0 = p(i, j, k) + (p(i, j, k) - p(i+1, j, k))*d0*(d0 - 0.5d0)/(d0 + 0.5d0)
			p0 = p(i, j, k) + (p(i, j, k) - p(i+1, j, k))*0.5d0*(d0 - 0.5d0)/(d0 + 0.5d0)
			p0 = p(i, j, k) + (p(i, j, k) - p(i+1, j, k))*d0*d0/(d0 + 0.5d0)
			p0 = p(i, j, k) + (p(i, j, k) - p(i+1, j, k))*d0
			p0 = p(i, j, k) + (p(i, j, k) - p(i+1, j, k))*0.5d0
			m0 = 1.0d0
		endif

		if( cidp1 /= 0 ) then
			p1 = p(i, j, k)
			p1 = p(i, j, k) + (p(i, j, k) - p(i-1, j, k))*d1*(d1 - 0.5d0)/(d1 + 0.5d0) 
			p1 = p(i, j, k) + (p(i, j, k) - p(i-1, j, k))*0.5d0*(d1 - 0.5d0)/(d1 + 0.5d0) 
			p1 = p(i, j, k) + (p(i, j, k) - p(i-1, j, k))*d1*d1/(d1 + 0.5d0)
			p1 = p(i, j, k) + (p(i, j, k) - p(i-1, j, k))*d1 
			p1 = p(i, j, k) + (p(i, j, k) - p(i-1, j, k))*0.5d0
			m1 = 1.0d0
		endif

		if( cidp2 /= 0 ) then
			p2 = p(i, j, k)
			p2 = p(i, j, k) + (p(i, j, k) - p(i, j+1, k))*d2*(d2 - 0.5d0)/(d2 + 0.5d0)
			p2 = p(i, j, k) + (p(i, j, k) - p(i, j+1, k))*0.5d0*(d2 - 0.5d0)/(d2 + 0.5d0)
			p2 = p(i, j, k) + (p(i, j, k) - p(i, j+1, k))*d2*d2/(d2 + 0.5d0)
			p2 = p(i, j, k) + (p(i, j, k) - p(i, j+1, k))*d2
			p2 = p(i, j, k) + (p(i, j, k) - p(i, j+1, k))*0.5d0
			m2 = 1.0d0
		endif

		if( cidp3 /= 0 ) then
			p3 = p(i, j, k)
			p3 = p(i, j, k) + (p(i, j, k) - p(i, j-1, k))*d3*(d3 - 0.5d0)/(d3 + 0.5d0)
			p3 = p(i, j, k) + (p(i, j, k) - p(i, j-1, k))*0.5d0*(d3 - 0.5d0)/(d3 + 0.5d0)
			p3 = p(i, j, k) + (p(i, j, k) - p(i, j-1, k))*d3*d3/(d3 + 0.5d0)
			p3 = p(i, j, k) + (p(i, j, k) - p(i, j-1, k))*d3
			p3 = p(i, j, k) + (p(i, j, k) - p(i, j-1, k))*0.5d0
			m3 = 1.0d0
		endif

		if( cidp4 /= 0 ) then
			p4 = p(i, j, k)
			p4 = p(i, j, k) + (p(i, j, k) - p(i, j, k+1))*d4*(d4 - 0.5d0)/(d4 + 0.5d0)
			p4 = p(i, j, k) + (p(i, j, k) - p(i, j, k+1))*0.5d0*(d4 - 0.5d0)/(d4 + 0.5d0)
			p4 = p(i, j, k) + (p(i, j, k) - p(i, j, k+1))*d4*d4/(d4 + 0.5d0)
			p4 = p(i, j, k) + (p(i, j, k) - p(i, j, k+1))*d4
			p4 = p(i, j, k) + (p(i, j, k) - p(i, j, k+1))*0.5d0
			m4 = 1.0d0
		endif

		if( cidp5 /= 0 ) then
			p5 = p(i, j, k)
			p5 = p(i, j, k) + (p(i, j, k) - p(i, j, k-1))*d5*(d5 - 0.5d0)/(d5 + 0.5d0)
			p5 = p(i, j, k) + (p(i, j, k) - p(i, j, k-1))*0.5d0*(d5 - 0.5d0)/(d5 + 0.5d0)
			p5 = p(i, j, k) + (p(i, j, k) - p(i, j, k-1))*d5*d5/(d5 + 0.5d0)
			p5 = p(i, j, k) + (p(i, j, k) - p(i, j, k-1))*d5
			p5 = p(i, j, k) + (p(i, j, k) - p(i, j, k-1))*0.5d0
			m5 = 1.0d0
		endif

		if( cidp0 /= 0 .and. cidp1 /= 0 ) then
			p0 = p(i, j, k)
			p1 = p(i, j, k)
			m0 = 1.0d0
			m1 = 1.0d0
		endif

		if( cidp2 /= 0 .and. cidp3 /= 0 ) then
			p2 = p(i, j, k)
			p3 = p(i, j, k)
			m2 = 1.0d0
			m3 = 1.0d0
		endif

		if( cidp4 /= 0 .and. cidp5 /= 0 ) then
			p4 = p(i, j, k)
			p5 = p(i, j, k)
			m4 = 1.0d0
			m5 = 1.0d0
		endif


		if( pidp /= 1 ) then
			m0 = 0.0d0
			m1 = 0.0d0
			m2 = 0.0d0
			m3 = 0.0d0
			m4 = 0.0d0
			m5 = 0.0d0
		endif

		fspx(i, j, k) = (m1*p1 - m0*p0)
		fspy(i, j, k) = (m3*p3 - m2*p2)
		fspz(i, j, k) = (m5*p5 - m4*p4)

!		fspz(i, j, k) = (m0 + m1 + m2 + m3 + m4 + m5)*dx*dx

		fsx = fsx + fspx(i, j, k)*dx*dx
		fsy = fsy + fspy(i, j, k)*dx*dx
		fsz = fsz + fspz(i, j, k)*dx*dx
	end do
	end do
	end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
	fsp(1) = fsx
	fsp(2) = fsy
	fsp(3) = fsz
end subroutine bcut_calc_f_p

subroutine bcut_calc_f_v( &
								fsvx, &
								fsvy, &
								fsvz, &
								fsv, &
								ux, &
								uy, &
								uz, &
								c0, c1, c2, c3, c4, c5, &
								cid0, cid1, cid2, cid3, cid4, cid5, &
								pid, &
								rho, &
								mu, &
								dx, dt, &
								Us, &
								sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: fsvx, fsvy, fsvz
	real, dimension(1:3)		:: fsv
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: ux, uy, uz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: c0, c1, c2, c3, c4, c5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: cid0, cid1, cid2, cid3, cid4, cid5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: pid
	real										:: rho
	real										:: mu
  real                    :: dx, dt
	real										:: Us
	integer									:: cidp
	integer									:: cidp0, cidp1, cidp2, cidp3, cidp4, cidp5
	integer									:: pidp, pidw, pide, pids, pidn, pidb, pidt
	real										:: d0, d1, d2, d3, d4, d5
	real										:: m0, m1, m2, m3, m4, m5
	real										:: uxp, uxw, uxe, uxs, uxn, uxb, uxt
	real										:: uyp, uyw, uye, uys, uyn, uyb, uyt
	real										:: uzp, uzw, uze, uzs, uzn, uzb, uzt
	real										:: fsx, fsy, fsz
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	fsx = 0.0
	fsy = 0.0
	fsz = 0.0
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) &
!$omp					 private(cidp) &
!$omp					 private(cidp0, cidp1, cidp2, cidp3, cidp4, cidp5) &
!$omp					 private(pidp, pidw, pide, pids, pidn, pidb, pidt) &
!$omp					 private(d0, d1, d2, d3, d4, d5) &
!$omp					 private(m0, m1, m2, m3, m4, m5) &
!$omp					 private(uxp, uxw, uxe, uxs, uxn, uxb, uxt) &
!$omp					 private(uyp, uyw, uye, uys, uyn, uyb, uyt) &
!$omp					 private(uzp, uzw, uze, uzs, uzn, uzb, uzt) 
!$omp do schedule(dynamic, 1), &
!$omp		 reduction(+:fsx, fsy, fsz)
#else
#endif
	do k=1, kx
	do j=1, jx
!ocl nouxsimd
	do i=1, ix
		uxp = ux(i, j, k)
		uxw = ux(i-1, j, k)
		uxe = ux(i+1, j, k)
		uxs = ux(i, j-1, k)
		uxn = ux(i, j+1, k)
		uxb = ux(i, j, k-1)
		uxt = ux(i, j, k+1)

		uyp = uy(i, j, k)
		uyw = uy(i-1, j, k)
		uye = uy(i+1, j, k)
		uys = uy(i, j-1, k)
		uyn = uy(i, j+1, k)
		uyb = uy(i, j, k-1)
		uyt = uy(i, j, k+1)

		uzp = uz(i, j, k)
		uzw = uz(i-1, j, k)
		uze = uz(i+1, j, k)
		uzs = uz(i, j-1, k)
		uzn = uz(i, j+1, k)
		uzb = uz(i, j, k-1)
		uzt = uz(i, j, k+1)

		m0 = 0.0d0
		m1 = 0.0d0
		m2 = 0.0d0
		m3 = 0.0d0
		m4 = 0.0d0
		m5 = 0.0d0

		d0 = c0(i, j, k)
		d1 = c1(i, j, k)
		d2 = c2(i, j, k)
		d3 = c3(i, j, k)
		d4 = c4(i, j, k)
		d5 = c5(i, j, k)

		cidp0 = cid0(i, j, k)
		cidp1 = cid1(i, j, k)
		cidp2 = cid2(i, j, k)
		cidp3 = cid3(i, j, k)
		cidp4 = cid4(i, j, k)
		cidp5 = cid5(i, j, k)

		pidp = pid(i, j, k)

		if( cidp0 /= 0 ) then
			uxw  = (1.0d0 - 1.0d0/d0)*uxp 
			uyw  = (1.0d0 - 1.0d0/d0)*uyp 
			uzw  = (1.0d0 - 1.0d0/d0)*uzp 
			m0 = 1.0d0
		endif
		if( cidp1 /= 0 ) then
			uxe  = (1.0d0 - 1.0d0/d1)*uxp
			uye  = (1.0d0 - 1.0d0/d1)*uyp
			uze  = (1.0d0 - 1.0d0/d1)*uzp
			m1 = 1.0d0
		endif
		if( cidp2 /= 0 ) then
			uxs  = (1.0d0 - 1.0d0/d2)*uxp 
			uys  = (1.0d0 - 1.0d0/d2)*uyp 
			uzs  = (1.0d0 - 1.0d0/d2)*uzp 
			m2 = 1.0d0
		endif
		if( cidp3 /= 0 ) then
			uxn  = (1.0d0 - 1.0d0/d3)*uxp
			uyn  = (1.0d0 - 1.0d0/d3)*uyp
			uzn  = (1.0d0 - 1.0d0/d3)*uzp
			m3 = 1.0d0
		endif
		if( cidp4 /= 0 ) then
			uxb  = (1.0d0 - 1.0d0/d4)*uxp
			uyb  = (1.0d0 - 1.0d0/d4)*uyp
			uzb  = (1.0d0 - 1.0d0/d4)*uzp
			m4 = 1.0d0
		endif
		if( cidp5 /= 0 ) then
			uxt  = (1.0d0 - 1.0d0/d5)*uxp
			uyt  = (1.0d0 - 1.0d0/d5)*uyp
			uzt  = (1.0d0 - 1.0d0/d5)*uzp
			m5 = 1.0d0
		endif

		if( pidp /= 1 ) then
			m0 = 0.0d0
			m1 = 0.0d0
			m2 = 0.0d0
			m3 = 0.0d0
			m4 = 0.0d0
			m5 = 0.0d0
		endif

!		fsvx(i, j, k) = + mu*((uxe - uxw)*0.5)/dx*m0*2.0 &
!										- mu*((uxe - uxw)*0.5)/dx*m1*2.0 &
!										+ mu*((uxn - uxs)*0.5 + (uye - uyw)*0.5)/dx*m2 &
!										- mu*((uxn - uxs)*0.5 + (uye - uyw)*0.5)/dx*m3 &
!										+ mu*((uxt - uxb)*0.5 + (uze - uzw)*0.5)/dx*m4 &
!										- mu*((uxt - uxb)*0.5 + (uze - uzw)*0.5)/dx*m5
!		fsvy(i, j, k) = + mu*((uye - uyw)*0.5 + (uxn - uxs)*0.5)/dx*m0 &
!										- mu*((uye - uyw)*0.5 + (uxn - uxs)*0.5)/dx*m1 &
!										+ mu*((uyn - uys)*0.5)/dx*m2*2.0 &
!										- mu*((uyn - uys)*0.5)/dx*m3*2.0 &
!										+ mu*((uyt - uyb)*0.5 + (uzn - uzs)*0.5)/dx*m4 &
!										- mu*((uyt - uyb)*0.5 + (uzn - uzs)*0.5)/dx*m5
!		fsvz(i, j, k) = + mu*((uze - uzw)*0.5 + (uxt - uxb)*0.5)/dx*m0 &
!										- mu*((uze - uzw)*0.5 + (uxt - uxb)*0.5)/dx*m1 &
!										+ mu*((uzn - uzs)*0.5 + (uyt - uyb)*0.5)/dx*m2 &
!										- mu*((uzn - uzs)*0.5 + (uyt - uyb)*0.5)/dx*m3 &
!										+ mu*((uzt - uzb)*0.5)/dx*m4*2.0 &
!										- mu*((uzt - uzb)*0.5)/dx*m5*2.0

		fsvx(i, j, k) = + mu*((uxp - uxw))/dx*m0*2.0 &
										- mu*((uxe - uxp))/dx*m1*2.0 &
										+ mu*((uxp - uxs) + (uye - uyw)*0.5)/dx*m2 &
										- mu*((uxn - uxp) + (uye - uyw)*0.5)/dx*m3 &
										+ mu*((uxp - uxb) + (uze - uzw)*0.5)/dx*m4 &
										- mu*((uxt - uxp) + (uze - uzw)*0.5)/dx*m5
		fsvy(i, j, k) = + mu*((uyp - uyw) + (uxn - uxs)*0.5)/dx*m0 &
										- mu*((uye - uyp) + (uxn - uxs)*0.5)/dx*m1 &
										+ mu*((uyp - uys))/dx*m2*2.0 &
										- mu*((uyn - uyp))/dx*m3*2.0 &
										+ mu*((uyp - uyb) + (uzn - uzs)*0.5)/dx*m4 &
										- mu*((uyt - uyp) + (uzn - uzs)*0.5)/dx*m5
		fsvz(i, j, k) = + mu*((uzp - uzw) + (uxt - uxb)*0.5)/dx*m0 &
										- mu*((uze - uzp) + (uxt - uxb)*0.5)/dx*m1 &
										+ mu*((uzp - uzs) + (uyt - uyb)*0.5)/dx*m2 &
										- mu*((uzn - uzp) + (uyt - uyb)*0.5)/dx*m3 &
										+ mu*((uzp - uzb))/dx*m4*2.0 &
										- mu*((uzt - uzp))/dx*m5*2.0

!		fsvz(i, j, k) = (m0 + m1 + m2 + m3 + m4 + m5)*dx*dx

		fsx = fsx + fsvx(i, j, k)*dx*dx
		fsy = fsy + fsvy(i, j, k)*dx*dx
		fsz = fsz + fsvz(i, j, k)*dx*dx
	end do
	end do
	end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
	fsv(1) = fsx
	fsv(2) = fsy
	fsv(3) = fsz
end subroutine bcut_calc_f_v

subroutine bcut_calc_f_v_2( &
								fsvx, &
								fsvy, &
								fsvz, &
								fsv, &
								ux, &
								uy, &
								uz, &
								c0, c1, c2, c3, c4, c5, &
								cid0, cid1, cid2, cid3, cid4, cid5, &
								pid, &
								rho, &
								mu, &
								dx, dt, &
								Us, &
								sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: fsvx, fsvy, fsvz
	real, dimension(1:3)		:: fsv
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: ux, uy, uz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: c0, c1, c2, c3, c4, c5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: cid0, cid1, cid2, cid3, cid4, cid5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: pid
	real										:: rho
	real										:: mu
  real                    :: dx, dt
	real										:: Us
	integer									:: cidp
	integer									:: cidp0, cidp1, cidp2, cidp3, cidp4, cidp5
	integer									:: pidp, pidw, pide, pids, pidn, pidb, pidt
	real										:: d0, d1, d2, d3, d4, d5
	real										:: mu0, mu1, mu2, mu3, mu4, mu5
	real										:: m0, m1, m2, m3, m4, m5
	real										:: fsx, fsy, fsz
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	fsx = 0.0
	fsy = 0.0
	fsz = 0.0
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) &
!$omp					 private(cidp) &
!$omp					 private(cidp0, cidp1, cidp2, cidp3, cidp4, cidp5) &
!$omp					 private(pidp, pidw, pide, pids, pidn, pidb, pidt) &
!$omp					 private(d0, d1, d2, d3, d4, d5) &
!$omp					 private(mu0, mu1, mu2, mu3, mu4, mu5) &
!$omp					 private(m0, m1, m2, m3, m4, m5) 
!$omp do schedule(dynamic, 1), &
!$omp		 reduction(+:fsx, fsy, fsz)
#else
#endif
	do k=1, kx
	do j=1, jx
!ocl nouxsimd
	do i=1, ix
		mu0 = mu
		mu1 = mu
		mu2 = mu
		mu3 = mu
		mu4 = mu
		mu5 = mu

		m0 = 0.0d0
		m1 = 0.0d0
		m2 = 0.0d0
		m3 = 0.0d0
		m4 = 0.0d0
		m5 = 0.0d0

		d0 = c0(i, j, k)
		d1 = c1(i, j, k)
		d2 = c2(i, j, k)
		d3 = c3(i, j, k)
		d4 = c4(i, j, k)
		d5 = c5(i, j, k)

		cidp0 = cid0(i, j, k)
		cidp1 = cid1(i, j, k)
		cidp2 = cid2(i, j, k)
		cidp3 = cid3(i, j, k)
		cidp4 = cid4(i, j, k)
		cidp5 = cid5(i, j, k)

		pidp = pid(i, j, k)

		if( cidp0 /= 0 ) then
			mu0 = mu/d0
			m0 = 1.0d0
		endif

		if( cidp1 /= 0 ) then
			mu1 = mu/d1
			m1 = 1.0d0
		endif

		if( cidp2 /= 0 ) then
			mu2 = mu/d2
			m2 = 1.0d0
		endif

		if( cidp3 /= 0 ) then
			mu3 = mu/d3
			m3 = 1.0d0
		endif

		if( cidp4 /= 0 ) then
			mu4 = mu/d4
			m4 = 1.0d0
		endif

		if( cidp5 /= 0 ) then
			mu5 = mu/d5
			m5 = 1.0d0
		endif

		if( pidp /= 1 ) then
			m0 = 0.0d0
			m1 = 0.0d0
			m2 = 0.0d0
			m3 = 0.0d0
			m4 = 0.0d0
			m5 = 0.0d0
		endif

		fsvx(i, j, k) = mu0*(ux(i, j, k))*dx*m0 &
									+ mu1*(ux(i, j, k))*dx*m1 &
									+ mu2*(ux(i, j, k))*dx*m2 &
									+ mu3*(ux(i, j, k))*dx*m3 &
									+ mu4*(ux(i, j, k))*dx*m4 &
									+ mu5*(ux(i, j, k))*dx*m5
		fsvy(i, j, k) = mu0*(uy(i, j, k))*dx*m0 &
									+ mu1*(uy(i, j, k))*dx*m1 &
									+ mu2*(uy(i, j, k))*dx*m2 &
									+ mu3*(uy(i, j, k))*dx*m3 &
									+ mu4*(uy(i, j, k))*dx*m4 &
									+ mu5*(uy(i, j, k))*dx*m5
		fsvz(i, j, k) = mu0*(uz(i, j, k))*dx*m0 &
									+ mu1*(uz(i, j, k))*dx*m1 &
									+ mu2*(uz(i, j, k))*dx*m2 &
									+ mu3*(uz(i, j, k))*dx*m3 &
									+ mu4*(uz(i, j, k))*dx*m4 &
									+ mu5*(uz(i, j, k))*dx*m5


		fsvx(i, j, k) = mu0*(ux(i, j, k))/dx*m0*2.0 &
									+ mu1*(ux(i, j, k))/dx*m1*2.0 &
									+ mu2*(ux(i, j, k))/dx*m2 &
									+ mu3*(ux(i, j, k))/dx*m3 &
									+ mu4*(ux(i, j, k))/dx*m4 &
									+ mu5*(ux(i, j, k))/dx*m5 
!									+ mu*(uy(i+1, j, k) - uy(i-1, j, k))*0.5/dx*m2 &
!									- mu*(uy(i+1, j, k) - uy(i-1, j, k))*0.5/dx*m3 &
!									+ mu*(uz(i+1, j, k) - uz(i-1, j, k))*0.5/dx*m4 &
!									- mu*(uz(i+1, j, k) - uz(i-1, j, k))*0.5/dx*m5
		fsvy(i, j, k) = mu0*(uy(i, j, k))/dx*m0 &
									+ mu1*(uy(i, j, k))/dx*m1 &
									+ mu2*(uy(i, j, k))/dx*m2*2.0 &
									+ mu3*(uy(i, j, k))/dx*m3*2.0 &
									+ mu4*(uy(i, j, k))/dx*m4 &
									+ mu5*(uy(i, j, k))/dx*m5 
!									+ mu*(ux(i, j+1, k) - ux(i, j-1, k))*0.5/dx*m0 &
!									- mu*(ux(i, j+1, k) - ux(i, j-1, k))*0.5/dx*m1 &
!									+ mu*(uz(i, j+1, k) - uz(i, j-1, k))*0.5/dx*m4 &
!									- mu*(uz(i, j+1, k) - uz(i, j-1, k))*0.5/dx*m5
		fsvz(i, j, k) = mu0*(uz(i, j, k))/dx*m0 &
									+ mu1*(uz(i, j, k))/dx*m1 &
									+ mu2*(uz(i, j, k))/dx*m2 &
									+ mu3*(uz(i, j, k))/dx*m3 &
									+ mu4*(uz(i, j, k))/dx*m4*2.0 &
									+ mu5*(uz(i, j, k))/dx*m5*2.0 
!									+ mu*(ux(i, j, k+1) - ux(i, j, k-1))*0.5/dx*m0 &
!									- mu*(ux(i, j, k+1) - ux(i, j, k-1))*0.5/dx*m1 &
!									+ mu*(uy(i, j, k+1) - uy(i, j, k-1))*0.5/dx*m2 &
!									- mu*(uy(i, j, k+1) - uy(i, j, k-1))*0.5/dx*m3

!		fsvz(i, j, k) = (m0 + m1 + m2 + m3 + m4 + m5)*dx*dx

		fsx = fsx + fsvx(i, j, k)*dx*dx
		fsy = fsy + fsvy(i, j, k)*dx*dx
		fsz = fsz + fsvz(i, j, k)*dx*dx
	end do
	end do
	end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
	fsv(1) = fsx
	fsv(2) = fsy
	fsv(3) = fsz
end subroutine bcut_calc_f_v_2

subroutine bcut_calc_q( &
								qx, &
								qy, &
								qz, &
								q, &
								t, &
								c0, c1, c2, c3, c4, c5, &
								cid0, cid1, cid2, cid3, cid4, cid5, &
								pid, &
								rhof, &
								cpf, &
								kf, &
								dx, dt, &
								Tc, &
								sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: qx, qy, qz
	real, dimension(1:3)		:: q
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: t
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: c0, c1, c2, c3, c4, c5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: cid0, cid1, cid2, cid3, cid4, cid5
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: pid
	real										:: rhof
	real										:: cpf
	real										:: kf
  real                    :: dx, dt
	real										:: Tc
	real										:: tp, tw, te, ts, tn, tb, tt
	integer									:: cidp
	integer									:: cidp0, cidp1, cidp2, cidp3, cidp4, cidp5
	integer									:: pidp, pidw, pide, pids, pidn, pidb, pidt
	real										:: d0, d1, d2, d3, d4, d5
	real										:: m0, m1, m2, m3, m4, m5
	real										:: qx0, qy0, qz0
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	qx0 = 0.0
	qy0 = 0.0
	qz0 = 0.0
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) 
!$omp do schedule(static, 1), &
!$omp		 reduction(+:qx0, qy0, qz0)
#else
#endif
	do k=1, kx
	do j=1, jx
!ocl nouxsimd
	do i=1, ix
		tp = t(i, j, k)
		tw = t(i-1, j, k)
		te = t(i+1, j, k)
		ts = t(i, j-1, k)
		tn = t(i, j+1, k)
		tb = t(i, j, k-1)
		tt = t(i, j, k+1)

		m0 = 0.0d0
		m1 = 0.0d0
		m2 = 0.0d0
		m3 = 0.0d0
		m4 = 0.0d0
		m5 = 0.0d0

		d0 = c0(i, j, k)
		d1 = c1(i, j, k)
		d2 = c2(i, j, k)
		d3 = c3(i, j, k)
		d4 = c4(i, j, k)
		d5 = c5(i, j, k)

		cidp0 = cid0(i, j, k)
		cidp1 = cid1(i, j, k)
		cidp2 = cid2(i, j, k)
		cidp3 = cid3(i, j, k)
		cidp4 = cid4(i, j, k)
		cidp5 = cid5(i, j, k)

		pidp = pid(i, j, k)

		if( cidp0 /= 0 ) then
			tw  = (1.0d0 - 1.0d0/d0)*tp + (1.0d0/d0)*Tc
			m0 = 1.0d0
		endif
		if( cidp1 /= 0 ) then
			te  = (1.0d0 - 1.0d0/d1)*tp + (1.0d0/d1)*Tc
			m1 = 1.0d0
		endif
		if( cidp2 /= 0 ) then
			ts  = (1.0d0 - 1.0d0/d2)*tp + (1.0d0/d2)*Tc
			m2 = 1.0d0
		endif
		if( cidp3 /= 0 ) then
			tn  = (1.0d0 - 1.0d0/d3)*tp + (1.0d0/d3)*Tc
			m3 = 1.0d0
		endif
		if( cidp4 /= 0 ) then
			tb  = (1.0d0 - 1.0d0/d4)*tp + (1.0d0/d4)*Tc
			m4 = 1.0d0
		endif
		if( cidp5 /= 0 ) then
			tt  = (1.0d0 - 1.0d0/d5)*tp + (1.0d0/d5)*Tc
			m5 = 1.0d0
		endif

		if( pidp /= 1 ) then
			m0 = 0.0d0
			m1 = 0.0d0
			m2 = 0.0d0
			m3 = 0.0d0
			m4 = 0.0d0
			m5 = 0.0d0
		endif

		qx(i, j, k) = - kf*(tp - tw)/dx*m0 &
									+ kf*(te - tp)/dx*m1
		qy(i, j, k) = - kf*(tp - ts)/dx*m2 &
									+ kf*(tn - tp)/dx*m3
		qz(i, j, k) = - kf*(tp - tb)/dx*m4 &
									+ kf*(tt - tp)/dx*m5

		qx0 = qx0 + qx(i, j, k)
		qy0 = qy0 + qy(i, j, k)
		qz0 = qz0 + qz(i, j, k)
	end do
	end do
	end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
	q(1) = qx0
	q(2) = qy0
	q(3) = qz0
end subroutine bcut_calc_q

subroutine bcut_set_a( &
								ux, uy, uz, &
								p, &
								pid, &
								rho, &
								mu, &
								a, &
								U, &
								dx, dt, &
								org, &
								sz, g)
  implicit none
  integer                 :: i, j, k
  integer                 :: ix, jx, kx
  integer                 :: g
  integer, dimension(3)   :: sz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: ux, uy, uz
  real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: p
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: pid
	real										:: rho
	real										:: mu
	real										:: a
	real										:: U
  real                    :: dx, dt
	real, dimension(3)			:: org
	real										:: x, y, z
	real										:: r
	integer									:: pidp
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) &
!$omp					 private(x, y, z) &
!$omp					 private(r) &
!$omp					 private(pidp)
!$omp do schedule(static, 1)
#else
#endif
  do k=1, kx
  do j=1, jx
!ocl nouxsimd
  do i=1, ix
		x = org(1) + (real(i) - 0.5)*dx
		y = org(2) + (real(j) - 0.5)*dx
		z = org(3) + (real(k) - 0.5)*dx
		r = sqrt(x*x + y*y + z*z)
		pidp = pid(i, j, k)
!		if( pidp == 1 ) then
		if( r > 0.3 ) then
			p(i, j, k) = -1.5*mu*a*U*x/(r*r*r)
			ux(i, j, k) = U - 0.25*a*U/r*(3.0 + a*a/(r*r)) - 0.75*a*U*x*x/(r*r*r)*(1.0 - a*a/(r*r))
			uy(i, j, k) =   - 0.75*a*U*x*y/(r*r*r)*(1.0 - a*a/(r*r))
			uz(i, j, k) =   - 0.75*a*U*x*z/(r*r*r)*(1.0 - a*a/(r*r))
		endif
  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bcut_set_a

subroutine bcut_set_fluidseed( &
								pid, &
								xs, ys, zs, &
								dx, &
								org, &
								sz, g)
  implicit none
  integer                 :: i, j, k
  integer                 :: ix, jx, kx
  integer                 :: g
  integer, dimension(3)   :: sz
  integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: pid
	real										:: xs, ys, zs
	real										:: dx
	real, dimension(3)			:: org
	real										:: x0, y0, z0
	real										:: x1, y1, z1
	integer									:: pidp
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) &
!$omp					 private(x0, y0, z0) &
!$omp					 private(x1, y1, z1) &
!$omp					 private(pidp)
!$omp do schedule(static, 1)
#else
#endif
  do k=1, kx
  do j=1, jx
!ocl nouxsimd
  do i=1, ix
		x0 = org(1) + (real(i-1))*dx
		y0 = org(2) + (real(j-1))*dx
		z0 = org(3) + (real(k-1))*dx
		x1 = org(1) + (real(i))*dx
		y1 = org(2) + (real(j))*dx
		z1 = org(3) + (real(k))*dx

		if( (x0 <= xs .and. xs <= x1) .and. &
				(y0 <= ys .and. ys <= y1) .and. &
				(z0 <= zs .and. zs <= z1) ) then
			pid(i, j, k) = 1
!			write(*, *) i, j, k
!			write(*, *) x0, xs, x1
!			write(*, *) y0, ys, y1
!			write(*, *) z0, zs, z1
		endif
  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bcut_set_fluidseed

