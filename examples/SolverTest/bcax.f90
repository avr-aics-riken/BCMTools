!
! BCMTools
!
! Copyright (C) 2011-2013 Institute of Industrial Science, The University of Tokyo.
! All rights reserved.
!
! Copyright (c) 2012-2013 Advanced Institute for Computational Science, RIKEN.
! All rights reserved.
!

function bc_get_u(x0, x1, um)
	implicit none
	real					:: bc_get_u
	real, dimension(3)	:: x0
	real, dimension(3)	:: x1
	real					:: rx, ry, rz, r2, r
	real					:: um
	real					:: u
	bc_get_u = 0.0

! for poiseuille.conf
!	x1 = x0
!	y1 = 0.0
!	z1 = 0.0

! for elbow-6a.conf
!	x1 = x0
!	y1 = 0.0
!	z1 = -3.0

	rx = x0(1) - x0(1)
	ry = x0(2) - x1(2)
	rz = x0(3) - x1(3)
	r2 = rx*rx + ry*ry + rz*rz
	r = sqrt(r2)

	u = um*2.0*(1.0 - r*r/(0.5*0.5))
	if( r > 0.5 ) then
		u = 0.0
	end if

	bc_get_u = u
	return
end function bc_get_u

function bc_get_gp(x0, x1, gp)
	implicit none
	real					:: bc_get_gp
	real, dimension(3)	:: x0
	real, dimension(3)	:: x1
	real					:: rx, ry, rz, r2, r
	real					:: gp
	bc_get_gp = 0.0

! for poiseuille.conf
!	x1 = x0
!	y1 = 0.0
!	z1 = 0.0

! for elbow-6a.conf
!	x1 = x0
!	y1 = 0.0
!	z1 = -3.0

	rx = x0(1) - x0(1)
	ry = x0(2) - x1(2)
	rz = x0(3) - x1(3)
	r2 = rx*rx + ry*ry + rz*rz
	r = sqrt(r2)

	bc_get_gp = gp
	if( r > 0.5 ) then
		bc_get_gp = 0.0d0
	end if

	return
end function bc_get_gp

subroutine bc_x3_poiseuille_u(x, xc, sz, g, center, org, bsize, csize)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(3)			:: center
	real, dimension(3)			:: org, bsize, csize
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real										:: xc
	real										:: dx
	real, dimension(3)			:: x0
	real, dimension(3)			:: x1
	real										:: u
	real										:: bc_get_u
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	i = 1
	dx = csize(1)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(j, k) &
!$omp					 private(x0) &
!$omp					 private(x1) &
!$omp					 private(u)
!$omp do
#else
#endif
	do k=1, kx
	do j=1, jx
		x0(1) = org(1) + (real(i) - 0.5)*dx
		x0(2) = org(2) + (real(j) - 0.5)*dx
		x0(3) = org(3) + (real(k) - 0.5)*dx
		x1(1) = center(1)
		x1(2) = center(2)
		x1(3) = center(3)
		u = bc_get_u(x0, x1, xc)
		x(i-1, j, k) = 2.0d0*u - x(i, j, k)
		x(i-2, j, k) = u
	end do
	end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bc_x3_poiseuille_u

subroutine bc_Aw_poiseuille_u(Ap, Aw, b, xc, sz, g, center, org, bsize, csize)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(3)			:: center
	real, dimension(3)			:: org, bsize, csize
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: Ap, Aw
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: b
	real										:: xc
	real										:: dx
	real, dimension(3)			:: x0
	real, dimension(3)			:: x1
	real										:: u
	real										:: bc_get_u
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	i = 1
	dx = csize(1)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(j, k) &
!$omp					 private(x0) &
!$omp					 private(x1) &
!$omp					 private(u)
!$omp do
#else
#endif
	do k=1, kx
	do j=1, jx
		x0(1) = org(1) + (real(i) - 0.5)*dx
		x0(2) = org(2) + (real(j) - 0.5)*dx
		x0(3) = org(3) + (real(k) - 0.5)*dx
		x1(1) = center(1)
		x1(2) = center(2)
		x1(3) = center(3)
		u = bc_get_u(x0, x1, xc)
		Ap(i, j, k) = Ap(i, j, k) - Aw(i, j, k)
		 b(i, j, k) = b(i, j, k) - Aw(i, j, k)*u*2.0d0
		Aw(i, j, k) = 0.0d0
	end do
	end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bc_Aw_poiseuille_u

subroutine bc_x3_poiseuille_p(x, xc, sz, g, center, org, bsize, csize)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(3)			:: center
	real, dimension(3)			:: org, bsize, csize
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real										:: xc
	real										:: dx
	real, dimension(3)			:: x0
	real, dimension(3)			:: x1
	real										:: gp
	real										:: bc_get_gp
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	i = 1
	dx = csize(1)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(j, k) &
!$omp					 private(x0) &
!$omp					 private(x1) &
!$omp					 private(gp) 
!$omp do
#else
#endif
	do k=1, kx
	do j=1, jx
		x0(1) = org(1) + (real(i) - 0.5)*dx
		x0(2) = org(2) + (real(j) - 0.5)*dx
		x0(3) = org(3) + (real(k) - 0.5)*dx
		x1(1) = center(1)
		x1(2) = center(2)
		x1(3) = center(3)
		gp = bc_get_gp(x0, x1, xc*dx)
		x(i-1, j, k) = x(i, j, k) - gp
		x(i-2, j, k) = x(i, j, k) - gp*2.0d0
	end do
	end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bc_x3_poiseuille_p

subroutine bc_Aw_poiseuille_p(Ap, Aw, b, xc, sz, g, center, org, bsize, csize)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(3)			:: center
	real, dimension(3)			:: org, bsize, csize
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: Ap, Aw
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: b
	real										:: xc
	real										:: dx
	real, dimension(3)			:: x0
	real, dimension(3)			:: x1
	real										:: gp
	real										:: bc_get_gp
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	i = 1
	dx = csize(1)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(j, k) &
!$omp					 private(x0) &
!$omp					 private(x1) &
!$omp					 private(gp) 
!$omp do
#else
#endif
	do k=1, kx
	do j=1, jx
		x0(1) = org(1) + (real(i) - 0.5)*dx
		x0(2) = org(2) + (real(j) - 0.5)*dx
		x0(3) = org(3) + (real(k) - 0.5)*dx
		x1(1) = center(1)
		x1(2) = center(2)
		x1(3) = center(3)
		gp = bc_get_gp(x0, x1, xc*dx)
		Ap(i, j, k) = Ap(i, j, k) + Aw(i, j, k)
		 b(i, j, k) = b(i, j, k) + Aw(i, j, k)*gp
		Aw(i, j, k) = 0.0d0
	end do
	end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bc_Aw_poiseuille_p

