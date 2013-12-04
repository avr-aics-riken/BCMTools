function bc_get_u(x0, y0, z0, um)
	implicit none
	real					:: bc_get_u
	real					:: x0, y0, z0
	real					:: x1, y1, z1
	real					:: rx, ry, rz, r2, r
	real					:: um
	real					:: u
	bc_get_u = 0.0

! for poiseuille.conf
	x1 = x0
	y1 = 0.0
	z1 = 0.0

! for elbow-6a.conf
	x1 = x0
	y1 = 0.0
	z1 = -3.0

	rx = x0 - x1
	ry = y0 - y1
	rz = z0 - z1
	r2 = rx*rx + ry*ry + rz*rz
	r = sqrt(r2)

	u = um*2.0*(1.0 - r*r/(0.5*0.5))
	if( r > 0.5 ) then
		u = 0.0
	end if

	bc_get_u = u
	return
end function bc_get_u

function bc_get_gp(x0, y0, z0, gp)
	implicit none
	real					:: bc_get_gp
	real					:: x0, y0, z0
	real					:: x1, y1, z1
	real					:: rx, ry, rz, r2, r
	real					:: gp
	bc_get_gp = 0.0

! for poiseuille.conf
	x1 = x0
	y1 = 0.0
	z1 = 0.0

! for elbow-6a.conf
	x1 = x0
	y1 = 0.0
	z1 = -3.0

	rx = x0 - x1
	ry = y0 - y1
	rz = z0 - z1
	r2 = rx*rx + ry*ry + rz*rz
	r = sqrt(r2)

	bc_get_gp = gp
	if( r > 0.5 ) then
		bc_get_gp = 0.0d0
	end if

	return
end function bc_get_gp

subroutine bc_x3_poiseuille_u(x, xc, sz, g, org, bsize, csize)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(3)			:: org, bsize, csize
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real										:: xc
	real										:: dx
	real										:: x0, y0, z0
	real										:: x1, y1, z1
	real										:: rx, ry, rz, r2, r
	real										:: u
	real										:: bc_get_u
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	i = 1
#ifdef _LARGE_BLOCK_
!$omp parallel private(j, k) &
!$omp					 private(dx) &
!$omp					 private(x0, y0, z0) &
!$omp					 private(x1, y1, z1) &
!$omp					 private(rx, ry, rz, r2, r) &
!$omp					 private(u)
!$omp do
#else
#endif
	do k=1, kx
	do j=1, jx
		dx = csize(1)
		x0 = org(1) + (real(i) - 0.5)*dx
		y0 = org(2) + (real(j) - 0.5)*dx
		z0 = org(3) + (real(k) - 0.5)*dx
		u = bc_get_u(x0, y0, z0, xc)
		x(i-1, j, k) = 2.0d0*u - x(i, j, k)
		x(i-2, j, k) = u
	end do
	end do
#ifdef _LARGE_BLOCK_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bc_x3_poiseuille_u

subroutine bc_Aw_poiseuille_u(Ap, Aw, b, xc, sz, g, org, bsize, csize)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(3)			:: org, bsize, csize
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: Ap, Aw
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: b
	real										:: xc
	real										:: dx
	real										:: x0, y0, z0
	real										:: x1, y1, z1
	real										:: rx, ry, rz, r2, r
	real										:: u
	real										:: bc_get_u
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	i = 1
#ifdef _LARGE_BLOCK_
!$omp parallel private(j, k) &
!$omp					 private(dx) &
!$omp					 private(x0, y0, z0) &
!$omp					 private(x1, y1, z1) &
!$omp					 private(rx, ry, rz, r2, r) &
!$omp					 private(u)
!$omp do
#else
#endif
	do k=1, kx
	do j=1, jx
		dx = csize(1)
		x0 = org(1) + (real(i) - 0.5)*dx
		y0 = org(2) + (real(j) - 0.5)*dx
		z0 = org(3) + (real(k) - 0.5)*dx
		u = bc_get_u(x0, y0, z0, xc)
		Ap(i, j, k) = Ap(i, j, k) - Aw(i, j, k)
		 b(i, j, k) = b(i, j, k) - Aw(i, j, k)*u*2.0d0
		Aw(i, j, k) = 0.0d0
	end do
	end do
#ifdef _LARGE_BLOCK_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bc_Aw_poiseuille_u

subroutine bc_x3_poiseuille_p(x, xc, sz, g, org, bsize, csize)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(3)			:: org, bsize, csize
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real										:: xc
	real										:: dx
	real										:: x0, y0, z0
	real										:: x1, y1, z1
	real										:: rx, ry, rz, r2, r
	real										:: gp
	real										:: bc_get_gp
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	i = 1
#ifdef _LARGE_BLOCK_
!$omp parallel private(j, k) &
!$omp					 private(dx) &
!$omp					 private(x0, y0, z0) &
!$omp					 private(x1, y1, z1) &
!$omp					 private(rx, ry, rz, r2, r) 
!$omp do
#else
#endif
	do k=1, kx
	do j=1, jx
		dx = csize(1)
		x0 = org(1) + (real(i) - 0.5)*dx
		y0 = org(2) + (real(j) - 0.5)*dx
		z0 = org(3) + (real(k) - 0.5)*dx
		gp = bc_get_gp(x0, y0, z0, xc*dx)
		x(i-1, j, k) = x(i, j, k) - gp
		x(i-2, j, k) = x(i, j, k) - gp*2.0d0
	end do
	end do
#ifdef _LARGE_BLOCK_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bc_x3_poiseuille_p

subroutine bc_Aw_poiseuille_p(Ap, Aw, b, xc, sz, g, org, bsize, csize)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(3)			:: org, bsize, csize
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: Ap, Aw
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: b
	real										:: xc
	real										:: dx
	real										:: x0, y0, z0
	real										:: x1, y1, z1
	real										:: rx, ry, rz, r2, r
	real										:: gp
	real										:: bc_get_gp
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	i = 1
#ifdef _LARGE_BLOCK_
!$omp parallel private(j, k) &
!$omp					 private(dx) &
!$omp					 private(x0, y0, z0) &
!$omp					 private(x1, y1, z1) &
!$omp					 private(rx, ry, rz, r2, r) 
!$omp do
#else
#endif
	do k=1, kx
	do j=1, jx
		dx = csize(1)
		x0 = org(1) + (real(i) - 0.5)*dx
		y0 = org(2) + (real(j) - 0.5)*dx
		z0 = org(3) + (real(k) - 0.5)*dx
		gp = bc_get_gp(x0, y0, z0, xc*dx)
		Ap(i, j, k) = Ap(i, j, k) + Aw(i, j, k)
		 b(i, j, k) = b(i, j, k) + Aw(i, j, k)*gp
		Aw(i, j, k) = 0.0d0
	end do
	end do
#ifdef _LARGE_BLOCK_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bc_Aw_poiseuille_p

