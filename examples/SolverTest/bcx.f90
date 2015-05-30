!
! BCMTools
!
! Copyright (C) 2011-2014 Institute of Industrial Science, The University of Tokyo.
! All rights reserved.
!
! Copyright (c) 2012-2015 Advanced Institute for Computational Science, RIKEN.
! All rights reserved.
!

subroutine bc_x1_d(x, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real										:: xc
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	i = sz(1)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(j, k)
!$omp do
#else
#endif
	do k=1, kx
	do j=1, jx
		x(i+1, j, k) = 2.0d0*xc - x(i, j, k)
		x(i+2, j, k) = xc
	end do
	end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bc_x1_d

subroutine bc_x3_d(x, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real										:: xc
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	i = 1
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(j, k)
!$omp do
#else
#endif
	do k=1, kx
	do j=1, jx
		x(i-1, j, k) = 2.0d0*xc - x(i, j, k)
		x(i-2, j, k) = xc
	end do
	end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bc_x3_d

subroutine bc_x2_d(x, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real										:: xc
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	j = sz(2)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, k)
!$omp do
#else
#endif
	do k=1, kx
	do i=1, ix
		x(i, j+1, k) = 2.0d0*xc - x(i, j, k)
		x(i, j+2, k) = xc
	end do
	end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bc_x2_d

subroutine bc_x4_d(x, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real										:: xc
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	j = 1
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, k)
!$omp do
#else
#endif
	do k=1, kx
	do i=1, ix
		x(i, j-1, k) = 2.0d0*xc - x(i, j, k)
		x(i, j-2, k) = xc
	end do
	end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bc_x4_d

subroutine bc_x5_d(x, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real										:: xc
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	k = kx
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j)
!$omp do
#else
#endif
	do j=1, jx
	do i=1, ix
		x(i, j, k+1) = 2.0d0*xc - x(i, j, k)
		x(i, j, k+2) = xc
	end do
	end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bc_x5_d

subroutine bc_x6_d(x, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real										:: xc
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	k = 1
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j)
!$omp do
#else
#endif
	do j=1, jx
	do i=1, ix
		x(i, j, k-1) = 2.0d0*xc - x(i, j, k)
		x(i, j, k-2) = xc
	end do
	end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bc_x6_d

subroutine bc_x1_n(x, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real										:: xc
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	i = sz(1)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(j, k)
!$omp do
#else
#endif
	do k=1, kx
	do j=1, jx
		x(i+1, j, k) = x(i, j, k) + xc
		x(i+2, j, k) = x(i, j, k) + xc*2.0d0
	end do
	end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bc_x1_n

subroutine bc_x3_n(x, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real										:: xc
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	i = 1
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(j, k)
!$omp do
#else
#endif
	do k=1, kx
	do j=1, jx
		x(i-1, j, k) = x(i, j, k) - xc
		x(i-2, j, k) = x(i, j, k) - xc*2.0d0
	end do
	end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bc_x3_n

subroutine bc_x2_n(x, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real										:: xc
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	j = sz(2)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, k)
!$omp do
#else
#endif
	do k=1, kx
	do i=1, ix
		x(i, j+1, k) = x(i, j, k) + xc
		x(i, j+2, k) = x(i, j, k) + xc*2.0d0
	end do
	end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bc_x2_n

subroutine bc_x4_n(x, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real										:: xc
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	j = 1
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, k)
!$omp do
#else
#endif
	do k=1, kx
	do i=1, ix
		x(i, j-1, k) = x(i, j, k) - xc
		x(i, j-2, k) = x(i, j, k) - xc*2.0d0
	end do
	end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bc_x4_n

subroutine bc_x5_n(x, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real										:: xc
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	k = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j)
!$omp do
#else
#endif
	do j=1, jx
	do i=1, ix
		x(i, j, k+1) = x(i, j, k) + xc
		x(i, j, k+2) = x(i, j, k) + xc*2.0d0
	end do
	end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bc_x5_n

subroutine bc_x6_n(x, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real										:: xc
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	k = 1
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j)
!$omp do
#else
#endif
	do j=1, jx
	do i=1, ix
		x(i, j, k-1) = x(i, j, k) - xc
		x(i, j, k-2) = x(i, j, k) - xc*2.0d0
	end do
	end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bc_x6_n

subroutine bc_x1_d_f(x, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real										:: xc
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	i = sz(1)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(j, k)
!$omp do
#else
#endif
	do k=1, kx
	do j=1, jx
		x(i  , j, k) = xc
		x(i+1, j, k) = xc
		x(i+2, j, k) = xc
	end do
	end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bc_x1_d_f

subroutine bc_x3_d_f(x, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real										:: xc
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	i = 1
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(j, k)
!$omp do
#else
#endif
	do k=1, kx
	do j=1, jx
		x(i-1, j, k) = xc
		x(i-2, j, k) = xc
	end do
	end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bc_x3_d_f

subroutine bc_x2_d_f(x, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real										:: xc
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	j = sz(2)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, k)
!$omp do
#else
#endif
	do k=1, kx
	do i=1, ix
		x(i, j  , k) = xc
		x(i, j+1, k) = xc
		x(i, j+2, k) = xc
	end do
	end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bc_x2_d_f

subroutine bc_x4_d_f(x, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real										:: xc
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	j = 1
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, k)
!$omp do
#else
#endif
	do k=1, kx
	do i=1, ix
		x(i, j-1, k) = xc
		x(i, j-2, k) = xc
	end do
	end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bc_x4_d_f

subroutine bc_x5_d_f(x, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real										:: xc
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	k = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j)
!$omp do
#else
#endif
	do j=1, jx
	do i=1, ix
		x(i, j, k  ) = xc
		x(i, j, k+1) = xc
		x(i, j, k+2) = xc
	end do
	end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bc_x5_d_f

subroutine bc_x6_d_f(x, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real										:: xc
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	k = 1
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j)
!$omp do
#else
#endif
	do j=1, jx
	do i=1, ix
		x(i, j, k-1) = xc
		x(i, j, k-2) = xc
	end do
	end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bc_x6_d_f

subroutine bc_x1_p(x, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real										:: xc
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(j, k)
!$omp do
#else
#endif
	do k=1-g, kx+g
	do j=1-g, jx+g
		x(ix+1, j, k) = x(1, j, k)
		x(ix+2, j, k) = x(2, j, k)
		x( 0, j, k) = x(ix  , j, k)
		x(-1, j, k) = x(ix-1, j, k)
	end do
	end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bc_x1_p

subroutine bc_x2_p(x, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real										:: xc
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, k)
!$omp do
#else
#endif
	do k=1-g, kx+g
	do i=1-g, ix+g
		x(i, jx+1, k) = x(i, 1, k)
		x(i, jx+2, k) = x(i, 2, k)
		x(i,  0, k) = x(i, jx  , k)
		x(i, -1, k) = x(i, jx-1, k)
	end do
	end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bc_x2_p

subroutine bc_x5_p(x, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real										:: xc
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j)
!$omp do
#else
#endif
	do j=1-g, jx+g
	do i=1-g, ix+g
		x(i, j, kx+1) = x(i, j, 1)
		x(i, j, kx+2) = x(i, j, 2)
		x(i, j,  0) = x(i, j, kx  )
		x(i, j, -1) = x(i, j, kx-1)
	end do
	end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bc_x5_p

