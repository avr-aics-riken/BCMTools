!
! BCMTools
!
! Copyright (C) 2011-2014 Institute of Industrial Science, The University of Tokyo.
! All rights reserved.
!
! Copyright (c) 2012-2015 Advanced Institute for Computational Science, RIKEN.
! All rights reserved.
!


subroutine bc_Aw_d(Ap, Aw, b, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: Ap, Aw
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: b
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
		Ap(i, j, k) = Ap(i, j, k) - Aw(i, j, k)
		 b(i, j, k) = b(i, j, k) - Aw(i, j, k)*xc*2.0d0
		Aw(i, j, k) = 0.0d0
	end do
	end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bc_Aw_d

subroutine bc_Ae_d(Ap, Ae, b, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: Ap, Ae
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: b
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
		Ap(i, j, k) = Ap(i, j, k) - Ae(i, j, k)
		 b(i, j, k) = b(i, j, k) - Ae(i, j, k)*xc*2.0d0
		Ae(i, j, k) = 0.0d0
	end do
	end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bc_Ae_d

subroutine bc_As_d(Ap, As, b, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: Ap, As
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: b
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
		Ap(i, j, k) = Ap(i, j, k) - As(i, j, k)
		 b(i, j, k) = b(i, j, k) - As(i, j, k)*xc*2.0d0
		As(i, j, k) = 0.0d0
	end do
	end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bc_As_d

subroutine bc_An_d(Ap, An, b, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: Ap, An
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: b
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
		Ap(i, j, k) = Ap(i, j, k) - An(i, j, k)
		 b(i, j, k) = b(i, j, k) - An(i, j, k)*xc*2.0d0
		An(i, j, k) = 0.0d0
	end do
	end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bc_An_d

subroutine bc_Ab_d(Ap, Ab, b, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: Ap, Ab
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: b
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
		Ap(i, j, k) = Ap(i, j, k) - Ab(i, j, k)
		 b(i, j, k) = b(i, j, k) - Ab(i, j, k)*xc*2.0d0
		Ab(i, j, k) = 0.0d0
	end do
	end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bc_Ab_d

subroutine bc_At_d(Ap, At, b, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: Ap, At
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: b
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
		Ap(i, j, k) = Ap(i, j, k) - At(i, j, k)
		 b(i, j, k) = b(i, j, k) - At(i, j, k)*xc*2.0d0
		At(i, j, k) = 0.0d0
	end do
	end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bc_At_d

subroutine bc_Aw_n(Ap, Aw, b, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: Ap, Aw
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: b
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
		Ap(i, j, k) = Ap(i, j, k) + Aw(i, j, k)
		 b(i, j, k) = b(i, j, k) + Aw(i, j, k)*xc
		Aw(i, j, k) = 0.0d0
	end do
	end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bc_Aw_n

subroutine bc_Ae_n(Ap, Ae, b, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: Ap, Ae
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: b
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
		Ap(i, j, k) = Ap(i, j, k) + Ae(i, j, k)
		 b(i, j, k) = b(i, j, k) - Ae(i, j, k)*xc
		Ae(i, j, k) = 0.0d0
	end do
	end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bc_Ae_n

subroutine bc_As_n(Ap, As, b, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: Ap, As
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: b
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
		Ap(i, j, k) = Ap(i, j, k) + As(i, j, k)
		 b(i, j, k) = b(i, j, k) + As(i, j, k)*xc
		As(i, j, k) = 0.0d0
	end do
	end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bc_As_n

subroutine bc_An_n(Ap, An, b, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: Ap, An
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: b
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
		Ap(i, j, k) = Ap(i, j, k) + An(i, j, k)
		 b(i, j, k) = b(i, j, k) - An(i, j, k)*xc
		An(i, j, k) = 0.0d0
	end do
	end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bc_An_n

subroutine bc_Ab_n(Ap, Ab, b, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: Ap, Ab
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: b
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
		Ap(i, j, k) = Ap(i, j, k) + Ab(i, j, k)
		 b(i, j, k) = b(i, j, k) + Ab(i, j, k)*xc
		Ab(i, j, k) = 0.0d0
	end do
	end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bc_Ab_n

subroutine bc_At_n(Ap, At, b, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: Ap, At
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: b
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
		Ap(i, j, k) = Ap(i, j, k) + At(i, j, k)
		 b(i, j, k) = b(i, j, k) - At(i, j, k)*xc
		At(i, j, k) = 0.0d0
	end do
	end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bc_At_n

