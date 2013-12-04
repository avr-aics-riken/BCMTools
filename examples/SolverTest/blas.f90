subroutine fill(x, a, sz, g)
  implicit none
  integer, dimension(3)   :: sz
  integer                 :: g
  integer                 :: i, j, k
  integer                 :: ix, jx, kx
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: x
  real                    :: a 
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _LARGE_BLOCK_
!$omp parallel private(i, j, k)
!$omp do
#else
#endif
  do k=1-g, kx+g
  do j=1-g, jx+g
  do i=1-g, ix+g
    x(i, j, k) = a
  end do
  end do
  end do
#ifdef _LARGE_BLOCK_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine fill

subroutine fill_vf3d(x, a, sz, g, ne)
  implicit none
  integer, dimension(3)   :: sz
  integer                 :: g
  integer                 :: i, j, k
  integer                 :: ix, jx, kx
  integer                 :: l
  integer                 :: ne
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g, 0:ne-1)  :: x
  real                    :: a
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _LARGE_BLOCK_
!$omp parallel private(l, i, j, k)
!$omp do
#else
#endif
  do l=0, ne-1
  do k=1-g, kx+g
  do j=1-g, jx+g
  do i=1-g, ix+g
    x(i, j, k, l) = a
  end do
  end do
  end do
  end do
#ifdef _LARGE_BLOCK_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine fill_vf3d

subroutine copy(y, x, sz, g)
	implicit none
	integer, dimension(3)		:: sz
	integer									:: g
	integer									:: i, j, k
	integer									:: ix, jx, kx
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: y
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
#ifdef _LARGE_BLOCK_
!$omp parallel private(i, j, k)
!$omp do
#else
#endif
	do k=1-g, kx+g
	do j=1-g, jx+g
	do i=1-g, ix+g
		y(i, j, k) = x(i, j, k)
	end do
	end do
	end do
#ifdef _LARGE_BLOCK_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine copy

subroutine add(C, A, B, sz, g)
	implicit none
	integer, dimension(3)		:: sz
	integer									:: g
	integer									:: i, j, k
	integer									:: ix, jx, kx
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: C
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: A, B
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
#ifdef _LARGE_BLOCK_
!$omp parallel private(i, j, k)
!$omp do
#else
!ocl nouxsimd
!ocl serial
!dir$ noparallel
#endif
	do k=1, kx
	do j=1, jx
	do i=1, ix
		C(i, j, k) = A(i, j, k) + B(i, j, k)
	end do
	end do
	end do
#ifdef _LARGE_BLOCK_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine add

subroutine triad(C, A, B, d, sz, g)
	implicit none
	integer, dimension(3)		:: sz
	integer									:: g
	integer									:: i, j, k
	integer									:: ix, jx, kx
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: C
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: A, B
	real										:: d
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
#ifdef _LARGE_BLOCK_
!$omp parallel private(i, j, k)
!$omp do
#else
!ocl nouxsimd
!ocl serial
!dir$ noparallel
#endif
	do k=1, kx
	do j=1, jx
	do i=1, ix
		C(i, j, k) = A(i, j, k) + d*B(i, j, k)
	end do
	end do
	end do
#ifdef _LARGE_BLOCK_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine triad

subroutine scal(y, a, sz, g)
	implicit none
	integer, dimension(3)		:: sz
	integer									:: g
	integer									:: i, j, k
	integer									:: ix, jx, kx
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: y
	real										:: a
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
#ifdef _LARGE_BLOCK_
!$omp parallel private(i, j, k)
!$omp do
#else
!ocl nouxsimd
!ocl serial
!dir$ noparallel
#endif
	do k=1, kx
	do j=1, jx
	do i=1, ix
		y(i, j, k) = a*y(i, j, k)
	end do
	end do
	end do
#ifdef _LARGE_BLOCK_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine scal

subroutine axpy(y, x, a, sz, g)
	implicit none
	integer, dimension(3)		:: sz
	integer									:: g
	integer									:: i, j, k
	integer									:: ix, jx, kx
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: y
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real										:: a
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
#ifdef _LARGE_BLOCK_
!$omp parallel private(i, j, k)
!$omp do
#else
!ocl nouxsimd
!ocl serial
!dir$ noparallel
#endif
	do k=1, kx
	do j=1, jx
	do i=1, ix
		y(i, j, k) = y(i, j, k) + a*x(i, j, k)
	end do
	end do
	end do
#ifdef _LARGE_BLOCK_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine axpy

subroutine xpay(y, x, a, sz, g)
	implicit none
	integer, dimension(3)		:: sz
	integer									:: g
	integer									:: i, j, k
	integer									:: ix, jx, kx
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: y
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real										:: a
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
#ifdef _LARGE_BLOCK_
!$omp parallel private(i, j, k)
!$omp do
#else
!ocl nouxsimd
!ocl serial
!dir$ noparallel
#endif
	do k=1, kx
	do j=1, jx
	do i=1, ix
		y(i, j, k) = x(i, j, k) + a*y(i, j, k)
	end do
	end do
	end do
#ifdef _LARGE_BLOCK_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine xpay

subroutine axpyz(z, x, y, a, sz, g)
	implicit none
	integer, dimension(3)		:: sz
	integer									:: g
	integer									:: i, j, k
	integer									:: ix, jx, kx
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: z
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: y
	real										:: a
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
#ifdef _LARGE_BLOCK_
!$omp parallel private(i, j, k)
!$omp do
#else
!ocl nouxsimd
!ocl serial
!dir$ noparallel
#endif
	do k=1, kx
	do j=1, jx
	do i=1, ix
		z(i, j, k) = a*x(i, j, k) + y(i, j, k)
	end do
	end do
	end do
#ifdef _LARGE_BLOCK_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine axpyz

subroutine axpbypz(z, x, y, a, b, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: z
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x, y
	real										:: a, b
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
#ifdef _LARGE_BLOCK_
!$omp parallel private(i, j, k)
!$omp do
#else
!ocl nouxsimd
!ocl serial
!dir$ noparallel
#endif
	do k=1, kx
	do j=1, jx
	do i=1, ix
		z(i, j, k) = a*x(i, j, k) + b*y(i, j, k) + z(i, j, k)
	end do
	end do
	end do
#ifdef _LARGE_BLOCK_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine axpbypz

subroutine dot(xy, x, y, sz, g)
	implicit none
	integer, dimension(3)		:: sz
	integer									:: g
	integer									:: i, j, k
	integer									:: ix, jx, kx
	real										:: xy
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: y
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	xy = 0.0
#ifdef _LARGE_BLOCK_
!$omp parallel private(i, j, k)
!$omp do reduction(+:xy), schedule(dynamic,1)
#else
!ocl nouxsimd
!ocl serial
!dir$ noparallel
#endif
	do k=1, kx
	do j=1, jx
	do i=1, ix
		xy = xy + x(i, j, k)*y(i, j, k)
	end do
	end do
	end do
#ifdef _LARGE_BLOCK_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine dot

subroutine dotx2(xy, xz, x, y, z, sz, g)
	implicit none
	integer, dimension(3)		:: sz
	integer									:: g
	integer									:: i, j, k
	integer									:: ix, jx, kx
	real										:: xy, xz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: y
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: z
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	xy = 0.0
	xz = 0.0
#ifdef _LARGE_BLOCK_
!$omp parallel private(i, j, k)
!$omp do reduction(+:xy,xz), schedule(dynamic,1)
#else
!ocl nouxsimd
!ocl serial
!dir$ noparallel
#endif
	do k=1, kx
	do j=1, jx
	do i=1, ix
		xy = xy + x(i, j, k)*y(i, j, k)
		xz = xz + x(i, j, k)*z(i, j, k)
	end do
	end do
	end do
#ifdef _LARGE_BLOCK_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine dotx2

subroutine jacobi_smoother( &
								x1, x0, &
								Ap, Aw, Ae, As, An, Ab, At, &
								b, &
								omega, &
								sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x1, x0
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: Ap, Aw, Ae, As, An, Ab, At
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: b
	real										:: omega
	real										:: r
	real										:: x_tmp
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
#ifdef _LARGE_BLOCK_
!$omp parallel private(i, j, k) &
!$omp					 private(r, x_tmp)
!$omp do schedule(dynamic,1)
#else
!ocl nouxsimd
!ocl serial
!dir$ noparallel
#endif
	do k=1, kx
!dir$ noparallel
	do j=1, jx
	do i=1, ix
		r = b(i, j, k) &
						- Aw(i, j, k)*x0(i-1, j, k) &
						- Ae(i, j, k)*x0(i+1, j, k) &
						- As(i, j, k)*x0(i, j-1, k) &
						- An(i, j, k)*x0(i, j+1, k) &
						- Ab(i, j, k)*x0(i, j, k-1) &
						- At(i, j, k)*x0(i, j, k+1) 
		x_tmp = r/Ap(i, j, k)
		x1(i, j, k) = x0(i, j, k) - omega*(x0(i, j, k) - x_tmp)
	end do
	end do
	end do
#ifdef _LARGE_BLOCK_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine jacobi_smoother

subroutine jacobi_smoother2( &
								x1, x0, &
								A, &
								b, &
								omega, &
								sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x1, x0
	real, dimension(0:6, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: A
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: b
	real										:: omega
	real										:: r
	real										:: x_tmp
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
#ifdef _LARGE_BLOCK_
!$omp parallel private(i, j, k) &
!$omp					 private(r, x_tmp)
!$omp do schedule(dynamic,1)
#else
!ocl nouxsimd
!ocl serial
!dir$ noparallel
#endif
	do k=1, kx
!dir$ noparallel
	do j=1, jx
	do i=1, ix
		r = b(i, j, k) &
						- A(1, i, j, k)*x0(i-1, j, k) &
						- A(2, i, j, k)*x0(i+1, j, k) &
						- A(3, i, j, k)*x0(i, j-1, k) &
						- A(4, i, j, k)*x0(i, j+1, k) &
						- A(5, i, j, k)*x0(i, j, k-1) &
						- A(6, i, j, k)*x0(i, j, k+1) 
		x_tmp = r/A(0, i, j, k)
		x1(i, j, k) = x0(i, j, k) - omega*(x0(i, j, k) - x_tmp)
	end do
	end do
	end do
#ifdef _LARGE_BLOCK_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine jacobi_smoother2

subroutine calc_ax( &
								Ax, &
								Ap, Aw, Ae, As, An, Ab, At, &
								x, &
								sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: Ax
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: Ap, Aw, Ae, As, An, Ab, At
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
#ifdef _LARGE_BLOCK_
!$omp parallel
!$omp do schedule(dynamic, 1) &
!$omp		,private(i, j, k) 
#else
!ocl nouxsimd
!ocl serial
!dir$ noparallel
#endif
	do k=1, kx
	do j=1, jx
	do i=1, ix
		Ax(i, j, k) = &
						+ Aw(i, j, k)*x(i-1, j, k) &
						+ Ae(i, j, k)*x(i+1, j, k) &
						+ As(i, j, k)*x(i, j-1, k) &
						+ An(i, j, k)*x(i, j+1, k) &
						+ Ab(i, j, k)*x(i, j, k-1) &
						+ At(i, j, k)*x(i, j, k+1) &
						+ Ap(i, j, k)*x(i, j, k)
	end do
	end do
	end do
#ifdef _LARGE_BLOCK_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine calc_ax

subroutine calc_r( &
								r, &
								Ap, Aw, Ae, As, An, Ab, At, &
								x, &
								b, &
								sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: r
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: Ap, Aw, Ae, As, An, Ab, At
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: b
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
#ifdef _LARGE_BLOCK_
!$omp parallel
!$omp do schedule(dynamic, 1) &
!$omp		,private(i, j, k) 
#else
!ocl nouxsimd
!ocl serial
!dir$ noparallel
#endif
	do k=1, kx
	do j=1, jx
	do i=1, ix
		r(i, j, k) = b(i, j, k) &
						- Aw(i, j, k)*x(i-1, j, k) &
						- Ae(i, j, k)*x(i+1, j, k) &
						- As(i, j, k)*x(i, j-1, k) &
						- An(i, j, k)*x(i, j+1, k) &
						- Ab(i, j, k)*x(i, j, k-1) &
						- At(i, j, k)*x(i, j, k+1) &
						- Ap(i, j, k)*x(i, j, k)
	end do
	end do
	end do
#ifdef _LARGE_BLOCK_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine calc_r

subroutine calc_r2( &
								rr, &
								Ap, Aw, Ae, As, An, Ab, At, &
								x, &
								b, &
								sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real										:: rr
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: Ap, Aw, Ae, As, An, Ab, At
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: b
	real										:: r
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	rr = 0.0
#ifdef _LARGE_BLOCK_
!$omp parallel
!$omp do schedule(dynamic, 1) &
!$omp		 private(i, j, k) &
!$omp		 private(r) &
!$omp		 reduction(+:rr)
#else
#endif
	do k=1, kx
	do j=1, jx
	do i=1, ix
		r = b(i, j, k) &
						- Aw(i, j, k)*x(i-1, j, k) &
						- Ae(i, j, k)*x(i+1, j, k) &
						- As(i, j, k)*x(i, j-1, k) &
						- An(i, j, k)*x(i, j+1, k) &
						- Ab(i, j, k)*x(i, j, k-1) &
						- At(i, j, k)*x(i, j, k+1) &
						- Ap(i, j, k)*x(i, j, k)
		rr = rr + r*r
	end do
	end do
	end do
#ifdef _LARGE_BLOCK_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine calc_r2

subroutine setup_mask( &
								m, &
								mask, &
								sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: m
	integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: mask
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
#ifdef _LARGE_BLOCK_
!$omp parallel private(i, j, k)
!$omp do
#else
#endif
	do k=1-g, kx+g
	do j=1-g, jx+g
	do i=1-g, ix+g
		m(i, j, k) = 0.0
		mask(i, j, k) = 0
	end do
	end do
	end do
#ifdef _LARGE_BLOCK_
!$omp end do
!$omp end parallel
#else
#endif

#ifdef _LARGE_BLOCK_
!$omp parallel private(i, j, k)
!$omp do
#else
#endif
	do k=1, kx
	do j=1, jx
	do i=1, ix
		m(i, j, k) = 1.0
		mask(i, j, k) = 1
	end do
	end do
	end do
#ifdef _LARGE_BLOCK_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine setup_mask

subroutine axpy_mask(y, x, a, mask, sz, g)
	implicit none
	integer									:: cx, cy, cz
	integer									:: ijk
	integer									:: ijk0, ijkx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: y
	real, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: x
	real										:: a
	integer, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: mask
  cx = sz(1) + 2*g
  cy = sz(2) + 2*g
  cz = sz(3) + 2*g
  ijk0 = cx*cy*g + cx*g + g + 1
  ijkx = cx*cy*cz - ijk0 + 1
!ocl serial
!ocl nouxsimd
!dir$ simd
	do ijk=ijk0, ijkx
		if( mask(ijk) == 1 ) then	
			y(ijk) = y(ijk) + a*x(ijk)
    end if
	end do
end subroutine axpy_mask

subroutine xpay_mask(y, x, a, mask, sz, g)
	implicit none
	integer									:: cx, cy, cz
	integer									:: ijk
	integer									:: ijk0, ijkx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: y
	real, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: x
	real										:: a
	integer, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: mask
  cx = sz(1) + 2*g
  cy = sz(2) + 2*g
  cz = sz(3) + 2*g
  ijk0 = cx*cy*g + cx*g + g + 1
  ijkx = cx*cy*cz - ijk0 + 1
!ocl serial
!ocl nouxsimd
!dir$ simd
	do ijk=ijk0, ijkx
		if( mask(ijk) == 1 ) then	
			y(ijk) = x(ijk) + a*y(ijk)
    end if
	end do
end subroutine xpay_mask

subroutine axpyz_mask(z, x, y, a, mask, sz, g)
	implicit none
	integer									:: cx, cy, cz
	integer									:: ijk
	integer									:: ijk0, ijkx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: z
	real, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: x
	real, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: y
	real										:: a
	integer, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: mask
  cx = sz(1) + 2*g
  cy = sz(2) + 2*g
  cz = sz(3) + 2*g
  ijk0 = cx*cy*g + cx*g + g + 1
  ijkx = cx*cy*cz - ijk0 + 1
!ocl serial
!ocl nouxsimd
!dir$ simd
	do ijk=ijk0, ijkx
		if( mask(ijk) == 1 ) then	
			z(ijk) = a*x(ijk) + y(ijk)
    end if
	end do
end subroutine axpyz_mask

subroutine axpbypz_mask(z, x, y, a, b, mask, sz, g)
	implicit none
	integer									:: cx, cy, cz
	integer									:: ijk
	integer									:: ijk0, ijkx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: z
	real, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: x, y
	real										:: a, b
	integer, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: mask
  cx = sz(1) + 2*g
  cy = sz(2) + 2*g
  cz = sz(3) + 2*g
  ijk0 = cx*cy*g + cx*g + g + 1
  ijkx = cx*cy*cz - ijk0 + 1
!ocl serial
!ocl nouxsimd
!dir$ simd
	do ijk=ijk0, ijkx
		if( mask(ijk) == 1 ) then	
			z(ijk) = a*x(ijk) + b*y(ijk) + z(ijk)
    end if
	end do
end subroutine axpbypz_mask

subroutine dot_mask(xy, x, y, mask, sz, g)
	implicit none
	integer									:: cx, cy, cz
	integer									:: ijk
	integer									:: ijk0, ijkx
	integer									:: g
	integer, dimension(3)		:: sz
	real										:: xy
	real, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: x
	real, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: y
	integer, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: mask
  cx = sz(1) + 2*g
  cy = sz(2) + 2*g
  cz = sz(3) + 2*g
  ijk0 = cx*cy*g + cx*g + g + 1
  ijkx = cx*cy*cz - ijk0 + 1
	xy = 0.0
!ocl serial
!ocl nouxsimd
!dir$ simd
	do ijk=ijk0, ijkx
		if( mask(ijk) == 1 ) then	
			xy = xy + x(ijk)*y(ijk)
    end if
	end do
end subroutine dot_mask

subroutine dotx2_mask(xy, xz, x, y, z, mask, sz, g)
	implicit none
	integer									:: cx, cy, cz
	integer									:: ijk
	integer									:: ijk0, ijkx
	integer									:: g
	integer, dimension(3)		:: sz
	real										:: xy, xz
	real, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: x
	real, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: y
	real, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: z
	integer, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: mask
  cx = sz(1) + 2*g
  cy = sz(2) + 2*g
  cz = sz(3) + 2*g
  ijk0 = cx*cy*g + cx*g + g + 1
  ijkx = cx*cy*cz - ijk0 + 1
	xy = 0.0
	xz = 0.0
!ocl serial
!ocl nouxsimd
!dir$ simd
	do ijk=ijk0, ijkx
		if( mask(ijk) == 1 ) then	
			xy = xy + x(ijk)*y(ijk)
			xz = xz + x(ijk)*z(ijk)
    end if
	end do
end subroutine dotx2_mask

subroutine jacobi_smoother_mask( &
								x1, x0, &
								Ap, Aw, Ae, As, An, Ab, At, &
								b, &
								mask, &
								omega, &
								sz, g)
	implicit none
	integer									:: cx, cy, cz
	integer									:: ijk
	integer									:: ijk0, ijkx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: x1, x0
	real, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: Ap, Aw, Ae, As, An, Ab, At
	real, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: b
	integer, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: mask
	real										:: omega
	real										:: r
	real										:: x_tmp
  cx = sz(1) + 2*g
  cy = sz(2) + 2*g
  cz = sz(3) + 2*g
  ijk0 = cx*cy*g + cx*g + g + 1
  ijkx = cx*cy*cz - ijk0 + 1
!ocl serial
!ocl nouxsimd
!dir$ simd
	do ijk=ijk0, ijkx
		if( mask(ijk) == 1 ) then	
			r = b(ijk) &
  						- Aw(ijk)*x0(ijk - 1) &
  						- Ae(ijk)*x0(ijk + 1) &
  						- As(ijk)*x0(ijk - cx) &
  						- An(ijk)*x0(ijk + cx) &
  						- Ab(ijk)*x0(ijk - cx*cy) &
  						- At(ijk)*x0(ijk + cx*cy)
			x_tmp = r/Ap(ijk)
			x1(ijk) = x0(ijk) - omega*(x0(ijk) - x_tmp)
    end if
	end do
end subroutine jacobi_smoother_mask

subroutine calc_ax_mask( &
								Ax, &
								Ap, Aw, Ae, As, An, Ab, At, &
								x, &
								mask, &
								sz, g)
	implicit none
	integer									:: cx, cy, cz
	integer									:: ijk
	integer									:: ijk0, ijkx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: Ax
	real, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: Ap, Aw, Ae, As, An, Ab, At
	real, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: x
	integer, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: mask
  cx = sz(1) + 2*g
  cy = sz(2) + 2*g
  cz = sz(3) + 2*g
  ijk0 = cx*cy*g + cx*g + g + 1
  ijkx = cx*cy*cz - ijk0 + 1
!ocl serial
!ocl nouxsimd
!dir$ simd
	do ijk=ijk0, ijkx
		if( mask(ijk) == 1 ) then	
			Ax(ijk) = &
							+ Aw(ijk)*x(ijk - 1) &
							+ Ae(ijk)*x(ijk + 1) &
							+ As(ijk)*x(ijk - cx) &
							+ An(ijk)*x(ijk + cx) &
							+ Ab(ijk)*x(ijk - cx*cy) &
							+ At(ijk)*x(ijk + cx*cy) &
							+ Ap(ijk)*x(ijk)
    end if
	end do
end subroutine calc_ax_mask

subroutine calc_r_mask( &
								r, &
								Ap, Aw, Ae, As, An, Ab, At, &
								x, &
								b, &
								mask, &
								sz, g)
	implicit none
	integer									:: cx, cy, cz
	integer									:: ijk
	integer									:: ijk0, ijkx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: r
	real, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: Ap, Aw, Ae, As, An, Ab, At
	real, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: x
	real, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: b
	integer, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: mask
  cx = sz(1) + 2*g
  cy = sz(2) + 2*g
  cz = sz(3) + 2*g
  ijk0 = cx*cy*g + cx*g + g + 1
  ijkx = cx*cy*cz - ijk0 + 1
!ocl serial
!ocl nouxsimd
!dir$ simd
	do ijk=ijk0, ijkx
		if( mask(ijk) == 1 ) then	
			r(ijk) = b(ijk) &
							- Aw(ijk)*x(ijk - 1) &
							- Ae(ijk)*x(ijk + 1) &
							- As(ijk)*x(ijk - cx) &
							- An(ijk)*x(ijk + cx) &
							- Ab(ijk)*x(ijk - cx*cy) &
							- At(ijk)*x(ijk + cx*cy) &
							- Ap(ijk)*x(ijk)
    end if
	end do
end subroutine calc_r_mask

subroutine calc_r2_mask( &
								rr, &
								Ap, Aw, Ae, As, An, Ab, At, &
								x, &
								b, &
								mask, &
								sz, g)
	implicit none
	integer									:: cx, cy, cz
	integer									:: ijk
	integer									:: ijk0, ijkx
	integer									:: g
	integer, dimension(3)		:: sz
	real										:: rr
	real, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: Ap, Aw, Ae, As, An, Ab, At
	real, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: x
	real, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: b
	integer, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: mask
	real										:: r
  cx = sz(1) + 2*g
  cy = sz(2) + 2*g
  cz = sz(3) + 2*g
  ijk0 = cx*cy*g + cx*g + g + 1
  ijkx = cx*cy*cz - ijk0 + 1
	rr = 0.0d0
!ocl serial
!ocl nouxsimd
!dir$ simd
	do ijk=ijk0, ijkx
		if( mask(ijk) == 1 ) then	
			r = b(ijk) &
							- Aw(ijk)*x(ijk - 1) &
							- Ae(ijk)*x(ijk + 1) &
							- As(ijk)*x(ijk - cx) &
							- An(ijk)*x(ijk + cx) &
							- Ab(ijk)*x(ijk - cx*cy) &
							- At(ijk)*x(ijk + cx*cy) &
							- Ap(ijk)*x(ijk)
			rr = rr + r*r
    end if
	end do
end subroutine calc_r2_mask

subroutine axpy_mask2(y, x, a, m, sz, g)
	implicit none
	integer									:: cx, cy, cz
	integer									:: ijk
	integer									:: ijk0, ijkx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: y
	real, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: x
	real										:: a
	real, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: m
  cx = sz(1) + 2*g
  cy = sz(2) + 2*g
  cz = sz(3) + 2*g
  ijk0 = cx*cy*g + cx*g + g + 1
  ijkx = cx*cy*cz - ijk0 + 1
!ocl serial
!ocl nouxsimd
	do ijk=ijk0, ijkx
		y(ijk) = y(ijk) + a*x(ijk)*m(ijk)
	end do
end subroutine axpy_mask2

subroutine xpay_mask2(y, x, a, m, sz, g)
	implicit none
	integer									:: cx, cy, cz
	integer									:: ijk
	integer									:: ijk0, ijkx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: y
	real, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: x
	real										:: a
	real, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: m
  cx = sz(1) + 2*g
  cy = sz(2) + 2*g
  cz = sz(3) + 2*g
  ijk0 = cx*cy*g + cx*g + g + 1
  ijkx = cx*cy*cz - ijk0 + 1
!ocl serial
!ocl nouxsimd
	do ijk=ijk0, ijkx
		y(ijk) = (x(ijk) + a*y(ijk))*m(ijk) + y(ijk)*(1.0 - m(ijk))
	end do
end subroutine xpay_mask2

subroutine axpyz_mask2(z, x, y, a, m, sz, g)
	implicit none
	integer									:: cx, cy, cz
	integer									:: ijk
	integer									:: ijk0, ijkx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: z
	real, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: x
	real, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: y
	real										:: a
	real, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: m
  cx = sz(1) + 2*g
  cy = sz(2) + 2*g
  cz = sz(3) + 2*g
  ijk0 = cx*cy*g + cx*g + g + 1
  ijkx = cx*cy*cz - ijk0 + 1
!ocl serial
!ocl nouxsimd
	do ijk=ijk0, ijkx
		z(ijk) = (a*x(ijk) + y(ijk))*m(ijk) + z(ijk)*(1.0 - m(ijk))
	end do
end subroutine axpyz_mask2

subroutine axpbypz_mask2(z, x, y, a, b, m, sz, g)
	implicit none
	integer									:: cx, cy, cz
	integer									:: ijk
	integer									:: ijk0, ijkx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: z
	real, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: x, y
	real										:: a, b
	real, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: m
  cx = sz(1) + 2*g
  cy = sz(2) + 2*g
  cz = sz(3) + 2*g
  ijk0 = cx*cy*g + cx*g + g + 1
  ijkx = cx*cy*cz - ijk0 + 1
!ocl serial
!ocl nouxsimd
	do ijk=ijk0, ijkx
		z(ijk) = (a*x(ijk) + b*y(ijk))*m(ijk) + z(ijk)
	end do
end subroutine axpbypz_mask2

subroutine dot_mask2(xy, x, y, m, sz, g)
	implicit none
	integer									:: cx, cy, cz
	integer									:: ijk
	integer									:: ijk0, ijkx
	integer									:: g
	integer, dimension(3)		:: sz
	real										:: xy
	real, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: x
	real, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: y
	real, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: m
  cx = sz(1) + 2*g
  cy = sz(2) + 2*g
  cz = sz(3) + 2*g
  ijk0 = cx*cy*g + cx*g + g + 1
  ijkx = cx*cy*cz - ijk0 + 1
	xy = 0.0
!ocl serial
!ocl nouxsimd
	do ijk=ijk0, ijkx
		xy = xy + x(ijk)*y(ijk)*m(ijk)
	end do
end subroutine dot_mask2

subroutine dotx2_mask2(xy, xz, x, y, z, m, sz, g)
	implicit none
	integer									:: cx, cy, cz
	integer									:: ijk
	integer									:: ijk0, ijkx
	integer									:: g
	integer, dimension(3)		:: sz
	real										:: xy, xz
	real, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: x
	real, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: y
	real, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: z
	real, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: m
  cx = sz(1) + 2*g
  cy = sz(2) + 2*g
  cz = sz(3) + 2*g
  ijk0 = cx*cy*g + cx*g + g + 1
  ijkx = cx*cy*cz - ijk0 + 1
	xy = 0.0
	xz = 0.0
!ocl serial
!ocl nouxsimd
	do ijk=ijk0, ijkx
		xy = xy + x(ijk)*y(ijk)*m(ijk)
		xz = xz + x(ijk)*z(ijk)*m(ijk)
	end do
end subroutine dotx2_mask2

subroutine jacobi_smoother_mask2( &
								x1, x0, &
								Ap, Aw, Ae, As, An, Ab, At, &
								b, &
								m, &
								omega, &
								sz, g)
	implicit none
	integer									:: cx, cy, cz
	integer									:: ijk
	integer									:: ijk0, ijkx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: x1, x0
	real, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: Ap, Aw, Ae, As, An, Ab, At
	real, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: b
	real, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: m
	real										:: omega
	real										:: r
	real										:: x_tmp
  cx = sz(1) + 2*g
  cy = sz(2) + 2*g
  cz = sz(3) + 2*g
  ijk0 = cx*cy*g + cx*g + g + 1
  ijkx = cx*cy*cz - ijk0 + 1
!ocl serial
!ocl nouxsimd
	do ijk=ijk0, ijkx
		r = b(ijk) &
 						- Aw(ijk)*x0(ijk - 1) &
 						- Ae(ijk)*x0(ijk + 1) &
 						- As(ijk)*x0(ijk - cx) &
 						- An(ijk)*x0(ijk + cx) &
 						- Ab(ijk)*x0(ijk - cx*cy) &
 						- At(ijk)*x0(ijk + cx*cy)
		x_tmp = r/Ap(ijk)
		x1(ijk) = x0(ijk) - omega*(x0(ijk) - x_tmp)*m(ijk)
	end do
end subroutine jacobi_smoother_mask2

subroutine calc_ax_mask2( &
								Ax, &
								Ap, Aw, Ae, As, An, Ab, At, &
								x, &
								m, &
								sz, g)
	implicit none
	integer									:: cx, cy, cz
	integer									:: ijk
	integer									:: ijk0, ijkx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: Ax
	real, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: Ap, Aw, Ae, As, An, Ab, At
	real, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: x
	real, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: m
  cx = sz(1) + 2*g
  cy = sz(2) + 2*g
  cz = sz(3) + 2*g
  ijk0 = cx*cy*g + cx*g + g + 1
  ijkx = cx*cy*cz - ijk0 + 1
!ocl serial
!ocl nouxsimd
	do ijk=ijk0, ijkx
		Ax(ijk) = &
						( Aw(ijk)*x(ijk - 1) &
						+ Ae(ijk)*x(ijk + 1) &
						+ As(ijk)*x(ijk - cx) &
						+ An(ijk)*x(ijk + cx) &
						+ Ab(ijk)*x(ijk - cx*cy) &
						+ At(ijk)*x(ijk + cx*cy) &
						+ Ap(ijk)*x(ijk) )*m(ijk)
	end do
end subroutine calc_ax_mask2

subroutine calc_r_mask2( &
								r, &
								Ap, Aw, Ae, As, An, Ab, At, &
								x, &
								b, &
								m, &
								sz, g)
	implicit none
	integer									:: cx, cy, cz
	integer									:: ijk
	integer									:: ijk0, ijkx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: r
	real, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: Ap, Aw, Ae, As, An, Ab, At
	real, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: x
	real, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: b
	real, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: m
  cx = sz(1) + 2*g
  cy = sz(2) + 2*g
  cz = sz(3) + 2*g
  ijk0 = cx*cy*g + cx*g + g + 1
  ijkx = cx*cy*cz - ijk0 + 1
!ocl serial
!ocl nouxsimd
	do ijk=ijk0, ijkx
		r(ijk) = (b(ijk) &
						- Aw(ijk)*x(ijk - 1) &
						- Ae(ijk)*x(ijk + 1) &
						- As(ijk)*x(ijk - cx) &
						- An(ijk)*x(ijk + cx) &
						- Ab(ijk)*x(ijk - cx*cy) &
						- At(ijk)*x(ijk + cx*cy) &
						- Ap(ijk)*x(ijk))*m(ijk)
	end do
end subroutine calc_r_mask2

subroutine calc_r2_mask2( &
								rr, &
								Ap, Aw, Ae, As, An, Ab, At, &
								x, &
								b, &
								m, &
								sz, g)
	implicit none
	integer									:: cx, cy, cz
	integer									:: ijk
	integer									:: ijk0, ijkx
	integer									:: g
	integer, dimension(3)		:: sz
	real										:: rr
	real, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: Ap, Aw, Ae, As, An, Ab, At
	real, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: x
	real, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: b
	real, dimension(1:(sz(1)+2*g)*(sz(2)+2*g)*(sz(3)+2*g))	:: m
	real										:: r
  cx = sz(1) + 2*g
  cy = sz(2) + 2*g
  cz = sz(3) + 2*g
  ijk0 = cx*cy*g + cx*g + g + 1
  ijkx = cx*cy*cz - ijk0 + 1
	rr = 0.0d0
!ocl serial
!ocl nouxsimd
	do ijk=ijk0, ijkx
		r = b(ijk) &
						- Aw(ijk)*x(ijk - 1) &
						- Ae(ijk)*x(ijk + 1) &
						- As(ijk)*x(ijk - cx) &
						- An(ijk)*x(ijk + cx) &
						- Ab(ijk)*x(ijk - cx*cy) &
						- At(ijk)*x(ijk + cx*cy) &
						- Ap(ijk)*x(ijk)
		rr = rr + r*r*m(ijk)
	end do
end subroutine calc_r2_mask2


