subroutine sf3d_copy_x2(x, xc, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real, dimension(1-g:sz(1)/2+g, 1-g:sz(2)/2+g, 1-g:sz(3)/2+g)	:: xc
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
#ifdef _LARGE_BLOCK_
!$omp parallel private(i, j, k) 
!$omp do 
#else
#endif
	do k=1, kx
	do j=1, jx
	do i=1, ix
		x(i, j, k) = xc((i+1)/2, (j+1)/2, (k+1)/2)
	end do
	end do
	end do
#ifdef _LARGE_BLOCK_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine sf3d_copy_x2

subroutine sf3d_calc_stats(dsum, dmax, dmin, dabsmax, dabsmin, ddata, sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: ddata
	real										:: ddata0
	real										:: dsum, dmax, dmin, dabsmax, dabsmin
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
#ifdef _LARGE_BLOCK_
!$omp parallel private(i, j, k) &
!$omp         ,private(ddata0)
!$omp do reduction(+: dsum) &
!$omp		,reduction(max: dmax, dabsmax) &
!$omp		,reduction(min: dmin, dabsmin)
#else
#endif
	do k=1, kx
	do j=1, jx
	do i=1, ix
		ddata0 = ddata(i, j, k)

		dsum = dsum + ddata0

		if( ddata0 > dmax ) then
			dmax = ddata0
		endif

		if( ddata0 < dmin ) then
			dmin = ddata0
		endif

		if( abs(ddata0) > dabsmax ) then
			dabsmax = abs(ddata0)
		endif

		if( abs(ddata0) < dabsmin ) then
			dabsmin = abs(ddata0)
		endif
	end do
	end do
	end do
#ifdef _LARGE_BLOCK_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine sf3d_calc_stats

