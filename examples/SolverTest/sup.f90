!
! BCMTools
!
! Copyright (C) 2011-2014 Institute of Industrial Science, The University of Tokyo.
! All rights reserved.
!
! Copyright (c) 2012-2014 Advanced Institute for Computational Science, RIKEN.
! All rights reserved.
!

subroutine sup_get_intp_index(ic, r, i, n)
	implicit none
	integer									:: ic
	real										:: r
	integer									:: i
	integer									:: n
	ic = i/2
	r = 0.25 + 0.5*mod(i, 2)
	if( i < 2 ) then
		ic = 1
		r = -0.25
	else if( i > 2*n-1 ) then
		ic = n-1
		r = 1.25
	end if
end subroutine sup_get_intp_index

subroutine sup_copy_from_neighbor( &
							data_dst, &
							i1_dst, &
							data_src, &
							i1_src, &
							sz_c, &
							sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	integer, dimension(3)		:: sz_c
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: data_dst
	integer, dimension(3)		:: i1_dst
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: data_src
	integer, dimension(3)		:: i1_src
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k)
!$omp do
#else
#endif
	do k=0, sz_c(3)-1
	do j=0, sz_c(2)-1
	do i=0, sz_c(1)-1
		data_dst(i + i1_dst(1), j + i1_dst(2), k + i1_dst(3)) = data_src(i + i1_src(1), j + i1_src(2), k + i1_src(3)) 
	end do
	end do
	end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine sup_copy_from_neighbor

subroutine sup_copy_from_neighbor_c2f( &
							data_dst, &
							i1_dst, &
							data_src, &
							i1_src, &
							sz_c, &
							sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	integer, dimension(3)		:: sz_c
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: data_dst
	integer, dimension(3)		:: i1_dst
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: data_src
	integer, dimension(3)		:: i1_src
	real										:: data_src_nnn, data_src_pnn, data_src_npn, data_src_nnp
	real										:: data_src_npp, data_src_pnp, data_src_ppn, data_src_ppp
	integer									:: i0, j0, k0
	integer									:: i1, j1, k1
	real										:: r, s, t
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) &
!$omp          private(data_src_nnn, data_src_pnn, data_src_npn, data_src_nnp) &
!$omp          private(data_src_npp, data_src_pnp, data_src_ppn, data_src_ppp) &
!$omp					 private(i0, j0, k0) &
!$omp					 private(i1, j1, k1) &
!$omp					 private(r, s, t) 
!$omp do
#else
#endif
	do k=0, sz_c(3)-1
	do j=0, sz_c(2)-1
	do i=0, sz_c(1)-1
		i0 = i + i1_src(1)
		j0 = j + i1_src(2)
		k0 = k + i1_src(3)
		call sup_get_intp_index(i1, r, i0, sz(1))
		call sup_get_intp_index(j1, s, j0, sz(2))
		call sup_get_intp_index(k1, t, k0, sz(3))
!		i1 = i0/2
!		j1 = j0/2
!		k1 = k0/2
!		r = 0.25 + 0.5*mod(i0, 2)
!		s = 0.25 + 0.5*mod(j0, 2)
!		t = 0.25 + 0.5*mod(k0, 2)
		data_src_nnn = data_src(i1  , j1  , k1  )
		data_src_pnn = data_src(i1+1, j1  , k1  )
		data_src_npn = data_src(i1  , j1+1, k1  )
		data_src_nnp = data_src(i1  , j1  , k1+1)
		data_src_npp = data_src(i1  , j1+1, k1+1)
		data_src_pnp = data_src(i1+1, j1  , k1+1)
		data_src_ppn = data_src(i1+1, j1+1, k1  )
		data_src_ppp = data_src(i1+1, j1+1, k1+1)
		data_dst(i + i1_dst(1), j + i1_dst(2), k + i1_dst(3)) = (1.0 - t)*( &
																																(1.0 - s)*( (1.0 - r)*data_src_nnn + r*data_src_pnn ) &
																															+         s*( (1.0 - r)*data_src_npn + r*data_src_ppn ) ) &
																														+       t*( &
																																(1.0 - s)*( (1.0 - r)*data_src_nnp + r*data_src_pnp ) &
																															+         s*( (1.0 - r)*data_src_npp + r*data_src_ppp ) ) 
	end do
	end do
	end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine sup_copy_from_neighbor_c2f

subroutine sup_copy_from_neighbor_f2c( &
							data_dst, &
							i1_dst, &
							data_src, &
							i1_src, &
							sz_c, &
							sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	integer, dimension(3)		:: sz_c
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: data_dst
	integer, dimension(3)		:: i1_dst
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: data_src
	integer, dimension(3)		:: i1_src
	real										:: data_src_nnn, data_src_pnn, data_src_npn, data_src_nnp
	real										:: data_src_npp, data_src_pnp, data_src_ppn, data_src_ppp
	integer									:: i0, j0, k0
	integer									:: i1, j1, k1
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) &
!$omp          private(data_src_nnn, data_src_pnn, data_src_npn, data_src_nnp) &
!$omp          private(data_src_npp, data_src_pnp, data_src_ppn, data_src_ppp) &
!$omp					 private(i0, j0, k0) &
!$omp					 private(i1, j1, k1) 
!$omp do
#else
#endif
	do k=0, sz_c(3)-1
	do j=0, sz_c(2)-1
	do i=0, sz_c(1)-1
		i0 = i + i1_src(1)
		j0 = j + i1_src(2)
		k0 = k + i1_src(3)
		i1 = 2*i0 - 1
		j1 = 2*j0 - 1
		k1 = 2*k0 - 1
		data_src_nnn = data_src(i1  , j1  , k1  )
		data_src_pnn = data_src(i1+1, j1  , k1  )
		data_src_npn = data_src(i1  , j1+1, k1  )
		data_src_nnp = data_src(i1  , j1  , k1+1)
		data_src_npp = data_src(i1  , j1+1, k1+1)
		data_src_pnp = data_src(i1+1, j1  , k1+1)
		data_src_ppn = data_src(i1+1, j1+1, k1  )
		data_src_ppp = data_src(i1+1, j1+1, k1+1)
		data_dst(i + i1_dst(1), j + i1_dst(2), k + i1_dst(3)) = 0.125*( data_src_nnn + data_src_pnn + data_src_npn + data_src_nnp &
																																	+ data_src_npp + data_src_pnp + data_src_ppn + data_src_ppp )
	end do
	end do
	end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine sup_copy_from_neighbor_f2c

subroutine sup_copy_to_buffer( &
							buffer, &
							data_src, &
							i1_src, &
							sz_c, &
							sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	integer, dimension(3)		:: sz_c
!	real, dimension(0:sz_c(1)*sz_c(2)*sz_c(3)-1)						:: buffer
	real, dimension(0:sz_c(1)-1, 0:sz_c(2)-1, 0:sz_c(3)-1)	:: buffer
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: data_src
	integer, dimension(3)		:: i1_src
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) 
!$omp do
#else
#endif
	do k=0, sz_c(3)-1
	do j=0, sz_c(2)-1
	do i=0, sz_c(1)-1
!		buffer(i + sz_c(1)*(j + sz_c(2)*k)) = data_src(i + i1_src(1), j + i1_src(2), k + i1_src(3)) 
		buffer(i, j, k) = data_src(i + i1_src(1), j + i1_src(2), k + i1_src(3)) 
	end do
	end do
	end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine sup_copy_to_buffer

subroutine sup_copy_to_buffer_c2f( &
							buffer, &
							data_src, &
							i1_src, &
							sz_c, &
							sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	integer, dimension(3)		:: sz_c
	real, dimension(0:sz_c(1)-1, 0:sz_c(2)-1, 0:sz_c(3)-1)	:: buffer
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: data_src
	integer, dimension(3)		:: i1_src
	real										:: data_src_nnn, data_src_pnn, data_src_npn, data_src_nnp
	real										:: data_src_npp, data_src_pnp, data_src_ppn, data_src_ppp
	integer									:: i0, j0, k0
	integer									:: i1, j1, k1
	real										:: r, s, t
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) &
!$omp          private(data_src_nnn, data_src_pnn, data_src_npn, data_src_nnp) &
!$omp          private(data_src_npp, data_src_pnp, data_src_ppn, data_src_ppp) &
!$omp					 private(i0, j0, k0) &
!$omp					 private(i1, j1, k1) &
!$omp					 private(r, s, t) 
!$omp do
#else
#endif
	do k=0, sz_c(3)-1
	do j=0, sz_c(2)-1
	do i=0, sz_c(1)-1
		i0 = i + i1_src(1)
		j0 = j + i1_src(2)
		k0 = k + i1_src(3)
		call sup_get_intp_index(i1, r, i0, sz(1))
		call sup_get_intp_index(j1, s, j0, sz(2))
		call sup_get_intp_index(k1, t, k0, sz(3))
!		i1 = i0/2
!		j1 = j0/2
!		k1 = k0/2
!		r = 0.25 + 0.5*mod(i0, 2)
!		s = 0.25 + 0.5*mod(j0, 2)
!		t = 0.25 + 0.5*mod(k0, 2)
		data_src_nnn = data_src(i1  , j1  , k1  )
		data_src_pnn = data_src(i1+1, j1  , k1  )
		data_src_npn = data_src(i1  , j1+1, k1  )
		data_src_nnp = data_src(i1  , j1  , k1+1)
		data_src_npp = data_src(i1  , j1+1, k1+1)
		data_src_pnp = data_src(i1+1, j1  , k1+1)
		data_src_ppn = data_src(i1+1, j1+1, k1  )
		data_src_ppp = data_src(i1+1, j1+1, k1+1)
		buffer(i, j, k) = (1.0 - t)*( &
													(1.0 - s)*( (1.0 - r)*data_src_nnn + r*data_src_pnn ) &
												+         s*( (1.0 - r)*data_src_npn + r*data_src_ppn ) ) &
													+       t*( &
													(1.0 - s)*( (1.0 - r)*data_src_nnp + r*data_src_pnp ) &
												+         s*( (1.0 - r)*data_src_npp + r*data_src_ppp ) ) 
	end do
	end do
	end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine sup_copy_to_buffer_c2f

subroutine sup_copy_to_buffer_f2c( &
							buffer, &
							data_src, &
							i1_src, &
							sz_c, &
							sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	integer, dimension(3)		:: sz_c
	real, dimension(0:sz_c(1)-1, 0:sz_c(2)-1, 0:sz_c(3)-1)	:: buffer
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: data_src
	integer, dimension(3)		:: i1_src
	real										:: data_src_nnn, data_src_pnn, data_src_npn, data_src_nnp
	real										:: data_src_npp, data_src_pnp, data_src_ppn, data_src_ppp
	integer									:: i0, j0, k0
	integer									:: i1, j1, k1
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) &
!$omp          private(data_src_nnn, data_src_pnn, data_src_npn, data_src_nnp) &
!$omp          private(data_src_npp, data_src_pnp, data_src_ppn, data_src_ppp) &
!$omp					 private(i0, j0, k0) &
!$omp					 private(i1, j1, k1) 
!$omp do
#else
#endif
	do k=0, sz_c(3)-1
	do j=0, sz_c(2)-1
	do i=0, sz_c(1)-1
		i0 = i + i1_src(1)
		j0 = j + i1_src(2)
		k0 = k + i1_src(3)
		i1 = 2*i0 - 1
		j1 = 2*j0 - 1
		k1 = 2*k0 - 1
		data_src_nnn = data_src(i1  , j1  , k1  )
		data_src_pnn = data_src(i1+1, j1  , k1  )
		data_src_npn = data_src(i1  , j1+1, k1  )
		data_src_nnp = data_src(i1  , j1  , k1+1)
		data_src_npp = data_src(i1  , j1+1, k1+1)
		data_src_pnp = data_src(i1+1, j1  , k1+1)
		data_src_ppn = data_src(i1+1, j1+1, k1  )
		data_src_ppp = data_src(i1+1, j1+1, k1+1)
		buffer(i, j, k) = 0.125*( data_src_nnn + data_src_pnn + data_src_npn + data_src_nnp &
														+ data_src_npp + data_src_pnp + data_src_ppn + data_src_ppp )
	end do
	end do
	end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine sup_copy_to_buffer_f2c

subroutine sup_copy_from_buffer( &
							data_dst, &
							i1_dst, &
							buffer, &
							sz_c, &
							sz, g)
	implicit none
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	integer, dimension(3)		:: sz_c
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: data_dst
	integer, dimension(3)		:: i1_dst
	real, dimension(0:sz_c(1)-1, 0:sz_c(2)-1, 0:sz_c(3)-1)	:: buffer
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) 
!$omp do
#else
#endif
	do k=0, sz_c(3)-1
	do j=0, sz_c(2)-1
	do i=0, sz_c(1)-1
		data_dst(i + i1_dst(1), j + i1_dst(2), k + i1_dst(3)) = buffer(i, j, k) 
	end do
	end do
	end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine sup_copy_from_buffer

