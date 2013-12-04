!
! BCMTools
!
! Copyright (C) 2011-2013 Institute of Industrial Science, The University of Tokyo.
! All rights reserved.
!
! Copyright (c) 2012-2013 Advanced Institute for Computational Science, RIKEN.
! All rights reserved.
!

subroutine allreduce(total, a)
	implicit none
	include 'mpif.h'
	real										:: total, a
	integer									:: ierror
	total = 0.0
#ifdef _REAL_IS_DOUBLE_
	call mpi_allreduce(a, total, 1, mpi_double_precision, mpi_sum, mpi_comm_world, ierror)
#else
	call mpi_allreduce(a, total, 1, mpi_real, mpi_sum, mpi_comm_world, ierror)
#endif
end subroutine allreduce

subroutine allreduce_max(total, a)
	implicit none
	include 'mpif.h'
	real										:: total, a
	integer									:: ierror
	total = 0.0
#ifdef _REAL_IS_DOUBLE_
	call mpi_allreduce(a, total, 1, mpi_double_precision, mpi_max, mpi_comm_world, ierror)
#else
	call mpi_allreduce(a, total, 1, mpi_real, mpi_max, mpi_comm_world, ierror)
#endif
end subroutine allreduce_max

subroutine allreduce_min(total, a)
	implicit none
	include 'mpif.h'
	real										:: total, a
	integer									:: ierror
	total = 0.0
#ifdef _REAL_IS_DOUBLE_
	call mpi_allreduce(a, total, 1, mpi_double_precision, mpi_min, mpi_comm_world, ierror)
#else
	call mpi_allreduce(a, total, 1, mpi_real, mpi_min, mpi_comm_world, ierror)
#endif
end subroutine allreduce_min

subroutine dot_all(xy, x, y, sz, g)
	implicit none
	include 'mpif.h'
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real										:: xy
	integer									:: ierror
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: y
	real										:: xy_
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
	xy_ = 0.0
!$omp parallel private(i, j, k)
!$omp do reduction(+:xy_)
	do k=1, kx
	do j=1, jx
	do i=1, ix
		xy_ = xy_ + x(i, j, k)*y(i, j, k)
	end do
	end do
	end do
!$omp end do
!$omp end parallel
	xy = 0.0
#ifdef _REAL_IS_DOUBLE_
	call mpi_allreduce(xy_, xy, 1, mpi_double_precision, mpi_sum, mpi_comm_world, ierror)
#else
	call mpi_allreduce(xy_, xy, 1, mpi_real, mpi_sum, mpi_comm_world, ierror)
#endif
end subroutine dot_all

subroutine comm_band_cells(x, &
													sz, g, &
													mx, &
													sx1, sx3, sx2, sx4, sx5, sx6, &
													rx1, rx3, rx2, rx4, rx5, rx6, &
													node)
	implicit none
	include 'mpif.h'
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	integer									:: m
	integer									:: mx
	real, dimension(1:mx, 1-g:sz(2)+g, 1-g:sz(3)+g) :: sx1, sx3, rx1, rx3
	real, dimension(1:mx, 1-g:sz(1)+g, 1-g:sz(3)+g) :: sx2, sx4, rx2, rx4
	real, dimension(1:mx, 1-g:sz(1)+g, 1-g:sz(2)+g) :: sx5, sx6, rx5, rx6
	integer, dimension(0:6)	:: node
	integer									:: ierror
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
!$omp parallel private(i, j, k)
!$omp do
	do k=1-g, kx+g
	do j=1-g, jx+g
	do m=1, mx
		sx1(m, j, k) = x(ix-mx+m, j, k)
		sx3(m, j, k) = x(m      , j, k)
	end do
	end do
	end do
!$omp end do
!$omp do
	do k=1-g, kx+g
	do i=1-g, ix+g
	do m=1, mx
		sx2(m, i, k) = x(i, jx-mx+m, k)
		sx4(m, i, k) = x(i, m      , k)
	end do
	end do
	end do
!$omp end do
!$omp do
	do j=1-g, jx+g
	do i=1-g, ix+g
	do m=1, mx
		sx5(m, i, j) = x(i, j, kx-mx+m)
		sx6(m, i, j) = x(i, j, m      )
	end do
	end do
	end do
!$omp end do
!$omp end parallel
#ifdef _REAL_IS_DOUBLE_
		call mpi_sendrecv(sx1, (jx+2*g)*(kx+2*g)*mx, mpi_double_precision, node(1), 1, &
											rx3, (jx+2*g)*(kx+2*g)*mx, mpi_double_precision, node(3), 1, &
											mpi_comm_world, mpi_statuses_ignore, ierror)

		call mpi_sendrecv(sx3, (jx+2*g)*(kx+2*g)*mx, mpi_double_precision, node(3), 3, &
											rx1, (jx+2*g)*(kx+2*g)*mx, mpi_double_precision, node(1), 3, &
											mpi_comm_world, mpi_statuses_ignore, ierror)

		call mpi_sendrecv(sx2, (ix+2*g)*(kx+2*g)*mx, mpi_double_precision, node(2), 2, &
											rx4, (ix+2*g)*(kx+2*g)*mx, mpi_double_precision, node(4), 2, &
											mpi_comm_world, mpi_statuses_ignore, ierror)

		call mpi_sendrecv(sx4, (ix+2*g)*(kx+2*g)*mx, mpi_double_precision, node(4), 4, &
											rx2, (ix+2*g)*(kx+2*g)*mx, mpi_double_precision, node(2), 4, &
											mpi_comm_world, mpi_statuses_ignore, ierror)

		call mpi_sendrecv(sx5, (ix+2*g)*(jx+2*g)*mx, mpi_double_precision, node(5), 5, &
											rx6, (ix+2*g)*(jx+2*g)*mx, mpi_double_precision, node(6), 5, &
											mpi_comm_world, mpi_statuses_ignore, ierror)

		call mpi_sendrecv(sx6, (ix+2*g)*(jx+2*g)*mx, mpi_double_precision, node(6), 6, &
											rx5, (ix+2*g)*(jx+2*g)*mx, mpi_double_precision, node(5), 6, &
											mpi_comm_world, mpi_statuses_ignore, ierror)
#else
		call mpi_sendrecv(sx1, (jx+2*g)*(kx+2*g)*mx, mpi_real, node(1), 1, &
											rx3, (jx+2*g)*(kx+2*g)*mx, mpi_real, node(3), 1, &
											mpi_comm_world, mpi_statuses_ignore, ierror)

		call mpi_sendrecv(sx3, (jx+2*g)*(kx+2*g)*mx, mpi_real, node(3), 3, &
											rx1, (jx+2*g)*(kx+2*g)*mx, mpi_real, node(1), 3, &
											mpi_comm_world, mpi_statuses_ignore, ierror)

		call mpi_sendrecv(sx2, (ix+2*g)*(kx+2*g)*mx, mpi_real, node(2), 2, &
											rx4, (ix+2*g)*(kx+2*g)*mx, mpi_real, node(4), 2, &
											mpi_comm_world, mpi_statuses_ignore, ierror)

		call mpi_sendrecv(sx4, (ix+2*g)*(kx+2*g)*mx, mpi_real, node(4), 4, &
											rx2, (ix+2*g)*(kx+2*g)*mx, mpi_real, node(2), 4, &
											mpi_comm_world, mpi_statuses_ignore, ierror)

		call mpi_sendrecv(sx5, (ix+2*g)*(jx+2*g)*mx, mpi_real, node(5), 5, &
											rx6, (ix+2*g)*(jx+2*g)*mx, mpi_real, node(6), 5, &
											mpi_comm_world, mpi_statuses_ignore, ierror)

		call mpi_sendrecv(sx6, (ix+2*g)*(jx+2*g)*mx, mpi_real, node(6), 6, &
											rx5, (ix+2*g)*(jx+2*g)*mx, mpi_real, node(5), 6, &
											mpi_comm_world, mpi_statuses_ignore, ierror)
#endif
!$omp parallel private(i, j, k)
!$omp do
	do k=1-g, kx+g
	do j=1-g, jx+g
	do m=1, mx
		x( ix+m, j, k) = rx1(m, j, k)
		x(-mx+m, j, k) = rx3(m, j, k)
	end do
	end do
	end do
!$omp end do
!$omp do
	do k=1-g, kx+g
	do i=1-g, ix+g
	do m=1, mx
		x(i,  jx+m, k) = rx2(m, i, k)
		x(i, -mx+m, k) = rx4(m, i, k)
	end do
	end do
	end do
!$omp end do
!$omp do
	do j=1-g, jx+g
	do i=1-g, ix+g
	do m=1, mx
		x(i, j,  kx+m) = rx5(m, i, j)
		x(i, j, -mx+m) = rx6(m, i, j)
	end do
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine comm_band_cells

subroutine comm_band_cells_nb(x, &
													sz, g, &
													mx, &
													sx1, sx3, sx2, sx4, sx5, sx6, &
													rx1, rx3, rx2, rx4, rx5, rx6, &
													node)
	implicit none
	include 'mpif.h'
	integer									:: i, j, k
	integer									:: ix, jx, kx
	integer									:: g
	integer, dimension(3)		:: sz
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)	:: x
	integer									:: m
	integer									:: mx
	real, dimension(1:mx, 1-g:sz(2)+g, 1-g:sz(3)+g) :: sx1, sx3, rx1, rx3
	real, dimension(1:mx, 1-g:sz(1)+g, 1-g:sz(3)+g) :: sx2, sx4, rx2, rx4
	real, dimension(1:mx, 1-g:sz(1)+g, 1-g:sz(2)+g) :: sx5, sx6, rx5, rx6
	integer, dimension(0:6)	:: node
	integer									:: ierror
	integer, dimension(4)		:: req
	ix = sz(1)
	jx = sz(2)
	kx = sz(3)
!$omp parallel private(i, j, k)
!$omp do
	do k=1-g, kx+g
	do j=1-g, jx+g
	do m=1, mx
		sx1(m, j, k) = x(ix-mx+m, j, k)
		sx3(m, j, k) = x(m      , j, k)
	end do
	end do
	end do
!$omp end do
!$omp do
	do k=1-g, kx+g
	do i=1-g, ix+g
	do m=1, mx
		sx2(m, i, k) = x(i, jx-mx+m, k)
		sx4(m, i, k) = x(i, m      , k)
	end do
	end do
	end do
!$omp end do
!$omp do
	do j=1-g, jx+g
	do i=1-g, ix+g
	do m=1, mx
		sx5(m, i, j) = x(i, j, kx-mx+m)
		sx6(m, i, j) = x(i, j, m      )
	end do
	end do
	end do
!$omp end do
!$omp end parallel
#ifdef _REAL_IS_DOUBLE_
		call mpi_irecv(rx1, (jx+2*g)*(kx+2*g)*mx, mpi_double_precision, node(1), 3, &
									 mpi_comm_world, req(4), ierror)
		call mpi_irecv(rx3, (jx+2*g)*(kx+2*g)*mx, mpi_double_precision, node(3), 1, &
									 mpi_comm_world, req(2), ierror)
		call mpi_isend(sx1, (jx+2*g)*(kx+2*g)*mx, mpi_double_precision, node(1), 1, &
									 mpi_comm_world, req(1), ierror)
		call mpi_isend(sx3, (jx+2*g)*(kx+2*g)*mx, mpi_double_precision, node(3), 3, &
									 mpi_comm_world, req(3), ierror)
		call mpi_waitall(4, req, mpi_statuses_ignore, ierror)

		call mpi_irecv(rx2, (ix+2*g)*(kx+2*g)*mx, mpi_double_precision, node(2), 4, &
									 mpi_comm_world, req(4), ierror)
		call mpi_irecv(rx4, (ix+2*g)*(kx+2*g)*mx, mpi_double_precision, node(4), 2, &
									 mpi_comm_world, req(2), ierror)
		call mpi_isend(sx2, (ix+2*g)*(kx+2*g)*mx, mpi_double_precision, node(2), 2, &
									 mpi_comm_world, req(1), ierror)
		call mpi_isend(sx4, (ix+2*g)*(kx+2*g)*mx, mpi_double_precision, node(4), 4, &
									 mpi_comm_world, req(3), ierror)
		call mpi_waitall(4, req, mpi_statuses_ignore, ierror)

		call mpi_irecv(rx5, (ix+2*g)*(jx+2*g)*mx, mpi_double_precision, node(5), 6, &
									 mpi_comm_world, req(4), ierror)
		call mpi_irecv(rx6, (ix+2*g)*(jx+2*g)*mx, mpi_double_precision, node(6), 5, &
									 mpi_comm_world, req(2), ierror)
		call mpi_isend(sx5, (ix+2*g)*(jx+2*g)*mx, mpi_double_precision, node(5), 5, &
									 mpi_comm_world, req(1), ierror)
		call mpi_isend(sx6, (ix+2*g)*(jx+2*g)*mx, mpi_double_precision, node(6), 6, &
									 mpi_comm_world, req(3), ierror)
		call mpi_waitall(4, req, mpi_statuses_ignore, ierror)
#else
		call mpi_irecv(rx1, (jx+2*g)*(kx+2*g)*mx, mpi_real, node(1), 3, &
									 mpi_comm_world, req(4), ierror)
		call mpi_irecv(rx3, (jx+2*g)*(kx+2*g)*mx, mpi_real, node(3), 1, &
									 mpi_comm_world, req(2), ierror)
		call mpi_isend(sx1, (jx+2*g)*(kx+2*g)*mx, mpi_real, node(1), 1, &
									 mpi_comm_world, req(1), ierror)
		call mpi_isend(sx3, (jx+2*g)*(kx+2*g)*mx, mpi_real, node(3), 3, &
									 mpi_comm_world, req(3), ierror)
		call mpi_waitall(4, req, mpi_statuses_ignore, ierror)

		call mpi_irecv(rx2, (ix+2*g)*(kx+2*g)*mx, mpi_real, node(2), 4, &
									 mpi_comm_world, req(4), ierror)
		call mpi_irecv(rx4, (ix+2*g)*(kx+2*g)*mx, mpi_real, node(4), 2, &
									 mpi_comm_world, req(2), ierror)
		call mpi_isend(sx2, (ix+2*g)*(kx+2*g)*mx, mpi_real, node(2), 2, &
									 mpi_comm_world, req(1), ierror)
		call mpi_isend(sx4, (ix+2*g)*(kx+2*g)*mx, mpi_real, node(4), 4, &
									 mpi_comm_world, req(3), ierror)
		call mpi_waitall(4, req, mpi_statuses_ignore, ierror)

		call mpi_irecv(rx5, (ix+2*g)*(jx+2*g)*mx, mpi_real, node(5), 6, &
									 mpi_comm_world, req(4), ierror)
		call mpi_irecv(rx6, (ix+2*g)*(jx+2*g)*mx, mpi_real, node(6), 5, &
									 mpi_comm_world, req(2), ierror)
		call mpi_isend(sx5, (ix+2*g)*(jx+2*g)*mx, mpi_real, node(5), 5, &
									 mpi_comm_world, req(1), ierror)
		call mpi_isend(sx6, (ix+2*g)*(jx+2*g)*mx, mpi_real, node(6), 6, &
									 mpi_comm_world, req(3), ierror)
		call mpi_waitall(4, req, mpi_statuses_ignore, ierror)
#endif
!$omp parallel private(i, j, k)
!$omp do
	do k=1-g, kx+g
	do j=1-g, jx+g
	do m=1, mx
		x( ix+m, j, k) = rx1(m, j, k)
		x(-mx+m, j, k) = rx3(m, j, k)
	end do
	end do
	end do
!$omp end do
!$omp do
	do k=1-g, kx+g
	do i=1-g, ix+g
	do m=1, mx
		x(i,  jx+m, k) = rx2(m, i, k)
		x(i, -mx+m, k) = rx4(m, i, k)
	end do
	end do
	end do
!$omp end do
!$omp do
	do j=1-g, jx+g
	do i=1-g, ix+g
	do m=1, mx
		x(i, j,  kx+m) = rx5(m, i, j)
		x(i, j, -mx+m) = rx6(m, i, j)
	end do
	end do
	end do
!$omp end do
!$omp end parallel
end subroutine comm_band_cells_nb

