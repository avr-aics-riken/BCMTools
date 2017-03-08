!
!##################################################################################
!
! BCMTools
!
! Copyright (c) 2011-2014 Institute of Industrial Science, The University of Tokyo.
! All rights reserved.
!
! Copyright (c) 2012-2016 Advanced Institute for Computational Science (AICS), RIKEN.
! All rights reserved.
!
! Copyright (c) 2017 Research Institute for Information Technology (RIIT), Kyushu University.
! All rights reserved.
!
!##################################################################################
!

subroutine bplt3d_open_file(filename, filenamelength, un)
	implicit none
	integer										:: filenamelength
	character(filenamelength)	:: filename
	integer										:: un
	open(unit=un, file=filename, status='unknown', form='unformatted')
end subroutine bplt3d_open_file

subroutine bplt3d_close_file(un)
	implicit none
	integer										:: un
	close(un)
end subroutine bplt3d_close_file

subroutine bplt3d_write_xyz_header(ix, jx, kx, ngrid, un)
	implicit none
  integer											:: ix, jx, kx
	integer											:: ngrid
	integer											:: un
	integer											:: n
	write(un) ngrid
	write(un) (ix, jx, kx, n=1, ngrid)
end subroutine bplt3d_write_xyz_header

subroutine bplt3d_write_xyz_block(x, y, z, ix, jx, kx, un)
	implicit none
  integer											:: ix, jx, kx
  real(4), dimension(ix, jx, kx)	:: x, y, z
	integer											:: un
	integer											:: i, j, k
	write(un)	(((x(i,j,k),i=1,ix),j=1,jx),k=1,kx), &
						(((y(i,j,k),i=1,ix),j=1,jx),k=1,kx), &
						(((z(i,j,k),i=1,ix),j=1,jx),k=1,kx)
end subroutine bplt3d_write_xyz_block

subroutine bplt3d_write_func_header(ix, jx, kx, nvar, ngrid, un)
	implicit none
  integer											:: ix, jx, kx
	integer											:: nvar
	integer											:: ngrid
	integer											:: un
	integer											:: n
	write(un) ngrid
	write(un) (ix, jx, kx, nvar, n=1, ngrid)
end subroutine bplt3d_write_func_header

subroutine bplt3d_write_func_block(p, ux, uy, uz, ix, jx, kx, un)
	implicit none
  integer											:: ix, jx, kx
  real(4), dimension(ix, jx, kx)	:: p, ux, uy, uz
	integer											:: un
	integer											:: i, j, k
	write(un)	(((p(i,j,k),i=1,ix),j=1,jx),k=1,kx), &
						(((ux(i,j,k),i=1,ix),j=1,jx),k=1,kx), &
						(((uy(i,j,k),i=1,ix),j=1,jx),k=1,kx), &
						(((uz(i,j,k),i=1,ix),j=1,jx),k=1,kx)
end subroutine bplt3d_write_func_block
