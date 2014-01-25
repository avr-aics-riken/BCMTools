!
! BCMTools
!
! Copyright (C) 2011-2014 Institute of Industrial Science, The University of Tokyo.
! All rights reserved.
!
! Copyright (c) 2012-2014 Advanced Institute for Computational Science, RIKEN.
! All rights reserved.
!

subroutine bstl_read_cut_1( &
								c0, c1, c2, c3, c4, c5, &
								cid, &
								cut, &
								bid, &
								sz, g)
  implicit none
  integer, dimension(3)   :: sz
  integer                 :: g
  integer                 :: i, j, k
  integer                 :: ix, jx, kx
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: c0, c1, c2, c3, c4, c5
	integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: cid
	real, dimension(0:5, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: cut
	integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: bid
	real										:: cut0, cut1, cut2, cut3, cut4, cut5
	integer									:: bidp
	integer									:: bidp0, bidp1, bidp2, bidp3, bidp4, bidp5
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) &
!$omp					 private(cut0, cut1, cut2, cut3, cut4, cut5) &
!$omp					 private(bidp) &
!$omp					 private(bidp0, bidp1, bidp2, bidp3, bidp4, bidp5) 
!$omp do 
#else
#endif
  do k=1, kx
  do j=1, jx
  do i=1, ix
		cut0 = (cut(0, i, j, k))
		cut1 = (cut(1, i, j, k))
		cut2 = (cut(2, i, j, k))
		cut3 = (cut(3, i, j, k))
		cut4 = (cut(4, i, j, k))
		cut5 = (cut(5, i, j, k))

    cut0 = max(0.0d0, cut0)
    cut1 = max(0.0d0, cut1)
    cut2 = max(0.0d0, cut2)
    cut3 = max(0.0d0, cut3)
    cut4 = max(0.0d0, cut4)
    cut5 = max(0.0d0, cut5)

    cut0 = min(1.0d0, cut0)
    cut1 = min(1.0d0, cut1)
    cut2 = min(1.0d0, cut2)
    cut3 = min(1.0d0, cut3)
    cut4 = min(1.0d0, cut4)
    cut5 = min(1.0d0, cut5)

		bidp = bid(i, j, k)
		bidp0 = ibits(bidp,  0, 5)
		bidp1 = ibits(bidp,  5, 5)
		bidp2 = ibits(bidp, 10, 5)
		bidp3 = ibits(bidp, 15, 5)
		bidp4 = ibits(bidp, 20, 5)
		bidp5 = ibits(bidp, 25, 5)

    c0(i, j, k) = 1.0d0
    c1(i, j, k) = 1.0d0
    c2(i, j, k) = 1.0d0
    c3(i, j, k) = 1.0d0
    c4(i, j, k) = 1.0d0
    c5(i, j, k) = 1.0d0

		if( bidp0 /= 0 ) then
			c0(i, j, k) = cut0
		endif
		if( bidp1 /= 0 ) then
			c1(i, j, k) = cut1
		endif
		if( bidp2 /= 0 ) then
			c2(i, j, k) = cut2
		endif
		if( bidp3 /= 0 ) then
			c3(i, j, k) = cut3
		endif
		if( bidp4 /= 0 ) then
			c4(i, j, k) = cut4
		endif
		if( bidp5 /= 0 ) then
			c5(i, j, k) = cut5
		endif

		cid(i, j, k) = bidp
  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bstl_read_cut_1

subroutine bstl_voxelize_1( &
								c0, c1, c2, c3, c4, c5, &
								cid, &
								sz, g)
  implicit none
  integer, dimension(3)   :: sz
  integer                 :: g
  integer                 :: i, j, k
  integer                 :: ix, jx, kx
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: c0, c1, c2, c3, c4, c5
	integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: cid
	integer									:: cidp
	integer									:: cidp0, cidp1, cidp2, cidp3, cidp4, cidp5
	integer									:: cidw
	integer									:: cidw0, cidw1, cidw2, cidw3, cidw4, cidw5
	integer									:: cide
	integer									:: cide0, cide1, cide2, cide3, cide4, cide5
	integer									:: cids
	integer									:: cids0, cids1, cids2, cids3, cids4, cids5
	integer									:: cidn
	integer									:: cidn0, cidn1, cidn2, cidn3, cidn4, cidn5
	integer									:: cidb
	integer									:: cidb0, cidb1, cidb2, cidb3, cidb4, cidb5
	integer									:: cidt
	integer									:: cidt0, cidt1, cidt2, cidt3, cidt4, cidt5
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) &
!$omp					 private(cidp) &
!$omp					 private(cidp0, cidp1, cidp2, cidp3, cidp4, cidp5) &
!$omp					 private(cidw) &
!$omp					 private(cidw0, cidw1, cidw2, cidw3, cidw4, cidw5) &
!$omp					 private(cide) &
!$omp					 private(cide0, cide1, cide2, cide3, cide4, cide5) &
!$omp					 private(cids) &
!$omp					 private(cids0, cids1, cids2, cids3, cids4, cids5) &
!$omp					 private(cidn) &
!$omp					 private(cidn0, cidn1, cidn2, cidn3, cidn4, cidn5) &
!$omp					 private(cidb) &
!$omp					 private(cidb0, cidb1, cidb2, cidb3, cidb4, cidb5) &
!$omp					 private(cidt) &
!$omp					 private(cidt0, cidt1, cidt2, cidt3, cidt4, cidt5) 
!$omp do 
#else
#endif
  do k=1, kx
  do j=1, jx
  do i=1, ix
		cidp = cid(i, j, k)
		cidp0 = ibits(cidp,  0, 5)
		cidp1 = ibits(cidp,  5, 5)
		cidp2 = ibits(cidp, 10, 5)
		cidp3 = ibits(cidp, 15, 5)
		cidp4 = ibits(cidp, 20, 5)
		cidp5 = ibits(cidp, 25, 5)

		cidw = cid(i-1, j, k)
		cidw0 = ibits(cidw,  0, 5)
		cidw1 = ibits(cidw,  5, 5)
		cidw2 = ibits(cidw, 10, 5)
		cidw3 = ibits(cidw, 15, 5)
		cidw4 = ibits(cidw, 20, 5)
		cidw5 = ibits(cidw, 25, 5)

		cide = cid(i+1, j, k)
		cide0 = ibits(cide,  0, 5)
		cide1 = ibits(cide,  5, 5)
		cide2 = ibits(cide, 10, 5)
		cide3 = ibits(cide, 15, 5)
		cide4 = ibits(cide, 20, 5)
		cide5 = ibits(cide, 25, 5)

		cids = cid(i, j-1, k)
		cids0 = ibits(cids,  0, 5)
		cids1 = ibits(cids,  5, 5)
		cids2 = ibits(cids, 10, 5)
		cids3 = ibits(cids, 15, 5)
		cids4 = ibits(cids, 20, 5)
		cids5 = ibits(cids, 25, 5)

		cidn = cid(i, j+1, k)
		cidn0 = ibits(cidn,  0, 5)
		cidn1 = ibits(cidn,  5, 5)
		cidn2 = ibits(cidn, 10, 5)
		cidn3 = ibits(cidn, 15, 5)
		cidn4 = ibits(cidn, 20, 5)
		cidn5 = ibits(cidn, 25, 5)

		cidb = cid(i, j, k-1)
		cidb0 = ibits(cidb,  0, 5)
		cidb1 = ibits(cidb,  5, 5)
		cidb2 = ibits(cidb, 10, 5)
		cidb3 = ibits(cidb, 15, 5)
		cidb4 = ibits(cidb, 20, 5)
		cidb5 = ibits(cidb, 25, 5)

		cidt = cid(i, j, k+1)
		cidt0 = ibits(cidt,  0, 5)
		cidt1 = ibits(cidt,  5, 5)
		cidt2 = ibits(cidt, 10, 5)
		cidt3 = ibits(cidt, 15, 5)
		cidt4 = ibits(cidt, 20, 5)
		cidt5 = ibits(cidt, 25, 5)

		if( cidp0 /= 0 ) then
			c0(i  , j, k) = 0.5
			c1(i-1, j, k) = 0.5
			cidw1 = cidp0
		endif
		if( cidp1 /= 0 ) then
			c1(i  , j, k) = 0.5
			c0(i+1, j, k) = 0.5
			cide0 = cidp1
		endif
		if( cidp2 /= 0 ) then
			c2(i, j  , k) = 0.5
			c3(i, j-1, k) = 0.5
			cids3 = cidp2
		endif
		if( cidp3 /= 0 ) then
			c3(i, j  , k) = 0.5
			c2(i, j+1, k) = 0.5
			cidn2 = cidp3
		endif
		if( cidp4 /= 0 ) then
			c4(i, j, k  ) = 0.5
			c5(i, j, k-1) = 0.5
			cidb5 = cidp4
		endif
		if( cidp5 /= 0 ) then
			c5(i, j, k  ) = 0.5
			c4(i, j, k+1) = 0.5
			cidt4 = cidp5
		endif

		cidw =  cidw0 &
					+ cidw1*32 &
					+ cidw2*32*32 &
					+ cidw3*32*32*32 &
					+ cidw4*32*32*32*32 &
					+ cidw5*32*32*32*32*32 
		cide =  cide0 &
					+ cide1*32 &
					+ cide2*32*32 &
					+ cide3*32*32*32 &
					+ cide4*32*32*32*32 &
					+ cide5*32*32*32*32*32 
		cids =  cids0 &
					+ cids1*32 &
					+ cids2*32*32 &
					+ cids3*32*32*32 &
					+ cids4*32*32*32*32 &
					+ cids5*32*32*32*32*32 
		cidn =  cidn0 &
					+ cidn1*32 &
					+ cidn2*32*32 &
					+ cidn3*32*32*32 &
					+ cidn4*32*32*32*32 &
					+ cidn5*32*32*32*32*32 
		cidb =  cidb0 &
					+ cidb1*32 &
					+ cidb2*32*32 &
					+ cidb3*32*32*32 &
					+ cidb4*32*32*32*32 &
					+ cidb5*32*32*32*32*32 
		cidt =  cidt0 &
					+ cidt1*32 &
					+ cidt2*32*32 &
					+ cidt3*32*32*32 &
					+ cidt4*32*32*32*32 &
					+ cidt5*32*32*32*32*32 

		cid(i-1, j, k) = cidw
		cid(i+1, j, k) = cide
		cid(i, j-1, k) = cids
		cid(i, j+1, k) = cidn
		cid(i, j, k-1) = cidb
		cid(i, j, k+1) = cidt
  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bstl_voxelize_1

subroutine bstl_cutoff_1( &
								c0, c1, c2, c3, c4, c5, &
								cid, &
								eps, &
								sz, g)
  implicit none
  integer, dimension(3)   :: sz
  integer                 :: g
  integer                 :: i, j, k
  integer                 :: ix, jx, kx
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: c0, c1, c2, c3, c4, c5
	integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: cid
	real										:: eps
	integer									:: cidp
	integer									:: cidp0, cidp1, cidp2, cidp3, cidp4, cidp5
	integer									:: cidw
	integer									:: cidw0, cidw1, cidw2, cidw3, cidw4, cidw5
	integer									:: cide
	integer									:: cide0, cide1, cide2, cide3, cide4, cide5
	integer									:: cids
	integer									:: cids0, cids1, cids2, cids3, cids4, cids5
	integer									:: cidn
	integer									:: cidn0, cidn1, cidn2, cidn3, cidn4, cidn5
	integer									:: cidb
	integer									:: cidb0, cidb1, cidb2, cidb3, cidb4, cidb5
	integer									:: cidt
	integer									:: cidt0, cidt1, cidt2, cidt3, cidt4, cidt5
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) &
!$omp					 private(cidp) &
!$omp					 private(cidp0, cidp1, cidp2, cidp3, cidp4, cidp5) &
!$omp					 private(cidw) &
!$omp					 private(cidw0, cidw1, cidw2, cidw3, cidw4, cidw5) &
!$omp					 private(cide) &
!$omp					 private(cide0, cide1, cide2, cide3, cide4, cide5) &
!$omp					 private(cids) &
!$omp					 private(cids0, cids1, cids2, cids3, cids4, cids5) &
!$omp					 private(cidn) &
!$omp					 private(cidn0, cidn1, cidn2, cidn3, cidn4, cidn5) &
!$omp					 private(cidb) &
!$omp					 private(cidb0, cidb1, cidb2, cidb3, cidb4, cidb5) &
!$omp					 private(cidt) &
!$omp					 private(cidt0, cidt1, cidt2, cidt3, cidt4, cidt5) 
!$omp do 
#else
#endif
  do k=1, kx
  do j=1, jx
  do i=1, ix
    c0(i, j, k) = max(eps, c0(i, j, k))
    c1(i, j, k) = max(eps, c1(i, j, k))
    c2(i, j, k) = max(eps, c2(i, j, k))
    c3(i, j, k) = max(eps, c3(i, j, k))
    c4(i, j, k) = max(eps, c4(i, j, k))
    c5(i, j, k) = max(eps, c5(i, j, k))
  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bstl_cutoff_1

subroutine bstl_symmetrize_1( &
								c0, c1, c2, c3, c4, c5, &
								cid, &
								sz, g)
  implicit none
  integer, dimension(3)   :: sz
  integer                 :: g
  integer                 :: i, j, k
  integer                 :: ix, jx, kx
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: c0, c1, c2, c3, c4, c5
	integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: cid
	integer									:: cidp
	integer									:: cidp0, cidp1, cidp2, cidp3, cidp4, cidp5
	integer									:: cidw
	integer									:: cidw0, cidw1, cidw2, cidw3, cidw4, cidw5
	integer									:: cide
	integer									:: cide0, cide1, cide2, cide3, cide4, cide5
	integer									:: cids
	integer									:: cids0, cids1, cids2, cids3, cids4, cids5
	integer									:: cidn
	integer									:: cidn0, cidn1, cidn2, cidn3, cidn4, cidn5
	integer									:: cidb
	integer									:: cidb0, cidb1, cidb2, cidb3, cidb4, cidb5
	integer									:: cidt
	integer									:: cidt0, cidt1, cidt2, cidt3, cidt4, cidt5
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) &
!$omp					 private(cidp) &
!$omp					 private(cidp0, cidp1, cidp2, cidp3, cidp4, cidp5) &
!$omp					 private(cidw) &
!$omp					 private(cidw0, cidw1, cidw2, cidw3, cidw4, cidw5) &
!$omp					 private(cide) &
!$omp					 private(cide0, cide1, cide2, cide3, cide4, cide5) &
!$omp					 private(cids) &
!$omp					 private(cids0, cids1, cids2, cids3, cids4, cids5) &
!$omp					 private(cidn) &
!$omp					 private(cidn0, cidn1, cidn2, cidn3, cidn4, cidn5) &
!$omp					 private(cidb) &
!$omp					 private(cidb0, cidb1, cidb2, cidb3, cidb4, cidb5) &
!$omp					 private(cidt) &
!$omp					 private(cidt0, cidt1, cidt2, cidt3, cidt4, cidt5) 
!$omp do 
#else
#endif
  do k=1, kx
  do j=1, jx
  do i=1, ix
		cidp = cid(i, j, k)
		cidp0 = ibits(cidp,  0, 5)
		cidp1 = ibits(cidp,  5, 5)
		cidp2 = ibits(cidp, 10, 5)
		cidp3 = ibits(cidp, 15, 5)
		cidp4 = ibits(cidp, 20, 5)
		cidp5 = ibits(cidp, 25, 5)

		cidw = cid(i-1, j, k)
		cidw0 = ibits(cidw,  0, 5)
		cidw1 = ibits(cidw,  5, 5)
		cidw2 = ibits(cidw, 10, 5)
		cidw3 = ibits(cidw, 15, 5)
		cidw4 = ibits(cidw, 20, 5)
		cidw5 = ibits(cidw, 25, 5)

		cide = cid(i+1, j, k)
		cide0 = ibits(cide,  0, 5)
		cide1 = ibits(cide,  5, 5)
		cide2 = ibits(cide, 10, 5)
		cide3 = ibits(cide, 15, 5)
		cide4 = ibits(cide, 20, 5)
		cide5 = ibits(cide, 25, 5)

		cids = cid(i, j-1, k)
		cids0 = ibits(cids,  0, 5)
		cids1 = ibits(cids,  5, 5)
		cids2 = ibits(cids, 10, 5)
		cids3 = ibits(cids, 15, 5)
		cids4 = ibits(cids, 20, 5)
		cids5 = ibits(cids, 25, 5)

		cidn = cid(i, j+1, k)
		cidn0 = ibits(cidn,  0, 5)
		cidn1 = ibits(cidn,  5, 5)
		cidn2 = ibits(cidn, 10, 5)
		cidn3 = ibits(cidn, 15, 5)
		cidn4 = ibits(cidn, 20, 5)
		cidn5 = ibits(cidn, 25, 5)

		cidb = cid(i, j, k-1)
		cidb0 = ibits(cidb,  0, 5)
		cidb1 = ibits(cidb,  5, 5)
		cidb2 = ibits(cidb, 10, 5)
		cidb3 = ibits(cidb, 15, 5)
		cidb4 = ibits(cidb, 20, 5)
		cidb5 = ibits(cidb, 25, 5)

		cidt = cid(i, j, k+1)
		cidt0 = ibits(cidt,  0, 5)
		cidt1 = ibits(cidt,  5, 5)
		cidt2 = ibits(cidt, 10, 5)
		cidt3 = ibits(cidt, 15, 5)
		cidt4 = ibits(cidt, 20, 5)
		cidt5 = ibits(cidt, 25, 5)

		if( cidp0 /= 0 ) then
			c1(i-1, j, k) = 1.0d0 - c0(i  , j, k) 
			cidw1 = cidp0
		endif
		if( cidp1 /= 0 ) then
			c0(i+1, j, k) = 1.0d0 - c1(i  , j, k) 
			cide0 = cidp1
		endif
		if( cidp2 /= 0 ) then
			c3(i, j-1, k) = 1.0d0 - c2(i, j  , k) 
			cids3 = cidp2
		endif
		if( cidp3 /= 0 ) then
			c2(i, j+1, k) = 1.0d0 - c3(i, j  , k) 
			cidn2 = cidp3
		endif
		if( cidp4 /= 0 ) then
			c5(i, j, k-1) = 1.0d0 - c4(i, j, k  ) 
			cidb5 = cidp4
		endif
		if( cidp5 /= 0 ) then
			c4(i, j, k+1) = 1.0d0 - c5(i, j, k  ) 
			cidt4 = cidp5
		endif

		cidw =  cidw0 &
					+ cidw1*32 &
					+ cidw2*32*32 &
					+ cidw3*32*32*32 &
					+ cidw4*32*32*32*32 &
					+ cidw5*32*32*32*32*32 
		cide =  cide0 &
					+ cide1*32 &
					+ cide2*32*32 &
					+ cide3*32*32*32 &
					+ cide4*32*32*32*32 &
					+ cide5*32*32*32*32*32 
		cids =  cids0 &
					+ cids1*32 &
					+ cids2*32*32 &
					+ cids3*32*32*32 &
					+ cids4*32*32*32*32 &
					+ cids5*32*32*32*32*32 
		cidn =  cidn0 &
					+ cidn1*32 &
					+ cidn2*32*32 &
					+ cidn3*32*32*32 &
					+ cidn4*32*32*32*32 &
					+ cidn5*32*32*32*32*32 
		cidb =  cidb0 &
					+ cidb1*32 &
					+ cidb2*32*32 &
					+ cidb3*32*32*32 &
					+ cidb4*32*32*32*32 &
					+ cidb5*32*32*32*32*32 
		cidt =  cidt0 &
					+ cidt1*32 &
					+ cidt2*32*32 &
					+ cidt3*32*32*32 &
					+ cidt4*32*32*32*32 &
					+ cidt5*32*32*32*32*32 

		cid(i-1, j, k) = cidw
		cid(i+1, j, k) = cide
		cid(i, j-1, k) = cids
		cid(i, j+1, k) = cidn
		cid(i, j, k-1) = cidb
		cid(i, j, k+1) = cidt
  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bstl_symmetrize_1

subroutine bstl_detect_zerocut_1( &
								c0, c1, c2, c3, c4, c5, &
								cid, &
                n, &
								sz, g)
  implicit none
  integer, dimension(3)   :: sz
  integer                 :: g
  integer                 :: i, j, k
  integer                 :: ix, jx, kx
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: c0, c1, c2, c3, c4, c5
	integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: cid
	integer									:: cidp
	integer									:: cidp0, cidp1, cidp2, cidp3, cidp4, cidp5
	integer									:: cidw
	integer									:: cidw0, cidw1, cidw2, cidw3, cidw4, cidw5
	integer									:: cide
	integer									:: cide0, cide1, cide2, cide3, cide4, cide5
	integer									:: cids
	integer									:: cids0, cids1, cids2, cids3, cids4, cids5
	integer									:: cidn
	integer									:: cidn0, cidn1, cidn2, cidn3, cidn4, cidn5
	integer									:: cidb
	integer									:: cidb0, cidb1, cidb2, cidb3, cidb4, cidb5
	integer									:: cidt
	integer									:: cidt0, cidt1, cidt2, cidt3, cidt4, cidt5
	integer									:: n
	real										:: cut_min
	integer									:: cid_min
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) &
!$omp					 private(cidp) &
!$omp					 private(cidp0, cidp1, cidp2, cidp3, cidp4, cidp5) &
!$omp					 private(cidw) &
!$omp					 private(cidw0, cidw1, cidw2, cidw3, cidw4, cidw5) &
!$omp					 private(cide) &
!$omp					 private(cide0, cide1, cide2, cide3, cide4, cide5) &
!$omp					 private(cids) &
!$omp					 private(cids0, cids1, cids2, cids3, cids4, cids5) &
!$omp					 private(cidn) &
!$omp					 private(cidn0, cidn1, cidn2, cidn3, cidn4, cidn5) &
!$omp					 private(cidb) &
!$omp					 private(cidb0, cidb1, cidb2, cidb3, cidb4, cidb5) &
!$omp					 private(cidt) &
!$omp					 private(cidt0, cidt1, cidt2, cidt3, cidt4, cidt5) &
!$omp					 private(cut_min) &
!$omp					 private(cid_min) 
!$omp do reduction(+: n)
#else
#endif
  do k=1, kx
  do j=1, jx
  do i=1, ix
		cidp = cid(i, j, k)
		cidp0 = ibits(cidp,  0, 5)
		cidp1 = ibits(cidp,  5, 5)
		cidp2 = ibits(cidp, 10, 5)
		cidp3 = ibits(cidp, 15, 5)
		cidp4 = ibits(cidp, 20, 5)
		cidp5 = ibits(cidp, 25, 5)

		cidw = cid(i-1, j, k)
		cidw0 = ibits(cidw,  0, 5)
		cidw1 = ibits(cidw,  5, 5)
		cidw2 = ibits(cidw, 10, 5)
		cidw3 = ibits(cidw, 15, 5)
		cidw4 = ibits(cidw, 20, 5)
		cidw5 = ibits(cidw, 25, 5)

		cide = cid(i+1, j, k)
		cide0 = ibits(cide,  0, 5)
		cide1 = ibits(cide,  5, 5)
		cide2 = ibits(cide, 10, 5)
		cide3 = ibits(cide, 15, 5)
		cide4 = ibits(cide, 20, 5)
		cide5 = ibits(cide, 25, 5)

		cids = cid(i, j-1, k)
		cids0 = ibits(cids,  0, 5)
		cids1 = ibits(cids,  5, 5)
		cids2 = ibits(cids, 10, 5)
		cids3 = ibits(cids, 15, 5)
		cids4 = ibits(cids, 20, 5)
		cids5 = ibits(cids, 25, 5)

		cidn = cid(i, j+1, k)
		cidn0 = ibits(cidn,  0, 5)
		cidn1 = ibits(cidn,  5, 5)
		cidn2 = ibits(cidn, 10, 5)
		cidn3 = ibits(cidn, 15, 5)
		cidn4 = ibits(cidn, 20, 5)
		cidn5 = ibits(cidn, 25, 5)

		cidb = cid(i, j, k-1)
		cidb0 = ibits(cidb,  0, 5)
		cidb1 = ibits(cidb,  5, 5)
		cidb2 = ibits(cidb, 10, 5)
		cidb3 = ibits(cidb, 15, 5)
		cidb4 = ibits(cidb, 20, 5)
		cidb5 = ibits(cidb, 25, 5)

		cidt = cid(i, j, k+1)
		cidt0 = ibits(cidt,  0, 5)
		cidt1 = ibits(cidt,  5, 5)
		cidt2 = ibits(cidt, 10, 5)
		cidt3 = ibits(cidt, 15, 5)
		cidt4 = ibits(cidt, 20, 5)
		cidt5 = ibits(cidt, 25, 5)

		cut_min = 1.0d0	
		cid_min = 0
		if( cut_min > c0(i, j, k) ) then
			cut_min = c0(i, j, k)
			cid_min = cidp0
		endif
		if( cut_min > c1(i, j, k) ) then
			cut_min = c1(i, j, k)
			cid_min = cidp1
		endif
		if( cut_min > c2(i, j, k) ) then
			cut_min = c2(i, j, k)
			cid_min = cidp2
		endif
		if( cut_min > c3(i, j, k) ) then
			cut_min = c3(i, j, k)
			cid_min = cidp3
		endif
		if( cut_min > c4(i, j, k) ) then
			cut_min = c4(i, j, k)
			cid_min = cidp4
		endif
		if( cut_min > c5(i, j, k) ) then
			cut_min = c5(i, j, k)
			cid_min = cidp5
		endif

		if( cut_min == 0.0d0 ) then
			c0(i, j, k) = 1.0d0
			c1(i, j, k) = 1.0d0
			c2(i, j, k) = 1.0d0
			c3(i, j, k) = 1.0d0
			c4(i, j, k) = 1.0d0
			c5(i, j, k) = 1.0d0
			cid(i, j, k) =  cid_min &
										+ cid_min*32 &
										+ cid_min*32*32 &
										+ cid_min*32*32*32 &
										+ cid_min*32*32*32*32 &
										+ cid_min*32*32*32*32*32 

      c1(i-1, j, k) = 1.0d0
      cidw1 = cid_min
      c0(i+1, j, k) = 1.0d0
      cide0 = cid_min
      c3(i, j-1, k) = 1.0d0
      cids3 = cid_min
      c2(i, j+1, k) = 1.0d0
      cidn2 = cid_min
      c5(i, j, k-1) = 1.0d0
      cidb5 = cid_min
      c4(i, j, k+1) = 1.0d0
      cidt4 = cid_min

			n = n + 1
		endif

		cidw =  cidw0 &
					+ cidw1*32 &
					+ cidw2*32*32 &
					+ cidw3*32*32*32 &
					+ cidw4*32*32*32*32 &
					+ cidw5*32*32*32*32*32 
		cide =  cide0 &
					+ cide1*32 &
					+ cide2*32*32 &
					+ cide3*32*32*32 &
					+ cide4*32*32*32*32 &
					+ cide5*32*32*32*32*32 
		cids =  cids0 &
					+ cids1*32 &
					+ cids2*32*32 &
					+ cids3*32*32*32 &
					+ cids4*32*32*32*32 &
					+ cids5*32*32*32*32*32 
		cidn =  cidn0 &
					+ cidn1*32 &
					+ cidn2*32*32 &
					+ cidn3*32*32*32 &
					+ cidn4*32*32*32*32 &
					+ cidn5*32*32*32*32*32 
		cidb =  cidb0 &
					+ cidb1*32 &
					+ cidb2*32*32 &
					+ cidb3*32*32*32 &
					+ cidb4*32*32*32*32 &
					+ cidb5*32*32*32*32*32 
		cidt =  cidt0 &
					+ cidt1*32 &
					+ cidt2*32*32 &
					+ cidt3*32*32*32 &
					+ cidt4*32*32*32*32 &
					+ cidt5*32*32*32*32*32 

		cid(i-1, j, k) = cidw
		cid(i+1, j, k) = cide
		cid(i, j-1, k) = cids
		cid(i, j+1, k) = cidn
		cid(i, j, k-1) = cidb
		cid(i, j, k+1) = cidt
  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bstl_detect_zerocut_1

subroutine bstl_fill_holes_1( &
								c0, c1, c2, c3, c4, c5, &
								cid, &
								n, &
								sz, g)
  implicit none
  integer, dimension(3)   :: sz
  integer                 :: g
  integer                 :: i, j, k
  integer                 :: ix, jx, kx
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: c0, c1, c2, c3, c4, c5
	integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: cid
	integer									:: n
	integer									:: cidp
	integer									:: cidp0, cidp1, cidp2, cidp3, cidp4, cidp5
	integer									:: cidw
	integer									:: cidw0, cidw1, cidw2, cidw3, cidw4, cidw5
	integer									:: cide
	integer									:: cide0, cide1, cide2, cide3, cide4, cide5
	integer									:: cids
	integer									:: cids0, cids1, cids2, cids3, cids4, cids5
	integer									:: cidn
	integer									:: cidn0, cidn1, cidn2, cidn3, cidn4, cidn5
	integer									:: cidb
	integer									:: cidb0, cidb1, cidb2, cidb3, cidb4, cidb5
	integer									:: cidt
	integer									:: cidt0, cidt1, cidt2, cidt3, cidt4, cidt5
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
	n = 0
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) &
!$omp					 private(cidp) &
!$omp					 private(cidp0, cidp1, cidp2, cidp3, cidp4, cidp5) &
!$omp					 private(cidw) &
!$omp					 private(cidw0, cidw1, cidw2, cidw3, cidw4, cidw5) &
!$omp					 private(cide) &
!$omp					 private(cide0, cide1, cide2, cide3, cide4, cide5) &
!$omp					 private(cids) &
!$omp					 private(cids0, cids1, cids2, cids3, cids4, cids5) &
!$omp					 private(cidn) &
!$omp					 private(cidn0, cidn1, cidn2, cidn3, cidn4, cidn5) &
!$omp					 private(cidb) &
!$omp					 private(cidb0, cidb1, cidb2, cidb3, cidb4, cidb5) &
!$omp					 private(cidt) &
!$omp					 private(cidt0, cidt1, cidt2, cidt3, cidt4, cidt5)
!$omp do reduction(+: n)
#else
#endif
  do k=1, kx
  do j=1, jx
  do i=1, ix
		cidp = cid(i, j, k)
		cidp0 = ibits(cidp,  0, 5)
		cidp1 = ibits(cidp,  5, 5)
		cidp2 = ibits(cidp, 10, 5)
		cidp3 = ibits(cidp, 15, 5)
		cidp4 = ibits(cidp, 20, 5)
		cidp5 = ibits(cidp, 25, 5)

		cidw = cid(i-1, j, k)
		cidw0 = ibits(cidw,  0, 5)
		cidw1 = ibits(cidw,  5, 5)
		cidw2 = ibits(cidw, 10, 5)
		cidw3 = ibits(cidw, 15, 5)
		cidw4 = ibits(cidw, 20, 5)
		cidw5 = ibits(cidw, 25, 5)

		cide = cid(i+1, j, k)
		cide0 = ibits(cide,  0, 5)
		cide1 = ibits(cide,  5, 5)
		cide2 = ibits(cide, 10, 5)
		cide3 = ibits(cide, 15, 5)
		cide4 = ibits(cide, 20, 5)
		cide5 = ibits(cide, 25, 5)

		cids = cid(i, j-1, k)
		cids0 = ibits(cids,  0, 5)
		cids1 = ibits(cids,  5, 5)
		cids2 = ibits(cids, 10, 5)
		cids3 = ibits(cids, 15, 5)
		cids4 = ibits(cids, 20, 5)
		cids5 = ibits(cids, 25, 5)

		cidn = cid(i, j+1, k)
		cidn0 = ibits(cidn,  0, 5)
		cidn1 = ibits(cidn,  5, 5)
		cidn2 = ibits(cidn, 10, 5)
		cidn3 = ibits(cidn, 15, 5)
		cidn4 = ibits(cidn, 20, 5)
		cidn5 = ibits(cidn, 25, 5)

		cidb = cid(i, j, k-1)
		cidb0 = ibits(cidb,  0, 5)
		cidb1 = ibits(cidb,  5, 5)
		cidb2 = ibits(cidb, 10, 5)
		cidb3 = ibits(cidb, 15, 5)
		cidb4 = ibits(cidb, 20, 5)
		cidb5 = ibits(cidb, 25, 5)

		cidt = cid(i, j, k+1)
		cidt0 = ibits(cidt,  0, 5)
		cidt1 = ibits(cidt,  5, 5)
		cidt2 = ibits(cidt, 10, 5)
		cidt3 = ibits(cidt, 15, 5)
		cidt4 = ibits(cidt, 20, 5)
		cidt5 = ibits(cidt, 25, 5)

		if( cidp0 == 0 ) then
			if( (cids0 /= 0 .or. cidp2 /=0 .or. cidw2 /=0) .and. &
					(cidn0 /= 0 .or. cidp3 /=0 .or. cidw3 /=0) .and. &
					(cidb0 /= 0 .or. cidp4 /=0 .or. cidw4 /=0) .and. &
					(cidt0 /= 0 .or. cidp5 /=0 .or. cidw5 /=0) ) then
				cidp0 = 1
				cidw1 = 1
				c0(i  , j, k) = 0.5
				c1(i-1, j, k) = 0.5
				n = n + 1
			endif
		endif

		if( cidp1 == 0 ) then
			if( (cids1 /= 0 .or. cidp2 /=0 .or. cide2 /=0) .and. &
					(cidn1 /= 0 .or. cidp3 /=0 .or. cide3 /=0) .and. &
					(cidb1 /= 0 .or. cidp4 /=0 .or. cide4 /=0) .and. &
					(cidt1 /= 0 .or. cidp5 /=0 .or. cide5 /=0) ) then
				cidp1 = 1
				cide0 = 1
				c1(i  , j, k) = 0.5
				c0(i+1, j, k) = 0.5
				n = n + 1
			endif
		endif

		if( cidp2 == 0 ) then
			if( (cidw2 /= 0 .or. cidp0 /=0 .or. cids0 /=0) .and. &
					(cide2 /= 0 .or. cidp1 /=0 .or. cids1 /=0) .and. &
					(cidb2 /= 0 .or. cidp4 /=0 .or. cids4 /=0) .and. &
					(cidt2 /= 0 .or. cidp5 /=0 .or. cids5 /=0) ) then
				cidp2 = 1
				cids3 = 1
				c2(i, j  , k) = 0.5
				c3(i, j-1, k) = 0.5
				n = n + 1
			endif
		endif

		if( cidp3 == 0 ) then
			if( (cidw3 /= 0 .or. cidp0 /=0 .or. cidn0 /=0) .and. &
					(cide3 /= 0 .or. cidp1 /=0 .or. cidn1 /=0) .and. &
					(cidb3 /= 0 .or. cidp4 /=0 .or. cidn4 /=0) .and. &
					(cidt3 /= 0 .or. cidp5 /=0 .or. cidn5 /=0) ) then
				cidp3 = 1
				cidn2 = 1
				c3(i, j  , k) = 0.5
				c2(i, j+1, k) = 0.5
				n = n + 1
			endif
		endif

		if( cidp4 == 0 ) then
			if( (cidw4 /= 0 .or. cidp0 /=0 .or. cidb0 /=0) .and. &
					(cide4 /= 0 .or. cidp1 /=0 .or. cidb1 /=0) .and. &
					(cids4 /= 0 .or. cidp2 /=0 .or. cidb2 /=0) .and. &
					(cidn4 /= 0 .or. cidp3 /=0 .or. cidb3 /=0) ) then
				cidp4 = 1
				cidb5 = 1
				c4(i, j, k  ) = 0.5
				c5(i, j, k-1) = 0.5
				n = n + 1
			endif
		endif

		if( cidp5 == 0 ) then
			if( (cidw5 /= 0 .or. cidp0 /=0 .or. cidt0 /=0) .and. &
					(cide5 /= 0 .or. cidp1 /=0 .or. cidt1 /=0) .and. &
					(cids5 /= 0 .or. cidp2 /=0 .or. cidt2 /=0) .and. &
					(cidn5 /= 0 .or. cidp3 /=0 .or. cidt3 /=0) ) then
				cidp5 = 1
				cidt4 = 1
				c5(i, j, k  ) = 0.5
				c4(i, j, k+1) = 0.5
				n = n + 1
			endif
		endif

		cidp =  cidp0 &
					+ cidp1*32 &
					+ cidp2*32*32 &
					+ cidp3*32*32*32 &
					+ cidp4*32*32*32*32 &
					+ cidp5*32*32*32*32*32 
		cidw =  cidw0 &
					+ cidw1*32 &
					+ cidw2*32*32 &
					+ cidw3*32*32*32 &
					+ cidw4*32*32*32*32 &
					+ cidw5*32*32*32*32*32 
		cide =  cide0 &
					+ cide1*32 &
					+ cide2*32*32 &
					+ cide3*32*32*32 &
					+ cide4*32*32*32*32 &
					+ cide5*32*32*32*32*32 
		cids =  cids0 &
					+ cids1*32 &
					+ cids2*32*32 &
					+ cids3*32*32*32 &
					+ cids4*32*32*32*32 &
					+ cids5*32*32*32*32*32 
		cidn =  cidn0 &
					+ cidn1*32 &
					+ cidn2*32*32 &
					+ cidn3*32*32*32 &
					+ cidn4*32*32*32*32 &
					+ cidn5*32*32*32*32*32 
		cidb =  cidb0 &
					+ cidb1*32 &
					+ cidb2*32*32 &
					+ cidb3*32*32*32 &
					+ cidb4*32*32*32*32 &
					+ cidb5*32*32*32*32*32 
		cidt =  cidt0 &
					+ cidt1*32 &
					+ cidt2*32*32 &
					+ cidt3*32*32*32 &
					+ cidt4*32*32*32*32 &
					+ cidt5*32*32*32*32*32 

		cid(i, j, k) = cidp
		cid(i-1, j, k) = cidw
		cid(i+1, j, k) = cide
		cid(i, j-1, k) = cids
		cid(i, j+1, k) = cidn
		cid(i, j, k-1) = cidb
		cid(i, j, k+1) = cidt
  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bstl_fill_holes_1

subroutine bstl_read_cut( &
								c0, c1, c2, c3, c4, c5, &
								cid0, cid1, cid2, cid3, cid4, cid5, &
								cut, &
								bid, &
								sz, g)
  implicit none
  integer, dimension(3)   :: sz
  integer                 :: g
  integer                 :: i, j, k
  integer                 :: ix, jx, kx
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: c0, c1, c2, c3, c4, c5
	integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: cid0, cid1, cid2, cid3, cid4, cid5
	real, dimension(0:5, 1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: cut
	integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: bid
	real										:: cut0, cut1, cut2, cut3, cut4, cut5
	integer									:: bidp
	integer									:: bidp0, bidp1, bidp2, bidp3, bidp4, bidp5
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) &
!$omp					 private(cut0, cut1, cut2, cut3, cut4, cut5) &
!$omp					 private(bidp) &
!$omp					 private(bidp0, bidp1, bidp2, bidp3, bidp4, bidp5) 
!$omp do 
#else
#endif
  do k=1, kx
  do j=1, jx
  do i=1, ix
		cut0 = (cut(0, i, j, k))
		cut1 = (cut(1, i, j, k))
		cut2 = (cut(2, i, j, k))
		cut3 = (cut(3, i, j, k))
		cut4 = (cut(4, i, j, k))
		cut5 = (cut(5, i, j, k))

    cut0 = max(0.0d0, cut0)
    cut1 = max(0.0d0, cut1)
    cut2 = max(0.0d0, cut2)
    cut3 = max(0.0d0, cut3)
    cut4 = max(0.0d0, cut4)
    cut5 = max(0.0d0, cut5)

    cut0 = min(1.0d0, cut0)
    cut1 = min(1.0d0, cut1)
    cut2 = min(1.0d0, cut2)
    cut3 = min(1.0d0, cut3)
    cut4 = min(1.0d0, cut4)
    cut5 = min(1.0d0, cut5)

		bidp = bid(i, j, k)
		bidp0 = ibits(bidp,  0, 5)
		bidp1 = ibits(bidp,  5, 5)
		bidp2 = ibits(bidp, 10, 5)
		bidp3 = ibits(bidp, 15, 5)
		bidp4 = ibits(bidp, 20, 5)
		bidp5 = ibits(bidp, 25, 5)

    c0(i, j, k) = 1.0d0
    c1(i, j, k) = 1.0d0
    c2(i, j, k) = 1.0d0
    c3(i, j, k) = 1.0d0
    c4(i, j, k) = 1.0d0
    c5(i, j, k) = 1.0d0

		if( bidp0 /= 0 ) then
			c0(i, j, k) = cut0
		endif
		if( bidp1 /= 0 ) then
			c1(i, j, k) = cut1
		endif
		if( bidp2 /= 0 ) then
			c2(i, j, k) = cut2
		endif
		if( bidp3 /= 0 ) then
			c3(i, j, k) = cut3
		endif
		if( bidp4 /= 0 ) then
			c4(i, j, k) = cut4
		endif
		if( bidp5 /= 0 ) then
			c5(i, j, k) = cut5
		endif

		cid0(i, j, k) = bidp0
		cid1(i, j, k) = bidp1
		cid2(i, j, k) = bidp2
		cid3(i, j, k) = bidp3
		cid4(i, j, k) = bidp4
		cid5(i, j, k) = bidp5
  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bstl_read_cut

subroutine bstl_voxelize( &
								c0, c1, c2, c3, c4, c5, &
								cid0, cid1, cid2, cid3, cid4, cid5, &
								sz, g)
  implicit none
  integer, dimension(3)   :: sz
  integer                 :: g
  integer                 :: i, j, k
  integer                 :: ix, jx, kx
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: c0, c1, c2, c3, c4, c5
	integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: cid0, cid1, cid2, cid3, cid4, cid5
	integer									:: cidp
	integer									:: cidp0, cidp1, cidp2, cidp3, cidp4, cidp5
	integer									:: cidw
	integer									:: cidw0, cidw1, cidw2, cidw3, cidw4, cidw5
	integer									:: cide
	integer									:: cide0, cide1, cide2, cide3, cide4, cide5
	integer									:: cids
	integer									:: cids0, cids1, cids2, cids3, cids4, cids5
	integer									:: cidn
	integer									:: cidn0, cidn1, cidn2, cidn3, cidn4, cidn5
	integer									:: cidb
	integer									:: cidb0, cidb1, cidb2, cidb3, cidb4, cidb5
	integer									:: cidt
	integer									:: cidt0, cidt1, cidt2, cidt3, cidt4, cidt5
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) &
!$omp					 private(cidp) &
!$omp					 private(cidp0, cidp1, cidp2, cidp3, cidp4, cidp5) &
!$omp					 private(cidw) &
!$omp					 private(cidw0, cidw1, cidw2, cidw3, cidw4, cidw5) &
!$omp					 private(cide) &
!$omp					 private(cide0, cide1, cide2, cide3, cide4, cide5) &
!$omp					 private(cids) &
!$omp					 private(cids0, cids1, cids2, cids3, cids4, cids5) &
!$omp					 private(cidn) &
!$omp					 private(cidn0, cidn1, cidn2, cidn3, cidn4, cidn5) &
!$omp					 private(cidb) &
!$omp					 private(cidb0, cidb1, cidb2, cidb3, cidb4, cidb5) &
!$omp					 private(cidt) &
!$omp					 private(cidt0, cidt1, cidt2, cidt3, cidt4, cidt5) 
!$omp do 
#else
#endif
  do k=1, kx
  do j=1, jx
  do i=1, ix
		cidp0 = cid0(i, j, k)
		cidp1 = cid1(i, j, k)
		cidp2 = cid2(i, j, k)
		cidp3 = cid3(i, j, k)
		cidp4 = cid4(i, j, k)
		cidp5 = cid5(i, j, k)

		cidw0 = cid0(i-1, j, k)
		cidw1 = cid1(i-1, j, k)
		cidw2 = cid2(i-1, j, k)
		cidw3 = cid3(i-1, j, k)
		cidw4 = cid4(i-1, j, k)
		cidw5 = cid5(i-1, j, k)

		cide0 = cid0(i+1, j, k)
		cide1 = cid1(i+1, j, k)
		cide2 = cid2(i+1, j, k)
		cide3 = cid3(i+1, j, k)
		cide4 = cid4(i+1, j, k)
		cide5 = cid5(i+1, j, k)

		cids0 = cid0(i, j-1, k)
		cids1 = cid1(i, j-1, k)
		cids2 = cid2(i, j-1, k)
		cids3 = cid3(i, j-1, k)
		cids4 = cid4(i, j-1, k)
		cids5 = cid5(i, j-1, k)

		cidn0 = cid0(i, j+1, k)
		cidn1 = cid1(i, j+1, k)
		cidn2 = cid2(i, j+1, k)
		cidn3 = cid3(i, j+1, k)
		cidn4 = cid4(i, j+1, k)
		cidn5 = cid5(i, j+1, k)

		cidb0 = cid0(i, j, k-1)
		cidb1 = cid1(i, j, k-1)
		cidb2 = cid2(i, j, k-1)
		cidb3 = cid3(i, j, k-1)
		cidb4 = cid4(i, j, k-1)
		cidb5 = cid5(i, j, k-1)

		cidt0 = cid0(i, j, k+1)
		cidt1 = cid1(i, j, k+1)
		cidt2 = cid2(i, j, k+1)
		cidt3 = cid3(i, j, k+1)
		cidt4 = cid4(i, j, k+1)
		cidt5 = cid5(i, j, k+1)

		if( cidp0 /= 0 ) then
			c0(i  , j, k) = 0.5
			c1(i-1, j, k) = 0.5
			cidw1 = cidp0
		endif
		if( cidp1 /= 0 ) then
			c1(i  , j, k) = 0.5
			c0(i+1, j, k) = 0.5
			cide0 = cidp1
		endif
		if( cidp2 /= 0 ) then
			c2(i, j  , k) = 0.5
			c3(i, j-1, k) = 0.5
			cids3 = cidp2
		endif
		if( cidp3 /= 0 ) then
			c3(i, j  , k) = 0.5
			c2(i, j+1, k) = 0.5
			cidn2 = cidp3
		endif
		if( cidp4 /= 0 ) then
			c4(i, j, k  ) = 0.5
			c5(i, j, k-1) = 0.5
			cidb5 = cidp4
		endif
		if( cidp5 /= 0 ) then
			c5(i, j, k  ) = 0.5
			c4(i, j, k+1) = 0.5
			cidt4 = cidp5
		endif

		cid0(i, j, k) = cidp0
		cid1(i, j, k) = cidp1
		cid2(i, j, k) = cidp2
		cid3(i, j, k) = cidp3
		cid4(i, j, k) = cidp4
		cid5(i, j, k) = cidp5

		cid0(i-1, j, k) = cidw0
		cid1(i-1, j, k) = cidw1
		cid2(i-1, j, k) = cidw2
		cid3(i-1, j, k) = cidw3
		cid4(i-1, j, k) = cidw4
		cid5(i-1, j, k) = cidw5

		cid0(i+1, j, k) = cide0
		cid1(i+1, j, k) = cide1
		cid2(i+1, j, k) = cide2
		cid3(i+1, j, k) = cide3
		cid4(i+1, j, k) = cide4
		cid5(i+1, j, k) = cide5

		cid0(i, j-1, k) = cids0
		cid1(i, j-1, k) = cids1
		cid2(i, j-1, k) = cids2
		cid3(i, j-1, k) = cids3
		cid4(i, j-1, k) = cids4
		cid5(i, j-1, k) = cids5

		cid0(i, j+1, k) = cidn0
		cid1(i, j+1, k) = cidn1
		cid2(i, j+1, k) = cidn2
		cid3(i, j+1, k) = cidn3
		cid4(i, j+1, k) = cidn4
		cid5(i, j+1, k) = cidn5

		cid0(i, j, k-1) = cidb0
		cid1(i, j, k-1) = cidb1
		cid2(i, j, k-1) = cidb2
		cid3(i, j, k-1) = cidb3
		cid4(i, j, k-1) = cidb4
		cid5(i, j, k-1) = cidb5

		cid0(i, j, k+1) = cidt0
		cid1(i, j, k+1) = cidt1
		cid2(i, j, k+1) = cidt2
		cid3(i, j, k+1) = cidt3
		cid4(i, j, k+1) = cidt4
		cid5(i, j, k+1) = cidt5
  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bstl_voxelize

subroutine bstl_cutoff( &
								c0, c1, c2, c3, c4, c5, &
								cid0, cid1, cid2, cid3, cid4, cid5, &
								eps, &
								sz, g)
  implicit none
  integer, dimension(3)   :: sz
  integer                 :: g
  integer                 :: i, j, k
  integer                 :: ix, jx, kx
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: c0, c1, c2, c3, c4, c5
	integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: cid0, cid1, cid2, cid3, cid4, cid5
	real										:: eps
	integer									:: cidp
	integer									:: cidp0, cidp1, cidp2, cidp3, cidp4, cidp5
	integer									:: cidw
	integer									:: cidw0, cidw1, cidw2, cidw3, cidw4, cidw5
	integer									:: cide
	integer									:: cide0, cide1, cide2, cide3, cide4, cide5
	integer									:: cids
	integer									:: cids0, cids1, cids2, cids3, cids4, cids5
	integer									:: cidn
	integer									:: cidn0, cidn1, cidn2, cidn3, cidn4, cidn5
	integer									:: cidb
	integer									:: cidb0, cidb1, cidb2, cidb3, cidb4, cidb5
	integer									:: cidt
	integer									:: cidt0, cidt1, cidt2, cidt3, cidt4, cidt5
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) &
!$omp					 private(cidp) &
!$omp					 private(cidp0, cidp1, cidp2, cidp3, cidp4, cidp5) &
!$omp					 private(cidw) &
!$omp					 private(cidw0, cidw1, cidw2, cidw3, cidw4, cidw5) &
!$omp					 private(cide) &
!$omp					 private(cide0, cide1, cide2, cide3, cide4, cide5) &
!$omp					 private(cids) &
!$omp					 private(cids0, cids1, cids2, cids3, cids4, cids5) &
!$omp					 private(cidn) &
!$omp					 private(cidn0, cidn1, cidn2, cidn3, cidn4, cidn5) &
!$omp					 private(cidb) &
!$omp					 private(cidb0, cidb1, cidb2, cidb3, cidb4, cidb5) &
!$omp					 private(cidt) &
!$omp					 private(cidt0, cidt1, cidt2, cidt3, cidt4, cidt5) 
!$omp do 
#else
#endif
  do k=1, kx
  do j=1, jx
  do i=1, ix
    c0(i, j, k) = max(eps, c0(i, j, k))
    c1(i, j, k) = max(eps, c1(i, j, k))
    c2(i, j, k) = max(eps, c2(i, j, k))
    c3(i, j, k) = max(eps, c3(i, j, k))
    c4(i, j, k) = max(eps, c4(i, j, k))
    c5(i, j, k) = max(eps, c5(i, j, k))
  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bstl_cutoff

subroutine bstl_symmetrize( &
								c0, c1, c2, c3, c4, c5, &
								cid0, cid1, cid2, cid3, cid4, cid5, &
								sz, g)
  implicit none
  integer, dimension(3)   :: sz
  integer                 :: g
  integer                 :: i, j, k
  integer                 :: ix, jx, kx
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: c0, c1, c2, c3, c4, c5
	integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: cid0, cid1, cid2, cid3, cid4, cid5
	integer									:: cidp
	integer									:: cidp0, cidp1, cidp2, cidp3, cidp4, cidp5
	integer									:: cidw
	integer									:: cidw0, cidw1, cidw2, cidw3, cidw4, cidw5
	integer									:: cide
	integer									:: cide0, cide1, cide2, cide3, cide4, cide5
	integer									:: cids
	integer									:: cids0, cids1, cids2, cids3, cids4, cids5
	integer									:: cidn
	integer									:: cidn0, cidn1, cidn2, cidn3, cidn4, cidn5
	integer									:: cidb
	integer									:: cidb0, cidb1, cidb2, cidb3, cidb4, cidb5
	integer									:: cidt
	integer									:: cidt0, cidt1, cidt2, cidt3, cidt4, cidt5
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) &
!$omp					 private(cidp) &
!$omp					 private(cidp0, cidp1, cidp2, cidp3, cidp4, cidp5) &
!$omp					 private(cidw) &
!$omp					 private(cidw0, cidw1, cidw2, cidw3, cidw4, cidw5) &
!$omp					 private(cide) &
!$omp					 private(cide0, cide1, cide2, cide3, cide4, cide5) &
!$omp					 private(cids) &
!$omp					 private(cids0, cids1, cids2, cids3, cids4, cids5) &
!$omp					 private(cidn) &
!$omp					 private(cidn0, cidn1, cidn2, cidn3, cidn4, cidn5) &
!$omp					 private(cidb) &
!$omp					 private(cidb0, cidb1, cidb2, cidb3, cidb4, cidb5) &
!$omp					 private(cidt) &
!$omp					 private(cidt0, cidt1, cidt2, cidt3, cidt4, cidt5) 
!$omp do 
#else
#endif
  do k=1, kx
  do j=1, jx
  do i=1, ix
		cidp0 = cid0(i, j, k)
		cidp1 = cid1(i, j, k)
		cidp2 = cid2(i, j, k)
		cidp3 = cid3(i, j, k)
		cidp4 = cid4(i, j, k)
		cidp5 = cid5(i, j, k)

		cidw0 = cid0(i-1, j, k)
		cidw1 = cid1(i-1, j, k)
		cidw2 = cid2(i-1, j, k)
		cidw3 = cid3(i-1, j, k)
		cidw4 = cid4(i-1, j, k)
		cidw5 = cid5(i-1, j, k)

		cide0 = cid0(i+1, j, k)
		cide1 = cid1(i+1, j, k)
		cide2 = cid2(i+1, j, k)
		cide3 = cid3(i+1, j, k)
		cide4 = cid4(i+1, j, k)
		cide5 = cid5(i+1, j, k)

		cids0 = cid0(i, j-1, k)
		cids1 = cid1(i, j-1, k)
		cids2 = cid2(i, j-1, k)
		cids3 = cid3(i, j-1, k)
		cids4 = cid4(i, j-1, k)
		cids5 = cid5(i, j-1, k)

		cidn0 = cid0(i, j+1, k)
		cidn1 = cid1(i, j+1, k)
		cidn2 = cid2(i, j+1, k)
		cidn3 = cid3(i, j+1, k)
		cidn4 = cid4(i, j+1, k)
		cidn5 = cid5(i, j+1, k)

		cidb0 = cid0(i, j, k-1)
		cidb1 = cid1(i, j, k-1)
		cidb2 = cid2(i, j, k-1)
		cidb3 = cid3(i, j, k-1)
		cidb4 = cid4(i, j, k-1)
		cidb5 = cid5(i, j, k-1)

		cidt0 = cid0(i, j, k+1)
		cidt1 = cid1(i, j, k+1)
		cidt2 = cid2(i, j, k+1)
		cidt3 = cid3(i, j, k+1)
		cidt4 = cid4(i, j, k+1)
		cidt5 = cid5(i, j, k+1)


		if( cidp0 /= 0 ) then
			c1(i-1, j, k) = 1.0d0 - c0(i  , j, k) 
			cidw1 = cidp0
		endif
		if( cidp1 /= 0 ) then
			c0(i+1, j, k) = 1.0d0 - c1(i  , j, k) 
			cide0 = cidp1
		endif
		if( cidp2 /= 0 ) then
			c3(i, j-1, k) = 1.0d0 - c2(i, j  , k) 
			cids3 = cidp2
		endif
		if( cidp3 /= 0 ) then
			c2(i, j+1, k) = 1.0d0 - c3(i, j  , k) 
			cidn2 = cidp3
		endif
		if( cidp4 /= 0 ) then
			c5(i, j, k-1) = 1.0d0 - c4(i, j, k  ) 
			cidb5 = cidp4
		endif
		if( cidp5 /= 0 ) then
			c4(i, j, k+1) = 1.0d0 - c5(i, j, k  ) 
			cidt4 = cidp5
		endif

		cid0(i, j, k) = cidp0
		cid1(i, j, k) = cidp1
		cid2(i, j, k) = cidp2
		cid3(i, j, k) = cidp3
		cid4(i, j, k) = cidp4
		cid5(i, j, k) = cidp5

		cid0(i-1, j, k) = cidw0
		cid1(i-1, j, k) = cidw1
		cid2(i-1, j, k) = cidw2
		cid3(i-1, j, k) = cidw3
		cid4(i-1, j, k) = cidw4
		cid5(i-1, j, k) = cidw5

		cid0(i+1, j, k) = cide0
		cid1(i+1, j, k) = cide1
		cid2(i+1, j, k) = cide2
		cid3(i+1, j, k) = cide3
		cid4(i+1, j, k) = cide4
		cid5(i+1, j, k) = cide5

		cid0(i, j-1, k) = cids0
		cid1(i, j-1, k) = cids1
		cid2(i, j-1, k) = cids2
		cid3(i, j-1, k) = cids3
		cid4(i, j-1, k) = cids4
		cid5(i, j-1, k) = cids5

		cid0(i, j+1, k) = cidn0
		cid1(i, j+1, k) = cidn1
		cid2(i, j+1, k) = cidn2
		cid3(i, j+1, k) = cidn3
		cid4(i, j+1, k) = cidn4
		cid5(i, j+1, k) = cidn5

		cid0(i, j, k-1) = cidb0
		cid1(i, j, k-1) = cidb1
		cid2(i, j, k-1) = cidb2
		cid3(i, j, k-1) = cidb3
		cid4(i, j, k-1) = cidb4
		cid5(i, j, k-1) = cidb5

		cid0(i, j, k+1) = cidt0
		cid1(i, j, k+1) = cidt1
		cid2(i, j, k+1) = cidt2
		cid3(i, j, k+1) = cidt3
		cid4(i, j, k+1) = cidt4
		cid5(i, j, k+1) = cidt5
  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bstl_symmetrize

subroutine bstl_detect_zerocut( &
								c0, c1, c2, c3, c4, c5, &
								cid0, cid1, cid2, cid3, cid4, cid5, &
                n, &
								sz, g)
  implicit none
  integer, dimension(3)   :: sz
  integer                 :: g
  integer                 :: i, j, k
  integer                 :: ix, jx, kx
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: c0, c1, c2, c3, c4, c5
	integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: cid0, cid1, cid2, cid3, cid4, cid5
	integer									:: cidp
	integer									:: cidp0, cidp1, cidp2, cidp3, cidp4, cidp5
	integer									:: cidw
	integer									:: cidw0, cidw1, cidw2, cidw3, cidw4, cidw5
	integer									:: cide
	integer									:: cide0, cide1, cide2, cide3, cide4, cide5
	integer									:: cids
	integer									:: cids0, cids1, cids2, cids3, cids4, cids5
	integer									:: cidn
	integer									:: cidn0, cidn1, cidn2, cidn3, cidn4, cidn5
	integer									:: cidb
	integer									:: cidb0, cidb1, cidb2, cidb3, cidb4, cidb5
	integer									:: cidt
	integer									:: cidt0, cidt1, cidt2, cidt3, cidt4, cidt5
	integer									:: n
	real										:: cut_min
	integer									:: cid_min
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) &
!$omp					 private(cidp) &
!$omp					 private(cidp0, cidp1, cidp2, cidp3, cidp4, cidp5) &
!$omp					 private(cidw) &
!$omp					 private(cidw0, cidw1, cidw2, cidw3, cidw4, cidw5) &
!$omp					 private(cide) &
!$omp					 private(cide0, cide1, cide2, cide3, cide4, cide5) &
!$omp					 private(cids) &
!$omp					 private(cids0, cids1, cids2, cids3, cids4, cids5) &
!$omp					 private(cidn) &
!$omp					 private(cidn0, cidn1, cidn2, cidn3, cidn4, cidn5) &
!$omp					 private(cidb) &
!$omp					 private(cidb0, cidb1, cidb2, cidb3, cidb4, cidb5) &
!$omp					 private(cidt) &
!$omp					 private(cidt0, cidt1, cidt2, cidt3, cidt4, cidt5) &
!$omp					 private(cut_min) &
!$omp					 private(cid_min) 
!$omp do reduction(+: n)
#else
#endif
  do k=1, kx
  do j=1, jx
  do i=1, ix
		cidp0 = cid0(i, j, k)
		cidp1 = cid1(i, j, k)
		cidp2 = cid2(i, j, k)
		cidp3 = cid3(i, j, k)
		cidp4 = cid4(i, j, k)
		cidp5 = cid5(i, j, k)

		cidw0 = cid0(i-1, j, k)
		cidw1 = cid1(i-1, j, k)
		cidw2 = cid2(i-1, j, k)
		cidw3 = cid3(i-1, j, k)
		cidw4 = cid4(i-1, j, k)
		cidw5 = cid5(i-1, j, k)

		cide0 = cid0(i+1, j, k)
		cide1 = cid1(i+1, j, k)
		cide2 = cid2(i+1, j, k)
		cide3 = cid3(i+1, j, k)
		cide4 = cid4(i+1, j, k)
		cide5 = cid5(i+1, j, k)

		cids0 = cid0(i, j-1, k)
		cids1 = cid1(i, j-1, k)
		cids2 = cid2(i, j-1, k)
		cids3 = cid3(i, j-1, k)
		cids4 = cid4(i, j-1, k)
		cids5 = cid5(i, j-1, k)

		cidn0 = cid0(i, j+1, k)
		cidn1 = cid1(i, j+1, k)
		cidn2 = cid2(i, j+1, k)
		cidn3 = cid3(i, j+1, k)
		cidn4 = cid4(i, j+1, k)
		cidn5 = cid5(i, j+1, k)

		cidb0 = cid0(i, j, k-1)
		cidb1 = cid1(i, j, k-1)
		cidb2 = cid2(i, j, k-1)
		cidb3 = cid3(i, j, k-1)
		cidb4 = cid4(i, j, k-1)
		cidb5 = cid5(i, j, k-1)

		cidt0 = cid0(i, j, k+1)
		cidt1 = cid1(i, j, k+1)
		cidt2 = cid2(i, j, k+1)
		cidt3 = cid3(i, j, k+1)
		cidt4 = cid4(i, j, k+1)
		cidt5 = cid5(i, j, k+1)


		cut_min = 1.0d0	
		cid_min = 0
		if( cut_min > c0(i, j, k) ) then
			cut_min = c0(i, j, k)
			cid_min = cidp0
		endif
		if( cut_min > c1(i, j, k) ) then
			cut_min = c1(i, j, k)
			cid_min = cidp1
		endif
		if( cut_min > c2(i, j, k) ) then
			cut_min = c2(i, j, k)
			cid_min = cidp2
		endif
		if( cut_min > c3(i, j, k) ) then
			cut_min = c3(i, j, k)
			cid_min = cidp3
		endif
		if( cut_min > c4(i, j, k) ) then
			cut_min = c4(i, j, k)
			cid_min = cidp4
		endif
		if( cut_min > c5(i, j, k) ) then
			cut_min = c5(i, j, k)
			cid_min = cidp5
		endif

		if( cut_min == 0.0d0 ) then
			c0(i, j, k) = 1.0d0
			c1(i, j, k) = 1.0d0
			c2(i, j, k) = 1.0d0
			c3(i, j, k) = 1.0d0
			c4(i, j, k) = 1.0d0
			c5(i, j, k) = 1.0d0
      c1(i-1, j, k) = 1.0d0

      cidw1 = cid_min
      c0(i+1, j, k) = 1.0d0
      cide0 = cid_min
      c3(i, j-1, k) = 1.0d0
      cids3 = cid_min
      c2(i, j+1, k) = 1.0d0
      cidn2 = cid_min
      c5(i, j, k-1) = 1.0d0
      cidb5 = cid_min
      c4(i, j, k+1) = 1.0d0
      cidt4 = cid_min

			n = n + 1
		endif

		cid0(i, j, k) = cidp0
		cid1(i, j, k) = cidp1
		cid2(i, j, k) = cidp2
		cid3(i, j, k) = cidp3
		cid4(i, j, k) = cidp4
		cid5(i, j, k) = cidp5

		cid0(i-1, j, k) = cidw0
		cid1(i-1, j, k) = cidw1
		cid2(i-1, j, k) = cidw2
		cid3(i-1, j, k) = cidw3
		cid4(i-1, j, k) = cidw4
		cid5(i-1, j, k) = cidw5

		cid0(i+1, j, k) = cide0
		cid1(i+1, j, k) = cide1
		cid2(i+1, j, k) = cide2
		cid3(i+1, j, k) = cide3
		cid4(i+1, j, k) = cide4
		cid5(i+1, j, k) = cide5

		cid0(i, j-1, k) = cids0
		cid1(i, j-1, k) = cids1
		cid2(i, j-1, k) = cids2
		cid3(i, j-1, k) = cids3
		cid4(i, j-1, k) = cids4
		cid5(i, j-1, k) = cids5

		cid0(i, j+1, k) = cidn0
		cid1(i, j+1, k) = cidn1
		cid2(i, j+1, k) = cidn2
		cid3(i, j+1, k) = cidn3
		cid4(i, j+1, k) = cidn4
		cid5(i, j+1, k) = cidn5

		cid0(i, j, k-1) = cidb0
		cid1(i, j, k-1) = cidb1
		cid2(i, j, k-1) = cidb2
		cid3(i, j, k-1) = cidb3
		cid4(i, j, k-1) = cidb4
		cid5(i, j, k-1) = cidb5

		cid0(i, j, k+1) = cidt0
		cid1(i, j, k+1) = cidt1
		cid2(i, j, k+1) = cidt2
		cid3(i, j, k+1) = cidt3
		cid4(i, j, k+1) = cidt4
		cid5(i, j, k+1) = cidt5
  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bstl_detect_zerocut

subroutine bstl_fill_holes( &
								c0, c1, c2, c3, c4, c5, &
								cid0, cid1, cid2, cid3, cid4, cid5, &
								n, &
								bClose, &
								sz, g)
  implicit none
  integer, dimension(3)   :: sz
  integer                 :: g
  integer                 :: i, j, k
  integer                 :: ix, jx, kx
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: c0, c1, c2, c3, c4, c5
	integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: cid0, cid1, cid2, cid3, cid4, cid5
	integer									:: n
	integer									:: bClose
	integer									:: cidp
	integer									:: cidp0, cidp1, cidp2, cidp3, cidp4, cidp5
	integer									:: cidw
	integer									:: cidw0, cidw1, cidw2, cidw3, cidw4, cidw5
	integer									:: cide
	integer									:: cide0, cide1, cide2, cide3, cide4, cide5
	integer									:: cids
	integer									:: cids0, cids1, cids2, cids3, cids4, cids5
	integer									:: cidn
	integer									:: cidn0, cidn1, cidn2, cidn3, cidn4, cidn5
	integer									:: cidb
	integer									:: cidb0, cidb1, cidb2, cidb3, cidb4, cidb5
	integer									:: cidt
	integer									:: cidt0, cidt1, cidt2, cidt3, cidt4, cidt5
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
	n = 0
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) &
!$omp					 private(cidp) &
!$omp					 private(cidp0, cidp1, cidp2, cidp3, cidp4, cidp5) &
!$omp					 private(cidw) &
!$omp					 private(cidw0, cidw1, cidw2, cidw3, cidw4, cidw5) &
!$omp					 private(cide) &
!$omp					 private(cide0, cide1, cide2, cide3, cide4, cide5) &
!$omp					 private(cids) &
!$omp					 private(cids0, cids1, cids2, cids3, cids4, cids5) &
!$omp					 private(cidn) &
!$omp					 private(cidn0, cidn1, cidn2, cidn3, cidn4, cidn5) &
!$omp					 private(cidb) &
!$omp					 private(cidb0, cidb1, cidb2, cidb3, cidb4, cidb5) &
!$omp					 private(cidt) &
!$omp					 private(cidt0, cidt1, cidt2, cidt3, cidt4, cidt5)
!$omp do reduction(+: n)
#else
#endif
  do k=1, kx
  do j=1, jx
  do i=1, ix
		cidp0 = cid0(i, j, k)
		cidp1 = cid1(i, j, k)
		cidp2 = cid2(i, j, k)
		cidp3 = cid3(i, j, k)
		cidp4 = cid4(i, j, k)
		cidp5 = cid5(i, j, k)

		cidw0 = cid0(i-1, j, k)
		cidw1 = cid1(i-1, j, k)
		cidw2 = cid2(i-1, j, k)
		cidw3 = cid3(i-1, j, k)
		cidw4 = cid4(i-1, j, k)
		cidw5 = cid5(i-1, j, k)

		cide0 = cid0(i+1, j, k)
		cide1 = cid1(i+1, j, k)
		cide2 = cid2(i+1, j, k)
		cide3 = cid3(i+1, j, k)
		cide4 = cid4(i+1, j, k)
		cide5 = cid5(i+1, j, k)

		cids0 = cid0(i, j-1, k)
		cids1 = cid1(i, j-1, k)
		cids2 = cid2(i, j-1, k)
		cids3 = cid3(i, j-1, k)
		cids4 = cid4(i, j-1, k)
		cids5 = cid5(i, j-1, k)

		cidn0 = cid0(i, j+1, k)
		cidn1 = cid1(i, j+1, k)
		cidn2 = cid2(i, j+1, k)
		cidn3 = cid3(i, j+1, k)
		cidn4 = cid4(i, j+1, k)
		cidn5 = cid5(i, j+1, k)

		cidb0 = cid0(i, j, k-1)
		cidb1 = cid1(i, j, k-1)
		cidb2 = cid2(i, j, k-1)
		cidb3 = cid3(i, j, k-1)
		cidb4 = cid4(i, j, k-1)
		cidb5 = cid5(i, j, k-1)

		cidt0 = cid0(i, j, k+1)
		cidt1 = cid1(i, j, k+1)
		cidt2 = cid2(i, j, k+1)
		cidt3 = cid3(i, j, k+1)
		cidt4 = cid4(i, j, k+1)
		cidt5 = cid5(i, j, k+1)

		if( cidp0 == 0 ) then
			if( (cids0 /= 0 .or. cidp2 /=0 .or. cidw2 /=0) .and. &
					(cidn0 /= 0 .or. cidp3 /=0 .or. cidw3 /=0) .and. &
					(cidb0 /= 0 .or. cidp4 /=0 .or. cidw4 /=0) .and. &
					(cidt0 /= 0 .or. cidp5 /=0 .or. cidw5 /=0) ) then
				if( bClose == 1 ) then
					cidp0 = 1
					cidw1 = 1
					c0(i  , j, k) = 0.5
					c1(i-1, j, k) = 0.5
				endif
				n = n + 1
			endif
		endif

		if( cidp1 == 0 ) then
			if( (cids1 /= 0 .or. cidp2 /=0 .or. cide2 /=0) .and. &
					(cidn1 /= 0 .or. cidp3 /=0 .or. cide3 /=0) .and. &
					(cidb1 /= 0 .or. cidp4 /=0 .or. cide4 /=0) .and. &
					(cidt1 /= 0 .or. cidp5 /=0 .or. cide5 /=0) ) then
				if( bClose == 1 ) then
					cidp1 = 1
					cide0 = 1
					c1(i  , j, k) = 0.5
					c0(i+1, j, k) = 0.5
				endif
				n = n + 1
			endif
		endif

		if( cidp2 == 0 ) then
			if( (cidw2 /= 0 .or. cidp0 /=0 .or. cids0 /=0) .and. &
					(cide2 /= 0 .or. cidp1 /=0 .or. cids1 /=0) .and. &
					(cidb2 /= 0 .or. cidp4 /=0 .or. cids4 /=0) .and. &
					(cidt2 /= 0 .or. cidp5 /=0 .or. cids5 /=0) ) then
				if( bClose == 1 ) then
					cidp2 = 1
					cids3 = 1
					c2(i, j  , k) = 0.5
					c3(i, j-1, k) = 0.5
				endif
				n = n + 1
			endif
		endif

		if( cidp3 == 0 ) then
			if( (cidw3 /= 0 .or. cidp0 /=0 .or. cidn0 /=0) .and. &
					(cide3 /= 0 .or. cidp1 /=0 .or. cidn1 /=0) .and. &
					(cidb3 /= 0 .or. cidp4 /=0 .or. cidn4 /=0) .and. &
					(cidt3 /= 0 .or. cidp5 /=0 .or. cidn5 /=0) ) then
				if( bClose == 1 ) then
					cidp3 = 1
					cidn2 = 1
					c3(i, j  , k) = 0.5
					c2(i, j+1, k) = 0.5
				endif
				n = n + 1
			endif
		endif

		if( cidp4 == 0 ) then
			if( (cidw4 /= 0 .or. cidp0 /=0 .or. cidb0 /=0) .and. &
					(cide4 /= 0 .or. cidp1 /=0 .or. cidb1 /=0) .and. &
					(cids4 /= 0 .or. cidp2 /=0 .or. cidb2 /=0) .and. &
					(cidn4 /= 0 .or. cidp3 /=0 .or. cidb3 /=0) ) then
				if( bClose == 1 ) then
					cidp4 = 1
					cidb5 = 1
					c4(i, j, k  ) = 0.5
					c5(i, j, k-1) = 0.5
				endif
				n = n + 1
			endif
		endif

		if( cidp5 == 0 ) then
			if( (cidw5 /= 0 .or. cidp0 /=0 .or. cidt0 /=0) .and. &
					(cide5 /= 0 .or. cidp1 /=0 .or. cidt1 /=0) .and. &
					(cids5 /= 0 .or. cidp2 /=0 .or. cidt2 /=0) .and. &
					(cidn5 /= 0 .or. cidp3 /=0 .or. cidt3 /=0) ) then
				if( bClose == 1 ) then
					cidp5 = 1
					cidt4 = 1
					c5(i, j, k  ) = 0.5
					c4(i, j, k+1) = 0.5
				endif
				n = n + 1
			endif
		endif

		cid0(i, j, k) = cidp0
		cid1(i, j, k) = cidp1
		cid2(i, j, k) = cidp2
		cid3(i, j, k) = cidp3
		cid4(i, j, k) = cidp4
		cid5(i, j, k) = cidp5

		cid0(i-1, j, k) = cidw0
		cid1(i-1, j, k) = cidw1
		cid2(i-1, j, k) = cidw2
		cid3(i-1, j, k) = cidw3
		cid4(i-1, j, k) = cidw4
		cid5(i-1, j, k) = cidw5

		cid0(i+1, j, k) = cide0
		cid1(i+1, j, k) = cide1
		cid2(i+1, j, k) = cide2
		cid3(i+1, j, k) = cide3
		cid4(i+1, j, k) = cide4
		cid5(i+1, j, k) = cide5

		cid0(i, j-1, k) = cids0
		cid1(i, j-1, k) = cids1
		cid2(i, j-1, k) = cids2
		cid3(i, j-1, k) = cids3
		cid4(i, j-1, k) = cids4
		cid5(i, j-1, k) = cids5

		cid0(i, j+1, k) = cidn0
		cid1(i, j+1, k) = cidn1
		cid2(i, j+1, k) = cidn2
		cid3(i, j+1, k) = cidn3
		cid4(i, j+1, k) = cidn4
		cid5(i, j+1, k) = cidn5

		cid0(i, j, k-1) = cidb0
		cid1(i, j, k-1) = cidb1
		cid2(i, j, k-1) = cidb2
		cid3(i, j, k-1) = cidb3
		cid4(i, j, k-1) = cidb4
		cid5(i, j, k-1) = cidb5

		cid0(i, j, k+1) = cidt0
		cid1(i, j, k+1) = cidt1
		cid2(i, j, k+1) = cidt2
		cid3(i, j, k+1) = cidt3
		cid4(i, j, k+1) = cidt4
		cid5(i, j, k+1) = cidt5
  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bstl_fill_holes

subroutine bstl_fill_holes_v2( &
								c0, c1, c2, c3, c4, c5, &
								cid0, cid1, cid2, cid3, cid4, cid5, &
								n, &
								bClose, &
								sz, g)
  implicit none
  integer, dimension(3)   :: sz
  integer                 :: g
  integer                 :: i, j, k
  integer                 :: ix, jx, kx
	real, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: c0, c1, c2, c3, c4, c5
	integer, dimension(1-g:sz(1)+g, 1-g:sz(2)+g, 1-g:sz(3)+g)  :: cid0, cid1, cid2, cid3, cid4, cid5
	integer									:: n
	integer									:: bClose
	integer									:: cidp
	integer									:: cidp0, cidp1, cidp2, cidp3, cidp4, cidp5
	integer									:: cidw
	integer									:: cidw0, cidw1, cidw2, cidw3, cidw4, cidw5
	integer									:: cide
	integer									:: cide0, cide1, cide2, cide3, cide4, cide5
	integer									:: cids
	integer									:: cids0, cids1, cids2, cids3, cids4, cids5
	integer									:: cidn
	integer									:: cidn0, cidn1, cidn2, cidn3, cidn4, cidn5
	integer									:: cidb
	integer									:: cidb0, cidb1, cidb2, cidb3, cidb4, cidb5
	integer									:: cidt
	integer									:: cidt0, cidt1, cidt2, cidt3, cidt4, cidt5
  ix = sz(1)
  jx = sz(2)
  kx = sz(3)
	n = 0
#ifdef _BLOCK_IS_LARGE_
!$omp parallel private(i, j, k) &
!$omp					 private(cidp) &
!$omp					 private(cidp0, cidp1, cidp2, cidp3, cidp4, cidp5) &
!$omp					 private(cidw) &
!$omp					 private(cidw0, cidw1, cidw2, cidw3, cidw4, cidw5) &
!$omp					 private(cide) &
!$omp					 private(cide0, cide1, cide2, cide3, cide4, cide5) &
!$omp					 private(cids) &
!$omp					 private(cids0, cids1, cids2, cids3, cids4, cids5) &
!$omp					 private(cidn) &
!$omp					 private(cidn0, cidn1, cidn2, cidn3, cidn4, cidn5) &
!$omp					 private(cidb) &
!$omp					 private(cidb0, cidb1, cidb2, cidb3, cidb4, cidb5) &
!$omp					 private(cidt) &
!$omp					 private(cidt0, cidt1, cidt2, cidt3, cidt4, cidt5)
!$omp do reduction(+: n)
#else
#endif
  do k=1, kx
  do j=1, jx
  do i=1, ix
		cidp0 = cid0(i, j, k)
		cidp1 = cid1(i, j, k)
		cidp2 = cid2(i, j, k)
		cidp3 = cid3(i, j, k)
		cidp4 = cid4(i, j, k)
		cidp5 = cid5(i, j, k)

		cidw0 = cid0(i-1, j, k)
		cidw1 = cid1(i-1, j, k)
		cidw2 = cid2(i-1, j, k)
		cidw3 = cid3(i-1, j, k)
		cidw4 = cid4(i-1, j, k)
		cidw5 = cid5(i-1, j, k)

		cide0 = cid0(i+1, j, k)
		cide1 = cid1(i+1, j, k)
		cide2 = cid2(i+1, j, k)
		cide3 = cid3(i+1, j, k)
		cide4 = cid4(i+1, j, k)
		cide5 = cid5(i+1, j, k)

		cids0 = cid0(i, j-1, k)
		cids1 = cid1(i, j-1, k)
		cids2 = cid2(i, j-1, k)
		cids3 = cid3(i, j-1, k)
		cids4 = cid4(i, j-1, k)
		cids5 = cid5(i, j-1, k)

		cidn0 = cid0(i, j+1, k)
		cidn1 = cid1(i, j+1, k)
		cidn2 = cid2(i, j+1, k)
		cidn3 = cid3(i, j+1, k)
		cidn4 = cid4(i, j+1, k)
		cidn5 = cid5(i, j+1, k)

		cidb0 = cid0(i, j, k-1)
		cidb1 = cid1(i, j, k-1)
		cidb2 = cid2(i, j, k-1)
		cidb3 = cid3(i, j, k-1)
		cidb4 = cid4(i, j, k-1)
		cidb5 = cid5(i, j, k-1)

		cidt0 = cid0(i, j, k+1)
		cidt1 = cid1(i, j, k+1)
		cidt2 = cid2(i, j, k+1)
		cidt3 = cid3(i, j, k+1)
		cidt4 = cid4(i, j, k+1)
		cidt5 = cid5(i, j, k+1)

		if( cidp0 == 0 ) then
			if( ( (cids0 /= 0 .or. cidp2 /=0 .or. cidw2 /=0) .and. &
					  (cidn0 /= 0 .or. cidp3 /=0 .or. cidw3 /=0) ) .or. &
					( (cidb0 /= 0 .or. cidp4 /=0 .or. cidw4 /=0) .and. &
					  (cidt0 /= 0 .or. cidp5 /=0 .or. cidw5 /=0) ) ) then
				if( bClose == 1 ) then
					cidp0 = 1
					cidw1 = 1
					c0(i  , j, k) = 0.5
					c1(i-1, j, k) = 0.5
				end if
				n = n + 1
			endif
		endif

		if( cidp1 == 0 ) then
			if( ( (cids1 /= 0 .or. cidp2 /=0 .or. cide2 /=0) .and. &
					  (cidn1 /= 0 .or. cidp3 /=0 .or. cide3 /=0) ) .or. &
					( (cidb1 /= 0 .or. cidp4 /=0 .or. cide4 /=0) .and. &
					  (cidt1 /= 0 .or. cidp5 /=0 .or. cide5 /=0) ) ) then
				if( bClose == 1 ) then
					cidp1 = 1
					cide0 = 1
					c1(i  , j, k) = 0.5
					c0(i+1, j, k) = 0.5
				end if
				n = n + 1
			endif
		endif

		if( cidp2 == 0 ) then
			if( ( (cidw2 /= 0 .or. cidp0 /=0 .or. cids0 /=0) .and. &
					  (cide2 /= 0 .or. cidp1 /=0 .or. cids1 /=0) ) .or. &
					( (cidb2 /= 0 .or. cidp4 /=0 .or. cids4 /=0) .and. &
					  (cidt2 /= 0 .or. cidp5 /=0 .or. cids5 /=0) ) ) then
				if( bClose == 1 ) then
					cidp2 = 1
					cids3 = 1
					c2(i, j  , k) = 0.5
					c3(i, j-1, k) = 0.5
				end if
				n = n + 1
			endif
		endif

		if( cidp3 == 0 ) then
			if( ( (cidw3 /= 0 .or. cidp0 /=0 .or. cidn0 /=0) .and. &
					  (cide3 /= 0 .or. cidp1 /=0 .or. cidn1 /=0) ) .or. &
					( (cidb3 /= 0 .or. cidp4 /=0 .or. cidn4 /=0) .and. &
					  (cidt3 /= 0 .or. cidp5 /=0 .or. cidn5 /=0) ) ) then
				if( bClose == 1 ) then
					cidp3 = 1
					cidn2 = 1
					c3(i, j  , k) = 0.5
					c2(i, j+1, k) = 0.5
				end if
				n = n + 1
			endif
		endif

		if( cidp4 == 0 ) then
			if( ( (cidw4 /= 0 .or. cidp0 /=0 .or. cidb0 /=0) .and. &
					  (cide4 /= 0 .or. cidp1 /=0 .or. cidb1 /=0) ) .or. &
					( (cids4 /= 0 .or. cidp2 /=0 .or. cidb2 /=0) .and. &
					  (cidn4 /= 0 .or. cidp3 /=0 .or. cidb3 /=0) ) ) then
				if( bClose == 1 ) then
					cidp4 = 1
					cidb5 = 1
					c4(i, j, k  ) = 0.5
					c5(i, j, k-1) = 0.5
				end if
				n = n + 1
			endif
		endif

		if( cidp5 == 0 ) then
			if( ( (cidw5 /= 0 .or. cidp0 /=0 .or. cidt0 /=0) .and. &
					  (cide5 /= 0 .or. cidp1 /=0 .or. cidt1 /=0) ) .or. &
					( (cids5 /= 0 .or. cidp2 /=0 .or. cidt2 /=0) .and. &
					  (cidn5 /= 0 .or. cidp3 /=0 .or. cidt3 /=0) ) ) then
				if( bClose == 1 ) then
					cidp5 = 1
					cidt4 = 1
					c5(i, j, k  ) = 0.5
					c4(i, j, k+1) = 0.5
				end if
				n = n + 1
			endif
		endif

		cid0(i, j, k) = cidp0
		cid1(i, j, k) = cidp1
		cid2(i, j, k) = cidp2
		cid3(i, j, k) = cidp3
		cid4(i, j, k) = cidp4
		cid5(i, j, k) = cidp5

		cid0(i-1, j, k) = cidw0
		cid1(i-1, j, k) = cidw1
		cid2(i-1, j, k) = cidw2
		cid3(i-1, j, k) = cidw3
		cid4(i-1, j, k) = cidw4
		cid5(i-1, j, k) = cidw5

		cid0(i+1, j, k) = cide0
		cid1(i+1, j, k) = cide1
		cid2(i+1, j, k) = cide2
		cid3(i+1, j, k) = cide3
		cid4(i+1, j, k) = cide4
		cid5(i+1, j, k) = cide5

		cid0(i, j-1, k) = cids0
		cid1(i, j-1, k) = cids1
		cid2(i, j-1, k) = cids2
		cid3(i, j-1, k) = cids3
		cid4(i, j-1, k) = cids4
		cid5(i, j-1, k) = cids5

		cid0(i, j+1, k) = cidn0
		cid1(i, j+1, k) = cidn1
		cid2(i, j+1, k) = cidn2
		cid3(i, j+1, k) = cidn3
		cid4(i, j+1, k) = cidn4
		cid5(i, j+1, k) = cidn5

		cid0(i, j, k-1) = cidb0
		cid1(i, j, k-1) = cidb1
		cid2(i, j, k-1) = cidb2
		cid3(i, j, k-1) = cidb3
		cid4(i, j, k-1) = cidb4
		cid5(i, j, k-1) = cidb5

		cid0(i, j, k+1) = cidt0
		cid1(i, j, k+1) = cidt1
		cid2(i, j, k+1) = cidt2
		cid3(i, j, k+1) = cidt3
		cid4(i, j, k+1) = cidt4
		cid5(i, j, k+1) = cidt5
  end do
  end do
  end do
#ifdef _BLOCK_IS_LARGE_
!$omp end do
!$omp end parallel
#else
#endif
end subroutine bstl_fill_holes_v2

