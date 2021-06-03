!-----------------------------------------------------------
! Test code to check the 1d-cufft:
!   1. Real to complex
!   2. Complex to real
!   3. complex to comples, backward and inverse
! several of instances of the samples together in batch
! 
! In all cases, transform is along X-direction and in place.
!-----------------------------------------------------------
program testr2c
implicit none
integer :: nx, ny, nz, i, j, k
real, dimension(:,:,:), allocatable :: a, ab
complex, dimension(:,:,:), allocatable :: ca, cab
real :: pi
nx = 16
ny = 1
nz = 1
allocate(a(nx+2,ny,nz),ab(nx,ny,nz),ca(nx,ny,nz),cab(nx,ny,nz))
pi = 4.*atan(1.)
a = 0.
ca = cmplx(0.,0.)
do k = 1, nz
   do j = 1, ny
      do i = 1, nx
          a(i,j,k) = sin(2.*pi*float(i-1)/float(nx))
          ca(i,j,k) = cmplx(sin(2.*pi*float(i-1)/float(nx)),cos(2.*pi*float(i-1)/float(nx)))
      enddo
   enddo
enddo
ab = a
cab = ca
!!$acc data copy(a,ca)
!!$acc host_data use_device(a,ca)
!call cud2z( nx, ny, nz, a )
!call fortcud2z( nx, 1, ny, 1, nz, (1./nx), a )
!call cuz2zf( nx, ny, nz, ca )
!!$acc end host_data
!!$acc end data
!a = a/nx
!ca = ca/nx
!!$acc data copy(a,ca)
!!$acc host_data use_device(a,ca)
!call cuz2d( nx, ny, nz, a )
call fortcuz2d( nx, 1, ny, 1, nz, a )
!call cuz2zb( nx, ny, nz, ca )
!!$acc end host_data
!!$acc end data
do i = 1, nx
  print '(A9,1x,i2,8(1x,f5.2))','r2dnback: ',i,ab(i,1,1),(a(i,j,1),j=1,ny)
enddo
do i = 1, nx
! print '(A9,1x,i2,8(1x,f5.2))','z2znback: ',i,cab(i,1,1),(ca(i,j,1),j=1,ny)
enddo
end
