program testr2c
implicit none
integer :: nx, ny, nz, i,j
!real, dimension(:,:,:), allocatable :: a, b
complex, dimension(:,:,:), allocatable :: a, b
real :: pi
nx = 16
ny = 1
nz = 1
allocate(a(nx,ny,nz),b(nx,ny,nz))
pi = 4.*atan(1.)
a = 0.
!a(1,1,1) = cmplx(1.0,0)
do i = 1, nx
   do j = 1, ny
!      a(i,j,1) = sin(2.*pi*float(i-1)/float(nx))
       a(i,j,1) = cmplx(sin(2.*pi*float(i-1)/float(nx)),cos(2.*pi*float(i-1)/float(nx)))
   enddo
enddo
b = a
!do i = 1, nx
!  print '(A7,1x,i2,4(1x,f12.5))','org: ',i,(a(i,j,1),j=1,ny)
!enddo
!print *,'org: ',a
!$acc data copy(a)
!$acc host_data use_device(a)
!call cud2z( nx, ny, nz, a )
call cuz2zf( nx, ny, nz, a )
!$acc end host_data
!$acc end data
a = a/nx
!do i = 1, nx
!  print '(A7,1x,i2,4(1x,f12.5))','trans: ',i,(a(i,j,1),j=1,ny)
!enddo
!stop
!$acc data copy(a)
!$acc host_data use_device(a)
!call cuz2d( nx, ny, nz, a )
call cuz2zb( nx, ny, nz, a )
!$acc end host_data
!$acc end data
do i = 1, nx
  print '(A6,1x,i2,4(1x,f12.5))','back: ',i,b(i,1,1),(a(i,j,1),j=1,ny)
enddo
end
