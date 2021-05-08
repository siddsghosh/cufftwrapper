program testr2c
implicit none
integer :: nx, ny, nz, i
real, dimension(:,:,:), allocatable :: a, b
!complex, dimension(:,:,:), allocatable :: a, b
real :: pi
nx = 16
ny = 2
nz = 1
allocate(a(nx+2,ny,nz),b(nx+2,ny,nz))
pi = 4.*atan(1.)
a = 0.
do i = 1, nx
   a(i,1,1) = sin(2.*pi*float(i-1)/float(nx))
   a(i,2,1) = sin(2.*pi*float(i-1)/float(nx))
!  a(i,1,1) = cmplx(sin(2.*pi*float(i-1)/float(nx)),cos(2.*pi*float(i-1)/float(nx)))
enddo
b = a
!print *,'org: ',a
!$acc data copy(a)
!$acc host_data use_device(a)
call cud2z( nx, ny, nz, a )
!$acc end host_data
!$acc end data
a = a/nx
do i = 1, nx
  print '(A7,1x,i2,2(1x,f12.5))','trans: ',i,a(i,1,1),a(i,2,1)
enddo
!$acc data copy(a)
!$acc host_data use_device(a)
call cuz2d( nx, ny, nz, a )
!$acc end host_data
!$acc end data
do i = 1, nx
  print '(A6,1x,i2,4(1x,f12.5))','back: ',i,b(i,1,1),a(i,1,1),b(i,2,1),a(i,2,1)
enddo
end
