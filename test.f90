program testr2c
implicit none
integer :: nx, ny, nz, i
real, dimension(:), allocatable :: a, b
real :: pi
nx = 16
ny = 1
nz = 1
allocate(a(nx+2),b(nx+2))
pi = 4.*atan(1.)
do i = 1, nx
   a(i) = sin(2.*pi*float(i-1)/float(nx))
enddo
a(nx+1) = 0.;a(nx+2) = 0.
b = a
print *,'org: ',a
!$acc data copy(a)
!$acc host_data use_device(a)
call cud2z( nx, ny, nz, a )
!$acc end host_data
!$acc end data
a = a/nx
do i = 1, nx
  print *,'trans: ',i,a(i)
enddo
!$acc data copy(a)
!$acc host_data use_device(a)
call cuz2d( nx, ny, nz, a )
!$acc end host_data
!$acc end data
do i = 1, nx
  print *,'back: ',i,b(i),a(i)
enddo
end
