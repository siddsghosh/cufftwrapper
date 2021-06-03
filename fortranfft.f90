module scratch_data
   implicit none
   integer :: init_state = 1, deljx, deliz
   real(8), allocatable :: a2d_c(:,:,:,:)
contains
   subroutine alloc_a2d_c( ny, jxs, jxe, iz1, iz2 )
   integer, intent(in) :: ny, jxs, jxe, iz1, iz2 
   deliz = (iz2 - iz1 + 1)
   deljx = (jxe - jxs + 1)
   if (mod(deljx,2) .eq. 0)then
      deljx = deljx/2
   else
      deljx = deljx/2 + 1
   endif
   deliz = (iz2 - iz1 + 1)
   allocate(a2d_c(2,ny,deljx,deliz))
   a2d_c = 0.
   init_state = 0
!!$acc enter data create(a2d_c)
   end subroutine alloc_a2d_c
end module scratch_data


subroutine fortcud2z( nx, iys, iye, iz1, iz2, fn, ax )
implicit none
integer :: nx, iys, iye, iz1, iz2, iz, iy, ix
real(8) :: ax(nx+2,iys:iye,iz1:iz2), fn
!!$acc data copy(ax), copyin(fn)
!!$acc kernels !collapse(3)
do iz=iz1,iz2
   do iy=iys,iye
      do ix=1,nx
         ax(ix,iy,iz) = ax(ix,iy,iz)*fn
      enddo
      ax(nx+1,iy,iz) = 0.; ax(nx+2,iy,iz) = 0.
   enddo
enddo
!!$acc end kernels
!!$acc end data
!$acc data copy(ax)
!$acc host_data use_device(ax)
call cud2z(nx, (iye-iys+1), (iz2-iz1+1), ax)
!$acc end host_data
!$acc end data
end subroutine fortcud2z


subroutine fortcuz2d( nx, iys, iye, iz1, iz2, ax )
implicit none
integer :: nx, iys, iye, iz1, iz2
real(8) :: ax(nx+2,iys:iye,iz1:iz2)
!$acc data copy(ax)
!$acc host_data use_device(ax)
call cuz2d(nx, (iye-iys+1), (iz2-iz1+1), ax)
!$acc end host_data
!$acc end data
end subroutine fortcuz2d


subroutine fortcuz2zf( ny, jxs, jxe, iz1, iz2, at )
use scratch_data
implicit none
integer :: ny, nx, nz, jxs, jxe, iz1, iz2, iz, ix, iy, ix1
real(8) :: at(ny,jxs:jxe,iz1:iz2)
if (init_state.eq.1) call alloc_a2d_c( ny, jxs, jxe, iz1, iz2 )
!!$acc data copy(at)
!!$acc kernels !collapse(3)
do iz=iz1,iz2
   ix1 = 1
   do ix=jxs,jxe,2
      do iy=1,ny
         a2d_c(1,iy,ix1,iz) = at(iy,ix,iz)
         a2d_c(2,iy,ix1,iz) = at(iy,ix+1,iz)
      enddo
      ix1 = ix1 + 1
   enddo
enddo
!!$acc end kernels
!!$acc end data
!$acc data copy(a2d_c)
!$acc host_data use_device(a2d_c)
call cuz2zf(ny, deljx, deliz, a2d_c)
!$acc end host_data
!$acc end data
!!$acc kernels !collapse(3)
do iz=iz1,iz2
   ix1 = 1
   do ix=jxs,jxe,2
      do iy=1,ny
         at(iy,ix,iz)   = a2d_c(1,iy,ix1,iz)
         at(iy,ix+1,iz) = a2d_c(2,iy,ix1,iz)
      enddo
      ix1 = ix1 + 1
   enddo
enddo
!!$acc end kernels
!!$acc end data
end subroutine fortcuz2zf


subroutine fortcuz2zb( ny, jxs, jxe, iz1, iz2, at )
use scratch_data
implicit none
integer :: ny, jxs, jxe, ix1, iz1, iz2, iz, ix, iy
real(8) :: at(ny,jxs:jxe,iz1:iz2)
if (init_state.eq.1) call alloc_a2d_c( ny, jxs, jxe, iz1, iz2 )
!!$acc data copy(at)
!!$acc kernels !collapse(3)
do iz=iz1,iz2
   ix1 = 1
   do ix=jxs,jxe,2
      do iy=1,ny
         a2d_c(1,iy,ix1,iz) = at(iy,ix,iz)
         a2d_c(2,iy,ix1,iz) = at(iy,ix+1,iz)
      enddo
      ix1 = ix1 + 1
   enddo
enddo
!!$acc end kernels
!!$acc end data
!$acc data copy(a2d_c)
!$acc host_data use_device(a2d_c)
call cuz2zb(ny, deljx, deliz, a2d_c)
!$acc end host_data
!$acc end data
!!$acc kernels !collapse(3)
do iz=iz1,iz2
   ix1 = 1
   do ix=jxs,jxe,2
      do iy=1,ny
         at(iy,ix,iz)   = a2d_c(1,iy,ix1,iz)
         at(iy,ix+1,iz) = a2d_c(2,iy,ix1,iz)
      enddo
      ix1 = ix1 + 1
   enddo
enddo
!!$acc end kernels
!!$acc end data
end subroutine fortcuz2zb
