module cufft2d_mpi_mod
   use mpi
   implicit none
   private
   public :: cufft2d_mpi
   integer :: init_a2d_c = 1, init_ftgt = 1, newcomm
   integer, allocatable :: ireqs(:), ireqr(:), status(:,:)
   real(8), allocatable :: a2d_c(:,:,:,:), ft(:,:), gt(:,:)
   integer :: deljx, deliz

contains
   subroutine cufft2d_mpi(ax,at,nx,ny,jxs,jxe,jx_s,jx_e,iys,iye, &
       iy_s,iy_e,iz1,iz2,myid,ncpu,isgn)
!
! -------- get 2d fft using fftpack routines and parallel mpi
!          use fftpack storage a0, (a1,b1), (a2,b2),...,
!
!         isgn = -1 do forward transform, get coefficients
!                   incoming array is ax(nx+2,iys:iye,iz1:iz2)
!                   outgoing array is ax(nx+2,iys:iye,iz1:iz2)
!
!         isgn = -2 do forward transform, get coefficients
!                   incoming array is ax(nx+2,iys:iye,iz1:iz2)
!                   outgoing array is at(ny,jxs:jxe,iz1:iz2)
!
!         isgn =  1 do inverse transform, move to physical space
!                   incoming array is ax(nx+2,iys:iye,iz1:iz2)
!                   outgoing array is ax(nx+2,iys:iye,iz1:iz2)
!
!         isgn =  2 do inverse transform, move to physical space
!                   incoming array is at(ny,jxs:jxe,iz1:iz2)
!                   outgoing array is ax(nx+2,iys:iye,iz1:iz2)
!
   implicit none
   integer :: nx,ny,jxs,jxe,iys,iye,iz1,iz2,myid,ncpu,isgn,nxp2,ix,iy,iz
   real :: ax(nx+2,iys:iye,iz1:iz2), at(ny,jxs:jxe,iz1:iz2), fn
   integer :: jx_s(0:ncpu-1), jx_e(0:ncpu-1),iy_s(0:ncpu-1), iy_e(0:ncpu-1)

   nxp2 = nx + 2
   if(isgn .lt. 0) then
      fn = 1.0/(float(nx)*float(ny))

! ------ 1d fft in x over [iys,iye] for all z

      !$acc data copyin(ax) copyout(at)
      call fortcud2z( nx, iys, iye, iz1, iz2, fn, ax )

      call xtoy_trans(ax,at,nxp2,ny,jxs,jxe,jx_s,jx_e, &
             iys,iye,iy_s,iy_e,iz1,iz2,myid,ncpu,ncpu)

! ------ 1d fft in y over [jxs,jxe] for all z

      call fortcuz2zf( ny, jxs, jxe, iz1, iz2, at )
 
! ---- decide whether to transpose back or leave as is
 
      if(isgn .eq. -1) then
         call ytox_trans(at,ax,nxp2,ny,jxs,jxe,jx_s,jx_e, &
                    iys,iye,iy_s,iy_e,iz1,iz2,myid,ncpu,ncpu)
      endif
      !$acc end data

   else

! ---- decide whether to first transpose or leave as is

      !$acc data copyin(at) copyout(ax)
      if(isgn .eq. 1) then
         call xtoy_trans(ax,at,nxp2,ny,jxs,jxe,jx_s,jx_e, &
                   iys,iye,iy_s,iy_e,iz1,iz2,myid,ncpu,ncpu)
      endif
 
! ------ 1d fft in y over [jxs,jxe] for all z
 
      call fortcuz2zb( ny, jxs, jxe, iz1, iz2, at )

      call ytox_trans(at,ax,nxp2,ny,jxs,jxe,jx_s,jx_e,  &
                 iys,iye,iy_s,iy_e,iz1,iz2,myid,ncpu,ncpu)

! ------  1d fft in x over [iys,iye] for all z

      call fortcuz2d( nx, iys, iye, iz1, iz2, ax )
      !$acc end data
   endif
   return
   end


   subroutine alloc_ftgt(nx,ny,ixs,ixe,iys,iye,iz2,iz1,ncpu_s)
   implicit none
   integer :: nx,ny,ixs,ixe,iys,iye,iz2,iz1,ncpu_s
   allocate(ft(nx*(iye+1-iys)*(iz2-iz1+1),ncpu_s))
   allocate(gt(ny*(ixe+1-ixs)*(iz2-iz1+1),ncpu_s))
   allocate(ireqs(ncpu_s),ireqr(ncpu_s),status(MPI_STATUS_SIZE,ncpu_s))
   !$acc enter data create(ft,gt)
   newcomm = mpi_comm_world
   init_ftgt = -1
   end
   

   subroutine xtoy_trans(f,g,nx,ny,ixs,ixe,ix_s,ix_e,iys,iye, &
       iy_s,iy_e,iz1,iz2,myid,ncpu_s,np)
! 
! ------- transpose array  f(nx,iys:iye,iz1:iz2) 
!                     ---> g(ny,ixs:ixe,iz1:iz2)
!         load balanced version, march 2017
!
   implicit none
   integer :: nx,ny,ixs,ixe,iys,iye,iz1,iz2,myid,ncpu_s,np,&
              jk,ik,iss,is,ir,i,nsend,nrecv,ierr,j

   real :: f(nx,iys:iye,iz1:iz2),g(ny,ixs:ixe,iz1:iz2)
   integer :: ix_s(0:np-1),ix_e(0:np-1),iy_s(0:np-1),iy_e(0:np-1)
 
   jk = (iye - iys + 1)*(iz2 - iz1 + 1)
   ik = (ixe - ixs + 1)*(iz2 - iz1 + 1)

   if(init_ftgt==1) call alloc_ftgt(nx,ny,ixs,ixe,iys,iye,iz2,iz1,ncpu_s)

! ----------- loop over cpus on a slab for given myid

   iss = (myid/ncpu_s)*ncpu_s
   j = 0
   do i=1,ncpu_s
      is    = mod(myid - iss + i,ncpu_s) + iss
      ir    = mod(myid - iss + (ncpu_s - i),ncpu_s) + iss
      nsend = (ix_e(is) - ix_s(is) + 1)*jk
      nrecv = (iy_e(ir) - iy_s(ir) + 1)*ik
      if(is == myid) then
         call send_xtoy(f,gt(1,i),nx,ix_s(is),ix_e(is),iy_s(myid),iy_e(myid),iz1,iz2)
      else
         call send_xtoy(f,ft(1,i),nx,ix_s(is),ix_e(is),iy_s(myid),iy_e(myid),iz1,iz2)
         j = j + 1
         !$acc host_data use_device(ft,gt)
         call mpi_isend(ft(1,i),nsend,mpi_real8,is,1,newcomm,ireqs(j),ierr)
         call mpi_irecv(gt(1,i),nrecv,mpi_real8,ir,1,newcomm,ireqr(j),ierr)
         !$acc end host_data
      endif
   enddo
   call MPI_Waitall((ncpu_s-1),ireqr,status,ierr)
   do i=1,ncpu_s
      ir    = mod(myid - iss + (ncpu_s - i),ncpu_s) + iss
      call recv_xtoy(g,gt(1,i),ny,ix_s(myid),ix_e(myid),iy_s(ir),iy_e(ir),iz1,iz2)
   enddo
   call MPI_Waitall((ncpu_s-1),ireqs,status,ierr)
   return
   end


   subroutine ytox_trans(g,f,nx,ny,ixs,ixe,ix_s,ix_e,iys,iye,iy_s,iy_e,iz1,iz2, &
        myid,ncpu_s,np)

! ------- transpose array g(ny,ixs:ixe,iz1:iz2) 
!                    ---> f(nx,iys:iye,iz1:iz2)
!         load balanced version, march 2017

   implicit none
   integer :: nx,ny,ixs,ixe,iys,iye,iz1,iz2,myid,ncpu_s,np, &
              jk,ik,i,is,ir,nsend,nrecv,iss,ierr,j
   real :: f(nx,iys:iye,iz1:iz2),g(ny,ixs:ixe,iz1:iz2)
   integer :: ix_s(0:np-1), ix_e(0:np-1), iy_s(0:np-1), iy_e(0:np-1)


   jk = (iye - iys + 1)*(iz2 - iz1 + 1)
   ik = (ixe - ixs + 1)*(iz2 - iz1 + 1)

! ----------- loop thru cpus on a slab
   if(init_ftgt==1) call alloc_ftgt(nx,ny,ixs,ixe,iys,iye,iz2,iz1,ncpu_s)

   iss   = (myid/ncpu_s)*ncpu_s
   j = 0
   do i=1,ncpu_s
      is    = mod(myid - iss + i,ncpu_s) + iss
      ir    = mod(myid - iss + (ncpu_s - i),ncpu_s) + iss
      nsend = (iy_e(is) - iy_s(is) + 1)*ik
      nrecv = (ix_e(ir) - ix_s(ir) + 1)*jk
      if(is == myid) then
          call send_ytox(g,ft(1,i),ny,ix_s(myid),ix_e(myid),iy_s(is),iy_e(is),iz1,iz2)
      else
          call send_ytox(g,gt(1,i),ny,ix_s(myid),ix_e(myid),iy_s(is),iy_e(is),iz1,iz2)
          j = j + 1
          !$acc host_data use_device(ft,gt)
          call mpi_isend(gt(1,i),nsend,mpi_real8,is,1,newcomm,ireqs(j),ierr)
          call mpi_irecv(ft(1,i),nrecv,mpi_real8,ir,1,newcomm,ireqr(j),ierr)
          !$acc end host_data
      endif
   enddo
   call MPI_Waitall((ncpu_s-1),ireqr,status,ierr)
   do i=1,ncpu_s
      ir    = mod(myid - iss + (ncpu_s - i),ncpu_s) + iss
      call recv_ytox(f,ft(1,i),nx,ix_s(ir),ix_e(ir),iy_s(myid),iy_e(myid),iz1,iz2)
   enddo
   call MPI_Waitall((ncpu_s-1),ireqs,status,ierr)
   return
   end

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
   init_a2d_c = 0
   !$acc enter data create(a2d_c)
   end subroutine alloc_a2d_c


   subroutine fortcud2z( nx, iys, iye, iz1, iz2, fn, ax )
   implicit none
   integer :: nx, iys, iye, iz1, iz2, iz, iy, ix
   real(8) :: ax(nx+2,iys:iye,iz1:iz2), fn
   !$acc data present(ax)
   !$acc parallel loop collapse(3)
   do iz=iz1,iz2
      do iy=iys,iye
         do ix=1,nx
            ax(ix,iy,iz) = ax(ix,iy,iz)*fn
         enddo
      enddo
   enddo
   !$acc end parallel
   !$acc host_data use_device(ax)
   call cud2z(nx, (iye-iys+1), (iz2-iz1+1), ax)
   !$acc end host_data
   !$acc end data
   end subroutine fortcud2z


   subroutine fortcuz2d( nx, iys, iye, iz1, iz2, ax )
   implicit none
   integer :: nx, iys, iye, iz1, iz2
   real(8) :: ax(nx+2,iys:iye,iz1:iz2)
   !$acc host_data use_device(ax)
   call cuz2d(nx, (iye-iys+1), (iz2-iz1+1), ax)
   !$acc end host_data
   end subroutine fortcuz2d


   subroutine fortcuz2zf( ny, jxs, jxe, iz1, iz2, at )
   implicit none
   integer :: ny, nx, nz, jxs, jxe, iz1, iz2, iz, ix, iy, ix1
   real(8) :: at(ny,jxs:jxe,iz1:iz2)
   if (init_a2d_c.eq.1) call alloc_a2d_c( ny, jxs, jxe, iz1, iz2 )
   !$acc data present(a2d_c,at)
   !$acc parallel loop collapse(3)
   do iz=iz1,iz2
      do ix=jxs,jxe,2
         do iy=1,ny
            ix1 = (ix-jxs)/2 + 1
            a2d_c(1,iy,ix1,iz) = at(iy,ix,iz)
            a2d_c(2,iy,ix1,iz) = at(iy,ix+1,iz)
         enddo
      enddo
   enddo
   !$acc end parallel
   !$acc host_data use_device(a2d_c)
   call cuz2zf(ny, deljx, deliz, a2d_c)
   !$acc end host_data
   !$acc parallel loop collapse(3)
   do iz=iz1,iz2
      do ix=jxs,jxe,2
         do iy=1,ny
            ix1 = (ix-jxs)/2 + 1
            at(iy,ix,iz)   = a2d_c(1,iy,ix1,iz)
            at(iy,ix+1,iz) = a2d_c(2,iy,ix1,iz)
         enddo
      enddo
   enddo
   !$acc end parallel
   !$acc end data
   end subroutine fortcuz2zf


   subroutine fortcuz2zb( ny, jxs, jxe, iz1, iz2, at )
   implicit none
   integer :: ny, jxs, jxe, ix1, iz1, iz2, iz, ix, iy
   real(8) :: at(ny,jxs:jxe,iz1:iz2)
   if (init_a2d_c.eq.1) call alloc_a2d_c( ny, jxs, jxe, iz1, iz2 )
   !$acc data present(a2d_c,at)
   !$acc parallel loop collapse(3)
   do iz=iz1,iz2
      do ix=jxs,jxe,2
         do iy=1,ny
            ix1 = (ix-jxs)/2 + 1
            a2d_c(1,iy,ix1,iz) = at(iy,ix,iz)
            a2d_c(2,iy,ix1,iz) = at(iy,ix+1,iz)
         enddo
      enddo
   enddo
   !$acc end parallel
   !$acc host_data use_device(a2d_c)
   call cuz2zb(ny, deljx, deliz, a2d_c)
   !$acc end host_data
   !$acc parallel loop collapse(3)
   do iz=iz1,iz2
      do ix=jxs,jxe,2
         do iy=1,ny
            ix1 = (ix-jxs)/2 + 1
            at(iy,ix,iz)   = a2d_c(1,iy,ix1,iz)
            at(iy,ix+1,iz) = a2d_c(2,iy,ix1,iz)
         enddo
      enddo
   enddo
   !$acc end parallel
   !$acc end data
   end subroutine fortcuz2zb
end module cufft2d_mpi_mod


subroutine send_xtoy(f,ft,nx,ixs,ixe,iys,iye,izs,ize)

! ------------- grab correct chunk of array to be sent

implicit none
integer :: i, j, k, nx,ixs,ixe,iys,iye,izs,ize
real :: f(nx,iys:iye,izs:ize), ft(ixs:ixe,iys:iye,izs:ize)
!$acc data present(ft,f)
!$acc parallel loop collapse(3)
do k=izs,ize
   do j=iys,iye
      do i=ixs,ixe
         ft(i,j,k) = f(i,j,k)
      enddo
   enddo
enddo
!$acc end parallel
!$acc end data
return
end

subroutine recv_xtoy(g,gt,ny,ixs,ixe,iys,iye,izs,ize)
implicit none
integer :: ny,ixs,ixe,iys,iye,izs,ize,i,j,k
real :: g(ny,ixs:ixe,izs:ize), gt(ixs:ixe,iys:iye,izs:ize)
!$acc data present(g,gt)
!$acc parallel loop collapse(3)
do k=izs,ize
   do j=iys,iye
      do i=ixs,ixe
         g(j,i,k) = gt(i,j,k)
      enddo
   enddo
enddo
!$acc end parallel
!$acc end data
return
end

subroutine send_ytox(g,gt,ny,ixs,ixe,iys,iye,izs,ize)

! ------------- grab correct chunk of array to be sent

implicit none
integer :: ny,ixs,ixe,iys,iye,izs,ize,i,j,k
real :: g(ny,ixs:ixe,izs:ize), gt(iys:iye,ixs:ixe,izs:ize)
!$acc data present(gt,g)
!$acc parallel loop collapse(3)
do k=izs,ize
   do i=ixs,ixe
      do j=iys,iye
         gt(j,i,k) = g(j,i,k)
      enddo
   enddo
enddo
!$acc end parallel
!$acc end data
return
end

subroutine recv_ytox(f,ft,nx,ixs,ixe,iys,iye,izs,ize)
implicit none
integer :: nx,ixs,ixe,iys,iye,izs,ize,i,j,k
real :: f(nx,iys:iye,izs:ize), ft(iys:iye,ixs:ixe,izs:ize)
!$acc data present(ft,f)
!$acc parallel loop collapse(3)
do k=izs,ize
   do i=ixs,ixe
      do j=iys,iye
         f(i,j,k) = ft(j,i,k)
      enddo
   enddo
enddo
!$acc end parallel
!$acc end data
return
end
