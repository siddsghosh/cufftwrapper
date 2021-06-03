module cufft2d_mpi_mod
   private
   public :: cufft2d_mpi

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
   integer :: nx,ny,jxs,jxe,iys,iye,iz1,iz2,myid,ncpu,isgn, nxp2
   real :: ax(nx+2,iys:iye,iz1:iz2), at(ny,jxs:jxe,iz1:iz2), fn
   integer :: jx_s(0:ncpu-1), jx_e(0:ncpu-1),iy_s(0:ncpu-1), iy_e(0:ncpu-1)

   nxp2 = nx + 2
   if(isgn .lt. 0) then
      fn = 1.0/(float(nx)*float(ny))

! ------ 1d fft in x over [iys,iye] for all z

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

   else

! ---- decide whether to first transpose or leave as is

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
   endif
   return
   end

   subroutine xtoy_trans(f,g,nx,ny,ixs,ixe,ix_s,ix_e,iys,iye, &
       iy_s,iy_e,iz1,iz2,myid,ncpu_s,np)
! 
! ------- transpose array  f(nx,iys:iye,iz1:iz2) 
!                     ---> g(ny,ixs:ixe,iz1:iz2)
!         load balanced version, march 2017
!
!     use pars, only : newcomm
   use mpi
   implicit none
   integer :: istatus(mpi_status_size),newcomm,nx,ny,ixs,ixe,iys,iye, &
              iz1,iz2,myid,ncpu_s,np,jk,ik,iss,is,ir,i,nsend,nrecv, &
              ireqs, ireqr, ierr

   real :: f(nx,iys:iye,iz1:iz2),g(ny,ixs:ixe,iz1:iz2)
   real :: ft(nx*(iye+1-iys)*(iz2 - iz1 + 1)), &
       gt(ny*(ixe+1-ixs)*(iz2 - iz1 + 1))
   integer :: ix_s(0:np-1),ix_e(0:np-1),iy_s(0:np-1),iy_e(0:np-1)
 
   newcomm = mpi_comm_world
   jk = (iye - iys + 1)*(iz2 - iz1 + 1)
   ik = (ixe - ixs + 1)*(iz2 - iz1 + 1)

! ----------- loop over cpus on a slab for given myid

   iss = (myid/ncpu_s)*ncpu_s

   do i=1,ncpu_s
      is    = mod(myid - iss + i,ncpu_s) + iss
      ir    = mod(myid - iss + (ncpu_s - i),ncpu_s) + iss
      nsend = (ix_e(is) - ix_s(is) + 1)*jk
      nrecv = (iy_e(ir) - iy_s(ir) + 1)*ik
      if(is == myid) then
         call send_xtoy(f,gt(1),nx,ix_s(is),ix_e(is), &
                          iy_s(myid),iy_e(myid),iz1,iz2)
      else
         call send_xtoy(f,ft(1),nx,ix_s(is),ix_e(is), &
                           iy_s(myid),iy_e(myid),iz1,iz2)
         call mpi_isend(ft(1),nsend,mpi_real8,is,1,newcomm,ireqs,ierr)
         call mpi_irecv(gt(1),nrecv,mpi_real8,ir,1,newcomm,ireqr,ierr)
         call mpi_wait(ireqs,istatus,ierr)
         call mpi_wait(ireqr,istatus,ierr)
      endif
      call recv_xtoy(g,gt(1),ny,ix_s(myid),ix_e(myid),iy_s(ir),iy_e(ir),iz1,iz2)
   enddo
   return
   end


   subroutine ytox_trans(g,f,nx,ny,ixs,ixe,ix_s,ix_e,iys,iye,iy_s,iy_e,iz1,iz2, &
        myid,ncpu_s,np)

! ------- transpose array g(ny,ixs:ixe,iz1:iz2) 
!                    ---> f(nx,iys:iye,iz1:iz2)
!         load balanced version, march 2017

   use mpi
   implicit none
   integer :: istatus(mpi_status_size),nx,ny,ixs,ixe,iys,iye,iz1,iz2,myid,ncpu_s,np, &
              jk,ik,i,is,ir,nsend,nrecv,newcomm,iss,ireqs,ireqr,ierr
   real :: f(nx,iys:iye,iz1:iz2),g(ny,ixs:ixe,iz1:iz2)
   real :: ft(nx*(iye+1-iys)*(iz2 - iz1 + 1)),gt(ny*(ixe+1-ixs)*(iz2 - iz1 + 1))
   integer :: ix_s(0:np-1), ix_e(0:np-1), iy_s(0:np-1), iy_e(0:np-1)

   newcomm = mpi_comm_world

   jk = (iye - iys + 1)*(iz2 - iz1 + 1)
   ik = (ixe - ixs + 1)*(iz2 - iz1 + 1)

! ----------- loop thru cpus on a slab

   iss   = (myid/ncpu_s)*ncpu_s
   do i=1,ncpu_s
      is    = mod(myid - iss + i,ncpu_s) + iss
      ir    = mod(myid - iss + (ncpu_s - i),ncpu_s) + iss
      nsend = (iy_e(is) - iy_s(is) + 1)*ik
      nrecv = (ix_e(ir) - ix_s(ir) + 1)*jk
      if(is == myid) then
          call send_ytox(g,ft(1),ny,ix_s(myid),ix_e(myid),iy_s(is),iy_e(is),iz1,iz2)
      else
          call send_ytox(g,gt(1),ny,ix_s(myid),ix_e(myid),iy_s(is),iy_e(is),iz1,iz2)
          call mpi_isend(gt(1),nsend,mpi_real8,is,1,newcomm,ireqs,ierr)
          call mpi_irecv(ft(1),nrecv,mpi_real8,ir,1,newcomm,ireqr,ierr)
          call mpi_wait(ireqs,istatus,ierr)
          call mpi_wait(ireqr,istatus,ierr)
      endif
      call recv_ytox(f,ft(1),nx,ix_s(ir),ix_e(ir),iy_s(myid),iy_e(myid),iz1,iz2)
   enddo
   return
   end

   subroutine xtoz_trans(f,g,nx,nz,ixs,ixe,ix_s,ix_e,iys,iye,izs,ize,iz_s,iz_e, &
           myid,ncpu_s,numprocs)

! ------- transpose array  f(nx,iys:iye,izs-1:ize+1) 
!                     ---> g(0:nz+1,iys:iye,ixs:ixe)
!         load balanced version, march 2017

   use mpi
   implicit none
   integer :: istatus(mpi_status_size),nx,nz,ixs,ixe,iys,iye,izs,ize, &
              myid,ncpu_s,numprocs,newcomm,i,is,ir,nsend,nrecv,jk,ij,nz_c, &
              ireqs, ireqr, ierr
   real :: f(nx,iys:iye,izs-1:ize+1), g(0:nz+1,iys:iye,ixs:ixe)
   real :: ft(nx*(iye+1-iys)*(ize-izs+1)), gt(nz*(ixe+1-ixs)*(iye-iys+1))
   integer :: ix_s(0:numprocs-1), ix_e(0:numprocs-1), iz_s(0:numprocs-1), &
             iz_e(0:numprocs-1)

   newcomm = mpi_comm_world
   jk = (ize - izs + 1)*(iye - iys + 1)
   ij = (ixe - ixs + 1)*(iye - iys + 1)

! ----------- find number of vertical slabs
!             build send and recv ids

   nz_c = numprocs/ncpu_s

   do i=1,nz_c
      is    = mod(myid + i*ncpu_s,numprocs)
      ir    = mod(myid + (nz_c - i)*ncpu_s,numprocs)
      nsend = (ix_e(is) - ix_s(is) + 1)*jk
      nrecv = (iz_e(ir) - iz_s(ir) + 1)*ij
      if(is == myid) then
         call send_xtoz(f,gt(1),nx,ix_s(is),ix_e(is),iys,iye,iz_s(myid),iz_e(myid))
      else
         call send_xtoz(f,ft(1),nx,ix_s(is),ix_e(is), iys,iye,iz_s(myid),iz_e(myid))
         call mpi_isend(ft(1),nsend,mpi_real8,is,1,newcomm,ireqs,ierr)
         call mpi_irecv(gt(1),nrecv,mpi_real8,ir,1,newcomm,ireqr,ierr)
         call mpi_wait(ireqs,istatus,ierr)
         call mpi_wait(ireqr,istatus,ierr)
      endif
      call recv_xtoz(g,gt(1),nz,ix_s(myid),ix_e(myid),iys,iye,iz_s(ir),iz_e(ir))
   enddo
   return
   end


   subroutine ztox_trans(g,f,nx,nz,ixs,ixe,ix_s,ix_e,iys,iye,izs,ize,iz_s,iz_e, &
        myid,ncpu_s,numprocs)

! ------- transpose array g(0:nz+1,iys:iye,ixs:ixe) 
!                    ---> f(nx,iys:iye,izs-1:ize+1)
   use mpi
   implicit none
   integer :: istatus(mpi_status_size),newcomm,jk,ij,i,is,ir,nsend,myid, &
              nrecv,nx,nz,ixs,ixe,iys,iye,izs,ize,ncpu_s,numprocs, nz_c, &
              ireqs, ireqr, ierr
   real :: f(nx,iys:iye,izs-1:ize+1), g(0:nz+1,iys:iye,ixs:ixe)
   real :: ft(nx*(iye+1-iys)*(ize-izs+3)), gt((nz+3)*(iye+1-iys)*(ixe-ixs+1))
   integer :: ix_s(0:numprocs-1), ix_e(0:numprocs-1), &
              iz_s(0:numprocs-1), iz_e(0:numprocs-1)

   newcomm = mpi_comm_world

   jk = (ize - izs + 3)*(iye - iys + 1)
   ij = (ixe - ixs + 1)*(iye - iys + 1)

! ------------- get number of vertical slabs
!               and build send and recv ids

   nz_c = numprocs/ncpu_s

   do i=1,nz_c
      is    = mod(myid + i*ncpu_s,numprocs)
      ir    = mod(myid + (nz_c - i)*ncpu_s,numprocs)
      nsend = (iz_e(is) - iz_s(is) + 3)*ij
      nrecv = (ix_e(ir) - ix_s(ir) + 1)*jk
      if(is == myid) then
         call send_ztox(g,ft(1),nz,ix_s(myid),ix_e(myid),iys,iye,iz_s(is),iz_e(is))
      else
         call send_ztox(g,gt(1),nz,ix_s(myid),ix_e(myid),iys,iye,iz_s(is),iz_e(is))
         call mpi_isend(gt(1),nsend,mpi_real8,is,1,newcomm,ireqs,ierr)
         call mpi_irecv(ft(1),nrecv,mpi_real8,ir,1,newcomm,ireqr,ierr)
         call mpi_wait(ireqs,istatus,ierr)
         call mpi_wait(ireqr,istatus,ierr)
      endif
      call recv_ztox(f,ft(1),nx,ix_s(ir),ix_e(ir),iys,iye,iz_s(myid),iz_e(myid))
   enddo
   return
   end


   subroutine send_xtoy(f,ft,nx,ixs,ixe,iys,iye,izs,ize)

! ------------- grab correct chunk of array to be sent

   implicit none
   integer :: i, j, k, nx,ixs,ixe,iys,iye,izs,ize
   real :: f(nx,iys:iye,izs:ize), ft(ixs:ixe,iys:iye,izs:ize)

   do k=izs,ize
      do j=iys,iye
         do i=ixs,ixe
            ft(i,j,k) = f(i,j,k)
         enddo
      enddo
   enddo
   return
   end

   subroutine recv_xtoy(g,gt,ny,ixs,ixe,iys,iye,izs,ize)
   implicit none
   integer :: ny,ixs,ixe,iys,iye,izs,ize,i,j,k
   real :: g(ny,ixs:ixe,izs:ize), gt(ixs:ixe,iys:iye,izs:ize)

   do k=izs,ize
      do j=iys,iye
         do i=ixs,ixe
            g(j,i,k) = gt(i,j,k)
         enddo
      enddo
   enddo
   return
   end

   subroutine send_ytox(g,gt,ny,ixs,ixe,iys,iye,izs,ize)

! ------------- grab correct chunk of array to be sent

   implicit none
   integer :: ny,ixs,ixe,iys,iye,izs,ize,i,j,k
   real :: g(ny,ixs:ixe,izs:ize), gt(iys:iye,ixs:ixe,izs:ize)

   do k=izs,ize
      do i=ixs,ixe
         do j=iys,iye
            gt(j,i,k) = g(j,i,k)
         enddo
      enddo
   enddo
   return
   end

   subroutine recv_ytox(f,ft,nx,ixs,ixe,iys,iye,izs,ize)
   implicit none
   integer :: nx,ixs,ixe,iys,iye,izs,ize,i,j,k
   real :: f(nx,iys:iye,izs:ize), ft(iys:iye,ixs:ixe,izs:ize)

   do k=izs,ize
      do i=ixs,ixe
         do j=iys,iye
            f(i,j,k) = ft(j,i,k)
         enddo
      enddo
   enddo
   return
   end


   subroutine send_xtoz(f,ft,nx,ixs,ixe,iys,iye,izs,ize)

! ------- grab correct chunk of array to be sent and skip ghost points

   implicit none
   integer :: nx,ixs,ixe,iys,iye,izs,ize,i,j,k
   real :: f(nx,iys:iye,izs-1:ize+1), ft(ixs:ixe,iys:iye,izs:ize)

   do k=izs,ize
      do j=iys,iye
         do i=ixs,ixe
            ft(i,j,k) = f(i,j,k)
         enddo
      enddo
   enddo
   return
   end

   subroutine recv_xtoz(g,gt,nz,ixs,ixe,iys,iye,izs,ize)
   implicit none
   integer :: i,k,j,nz,ixs,ixe,iys,iye,izs,ize
   real :: g(0:nz+1,iys:iye,ixs:ixe), gt(ixs:ixe,iys:iye,izs:ize)

   do k=izs,ize
      do j=iys,iye
         do i=ixs,ixe
            g(k,j,i) = gt(i,j,k)
         enddo
      enddo
   enddo
   return
   end


   subroutine send_ztox(g,gt,nz,ixs,ixe,iys,iye,izs,ize)

! ------------- grab correct chunk of array to be sent,
!               account for ghost points

   implicit none
   integer :: nz,ixs,ixe,iys,iye,izs,ize, i, j, k
   real :: g(0:nz+1,iys:iye,ixs:ixe), gt(izs-1:ize+1,iys:iye,ixs:ixe)

   do j=iys,iye
      do i=ixs,ixe
         do k=izs-1,ize+1
            gt(k,j,i) = g(k,j,i)
         enddo
      enddo
   enddo
   return
   end

   subroutine recv_ztox(f,ft,nx,ixs,ixe,iys,iye,izs,ize)
   implicit none
   integer :: nx,ixs,ixe,iys,iye,izs,ize,i,j,k
   real :: f(nx,iys:iye,izs-1:ize+1), ft(izs-1:ize+1,iys:iye,ixs:ixe)

   do i=ixs,ixe
      do j=iys,iye
         do k=izs-1,ize+1
            f(i,j,k) = ft(k,j,i)
         enddo
      enddo
   enddo
   return
   end

end module cufft2d_mpi_mod
