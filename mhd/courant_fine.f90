subroutine courant_fine(ilevel)
  use amr_commons
  use hydro_commons
  use hydro_parameters,only:cr_vmax, cr_c_fraction, iCRu
  use poisson_commons
  use mpi_mod
#if USE_TURB==1
  use turb_commons
#endif
  implicit none
#ifdef CRFLX
  integer,parameter::ncomm=5
  real(kind=8)::ecrs_loc,ecrs_all
#else
  integer,parameter::ncomm=4
#endif 
#ifndef WITHOUTMPI
  integer::info
  real(kind=8),dimension(ncomm)::comm_buffin,comm_buffout
#endif
  integer::ilevel
  !----------------------------------------------------------------------
  ! Using the Courant-Friedrich-Levy stability condition,               !
  ! this routine computes the maximum allowed time-step.                !
  !----------------------------------------------------------------------
  integer::i,ivar,idim,ind,ncache,igrid,iskip
  integer::nleaf,ngrid,nx_loc
  integer,dimension(1:nvector),save::ind_grid,ind_cell,ind_leaf

  real(dp)::dt_lev,dx,vol,scale
  real(kind=8)::mass_loc,ekin_loc,eint_loc,emag_loc,dt_loc
  real(kind=8)::mass_all,ekin_all,eint_all,emag_all,dt_all
  real(dp),dimension(1:nvector,1:nvar+3+ncrvars),save::uu
  real(dp),dimension(1:nvector,1:ndim),save::gg
#ifdef CRFLX
  real(dp)::dt_cr, cr_vgas_max_all
#endif


  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  mass_all=0.0d0; mass_loc=0.0d0
  ekin_all=0.0d0; ekin_loc=0.0d0
  emag_all=0.0d0; emag_loc=0.0d0
  eint_all=0.0d0; eint_loc=0.0d0
#ifdef CRFLX
  ecrs_all=0.0d0; ecrs_loc=0.0d0
  cr_vgas_max = 0.0d0
#endif
  dt_all=dtnew(ilevel); dt_loc=dt_all

  ! Mesh spacing at that level
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx=0.5D0**ilevel*scale
  vol=dx**ndim

  if (ischeme .eq. 1) then
     CALL velocity_fine(ilevel)
     do ivar=1,nvar+3
        call make_virtual_fine_dp(uold(1,ivar),ilevel)
     end do
     if(simple_boundary)call make_boundary_hydro(ilevel)
  endif

  ! Loop over active grids by vector sweeps
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do

     ! Loop over cells
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=ind_grid(i)+iskip
        end do

        ! Gather leaf cells
        nleaf=0
        do i=1,ngrid
           if(son(ind_cell(i))==0)then
              nleaf=nleaf+1
              ind_leaf(nleaf)=ind_cell(i)
           end if
        end do

        ! Gather hydro variables
        do ivar=1,nvar+3+ncrvars
           do i=1,nleaf
              uu(i,ivar)=uold(ind_leaf(i),ivar)
           end do
        end do

        ! Gather gravitational acceleration
        gg=0.0d0
        if(poisson)then
           do idim=1,ndim
              do i=1,nleaf
                 gg(i,idim)=f(ind_leaf(i),idim)
              end do
           end do
        end if


#if USE_TURB==1
        if (turb .AND. turb_type/=3) then
           do idim=1,ndim
              do i=1,nleaf
                 gg(i,idim)=gg(i,idim)+fturb(ind_leaf(i),idim)
              end do
           end do
        end if
#endif

        ! Compute total mass
        do i=1,nleaf
           mass_loc=mass_loc+uu(i,1)*vol
        end do

        ! Compute total energy
        do i=1,nleaf
           ekin_loc=ekin_loc+uu(i,5)*vol
        end do

        ! Compute total magnetic energy
        do ivar=1,3
           do i=1,nleaf
              emag_loc=emag_loc+0.125d0*(uu(i,5+ivar)+uu(i,nvar+ivar))**2*vol
           end do
        end do

        ! Compute total internal energy
        do i=1,nleaf
           eint_loc=eint_loc+uu(i,5)*vol
        end do
        do ivar=1,3
           do i=1,nleaf
              eint_loc=eint_loc-0.5d0*uu(i,1+ivar)**2/uu(i,1)*vol &
                   & -0.125d0*(uu(i,5+ivar)+uu(i,nvar+ivar))**2*vol
           end do
        end do
#if NENER>0
        do ivar=1,nener
           do i=1,nleaf
              eint_loc=eint_loc-uu(i,8+ivar)*vol
           end do
        end do
#endif

        ! Compute cosmic rays energy
#ifdef CRFLX
        do ivar=1,ncr
           do i=1,nleaf
              ecrs_loc=ecrs_loc+uold(ind_leaf(i),icru+(ndim+1)*(ivar-1)) * vol
           end do
        end do
#endif
        ! Compute CFL time-step
        if(nleaf>0)then
           call cmpdt(uu,gg,dx,dt_lev,nleaf)
           dt_loc=min(dt_loc,dt_lev)
        end if

     end do
     ! End loop over cells

  end do
  ! End loop over grids

  ! Compute global quantities
#ifndef WITHOUTMPI
  comm_buffin(1)=mass_loc
  comm_buffin(2)=ekin_loc
  comm_buffin(3)=eint_loc
  comm_buffin(4)=emag_loc
#ifdef CRFLX
  comm_buffin(5)=ecrs_loc
#endif
  call MPI_ALLREDUCE(comm_buffin,comm_buffout,ncomm,MPI_DOUBLE_PRECISION,MPI_SUM,&
       &MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(dt_loc     ,dt_all      ,1,MPI_DOUBLE_PRECISION,MPI_MIN,&
       &MPI_COMM_WORLD,info)
  mass_all=comm_buffout(1)
  ekin_all=comm_buffout(2)
  eint_all=comm_buffout(3)
  emag_all=comm_buffout(4)
#ifdef CRFLX
  ecrs_all=comm_buffout(5)
#endif
#endif
#ifdef WITHOUTMPI
  mass_all=mass_loc
  ekin_all=ekin_loc
  eint_all=eint_loc
  emag_all=emag_loc
#ifdef CRFLX
  ecrs_all=ecrs_loc
#endif
  dt_all=dt_loc
#endif

  mass_tot=mass_tot+mass_all
  ekin_tot=ekin_tot+ekin_all
  eint_tot=eint_tot+eint_all
  emag_tot=emag_tot+emag_all
#ifdef CRFLX
  ecrs_tot=ecrs_tot+ecrs_all
#endif
  dtnew(ilevel)=MIN(dtnew(ilevel),dt_all)

#ifdef CRFLX
  call update_cr_vmax(cr_vmax(ilevel))
  if(variable_cr_vmax) then
      ! Adaptively set cr_vmax to be at least N=6 times larger than gas vmax
      cr_vgas_max_all=0.0d0
#ifndef WITHOUTMPI
      call MPI_ALLREDUCE(cr_vgas_max,cr_vgas_max_all,1,MPI_DOUBLE_PRECISION,MPI_MAX,&
          & MPI_COMM_WORLD,info)
      cr_vgas_max=cr_vgas_max_all
#endif
      !cr_vmax(ilevel) = max(cr_vmax(ilevel), cr_vgas_max_all*3)
      cr_vmax(ilevel) = max(cr_vmax(ilevel), dx/3./dtnew(ilevel)*6. )
  endif
  
  ! Maximum time step for cosmic ray moment flux, from Courant condition on cr_vmax
  call get_crmom_courant(dt_cr, ilevel)
  if(mom_cr_advect) dtnew(ilevel) = MIN(dtnew(ilevel), dt_cr * cr_nsubcycle*0.99999)
#endif


111 format('   Entering courant_fine for level ',I2)

end subroutine courant_fine

!-------------------------------
! Matts edit to add energy_fine
!-------------------------------
subroutine energy_fine(ilevel)
  use amr_commons
  use hydro_commons
  use mpi_mod

  implicit none

  integer :: ilevel
  !----------------------------------------------------------
  ! This routine add heating & cooling at the bottom & top
  ! of the convection zone
  !----------------------------------------------------------
  integer::igrid,ngrid,ncache,i,ind,iskip,ix,iy,iz
  integer::nx_loc,idim,ivar
  integer::neul=5
  real(dp)::dx,dx_loc,scale,d,u,v,w,A,B,C
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:3)::skip_loc
  real(dp)::t_inject=0.15

  integer ,dimension(1:nvector),save::ind_grid,ind_cell
  real(dp),dimension(1:nvector,1:ndim),save::xx
  real(dp),dimension(1:nvector)::ee
  real(dp),dimension(1:nvector,1:nvar+3)::uut 
  real(kind=8)::uud,ekin,uuv
  real(dp),save::vol_heat=0.0d0
  real(dp),save::vol_cool=0.0d0
  logical,save::inject_cr=.true.
#ifndef WITHOUTMPI
  integer::info
  real(dp),dimension(2)::comm_buffin,comm_buffout
#endif

  if(numbtot(1,ilevel)==0)return
  if(.not. inject_cr)return
  
  ! Mesh size at level ilevel in coarse cell units
  dx=0.5D0**ilevel

  ! Rescaling factors
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=dble(nx_loc)/boxlen
  dx_loc=dx/scale

  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do

  !-------------------------------------
  ! Add analytical heating/cooling
  !-------------------------------------
  ncache=active(ilevel)%ngrid

  ! Loop over grids by vector sweeps
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do

     ! Loop over cells
     do ind=1,twotondim

        ! Gather cell indices
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do

        ! Gather cell centre positions
        do idim=1,ndim
           do i=1,ngrid
              xx(i,idim)=xg(ind_grid(i),idim)+xc(ind,idim)
           end do
        end do
        ! Rescale position from code units to user units
        do idim=1,ndim
           do i=1,ngrid
              xx(i,idim)=(xx(i,idim)-skip_loc(idim))/scale
           end do
        end do

        ! Impose analytical energy field
        if(t>t_inject)then
           call eneana(xx,ee,dx_loc,ngrid) ! Here need to ensure that modify the energy on the CR particle
        !endif
        
        ! Update total energy
          do i=1,ngrid
            uold(ind_cell(i),iCRu)=uold(ind_cell(i),iCRu)+ee(i)
          end do
        endif

     end do
     ! End loop over cells
  end do
  ! End loop over grids
  ! If first call, then set inject_cr to false
  if(t>t_inject)inject_cr=.false.
end subroutine energy_fine 
!####################################################
!####################################################
!####################################################
!####################################################


!*************************************************************************
SUBROUTINE update_cr_vmax(cr_vmax)

! Update the maximum CR speed, in code units.
! This cannot be just a constant, since scale_v changes with time in
! cosmological simulations.
! This routine probably belongs in some other .f90 though 
!-------------------------------------------------------------------------
   use constants
   use hydro_parameters,only:cr_c_fraction
   implicit none
   real(dp),intent(out)::cr_vmax
   real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
!-------------------------------------------------------------------------
   call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
   cr_vmax=c_cgs/scale_v * cr_c_fraction
END SUBROUTINE update_cr_vmax


!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine velocity_fine(ilevel)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel
  !----------------------------------------------------------
  ! This routine computes the gravitational acceleration,
  ! the maximum density rho_max, and the potential energy
  !----------------------------------------------------------
  integer::igrid,ngrid,ncache,i,ind,iskip,ix,iy,iz
  integer::nx_loc,idim,neul=5
  real(dp)::dx,dx_loc,scale,d,u,v,w,A,B,C
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:3)::skip_loc

  integer ,dimension(1:nvector),save::ind_grid,ind_cell
  real(dp),dimension(1:nvector,1:ndim),save::xx
  real(dp),dimension(1:nvector,1:3),save::vv

  if(numbtot(1,ilevel)==0)return

  ! Mesh size at level ilevel in coarse cell units
  dx=0.5D0**ilevel

  ! Rescaling factors
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=dble(nx_loc)/boxlen
  dx_loc=dx/scale

  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do

  !-------------------------------------
  ! Compute analytical velocity field
  !-------------------------------------
  ncache=active(ilevel)%ngrid

  ! Loop over grids by vector sweeps
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do

     ! Loop over cells
     do ind=1,twotondim

        ! Gather cell indices
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do

        ! Gather cell centre positions
        do idim=1,ndim
           do i=1,ngrid
              xx(i,idim)=xg(ind_grid(i),idim)+xc(ind,idim)
           end do
        end do
        ! Rescale position from code units to user units
        do idim=1,ndim
           do i=1,ngrid
              xx(i,idim)=(xx(i,idim)-skip_loc(idim))/scale
           end do
        end do

        ! Impose analytical velocity field
        call velana(xx,vv,dx_loc,t,ngrid)

        ! Impose induction variables
        do i=1,ngrid
           uold(ind_cell(i),1)=1.0
        end do
        do idim=1,3
           do i=1,ngrid
              uold(ind_cell(i),idim+1)=vv(i,idim)
           end do
        end do
        ! Update total energy
        do i=1,ngrid
           d=uold(ind_cell(i),1)
           u=uold(ind_cell(i),2)/d
           v=uold(ind_cell(i),3)/d
           w=uold(ind_cell(i),4)/d
           A=0.5*(uold(ind_cell(i),6)+uold(ind_cell(i),nvar+1))
           B=0.5*(uold(ind_cell(i),7)+uold(ind_cell(i),nvar+2))
           C=0.5*(uold(ind_cell(i),8)+uold(ind_cell(i),nvar+3))
           uold(ind_cell(i),neul)=1.0+0.5*d*(u**2+v**2+w**2)+0.5*(A**2+B**2+C**2)
        end do

     end do
     ! End loop over cells

  end do
  ! End loop over grids

end subroutine velocity_fine
!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine eneana(x,e,dx,ncell)
  use amr_parameters
  use hydro_parameters
  use poisson_parameters
  
  implicit none
  integer ::ncell                         ! Size of input arrays
  real(dp)::dx                            ! Cell size
  real(dp)::t                             ! Current time
  real(dp),dimension(1:nvector)::e        ! Energy
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.

  integer :: i, idim
  real(dp):: rx, ry, rz, x_mid, y_mid, z_mid, rr, r1  ! Radial coordinates

  ! choose value of sphere
  x_mid=0.5! ncell * dx / 2 
  y_mid=0.5! x_mid 
  z_mid=0.5! x_mid
  r1=0.10 ! radius of the sphere containing CRs

  ! Initialize
  do i=1,ncell
    e(i) = 0.0d0
  end do

  ! Loop over cells and check if inside sphere
  do i=1,ncell
     rx=x(i,1)-x_mid
#if NDIM>1
     ry=x(i,2)-y_mid
#endif
#if NDIM>2
     rz=x(i,3)-z_mid
#endif
     rr=sqrt(rx**2+ry**2+rz**2)
     if ((rr .lt. r1)) then
        e(i) = e(i) + 500! note here will eed to determine amount erg/cm^3
     end if
  end do
  
end subroutine eneana

