!************************************************************************
SUBROUTINE cr_godunov_fine(ilevel)

! This routine is a wrapper to the grid solver for radiative transfer.
! Small grids (2x2x2) are gathered from level ilevel and sent to the
! CR solver. On entry, CR variables are gathered from array uold.
! On exit, unew has been updated.
!------------------------------------------------------------------------
   use amr_commons
   use hydro_commons
   implicit none
   integer::ilevel
   integer::i,igrid,ncache,ngrid
   integer,dimension(1:nvector),save::ind_grid
!------------------------------------------------------------------------
   if(numbtot(1,ilevel)==0)return  ! # of grids at ilevel
   if(verbose)write(*,111)ilevel

   ncache=active(ilevel)%ngrid  ! total # of grids at level ilevel
   do igrid=1,ncache,nvector    ! take steps of 500 grids up to ncache
      ngrid=MIN(nvector,ncache-igrid+1) ! # of grids in each sweep
      do i=1,ngrid              ! collect grid indices for one sweep
         ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
      end do
      call cr_godfine1(ind_grid,ngrid,ilevel)
   end do

   if(verbose)write(*,112)ilevel

111 format('   Entering cr_godunov_fine for level ',i2)
112 format('   Exiting cr_godunov_fine for level ',i2)

END SUBROUTINE cr_godunov_fine

!###########################################################
!###########################################################
!###########################################################
!###########################################################
SUBROUTINE cr_set_unew(ilevel)
   use amr_commons
   use hydro_commons
   implicit none
   integer::ilevel
   !--------------------------------------------------------------------------
   ! This routine sets array unew to its initial value uold before calling
   ! the CR advection scheme. unew is set to zero in virtual boundaries.
   !--------------------------------------------------------------------------
   integer::i,ivar,ind,icpu,iskip
 
   if(numbtot(1,ilevel)==0)return
   if(verbose)write(*,111)ilevel
 
   ! Set unew to uold for myid cells
   do ind=1,twotondim
      iskip=ncoarse+(ind-1)*ngridmax
      do ivar=nvar+3+1,nvar+3+ncrvars
         do i=1,active(ilevel)%ngrid
            unew(active(ilevel)%igrid(i)+iskip,ivar) = uold(active(ilevel)%igrid(i)+iskip,ivar)
         end do
      end do
   end do
 
   ! Set unew to 0 for virtual boundary cells
   do icpu=1,ncpu
   do ind=1,twotondim
      iskip=ncoarse+(ind-1)*ngridmax
      do ivar=nvar+3+1,nvar+3+ncrvars
         do i=1,reception(icpu,ilevel)%ngrid
            unew(reception(icpu,ilevel)%igrid(i)+iskip,ivar)=0.0
         end do
      end do
   end do
   end do
 
 111 format('   Entering cr_set_unew for level ',i2)
 
END SUBROUTINE cr_set_unew
 
 !###########################################################
!###########################################################
!###########################################################
!###########################################################
SUBROUTINE cr_set_uold(ilevel)
   use amr_commons
   use hydro_commons
   use poisson_commons
   implicit none
   integer::ilevel
   !--------------------------------------------------------------------------
   ! This routine sets array uold to its new value unew after the
   ! CR step.
   !--------------------------------------------------------------------------
   integer::i,ivar,ind,iskip
 
   if(numbtot(1,ilevel)==0)return
   if(verbose)write(*,111)ilevel
  
   ! Set uold to unew for myid cells
   do ind=1,twotondim
      iskip=ncoarse+(ind-1)*ngridmax
      do ivar=nvar+3+1,nvar+3+ncrvars
         do i=1,active(ilevel)%ngrid
            uold(active(ilevel)%igrid(i)+iskip,ivar) = unew(active(ilevel)%igrid(i)+iskip,ivar)
         end do
      end do
   end do
 
 111 format('   Entering cr_set_uold for level ',i2)
 
END SUBROUTINE cr_set_uold

!*************************************************************************
SUBROUTINE cr_godfine1(ind_grid, ncache, ilevel)

! This routine first gathers CR variables from neighboring grids
! to set initial conditions in a 6x6x6 grid. It interpolates from
! coarser level to missing grid variables. It then calls the solver
! that computes fluxes. These fluxes are zeroed at coarse-fine boundaries,
! since contribution from finer levels has already been taken into
! account. Conservative variables are updated and stored in array unew(:),
! both at the current level and at the coarser level if necessary.
!
! in ind_grid: Indexes of grids/octs to solve in
! in ncache:   Length of ind_grid (i.e. number of grids)
! in ilevel:   Level at which the grids are
! in dt:       Timestep length
!-------------------------------------------------------------------------
   use amr_commons
   use hydro_commons
   use cr_flux_module
   implicit none
   integer::ilevel,ncache
   real(dp)::dt
   integer   ,dimension(1:nvector)::ind_grid
   integer   ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
   integer   ,dimension(1:nvector,1:twotondim  ),save::nbors_father_grids
   integer   ,dimension(1:nvector,0:twondim    ),save::ibuffer_father
   integer   ,dimension(1:nvector,0:twondim    ),save::ind1
   real(dp),dimension(1:nvector,0:twondim  ,1:nvar+3+ncr*(1+ndim)),save::u1
   real(dp),dimension(1:nvector,1:twotondim,1:nvar+3+ncr*(1+ndim)),save::u2

   ! 500*6*6*6*4:
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3+ncr*(1+ndim)),save::uloc
   real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:ncr*(ndim+1),1:ndim),save::flux
   logical,dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),save::ok
   integer,dimension(1:nvector),save::igrid_nbor,ind_cell,ind_buffer
   integer,dimension(1:nvector),save::ind_exist,ind_nexist

   integer::i,j,idim,ind_son,ind_father,iskip,nbuffer
   integer::i0,j0,k0,i1,j1,k1,i2,j2,k2,i3,j3,k3,nx_loc,nb_noneigh,nexist
   integer::i1min,i1max,j1min,j1max,k1min,k1max
   integer::i2min,i2max,j2min,j2max,k2min,k2max
   integer::i3min,i3max,j3min,j3max,k3min,k3max
   real(dp)::dx,scale,oneontwotondim
   integer::icr,icrE
!------------------------------------------------------------------------
   oneontwotondim = 1d0/dble(twotondim) ! 1/8 in 3D
   ! Mesh spacing
   nx_loc=icoarse_max-icoarse_min+1    ! =1
   scale=boxlen/dble(nx_loc)           ! length per coarse oct (=boxlen)
   dx=0.5D0**ilevel*scale              ! length per oct/grid at ilevel
   dt=dtnew(ilevel)

   ! Integer constants
   i1min=0; i1max=0; i2min=0; i2max=0; i3min=1; i3max=1
   j1min=0; j1max=0; j2min=0; j2max=0; j3min=1; j3max=1
   k1min=0; k1max=0; k2min=0; k2max=0; k3min=1; k3max=1
   if(ndim>0)then
      i1max=2; i2max=1; i3max=2
   end if
   if(ndim>1)then
      j1max=2; j2max=1; j3max=2
   end if
   if(ndim>2)then
      k1max=2; k2max=1; k3max=2
   end if
   ! in 3D:
   !            min      max    tot
   ! -----------------------------------------
   ! i,j,k1      0        2      3
   ! i,j,k2      0        1      2
   ! i,j,k3      1        2      2

   !------------------------------------------
   ! Gather 3^ndim neighboring father cells
   ! of grid (ilevel-1). ind_cell are indexes
   ! in uold.
   !------------------------------------------
   do i=1,ncache
      ind_cell(i)=father(ind_grid(i))
   end do
   ! ..and father cells of neighbor grids:
   call get3cubefather(ind_cell,nbors_father_cells,nbors_father_grids,ncache,ilevel)
   ! now for the parent cell (ind_cell(i)) of each grid i in cache:
   !
   ! nbors_father_cells contains indexes of all it's neighbor cells
   ! (ilevel-1), plus itself, total 3^ndim.
   !
   ! nbors_father_grids contains indexes of all the neighboring
   ! (and containing) grids of the father cell, total 2^ndim
   ! (in case interpolation is needed I guess)

   !---------------------------
   ! Gather 6x6x6 cells stencil
   !---------------------------
   ! Loop over 3x3x3 neighboring father cells (ilevel-1)
   do k1=k1min,k1max
   do j1=j1min,j1max
   do i1=i1min,i1max
      ! Check if neighbor has a grid (or is just a cell)

      ! # of neighbors (out of all the ncache cells) at k1, j1, i1
      ! that are only a cell (not a grid):
      nbuffer=0
      ! # of neighbors (out of all the ncache cells) at k1, j1, i1
      ! that have subgrids:
      nexist=0

      ind_father=1+i1+3*j1+9*k1
      do i=1,ncache
         igrid_nbor(i)=son(nbors_father_cells(i,ind_father))
         if(igrid_nbor(i)>0) then
            nexist=nexist+1
            ind_exist(nexist)=i
         else
            nbuffer=nbuffer+1
            ind_nexist(nbuffer)=i
            ind_buffer(nbuffer)=nbors_father_cells(i,ind_father)
         end if
      end do
      !--------------------------------------------------------------------
      ! we should now have: nexist+nbuffer==ncache
      !
      ! ind_exist has size nexist and contains indexes of cells in cache
      ! that have a neighboring father grid at k1, j1, i1
      !
      ! ind_nexist has size nbuffer and contains indexes of cells in cache
      ! that only have a neighboring father cell at k1, j1, i1
      !
      ! ind_buffer has size nbuffer and contains tree indexes of these
      ! father nongrid-cells
      !--------------------------------------------------------------------

      if(nbuffer>0)then
         ! For those nongrid cells we interpolate variables from parent cells:
         call getnborfather(ind_buffer,ibuffer_father,nbuffer,ilevel)
         do j=0,twondim
            do icr=1,nvar+3+ncr*(1+ndim)
               do i=1,nbuffer
                  u1(i,j,icr)= uold(ibuffer_father(i,j),icr)
               end do
            end do
            do i=1,nbuffer
               ind1(i,j)=son(ibuffer_father(i,j))
            end do
         end do
         call interpol_hydro(u1,ind1,u2,nbuffer)
      endif

      ! Loop 2x2x2 cells within father cell and add them to stencil (uloc)
      do k2=k2min,k2max
      do j2=j2min,j2max
      do i2=i2min,i2max

         ind_son=1+i2+2*j2+4*k2
         iskip=ncoarse+(ind_son-1)*ngridmax
         do i=1,nexist
            ind_cell(i)=iskip+igrid_nbor(ind_exist(i))
         end do

         i3=1; j3=1; k3=1
         if(ndim>0)i3=1+2*(i1-1)+i2 ! From -1 to 4 over outer loop, but
         if(ndim>1)j3=1+2*(j1-1)+j2 ! only over 2x2x2 indexes in inner loop
         if(ndim>2)k3=1+2*(k1-1)+k2

         ! Gather hydro variables

         do icr=1,nvar+3+ncr*(1+ndim)
            do i=1,nexist
               uloc(ind_exist(i),i3,j3,k3,icr) = uold(ind_cell(i),icr)
            end do
            do i=1,nbuffer
               uloc(ind_nexist(i),i3,j3,k3,icr) = u2(i,ind_son,icr)
            end do
         end do

         ! Gather refinement flag
         do i=1,nexist
            ok(ind_exist(i),i3,j3,k3)=son(ind_cell(i))>0
         end do
         do i=1,nbuffer
            ok(ind_nexist(i),i3,j3,k3)=.false.
         end do

      end do
      end do
      end do
      ! End loop over cells

   end do
   end do
   end do
   ! End loop over neighboring grids

   !----------------------------------------------------------------------
   ! Compute fluxes of each photon group, using Eddington tensor
   !----------------------------------------------------------------------
   do i = 1,ncr
      call cmp_cr_faces(uloc(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3+ncr*(1+ndim)) &
              ,flux,dx,dt,i,ncache,ilevel)
   end do

   !----------------------------------------------------------------------
   ! Reset flux along direction at refined interface, if not cr-subcycling
   !----------------------------------------------------------------------
   if (cr_nsubcycle == 1)then
   do idim=1,ndim
      i0=0; j0=0; k0=0
      if(idim==1)i0=1
      if(idim==2)j0=1
      if(idim==3)k0=1
      do k3=k3min,k3max+k0
      do j3=j3min,j3max+j0
      do i3=i3min,i3max+i0
         do icr=1,ncr
            icrE = 1+(icr-1)*(ndim+1)
            do i=1,ncache
               if(ok(i,i3-i0,j3-j0,k3-k0) .or. ok(i,i3,j3,k3))then
                  flux(i,i3,j3,k3,icrE:icrE+ndim,idim)=0.0d0
               end if
            end do
         end do
      end do
      end do
      end do
   end do
   endif
   !--------------------------------------
   ! Conservative update at level ilevel
   !--------------------------------------
   do idim=1,ndim
      i0=0; j0=0; k0=0
      if(idim==1)i0=1
      if(idim==2)j0=1
      if(idim==3)k0=1
      do k2=k2min,k2max      ! all from 0 to 1
      do j2=j2min,j2max      ! => update 2x2x2 cells
      do i2=i2min,i2max      ! i.e. update one grid at level ilevel
         ind_son=1+i2+2*j2+4*k2
         iskip=ncoarse+(ind_son-1)*ngridmax
         do i=1,ncache
            ind_cell(i)=iskip+ind_grid(i)
         end do
         i3=1+i2
         j3=1+j2   ! just because the flux indexes are (1:3), not (0:2)
         k3=1+k2
         ! Update conservative variables new state vector
         do icr=1,ncr
            icrE = 1+(icr-1)*(ndim+1) ! starting index of cr group in u1,u2,uloc,flux
            do i=1,ncache
               unew(ind_cell(i),icrU:icrU+ncr*(ndim+1)-1)= &
                  unew(ind_cell(i),icrU:icrU+ncr*(ndim+1)-1)   &
                  & +(flux(i,i3   ,j3   ,k3   ,icrE:icrE+ndim,idim)  &
                  & - flux(i,i3+i0,j3+j0,k3+k0,icrE:icrE+ndim,idim))
            end do
         end do

      end do
      end do
      end do
   end do

   !--------------------------------------
   ! Conservative update at level ilevel-1
   !--------------------------------------
   if (cr_nsubcycle == 1)then
      ! Loop over dimensions
      do idim=1,ndim
         i0=0; j0=0; k0=0
         if(idim==1)i0=1
         if(idim==2)j0=1
         if(idim==3)k0=1

         !----------------------
         ! Left flux at boundary
         !----------------------
         ! Check if grids sits near left boundary
         ! and gather neighbor father cells index
         nb_noneigh=0
         do i=1,ncache
            if (son(nbor(ind_grid(i),2*idim-1))==0) then
               nb_noneigh = nb_noneigh + 1
               ind_buffer(nb_noneigh) = nbor(ind_grid(i),2*idim-1)
               ind_cell(nb_noneigh) = i
            end if
         end do
         ! Conservative update of new state variables
         do icr=1,ncr
            icrE = 1+(icr-1)*(ndim+1) ! starting index of cr group in u1,u2,uloc,flux
            ! Loop over boundary cells
            do k3=k3min,k3max-k0 ! 1 to 1 if dim=3, 1 to 2 otherwise
               do j3=j3min,j3max-j0 ! 1 to 1 if dim=2, 1 to 2 otherwise
                  do i3=i3min,i3max-i0 ! 1 to 1 if dim=1, 1 to 2 otherwise
                     do i=1,nb_noneigh
                        unew(ind_buffer(i),icrU:icrU+ncr*(ndim+1)-1)= &
                           unew(ind_buffer(i),icrU:icrU+ncr*(ndim+1)-1)   &
                           & -flux(ind_cell(i),i3   ,j3   ,k3   ,icrE:icrE+ndim,idim)  &
                           & * oneontwotondim
                    end do
                  end do
               end do
            end do
         end do

         !-----------------------
         ! Right flux at boundary
         !-----------------------
         ! Check if grids sits near right boundary
         ! and gather neighbor father cells index
         nb_noneigh=0
         do i=1,ncache
            if (son(nbor(ind_grid(i),2*idim))==0) then
               nb_noneigh = nb_noneigh + 1
               ind_buffer(nb_noneigh) = nbor(ind_grid(i),2*idim)
               ind_cell(nb_noneigh) = i
            end if
         end do
         ! Conservative update of new state variables
         do icr=1,ncr
            icrE = 1+(icr-1)*(ndim+1) ! starting index of cr group in u1,u2,uloc,flux
            ! Loop over boundary cells
            do k3=k3min+k0,k3max
               do j3=j3min+j0,j3max
                  do i3=i3min+i0,i3max
                     do i=1,nb_noneigh
                        unew(ind_buffer(i),icrU:icrU+ncr*(ndim+1)-1)= &
                           unew(ind_buffer(i),icrU:icrU+ncr*(ndim+1)-1)   &
                           & +flux(ind_cell(i),i3+i0   ,j3+j0   ,k3+k0   ,icrE:icrE+ndim,idim)  &
                           & * oneontwotondim
                     end do
                  end do
               end do
            end do
         end do

      end do
      ! End loop over dimensions
   end if
   ! End if-clause for cr-subcycling

END SUBROUTINE cr_godfine1

 
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine add_cr_source_terms(ilevel)
  use amr_commons
  use hydro_commons
  use cr_flux_module
  implicit none
  integer::ilevel
  !---------------------------------------------------------
  ! This routine adds the cosmic ray source term .
  !---------------------------------------------------------
  integer::i,ind,iskip,nx_loc
  integer::ncache,igrid,ngrid,idim,id1,ig1,ih1,id2,ig2,ih2
  integer,dimension(1:3,1:2,1:8)::iii,jjj
  real(dp)::scale,dx,dx_loc,B_field(3),dt

  integer ,dimension(1:nvector),save::ind_grid,ind_cell
  integer ,dimension(1:nvector,0:twondim),save::igridn
  integer ,dimension(1:nvector,1:ndim),save::ind_left,ind_right
  real(dp),dimension(1:nvector,1:ndim),save::dx_g,dx_d
  real(dp),dimension(1:nvector,1:ndim,1:ncr),save::gradpcr_loc
  real(dp),dimension(1:nvector,1:3),save::B_field_loc
  real(dp),dimension(1:nvector),save::bdotgradp_loc
  real(dp),dimension(1:nvector,1:3),save::vs_loc
  real(dp),dimension(1:nvector),save::va_loc
  real(dp)::norm,frotx,froty,frotz,bxby,cosp,sinp,cost,sint
  real(dp)::f1,f2,f3
  integer::j,icr,icrE
  real(dp),dimension(1:nvector,1:ndim,1:ncr),save::pcrg,pcrd
  real(dp)::coef_11, coef_12, coef_13, coef_14, coef_21, coef_22
  real(dp)::coef_31, coef_33, coef_41, coef_44
  real(dp)::e_coef, new_ec, old_ec, sigma_x, sigma_y, sigma_z, sigma_stream
  real(dp)::rhs1, rhs2, rhs3, rhs4, fred
  real(dp)::v1, v2, v3, vtot1, vtot2, vtot3
  real(dp)::scale_l, scale_t, scale_d, scale_v, scale_nH, scale_T2
  real(dp)::scale_kappa, DCR_code, DCRmax_code, mom_change
  real(dp)::etherm, ekin, emag, err, f_decouple, smallp

#if NENER>0
  integer::irad
#endif

  dt=dtnew(ilevel)

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  smallp = smallc**2/gamma

  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx=0.5d0**ilevel
  dx_loc=dx*scale

  iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/); jjj(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(1,2,1:8)=(/0,2,0,2,0,2,0,2/); jjj(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(2,1,1:8)=(/3,3,0,0,3,3,0,0/); jjj(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(2,2,1:8)=(/0,0,4,4,0,0,4,4/); jjj(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(3,1,1:8)=(/5,5,5,5,0,0,0,0/); jjj(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  iii(3,2,1:8)=(/0,0,0,0,6,6,6,6/); jjj(3,2,1:8)=(/5,6,7,8,1,2,3,4/)

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_kappa = scale_l**2/scale_t
  DCRmax_code=DCRmax/scale_kappa
  DCR_code=DCR/scale_kappa

  ! Loop over myid grids by vector sweeps
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector

     ! Gather nvector grids
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do

     ! Gather neighboring grids
     do i=1,ngrid
        igridn(i,0)=ind_grid(i)
     end do
     do idim=1,ndim
        do i=1,ngrid
           ind_left (i,idim)=nbor(ind_grid(i),2*idim-1)
           ind_right(i,idim)=nbor(ind_grid(i),2*idim  )
           igridn(i,2*idim-1)=son(ind_left (i,idim))
           igridn(i,2*idim  )=son(ind_right(i,idim))
        end do
     end do

     ! Loop over cells
     do ind=1,twotondim

        ! Compute central cell index
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do

        ! Gather all neighboring CR energies
        do idim=1,ndim
           id1=jjj(idim,1,ind); ig1=iii(idim,1,ind)
           ih1=ncoarse+(id1-1)*ngridmax
           do i=1,ngrid
           do icr=1,ncr
              icrE = icru+(ndim+1)*(icr-1)  ! starting index of cr variables
              if(igridn(i,ig1)>0)then
                  pcrg(i,idim,icr)   = max(uold(igridn(i,ig1)+ih1,icrE),smallcr)
                  dx_g(i,idim) = dx_loc
              else
                  pcrg(i,idim,icr)   = max(uold(ind_left(i,idim),icrE),smallcr)
                  dx_g(i,idim) = dx_loc*1.5_dp
              end if
           enddo
           enddo
           id2=jjj(idim,2,ind); ig2=iii(idim,2,ind)
           ih2=ncoarse+(id2-1)*ngridmax
           do i=1,ngrid
           do icr=1,ncr
              icrE = icru+(ndim+1)*(icr-1)  ! starting index of cr variables
              if(igridn(i,ig2)>0)then
                  pcrd(i,idim,icr)  = max(uold(igridn(i,ig2)+ih2,icrE),smallcr)
                  dx_d(i,idim)=dx_loc
              else
                  pcrd(i,idim,icr)  = max(uold(ind_right(i,idim),icrE),smallcr)
                  dx_d(i,idim)=dx_loc*1.5_dp
              end if
           enddo
           enddo
        end do
        ! End loop over dimensions
        
        do i=1,ngrid
        do icr=1,ncr
           do idim=1,ndim
              gradpcr_loc(i,idim,icr) = (pcrd(i,idim,icr)-pcrg(i,idim,icr)) &
              &                    / (dx_g(i,idim)     +dx_d(i,idim)) * (gamma_cr(icr)-1.)
           enddo
           vs_loc(i,:)=0.
           va_loc(i)=0.
           do idim=1,3
              B_field_loc(i,idim) = 0.5 * (uold(ind_cell(i), 5+idim) + uold(ind_cell(i), nvar+idim) )
              vs_loc(i,idim) = B_field_loc(i,idim)/sqrt(uold(ind_cell(i),1))
              va_loc(i) = va_loc(i) + vs_loc(i,idim)**2
           enddo
           va_loc(i) = sqrt(va_loc(i))

           if(v_alfven.gt.0.) then
              vs_loc(i,:) = 0. 
              vs_loc(i,1) = v_alfven
              va_loc(i) = v_alfven
           endif

           bdotgradp_loc(i)=0.
           ! bdotgradp is needed for eq 3 in Jiang & Oh
           do idim=1,ndim
              bdotgradp_loc(i) = bdotgradp_loc(i) + B_field_loc(i,idim)*gradpcr_loc(i,idim,icr)
           enddo
           ! Local streaming velocity, needed for eq 18 in Jiang & Oh
           vs_loc(i,:) = -vs_loc(i,:)  * bdotgradp_loc(i)/max(1d-50,abs(bdotgradp_loc(i))) !! eq 3 in PO17
        end do ! ncr
        end do ! ngrid

        ! Source term, eq 18 in Jiang & Oh 2017
        do i=1,ngrid
        do icr=1,ncr
            icrE = icru+(ndim+1)*(icr-1)  ! starting index of cr variables
            norm=0.
            do j=1,3
               B_field(j) = 0.5 * (uold(ind_cell(i), 5+j) + uold(ind_cell(i), nvar+j) )
               norm = norm + B_field(j)**2
            end do
            norm = max(sqrt(norm), 1d-30)

            bxby = sqrt(B_field(1)**2+B_field(2)**2)
            if(norm.gt.1e-10) then
               sint = bxby/norm     
               cost = B_field(3)/norm
            else
               sint = 1d0
               cost = 0d0
            endif
            if(bxby.gt.1e-10) then
               sinp = B_field(2)/bxby     
               cosp = B_field(1)/bxby
            else
               sinp = 0d0
               cosp = 1d0
            endif
            !B_field = B_field/norm
 
            rhs1 = max(unew(ind_cell(i),icrE),smallcr)
  
            ! Flux update ... first rotate flux vector onto B-field,
            ! i.e. describing the flux in the B coordinate system
            f2=0. ; f3=0.
            f1 = unew(ind_cell(i),icrE+1)
#if NDIM>1
            f2 = unew(ind_cell(i),icrE+2)
#endif
#if NDIM>2
            f3 = unew(ind_cell(i),icrE+3)
#endif

            ! Make sure that |F|<=cE
            fred = sqrt(f1**2+f2**2+f3**2)/(cr_vmax(ilevel)*rhs1)
            if(fred .gt. 1d0) then
               f1 = f1/fred ; f2 = f2/fred ; f3 = f3/fred
            endif

            frotx=f1; froty=f2; frotz=f3
            call rotatevec(sint, cost, sinp, cosp, frotx, froty, frotz)
            v1=0. ; v2=0. ; v3=0.
            v1 = uold(ind_cell(i),2)/uold(ind_cell(i),1)
#if NDIM>1
            v2 = uold(ind_cell(i),3)/uold(ind_cell(i),1)
#endif
#if NDIM>2
            v3 = uold(ind_cell(i),4)/uold(ind_cell(i),1)
#endif
            vtot1 = v1 ; vtot2 = v2 ; vtot3 = v3
            if(mom_streaming_heating) then
               vtot1 = vtot1+vs_loc(i,1)
               vtot2 = vtot2+vs_loc(i,2)
               vtot3 = vtot3+vs_loc(i,3)
            endif

            ! Rotate velocity
            call rotatevec(sint, cost, sinp, cosp, v1, v2, v3)          
            ! Rotate streaming velocity
            call rotatevec(sint, cost, sinp, cosp, vtot1, vtot2, vtot3)

            rhs2 = frotx
            rhs3 = froty
            rhs4 = frotz
            if(abs(rhs2).lt.1e-20) rhs2=0d0
            if(abs(rhs3).lt.1e-20) rhs3=0d0
            if(abs(rhs4).lt.1e-20) rhs4=0d0

            ! Factor for decoupling CRs from gas at low densities
            f_decouple = MAX(exp(-smallr*cr_smallr_decouple/uold(ind_cell(i),1)),1d-10)

            sigma_stream = max(1./DCRmax_code &
               & ,abs(bdotgradp_loc(i)) / norm / va_loc(i) / gamma_cr(icr) / max(unew(ind_cell(i),icrE),smallcr))
            if(mom_streaming_diffusion) then
               sigma_x = 1./(Dcr_code + 1./sigma_stream)
            else
               sigma_x = 1./Dcr_code
            endif
            sigma_y = 1./Dcr_code*1e6
            sigma_z = 1./Dcr_code*1e6

            coef_11 = 1.0 - dt * sigma_x * vtot1 * v1 * gamma_cr(icr)   &
                         - dt * sigma_y * vtot2 * v2 * gamma_cr(icr)   &
                         - dt * sigma_z * vtot3 * v3 * gamma_cr(icr) 
            coef_12 = dt * sigma_x * vtot1
            coef_13 = dt * sigma_y * vtot2
            coef_14 = dt * sigma_z * vtot3
   
            coef_21 = -dt * v1 * sigma_x * gamma_cr(icr) * cr_vmax(ilevel)**2
            coef_22 = 1.0 + dt * sigma_x * cr_vmax(ilevel)**2

            coef_31 = -dt * v2 * sigma_y * gamma_cr(icr) * cr_vmax(ilevel)**2
            coef_33 = 1.0 + dt * sigma_y * cr_vmax(ilevel)**2
   
            coef_41 = -dt * v3 * sigma_z * gamma_cr(icr) * cr_vmax(ilevel)**2
            coef_44 = 1.0 + dt * sigma_z * cr_vmax(ilevel)**2
   
            ! newfr1 = (rhs2 - coef21 * newEc)/coef22
            ! newfr2= (rhs3 - coef31 * newEc)/coef33
            ! newfr3 = (rhs4 - coef41 * newEc)/coef44
            ! coef11 - coef21 * coef12 /coef22 - coef13 * coef31 /coef33 - coef41 * coef14 /coef44)* newec
            !    =rhs1 - coef12 *rhs2/coef22 - coef13 * rhs3/coef33 - coef14 * rhs4/coef44
   
            e_coef = coef_11 - coef_12 * coef_21/coef_22 - coef_13 * coef_31/coef_33 &
                            - coef_14 * coef_41/coef_44
            new_ec = rhs1 - coef_12 * rhs2/coef_22 - coef_13 * rhs3/coef_33 &
                         - coef_14 * rhs4/coef_44
            new_ec = new_ec / e_coef

            old_ec = unew(ind_cell(i),icrE)
            unew(ind_cell(i),icrE) = unew(ind_cell(i),icrE) + (new_ec-old_ec) * f_decouple

            ! Floor the CR energy and update total energy if necessary
            if ( unew(ind_cell(i),icrE) .lt. smallcr ) unew(ind_cell(i),icrE) = smallcr
            ! Thermal energy update:
            if(.not. static_gas) then
               unew(ind_cell(i),5) = unew(ind_cell(i),5) - (unew(ind_cell(i),icrE) - old_ec)*f_decouple
               unew(ind_cell(i),5) = max(smallp*uold(ind_cell(i),1), unew(ind_cell(i),5))
            endif
            frotx = (rhs2 - coef_21 * new_ec)/coef_22
            froty = (rhs3 - coef_31 * new_ec)/coef_33
            frotz = (rhs4 - coef_41 * new_ec)/coef_44

            ! Make sure that |F|<=cE
            fred = sqrt(frotx**2+froty**2+frotz**2)/(cr_vmax(ilevel)*new_ec)
            if(fred .gt. 1d0) then
               !print*,'maybe a problem with fred'
               frotx = frotx/fred ; froty = froty/fred ; frotz = frotz/fred
            endif

            ! We are missing the perpendicular energy source term which is done here in Athena!!!
            ! Rotate the flux back to the simulation coordinate system
            call invrotatevec(sint, cost, sinp, cosp, frotx, froty, frotz)

            unew(ind_cell(i),icrE+1) = frotx
            ! Momentum update
            mom_change = ( f1 - unew(ind_cell(i),icrE+1) ) / cr_vmax(ilevel)**2
            if(.not. static_gas) &
               unew(ind_cell(i),2) = unew(ind_cell(i),2) + mom_change*f_decouple
#if NDIM>1
            unew(ind_cell(i),icrE+2) = froty
            mom_change = ( f2 - unew(ind_cell(i),icrE+2) ) / cr_vmax(ilevel)**2
            if(.not. static_gas) &
               unew(ind_cell(i),3) = unew(ind_cell(i),3) + mom_change*f_decouple
#endif
#if NDIM>2
            unew(ind_cell(i),icrE+3) = frotz
            mom_change = ( f3 - unew(ind_cell(i),icrE+3) ) / cr_vmax(ilevel)**2
            if(.not. static_gas) &
               unew(ind_cell(i),4) = unew(ind_cell(i),4) + mom_change*f_decouple
#endif

! here add the cr routine from cooling_fine
! update this before next loop

        end do ! ncr 

            ekin=0d0
            do idim=1,ndim
               ekin = ekin + 0.5d0*unew(ind_cell(i),idim+1)**2/unew(ind_cell(i),1)
            end do
            err=0d0
#if NENER>0
            do irad=0,nener-1
               err = err + unew(ind_cell(i),inener+irad)
            end do
#endif
            emag=0d0
            do idim=1,ndim
               emag = emag + 0.125d0*(unew(ind_cell(i),idim+5)+unew(ind_cell(i),idim+nvar))**2
            end do

            etherm = unew(ind_cell(i),5) - ekin - emag - err
            if(etherm .lt. smallp*uold(ind_cell(i),1)) then
               unew(ind_cell(i),5) = ekin + emag + err + smallp*uold(ind_cell(i),1)
               if(myid.eq.1) print*,'Noooooo negative energy!!'
               !call clean_stop
            endif

        end do ! ncache

     enddo
     ! End loop over cells
  end do
  ! End loop over grids

111 format('   Entering add_cr_source_terms for level ',i2)

end subroutine add_cr_source_terms

!##########################################################################
!##########################################################################
!##########################################################################
!##########################################################################

subroutine crmom_step(ilevel)
  use amr_parameters, only: dp, simple_boundary, ncontrol
  use amr_commons,    only: t, dtnew, myid, nstep_coarse
  use hydro_parameters, only:nvar,ncrvars,cr_nsubcycle,cr_vmax,cr_vgas_max
  use hydro_commons, only: unew, uold
  use constants
  use mpi_mod
  implicit none
  integer, intent(in) :: ilevel
  real(dp) :: dt_hydro, t_left, dt_cr, t_save
  integer  :: i_substep, ivar
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  !------------------------------------------------------------------------
  !  CR advection step. Either do one step on ilevel,
  !  with CR field updates in coarser level neighbours, or, if
  !  cr_nsubcycle>1, do many substeps in ilevel only, using Dirichlet
  !  boundary conditions for the level boundaries.
  !------------------------------------------------------------------------

  dt_hydro = dtnew(ilevel)                   ! Store hydro timestep length
  t_left = dt_hydro
  ! We shift the time backwards one hydro-dt, in the case of rt subcycling.
  ! This is not really necessary as is though, since there is no dependency 
  ! in CR advection on the simulation time.
  t_save=t ; t=t-t_left

  call get_crmom_courant(dt_cr,ilevel)
  i_substep = 0
  do while (t_left > 0)                      !                RT sub-cycle
     i_substep = i_substep + 1
     ! Temporarily change timestep length to cr step:
     dtnew(ilevel) = MIN(t_left, dt_cr)
     if(i_substep .gt. cr_nsubcycle) then
         print*,'Doing CR substeps but should not!! ',i_substep, cr_nsubcycle, dt_hydro, t_left
         call clean_stop
     endif
     t = t + dtnew(ilevel) ! Shift the time forwards one dt_rt

     if (i_substep > 1) call cr_set_unew(ilevel)

     call cr_godunov_fine(ilevel)
     ! Pull in updates from finer level on other cpus -- only needed when updating coarser level
     if(cr_nsubcycle==1) then
         do ivar=nvar+3+1,nvar+3+ncrvars
            call make_virtual_reverse_dp(unew(1,ivar),ilevel)
         end do
     endif
     call add_cr_source_terms(ilevel)

     ! Set uold equal to unew, but only for CR vars
     call cr_set_uold(ilevel)

     do ivar=1,nvar+3+ncrvars
        call make_virtual_fine_dp(uold(1,ivar),ilevel)
     end do

     ! Should update this to cr_make_boundary_hydro
     if(simple_boundary)call make_boundary_hydro(ilevel)

     t_left = t_left - dtnew(ilevel)
  end do                                   !          End CR subcycle loop
  dtnew(ilevel) = dt_hydro                 ! Restore hydro timestep length
  t = t_save       ! Restore original time (otherwise tiny roundoff error)

  ! Restriction operator to update this level split cells
  call upload_fine(ilevel)

  if (myid==1 .and. mod(nstep_coarse,ncontrol)==0) then
      call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
      write(*,901) ilevel, i_substep, cr_vmax(ilevel)*scale_v/1e5, cr_vmax(ilevel)*scale_v/c_cgs, dt_cr
      !print*,'The max gas velocity in km/s is ',cr_vgas_max_tmp*65.59266058737735
   endif

901 format (' Performed level', I3, ' CR-step with ', I5, ' subcycles, cr_vmax(km/s)=', 1pe9.2, ' cr_c_fraction=', 1pe9.2, ',  dt_cr=', 1pe9.2)

end subroutine crmom_step

!###########################################################
!###########################################################
!###########################################################                                     
SUBROUTINE get_crmom_courant(dt,ilevel)

   ! Determine the coarse CR timestep length set by the Courant condition
   !-------------------------------------------------------------------------
     use amr_parameters
     use hydro_parameters
     implicit none
     integer:: nx_loc,ilevel
     real(dp):: dt, scale, dx
   !-------------------------------------------------------------------------
     ! Mesh spacing at coarse level
     nx_loc=icoarse_max-icoarse_min+1
     scale=boxlen/dble(nx_loc)
     dx=0.5D0**levelmin*scale
     dt=courant_factor*dx/3d0/cr_vmax(ilevel)/2.0**(ilevel-levelmin)
END SUBROUTINE get_crmom_courant
