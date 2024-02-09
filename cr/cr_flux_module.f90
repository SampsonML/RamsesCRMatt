MODULE cr_flux_module

  use amr_commons
  use hydro_commons
  use hydro_parameters
  implicit none

  private   ! default

  public cmp_cr_faces, rotatevec, invrotatevec

CONTAINS

!************************************************************************
subroutine cmp_cr_flux_tensors(uin, iGrp, nGrid, ftens, vmax)
  
! Compute central fluxes for a CR group, for each cell in a vector 
! of grids. 
! The flux tensor is a three by four tensor (2*3 and 1*2 in 1D and 2D, 
! respectively) where the first row is CR flux (x,y,z) and 
! the other three rows compose the Eddington tensor (see Yiang&Peng 2017)
! input/output:
! uin       => uold variables of all cells in a vector of grids
! igrp      => CR group number
! ngrid     => Number of 'valid' grids in uin.
! ftens     <=  Group flux tensors for all the cells.
!------------------------------------------------------------------------
  real(dp),dimension(nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3+ncrvars)::uin
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nDim+1,1:ndim)::ftens
  integer::iGrp, nGrid!---------------------------------------------------
  real(dp),dimension(1:ndim)::crflux
  real(dp)::Ecr, norm, vmax
  integer::i, j, k, idim, jdim, n
!------------------------------------------------------------------------

  ! Loop 6X6X6 cells in grid, from -1 to 4.
  do k = ku1, ku2
  do j = ju1, ju2
  do i = iu1, iu2
    
     do n=1,ngrid
        Ecr =   uin(n, i, j, k, icrU+(iGrp-1)*(ndim+1))  ! CR density in cell
        crflux = uin(n,i,j,k,icrU+1+(iGrp-1)*(ndim+1) : icrU+1+(iGrp-1)*(ndim+1)+ndim-1) ! CR flux vector
        norm=0.

        if(Ecr .lt. 0d0) then
          write(*,*)'negative CR density in cmp_flux_tensors. -EXITING-'
          call clean_stop
        endif
        ftens(n,i,j,k,1,1:nDim)= crflux  !   First row is CR flux
        ! Rest is Eddington tensor
        if(isotropic_pressure) then
          ftens(n,i,j,k,2:ndim+1,1:nDim) = 0d0
          do idim = 1, ndim
            ftens(n,i,j,k,idim+1,idim) = &
              (gamma_cr(iGrp)-1.)*Ecr*vmax**2
          enddo
        else
          do idim = 1, ndim
            do jdim = 1, ndim
              ftens(n,i,j,k,idim+1,jdim) =          &
                (gamma_cr(iGrp)-1.)* Ecr * vmax**2
            end do
          end do
        endif
    enddo
  enddo
  enddo
  enddo

end subroutine cmp_cr_flux_tensors

!************************************************************************
SUBROUTINE cmp_cr_wavespeeds(uin, iGrp, ngrid, lmax, ilevel)

  !  Compute CR wavespeeds for given vector of sub-grids.
  !
  !  inputs/outputs
  !  uin         => input cell states
  !  iGrp        => CR group number
  !  ngrid       => number of sub-grids of 3^ndim cells
  !  lmax       <=  return maximum cell wavespeeds
  !
  !  other vars
  !  iu1,iu2     |first and last index of input array,
  !  ju1,ju2     |cell centered,
  !  ku1,ku2     |including buffer cells.
  !------------------------------------------------------------------------
    real(dp),dimension(nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3+ncrvars),  &
                                                            intent(in)::uin 
    real(dp),dimension(nvector,iu1:iu2,ju1:ju2,ku1:ku2,ndim)::   lmax
    integer,intent(in)::iGrp, ngrid!---------------------------------------
    real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
    real(dp)::scale_kappa, dx, dx_loc, scale, Ecr, va, Dcrmax_code
    integer::iP0, i, j, k, n, ilevel, nx_loc, idim
    real(dp),dimension(1:3)::skip_loc, B_field, gradpcr
    real(dp)::norm,bdotgradp,cosp,sinp,cost,sint,bxby,Dcr_dir,Dcr_vec(3)
  !------------------------------------------------------------------------
    iP0 = icrU+(iGrp-1)*(ndim+1) ! Index of CR group energy density
    ! Conversion factor from user units to cgs units
    call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
    scale_kappa = scale_l**2/scale_t
    DCRmax_code = DCRmax/scale_kappa
  
    ! Cell width in ilevel
    dx=0.5D0**ilevel
    nx_loc=(icoarse_max-icoarse_min+1)
    skip_loc=(/0.0d0,0.0d0,0.0d0/)
    if(ndim>0)skip_loc(1)=dble(icoarse_min)
    if(ndim>1)skip_loc(2)=dble(jcoarse_min)
    if(ndim>2)skip_loc(3)=dble(kcoarse_min)
    scale=boxlen/dble(nx_loc)
    dx_loc=dx*scale
    do k=kfcr1,kf2                                !
    do j=jfcr1,jf2                                !  Loop each cell in grid
    do i=ifcr1,if2                                !
      
    do n=1,ngrid                                     ! Loop buffer of grids
  
       Ecr    = uin(n, i, j, k, iP0)
  
       ! Magnetic field, needed to rotate Dcr and for bdotgradPcr
       norm=0.
       do idim=1,3
          B_field(idim) = 0.5 * &
            (uin(n, i, j, k, 5+idim) + uin(n, i, j, k, nvar+idim) )
          norm = norm + B_field(idim)**2
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
       B_field = B_field/norm

       va=0.
       if(mom_streaming_diffusion) va=norm/sqrt(uin(n,i,j,k,1))
       if(mom_streaming_diffusion .and. v_alfven.gt.0.0) va = v_alfven

       ! Calculate grad Pcr
       gradpcr(1) = (uin(n, i+1, j, k, iP0)-uin(n, i-1, j, k, iP0)) &
                      * (gamma_cr(iGrp)-1.)  / (2d0*dx_loc)
#if NDIM>1
       gradpcr(2) = (uin(n, i, j+1, k, iP0)-uin(n, i, j-1, k, iP0)) &
                      * (gamma_cr(iGrp)-1.)  / (2d0*dx_loc)
#endif
#if NDIM>2
       gradpcr(3) = (uin(n, i, j, k+1, iP0)-uin(n, i, j, k-1, iP0)) &
                      * (gamma_cr(iGrp)-1.)  / (2d0*dx_loc)
#endif

       ! Calculate B dot grad Pcr
       bdotgradp = 1e-30
       do idim=1,ndim
          bdotgradp = bdotgradp + B_field(idim) * gradpcr(idim)
       enddo

       ! Diffusion, eq 10 in JO17
       Dcr_vec = (/ Dcr, Dcr/1d6, Dcr/1d6 /)
       Dcr_vec = Dcr_vec/scale_kappa ! convert to code units
       if(mom_streaming_diffusion) &
          Dcr_vec(1) = Dcr_vec(1) + &
              min(DCRmax_code, 1./abs(bdotgradp) * va * gamma_cr(iGrp) * max(Ecr,smallcr))

       ! Rotate Dcr_vec so it is parallel with B, hence
       ! describing Dcr in the simulation coordinate system
       ! (instead of the B coordinate system)
       call invrotatevec(sint, cost, sinp, cosp, Dcr_vec(1), Dcr_vec(2), Dcr_vec(3))
  
       ! Calculate wavespeeds
       Dcr_dir = abs(Dcr_vec(1)) ! x component of rotated Dcr
       lmax(n,i,j,k,1) = cmp_cr_lmax(dx_loc, gamma_cr(iGrp), Dcr_dir, cr_vmax(ilevel))

#if NDIM>1
       Dcr_dir = abs(Dcr_vec(2)) ! y component of rotated Dcr
       lmax(n,i,j,k,2) = cmp_cr_lmax(dx_loc, gamma_cr(iGrp), Dcr_dir, cr_vmax(ilevel))
#endif
  
#if NDIM>2
       Dcr_dir = abs(Dcr_vec(3)) ! z component of rotated Dcr
       lmax(n,i,j,k,3) = cmp_cr_lmax(dx_loc, gamma_cr(iGrp), Dcr_dir, cr_vmax(ilevel))
#endif
  
    end do
    end do
    end do
    end do
  
END SUBROUTINE cmp_cr_wavespeeds
  
!************************************************************************
FUNCTION cmp_cr_lmax(dx_loc, gamma_rad_loc, dcoeff, vmax)
  
! Compute maximum local wavespeed
!------------------------------------------------------------------------
  real(dp)::dx_loc, gamma_rad_loc, cmp_cr_lmax, vmax
  real(dp)::tau, dcoeff, r_factor
!------------------------------------------------------------------------
  tau = cr_f_taucell * dx_loc / dcoeff * vmax * sqrt(1.5)
  if(tau.lt.1e-3) then
    r_factor = sqrt((1.0 - 0.5*tau**2))
  else
    r_factor = sqrt((1.-exp(-min(tau,10.)**2))/min(tau,1e8)**2) ! Capital R on p 6 in YP17
  endif
  cmp_cr_lmax = r_factor * vmax * sqrt(gamma_rad_loc-1.)

END FUNCTION cmp_cr_lmax

!************************************************************************
FUNCTION cmp_cr_face(fdn, fup, udn, uup, lminus, lplus)
  
! Compute HLLE intercell fluxes for all (four) CR variables.
! fdn    => flux function in the cell downwards from the border
! fup    => flux function in the cell upwards from the border
! udn    => state (cr density and flux) downwards from the border
! uup    => state (cr density and flux) upwards from the border
! lminus => minimum intercell wavespeed
! lplus  => maximum intercell wavespeed
! returns      flux vector for the given state variables, i.e. line nr dim
!              in the 3*4 flux function tensor
!------------------------------------------------------------------------
  real(dp),dimension(nDim+1)::fdn, fup, udn, uup, cmp_cr_face
  real(dp)::lminus, lplus, coeff
!------------------------------------------------------------------------
  coeff=0d0
  if(abs(lplus-lminus).gt.1d-20) coeff = 0.5d0*(lplus+lminus)/(lplus-lminus)
  cmp_cr_face = 0.5d0 * (fdn+fup - lminus*udn - lplus*uup) &
              +   coeff * (fdn-fup - lminus*udn + lplus*uup  )
END FUNCTION cmp_cr_face

!************************************************************************
FUNCTION minmod(left, right)

! Minmod interpolation function for intercell fluxes  
  real(dp),dimension(ndim+1)::left,right,minmod
  integer::i
!------------------------------------------------------------------------
  do i=1,ndim+1
    if (abs(left(i) ) .gt. abs(right(i))) minmod(i)=right(i)
    if (abs(right(i)) .ge. abs(left(i) )) minmod(i)=left(i)
    if (left(i)*right(i) .le. 0.0d0     ) minmod(i)=0.0d0
  end do
END FUNCTION minmod

!************************************************************************
FUNCTION vminmod(left, right)

! Minmod interpolation function for velocity 
  real(dp)::left,right,vminmod
!------------------------------------------------------------------------
  if (abs(left ) .gt. abs(right)) vminmod=right
  if (abs(right) .ge. abs(left )) vminmod=left
  if (left*right .le. 0.0d0     ) vminmod=0.0d0
END FUNCTION vminmod

!************************************************************************
SUBROUTINE cmp_cr_faces(uin,iFlx,dx,dt,iGrp,ngrid,ilevel)

!  Compute intercell fluxes for one CR group in all dimensions,
!  using the Eddington tensor with the Yiang+Peng'17 closure relation.
!  The intercell fluxes are the right-hand sides of the equation:
!      dq/dt = - nabla dot f   (eq 12 in YP17)
!  where q=[Ecr, Fx/ccr^2, Fy/ccr^2, Fz/ccr^2], ccr the reduced wavespeed
!  and f the Eddington pressure tensor. A flux at index i,j,k represents
!  flux across the lower faces of that cell, i.e. at i-1/2 etc.
!
!  inputs/outputs
!  uin         => input states
!  iFlx       <=  return fluxes in the 3 coord directions.
!  dx          => cell width
!  dt          => time step
!  iGrp        => CR group number
!  ngrid       => number of sub-grids
!  ilevel      => level being updated
!
!  other vars
!  iu1,iu2     |First and last index of input array,
!  ju1,ju2     |cell centered,
!  ku1,ku2     |including buffer cells (6x6x6).
!  if1,if2     |First and last index of output array,
!  jf1,jf2     |edge centered, for active
!  kf1,kf2     |cells only (3x3x3).
!------------------------------------------------------------------------
  use amr_parameters
  use amr_commons
  use const
  implicit none
  real(dp),dimension(nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3+ncr*(1+ndim))::uin 
  real(dp),dimension(nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ncr*(ndim+1))::cru 
  real(dp),dimension(nvector,if1:if2,jf1:jf2,kf1:kf2,1:ncr*(ndim+1),1:ndim)::iFlx
  real(dp)::dx, dt
  integer ::iGrp, iP0, iP1, icrE, nGrid, ilevel !------------------------
  real(dp),save, &                                     !   Central fluxes
           dimension(nvector,iu1:iu2,ju1:ju2,ku1:ku2, ndim+1, ndim)::cFlx
  real(dp),save, &                                     ! Cell wavespeeds
           dimension(nvector,iu1:iu2,ju1:ju2,ku1:ku2,ndim)::  lmax
  ! Upwards and downwards fluxes and states of the group
  real(dp),dimension(nDim+1),save:: fdn, fup, udn, uup
  real(dp):: lminus, lplus                        ! Intercell wavespeeds
  real(dp)::dtdx, prod(ndim+1)
  integer ::i, j, k, n
  real(dp),dimension(ndim+1),save::slopeLM,slopeRM,slopeM
  real(dp),dimension(ndim+1),save::slopeLL,slopeL
  real(dp),save::vslopeLM,vslopeRM,vslopeM
  real(dp),save::vslopeLL,vslopeL,vprod
  real(dp):: vup,vdn,meanadv,meandiffv,aup,adn
  logical::interpolation=.true.,use_minmod=.false.
!------------------------------------------------------------------------
  iP0 = 1+(iGrp-1)*(ndim+1)            ! Index of CR group energy density
  iP1 = iP0+nDim     
  icrE = icrU+(ndim+1)*(iGrp-1)                    ! end index of CR group
  cru(:,:,:,:,iP0:iP1)   = uin(:,:,:,:,icrE:icrE+ndim)
  ! compute flux tensors for all the cells with correction
  call cmp_cr_flux_tensors(uin, iGrp, ngrid, cFlx, cr_vmax(ilevel))!  flux tensors

  ! Wavespeeds in each cell
  call cmp_cr_wavespeeds(uin, iGrp, ngrid, lmax, ilevel)

  ! Solve for 1D flux in X direction
  !----------------------------------------------------------------------
  dtdx=dt/dx
  do i=if1,if2                                 !
  do j=jf1,jf2                                 !        each cell in grid
  do k=kf1,kf2                                 !
     do n=1,ngrid                              ! <- buffer of grids
        fdn = cFlx(n,  i-1, j, k, :, 1    )    !
        fup = cFlx(n,  i,   j, k, :, 1    )   !  upwards and downwards
        udn = cru( n,  i-1, j, k, iP0:iP1 )   !  conditions
        uup = cru( n,  i,   j, k, iP0:iP1 )    !
        vdn  = uin( n,  i-1, j, k, 2) / uin(n,i-1,j,k,1) ! left velocity
        vup  = uin( n,  i,   j, k, 2) / uin(n,i  ,j,k,1) ! right velocity
        
        if(interpolation) then
        ! interpolation of U
        slopeLM = (fup-fdn)/dx
        slopeRM = (cFlx(n, i+1, j, k, :, 1) - fup)/dx
        prod = slopeLM*slopeRM
        slopeM=0.
        where(prod.gt.0.) slopeM=2.*prod/(slopeLM+slopeRM)
        slopeLL = (fdn - cFlx(n, i-2, j, k, :, 1))/dx
        prod = slopeLL*slopeLM
        slopeL=0.
        where(prod.gt.0) slopeL=2.*prod/(slopeLL+slopeLM)
        fdn = fdn+slopeL*0.5d0*dx
        fup = fup-slopeM*0.5d0*dx
        
        ! interpolation of F
        slopeLM = (uup-udn)/dx
        slopeRM = (cru(n, i+1, j, k, iP0:iP1) - uup)/dx
        prod = slopeLM*slopeRM
        slopeM=0.
        where(prod.gt.0) slopeM=2.*prod/(slopeLM+slopeRM)
        slopeLL = (udn - cru(n, i-2, j, k, iP0:iP1 ))/dx
        prod = slopeLL*slopeLM
        slopeL=0.
        where(prod.gt.0.) slopeL=2.*prod/(slopeLL+slopeLM)
        udn = udn+slopeL*0.5d0*dx
        uup = uup-slopeM*0.5d0*dx            
        
        ! interpolation of velocities
        vslopeLM = (vup-vdn)/dx
        vslopeRM = (uin(n,i+1,j,k,2)/uin(n,i+1,j,k,1) - vup)/dx
        vprod = vslopeLM*vslopeRM
        vslopeM=0.
        if(vprod.gt.0) vslopeM=2.*vprod/(vslopeLM+vslopeRM)
        vslopeLL = (vdn - uin(n,i-2,j,k,2)/uin(n,i-2,j,k,1))/dx
        vprod = vslopeLL*vslopeLM
        vslopeL=0.
        if(vprod.gt.0.) vslopeL=2.*vprod/(vslopeLL+vslopeLM)
        !vdn = vdn+vslopeL*0.5d0*dx
        !vup = vup-vslopeM*0.5d0*dx
        
        else if(use_minmod) then
            slopeLM = (fup-fdn)/dx
            slopeRM = (cFlx(n, i+1, j, k, :, 1) - fup)/dx
            slopeM  = minmod(slopeLM,slopeRM)
            slopeLL = (fdn - cFlx(n, i-2, j, k, :, 1))/dx
            slopeL  = minmod(slopeLL,slopeLM)
            fdn = fdn+slopeL*0.5d0*dx
            fup = fup-slopeM*0.5d0*dx

            slopeLM = (uup-udn)/dx
            slopeRM = (cru(n, i+1, j, k, iP0:iP1) - uup)/dx
            slopeM  = minmod(slopeLM,slopeRM)
            slopeLL = (udn - cru(n, i-2, j, k, iP0:iP1 ))/dx
            slopeL  = minmod(slopeLL,slopeLM)
            udn = udn+slopeL*0.5d0*dx
            uup = uup-slopeM*0.5d0*dx

            vslopeLM = (vup-vdn)/dx
            vslopeRM = (uin(n,i+1,j,k,2)/uin(n,i+1,j,k,1) - vup)/dx
            vslopeM  = vminmod(vslopeLM,vslopeRM)
            vslopeLL = (vdn - uin(n,i-2,j,k,2)/uin(n,i-2,j,k,1))/dx
            vslopeL  = vminmod(vslopeLL,vslopeLM)
            !vdn = vdn+vslopeL*0.5d0*dx
            !vup = vup-vslopeM*0.5d0*dx
        endif
        meanadv = 0.5*(vdn+vup)
        meandiffv = 0.5*( lmax(n,i-1,j,k,1) + lmax(n,i,j,k,1) )
        adn = min(meanadv-meandiffv, vdn-lmax(n,i-1,j,k,1))
        adn = max(adn,-cr_vmax(ilevel)*sqrt(1./3.))
        aup = max(meanadv+meandiffv, vup+lmax(n,i,j,k,1))
        aup = min(aup,cr_vmax(ilevel)*sqrt(1./3.))
        lminus = min(adn,0.)
        lplus = max(aup,0.)

        iFlx(n, i, j, k, iP0:iP1, 1)=&
             cmp_cr_face( fdn, fup, udn, uup, lminus, lplus)*dtdx
     end do
  end do
  end do
  end do

  ! Solve for 1D flux in Y direction
  !----------------------------------------------------------------------
#if NDIM>1
  dtdx=dt/dx
  do i=if1,if2
  do j=jf1,jf2
  do k=kf1,kf2
     do n=1,ngrid
           fdn = cFlx(n,  i, j-1, k, :, 2    )
           fup = cFlx(n,  i, j,   k, :, 2    )
           udn = cru( n,  i, j-1, k, iP0:iP1 )
           uup = cru( n,  i, j,   k, iP0:iP1 )
           vdn  = uin( n,  i, j-1, k,3) / uin(n,i,j-1,k,1) ! left velocity
           vup  = uin( n,  i ,j,   k,3) / uin(n,i,j,  k,1) ! right velocity
           
           if(interpolation) then
           ! interpolation of U
           slopeLM = (fup-fdn)/dx
           slopeRM = (cFlx(n, i, j+1, k, :, 2) - fup)/dx
           prod = slopeLM*slopeRM
           slopeM=0.
           where(prod.gt.0.) slopeM=2.*prod/(slopeLM+slopeRM)
           slopeLL = (fdn - cFlx(n, i, j-2, k, :, 2))/dx
           prod = slopeLL*slopeLM
           slopeL=0.
           where(prod.gt.0) slopeL=2.*prod/(slopeLL+slopeLM)
           fdn = fdn+slopeL*0.5d0*dx
           fup = fup-slopeM*0.5d0*dx
           
           ! interpolation of F
           slopeLM = (uup-udn)/dx
           slopeRM = (cru(n, i, j+1, k, iP0:iP1) - uup)/dx
           prod = slopeLM*slopeRM
           slopeM=0.
           where(prod.gt.0) slopeM=2.*prod/(slopeLM+slopeRM)
           slopeLL = (udn - cru(n, i, j-2, k, iP0:iP1 ))/dx
           prod = slopeLL*slopeLM
           slopeL=0.
           where(prod.gt.0.) slopeL=2.*prod/(slopeLL+slopeLM)
           udn = udn+slopeL*0.5d0*dx
           uup = uup-slopeM*0.5d0*dx            
           
           ! interpolation of velocities
           vslopeLM = (vup-vdn)/dx
           vslopeRM = (uin(n,i,j+1,k,3)/uin(n,i,j+1,k,1) - vup)/dx
           vprod = vslopeLM*vslopeRM
           vslopeM=0.
           if(vprod.gt.0) vslopeM=2.*vprod/(vslopeLM+vslopeRM)
           vslopeLL = (vdn - uin(n,i,j-2,k,3)/uin(n,i,j-2,k,1))/dx
           vprod = vslopeLL*vslopeLM
           vslopeL=0.
           if(vprod.gt.0.) vslopeL=2.*vprod/(vslopeLL+vslopeLM)
           !vdn = vdn+vslopeL*0.5d0*dx
           !vup = vup-vslopeM*0.5d0*dx
           
           else if(use_minmod) then
              slopeLM = (fup-fdn)/dx
              slopeRM = (cFlx(n, i, j+1, k, :, 2) - fup)/dx
              slopeM  = minmod(slopeLM,slopeRM)
              slopeLL = (fdn - cFlx(n, i, j-2, k, :, 2))/dx
              slopeL  = minmod(slopeLL,slopeLM)
              fdn = fdn+slopeL*0.5d0*dx
              fup = fup-slopeM*0.5d0*dx

              slopeLM = (uup-udn)/dx
              slopeRM = (cru(n, i, j+1, k, iP0:iP1) - uup)/dx
              slopeM  = minmod(slopeLM,slopeRM)
              slopeLL = (udn - cru(n, i, j-2, k, iP0:iP1))/dx
              slopeL  = minmod(slopeLL,slopeLM)
              udn = udn+slopeL*0.5d0*dx
              uup = uup-slopeM*0.5d0*dx

              vslopeLM = (vup-vdn)/dx
              vslopeRM = (uin(n,i,j+1,k,3)/uin(n,i,j+1,k,1) - vup)/dx
              vslopeM  = vminmod(vslopeLM,vslopeRM)
              vslopeLL = (vdn - uin(n,i,j-2,k,3)/uin(n,i,j-2,k,1))/dx
              vslopeL  = vminmod(vslopeLL,vslopeLM)
              !vdn = vdn+vslopeL*0.5d0*dx
              !vup = vup-vslopeM*0.5d0*dx
           endif
           meanadv = 0.5*(vdn+vup)
           meandiffv = 0.5*( lmax(n,i,j-1,k,2) + lmax(n,i,j,k,2) )
           adn = min(meanadv-meandiffv, vdn-lmax(n,i,j-1,k,2))
           adn = max(adn,-cr_vmax(ilevel)*sqrt(1./3.))
           aup = max(meanadv+meandiffv, vup+lmax(n,i,j,k,2))
           aup = min(aup,cr_vmax(ilevel)*sqrt(1./3.))
           lminus = min(adn,0.)
           lplus = max(aup,0.)
           
           iFlx(n, i, j, k, iP0:iP1, 2)=&
                cmp_cr_face( fdn, fup, udn, uup, lminus, lplus)*dtdx
     end do
  end do
  end do
  end do
#endif

  ! Solve for 1D flux in Z direction
  !----------------------------------------------------------------------
#if NDIM>2
  dtdx=dt/dx
  do i=if1,if2
  do j=jf1,jf2
  do k=kf1,kf2
     do n=1,ngrid
           fdn = cFlx(n,  i, j, k-1, :, 3    )
           fup = cFlx(n,  i, j, k,   :, 3    )
           udn = cru( n,  i, j, k-1, iP0:iP1 )
           uup = cru( n,  i, j, k,   iP0:iP1 )
           vdn  = uin( n,  i, j, k-1, 4) / uin(n,i,  j,k-1,1) ! left velocity
           vup  = uin( n,  i ,j, k,   4) / uin(n,i  ,j,k,  1) ! right velocity
           
           if(interpolation) then
           ! interpolation of U
           slopeLM = (fup-fdn)/dx
           slopeRM = (cFlx(n, i, j, k+1, :, 3) - fup)/dx
           prod = slopeLM*slopeRM
           slopeM=0.
           where(prod.gt.0.) slopeM=2.*prod/(slopeLM+slopeRM)
           slopeLL = (fdn - cFlx(n, i, j, k-2, :, 3))/dx
           prod = slopeLL*slopeLM
           slopeL=0.
           where(prod.gt.0) slopeL=2.*prod/(slopeLL+slopeLM)
           fdn = fdn+slopeL*0.5d0*dx
           fup = fup-slopeM*0.5d0*dx
           
           ! interpolation of F
           slopeLM = (uup-udn)/dx
           slopeRM = (cru(n, i, j, k+1, iP0:iP1) - uup)/dx
           prod = slopeLM*slopeRM
           slopeM=0.
           where(prod.gt.0) slopeM=2.*prod/(slopeLM+slopeRM)
           slopeLL = (udn - cru(n, i, j, k-2, iP0:iP1 ))/dx
           prod = slopeLL*slopeLM
           slopeL=0.
           where(prod.gt.0.) slopeL=2.*prod/(slopeLL+slopeLM)
           udn = udn+slopeL*0.5d0*dx
           uup = uup-slopeM*0.5d0*dx            
           
           ! interpolation of velocities
           vslopeLM = (vup-vdn)/dx
           vslopeRM = (uin(n,i,j,k+1,4)/uin(n,i,j,k+1,1) - vup)/dx
           vprod = vslopeLM*vslopeRM
           vslopeM=0.
           if(vprod.gt.0) vslopeM=2.*vprod/(vslopeLM+vslopeRM)
           vslopeLL = (vdn - uin(n,i,j,k-2,4)/uin(n,i,j,k-2,1))/dx
           vprod = vslopeLL*vslopeLM
           vslopeL=0.
           if(vprod.gt.0.) vslopeL=2.*vprod/(vslopeLL+vslopeLM)
           !vdn = vdn+vslopeL*0.5d0*dx
           !vup = vup-vslopeM*0.5d0*dx
           
           else if(use_minmod) then
              slopeLM = (fup-fdn)/dx
              slopeRM = (cFlx(n, i, j, k+1, :, 3) - fup)/dx
              slopeM  = minmod(slopeLM,slopeRM)
              slopeLL = (fdn - cFlx(n, i, j, k-2, :, 3))/dx
              slopeL  = minmod(slopeLL,slopeLM)
              fdn = fdn+slopeL*0.5d0*dx
              fup = fup-slopeM*0.5d0*dx

              slopeLM = (uup-udn)/dx
              slopeRM = (cru(n, i, j, k+1, iP0:iP1) - uup)/dx
              slopeM  = minmod(slopeLM,slopeRM)
              slopeLL = (udn - cru(n, i, j, k-2, iP0:iP1))/dx
              slopeL  = minmod(slopeLL,slopeLM)
              udn = udn+slopeL*0.5d0*dx
              uup = uup-slopeM*0.5d0*dx

              vslopeLM = (vup-vdn)/dx
              vslopeRM = (uin(n,i,j,k+1,4)/uin(n,i,j,k+1,1) - vup)/dx
              vslopeM  = vminmod(vslopeLM,vslopeRM)
              vslopeLL = (vdn - uin(n,i,j,k-2,4)/uin(n,i,j,k-2,1))/dx
              vslopeL  = vminmod(vslopeLL,vslopeLM)
              !vdn = vdn+vslopeL*0.5d0*dx
              !vup = vup-vslopeM*0.5d0*dx
           endif

           meanadv = 0.5*(vdn+vup)
           meandiffv = 0.5*( lmax(n,i,j,k-1,3) + lmax(n,i,j,k,3) )
           adn = min(meanadv-meandiffv, vdn-lmax(n,i,j,k-1,3))
           adn = max(adn,-cr_vmax(ilevel)*sqrt(1./3.))
           aup = max(meanadv+meandiffv, vup+lmax(n,i,j,k,3))
           aup = min(aup,cr_vmax(ilevel)*sqrt(1./3.))
           lminus = min(adn,0.)
           lplus  = max(aup,0.)

           iFlx(n, i, j, k, iP0:iP1, 3)=&
                cmp_cr_face( fdn, fup, udn, uup, lminus, lplus)*dtdx
      end do
  end do
  end do
  end do
#endif

end subroutine cmp_cr_faces

!************************************************************************
SUBROUTINE rotatevec(sint, cost, sinp, cosp, v1, v2, v3)
  !  Rotate vector v by t=theta and p=phi
  !  i.e. rotate to the local coordinate system from theta, phi.
  !  Hence the x-component of the result is the component of v parallel 
  !  with the theta,phi vector.
  !------------------------------------------------------------------------
    implicit none
    real(dp),intent(in):: sint, cost, sinp, cosp
    real(dp),intent(inout)::v1,v2,v3
    real(dp)::newv1, newv3
  !------------------------------------------------------------------------
    ! First apply R1, then apply R2
    newv1 =  cosp * v1 + sinp * v2
    v2 = -sinp * v1 + cosp * v2
  
    ! Now apply R2
    v1 =  sint * newv1 + cost * v3
    newv3 = -cost * newv1 + sint * v3
    v3 = newv3
END SUBROUTINE rotatevec
  
!************************************************************************
SUBROUTINE invrotatevec(sint, cost, sinp, cosp, v1, v2, v3)
  !  Inverse-rotate vector v by t=theta and p=phi,
  !  i.e. rotate v onto theta, pi
  !
  !------------------------------------------------------------------------
    implicit none
    real(dp),intent(in):: sint, cost, sinp, cosp
    real(dp),intent(inout)::v1,v2,v3
    real(dp)::newv1, newv2
  !------------------------------------------------------------------------
    ! First apply R2^-1, then apply R1^-1
    newv1 = sint * v1 - cost * v3
    v3 = cost * v1 + sint * v3
  
    ! now apply R1^-1
    v1 = cosp * newv1 - sinp * v2
    newv2 = sinp * newv1 + cosp * v2
    v2 = newv2
END SUBROUTINE invrotatevec

END MODULE cr_flux_module
