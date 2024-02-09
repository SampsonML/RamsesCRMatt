!================================================================
!================================================================
!================================================================
!================================================================
subroutine condinit(x,u,dx,nn)
  use amr_parameters
  use hydro_parameters
  implicit none
  integer ::nn                            ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:nvar+3+ncrvars)::u ! Conservative variables
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine generates initial conditions for RAMSES.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! U is the conservative variable vector. Conventions are here:
  ! U(i,1): d, U(i,2:4): d.u,d.v,d.w, U(i,5): E, U(i,6:8): Bleft, 
  ! U(i,nvar+1:nvar+3): Bright
  ! Q is the primitive variable vector. Conventions are here:
  ! Q(i,1): d, Q(i,2:4):u,v,w, Q(i,5): P, Q(i,6:8): Bleft, 
  ! Q(i,nvar+1:nvar+3): Bright, Q(i,9:8+nrad+1): Er (if FLD or electronic conduction)
  ! If nvar > 8+nrad, remaining variables (9+nrad:nvar) are treated as passive
  ! scalars in the hydro solver.
  ! U(:,:) and Q(:,:) are in user units.
  !================================================================
  integer::ivar,i
  real(dp),dimension(1:nvector,1:nvar+3+ncrvars),save::q   ! Primitive variables
  real(dp)::xl,yl,xr,yr,rl,rr,xx,yy,zz,ttmin,ttmax,pi,xcenter,Al,Ar,R0,A0,theta

  ! Call built-in initial condition generator
  call region_condinit(x,q,dx,nn)

  ! Add here, if you wish, some user-defined initial conditions
  ! ........

  ! Density
  q(1:nn,1)=1d0
  ! Velocities
  q(1:nn,2:4)=0d0
  pi=acos(-1d0)
  !tmin=11d0/12d0*pi
  !tmax=13d0/12d0*pi

  ! Left
  ttmin=-pi/12d0
  ttmax= pi/12d0
  ! Top
!!$  ttmin=5d0*pi/12d0
!!$  ttmax=7d0*pi/12d0
  xcenter=boxlen/2d0 
  R0=0.3
  A0=1d-3
  q(1:nn,5)=1d-2
  zz=0.0d0
  do i=1,nn
     xx=x(i,1)-xcenter
     yy=x(i,2)-xcenter
#if NDIM>2
     zz=x(i,3)-xcenter
#endif

     !if(xx.gt.-0.3d0*boxlen .and. xx.lt.-0.2d0*boxlen .and. yy.gt.-0.025*boxlen .and. yy.lt.0.025*boxlen .and. zz.gt.-0.1 .and. zz.lt.0.1)then

     ! Energy patch for loop around z-axis
     rr=sqrt(xx**2+yy**2)
     ! Left
     theta=atan(yy/xx)
     q(i,5)=1d-2
     if(rr.gt.0.25d0*boxlen .and. rr.lt.0.35d0*boxlen .and. theta.gt.ttmin .and. theta.lt.ttmax .and. xx.gt.0d0 .and. zz.gt.-0.1*boxlen .and. zz.lt.0.1*boxlen)then
     ! Top
!!$     theta=atan(abs(yy/xx))
!!$     if(rr.gt.0.2d0*boxlen .and. rr.lt.0.3d0*boxlen .and. theta.gt.ttmin .and. theta.lt.ttmax .and. yy.gt.0d0 .and. zz.gt.-0.1*boxlen .and. zz.lt.0.1*boxlen)then
!!$
!!$     ! Energy patch for loop around x-axis
!!$     rr=sqrt(zz**2+yy**2)
!!$     theta=atan(yy/zz)
!!$     if(rr.gt.0.2d0*boxlen .and. rr.lt.0.3d0*boxlen .and. theta.gt.ttmin .and. theta.lt.ttmax .and. zz.lt.0d0 .and. xx.gt.-0.1*boxlen .and. xx.lt.0.1*boxlen)then
           !q(i,8+ivar)=1.2d1
           q(i,icrU)=1.2d1
     else
           !q(i,8+ivar)=1.0d1
           q(i,icrU)=1.0d1
     endif

!!$     xl=x(i,1)-0.5*dx-xcenter
!!$     xr=x(i,1)+0.5*dx-xcenter
!!$     yl=x(i,2)-0.5*dx-xcenter
!!$     yr=x(i,2)+0.5*dx-xcenter
!!$     rl=sqrt(xl*xl+yl*yl)
!!$     rr=sqrt(xr*xr+yr*yr)

     ! loop around x-axis
!!$     xl=x(i,3)-0.5*dx-boxlen/2.0
!!$     xr=x(i,3)+0.5*dx-boxlen/2.0
!!$     yl=x(i,2)-0.5*dx-boxlen/2.0
!!$     yr=x(i,2)+0.5*dx-boxlen/2.0
!!$
!!$     Ar = A0*max(R0-sqrt(xl**2+yr**2),-boxlen)
!!$     Al = A0*max(R0-sqrt(xl**2+yl**2),-boxlen)
!!$     q(i,8)=(Ar-Al)/dx
!!$     Ar = A0*max(R0-sqrt(xr**2+yr**2),-boxlen)
!!$     Al = A0*max(R0-sqrt(xr**2+yl**2),-boxlen)
!!$     q(i,nvar+3)=(Ar-Al)/dx
!!$     Ar = A0*max(R0-sqrt(xr**2+yl**2),-boxlen)
!!$     Al = A0*max(R0-sqrt(xl**2+yl**2),-boxlen)
!!$     q(i,7)=(Al-Ar)/dx
!!$     Ar = A0*max(R0-sqrt(xr**2+yr**2),-boxlen)
!!$     Al = A0*max(R0-sqrt(xl**2+yr**2),-boxlen)
!!$     q(i,nvar+2)=(Al-Ar)/dx
!!$     q(i,6)=0.0
!!$     q(i,nvar+1)=0.0

     ! loop around z-axis
     xl=x(i,1)-0.5*dx-boxlen/2.0
     xr=x(i,1)+0.5*dx-boxlen/2.0
     yl=x(i,2)-0.5*dx-boxlen/2.0
     yr=x(i,2)+0.5*dx-boxlen/2.0

     Ar = A0*max(R0-sqrt(xl**2+yr**2),-boxlen)
     Al = A0*max(R0-sqrt(xl**2+yl**2),-boxlen)
     q(i,6)=(Ar-Al)/dx
     Ar = A0*max(R0-sqrt(xr**2+yr**2),-boxlen)
     Al = A0*max(R0-sqrt(xr**2+yl**2),-boxlen)
     q(i,nvar+1)=(Ar-Al)/dx
     Ar = A0*max(R0-sqrt(xr**2+yl**2),-boxlen)
     Al = A0*max(R0-sqrt(xl**2+yl**2),-boxlen)
     q(i,7)=(Al-Ar)/dx
     Ar = A0*max(R0-sqrt(xr**2+yr**2),-boxlen)
     Al = A0*max(R0-sqrt(xl**2+yr**2),-boxlen)
     q(i,nvar+2)=(Al-Ar)/dx
     q(i,8)=0.0
     q(i,nvar+3)=0.0

!!$     q(i,6)=1.
!!$     q(i,nvar+1)=1.
!!$     q(i,7)=0.
!!$     q(i,nvar+2)=0.

!!$     q(i,6     )= yy/sqrt(xl*xl+yy*yy)
!!$     q(i,7     )=-xx/sqrt(xx*xx+yl*yl)
!!$     q(i,nvar+1)= yy/sqrt(xr*xr+yy*yy)
!!$     q(i,nvar+2)=-xx/sqrt(xx*xx+yr*yr)
!!$     write(*,'(2f10.5,5ES13.5)')xx,yy,q(i,6),q(i,7),q(i,nvar+1),q(i,nvar+2),(q(i,nvar+2)-q(i,7))+(q(i,nvar+1)-q(i,6))

!!$     if (rl==0.)then
!!$        q(i,6     )=0d0
!!$        q(i,7     )=0d0
!!$     else
!!$        q(i,6     )= yy/sqrt(xl*xl+yy*yy)
!!$        q(i,7     )=-xx/sqrt(xx*xx+yl*yl)
!!$     endif
!!$     if (rr==0.)then
!!$        q(i,nvar+1)=0d0
!!$        q(i,nvar+2)=0d0
!!$     else
!!$        q(i,nvar+1)= yy/sqrt(xr*xr+yy*yy)
!!$        q(i,nvar+2)=-xx/sqrt(xx*xx+yr*yr)
!!$     endif
  enddo

!!$  q(1:nn,6     )=1d0
!!$  q(1:nn,nvar+1)=1d0
!!$  q(1:nn,7     )=0d0
!!$  q(1:nn,nvar+2)=0d0
!!$  q(1:nn,8     )=1d0
!!$  q(1:nn,nvar+3)=1d0
     
  ! Convert primitive to conservative variables
  ! density -> density
  u(1:nn,1)=q(1:nn,1)
  ! velocity -> momentum
  u(1:nn,2)=q(1:nn,1)*q(1:nn,2)
  u(1:nn,3)=q(1:nn,1)*q(1:nn,3)
  u(1:nn,4)=q(1:nn,1)*q(1:nn,4)
  ! kinetic energy
  u(1:nn,5)=0.0d0
  u(1:nn,5)=u(1:nn,5)+0.5*q(1:nn,1)*q(1:nn,2)**2
  u(1:nn,5)=u(1:nn,5)+0.5*q(1:nn,1)*q(1:nn,3)**2
  u(1:nn,5)=u(1:nn,5)+0.5*q(1:nn,1)*q(1:nn,4)**2
  ! pressure -> total fluid energy
  u(1:nn,5)=u(1:nn,5)+q(1:nn,5)/(gamma-1.0d0)
  ! magnetic energy -> total fluid energy
  u(1:nn,5)=u(1:nn,5)+0.125d0*(q(1:nn,6)+q(1:nn,nvar+1))**2
  u(1:nn,5)=u(1:nn,5)+0.125d0*(q(1:nn,7)+q(1:nn,nvar+2))**2
  u(1:nn,5)=u(1:nn,5)+0.125d0*(q(1:nn,8)+q(1:nn,nvar+3))**2
  u(1:nn,6:8)=q(1:nn,6:8)
  u(1:nn,nvar+1:nvar+3)=q(1:nn,nvar+1:nvar+3)
  !u(1:nn,iCRflx:iCRflx-1+ncrflx)=q(1:nn,iCRflx:iCRflx-1+ncrflx)
  u(1:nn,icrU:icrU+ncrvars-1)=q(1:nn,icrU:icrU+ncrvars-1)

#if NENER>0
  ! radiative pressure -> radiative energy
  ! radiative energy -> total fluid energy
  do ivar=1,nener
     u(1:nn,8+ivar)=q(1:nn,8+ivar)/(gamma_rad(ivar)-1.0d0)
     u(1:nn,5)=u(1:nn,5)+q(1:nn,8+ivar)/(gamma_rad(ivar)-1.0d0)
!!$     u(1:nn,ndim+2)=u(1:nn,ndim+2)+u(1:nn,ndim+2+ivar) !WARNING HERE!!!! DO  NOT WORK FOR NON-THERMAL ENERGIES
  enddo
#endif
#if NVAR>8+NENER
  ! passive scalars
  do ivar=9+nener,nvar
     u(1:nn,ivar)=q(1:nn,1)*q(1:nn,ivar)
  end do
#endif
  
end subroutine condinit
!================================================================
!================================================================
!================================================================
!================================================================
subroutine velana(x,v,dx,t,ncell)
  use amr_parameters
  use hydro_parameters  
  implicit none
  integer ::ncell                         ! Size of input arrays
  real(dp)::dx                            ! Cell size
  real(dp)::t                             ! Current time
  real(dp),dimension(1:nvector,1:3)::v    ! Velocity field
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine computes the user defined velocity fields.
  ! x(i,1:ndim) are cell center position in [0,boxlen] (user units).
  ! v(i,1:3) is the imposed 3-velocity in user units.
  !================================================================
  integer::i
  real(dp)::xx,yy,zz,vx,vy,vz,rr,tt,omega,aa,twopi

  ! Add here, if you wish, some user-defined initial conditions
  aa=1.0
  twopi=2d0*ACOS(-1d0)
  do i=1,ncell

!!$     xx=x(i,1)
!!$     yy=x(i,2)
!!$     zz=x(i,3)
!!$
!!$     ! ABC
!!$     vx=aa*(cos(twopi*yy)+sin(twopi*zz))
!!$     vy=aa*(sin(twopi*xx)+cos(twopi*zz))
!!$     vz=aa*(cos(twopi*xx)+sin(twopi*yy))

!!$     ! 1D advection test
!!$     vx=1.0_dp
!!$     vy=0.0_dp
!!$     vz=0.0_dp

!!$     ! Ponomarenko
!!$     xx=xx-boxlen/2.0
!!$     yy=yy-boxlen/2.0
!!$     rr=sqrt(xx**2+yy**2)
!!$     if(yy>0)then
!!$        tt=acos(xx/rr)
!!$     else
!!$        tt=-acos(xx/rr)+twopi
!!$     endif
!!$     if(rr<1.0)then
!!$        omega=0.609711
!!$        vz=0.792624
!!$     else
!!$        omega=0.0
!!$        vz=0.0
!!$     endif
!!$     vx=-sin(tt)*rr*omega
!!$     vy=+cos(tt)*rr*omega
     
!!$     v(i,1)=vx
!!$     v(i,2)=vy
!!$     v(i,3)=vz

  end do


end subroutine velana
