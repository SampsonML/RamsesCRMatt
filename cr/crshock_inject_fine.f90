!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine crshock_inject_fine(ilevel)
  use amr_commons
  use hydro_commons
#ifndef WITHOUTMPI
  use mpi_mod
#endif
  implicit none
  integer::ilevel
  !--------------------------------------------------------------------------
  ! 
  ! 
  ! 
  ! 
  !--------------------------------------------------------------------------
  integer,dimension(1:nvector),save::ind_grid,ind_cell,ind_shock2
  integer,dimension(1:nvector),save::grid_index
  integer,dimension(1:nvector),save::cell_lev
  integer,dimension(1:nvector)::cell_index,cell_index_father
  integer,dimension(1:ncpu   )::nshock_cpu,nshock_cpu_all,nshock_tot

  logical,dimension(1:ndim)::period

  real(dp),dimension(1:nvector,1:ndim),save::xtest

  integer,allocatable,dimension(:)::ind_shock,ind_shock3,ind_shock3_all
  integer,allocatable,dimension(:)::cpu_map_loc,cpu_map_loc_all
  integer,allocatable,dimension(:)::idist_inject,idist_inject_all,idist_inject_tosend

  real(dp),allocatable,dimension(:,:)::x_shock,x_shock_big,x_shock_all,x_big_all
  real(dp),allocatable,dimension(:,:)::shock_dir,shock_dir_big,shock_dir_all
  real(dp),allocatable,dimension(:)::divv_shock,divv_shock_all
  real(dp),allocatable,dimension(:)::normdivv2,normdivv2_all
  real(dp),allocatable,dimension(:)::etherm_shock,etherm_shock_all
  real(dp),allocatable,dimension(:)::eint_shock,eint_shock_all
  real(dp),allocatable,dimension(:)::ecrdiss,ecrdiss_all
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:ndim)::delta

  integer::info
  integer::i,ii,ind,iskip,igrid,ncache,ngrid,ngrid_tot,imin,imax,icpu,iener
  integer::ishock,ishock2,ishock_myid,nshock,nshock_all,nshock_myid,nx_loc,idist,mydist
  integer::icell,isweep,nsweep,ncell,ind_cr,idim

  real(dp)::dx,dx_loc,scale,d,oneoverd,u,v,w,A,B,C,e,eo,Ptot,etot,mach_loc
  real(dp)::oneoverdx,gg,g1,g2,deltae,decr,theta,ctheta
  real(dp)::pi,theta_crit,oneoverdtheta,xi,xi0,divv,weight

  if(verbose)write(*,111)ilevel

  period(1)=(nx==1)
#if NDIM>1
  if(ndim>1)period(2)=(ny==1)
#endif
#if NDIM>2
  if(ndim>2)period(3)=(nz==1)  
#endif

  ind_cr=8
  if(twotemp)ind_cr=9
  if(mom_cr_advect)ind_cr=icru-1
  
  pi=acos(-1d0)
  theta_crit=0.25d0*pi
  oneoverdtheta=18d0/pi
  xi0=0.5d0*eta_mach

  ! Mesh spacing in that level
  dx=0.5D0**ilevel 
  oneoverdx=1d0/dx
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale

  ! ========================================================
  ! Do the shock injection in the post-shock region
  ! ========================================================
  if(t.gt.tmin_acc)then

     ! ========================================================
     ! Count the number of shock cells 
     ! ========================================================
     ncache=active(ilevel)%ngrid
     nshock=0
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
        end do
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              icell=iskip+ind_grid(i)
              mach_loc=uold(icell,imach)/uold(icell,1)
              if(mach_loc.gt.0.0d0.and.ediss_cr(icell).gt.0.0d0)then
                 nshock=nshock+1
              endif
           enddo
        enddo
     end do
     nshock_myid=nshock
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(nshock,nshock_all,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
     nshock=nshock_all
#endif
     if(myid==1)write(*,*)nshock,' shock cells in level ',ilevel

     ! Store all the upstream/downstream positions in one big array 
     nshock_cpu=0
     nshock_cpu(myid)=nshock_myid
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(nshock_cpu,nshock_cpu_all,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
     nshock_tot(1)=nshock_cpu_all(1)
#endif
     do icpu=2,ncpu
        nshock_tot(icpu)=nshock_tot(icpu-1)+nshock_cpu_all(icpu)
     end do
     nshock_cpu=nshock_cpu_all
     ! Store all the upstream/downstream positions in one big array 
     nshock_cpu=0
     nshock_cpu(myid)=nshock_myid
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(nshock_cpu,nshock_cpu_all,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
     nshock_tot(1)=nshock_cpu_all(1)
#endif
     do icpu=2,ncpu
        nshock_tot(icpu)=nshock_tot(icpu-1)+nshock_cpu_all(icpu)
     end do
     nshock_cpu=nshock_cpu_all

     if(myid==1)then
        imin=1
     else
        imin=nshock_tot(myid-1)+1
     endif
     imax=nshock_tot(myid)
     if(nshock_cpu(myid).ne.(imax-imin+1).and.nshock_cpu(myid).ne.nshock_myid)then
        write(*,*)'something bad in computing number of shock cells'
#ifndef WITHOUTMPI
        call MPI_ABORT(MPI_COMM_WORLD,1,info)
#else
        stop
#endif
     endif

     allocate(x_shock  (1:nshock_myid,1:ndim))
     allocate(shock_dir(1:nshock_myid,1:ndim)) ! This time it runs over ndim (since it is not used to project B)
     allocate(x_shock_big  (1:nshock,1:ndim))
     allocate(shock_dir_big(1:nshock,1:ndim))
     allocate(x_shock_all    (1:nshock,1:ndim))
     allocate(shock_dir_all(1:nshock,1:ndim))
     allocate(divv_shock    (1:nshock))
     allocate(divv_shock_all(1:nshock))
     allocate(normdivv2    (1:nshock))
     allocate(normdivv2_all(1:nshock))
     allocate(etherm_shock    (1:nshock))
     allocate(etherm_shock_all(1:nshock))
     allocate(eint_shock    (1:nshock))
     allocate(eint_shock_all(1:nshock))
     allocate(ecrdiss    (1:nshock))
     allocate(ecrdiss_all(1:nshock))
     allocate(idist_inject    (1:nshock))
     allocate(idist_inject_all(1:nshock))
     x_shock  =0.0d0;x_shock_big  =0.0d0;x_shock_all  =0.0d0
     shock_dir=0.0d0;shock_dir_big=0.0d0;shock_dir_all=0.0d0
     divv_shock=0.0d0;divv_shock_all=0.0d0
     normdivv2=0.0d0;normdivv2_all=0.0d0
     etherm_shock=0.0d0;etherm_shock_all=0.0d0
     eint_shock=0.0d0;eint_shock_all=0.0d0
     ecrdiss   =0.0d0;ecrdiss_all   =0.0d0


     ! ========================================================
     ! Get the shock unit vector for shock cells only.
     ! Get the cell position for shock cells only.
     ! Compute the divu in the shock cell only
     ! ========================================================
     ! This part is to scan only octs having active shock cells
     ngrid=0;ishock=0;ishock2=0;ngrid_tot=0
     do igrid=1,active(ilevel)%ngrid
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           icell=iskip+active(ilevel)%igrid(igrid)
           mach_loc=uold(icell,imach)/uold(icell,1)
           if(mach_loc.gt.0.0d0.and.ediss_cr(icell).gt.0.0d0)then
              ngrid=ngrid+1
              ngrid_tot=ngrid_tot+1
              ind_grid  (ngrid)=active(ilevel)%igrid(igrid)
              ind_shock2(ngrid)=icell
              ecrdiss(imin+ishock2)=ediss_cr(icell)

              call cmp_divvold(divv,icell,ilevel)
              divv_shock(imin+ishock2)=abs(divv)

              call cmp_ethermold(eo,icell)
              etherm_shock(imin+ishock2)=eo
              call cmp_eintold(eo,icell)
              eint_shock(imin+ishock2)=eo

              ishock2=ishock2+1

           endif
           if(ngrid.eq.nvector)then
              call get_shock_unit_vector2(ind_grid,ngrid,ind_shock2,nshock_myid,ishock,x_shock,shock_dir,ilevel)
              ngrid=0
           endif
        end do
     end do

     if(ngrid.gt.0)call get_shock_unit_vector2(ind_grid,ngrid,ind_shock2,nshock_myid,ishock,x_shock,shock_dir,ilevel)

     if(nshock_cpu(myid).gt.0)then
        x_shock_big  (imin:imax,1:ndim)=x_shock  (1:nshock_myid,1:ndim)
        shock_dir_big(imin:imax,1:ndim)=shock_dir(1:nshock_myid,1:ndim)
     endif
     
     ! Share the shock positions, shock normal and energy to dissipate amongst all CPUs
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(x_shock_big  ,x_shock_all     ,ndim*nshock,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     x_shock_big  =x_shock_all
     call MPI_ALLREDUCE(shock_dir_big,shock_dir_all   ,ndim*nshock,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     shock_dir_big=shock_dir_all
     call MPI_ALLREDUCE(divv_shock   ,divv_shock_all  ,nshock     ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     divv_shock=divv_shock_all
     call MPI_ALLREDUCE(etherm_shock ,etherm_shock_all,nshock     ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     etherm_shock=etherm_shock_all
     call MPI_ALLREDUCE(eint_shock   ,eint_shock_all ,nshock      ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     eint_shock=eint_shock_all
     call MPI_ALLREDUCE(ecrdiss      ,ecrdiss_all    ,nshock      ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     ecrdiss=ecrdiss_all
#endif
     ! ========================================================

  
     ! ========================================================
     ! Compute divu along the post-shock direction up to
     ! naway_inject cells away from the shock and 
     ! the cell with the lowest value of |divu|
     ! ========================================================
     idist_inject=0 ! by defaul inject in the shock cell
     do idist=1,naway_inject*10
        do ishock=1,nshock,nvector
           nsweep=MIN(nvector,nshock-ishock+1)

           ! Now we search for the post-shock cells in all CPUs
           ! The cpu_map(icell)==myid condition tests if the cell 
           ! actually belongs to the current CPU
           xtest=0.0d0
           do isweep=1,nsweep
              ii=isweep+ishock-1
              delta(1:ndim)=shock_dir_big(ii,1:ndim)*(0.55d0+dble(idist-1)/10d0)*dx ! 1/10*dx steps, starting at 0.55dx
              xtest(isweep,1:ndim)=(x_shock_big(ii,1:ndim)+delta(1:ndim))*scale
              do idim=1,ndim
                 if(period(idim).and.xtest(isweep,idim)>boxlen)xtest(isweep,idim)=xtest(isweep,idim)-boxlen
                 if(period(idim).and.xtest(isweep,idim)<0.    )xtest(isweep,idim)=xtest(isweep,idim)+boxlen
              enddo
           enddo
           ! grid_index might end up not being necessary: remove everywhere?
           call get_mycell_index(cell_index,cell_index_father,grid_index,cell_lev,xtest,nlevelmax,nsweep)

           do isweep=1,nsweep
              if(cpu_map(cell_index_father(isweep))==myid)then              
                 icell=cell_index(isweep)
                 call cmp_divvold(divv,icell,ilevel)
                 divv=abs(divv)
                 ii=isweep+ishock-1
                 call cmp_ethermold(eo,icell)
                 if(divv.lt.divv_shock(ii))then
                    idist_inject(ii) = idist
                    divv_shock  (ii) = divv
                 endif
                 if(idist.ge.1.and.idist.le.3)then
                    weight=1d0/MAX(divv*divv,1d-10)
                 endif
                 call cmp_eintold(eo,icell)

              endif
           enddo

        enddo

     enddo

     ! Share the min(|divu|) and the corresponding max distance amongst all CPUs
#ifndef WITHOUTMPI
     ! divv is not supposed to be zero, since all CPUs have at least the value of divv for the shock cell (strongest compression)
     call MPI_ALLREDUCE(divv_shock  ,divv_shock_all  ,nshock,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,info)
     allocate(idist_inject_tosend(1:nshock))
     idist_inject_tosend=0
     do ishock=1,nshock
        if(divv_shock(ishock)==divv_shock_all(ishock))then
           idist_inject_tosend(ishock)=idist_inject(ishock)
        endif
     enddo
     divv_shock=divv_shock_all
     ! everyone's value is zero except for the CPU that has min_CPU(divv)=min_all(divv)
     call MPI_ALLREDUCE(idist_inject_tosend,idist_inject_all,nshock,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
     idist_inject=idist_inject_all
     deallocate(idist_inject_tosend)

     call MPI_ALLREDUCE(etherm_shock,etherm_shock_all,nshock,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,info)
     etherm_shock=etherm_shock_all
     call MPI_ALLREDUCE(eint_shock  ,eint_shock_all  ,nshock,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,info)
     eint_shock=eint_shock_all
#endif
     ! ========================================================


     ! ========================================================
     ! Repeat to get the 1/divv**2 normalisation
     ! ========================================================
!!$     do idist=1,4
     do idist=1,twotondim
        do ishock=1,nshock,nvector
           nsweep=MIN(nvector,nshock-ishock+1)

           ! Now we search for the post-shock cells in all CPUs
           ! The cpu_map(icell)==myid condition tests if the cell 
           ! actually belongs to the current CPU
           xtest=0.0d0
           do isweep=1,nsweep
              ii=isweep+ishock-1
              mydist=idist_inject(ii)

              delta(1:ndim)=shock_dir_big(ii,1:ndim)*(0.55d0+dble(mydist-1)/10d0)*dx

#if NDIM==1
              ! 1D poor-man CIC
              if(idist==1)then
                 delta(1)=delta(1)-0.50d0*dx
              endif
              if(idist==2)then
                 delta(1)=delta(1)+0.50d0*dx
              endif
#endif
#if NDIM==2
              ! 2D poor-man CIC
              if(idist==1)then
                 delta(1)=delta(1)-0.50d0*dx
                 delta(2)=delta(2)-0.50d0*dx
              endif
              if(idist==2)then
                 delta(1)=delta(1)+0.50d0*dx
                 delta(2)=delta(2)-0.50d0*dx
              endif
              if(idist==3)then
                 delta(1)=delta(1)-0.50d0*dx
                 delta(2)=delta(2)+0.50d0*dx
              endif
              if(idist==4)then
                 delta(1)=delta(1)+0.50d0*dx
                 delta(2)=delta(2)+0.50d0*dx
              endif
#endif
#if NDIM==3
              ! 3D poor-man CIC
              if(idist==1)then
                 delta(1)=delta(1)-0.50d0*dx
                 delta(2)=delta(2)-0.50d0*dx
                 delta(3)=delta(3)-0.50d0*dx
              endif
              if(idist==2)then
                 delta(1)=delta(1)+0.50d0*dx
                 delta(2)=delta(2)-0.50d0*dx
                 delta(3)=delta(3)-0.50d0*dx
              endif
              if(idist==3)then
                 delta(1)=delta(1)-0.50d0*dx
                 delta(2)=delta(2)+0.50d0*dx
                 delta(3)=delta(3)-0.50d0*dx
              endif
              if(idist==4)then
                 delta(1)=delta(1)+0.50d0*dx
                 delta(2)=delta(2)+0.50d0*dx
                 delta(3)=delta(3)-0.50d0*dx
              endif
              if(idist==5)then
                 delta(1)=delta(1)-0.50d0*dx
                 delta(2)=delta(2)-0.50d0*dx
                 delta(3)=delta(3)+0.50d0*dx
              endif
              if(idist==6)then
                 delta(1)=delta(1)+0.50d0*dx
                 delta(2)=delta(2)-0.50d0*dx
                 delta(3)=delta(3)+0.50d0*dx
              endif
              if(idist==7)then
                 delta(1)=delta(1)-0.50d0*dx
                 delta(2)=delta(2)+0.50d0*dx
                 delta(3)=delta(3)+0.50d0*dx
              endif
              if(idist==8)then
                 delta(1)=delta(1)+0.50d0*dx
                 delta(2)=delta(2)+0.50d0*dx
                 delta(3)=delta(3)+0.50d0*dx
              endif
#endif

              xtest(isweep,1:ndim)=(x_shock_big(ii,1:ndim)+delta(1:ndim))*scale
              do idim=1,ndim
                 if(period(idim).and.xtest(isweep,idim)>boxlen)xtest(isweep,idim)=xtest(isweep,idim)-boxlen
                 if(period(idim).and.xtest(isweep,idim)<0.    )xtest(isweep,idim)=xtest(isweep,idim)+boxlen
              enddo
           enddo
           ! grid_index might end up not being necessary: remove everywhere?
           call get_mycell_index(cell_index,cell_index_father,grid_index,cell_lev,xtest,nlevelmax,nsweep)

           do isweep=1,nsweep
              if(cpu_map(cell_index_father(isweep))==myid)then              
                 icell=cell_index(isweep)
                 call cmp_divvold(divv,icell,ilevel)
                 divv=abs(divv)
                 ii=isweep+ishock-1
                 weight=1d0/MAX(abs(divv),1d-10)
                 normdivv2(ii)=normdivv2(ii)+weight
              endif
           enddo
        enddo
     enddo
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(normdivv2   ,normdivv2_all   ,nshock,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     normdivv2=normdivv2_all
#endif



     ! ========================================================
     ! Now do the shock injection in the cell an the -1,+1 cell
     ! identified as the max of the energy (post-shock)
     ! ========================================================
     do idist=1,1
        do ishock=1,nshock,nvector
           nsweep=MIN(nvector,nshock-ishock+1)

           ! Now we search for the post-shock cells in all CPUs
           ! The cpu_map(icell)==myid condition tests if the cell 
           ! actually belongs to the current CPU
           do isweep=1,nsweep
              ii=isweep+ishock-1
              mydist=idist_inject(ii)

              delta(1:ndim)=shock_dir_big(ii,1:ndim)*(0.55d0+dble(mydist-1)/10d0)*dx
              xtest(isweep,1:ndim)=(x_shock_big(ii,1:ndim)+delta(1:ndim))*scale
              do idim=1,ndim
                 if(period(idim).and.xtest(isweep,idim)>boxlen)xtest(isweep,idim)=xtest(isweep,idim)-boxlen
                 if(period(idim).and.xtest(isweep,idim)<0.    )xtest(isweep,idim)=xtest(isweep,idim)+boxlen
              enddo
           enddo
           call get_mycell_index(cell_index,cell_index_father,grid_index,cell_lev,xtest,nlevelmax,nsweep)
           do isweep=1,nsweep
              icell=cell_index(isweep)
              if(cpu_map(cell_index_father(isweep))==myid)then
                 
                 ii=isweep+ishock-1
                 call cmp_ethermold(e,icell)
                 call cmp_divvold(divv,icell,ilevel)

                 weight=1d0/MAX(abs(divv),1d-10)
                 decr=ecrdiss(ii)
                 decr=decr* 2d0**(ndim*(ilevel-cell_lev(isweep)))
                 deltae=MIN(decr,0.5d0*e)
                 if(decr.gt.0.5d0*e)then
                    if(cell_lev(isweep).lt.ilevel)then
                       if(e.lt.0.0d0)write(*,'(A,2e12.4,4I9)')'negative e',unew(icell,ind_cr+1),e,icell,cell_lev(isweep),ilevel,myid
                       unew(icell,ind_cr+1)=unew(icell,ind_cr+1)+0.5d0*e
                    else
                       if(e.lt.0.0d0)write(*,'(A,2e12.4,4I9)')'negative e',uold(icell,ind_cr+1),e,icell,cell_lev(isweep),ilevel,myid
                       uold(icell,ind_cr+1)=uold(icell,ind_cr+1)+0.5d0*e
                    endif
                 else
                    if(cell_lev(isweep).lt.ilevel)then
                       if(deltae.lt.0.0d0)write(*,'(A,2e12.4,4I9)')'negative deltae',unew(icell,ind_cr+1),deltae,icell,cell_lev(isweep),ilevel,myid
                       unew(icell,ind_cr+1)=unew(icell,ind_cr+1)+deltae
                    else
                       if(deltae.lt.0.0d0)write(*,'(A,2e12.4,4I9)')'negative deltae',uold(icell,ind_cr+1),deltae,icell,cell_lev(isweep),ilevel,myid
                       uold(icell,ind_cr+1)=uold(icell,ind_cr+1)+deltae
                    endif
                 endif
                 
              endif
           enddo

        enddo
     enddo

     deallocate(x_shock  ,x_shock_big  ,x_shock_all  )
     deallocate(shock_dir,shock_dir_big,shock_dir_all)
     deallocate(divv_shock,divv_shock_all)
     deallocate(normdivv2,normdivv2_all)
     deallocate(etherm_shock,etherm_shock_all)
     deallocate(eint_shock,eint_shock_all)
     deallocate(ecrdiss   ,ecrdiss_all   )
     deallocate(idist_inject,idist_inject_all)

  endif


111 format('   Entering shock_fine for level ',i2)

end subroutine crshock_inject_fine
