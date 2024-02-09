!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine shock_fine(ilevel)
  use amr_commons
  use hydro_commons
  use mpi_mod
  implicit none
  integer::ilevel
  !--------------------------------------------------------------------------
  ! 
  ! 
  ! 
  ! 
  !--------------------------------------------------------------------------
  integer,dimension(1:nvector),save::ind_grid,ind_cell,ind_shock2
  integer,dimension(1:nvector),save::cell_index,grid_index
  integer,dimension(1:nvector),save::cell_lev
  integer,dimension(1:nvector)::cell_index_pre,cell_index_pos,cell_index_father
  integer,dimension(1:ncpu   )     ::nshock_cpu,nshock_cpu_all,nshock_tot

  logical,dimension(1:ndim)::period

  real(dp),dimension(1:nvector,1:ndim),save::xtest
  real(dp),dimension(1:nvector,1:ndim),save::xgrid
  real(dp),dimension(1:nvector)::dens_pre,pres_pre,ethe_pre,geff_pre,gamm_pre,csou_pre
  real(dp),dimension(1:nvector)::dens_pos,pres_pos,ethe_pos,geff_pos,gamm_pos
  real(dp),dimension(1:nvector)::dens_pre_all,pres_pre_all,ethe_pre_all,geff_pre_all,gamm_pre_all,csou_pre_all
  real(dp),dimension(1:nvector)::dens_pos_all,pres_pos_all,ethe_pos_all,geff_pos_all,gamm_pos_all
  real(dp),dimension(1:nvector,1:3)::Bmag_pre,Bmag_pre_all
  real(dp),dimension(1:nvector,1:ncr)::ecr_pre,ecr_pre_all

  integer,allocatable,dimension(:)::ind_shock,ind_shock3,ind_shock3_all
  integer,allocatable,dimension(:)::cpu_map_loc,cpu_map_loc_all
  integer,allocatable,dimension(:)::idist_inject,idist_inject_all

  real(dp),allocatable,dimension(:,:)::x_pos,x_pre,x_pre_big,x_pos_big
  real(dp),allocatable,dimension(:,:)::x_shock,x_shock_big,x_shock_all,x_big_all
  real(dp),allocatable,dimension(:,:)::shock_dir,shock_dir_big,shock_dir_all
  real(dp),allocatable,dimension(:,:)::Bpre,ecpre
  real(dp),allocatable,dimension(:)::dpre,dpos,ppre,ppos,epre,epos,cpre,depos
  real(dp),allocatable,dimension(:)::gpre,gpos,gapre,gapos
  real(dp),allocatable,dimension(:)::divv_shock,divv_shock_all
  real(dp),allocatable,dimension(:)::etherm_shock,etherm_shock_all
  real(dp),allocatable,dimension(:)::eint_shock,eint_shock_all
  real(dp),allocatable,dimension(:)::ecrdiss,ecrdiss_all
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:3)::Bpre_proj
  real(dp),dimension(1:ndim)::delta,xx
  real(dp),dimension(1:twotondim,1:3)::xc

  integer::info
  integer::i,ii,ind,iskip,igrid,ncache,ngrid,ngrid_tot,imin,imax,icpu,iener
  integer::ishock,ishock2,ishock_myid,nshock,nshock_all,nshock_myid,nx_loc,idist
  integer::icell,isweep,nsweep,ncell,ind_cr,ix,iy,iz,idim

  real(dp)::dx,dx_loc,scale,d,oneoverd,u,v,w,A,B,C,e,eo,Ptot,etot,pratio,mach,mach2,mach_loc
  real(dp)::oneoverdx,gg,g1,g2,CC,deltae,decr,Rcomp,ediss_th,Bnorm2,Bnorm_pre,theta,ctheta
  real(dp)::pi,theta_crit,oneoverdtheta,xi,xi0,divv
  real(dp)::pcr_tot,fcr,creff

  logical::okoct

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

  ! =============== Step 1 =================================
  ! Detect cells with local minimum of divu and get their
  ! Mach number using at most 2 cells distant from the shock
  ! ========================================================
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     call get_shock_cells(ind_grid,ngrid,ilevel)
  end do
  ! =============== Step 1 =================================


#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(decr_lost_mpi,decr_lost_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  decr_lost=decr_lost+decr_lost_all
  decr_lost_mpi=0.0d0
  call MPI_ALLREDUCE(decr_gain_mpi,decr_gain_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  decr_gain=decr_gain+decr_gain_all
  decr_gain_mpi=0.0d0
#else
  decr_lost=decr_lost_mpi
  decr_gain=decr_gain_mpi
#endif

  if(ndist>2)then ! Search for cells further than one oct (2 cells) away from shock cells
     ! =============== Step 2 =================================
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
              if(uold(icell,imach)/uold(icell,1).gt.mach_strong)then
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
     ! =============== Step 2 =================================

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
     ! Now we have:
     ! nshock_cpu stores the number of shock cells in each CPU (know by all CPUs)
     ! nshock_tot stores the number of shock cells from CPU 1 to CPU n (know by all CPUs)

     allocate(shock_dir(1:nshock_myid,1:3)) ! This has to run over the 3 possible dimensions
     allocate(Bpre(1:nshock_myid,1:3))      ! since B field is 3D even for ndim<3 runs
     allocate(x_shock  (1:nshock_myid,1:ndim))
     allocate(x_pre    (1:nshock_myid,1:ndim))
     allocate(x_pos    (1:nshock_myid,1:ndim))
     allocate(dpre(1:nshock_myid),dpos(1:nshock_myid))
     allocate(ppre(1:nshock_myid),ppos(1:nshock_myid))
     allocate(gpre(1:nshock_myid),gpos(1:nshock_myid))
     allocate(gapre(1:nshock_myid),gapos(1:nshock_myid))
     allocate(epre(1:nshock_myid),epos(1:nshock_myid),cpre(1:nshock_myid))
     allocate(depos(1:nshock_myid))
     allocate(ecpre(1:nshock_myid,1:ncr))
     allocate(ind_shock(1:nshock_myid))
     allocate(ind_shock3(1:nshock),ind_shock3_all(1:nshock))
     allocate(cpu_map_loc(1:nshock),cpu_map_loc_all(1:nshock))
     allocate(x_pre_big(1:nshock,1:ndim))
     allocate(x_pos_big(1:nshock,1:ndim))
     allocate(x_big_all(1:nshock,1:ndim))

     ! =============== Step 3 =================================
     ! Get the shock unit vector for shock cells only.
     ! Get the cell position for shock cells only.
     ! IMPORTANT: Cell index are stored in ind_shock. It will 
     ! be used to get back the corresponding shock cells
     ! ========================================================

     ! This part is to scan only octs having active shock cells
     ngrid=0;ishock=0;ngrid_tot=0
     shock_dir=0.0d0 ! This MUST be set to zero for all dimensions >ndim+1
     do igrid=1,active(ilevel)%ngrid
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           icell=iskip+active(ilevel)%igrid(igrid)
           mach_loc=uold(icell,imach)/uold(icell,1)
           if(mach_loc.gt.mach_strong)then
              ngrid=ngrid+1
              ngrid_tot=ngrid_tot+1
              ind_grid  (ngrid)=active(ilevel)%igrid(igrid)
              ind_shock2(ngrid)=icell
           endif
           if(ngrid.eq.nvector)then
              call get_shock_unit_vector(ind_grid,ngrid,ind_shock,ind_shock2,nshock_myid,ishock,x_shock,shock_dir,dpre,dpos,ppre,ppos,epre,epos,ecpre,gpre,gpos,gapre,gapos,cpre,Bpre,depos,ilevel)
              ngrid=0
           endif
        end do
     end do
     if(ngrid.gt.0)call get_shock_unit_vector(ind_grid,ngrid,ind_shock,ind_shock2,nshock_myid,ishock,x_shock,shock_dir,dpre,dpos,ppre,ppos,epre,epos,ecpre,gpre,gpos,gapre,gapos,cpre,Bpre,depos,ilevel)
     ! =============== Step 3 =================================

     ! With a while loop it should be possible to loop over chunks of shock cells 
     ! (and get rid of some of the dynamical allocated arrays)...

     ! Loop over distance from shock cells by dx increment
     do idist=3,ndist

        do ishock=1,nshock_myid
           delta(1:ndim)=shock_dir(ishock,1:ndim)*dx*dble(idist)
           x_pre(ishock,1:ndim)=(x_shock(ishock,1:ndim)-delta)*scale
           x_pos(ishock,1:ndim)=(x_shock(ishock,1:ndim)+delta)*scale
        enddo
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

        x_pre_big=0.0d0
        x_pos_big=0.0d0
        ind_shock3=0
        cpu_map_loc=0
        if(nshock_cpu(myid).gt.0)then
           ind_shock3(imin:imax)=ind_shock(1:nshock_myid)
           cpu_map_loc(imin:imax)=myid ! use my own cpu map (pb with cpu_map...)
           x_pre_big(imin:imax,1:ndim)=x_pre(1:nshock_myid,1:ndim)
           x_pos_big(imin:imax,1:ndim)=x_pos(1:nshock_myid,1:ndim)
        endif
     
        ! Share the upstream/downstream positions amongst all CPUs
#ifndef WITHOUTMPI
        call MPI_ALLREDUCE(ind_shock3,ind_shock3_all,nshock,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
        ind_shock3=ind_shock3_all
        call MPI_ALLREDUCE(cpu_map_loc,cpu_map_loc_all,nshock,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
        cpu_map_loc=cpu_map_loc_all
        call MPI_ALLREDUCE(x_pre_big,x_big_all,ndim*nshock,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
        x_pre_big=x_big_all
        call MPI_ALLREDUCE(x_pos_big,x_big_all,ndim*nshock,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
        x_pos_big=x_big_all
#endif
  
        ishock_myid=0
        do ishock=1,nshock,nvector
           nsweep=MIN(nvector,nshock-ishock+1)

           ! Now we search for pre- and post-shock cells in all CPUs
           ! The cpu_map(icell)==myid condition tests if the cell 
           ! actually belongs to the current CPU
           ! Post-shock region
           do isweep=1,nsweep
              xtest(isweep,1:ndim)=x_pos_big(isweep+ishock-1,1:ndim)
              do idim=1,ndim
                 if(period(idim).and.xtest(isweep,idim)>boxlen)xtest(isweep,idim)=xtest(isweep,idim)-boxlen
                 if(period(idim).and.xtest(isweep,idim)<0.    )xtest(isweep,idim)=xtest(isweep,idim)+boxlen
              enddo
           enddo
           call get_mycell_index(cell_index_pos,cell_index_father,grid_index,cell_lev,xtest,ilevel,nsweep)
           pres_pos=0.0d0; geff_pos=0.0d0; gamm_pos=0.0d0
           ethe_pos=0.0d0; dens_pos=0.0d0
           do i=1,nsweep
              icell=cell_index_pos(i)
              if(cpu_map(cell_index_father(i))==myid)then
                 d=MAX(uold(icell,1),smallr)
                 oneoverd=1d0/d
                 u=uold(icell,2)*oneoverd
                 v=uold(icell,3)*oneoverd
                 w=uold(icell,4)*oneoverd
                 A=0.5*(uold(icell,6)+uold(icell,nvar+1))
                 B=0.5*(uold(icell,7)+uold(icell,nvar+2))
                 C=0.5*(uold(icell,8)+uold(icell,nvar+3))
                 e=uold(icell,5)-0.5*d*(u**2+v**2+w**2)-0.5*(A**2+B**2+C**2)
#if NENER>0
                 do iener=1,nener
                    e=e-uold(icell,8+iener)
                 end do
#endif
                 pres_pos(i)=e*(gamma-1d0)
                 geff_pos(i)=gamma* e*(gamma-1d0)
                 Ptot=e*(gamma-1d0)
                 ethe_pos(i)=e
                 gamm_pos(i)=gamma*e
                 etot=e
#if NENER>0
                 do iener=1,nener
                    pres_pos(i)=pres_pos(i)+uold(icell,8+iener)*(gamma_rad(iener)-1d0)
                    geff_pos(i)=geff_pos(i)+gamma_rad(iener) *uold(icell,8+iener)*(gamma_rad(iener)-1d0)
                    gamm_pos(i)=gamm_pos(i)+gamma_rad(iener) *uold(icell,8+iener)
                    Ptot=Ptot+uold(icell,8+iener)*(gamma_rad(iener)-1d0)
                    etot=etot+uold(icell,8+iener)
                 end do
#endif
#ifdef CRFLX
                 do iener=1,ncr
                    pres_pos(i)=pres_pos(i)+uold(icell,ind_cr+iener)*(gamma_rad(nener+iener)-1d0)
                    geff_pos(i)=geff_pos(i)+gamma_rad(nener+iener) *uold(icell,ind_cr+iener)*(gamma_rad(nener+iener)-1d0)
                    gamm_pos(i)=gamm_pos(i)+gamma_rad(nener+iener) *uold(icell,ind_cr+iener)
                    Ptot=Ptot+uold(icell,ind_cr+iener)*(gamma_rad(nener+iener)-1d0)
                    etot=etot+uold(icell,ind_cr+iener)
                 end do
#endif
                 ethe_pos(i)=etot ! YD WARNING: testing with total energy, should be called eint
                 geff_pos(i)=geff_pos(i)/Ptot
                 gamm_pos(i)=gamm_pos(i)/etot
                 dens_pos(i)=d
              endif
           enddo

           ! Pre-shock region
           do isweep=1,nsweep
              xtest(isweep,1:ndim)=x_pre_big(isweep+ishock-1,1:ndim)
              do idim=1,ndim
                 if(period(idim).and.xtest(isweep,idim)>boxlen)xtest(isweep,idim)=xtest(isweep,idim)-boxlen
                 if(period(idim).and.xtest(isweep,idim)<0.    )xtest(isweep,idim)=xtest(isweep,idim)+boxlen
              enddo
           enddo
           call get_mycell_index(cell_index_pre,cell_index_father,grid_index,cell_lev,xtest,ilevel,nsweep)
           pres_pre=0.0d0; geff_pre=0.0d0; gamm_pre=0.0d0
           ethe_pre=0.0d0; ecr_pre=0.0d0 ; dens_pre=0.0d0
           csou_pre=0.0d0; Bmag_pre=0.0d0
           do i=1,nsweep
              icell=cell_index_pre(i)
              if(cpu_map(cell_index_father(i))==myid)then
                 d=MAX(uold(icell,1),smallr)
                 oneoverd=1d0/d
                 u=uold(icell,2)*oneoverd
                 v=uold(icell,3)*oneoverd
                 w=uold(icell,4)*oneoverd
                 A=0.5*(uold(icell,6)+uold(icell,nvar+1))
                 B=0.5*(uold(icell,7)+uold(icell,nvar+2))
                 C=0.5*(uold(icell,8)+uold(icell,nvar+3))
                 Bmag_pre(i,1)=A
                 Bmag_pre(i,2)=B
                 Bmag_pre(i,3)=C
                 e=uold(icell,5)-0.5*d*(u**2+v**2+w**2)-0.5*(A**2+B**2+C**2)
#if NENER>0
                 do iener=1,nener
                    e=e-uold(icell,8+iener)
                 end do
#endif
                 pres_pre(i)=e*(gamma-1d0)
                 geff_pre(i)=gamma* e*(gamma-1d0)
                 Ptot=e*(gamma-1d0)
                 ethe_pre(i)=e
                 gamm_pre(i)=gamma*e
                 etot=e
                 csou_pre(i)=sqrt(gamma*(gamma-1d0)*e*oneoverd)
#if NENER>0
                 do iener=1,nener
                    pres_pre(i)=pres_pre(i)+uold(icell,8+iener)*(gamma_rad(iener)-1d0)
                    geff_pre(i)=geff_pre(i)+gamma_rad(iener) *uold(icell,8+iener)*(gamma_rad(iener)-1d0)
                    gamm_pre(i)=gamm_pre(i)+gamma_rad(iener) *uold(icell,8+iener)
                    Ptot=Ptot+uold(icell,8+iener)*(gamma_rad(iener)-1d0)
                    etot=etot+uold(icell,8+iener)
                    ecr_pre(i,iener)=uold(icell,8+iener)
                 end do
#endif
#ifdef CRFLX
                 do iener=1,ncr
                    pres_pre(i)=pres_pre(i)+uold(icell,ind_cr+iener)*(gamma_rad(nener+iener)-1d0)
                    geff_pre(i)=geff_pre(i)+gamma_rad(nener+iener) *uold(icell,ind_cr+iener)*(gamma_rad(nener+iener)-1d0)
                    gamm_pre(i)=gamm_pre(i)+gamma_rad(nener+iener) *uold(icell,ind_cr+iener)
                    Ptot=Ptot+uold(icell,ind_cr+iener)*(gamma_rad(nener+iener)-1d0)
                    etot=etot+uold(icell,ind_cr+iener)
                    ecr_pre(i,iener)=uold(icell,ind_cr+iener)
                 end do
#endif
                 geff_pre(i)=geff_pre(i)/Ptot
                 gamm_pre(i)=gamm_pre(i)/etot
                 dens_pre(i)=d
              endif
           enddo

           ! Scatter the pre-/post-shock gas properties to all CPUs
#ifndef WITHOUTMPI
           call MPI_ALLREDUCE(dens_pre,dens_pre_all,nvector,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
           call MPI_ALLREDUCE(dens_pos,dens_pos_all,nvector,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
           call MPI_ALLREDUCE(pres_pre,pres_pre_all,nvector,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
           call MPI_ALLREDUCE(pres_pos,pres_pos_all,nvector,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
           call MPI_ALLREDUCE(geff_pre,geff_pre_all,nvector,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
           call MPI_ALLREDUCE(geff_pos,geff_pos_all,nvector,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
           call MPI_ALLREDUCE(gamm_pre,gamm_pre_all,nvector,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
           call MPI_ALLREDUCE(gamm_pos,gamm_pos_all,nvector,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
           dens_pre=dens_pre_all
           dens_pos=dens_pos_all
           pres_pre=pres_pre_all
           pres_pos=pres_pos_all
           geff_pre=geff_pre_all
           geff_pos=geff_pos_all
           gamm_pre=gamm_pre_all
           gamm_pos=gamm_pos_all

           call MPI_ALLREDUCE(ethe_pre,ethe_pre_all,nvector      ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
           call MPI_ALLREDUCE(ecr_pre ,ecr_pre_all ,nvector*ncr  ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
           call MPI_ALLREDUCE(ethe_pos,ethe_pos_all,nvector      ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
           call MPI_ALLREDUCE(csou_pre,csou_pre_all,nvector      ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
           call MPI_ALLREDUCE(Bmag_pre,Bmag_pre_all,nvector*3    ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
           ethe_pre=ethe_pre_all
           ecr_pre =ecr_pre_all
           ethe_pos=ethe_pos_all
           csou_pre=csou_pre_all
           Bmag_pre=Bmag_pre_all
#endif

           do i=1,nsweep
              icell=ind_shock3(ishock+i-1)
              if(icell.le.0)then
#ifndef WITHOUTMPI
                 call MPI_ABORT(MPI_COMM_WORLD,1,info)
#else
                 stop
#endif                 
              endif
              ! for a given shared ind_cell, cpu_map is not unique
              ! Test which CPU it belongs with the custom cpu_map_loc
              if(cpu_map_loc(ishock+i-1)==myid)then
                 ishock_myid=ishock_myid+1
                 if(ndim==1)write(*,'(A,4e12.5,I3)')'extra check:',(gamma-1d0)*epos(ishock_myid),(gamma-1d0)*ethe_pos(i),ppos(ishock_myid),pres_pos(i),ishock_myid
                 if(ndim==1)write(*,'(A,6e12.4,2I4)')'First',Bpre(ishock_myid,1:3),Bmag_pre(i,1:3),ishock_myid,idist
                 ! check that the pressure and density are lower
                 if(pres_pre(i).lt.ppre(ishock_myid).and.dens_pre(i).lt.dpre(ishock_myid))then
                    dpre (ishock_myid)=dens_pre(i)
                    ppre (ishock_myid)=pres_pre(i)
                    epre (ishock_myid)=ethe_pre(i)
                    ecpre(ishock_myid,1:ncr)=ecr_pre(i,1:ncr)
                    gpre (ishock_myid)=geff_pre(i)
                    gapre(ishock_myid)=gamm_pre(i)
                    cpre (ishock_myid)=csou_pre(i)
                    Bpre (ishock_myid,1:3)=Bmag_pre(i,1:3)
                    if(ndim==1)write(*,'(A,6e12.4,2I4)')'<<<',Bpre(ishock_myid,1:3),Bmag_pre(i,1:3),ishock_myid,idist
                 endif
                 if((ethe_pos(i)-epos(ishock_myid)).lt.depos(ishock_myid))then
                    if(pres_pos(i).gt.ppos(ishock_myid))then
!!$                    ! check that the pressure and density are higher
                       dpos(ishock_myid)=dens_pos(i)
                       ppos(ishock_myid)=pres_pos(i)
                       epos(ishock_myid)=ethe_pos(i)
                       gpos(ishock_myid)=geff_pos(i)
                       gapos(ishock_myid)=gamm_pos(i)
                    endif
                 else
                    depos(ishock_myid)=0.0d0
                 endif
                 pratio=ppos(ishock_myid)/ppre(ishock_myid)

                 g1=gapre(ishock_myid);g2=gapos(ishock_myid)
                 gg=gpre(ishock_myid)
                 CC=(g1-1d0)*((g2+1d0)*pratio+g2-1d0)
                 mach=sqrt( (pratio-1d0)*CC/gg /(CC-(g1+1d0+(g1-1d0)*pratio)*(g2-1d0)) )
                 if(mach.gt.uold(icell,imach)/uold(icell,1))then
                    uold(icell,imach)=mach*uold(icell,1)

                    if(cr_shock)then
                       mach2=mach*mach
                       gg=gpre(ishock_myid)
                       Rcomp=dpos(ishock_myid)/dpre(ishock_myid)
                       ediss_th=epos(ishock_myid)-epre(ishock_myid)*Rcomp**gamma ! Dissipated thermal energy
#ifdef CRFLX
                       do iener=1,ncr
                          ediss_th=ediss_th-ecpre(ishock_myid,iener)*Rcomp**gamma_rad(nener+iener)
#else                          
                       do iener=1,nener
                          ediss_th=ediss_th-ecpre(ishock_myid,iener)*Rcomp**gamma_rad(iener)
#endif
                       enddo
                       if(cr_shock_mag)then
                          Bnorm2=0.0d0
                          do idim=1,3
                             Bnorm2=Bnorm2+Bpre(ishock_myid,idim)**2
                          enddo
                          Bpre_proj(1:3)=Bpre(ishock_myid,1:3)*shock_dir(ishock_myid,1:3)
                          Bnorm_pre=0.0d0
                          do idim=1,3
                             Bnorm_pre=Bnorm_pre+Bpre_proj(idim)
                          enddo
                          ctheta=ABS(Bnorm_pre)/sqrt(Bnorm2)
                          theta=acos(ctheta)           
                          xi=xi0 * ( tanh( (theta_crit-theta)*oneoverdtheta) + 1d0 )
                          if(xi.lt.0.0d0.or.xi.gt.1.0d0)then
                             write(*,'(A,4e10.3)')'xi=',xi,theta,theta_crit,1d0/oneoverdtheta
                             call clean_stop
                          endif
                       else
                          xi=eta_mach
                       endif
                       if(cr_shock_mach)then
                          pcr_tot=0.0d0
#ifdef CRFLX
                          do iener=1,ncr
                             pcr_tot=pcr_tot+ecpre(ishock_myid,iener)*(gamma_rad(nener+iener)-1d0)
#else                          
                          do iener=1,nener
                             pcr_tot=pcr_tot+ecpre(ishock_myid,iener)*(gamma_rad(iener)-1d0)
#endif
                          enddo
                          fcr=pcr_tot/(epre(ishock_myid)*(gamma-1d0))
                          call interpol_creff(mach,fcr,creff)
                          xi=xi*creff
                       endif
                       if(ediss_th>0.0d0)then
                          ediss_cr(icell)=xi*ediss_th*mach*cpre(ishock_myid)/Rcomp*dtnew(ilevel)/dx_loc
#if NVAR>9+NENER
                          uold(icell,imach+1)=ediss_th*mach*cpre(ishock_myid)/Rcomp*MAX(uold(icell,1),smallr)
#endif
#if NVAR>10+NENER
                          if(cr_shock_mag)uold(icell,imach+2)=theta*MAX(uold(icell,1),smallr)
#endif
#if NVAR>11+NENER
                          if(cr_shock_mach)uold(ind_cell(i),imach+3)=fcr*MAX(uold(ind_cell(i),1),smallr)
#endif
#if NVAR>12+NENER
                          if(cr_shock_mach)uold(ind_cell(i),imach+4)=creff*MAX(uold(ind_cell(i),1),smallr)
#endif
                       endif
!!$                       uold(icell,imach+5)=ediss_cr(icell)*MAX(uold(icell,1),smallr)
                    endif
                 endif

              endif
           enddo
           

        enddo
        
     enddo
     ! End loop over distance steps

     deallocate(x_shock,x_pre,x_pos)
     deallocate(shock_dir)
     deallocate(x_pre_big,x_pos_big,x_big_all)
     deallocate(ind_shock,ind_shock3,ind_shock3_all)
     deallocate(cpu_map_loc,cpu_map_loc_all)
     deallocate(dpre,dpos)
     deallocate(ppre,ppos)
     deallocate(gpre,gpos)
     deallocate(gapre,gapos)
     deallocate(epre,epos,ecpre,cpre,depos,Bpre)
  endif


#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(decr_lost_mpi,decr_lost_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  decr_lost=decr_lost+decr_lost_all
  decr_lost_mpi=0.0d0
  call MPI_ALLREDUCE(decr_gain_mpi,decr_gain_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  decr_gain=decr_gain+decr_gain_all
  decr_gain_mpi=0.0d0
#else
  decr_lost=decr_lost_mpi
  decr_gain=decr_gain_mpi
#endif


111 format('   Entering shock_fine for level ',i2)

end subroutine shock_fine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine get_shock_cells(ind_grid,ncache,ilevel)

  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel,ncache
  integer,dimension(1:nvector)::ind_grid
  !-------------------------------------------------------------------
  !-------------------------------------------------------------------
  integer ,dimension(1:twotondim)::icic,jcic,kcic
  real(dp),dimension(1:twotondim)::vol
  integer ,dimension(1:nvector,1:threetondim     ),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim       ),save::nbors_father_grids
  integer ,dimension(1:nvector,0:twondim         ),save::ibuffer_father
  integer ,dimension(1:nvector,0:twondim         ),save::ind1
  real(dp),dimension(1:nvector,0:twondim  ,1:nvar+3),save::u1
  real(dp),dimension(1:nvector,1:twotondim,1:nvar+3),save::u2

  integer ,dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),save::ind_cell_loc

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3),save::uloc
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),save::divuloc,ddx
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),save::pres,csou,dens,geff,gamm,ethe
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ncr),save::ecr
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndim),save::gradd,gradT,gradS
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3),save::uvel,Bmag
  real(dp),dimension(1:ndim),save::x,dd,dg,position
  real(dp),dimension(1:ndim),save::uvel_pre,uvel_pos,uvel_pre_far,uvel_pos_far
  real(dp),dimension(1:3   ),save::Bmag_pre,Bmag_pos,Bmag_pre_far,Bmag_pos_far
  real(dp),dimension(1:3   ),save::gradT_extended
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2),save::quantity
  real(dp),dimension(1:ncr),save::ecr_pre,ecr_pos

  logical ,dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),save::okshock

  integer,dimension(1:nvector),save::igrid_nbor,ind_cell,ind_buffer,ind_exist,ind_nexist

  integer::ind_buffer1,ind_buffer2,ind_buffer3
  integer::ind_father1,ind_father2,ind_father3
  integer::i,j,ivar,idim,ind_son,ind_father,iskip,nbuffer,ibuffer
  integer::i0,j0,k0,i1,j1,k1,i2,j2,k2,i3,j3,k3,nx_loc,nb_noneigh,nexist
  integer::i1min,i1max,j1min,j1max,k1min,k1max
  integer::i2min,i2max,j2min,j2max,k2min,k2max
  integer::i3min,i3max,j3min,j3max,k3min,k3max
  integer::i4min,i4max,j4min,j4max,k4min,k4max
  integer::icr,irad,ii,jj,kk,ig,id,jg,jd,kg,kd,icell
  integer::ipre,ipos,jpre,jpos,idim_max,ind_cr

  real(dp)::dx,dx_loc,scale,oneoverd,oneoverdeltax
  real(dp)::d,u,v,w,A,B,C,e,mach,mach2,d_up,csound,cr_flux,deltae,machold,machloc,mymax
  real(dp)::Ptot,etot,pratio,Rcomp,ediss,eth_pos,eth_pre
  real(dp)::max_loc,factor,divu_pre,divu_pos,u_up
  real(dp)::dens_pre,dens_pos,dens_pre_far,dens_pos_far
  real(dp)::pres_pre,pres_pos,pres_pre_far,pres_pos_far
  real(dp)::csou_pre,csou_pos,csou_pre_far,csou_pos_far
  real(dp)::geff_pre,geff_pos,geff_pre_far,geff_pos_far
  real(dp)::gamm_pre,gamm_pos,gamm_pre_far,gamm_pos_far
  real(dp)::unorm_pre,unorm_pos,Bnorm2,Bnorm_pre,ctheta,theta,den
  real(dp)::densg,densd,tempg,tempd,entrg,entrd,dTdd,dTdS,normgradT,gamma_eff
  real(dp)::ul,ur,pi,theta_crit,oneoverdtheta,xi,xi0,newde,gg,g1,g2,CC,minmach2
  real(dp)::nearlytwo=1.999d0
  real(dp)::cic_interpol,cic_value
  real(dp)::pcr_tot,fcr,creff

  logical::okdivu

  ind_cr=8
  if(twotemp)ind_cr=9
  if(mom_cr_advect)ind_cr=icru-1

  pi=acos(-1d0)
  theta_crit=0.25d0*pi
  oneoverdtheta=18d0/pi
  xi0=0.5d0*eta_mach
  minmach2=minmach*minmach

  okshock=.true.

  ! Mesh spacing in that level
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx=0.5D0**ilevel
  dx_loc=dx*scale

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

  !------------------------------------------
  ! Gather 3^ndim neighboring father cells
  !------------------------------------------
  do i=1,ncache
     ind_cell(i)=father(ind_grid(i))
  end do
  call get3cubefather(ind_cell,nbors_father_cells,nbors_father_grids,ncache,ilevel)
  
  !---------------------------
  ! Gather 6x6x6 cells stencil
  !---------------------------
  ! Loop over 3x3x3 neighboring father cells
  do k1=k1min,k1max
  do j1=j1min,j1max
  do i1=i1min,i1max
     
     ! Check if neighboring grid exists
     nbuffer=0
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
     
     ! If not, interpolate variables from parent cells
     if(nbuffer>0)then
        call getnborfather(ind_buffer,ibuffer_father,nbuffer,ilevel)
        do j=0,twondim
           do ivar=1,nvar+3
              do i=1,nbuffer
                 u1(i,j,ivar)=uold(ibuffer_father(i,j),ivar)
              end do
           end do
           do i=1,nbuffer
              ind1(i,j)=son(ibuffer_father(i,j))
           end do
        end do
        call interpol_hydro(u1,ind1,u2,nbuffer)
     endif

     ! Loop over 2x2x2 cells
     do k2=k2min,k2max
     do j2=j2min,j2max
     do i2=i2min,i2max

        ind_son=1+i2+2*j2+4*k2
        iskip=ncoarse+(ind_son-1)*ngridmax
        do i=1,nexist
           ind_cell(i)=iskip+igrid_nbor(ind_exist(i))
        end do
        
        i3=1; j3=1; k3=1
        if(ndim>0)i3=1+2*(i1-1)+i2
        if(ndim>1)j3=1+2*(j1-1)+j2
        if(ndim>2)k3=1+2*(k1-1)+k2
        
        ! Gather hydro variables
        do ivar=1,nvar+3
           do i=1,nexist
              uloc(ind_exist(i),i3,j3,k3,ivar)=uold(ind_cell(i),ivar)
           end do
           do i=1,nbuffer
              uloc(ind_nexist(i),i3,j3,k3,ivar)=u2(i,ind_son,ivar)
           end do
        end do

        ! Gather dx
        do i=1,nexist
           ddx(ind_exist(i),i3,j3,k3)=dx_loc
        end do
        do i=1,nbuffer
           if(interpol_type.gt.0)then
              ddx(ind_nexist(i),i3,j3,k3)=dx_loc
           else
              ddx(ind_nexist(i),i3,j3,k3)=1.5d0*dx_loc
           endif
        end do
        
     end do
     end do
     end do
     ! End loop over cells

  end do
  end do
  end do
  ! End loop over neighboring grids


  !Warning: replace ijk2minmax=[0,1] by [1,2]
  ! Integer constants
  i2min=1; i2max=1
  j2min=1; j2max=1
  k2min=1; k2max=1
  if(ndim>0)then
     i2min=1;i2max=2
  end if
  if(ndim>1)then
     j2min=1;j2max=2
  end if
  if(ndim>2)then
     k2min=1;k2max=2
  end if

#if NVAR>8
  ! Set the mach number (and other analysis quantities) to zero everywhere
  do ivar=9,nvar-4*ncr
  do k2=k2min,k2max
  do j2=j2min,j2max
  do i2=i2min,i2max
     do i=1,ncache
        if(ndim==1)ind_son=1+(i2-1)
        if(ndim==2)ind_son=1+(i2-1)+2*(j2-1)
        if(ndim==3)ind_son=1+(i2-1)+2*(j2-1)+4*(k2-1)
        iskip=ncoarse+(ind_son-1)*ngridmax
        ind_cell(i)=iskip+ind_grid(i)
        uold(ind_cell(i),ivar)=0.0d0
     enddo
  enddo
  enddo
  enddo
  enddo
#endif

  ! =====================================================================
  ! Compute gradT and gradS
  ! =====================================================================
  do idim=1,ndim
  if(idim==1)then
     ii=1;jj=0;kk=0 ! jg=jd, kg=kd
  endif
  if(idim==2)then
     ii=0;jj=1;kk=0 ! ig=id, kg=kd
  endif
  if(idim==3)then
     ii=0;jj=0;kk=1 ! ig=id, jg=jd
  endif
  do k2=k2min,k2max
  kg=k2-kk;kd=k2+kk
  do j2=j2min,j2max
  jg=j2-jj;jd=j2+jj
  do i2=i2min,i2max
  ig=i2-ii;id=i2+ii
     do i=1,ncache
        ! Right; Top; Front
        d=max(uloc(i,id,jd,kd,1),smallr)
        oneoverd=1d0/d
        u=uloc(i,id,jd,kd,2)*oneoverd
        v=uloc(i,id,jd,kd,3)*oneoverd
        w=uloc(i,id,jd,kd,4)*oneoverd
        A=0.5*(uloc(i,id,jd,kd,6)+uloc(i,id,jd,kd,nvar+1))
        B=0.5*(uloc(i,id,jd,kd,7)+uloc(i,id,jd,kd,nvar+2))
        C=0.5*(uloc(i,id,jd,kd,8)+uloc(i,id,jd,kd,nvar+3))
        e=uloc(i,id,jd,kd,5)-0.5*d*(u**2+v**2+w**2)-0.5*(A**2+B**2+C**2)
#if NENER>0
        do irad=1,nener
           e=e-uloc(i,id,jd,kd,8+irad)
        end do
#endif
        densd=d
        tempd=e*(gamma-1d0)*oneoverd
        entrd=e*oneoverd**(gamma-1d0)
#if NENER>0
        do irad=1,nener
           tempd=tempd+uloc(i,id,jd,kd,8+irad)*(gamma_rad(irad)-1d0)*oneoverd
           entrd=entrd+uloc(i,id,jd,kd,8+irad)*oneoverd**(gamma_rad(irad)-1d0)
        end do
#endif
#ifdef CRFLX
        do icr=1,ncr
           tempd=tempd+uloc(i,id,jd,kd,ind_cr+irad)*(gamma_cr(icr)-1d0)*oneoverd
           entrd=entrd+uloc(i,id,jd,kd,ind_cr+irad)*oneoverd**(gamma_cr(icr)-1d0)
        end do
#endif

        ! Left; Bottom; Back
        d=max(uloc(i,ig,jg,kg,1),smallr)
        oneoverd=1d0/d
        u=uloc(i,ig,jg,kg,2)*oneoverd
        v=uloc(i,ig,jg,kg,3)*oneoverd
        w=uloc(i,ig,jg,kg,4)*oneoverd
        A=0.5*(uloc(i,ig,jg,kg,6)+uloc(i,ig,jg,kg,nvar+1))
        B=0.5*(uloc(i,ig,jg,kg,7)+uloc(i,ig,jg,kg,nvar+2))
        C=0.5*(uloc(i,ig,jg,kg,8)+uloc(i,ig,jg,kg,nvar+3))
        e=uloc(i,ig,jg,kg,5)-0.5*d*(u**2+v**2+w**2)-0.5*(A**2+B**2+C**2)
#if NENER>0
        do irad=1,nener
           e=e-uloc(i,ig,jg,kg,8+irad) ! Here we keep only the ion temperature (should we keep the mean temperature? with CR?)
        end do
#endif
        densg=d
        tempg=e*(gamma-1d0)*oneoverd
        entrg=e*oneoverd**(gamma-1d0)
#if NENER>0
        do irad=1,nener
           tempg=tempg+uloc(i,ig,jg,kg,8+irad)*(gamma_rad(irad)-1d0)*oneoverd
           entrg=entrg+uloc(i,ig,jg,kg,8+irad)*oneoverd**(gamma_rad(irad)-1d0)
        end do
#endif
#ifdef CRFLX
        do icr=1,ncr
           tempg=tempg+uloc(i,ig,jg,kg,ind_cr+irad)*(gamma_cr(icr)-1d0)*oneoverd
           entrg=entrg+uloc(i,ig,jg,kg,ind_cr+irad)*oneoverd**(gamma_cr(icr)-1d0)
        end do
#endif

        oneoverdeltax=1d0/(ddx(i,id,jd,kd)+ddx(i,ig,jg,kg))

        gradd(i,i2,j2,k2,idim)=(densd-densg)*oneoverdeltax ! Grad range from i=[i2min,i2max]=[1,2], j=[j2min,j2max], k=[k2min,k2max]
        gradT(i,i2,j2,k2,idim)=(tempd-tempg)*oneoverdeltax
        gradS(i,i2,j2,k2,idim)=(entrd-entrg)*oneoverdeltax
     enddo
  enddo
  enddo
  enddo
  enddo   

  ! =====================================================================
  ! Tag cells with div.u <0 (divu>0) where divu is computed here
  ! i.e. not using stored divu
  ! =====================================================================
  divuloc=0.0d0
  ! Change this loop for indexes running over [0,3]
  i4min=1; i4max=1; j4min=1; j4max=1; k4min=1; k4max=1
  if(ndim>0)then
     i4min=0; i4max=3
  endif
  if(ndim>1)then
     j4min=0; j4max=3
  endif
  if(ndim>2)then
     k4min=0; k4max=3
  endif
  do k2=k4min,k4max
  do j2=j4min,j4max
  do i2=i4min,i4max
     do i=1,ncache
        ul=uloc(i,i2-1,j2  ,k2  ,2)/max(uloc(i,i2-1,j2  ,k2  ,1),smallr)
        ur=uloc(i,i2+1,j2  ,k2  ,2)/max(uloc(i,i2+1,j2  ,k2  ,1),smallr)
        divuloc(i,i2,j2,k2)=divuloc(i,i2,j2,k2) &
             & - (ur-ul)/(ddx(i,i2+1,j2  ,k2  )+ddx(i,i2-1,j2  ,k2  ))
#if NDIM>1
        ul=uloc(i,i2  ,j2-1,k2  ,3)/max(uloc(i,i2  ,j2-1,k2  ,1),smallr)
        ur=uloc(i,i2  ,j2+1,k2  ,3)/max(uloc(i,i2  ,j2+1,k2  ,1),smallr)
        divuloc(i,i2,j2,k2)=divuloc(i,i2,j2,k2) &
             & - (ur-ul)/(ddx(i,i2  ,j2+1,k2  )+ddx(i,i2  ,j2-1,k2  ))
#endif
#if NDIM>2
        ul=uloc(i,i2  ,j2  ,k2-1,4)/max(uloc(i,i2  ,j2,k2-1,1),smallr)
        ur=uloc(i,i2  ,j2  ,k2+1,4)/max(uloc(i,i2  ,j2,k2+1,1),smallr)
        divuloc(i,i2,j2,k2)=divuloc(i,i2,j2,k2) &
             & - (ur-ul)/(ddx(i,i2  ,j2  ,k2+1)+ddx(i,i2  ,j2  ,k2-1))
#endif
     enddo
  enddo
  enddo
  enddo

  ! =====================================================================
  ! Tag leaf cells only
  ! =====================================================================
  do k2=k2min,k2max
  do j2=j2min,j2max
  do i2=i2min,i2max
     do i=1,ncache
        ! Get cell index
        if(ndim==1)ind_son=1+(i2-1)
        if(ndim==2)ind_son=1+(i2-1)+2*(j2-1)
        if(ndim==3)ind_son=1+(i2-1)+2*(j2-1)+4*(k2-1)
        iskip=ncoarse+(ind_son-1)*ngridmax
        ind_cell(i)=iskip+ind_grid(i)
        okshock(i,i2,j2,k2)=okshock(i,i2,j2,k2).and.son(ind_cell(i))==0
     enddo
  enddo
  enddo
  enddo   

  ! =====================================================================
  ! Tag cells with div.u <0 (divu>0): regions of compression
  ! =====================================================================
  do k2=k2min,k2max
  do j2=j2min,j2max
  do i2=i2min,i2max
     do i=1,ncache
        okshock(i,i2,j2,k2)=okshock(i,i2,j2,k2).and.divuloc(i,i2,j2,k2)>0.0d0
     enddo
  enddo
  enddo
  enddo   

  ! =====================================================================
  ! Tag cells with gradT.gradS>0 and gradT.gradd>0
  ! =====================================================================
  do k2=k2min,k2max
  do j2=j2min,j2max
  do i2=i2min,i2max
     do i=1,ncache
        dTdS=0.0d0
        do idim=1,ndim
           dTdS=dTdS+gradT(i,i2,j2,k2,idim)*gradS(i,i2,j2,k2,idim)
        enddo
        dTdd=0.0d0
        do idim=1,ndim
           dTdd=dTdd+gradT(i,i2,j2,k2,idim)*gradd(i,i2,j2,k2,idim)
        enddo
        okshock(i,i2,j2,k2)=okshock(i,i2,j2,k2).and.dTdS>0.0d0.and.dTdd>0.0d0
     enddo
  enddo
  enddo
  enddo   


  ! =====================================================================
  ! Find the normal to the shock
  ! if it doesn't exist (gradT=0) prevent cell to be candidate for shock
  ! =====================================================================
  do k2=k2min,k2max
  do j2=j2min,j2max
  do i2=i2min,i2max
     do i=1,ncache
        normgradT=0.0d0
        do idim=1,ndim
           normgradT=normgradT+gradT(i,i2,j2,k2,idim)**2
        enddo
        if(normgradT>0.0d0)then
           normgradT=sqrt(normgradT)
           do idim=1,ndim
              gradT(i,i2,j2,k2,idim)=gradT(i,i2,j2,k2,idim)/normgradT
           enddo           
        else
           okshock(i,i2,j2,k2)=.false.
        endif
     enddo
  enddo
  enddo
  enddo   

  ! =====================================================================
  ! Tag cells with local minimum of div.u along the normal to the shock
  ! =====================================================================
  do k2=k2min,k2max
  do j2=j2min,j2max
  do i2=i2min,i2max
     do i=1,ncache
        normgradT=0.0d0
        do idim=1,ndim
           normgradT=normgradT+gradT(i,i2,j2,k2,idim)**2
        enddo
        if(normgradT>0.0d0)then
           quantity(iu1:iu2,ju1:ju2,ku1:ku2)=divuloc(i,iu1:iu2,ju1:ju2,ku1:ku2)
           position(1:ndim)=-gradT(i,i2,j2,k2,1:ndim)
           divu_pre=cic_interpol(position,quantity,i2,j2,k2)
           position(1:ndim)= gradT(i,i2,j2,k2,1:ndim)
           divu_pos=cic_interpol(position,quantity,i2,j2,k2)

           okdivu=(divuloc(i,i2,j2,k2)>divu_pre) .and. (divuloc(i,i2,j2,k2)>divu_pos)
           okshock(i,i2,j2,k2)=okshock(i,i2,j2,k2).and.okdivu
        endif
     enddo
  enddo
  enddo
  enddo   

  ! =====================================================================
  ! Store some useful quantities
  ! dens: gas density
  ! pres: total pressure (thermal + non-thermal)
  ! ethe: thermal energy
  ! csou: effective sound speed (thermal + non-thermal)
  ! geff: effective adiabatic index
  ! =====================================================================
  do k2=ku1,ku2
  do j2=ju1,ju2
  do i2=iu1,iu2
     do i=1,ncache
        d=max(uloc(i,i2,j2,k2,1),smallr)
        oneoverd=1d0/d
        u=uloc(i,i2,j2,k2,2)*oneoverd
        v=uloc(i,i2,j2,k2,3)*oneoverd
        w=uloc(i,i2,j2,k2,4)*oneoverd
        A=0.5*(uloc(i,i2,j2,k2,6)+uloc(i,i2,j2,k2,nvar+1))
        B=0.5*(uloc(i,i2,j2,k2,7)+uloc(i,i2,j2,k2,nvar+2))
        C=0.5*(uloc(i,i2,j2,k2,8)+uloc(i,i2,j2,k2,nvar+3))
        e=uloc(i,i2,j2,k2,5)-0.5*d*(u**2+v**2+w**2)-0.5*(A**2+B**2+C**2)
#if NENER>0
        do irad=1,nener
           e=e-uloc(i,i2,j2,k2,8+irad)
        end do
#endif
        ethe(i,i2,j2,k2)=e
        pres(i,i2,j2,k2)=e*(gamma-1d0)
        csou(i,i2,j2,k2)=gamma*(gamma-1d0)*e*oneoverd ! cs^2
        geff(i,i2,j2,k2)=gamma* e*(gamma-1d0)
        gamm(i,i2,j2,k2)=gamma* e
        uvel(i,i2,j2,k2,1)=u
        uvel(i,i2,j2,k2,2)=v
        uvel(i,i2,j2,k2,3)=w
        Bmag(i,i2,j2,k2,1)=A
        Bmag(i,i2,j2,k2,2)=B
        Bmag(i,i2,j2,k2,3)=C
        Ptot=e*(gamma-1d0)
        etot=e
#if NENER>0
        do irad=1,nener
           pres(i,i2,j2,k2)=pres(i,i2,j2,k2)+uloc(i,i2,j2,k2,8+irad)*(gamma_rad(irad)-1d0)
           csou(i,i2,j2,k2)=csou(i,i2,j2,k2)+gamma_rad(irad)*(gamma_rad(irad)-1d0)*uloc(i,i2,j2,k2,8+irad)*oneoverd
           geff(i,i2,j2,k2)=geff(i,i2,j2,k2)+gamma_rad(irad) *uloc(i,i2,j2,k2,8+irad)*(gamma_rad(irad)-1d0)
           gamm(i,i2,j2,k2)=gamm(i,i2,j2,k2)+gamma_rad(irad) *uloc(i,i2,j2,k2,8+irad)
           Ptot=Ptot+uloc(i,i2,j2,k2,8+irad)*(gamma_rad(irad)-1d0)
           etot=etot+uloc(i,i2,j2,k2,8+irad)
           ecr(i,i2,j2,k2,irad)=uloc(i,i2,j2,k2,8+irad)
        end do
#endif
#ifdef CRFLX
        do icr=1,ncr
           pres(i,i2,j2,k2)=pres(i,i2,j2,k2)+uloc(i,i2,j2,k2,ind_cr+irad)*(gamma_cr(icr)-1d0)
           csou(i,i2,j2,k2)=csou(i,i2,j2,k2)+gamma_cr(icr)*(gamma_cr(icr)-1d0)*uloc(i,i2,j2,k2,ind_cr+irad)*oneoverd
           geff(i,i2,j2,k2)=geff(i,i2,j2,k2)+gamma_cr(icr) *uloc(i,i2,j2,k2,ind_cr+irad)*(gamma_cr(icr)-1d0)
           gamm(i,i2,j2,k2)=gamm(i,i2,j2,k2)+gamma_cr(icr) *uloc(i,i2,j2,k2,ind_cr+irad)
           Ptot=Ptot+uloc(i,i2,j2,k2,ind_cr+irad)*(gamma_cr(icr)-1d0)
           etot=etot+uloc(i,i2,j2,k2,ind_cr+irad)
           ecr(i,i2,j2,k2,irad)=uloc(i,i2,j2,k2,ind_cr+irad)
        end do
#endif
        ethe(i,i2,j2,k2)=etot ! change for internal energy
        csou(i,i2,j2,k2)=sqrt(csou(i,i2,j2,k2))
        geff(i,i2,j2,k2)=geff(i,i2,j2,k2)/Ptot
        gamm(i,i2,j2,k2)=gamm(i,i2,j2,k2)/etot
        dens(i,i2,j2,k2)=d
     enddo
  enddo
  enddo
  enddo   

  ! =====================================================================
  ! Compute the Mach number from pressure jump
  ! =====================================================================
  do k2=k2min,k2max
  do j2=j2min,j2max
  do i2=i2min,i2max
     do i=1,ncache
        if(okshock(i,i2,j2,k2))then
        divu_pre=0.0d0;divu_pos=0.0d0

        normgradT=0.0d0
        do idim=1,ndim
           normgradT=normgradT+gradT(i,i2,j2,k2,idim)**2
        enddo
        if(normgradT>0.0d0)then

           position(1:ndim)=-gradT(i,i2,j2,k2,1:ndim)
           call get_cloud_cic(position,icic,jcic,kcic,vol)
           icic=icic+i2;jcic=jcic+j2;kcic=kcic+k2
           quantity(iu1:iu2,ju1:ju2,ku1:ku2)=dens(i,iu1:iu2,ju1:ju2,ku1:ku2)
           dens_pre=cic_value(quantity,icic,jcic,kcic,vol)
           quantity(iu1:iu2,ju1:ju2,ku1:ku2)=pres(i,iu1:iu2,ju1:ju2,ku1:ku2)
           pres_pre=cic_value(quantity,icic,jcic,kcic,vol)
           quantity(iu1:iu2,ju1:ju2,ku1:ku2)=csou(i,iu1:iu2,ju1:ju2,ku1:ku2)
           csou_pre=cic_value(quantity,icic,jcic,kcic,vol)
           quantity(iu1:iu2,ju1:ju2,ku1:ku2)=geff(i,iu1:iu2,ju1:ju2,ku1:ku2)
           geff_pre=cic_value(quantity,icic,jcic,kcic,vol)
           quantity(iu1:iu2,ju1:ju2,ku1:ku2)=gamm(i,iu1:iu2,ju1:ju2,ku1:ku2)
           gamm_pre=cic_value(quantity,icic,jcic,kcic,vol)
           do idim=1,ndim
              quantity(iu1:iu2,ju1:ju2,ku1:ku2)=uvel(i,iu1:iu2,ju1:ju2,ku1:ku2,idim)
              uvel_pre(idim)=cic_value(quantity,icic,jcic,kcic,vol)
           enddo
           do idim=1,3
              quantity(iu1:iu2,ju1:ju2,ku1:ku2)=Bmag(i,iu1:iu2,ju1:ju2,ku1:ku2,idim)
              Bmag_pre(idim)=cic_value(quantity,icic,jcic,kcic,vol)
           enddo

           position(1:ndim)= gradT(i,i2,j2,k2,1:ndim)
           call get_cloud_cic(position,icic,jcic,kcic,vol)
           icic=icic+i2;jcic=jcic+j2;kcic=kcic+k2
           quantity(iu1:iu2,ju1:ju2,ku1:ku2)=dens(i,iu1:iu2,ju1:ju2,ku1:ku2)
           dens_pos=cic_value(quantity,icic,jcic,kcic,vol)
           quantity(iu1:iu2,ju1:ju2,ku1:ku2)=pres(i,iu1:iu2,ju1:ju2,ku1:ku2)
           pres_pos=cic_value(quantity,icic,jcic,kcic,vol)
           quantity(iu1:iu2,ju1:ju2,ku1:ku2)=csou(i,iu1:iu2,ju1:ju2,ku1:ku2)
           csou_pos=cic_value(quantity,icic,jcic,kcic,vol)
           quantity(iu1:iu2,ju1:ju2,ku1:ku2)=geff(i,iu1:iu2,ju1:ju2,ku1:ku2)
           geff_pos=cic_value(quantity,icic,jcic,kcic,vol)
           quantity(iu1:iu2,ju1:ju2,ku1:ku2)=gamm(i,iu1:iu2,ju1:ju2,ku1:ku2)
           gamm_pos=cic_value(quantity,icic,jcic,kcic,vol)
           do idim=1,ndim
              quantity(iu1:iu2,ju1:ju2,ku1:ku2)=uvel(i,iu1:iu2,ju1:ju2,ku1:ku2,idim)
              uvel_pos(idim)=cic_value(quantity,icic,jcic,kcic,vol)
           enddo
           do idim=1,3
              quantity(iu1:iu2,ju1:ju2,ku1:ku2)=Bmag(i,iu1:iu2,ju1:ju2,ku1:ku2,idim)
              Bmag_pos(idim)=cic_value(quantity,icic,jcic,kcic,vol)
           enddo

           position(1:ndim)=-gradT(i,i2,j2,k2,1:ndim)*nearlytwo
           call get_cloud_cic(position,icic,jcic,kcic,vol)
           icic=icic+i2;jcic=jcic+j2;kcic=kcic+k2
           quantity(iu1:iu2,ju1:ju2,ku1:ku2)=dens(i,iu1:iu2,ju1:ju2,ku1:ku2)
           dens_pre_far=cic_value(quantity,icic,jcic,kcic,vol)
           quantity(iu1:iu2,ju1:ju2,ku1:ku2)=pres(i,iu1:iu2,ju1:ju2,ku1:ku2)
           pres_pre_far=cic_value(quantity,icic,jcic,kcic,vol)
           quantity(iu1:iu2,ju1:ju2,ku1:ku2)=csou(i,iu1:iu2,ju1:ju2,ku1:ku2)
           csou_pre_far=cic_value(quantity,icic,jcic,kcic,vol)
           quantity(iu1:iu2,ju1:ju2,ku1:ku2)=geff(i,iu1:iu2,ju1:ju2,ku1:ku2)
           geff_pre_far=cic_value(quantity,icic,jcic,kcic,vol)
           quantity(iu1:iu2,ju1:ju2,ku1:ku2)=gamm(i,iu1:iu2,ju1:ju2,ku1:ku2)
           gamm_pre_far=cic_value(quantity,icic,jcic,kcic,vol)
           do idim=1,ndim
              quantity(iu1:iu2,ju1:ju2,ku1:ku2)=uvel(i,iu1:iu2,ju1:ju2,ku1:ku2,idim)
              uvel_pre_far(idim)=cic_value(quantity,icic,jcic,kcic,vol)
           enddo
           do idim=1,3
              quantity(iu1:iu2,ju1:ju2,ku1:ku2)=Bmag(i,iu1:iu2,ju1:ju2,ku1:ku2,idim)
              Bmag_pre_far(idim)=cic_value(quantity,icic,jcic,kcic,vol)
           enddo

           position(1:ndim)= gradT(i,i2,j2,k2,1:ndim)*nearlytwo
           call get_cloud_cic(position,icic,jcic,kcic,vol)
           icic=icic+i2;jcic=jcic+j2;kcic=kcic+k2
           quantity(iu1:iu2,ju1:ju2,ku1:ku2)=dens(i,iu1:iu2,ju1:ju2,ku1:ku2)
           dens_pos_far=cic_value(quantity,icic,jcic,kcic,vol)
           quantity(iu1:iu2,ju1:ju2,ku1:ku2)=pres(i,iu1:iu2,ju1:ju2,ku1:ku2)
           pres_pos_far=cic_value(quantity,icic,jcic,kcic,vol)
           quantity(iu1:iu2,ju1:ju2,ku1:ku2)=csou(i,iu1:iu2,ju1:ju2,ku1:ku2)
           csou_pos_far=cic_value(quantity,icic,jcic,kcic,vol)
           quantity(iu1:iu2,ju1:ju2,ku1:ku2)=geff(i,iu1:iu2,ju1:ju2,ku1:ku2)
           geff_pos_far=cic_value(quantity,icic,jcic,kcic,vol)
           quantity(iu1:iu2,ju1:ju2,ku1:ku2)=gamm(i,iu1:iu2,ju1:ju2,ku1:ku2)
           gamm_pos_far=cic_value(quantity,icic,jcic,kcic,vol)
           do idim=1,ndim
              quantity(iu1:iu2,ju1:ju2,ku1:ku2)=uvel(i,iu1:iu2,ju1:ju2,ku1:ku2,idim)
              uvel_pos_far(idim)=cic_value(quantity,icic,jcic,kcic,vol)
           enddo
           do idim=1,3
              quantity(iu1:iu2,ju1:ju2,ku1:ku2)=Bmag(i,iu1:iu2,ju1:ju2,ku1:ku2,idim)
              Bmag_pos_far(idim)=cic_value(quantity,icic,jcic,kcic,vol)
           enddo

           if (pres_pre .lt. pres_pre_far)then
              gamma_eff=geff_pre
              g1=gamm_pre
           else
              gamma_eff=geff_pre_far
              g1=gamm_pre_far
              dens_pre=dens_pre_far
              pres_pre=pres_pre_far
              csou_pre=csou_pre_far
              uvel_pre(1:ndim)=uvel_pre_far(1:ndim)
              Bmag_pre(1:3   )=Bmag_pre_far(1:3   )
           endif
           if (pres_pos .lt. pres_pos_far)then
              g2=gamm_pos_far
              dens_pos=dens_pos_far
              pres_pos=pres_pos_far
              csou_pos=csou_pos_far
              uvel_pos(1:ndim)=uvel_pos_far(1:ndim)
              Bmag_pos(1:3   )=Bmag_pos_far(1:3   )
           else
              g2=gamm_pos
           endif
           
           if(pres_pre.ne.0.0d0.and.pres_pos.ne.0.0d0)then
              pratio=MAX(pres_pos,pres(i,i2,j2,k2))/MIN(pres_pre,pres(i,i2,j2,k2))
           else
              pratio=1.0d0
           endif
           gg=gamma_eff
           CC=(g1-1d0)*((g2+1d0)*pratio+g2-1d0)
           den=CC-(g1+1d0+(g1-1d0)*pratio)*(g2-1d0)
           
           if(den.gt.0.0d0)then

           mach2=(pratio-1d0)*CC/gg /den
           if(mach2.gt.minmach2)then

           mach=sqrt(mach2)

           Rcomp=dens_pos/dens_pre

           if(pres_pos.gt.pres_pos_far)then
              position(1:ndim)= gradT(i,i2,j2,k2,1:ndim)
           else
              position(1:ndim)= gradT(i,i2,j2,k2,1:ndim)*nearlytwo
           endif
           quantity(iu1:iu2,ju1:ju2,ku1:ku2)=ethe(i,iu1:iu2,ju1:ju2,ku1:ku2)
           eth_pos=cic_interpol(position,quantity,i2,j2,k2)
           ediss=eth_pos

           if(pres_pre.lt.pres_pre_far)then
              position(1:ndim)=-gradT(i,i2,j2,k2,1:ndim)
           else
              position(1:ndim)=-gradT(i,i2,j2,k2,1:ndim)*nearlytwo
           endif
           quantity(iu1:iu2,ju1:ju2,ku1:ku2)=ethe(i,iu1:iu2,ju1:ju2,ku1:ku2)
           eth_pre=cic_interpol(position,quantity,i2,j2,k2)
           ediss=ediss-eth_pre*Rcomp**gamma
#if NENER>0
           do irad=1,nener
              quantity(iu1:iu2,ju1:ju2,ku1:ku2)=ecr(i,iu1:iu2,ju1:ju2,ku1:ku2,irad)
              ecr_pre(irad)=cic_interpol(position,quantity,i2,j2,k2)
              ediss=ediss-ecr_pre(irad)*Rcomp**gamma_rad(irad)
           enddo
#endif

           csound=MIN(csou(i,i2,j2,k2),csou_pre)
!!$           if(mach.gt.mach_strong)then
!!$              uvel_pre(1:ndim)=(uvel_pre(1:ndim)*gradT(i,i2,j2,k2,1:ndim))**2
!!$              uvel_pos(1:ndim)=(uvel_pos(1:ndim)*gradT(i,i2,j2,k2,1:ndim))**2
!!$              unorm_pre=sqrt(SUM(uvel_pre))
!!$              unorm_pos=sqrt(SUM(uvel_pos))
!!$              u_up=ABS(4d0*(unorm_pre-unorm_pos)/3d0)
!!$           else
!!$              u_up=mach*csound
!!$           endif
           u_up=mach*csound

           if(cr_shock_mag)then
              Bnorm2=0.0d0
              do idim=1,3
                 Bnorm2=Bnorm2+Bmag_pre(idim)**2
              enddo
              gradT_extended=0.0d0
              gradT_extended(1:ndim)=gradT(i,i2,j2,k2,1:ndim)
              Bmag_pre(1:3)=Bmag_pre(1:3)*gradT_extended(1:3)
              Bnorm_pre=ABS(SUM(Bmag_pre))
              ctheta=Bnorm_pre/sqrt(Bnorm2)
              theta=acos(ctheta)           
              xi=xi0 * ( tanh( (theta_crit-theta)*oneoverdtheta) + 1d0 )
              if(xi.lt.0.0d0.or.xi.gt.1.0d0)then
                 write(*,'(A,4e10.3)')'xi=',xi,theta,theta_crit,1d0/oneoverdtheta
                 call clean_stop
              endif
           else
              xi=eta_mach
           endif
           if(cr_shock_mach)then
              pcr_tot=0.0d0
              do icr=1,ncr
                 pcr_tot=pcr_tot+ecr_pre(irad)*(gamma_rad(irad)-1d0)
              enddo
              fcr=pcr_tot/(eth_pre*(gamma-1d0))
              call interpol_creff(mach,fcr,creff)
              xi=xi*creff
           endif

!!$           cr_flux=xi*dens_pre*u_up**3*0.5d0
           cr_flux=xi*ediss*u_up/Rcomp
           newde=cr_flux*dtnew(ilevel)/dx_loc

           ! Get cell index
           if(ndim==1)ind_son=1+(i2-1)
           if(ndim==2)ind_son=1+(i2-1)+2*(j2-1)
           if(ndim==3)ind_son=1+(i2-1)+2*(j2-1)+4*(k2-1)
           iskip=ncoarse+(ind_son-1)*ngridmax
           ind_cell(i)=iskip+ind_grid(i)

           if(ediss.gt.0.0d0)then
              okshock(i,i2,j2,k2)=.true.
              uold(ind_cell(i),imach)=mach*MAX(uold(ind_cell(i),1),smallr)
#if NVAR>9+NENER
              uold(ind_cell(i),imach+1)=ediss*u_up/Rcomp*MAX(uold(ind_cell(i),1),smallr)
#endif
#if NVAR>10+NENER
              if(cr_shock_mag)uold(ind_cell(i),imach+2)=theta*MAX(uold(ind_cell(i),1),smallr)
#endif
#if NVAR>11+NENER
              if(cr_shock_mach)uold(ind_cell(i),imach+3)=fcr*MAX(uold(ind_cell(i),1),smallr)
#endif
#if NVAR>12+NENER
              if(cr_shock_mach)uold(ind_cell(i),imach+4)=creff*MAX(uold(ind_cell(i),1),smallr)
#endif

              if(cr_shock.and.t.gt.tmin_acc)then 
                 if(mach.gt.mach_strong.and.ndist.gt.2)then
                    ! In that case the actual ediss_cr will be computed later
                    ! with 3 cells distance or more
                    ediss_cr(ind_cell(i))=0.0d0
!!$                    uold(icell,imach+5)=0.0d0
                 else
                    decr_gain_mpi=decr_gain_mpi+newde*dx_loc**ndim
                    ediss_cr(ind_cell(i))=newde
!!$                    uold(icell,imach+5)=newde*MAX(uold(icell,1),smallr)
                 endif
              endif


           endif
           else
              okshock(i,i2,j2,k2)=.false.
           endif
           endif
        endif
     endif
     enddo
  enddo
  enddo
  enddo

  ! =====================================================================
  ! 
  ! =====================================================================
  if(cr_shock)then
  do k2=k2min,k2max
  do j2=j2min,j2max
  do i2=i2min,i2max
     do i=1,ncache
        if(.not.okshock(i,i2,j2,k2))then
           if(ndim==1)ind_son=1+(i2-1)
           if(ndim==2)ind_son=1+(i2-1)+2*(j2-1)
           if(ndim==3)ind_son=1+(i2-1)+2*(j2-1)+4*(k2-1)
           iskip=ncoarse+(ind_son-1)*ngridmax
           ind_cell(i)=iskip+ind_grid(i)
           ediss_cr(ind_cell(i))=0.0d0
        endif
     enddo
  enddo
  enddo
  enddo
  endif

111 format('   Entering get_shock_cells for level ',i2)

end subroutine get_shock_cells
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine get_shock_unit_vector(ind_grid,ncache,ind_shock,ind_shock2,nshock,ishock,x_shock,shock_dir,dpre,dpos,ppre,ppos,epre,epos,ecpre,gpre,gpos,gapre,gapos,cpre,Bpre,depos,ilevel)

  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel,ncache,nshock,ishock
  integer ,dimension(1:nvector                )::ind_grid,ind_shock2
  integer ,dimension(1:nshock                 )::ind_shock
  real(dp),dimension(1:nshock,1:ndim)::x_shock
  real(dp),dimension(1:nshock,1:3)::shock_dir,Bpre
  real(dp),dimension(1:nshock)::dpre,dpos,ppre,ppos,epre,epos
  real(dp),dimension(1:nshock,1:ncr)::ecpre
  real(dp),dimension(1:nshock)::gpre,gpos,gapre,gapos,cpre,depos
  !-------------------------------------------------------------------
  !-------------------------------------------------------------------
  integer ,dimension(1:nvector,1:threetondim     ),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim       ),save::nbors_father_grids
  integer ,dimension(1:nvector,0:twondim         ),save::ibuffer_father
  integer ,dimension(1:nvector,0:twondim         ),save::ind1

  real(dp),dimension(1:nvector,0:twondim  ,1:nvar+3),save::u1
  real(dp),dimension(1:nvector,1:twotondim,1:nvar+3),save::u2

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3),save::uloc
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),save::ddx,dens,pres,geff,gamm
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),save::ethe,csou
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ncr),save::ecr
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndim),save::gradT
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3   ),save::Bmag
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2),save::quantity
  real(dp),dimension(1:ndim),save::position
  real(dp),dimension(1:3)::skip_loc,xcell
  real(dp),dimension(1:3)::Bmag_pre,Bmag_pre_far
  real(dp),dimension(1:ncr)::ecr_pre,ecr_pre_far

  logical ,dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),save::okshock

  integer,dimension(1:nvector),save::igrid_nbor,ind_cell,ind_buffer,ind_exist,ind_nexist

  integer::ind_buffer1,ind_buffer2,ind_buffer3
  integer::ind_father1,ind_father2,ind_father3
  integer::i,j,ivar,idim,ind_son,ind_father,iskip,nbuffer,ibuffer
  integer::i0,j0,k0,i1,j1,k1,i2,j2,k2,i3,j3,k3,nx_loc,nb_noneigh,nexist
  integer::i1min,i1max,j1min,j1max,k1min,k1max
  integer::i2min,i2max,j2min,j2max,k2min,k2max
  integer::i3min,i3max,j3min,j3max,k3min,k3max
  integer::icr,irad,ii,jj,kk,ig,id,jg,jd,kg,kd,icell,ind_cr

  real(dp)::dx,dx_loc,scale,oneoverd,oneoverdeltax
  real(dp)::d,u,v,w,A,B,C,e,etot
  real(dp)::tempg,tempd,normgradT,delta
  real(dp)::dens_pre,dens_pos,dens_pre_far,dens_pos_far
  real(dp)::pres_pre,pres_pos,pres_pre_far,pres_pos_far
  real(dp)::geff_pre,geff_pre_far,geff_pos,geff_pos_far
  real(dp)::gamm_pre,gamm_pre_far,gamm_pos,gamm_pos_far
  real(dp)::ethe_pre,ethe_pos,ethe_pre_far,ethe_pos_far
  real(dp)::csou_pre,csou_pre_far
  real(dp)::cic_interpol,nearlytwo=1.999d0

  ind_cr=8
  if(twotemp)ind_cr=9
  if(mom_cr_advect)ind_cr=icru-1

  ! Mesh spacing in that level
  nx_loc=icoarse_max-icoarse_min+1
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx=0.5D0**ilevel
  dx_loc=dx*scale

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

  !------------------------------------------
  ! Gather 3^ndim neighboring father cells
  !------------------------------------------
  do i=1,ncache
     ind_cell(i)=father(ind_grid(i))
  end do
  call get3cubefather(ind_cell,nbors_father_cells,nbors_father_grids,ncache,ilevel)
  
  !---------------------------
  ! Gather 6x6x6 cells stencil
  !---------------------------
  ! Loop over 3x3x3 neighboring father cells
  do k1=k1min,k1max
  do j1=j1min,j1max
  do i1=i1min,i1max
     
     ! Check if neighboring grid exists
     nbuffer=0
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
     
     ! If not, interpolate variables from parent cells
     if(nbuffer>0)then
        call getnborfather(ind_buffer,ibuffer_father,nbuffer,ilevel)
        do j=0,twondim
           do ivar=1,nvar+3
              do i=1,nbuffer
                 u1(i,j,ivar)=uold(ibuffer_father(i,j),ivar)
              end do
           end do
           do i=1,nbuffer
              ind1(i,j)=son(ibuffer_father(i,j))
           end do
        end do
        call interpol_hydro(u1,ind1,u2,nbuffer)
     endif

     ! Loop over 2x2x2 cells
     do k2=k2min,k2max
     do j2=j2min,j2max
     do i2=i2min,i2max

        ind_son=1+i2+2*j2+4*k2
        iskip=ncoarse+(ind_son-1)*ngridmax
        do i=1,nexist
           ind_cell(i)=iskip+igrid_nbor(ind_exist(i))
        end do
        
        i3=1; j3=1; k3=1
        if(ndim>0)i3=1+2*(i1-1)+i2
        if(ndim>1)j3=1+2*(j1-1)+j2
        if(ndim>2)k3=1+2*(k1-1)+k2
        
        ! Gather hydro variables
        do ivar=1,nvar+3
           do i=1,nexist
              uloc(ind_exist(i),i3,j3,k3,ivar)=uold(ind_cell(i),ivar)
           end do
           do i=1,nbuffer
              uloc(ind_nexist(i),i3,j3,k3,ivar)=u2(i,ind_son,ivar)
           end do
        end do

        ! Gather divu
        do i=1,nexist
           ddx(ind_exist(i),i3,j3,k3)=dx_loc
        end do
        do i=1,nbuffer
           if(interpol_type.gt.0)then
              ddx(ind_nexist(i),i3,j3,k3)=dx_loc
           else
              ddx(ind_nexist(i),i3,j3,k3)=1.5d0*dx_loc
           endif
        end do
        
     end do
     end do
     end do
     ! End loop over cells

  end do
  end do
  end do
  ! End loop over neighboring grids


  ! Warning: replace ijk2minmax=[0,1] by [1,2]
  ! Integer constants
  i2min=1; i2max=1
  j2min=1; j2max=1
  k2min=1; k2max=1
  if(ndim>0)then
     i1min=1;i2max=2
  end if
  if(ndim>1)then
     j2min=1;j2max=2
  end if
  if(ndim>2)then
     k2min=1;k2max=2
  end if

  okshock=.false.
  do k2=k2min,k2max
  do j2=j2min,j2max
  do i2=i2min,i2max
     do i=1,ncache
        if(ndim==1)ind_son=1+(i2-1)
        if(ndim==2)ind_son=1+(i2-1)+2*(j2-1)
        if(ndim==3)ind_son=1+(i2-1)+2*(j2-1)+4*(k2-1)
        iskip=ncoarse+(ind_son-1)*ngridmax
        icell=iskip+ind_grid(i)
        if(icell==ind_shock2(i))okshock(i,i2,j2,k2)=.true.
     enddo
  enddo
  enddo
  enddo

  ! =====================================================================
  ! Compute gradT and gradS
  ! =====================================================================
  do idim=1,ndim
  if(idim==1)then
     ii=1;jj=0;kk=0 ! jg=jd, kg=kd
  endif
  if(idim==2)then
     ii=0;jj=1;kk=0 ! ig=id, kg=kd
  endif
  if(idim==3)then
     ii=0;jj=0;kk=1 ! ig=id, jg=jd
  endif
  do k2=k2min,k2max
  kg=k2-kk;kd=k2+kk
  do j2=j2min,j2max
  jg=j2-jj;jd=j2+jj
  do i2=i2min,i2max
  ig=i2-ii;id=i2+ii
     do i=1,ncache
        if(okshock(i,i2,j2,k2))then
        ! Right; Top; Front
        d=max(uloc(i,id,jd,kd,1),smallr)
        oneoverd=1d0/d
        u=uloc(i,id,jd,kd,2)*oneoverd
        v=uloc(i,id,jd,kd,3)*oneoverd
        w=uloc(i,id,jd,kd,4)*oneoverd
        A=0.5*(uloc(i,id,jd,kd,6)+uloc(i,id,jd,kd,nvar+1))
        B=0.5*(uloc(i,id,jd,kd,7)+uloc(i,id,jd,kd,nvar+2))
        C=0.5*(uloc(i,id,jd,kd,8)+uloc(i,id,jd,kd,nvar+3))
        e=uloc(i,id,jd,kd,5)-0.5*d*(u**2+v**2+w**2)-0.5*(A**2+B**2+C**2)
#if NENER>0
        do irad=1,nener
           e=e-uloc(i,id,jd,kd,8+irad)
        end do
#endif
        tempd=e*(gamma-1d0)*oneoverd
#if NENER>0
        do irad=1,nener
           tempd=tempd+uloc(i,id,jd,kd,8+irad)*(gamma_rad(irad)-1d0)*oneoverd
        end do
#endif
#ifdef CRFLX
        do icr=1,ncr
           tempd=tempd+uloc(i,id,jd,kd,ind_cr+irad)*(gamma_cr(icr)-1d0)*oneoverd
        end do
#endif

        ! Left; Bottom; Back
        d=max(uloc(i,ig,jg,kg,1),smallr)
        oneoverd=1d0/d
        u=uloc(i,ig,jg,kg,2)*oneoverd
        v=uloc(i,ig,jg,kg,3)*oneoverd
        w=uloc(i,ig,jg,kg,4)*oneoverd
        A=0.5*(uloc(i,ig,jg,kg,6)+uloc(i,ig,jg,kg,nvar+1))
        B=0.5*(uloc(i,ig,jg,kg,7)+uloc(i,ig,jg,kg,nvar+2))
        C=0.5*(uloc(i,ig,jg,kg,8)+uloc(i,ig,jg,kg,nvar+3))
        e=uloc(i,ig,jg,kg,5)-0.5*d*(u**2+v**2+w**2)-0.5*(A**2+B**2+C**2)
#if NENER>0
        do irad=1,nener
           e=e-uloc(i,ig,jg,kg,8+irad)
        end do
#endif
        tempg=e*(gamma-1d0)*oneoverd
#if NENER>0
        do irad=1,nener
           tempg=tempg+uloc(i,ig,jg,kg,8+irad)*(gamma_rad(irad)-1d0)*oneoverd
        end do
#endif
#ifdef CRFLX
        do icr=1,ncr
           tempg=tempg+uloc(i,ig,jg,kg,ind_cr+irad)*(gamma_cr(icr)-1d0)*oneoverd
        end do
#endif

        oneoverdeltax=1d0/(ddx(i,id,jd,kd)+ddx(i,ig,jg,kg))
        gradT(i,i2,j2,k2,idim)=(tempd-tempg)*oneoverdeltax

     endif
     enddo
  enddo
  enddo
  enddo
  enddo   

  ! =====================================================================
  ! Store some useful quantities
  ! pres: total pressure (thermal + non-thermal)
  ! It has to be done for all cells because of CIC interpolation
  ! =====================================================================
  do k2=ku1,ku2
  do j2=ju1,ju2
  do i2=iu1,iu2
     do i=1,ncache
        d=max(uloc(i,i2,j2,k2,1),smallr)
        oneoverd=1d0/d
        u=uloc(i,i2,j2,k2,2)*oneoverd
        v=uloc(i,i2,j2,k2,3)*oneoverd
        w=uloc(i,i2,j2,k2,4)*oneoverd
        A=0.5*(uloc(i,i2,j2,k2,6)+uloc(i,i2,j2,k2,nvar+1))
        B=0.5*(uloc(i,i2,j2,k2,7)+uloc(i,i2,j2,k2,nvar+2))
        C=0.5*(uloc(i,i2,j2,k2,8)+uloc(i,i2,j2,k2,nvar+3))
        e=uloc(i,i2,j2,k2,5)-0.5*d*(u**2+v**2+w**2)-0.5*(A**2+B**2+C**2)
#if NENER>0
        do irad=1,nener
           e=e-uloc(i,i2,j2,k2,8+irad)
        end do
#endif
        ethe(i,i2,j2,k2)=e
        pres(i,i2,j2,k2)=e*(gamma-1d0)
        csou(i,i2,j2,k2)=gamma*(gamma-1d0)*e*oneoverd ! cs^2
        geff(i,i2,j2,k2)=gamma* e*(gamma-1d0)
        gamm(i,i2,j2,k2)=gamma* e
        etot=e
        Bmag(i,i2,j2,k2,1)=A
        Bmag(i,i2,j2,k2,2)=B
        Bmag(i,i2,j2,k2,3)=C
#if NENER>0
        do irad=1,nener
           pres(i,i2,j2,k2)=pres(i,i2,j2,k2)+uloc(i,i2,j2,k2,8+irad)*(gamma_rad(irad)-1d0)
           csou(i,i2,j2,k2)=csou(i,i2,j2,k2)+gamma_rad(irad)*(gamma_rad(irad)-1d0)*uloc(i,i2,j2,k2,8+irad)*oneoverd
           geff(i,i2,j2,k2)=geff(i,i2,j2,k2)+gamma_rad(irad) *uloc(i,i2,j2,k2,8+irad)*(gamma_rad(irad)-1d0)
           gamm(i,i2,j2,k2)=gamm(i,i2,j2,k2)+gamma_rad(irad) *uloc(i,i2,j2,k2,8+irad)
           etot=etot+uloc(i,i2,j2,k2,8+irad)
           ecr(i,i2,j2,k2,irad)=uloc(i,i2,j2,k2,8+irad)
        end do
#endif
#ifdef CRFLX
        do icr=1,ncr
           pres(i,i2,j2,k2)=pres(i,i2,j2,k2)+uloc(i,i2,j2,k2,ind_cr+irad)*(gamma_cr(icr)-1d0)
           csou(i,i2,j2,k2)=csou(i,i2,j2,k2)+gamma_cr(icr)*(gamma_cr(icr)-1d0)*uloc(i,i2,j2,k2,ind_cr+irad)*oneoverd
           geff(i,i2,j2,k2)=geff(i,i2,j2,k2)+gamma_cr(icr) *uloc(i,i2,j2,k2,ind_cr+irad)*(gamma_cr(icr)-1d0)
           gamm(i,i2,j2,k2)=gamm(i,i2,j2,k2)+gamma_cr(icr) *uloc(i,i2,j2,k2,ind_cr+irad)
           etot=etot+uloc(i,i2,j2,k2,ind_cr+irad)
           ecr(i,i2,j2,k2,irad)=uloc(i,i2,j2,k2,ind_cr+irad)
        end do
#endif
        ethe(i,i2,j2,k2)=etot ! YD WARNING: testing with total energy, should be called eint
        geff(i,i2,j2,k2)=geff(i,i2,j2,k2)/pres(i,i2,j2,k2)
        gamm(i,i2,j2,k2)=gamm(i,i2,j2,k2)/etot
        csou(i,i2,j2,k2)=sqrt(csou(i,i2,j2,k2))
        dens(i,i2,j2,k2)=d
     enddo
  enddo
  enddo
  enddo   

  
  ! =====================================================================
  ! Store the shock cell position, its cell index, 
  ! the shock normal direction, the pre-shock and post-shock pressures.
  ! =====================================================================  
     do i=1,ncache
  do k2=k2min,k2max
  do j2=j2min,j2max
  do i2=i2min,i2max
        if(okshock(i,i2,j2,k2))then
        normgradT=0.0d0
        do idim=1,ndim
           normgradT=normgradT+gradT(i,i2,j2,k2,idim)**2
        enddo
        normgradT=sqrt(normgradT)
        do idim=1,ndim
           gradT(i,i2,j2,k2,idim)=gradT(i,i2,j2,k2,idim)/normgradT
        enddo

        ishock=ishock+1

        shock_dir(ishock,1:ndim)=gradT(i,i2,j2,k2,1:ndim)
        ind_shock(ishock)=ind_shock2(i)

        x_shock(ishock,1)=xg(ind_grid(i),1)+(dble(i2)-1.5d0)*dx-skip_loc(1)
#if NDIM>1
        x_shock(ishock,2)=xg(ind_grid(i),2)+(dble(j2)-1.5d0)*dx-skip_loc(2)
#endif
#if NDIM>2
        x_shock(ishock,3)=xg(ind_grid(i),3)+(dble(k2)-1.5d0)*dx-skip_loc(3)
#endif

        ! Get cell index to store the shock_dir for post-processing use (Remove this after tests)
        ! acting on uold because set_uold have been called just before
        if(ndim==1)ind_son=1+(i2-1)
        if(ndim==2)ind_son=1+(i2-1)+2*(j2-1)
        if(ndim==3)ind_son=1+(i2-1)+2*(j2-1)+4*(k2-1)
        iskip=ncoarse+(ind_son-1)*ngridmax
        ind_cell(i)=iskip+ind_grid(i)
!!$        uold(ind_cell(i),imach+1)=shock_dir(ishock,1)*MAX(uold(ind_cell(i),1),smallr)
!!$        uold(ind_cell(i),imach+2)=shock_dir(ishock,2)*MAX(uold(ind_cell(i),1),smallr)

        ! Pre-shock 1 dx
        position(1:ndim)=-gradT(i,i2,j2,k2,1:ndim)
        quantity(iu1:iu2,ju1:ju2,ku1:ku2)=dens(i,iu1:iu2,ju1:ju2,ku1:ku2)           
        dens_pre=cic_interpol(position,quantity,i2,j2,k2)
        quantity(iu1:iu2,ju1:ju2,ku1:ku2)=pres(i,iu1:iu2,ju1:ju2,ku1:ku2)           
        pres_pre=cic_interpol(position,quantity,i2,j2,k2)
        quantity(iu1:iu2,ju1:ju2,ku1:ku2)=ethe(i,iu1:iu2,ju1:ju2,ku1:ku2)           
        ethe_pre=cic_interpol(position,quantity,i2,j2,k2)
        do icr=1,ncr
           quantity(iu1:iu2,ju1:ju2,ku1:ku2)=ecr(i,iu1:iu2,ju1:ju2,ku1:ku2,irad)
           ecr_pre(irad)=cic_interpol(position,quantity,i2,j2,k2)
        enddo
        quantity(iu1:iu2,ju1:ju2,ku1:ku2)=geff(i,iu1:iu2,ju1:ju2,ku1:ku2)           
        geff_pre=cic_interpol(position,quantity,i2,j2,k2)
        quantity(iu1:iu2,ju1:ju2,ku1:ku2)=gamm(i,iu1:iu2,ju1:ju2,ku1:ku2)           
        gamm_pre=cic_interpol(position,quantity,i2,j2,k2)
        quantity(iu1:iu2,ju1:ju2,ku1:ku2)=csou(i,iu1:iu2,ju1:ju2,ku1:ku2)           
        csou_pre=cic_interpol(position,quantity,i2,j2,k2)
        do idim=1,3
           quantity(iu1:iu2,ju1:ju2,ku1:ku2)=Bmag(i,iu1:iu2,ju1:ju2,ku1:ku2,idim)
           Bmag_pre(idim)=cic_interpol(position,quantity,i2,j2,k2)
        enddo
        ! Post-shock 1 dx
        position(1:ndim)= gradT(i,i2,j2,k2,1:ndim)
        quantity(iu1:iu2,ju1:ju2,ku1:ku2)=dens(i,iu1:iu2,ju1:ju2,ku1:ku2)           
        dens_pos=cic_interpol(position,quantity,i2,j2,k2)
        quantity(iu1:iu2,ju1:ju2,ku1:ku2)=pres(i,iu1:iu2,ju1:ju2,ku1:ku2)           
        pres_pos=cic_interpol(position,quantity,i2,j2,k2)
        quantity(iu1:iu2,ju1:ju2,ku1:ku2)=ethe(i,iu1:iu2,ju1:ju2,ku1:ku2)           
        ethe_pos=cic_interpol(position,quantity,i2,j2,k2)
        quantity(iu1:iu2,ju1:ju2,ku1:ku2)=geff(i,iu1:iu2,ju1:ju2,ku1:ku2)           
        geff_pos=cic_interpol(position,quantity,i2,j2,k2)
        quantity(iu1:iu2,ju1:ju2,ku1:ku2)=gamm(i,iu1:iu2,ju1:ju2,ku1:ku2)           
        gamm_pos=cic_interpol(position,quantity,i2,j2,k2)
        ! Pre-shock 2 dx
        position(1:ndim)=-gradT(i,i2,j2,k2,1:ndim)*nearlytwo
        quantity(iu1:iu2,ju1:ju2,ku1:ku2)=dens(i,iu1:iu2,ju1:ju2,ku1:ku2)           
        dens_pre_far=cic_interpol(position,quantity,i2,j2,k2)
        quantity(iu1:iu2,ju1:ju2,ku1:ku2)=pres(i,iu1:iu2,ju1:ju2,ku1:ku2)           
        pres_pre_far=cic_interpol(position,quantity,i2,j2,k2)
        quantity(iu1:iu2,ju1:ju2,ku1:ku2)=ethe(i,iu1:iu2,ju1:ju2,ku1:ku2)           
        ethe_pre_far=cic_interpol(position,quantity,i2,j2,k2)
        do icr=1,ncr
           quantity(iu1:iu2,ju1:ju2,ku1:ku2)=ecr(i,iu1:iu2,ju1:ju2,ku1:ku2,irad)
           ecr_pre_far(irad)=cic_interpol(position,quantity,i2,j2,k2)
        enddo
        quantity(iu1:iu2,ju1:ju2,ku1:ku2)=geff(i,iu1:iu2,ju1:ju2,ku1:ku2)           
        geff_pre_far=cic_interpol(position,quantity,i2,j2,k2)
        quantity(iu1:iu2,ju1:ju2,ku1:ku2)=gamm(i,iu1:iu2,ju1:ju2,ku1:ku2)           
        gamm_pre_far=cic_interpol(position,quantity,i2,j2,k2)
        quantity(iu1:iu2,ju1:ju2,ku1:ku2)=csou(i,iu1:iu2,ju1:ju2,ku1:ku2)           
        csou_pre_far=cic_interpol(position,quantity,i2,j2,k2)
        do idim=1,3
           quantity(iu1:iu2,ju1:ju2,ku1:ku2)=Bmag(i,iu1:iu2,ju1:ju2,ku1:ku2,idim)
           Bmag_pre_far(idim)=cic_interpol(position,quantity,i2,j2,k2)
        enddo
        ! Post-shock 2 dx
        position(1:ndim)= gradT(i,i2,j2,k2,1:ndim)*nearlytwo
        quantity(iu1:iu2,ju1:ju2,ku1:ku2)=dens(i,iu1:iu2,ju1:ju2,ku1:ku2)           
        dens_pos_far=cic_interpol(position,quantity,i2,j2,k2)
        quantity(iu1:iu2,ju1:ju2,ku1:ku2)=pres(i,iu1:iu2,ju1:ju2,ku1:ku2)           
        pres_pos_far=cic_interpol(position,quantity,i2,j2,k2)
        quantity(iu1:iu2,ju1:ju2,ku1:ku2)=ethe(i,iu1:iu2,ju1:ju2,ku1:ku2)           
        ethe_pos_far=cic_interpol(position,quantity,i2,j2,k2)
        quantity(iu1:iu2,ju1:ju2,ku1:ku2)=geff(i,iu1:iu2,ju1:ju2,ku1:ku2)           
        geff_pos_far=cic_interpol(position,quantity,i2,j2,k2)
        quantity(iu1:iu2,ju1:ju2,ku1:ku2)=gamm(i,iu1:iu2,ju1:ju2,ku1:ku2)           
        gamm_pos_far=cic_interpol(position,quantity,i2,j2,k2)


        if (pres_pre .lt. pres_pre_far)then
           dpre(ishock)=dens_pre
           ppre(ishock)=pres_pre
           epre(ishock)=ethe_pre
           ecpre(ishock,1:ncr)=ecr_pre(1:ncr)
           gpre(ishock)=geff_pre
           gapre(ishock)=gamm_pre
           cpre(ishock)=csou_pre
           Bpre(ishock,1:3)=Bmag_pre(1:3)
        else
           dpre(ishock)=dens_pre_far
           ppre(ishock)=pres_pre_far
           epre(ishock)=ethe_pre_far
           ecpre(ishock,1:ncr)=ecr_pre_far(1:ncr)
           gpre(ishock)=geff_pre_far
           gapre(ishock)=gamm_pre_far
           cpre(ishock)=csou_pre_far
           Bpre(ishock,1:3)=Bmag_pre_far(1:3)
        endif
        if (pres_pos .gt. pres_pos_far)then
           dpos(ishock)=dens_pos
           ppos(ishock)=pres_pos
           epos(ishock)=ethe_pos
           gpos(ishock)=geff_pos
           gapos(ishock)=gamm_pos
        else
           dpos(ishock)=dens_pos_far
           ppos(ishock)=pres_pos_far
           epos(ishock)=ethe_pos_far
           gpos(ishock)=geff_pos_far
           gapos(ishock)=gamm_pos_far
        endif
        depos(ishock)=ethe_pos_far-ethe_pos

     endif
  enddo
  enddo
  enddo
     enddo

111 format('   Entering get_shock_unit_vector for level ',i2)

end subroutine get_shock_unit_vector
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine get_shock_unit_vector2(ind_grid,ncache,ind_shock2,nshock,ishock,x_shock,shock_dir,ilevel)

  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel,ncache,nshock,ishock
  integer ,dimension(1:nvector                )::ind_grid,ind_shock2
  real(dp),dimension(1:nshock,1:ndim)::x_shock
  real(dp),dimension(1:nshock,1:ndim)::shock_dir
  !-------------------------------------------------------------------
  !-------------------------------------------------------------------
  integer ,dimension(1:nvector,1:threetondim     ),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim       ),save::nbors_father_grids
  integer ,dimension(1:nvector,0:twondim         ),save::ibuffer_father
  integer ,dimension(1:nvector,0:twondim         ),save::ind1

  real(dp),dimension(1:nvector,0:twondim  ,1:nvar+3),save::u1
  real(dp),dimension(1:nvector,1:twotondim,1:nvar+3),save::u2

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3),save::uloc
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),save::ddx
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndim),save::gradT
  real(dp),dimension(1:3)::skip_loc,xcell

  logical ,dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),save::okshock

  integer,dimension(1:nvector),save::igrid_nbor,ind_cell,ind_buffer,ind_exist,ind_nexist

  integer::ind_buffer1,ind_buffer2,ind_buffer3
  integer::ind_father1,ind_father2,ind_father3
  integer::i,j,ivar,idim,ind_son,ind_father,iskip,nbuffer,ibuffer
  integer::i0,j0,k0,i1,j1,k1,i2,j2,k2,i3,j3,k3,nx_loc,nb_noneigh,nexist
  integer::i1min,i1max,j1min,j1max,k1min,k1max
  integer::i2min,i2max,j2min,j2max,k2min,k2max
  integer::i3min,i3max,j3min,j3max,k3min,k3max
  integer::icr,irad,ii,jj,kk,ig,id,jg,jd,kg,kd,icell,ind_cr

  real(dp)::dx,dx_loc,scale,oneoverd,oneoverdeltax
  real(dp)::d,u,v,w,A,B,C,e
  real(dp)::tempg,tempd,normgradT,delta

  ind_cr=8
  if(twotemp)ind_cr=9
  if(mom_cr_advect)ind_cr=icru-1

  ! Mesh spacing in that level
  nx_loc=icoarse_max-icoarse_min+1
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx=0.5D0**ilevel
  dx_loc=dx*scale

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

  !------------------------------------------
  ! Gather 3^ndim neighboring father cells
  !------------------------------------------
  do i=1,ncache
     ind_cell(i)=father(ind_grid(i))
  end do
  call get3cubefather(ind_cell,nbors_father_cells,nbors_father_grids,ncache,ilevel)
  
  !---------------------------
  ! Gather 6x6x6 cells stencil
  !---------------------------
  ! Loop over 3x3x3 neighboring father cells
  do k1=k1min,k1max
  do j1=j1min,j1max
  do i1=i1min,i1max
     
     ! Check if neighboring grid exists
     nbuffer=0
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
     
     ! If not, interpolate variables from parent cells
     if(nbuffer>0)then
        call getnborfather(ind_buffer,ibuffer_father,nbuffer,ilevel)
        do j=0,twondim
           do ivar=1,nvar+3
              do i=1,nbuffer
                 u1(i,j,ivar)=uold(ibuffer_father(i,j),ivar)
              end do
           end do
           do i=1,nbuffer
              ind1(i,j)=son(ibuffer_father(i,j))
           end do
        end do
        call interpol_hydro(u1,ind1,u2,nbuffer)
     endif

     ! Loop over 2x2x2 cells
     do k2=k2min,k2max
     do j2=j2min,j2max
     do i2=i2min,i2max

        ind_son=1+i2+2*j2+4*k2
        iskip=ncoarse+(ind_son-1)*ngridmax
        do i=1,nexist
           ind_cell(i)=iskip+igrid_nbor(ind_exist(i))
        end do
        
        i3=1; j3=1; k3=1
        if(ndim>0)i3=1+2*(i1-1)+i2
        if(ndim>1)j3=1+2*(j1-1)+j2
        if(ndim>2)k3=1+2*(k1-1)+k2
        
        ! Gather hydro variables
        do ivar=1,nvar+3
           do i=1,nexist
              uloc(ind_exist(i),i3,j3,k3,ivar)=uold(ind_cell(i),ivar)
           end do
           do i=1,nbuffer
              uloc(ind_nexist(i),i3,j3,k3,ivar)=u2(i,ind_son,ivar)
           end do
        end do

        ! Gather divu
        do i=1,nexist
           ddx(ind_exist(i),i3,j3,k3)=dx_loc
        end do
        do i=1,nbuffer
           if(interpol_type.gt.0)then
              ddx(ind_nexist(i),i3,j3,k3)=dx_loc
           else
              ddx(ind_nexist(i),i3,j3,k3)=1.5d0*dx_loc
           endif
        end do
        
     end do
     end do
     end do
     ! End loop over cells

  end do
  end do
  end do
  ! End loop over neighboring grids


  ! Warning: replace ijk2minmax=[0,1] by [1,2]
  ! Integer constants
  i2min=1; i2max=1
  j2min=1; j2max=1
  k2min=1; k2max=1
  if(ndim>0)then
     i1min=1;i2max=2
  end if
  if(ndim>1)then
     j2min=1;j2max=2
  end if
  if(ndim>2)then
     k2min=1;k2max=2
  end if

  okshock=.false.
  do k2=k2min,k2max
  do j2=j2min,j2max
  do i2=i2min,i2max
     do i=1,ncache
        if(ndim==1)ind_son=1+(i2-1)
        if(ndim==2)ind_son=1+(i2-1)+2*(j2-1)
        if(ndim==3)ind_son=1+(i2-1)+2*(j2-1)+4*(k2-1)
        iskip=ncoarse+(ind_son-1)*ngridmax
        icell=iskip+ind_grid(i)
        if(icell==ind_shock2(i))okshock(i,i2,j2,k2)=.true.
     enddo
  enddo
  enddo
  enddo

  ! =====================================================================
  ! Compute gradT and gradS
  ! =====================================================================
  do idim=1,ndim
  if(idim==1)then
     ii=1;jj=0;kk=0 ! jg=jd, kg=kd
  endif
  if(idim==2)then
     ii=0;jj=1;kk=0 ! ig=id, kg=kd
  endif
  if(idim==3)then
     ii=0;jj=0;kk=1 ! ig=id, jg=jd
  endif
  do k2=k2min,k2max
  kg=k2-kk;kd=k2+kk
  do j2=j2min,j2max
  jg=j2-jj;jd=j2+jj
  do i2=i2min,i2max
  ig=i2-ii;id=i2+ii
     do i=1,ncache
        if(okshock(i,i2,j2,k2))then
        ! Right; Top; Front
        d=max(uloc(i,id,jd,kd,1),smallr)
        oneoverd=1d0/d
        u=uloc(i,id,jd,kd,2)*oneoverd
        v=uloc(i,id,jd,kd,3)*oneoverd
        w=uloc(i,id,jd,kd,4)*oneoverd
        A=0.5*(uloc(i,id,jd,kd,6)+uloc(i,id,jd,kd,nvar+1))
        B=0.5*(uloc(i,id,jd,kd,7)+uloc(i,id,jd,kd,nvar+2))
        C=0.5*(uloc(i,id,jd,kd,8)+uloc(i,id,jd,kd,nvar+3))
        e=uloc(i,id,jd,kd,5)-0.5*d*(u**2+v**2+w**2)-0.5*(A**2+B**2+C**2)
#if NENER>0
        do irad=1,nener
           e=e-uloc(i,id,jd,kd,8+irad)
        end do
#endif
        tempd=e*(gamma-1d0)*oneoverd
#if NENER>0
        do irad=1,nener
           tempd=tempd+uloc(i,id,jd,kd,8+irad)*(gamma_rad(irad)-1d0)*oneoverd
        end do
#endif
#ifdef CRFLX
        do icr=1,ncr
           tempd=tempd+uloc(i,id,jd,kd,ind_cr+irad)*(gamma_cr(icr)-1d0)*oneoverd
        end do
#endif

        ! Left; Bottom; Back
        d=max(uloc(i,ig,jg,kg,1),smallr)
        oneoverd=1d0/d
        u=uloc(i,ig,jg,kg,2)*oneoverd
        v=uloc(i,ig,jg,kg,3)*oneoverd
        w=uloc(i,ig,jg,kg,4)*oneoverd
        A=0.5*(uloc(i,ig,jg,kg,6)+uloc(i,ig,jg,kg,nvar+1))
        B=0.5*(uloc(i,ig,jg,kg,7)+uloc(i,ig,jg,kg,nvar+2))
        C=0.5*(uloc(i,ig,jg,kg,8)+uloc(i,ig,jg,kg,nvar+3))
        e=uloc(i,ig,jg,kg,5)-0.5*d*(u**2+v**2+w**2)-0.5*(A**2+B**2+C**2)
#if NENER>0
        do irad=1,nener
           e=e-uloc(i,ig,jg,kg,8+irad)
        end do
#endif
        tempg=e*(gamma-1d0)*oneoverd
#if NENER>0
        do irad=1,nener
           tempg=tempg+uloc(i,ig,jg,kg,8+irad)*(gamma_rad(irad)-1d0)*oneoverd
        end do
#endif
#ifdef CRFLX
        do icr=1,ncr
           tempg=tempg+uloc(i,ig,jg,kg,ind_cr+irad)*(gamma_cr(icr)-1d0)*oneoverd
        end do
#endif

        oneoverdeltax=1d0/(ddx(i,id,jd,kd)+ddx(i,ig,jg,kg))
        gradT(i,i2,j2,k2,idim)=(tempd-tempg)*oneoverdeltax

     endif
     enddo
  enddo
  enddo
  enddo
  enddo   

  ! =====================================================================
  ! Store the shock cell position, its cell index, 
  ! the shock normal direction, the pre-shock and post-shock pressures.
  ! =====================================================================
  ! WARNING: weird ordering of loops (check that...)
     do i=1,ncache
  do k2=k2min,k2max
  do j2=j2min,j2max
  do i2=i2min,i2max
     if(okshock(i,i2,j2,k2))then
        normgradT=0.0d0
        do idim=1,ndim
           normgradT=normgradT+gradT(i,i2,j2,k2,idim)**2
        enddo
        normgradT=sqrt(normgradT)
        do idim=1,ndim
           gradT(i,i2,j2,k2,idim)=gradT(i,i2,j2,k2,idim)/normgradT
        enddo

        ishock=ishock+1

        shock_dir(ishock,1:ndim)=gradT(i,i2,j2,k2,1:ndim)

        x_shock(ishock,1)=xg(ind_grid(i),1)+(dble(i2)-1.5d0)*dx-skip_loc(1)
#if NDIM>1
        x_shock(ishock,2)=xg(ind_grid(i),2)+(dble(j2)-1.5d0)*dx-skip_loc(2)
#endif
#if NDIM>2
        x_shock(ishock,3)=xg(ind_grid(i),3)+(dble(k2)-1.5d0)*dx-skip_loc(3)
#endif

     endif
  enddo
  enddo
  enddo
     enddo

111 format('   Entering get_shock_unit_vector2 for level ',i2)

end subroutine get_shock_unit_vector2
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine get_mycell_index(cell_index,cell_index_father,grid_index,cell_lev,xpart,ilevel,n)
  use amr_commons
  implicit none

  integer::n,ilevel
  integer ,dimension(1:nvector)::cell_index,cell_index_father
  integer ,dimension(1:nvector)::cell_lev
  integer ,dimension(1:nvector)::grid_index
  real(dp),dimension(1:nvector,1:ndim)::xpart

  !----------------------------------------------------------------------------
  ! This routine returns the index and level of the cell, (at maximum level
  ! ilevel), in which the input the position specified by xpart lies
  !----------------------------------------------------------------------------

  real(dp)::xx,yy,zz
  integer::i,j,ii,jj,kk,ind,iskip,igrid,ind_cell,ind_cell_father,igrid0

#if NDIM==1
  igrid0=son(1+icoarse_min)
#endif
#if NDIM==2
  igrid0=son(1+icoarse_min+jcoarse_min*nx)
#endif
#if NDIM==3
  igrid0=son(1+icoarse_min+jcoarse_min*nx+kcoarse_min*nx*ny)
#endif

  ii=0;jj=0;kk=0
  do i=1,n
     xx = xpart(i,1)/boxlen + (nx-1)/2.0
#if NDIM>1
     yy = xpart(i,2)/boxlen + (ny-1)/2.0
#endif
#if NDIM>2
     zz = xpart(i,3)/boxlen + (nz-1)/2.0
#endif

     if(xx<0.)xx=xx+dble(nx)
     if(xx>dble(nx))xx=xx-dble(nx)
#if NDIM>1
     if(yy<0.)yy=yy+dble(ny)
     if(yy>dble(ny))yy=yy-dble(ny)
#endif
#if NDIM>2
     if(zz<0.)zz=zz+dble(nz)
     if(zz>dble(nz))zz=zz-dble(nz)
#endif

     ind_cell=0
     ind_cell_father=0

     igrid=igrid0
     do j=1,ilevel 
        ii=1
        if(xx<xg(igrid,1))ii=0
#if NDIM>1
        jj=1
        if(yy<xg(igrid,2))jj=0
#endif
#if NDIM>2
        kk=1
        if(zz<xg(igrid,3))kk=0
#endif
        ind=1+ii+2*jj+4*kk
        iskip=ncoarse+(ind-1)*ngridmax
        ind_cell=iskip+igrid
        grid_index(i)=igrid
        igrid=son(ind_cell)
        cell_lev(i)=j
        if(igrid==0.or.j==ilevel)exit
        ind_cell_father=ind_cell
     end do
     cell_index(i)=ind_cell
     cell_index_father(i)=ind_cell_father
  end do
end subroutine get_mycell_index
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine get_mycell_index2(cell_index,xpart,xgrid,ilevel,n)
  use amr_commons
  implicit none

  integer::n,ilevel
  integer ,dimension(1:nvector)::cell_index
  real(dp),dimension(1:nvector,1:ndim)::xpart
  real(dp),dimension(1:nvector,1:ndim)::xgrid

  !----------------------------------------------------------------------------
  ! This routine returns the index and level of the cell, (at maximum level
  ! ilevel), in which the input the position specified by xpart lies.
  ! Difference with get_mycell_index is that it outputs the grid position
  !----------------------------------------------------------------------------

  real(dp)::xx,yy,zz
  integer::i,j,ii,jj,kk,ind,iskip,igrid,igrid_prev,ind_cell,igrid0

  ind_cell=0
#if NDIM==1
  igrid0=son(1+icoarse_min)
#endif
#if NDIM==2
  igrid0=son(1+icoarse_min+jcoarse_min*nx)
#endif
#if NDIM==3
  igrid0=son(1+icoarse_min+jcoarse_min*nx+kcoarse_min*nx*ny)
#endif

  ii=0;jj=0;kk=0
  do i=1,n
     xx = xpart(i,1)/boxlen + (nx-1)/2.0
#if NDIM>1
     yy = xpart(i,2)/boxlen + (ny-1)/2.0
#endif
#if NDIM>2
     zz = xpart(i,3)/boxlen + (nz-1)/2.0
#endif

     if(xx<0.)xx=xx+dble(nx)
     if(xx>dble(nx))xx=xx-dble(nx)
#if NDIM>1
     if(yy<0.)yy=yy+dble(ny)
     if(yy>dble(ny))yy=yy-dble(ny)
#endif
#if NDIM>2
     if(zz<0.)zz=zz+dble(nz)
     if(zz>dble(nz))zz=zz-dble(nz)
#endif

     igrid=igrid0
     do j=1,ilevel 
        ii=1
        if(xx<xg(igrid,1))ii=0
#if NDIM>1
        jj=1
        if(yy<xg(igrid,2))jj=0
#endif
#if NDIM>2
        kk=1
        if(zz<xg(igrid,3))kk=0
#endif
        ind=1+ii+2*jj+4*kk
        iskip=ncoarse+(ind-1)*ngridmax
        ind_cell=iskip+igrid
        igrid_prev=igrid
        igrid=son(ind_cell)
        if(igrid==0.or.j==ilevel)exit
     end do
     cell_index(i)=ind_cell
     xgrid(i,1:ndim)=xg(igrid_prev,1:ndim)
  end do
end subroutine get_mycell_index2
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
function cic_interpol(position,quantity,i2,j2,k2)
  use amr_commons
  use hydro_commons
  implicit none
  real(dp),dimension(1:ndim)::position
  real(dp),dimension(1:ndim)::x,dd,dg
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2)::quantity
  real(dp),dimension(1:twotondim)::vol

  integer,dimension(1:ndim)::itarget
  integer,dimension(1:ndim)::iig,iid
  integer,dimension(1:twotondim)::icic,jcic,kcic

  integer::idim,icell,i2,j2,k2

  real(dp)::cic_interpol

  cic_interpol=0.0d0

  itarget(1:ndim)=nint(position(1:ndim))
  x(1:ndim)=position(1:ndim)-dble(itarget(1:ndim))+0.5D0 ! Rescale position to the cell interval [0,1]
  do idim=1,ndim
     dd(idim)=x(idim)+0.5D0
     iid(idim)=dd(idim)
     dd(idim)=dd(idim)-iid(idim)
     dg(idim)=1.0D0-dd(idim)
     iig(idim)=iid(idim)-1
  enddo
#if NDIM==1
  icic(1)=iig(1)+itarget(1)
  icic(2)=iid(1)+itarget(1)
  vol(1)=dg(1)
  vol(2)=dd(1)
  do icell=1,twotondim
     cic_interpol=cic_interpol+vol(icell)*quantity(i2+icic(icell),j2,k2)
  enddo
#endif
#if NDIM==2
  icic(1)=iig(1)+itarget(1)
  icic(2)=iid(1)+itarget(1)
  icic(3)=iig(1)+itarget(1)
  icic(4)=iid(1)+itarget(1)
  jcic(1)=iig(2)+itarget(2)
  jcic(2)=iig(2)+itarget(2)
  jcic(3)=iid(2)+itarget(2)
  jcic(4)=iid(2)+itarget(2)
  vol(1)=dg(1)*dg(2)
  vol(2)=dd(1)*dg(2)
  vol(3)=dg(1)*dd(2)
  vol(4)=dd(1)*dd(2)
  do icell=1,twotondim
     cic_interpol=cic_interpol+vol(icell)*quantity(i2+icic(icell),j2+jcic(icell),k2)
  enddo
#endif
#if NDIM==3
  icic(1)=iig(1)+itarget(1)
  icic(2)=iid(1)+itarget(1)
  icic(3)=iig(1)+itarget(1)
  icic(4)=iid(1)+itarget(1)
  icic(5)=iig(1)+itarget(1)
  icic(6)=iid(1)+itarget(1)
  icic(7)=iig(1)+itarget(1)
  icic(8)=iid(1)+itarget(1)
  jcic(1)=iig(2)+itarget(2)
  jcic(2)=iig(2)+itarget(2)
  jcic(3)=iid(2)+itarget(2)
  jcic(4)=iid(2)+itarget(2)
  jcic(5)=iig(2)+itarget(2)
  jcic(6)=iig(2)+itarget(2)
  jcic(7)=iid(2)+itarget(2)
  jcic(8)=iid(2)+itarget(2)
  kcic(1)=iig(3)+itarget(3)
  kcic(2)=iig(3)+itarget(3)
  kcic(3)=iig(3)+itarget(3)
  kcic(4)=iig(3)+itarget(3)
  kcic(5)=iid(3)+itarget(3)
  kcic(6)=iid(3)+itarget(3)
  kcic(7)=iid(3)+itarget(3)
  kcic(8)=iid(3)+itarget(3)
  vol(1)=dg(1)*dg(2)*dg(3)
  vol(2)=dd(1)*dg(2)*dg(3)
  vol(3)=dg(1)*dd(2)*dg(3)
  vol(4)=dd(1)*dd(2)*dg(3)
  vol(5)=dg(1)*dg(2)*dd(3)
  vol(6)=dd(1)*dg(2)*dd(3)
  vol(7)=dg(1)*dd(2)*dd(3)
  vol(8)=dd(1)*dd(2)*dd(3)
  do icell=1,twotondim
     cic_interpol=cic_interpol+vol(icell)*quantity(i2+icic(icell),j2+jcic(icell),k2+kcic(icell))
  enddo
#endif

end function cic_interpol
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
function cic_value(quantity,icic,jcic,kcic,vol)
  use amr_commons  ,only:twotondim,dp
  use hydro_commons,only:iu1,iu2,ju1,ju2,ku1,ku2
  implicit none
  integer,dimension(1:twotondim)::icic,jcic,kcic
  real(dp),dimension(iu1:iu2,ju1:ju2,ku1:ku2)::quantity
  real(dp),dimension(1:twotondim)::vol
  integer::icell

  real(dp)::cic_value

  cic_value=0.0d0
  do icell=1,twotondim
     cic_value=cic_value+vol(icell)*quantity(icic(icell),jcic(icell),kcic(icell))
  enddo

end function cic_value
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine get_cloud_cic(position,icic,jcic,kcic,vol)
  use amr_commons, only:ndim,twotondim,dp
  implicit none
  real(dp),dimension(1:ndim)::position
  real(dp),dimension(1:ndim)::x,dd,dg
  real(dp),dimension(1:twotondim)::vol

  integer,dimension(1:ndim)::itarget
  integer,dimension(1:ndim)::iig,iid
  integer,dimension(1:twotondim)::icic,jcic,kcic

  integer::idim,icell,i2,j2,k2

  icic=0;jcic=0;kcic=0

  itarget(1:ndim)=nint(position(1:ndim))
  x(1:ndim)=position(1:ndim)-dble(itarget(1:ndim))+0.5D0 ! Rescale position to the cell interval [0,1]
  do idim=1,ndim
     dd(idim)=x(idim)+0.5D0
     iid(idim)=dd(idim)
     dd(idim)=dd(idim)-iid(idim)
     dg(idim)=1.0D0-dd(idim)
     iig(idim)=iid(idim)-1
  enddo
#if NDIM==1
  icic(1)=iig(1)+itarget(1)
  icic(2)=iid(1)+itarget(1)
  vol(1)=dg(1)
  vol(2)=dd(1)
#endif
#if NDIM==2
  icic(1)=iig(1)+itarget(1)
  icic(2)=iid(1)+itarget(1)
  icic(3)=iig(1)+itarget(1)
  icic(4)=iid(1)+itarget(1)
  jcic(1)=iig(2)+itarget(2)
  jcic(2)=iig(2)+itarget(2)
  jcic(3)=iid(2)+itarget(2)
  jcic(4)=iid(2)+itarget(2)
  vol(1)=dg(1)*dg(2)
  vol(2)=dd(1)*dg(2)
  vol(3)=dg(1)*dd(2)
  vol(4)=dd(1)*dd(2)
#endif
#if NDIM==3
  icic(1)=iig(1)+itarget(1)
  icic(2)=iid(1)+itarget(1)
  icic(3)=iig(1)+itarget(1)
  icic(4)=iid(1)+itarget(1)
  icic(5)=iig(1)+itarget(1)
  icic(6)=iid(1)+itarget(1)
  icic(7)=iig(1)+itarget(1)
  icic(8)=iid(1)+itarget(1)
  jcic(1)=iig(2)+itarget(2)
  jcic(2)=iig(2)+itarget(2)
  jcic(3)=iid(2)+itarget(2)
  jcic(4)=iid(2)+itarget(2)
  jcic(5)=iig(2)+itarget(2)
  jcic(6)=iig(2)+itarget(2)
  jcic(7)=iid(2)+itarget(2)
  jcic(8)=iid(2)+itarget(2)
  kcic(1)=iig(3)+itarget(3)
  kcic(2)=iig(3)+itarget(3)
  kcic(3)=iig(3)+itarget(3)
  kcic(4)=iig(3)+itarget(3)
  kcic(5)=iid(3)+itarget(3)
  kcic(6)=iid(3)+itarget(3)
  kcic(7)=iid(3)+itarget(3)
  kcic(8)=iid(3)+itarget(3)
  vol(1)=dg(1)*dg(2)*dg(3)
  vol(2)=dd(1)*dg(2)*dg(3)
  vol(3)=dg(1)*dd(2)*dg(3)
  vol(4)=dd(1)*dd(2)*dg(3)
  vol(5)=dg(1)*dg(2)*dd(3)
  vol(6)=dd(1)*dg(2)*dd(3)
  vol(7)=dg(1)*dd(2)*dd(3)
  vol(8)=dd(1)*dd(2)*dd(3)
#endif

end subroutine get_cloud_cic
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine cmp_divvold(divv,icell,ilevel)
  use amr_commons  ,only:dp,nvector,twondim
  use hydro_commons,only:uold,smallr
  implicit none
  !----------------------------------------------------------------------------
  ! This routine computes divv from the uold values
  !----------------------------------------------------------------------------
  integer::icell,ilevel,ncell
  integer ,dimension(1:nvector),save::ind_cell2
  integer ,dimension(1:nvector,0:twondim)::ind_nbor
  real(dp)::d1,d2,d3,d4,d5,d6,divv,ur,ul

  ncell = 1 ! we just want the neighbors of that cell
  ind_cell2(1) = icell
  call getnbor(ind_cell2,ind_nbor,ncell,ilevel)
  
  d1 = MAX(uold(ind_nbor(1,1),1),smallr) ; d2 = MAX(uold(ind_nbor(1,2),1),smallr)
  ul   = uold(ind_nbor(1,2),2)/d2
  ur   = uold(ind_nbor(1,1),2)/d1
  divv = (ur-ul)
#if NDIM>1
  d3 = MAX(uold(ind_nbor(1,3),1),smallr) ; d4 = MAX(uold(ind_nbor(1,4),1),smallr)
  ul   = uold(ind_nbor(1,4),3)/d4
  ur   = uold(ind_nbor(1,3),3)/d3
  divv = divv + (ur-ul)
#endif
#if NDIM>2
  d5 = MAX(uold(ind_nbor(1,5),1),smallr) ; d6 = MAX(uold(ind_nbor(1,6),1),smallr)
  ul   = uold(ind_nbor(1,6),4)/d6
  ur   = uold(ind_nbor(1,5),4)/d5
  divv = divv + (ur-ul)
#endif

end subroutine cmp_divvold
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine cmp_ethermnew(e,icell)
  use amr_commons
  use hydro_commons
  implicit none
  !----------------------------------------------------------------------------
  ! This routine computes etherm from the unew values
  !----------------------------------------------------------------------------
  integer::icell,i
  real(dp)::d,u,v,w,A,B,C,e,oneoverd

  d=MAX(unew(icell,1),smallr)
  oneoverd=1d0/d
  u=unew(icell,2)*oneoverd
  v=unew(icell,3)*oneoverd
  w=unew(icell,4)*oneoverd
  A=0.5*(unew(icell,6)+unew(icell,nvar+1))
  B=0.5*(unew(icell,7)+unew(icell,nvar+2))
  C=0.5*(unew(icell,8)+unew(icell,nvar+3))
  e=unew(icell,5)-0.5*d*(u**2+v**2+w**2)-0.5*(A**2+B**2+C**2)
#if NENER>0
  do i=1,nener
     e=e-unew(icell,8+i)
  end do
#endif

end subroutine cmp_ethermnew
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine cmp_ethermold(e,icell)
  use amr_commons
  use hydro_commons
  implicit none

  integer::icell,i

  !----------------------------------------------------------------------------
  ! This routine computes etherm from the uold values
  !----------------------------------------------------------------------------

  real(dp)::d,u,v,w,A,B,C,e,oneoverd

  d=MAX(uold(icell,1),smallr)
  oneoverd=1d0/d
  u=uold(icell,2)*oneoverd
  v=uold(icell,3)*oneoverd
  w=uold(icell,4)*oneoverd
  A=0.5*(uold(icell,6)+uold(icell,nvar+1))
  B=0.5*(uold(icell,7)+uold(icell,nvar+2))
  C=0.5*(uold(icell,8)+uold(icell,nvar+3))
  e=uold(icell,5)-0.5*d*(u**2+v**2+w**2)-0.5*(A**2+B**2+C**2)
#if NENER>0
  do i=1,nener
     e=e-uold(icell,8+i)
  end do
#endif

end subroutine cmp_ethermold
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine cmp_eintnew(e,icell)
  use amr_commons
  use hydro_commons
  implicit none

  integer::icell,i

  !----------------------------------------------------------------------------
  ! This routine computes eint from the unew values
  !----------------------------------------------------------------------------

  real(dp)::d,u,v,w,A,B,C,e,oneoverd

  d=MAX(unew(icell,1),smallr)
  oneoverd=1d0/d
  u=unew(icell,2)*oneoverd
  v=unew(icell,3)*oneoverd
  w=unew(icell,4)*oneoverd
  A=0.5*(unew(icell,6)+unew(icell,nvar+1))
  B=0.5*(unew(icell,7)+unew(icell,nvar+2))
  C=0.5*(unew(icell,8)+unew(icell,nvar+3))
  e=unew(icell,5)-0.5*d*(u**2+v**2+w**2)-0.5*(A**2+B**2+C**2)

end subroutine cmp_eintnew
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine cmp_eintold(e,icell)
  use amr_commons
  use hydro_commons
  implicit none

  integer::icell,i

  !----------------------------------------------------------------------------
  ! This routine computes eint from the uold values
  !----------------------------------------------------------------------------

  real(dp)::d,u,v,w,A,B,C,e,oneoverd

  d=MAX(uold(icell,1),smallr)
  oneoverd=1d0/d
  u=uold(icell,2)*oneoverd
  v=uold(icell,3)*oneoverd
  w=uold(icell,4)*oneoverd
  A=0.5*(uold(icell,6)+uold(icell,nvar+1))
  B=0.5*(uold(icell,7)+uold(icell,nvar+2))
  C=0.5*(uold(icell,8)+uold(icell,nvar+3))
  e=uold(icell,5)-0.5*d*(u**2+v**2+w**2)-0.5*(A**2+B**2+C**2)

end subroutine cmp_eintold
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine interpol_creff(mach,fcr,creff)
  use amr_commons,only:dp
  use cr_commons
  implicit none

  logical::continuemach,continuefcr
  integer::indmach,indfcr
  real(dp)::mach,fcr,creff
  real(dp)::x1,x2,y1,y2,dlfx,dlfy,dlfxy,delx,dely

  indmach=0
  continuemach=.true.
  do while(continuemach)
     if(mach .lt. mach_sampling(indmach+1))then
        continuemach=.false.
     else
        if(indmach.lt.nbinm)then
           indmach=indmach+1
        else
           continuemach=.false.
        endif
     endif
  enddo
  indfcr=0
  continuefcr=.true.
  do while(continuefcr)
     if(fcr .lt. fcr_sampling(indfcr+1))then
        continuefcr=.false.
     else
        if(indfcr.lt.nbinm)then
           indfcr=indfcr+1
        else
           continuefcr=.false.
        endif
     endif
  enddo
  
  if    (indmach.ge.1 .and. indmach.lt.nbinm .and.&
       & indfcr .ge.1 .and. indfcr .lt.nbinc)then
     ! bilinear interpolation
     x1=mach_sampling(indmach  )
     x2=mach_sampling(indmach+1)
     y1=fcr_sampling (indfcr   )
     y2=fcr_sampling (indfcr +1)
     dlfx = lcref_sampling(indmach+1,indfcr  )-lcref_sampling(indmach,indfcr)
     dlfy = lcref_sampling(indmach  ,indfcr+1)-lcref_sampling(indmach,indfcr)
     dlfxy= lcref_sampling(indmach+1,indfcr+1)+lcref_sampling(indmach,indfcr  ) &
          &-lcref_sampling(indmach+1,indfcr  )-lcref_sampling(indmach,indfcr+1)
     delx=(mach-x1)/(x2-x1)
     dely=(fcr -y1)/(y2-y1)
     creff= 10d0**(lcref_sampling(indmach,indfcr)+delx*dlfx+dely*dlfy+delx*dely*dlfxy)
  else if(indmach.lt.1)then
     creff=0.0d0
  else if(indmach.ge.nbinm)then
     creff=1.0d0
  else if(indfcr .ge.nbinc)then
     creff=1.0d0
  endif
  
end subroutine interpol_creff
