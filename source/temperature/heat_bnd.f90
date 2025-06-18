!**************************************************************
    subroutine heat_bc_init
! Initializes the heat boundary conditions
!**************************************************************
    use heat_def
    use bnd_def
    use diag_lib
    implicit none
    integer :: i,iheat
    
!=== Find all temperature boundaries ================
    !Assumes that temperature is specified at all boundaries
    nheatstr = nQstr + nTHstr + nHstr + nMHstr + nMHVstr + nCSstr   &
             + nNHstr + nNHVstr + nNTHstr + nNTHVstr
    allocate(heat_str(nheatstr))
    heat_str%inc = 1 !Initialize
    iheat = 0
    
!--- River BC -------------------------------------------
    do i=1,nQstr        
      call heat_bnd_chk(Q_str(i)%bidfile,Q_str(i)%bidpath,&
         Q_str(i)%ncells,Q_str(i)%cells,Q_str(i)%faces,iheat)  
    enddo

!--- Tidal/Harmonic BC ---------------------------------
    do i=1,nTHstr
      call heat_bnd_chk(TH_str(i)%bidfile,TH_str(i)%bidpath,&
         TH_str(i)%ncells,TH_str(i)%cells,TH_str(i)%faces,iheat)
    enddo
    
!--- Single Water Level BC --------------------------
    do i=1,nHstr  
      call heat_bnd_chk(H_str(i)%bidfile,H_str(i)%bidpath,&
         H_str(i)%ncells,H_str(i)%cells,H_str(i)%faces,iheat)
    enddo

!--- Multiple Water Level BC ------------------------------
    do i=1,nMHstr
      call heat_bnd_chk(MH_str(i)%bidfile,MH_str(i)%bidpath,&
         MH_str(i)%ncells,MH_str(i)%cells,MH_str(i)%faces,iheat)
    enddo
    
!--- Multiple Water Level and Velocity BC -------------------------
    do i=1,nMHVstr
      call heat_bnd_chk(MHV_str(i)%bidfile,MHV_str(i)%bidpath,&
         MHV_str(i)%ncells,MHV_str(i)%cells,MHV_str(i)%faces,iheat)
    enddo 

!--- Cross-shore BC (Type CS=6) --------------------------
    do i=1,nCSstr 
      call heat_bnd_chk(CS_str(i)%bidfile,CS_str(i)%bidpath,&
         CS_str(i)%ncells,CS_str(i)%cells,CS_str(i)%faces,iheat)  
    enddo
    
!--- Nested Water Level BC (Type 7=NH) ----------------------------
    do i=1,nNHstr    
      call heat_bnd_chk(NH_str(i)%bidfile,NH_str(i)%bidpath,&
         NH_str(i)%ncells,NH_str(i)%cells,NH_str(i)%faces,iheat)  
    enddo        
        
!--- Nested Water Level and Velocity BC (Type 8=NHV) ---------------
    do i=1,nNHVstr    
      call heat_bnd_chk(NHV_str(i)%bidfile,NHV_str(i)%bidpath,&
         NHV_str(i)%ncells,NHV_str(i)%cells,NHV_str(i)%faces,iheat)  
    enddo 
    
!--- Nested Tidal Database WSE BC (Type 9=NTH) -----------------------
    do i=1,nNTHstr    
      call heat_bnd_chk(NTH_str(i)%bidfile,NTH_str(i)%bidpath,&
         NTH_str(i)%ncells,NTH_str(i)%cells,NTH_str(i)%faces,iheat)  
    enddo 
    
!--- Nested Tidal Database WSE and Velocity BC (Type 10=NTHV) ----------  
    do i=1,nNTHVstr    
      call heat_bnd_chk(NTHV_str(i)%bidfile,NTHV_str(i)%bidpath,&
         NTHV_str(i)%ncells,NTHV_str(i)%cells,NTHV_str(i)%faces,iheat)  
    enddo 
    
    !Check
    if(iheat/=nheatstr)then
      call diag_print_error('Problem determining the temperature boundary conditions')
    endif
    
    return
    end subroutine heat_bc_init
    
!************************************************************************
    subroutine heat_bnd_chk(bidfile,bidpath,nbndcells,icells,kfaces,iheat)
! Checks the flow boundary condition file and paths 
! for temperature information
!************************************************************************
#include "CMS_cpp.h"
    use bnd_def
    use diag_lib
    use heat_def
    use size_def
    use geo_def,   only: cell2cell,ncface
    use comvarbl,  only: mpfile
    use const_def, only: READONLY
#ifdef XMDF_IO
    use xmdf
#endif
    implicit none
    !Input/Output
    character(len=*),intent(in) :: bidfile,bidpath
    integer, intent(in) :: nbndcells
    integer, intent(in) :: icells(nbndcells),kfaces(nbndcells)
    integer, intent(inout) :: iheat
    !Internal Variables
    integer(XID) :: PID,GID
    integer :: nstimes,ierr
    
#ifdef XMDF_IO
    call XF_OPEN_FILE(bidfile,READONLY,PID,ierr)
    call XF_OPEN_GROUP(PID,bidpath,GID,ierr)
    call XF_GET_PROPERTY_NUMBER(GID,'Temp_Times',nstimes,ierr)
    if(ierr<0)then
      call diag_print_error('Heat must be specified at multi-wse boundary conditions')
    endif
    iheat = iheat + 1
    heat_str(iheat)%ntimes = nstimes
    allocate(heat_str(iheat)%val(nstimes))
    allocate(heat_str(iheat)%timeheat(nstimes))
    call XF_READ_PROPERTY_FLOAT(GID,'Temp_Times',nstimes,heat_str(iheat)%timeheat(1),ierr)
    call XF_READ_PROPERTY_FLOAT(GID,'Temperature',nstimes,heat_str(iheat)%val(1),ierr)
    call XF_CLOSE_GROUP(GID,ierr)
    call XF_CLOSE_FILE(PID,ierr)
    
    !Cell id's and faces
    heat_str(iheat)%ncells = nbndcells
    allocate(heat_str(iheat)%cells(nbndcells))
    allocate(heat_str(iheat)%faces(nbndcells))
    heat_str(iheat)%cells = icells
    heat_str(iheat)%faces = kfaces
#endif

    return
    end subroutine heat_bnd_chk
    
!**************************************************************
    subroutine bndheateval
! Evaluate heat at boundary cell strings
!**************************************************************    
    use size_def
    use flow_def
    use geo_def, only: cell2cell
    use bnd_def
    use comvarbl, only: timehrs,ramp
    use heat_def
    use const_def, only: eps
    use prec_def
    implicit none
    integer :: i,ntimes,inc,iriv,iwse,j,k,nck  
    real(ikind) :: fac, heatconst

    do i=1,nheatstr
      if(heat_str(i)%ntimes==0) cycle !Boundary value is constant
      !find out where we are in the time/value arrays
      inc = heat_str(i)%inc
      ntimes = heat_str(i)%ntimes
      do while((timehrs+eps)>=heat_str(i)%timeheat(inc+1) .and. inc<ntimes)
        inc = inc + 1
        heat_str(i)%inc = inc
      enddo
      fac = (timehrs-heat_str(i)%timeheat(inc))/                    &
            (heat_str(i)%timeheat(inc+1)-heat_str(i)%timeheat(inc))
      heat_str(i)%heatbnd = (1.0-fac)*heat_str(i)%val(inc) +        &
            fac*heat_str(i)%val(inc+1)
    enddo
 
    !do iriv=1,nQstr
    !   if(iriv.eq.1) then     !Temporary
    !      heatconst = 15
    !   elseif(iriv.eq.2) then
    !      heatconst = 5
    !   else
    !      heatconst = 10
    !   endif
    !   do j=1,Q_str(iriv)%ncells
    !      i = Q_str(iriv)%cells(j)
    !      k = Q_str(iriv)%faces(j)
    !      nck = cell2cell(k,i)
    !      heat(nck) = heatconst
    !   enddo
    !enddo   

!    do iwse=1,nHstr
!       do j=1,H_str(iwse)%ncells
!          i = H_str(iwse)%cells(j)
!          k = H_str(iwse)%faces(j)
!          nck = cell2cell(k,i)
!          heat(nck) = 20        !Temporary
!       enddo
!    enddo      

!    do iwse=1,nTHstr
!       do j=1,TH_str(iwse)%ncells
!          i = TH_str(iwse)%cells(j)
!          k = TH_str(iwse)%faces(j)
!          nck = cell2cell(k,i)
!          heat(nck) = 20        !Temporary
!       enddo
!    enddo      
    
    return
    end subroutine bndheateval
    
!***********************************************************************
    subroutine bound_heat
! Applies boundary condition to heat transfer equation

!***********************************************************************
    use size_def
    use geo_def
    use bnd_def
    use flow_def
    use heat_def
    use comvarbl, only: ramp
    use prec_def
    implicit none
    integer :: i,j,k,iht,nck,idrym,iriv,iwse

!--- Dry nodes -------------------------------
!$OMP PARALLEL DO PRIVATE(i,k,idrym,nck)
    do i=1,ncells
      if(iwet(i)==1)then
        do k=1,ncface(i)
          acoef(k,i)=acoef(k,i)*iwet(cell2cell(k,i))  
        enddo             
      elseif(iwet(i)==0)then
        idrym=0
        do k=1,ncface(i)
          nck=cell2cell(k,i)
          if(nck<=ncells) idrym=idrym+iwet(nck)
        enddo
!        apheat=0.0
        suheat0(i)=0.0
        do k=1,ncface(i)
          nck=cell2cell(k,i)
          if(idrym>=1)then
            acoef(k,i)=iwet(nck)
          else
            acoef(k,i)=dsxy(k,i)
          endif
          if(nck>ncells) acoef(k,i)=0.0
          !apheat=apheat+acoef(k,i)
          !suheat0(i)=suheat0(i)+acoef(k,i)*heat(nck)
        enddo
        !suheat0(i)=suheat0(i)-apheat*heat(i)
        sp(i)=0.0
      endif
    enddo
 !$OMP END PARALLEL DO
 
!--- heat boundary ---------------------
  !  do iht=1,nheatstr
  !    do j=1,heat_str(iht)%ncells
  !      i=heat_str(iht)%cells(j)
  !      !if(i==0) cycle !Land cell   
  !       k=heat_str(iht)%faces(j)
  !      !do k=1,ncface(i)
  !        nck=cell2cell(k,i)
  !        !if(nck>ncells)then
  !          if(flux(k,i)<0.0)then !Inflow
  !             heat(nck)=ramp*heat_str(iht)%heatbnd+(1-ramp)*heat_str(iht)%heatbnd0(j)
  !             suheat0(i)=suheat0(i)+acoef(k,i)*heat(nck)
  !             sp(i)=sp(i)-acoef(k,i)
  !          else  !Outflow
  !            heat(nck)=heat(i)
  !          endif
  !          acoef(k,i)=0.0
  !        !endif
  !      !enddo
  !    enddo
  !  enddo   

     do iriv=1,nQstr
       do j=1,Q_str(iriv)%ncells
         i = Q_str(iriv)%cells(j)
         k = Q_str(iriv)%faces(j)
         nck = cell2cell(k,i)
         if(flux(k,i)<0.0)then !Inflow
           suheat0(i)=suheat0(i)+acoef(k,i)*heat(nck)
           sp(i)=sp(i)-acoef(k,i)
         else  !Outflow
           heat(nck)=heat(i)
         endif
         acoef(k,i)=0.0
       enddo    
     enddo    

     do iwse=1,nHstr
       do j=1,H_str(iwse)%ncells
         i = H_str(iwse)%cells(j)
         k = H_str(iwse)%faces(j)
         nck = cell2cell(k,i)
         if(flux(k,i)<0.0)then !Inflow
           suheat0(i)=suheat0(i)+acoef(k,i)*heat(nck)
           sp(i)=sp(i)-acoef(k,i)
         else  !Outflow
           heat(nck)=heat(i)
         endif
         acoef(k,i)=0.0
       enddo
     enddo      

     do iwse=1,nTHstr
       do j=1,TH_str(iwse)%ncells
         i = TH_str(iwse)%cells(j)
         k = TH_str(iwse)%faces(j)
         nck = cell2cell(k,i)
         if(flux(k,i)<0.0)then !Inflow
           suheat0(i)=suheat0(i)+acoef(k,i)*heat(nck)
           sp(i)=sp(i)-acoef(k,i)
         else  !Outflow
           heat(nck)=heat(i)
         endif
         acoef(k,i)=0.0
       enddo
     enddo      
     
!--- Structures --------------------
!    call struct_heat
    
    return
    end subroutine bound_heat
    