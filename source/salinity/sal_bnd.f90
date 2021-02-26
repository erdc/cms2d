!**************************************************************
    subroutine sal_bc_init
! Initializes the salinity boundary conditions
! Alex Sanchez, USACE-CHL; Weiming Wu, NCCHE
!**************************************************************
    use sal_def
    use bnd_def
    use diag_lib
    implicit none
    integer :: i,isal
    
!=== Find all salinity boundaries ================
    !Assumes that salinity is specified at all boundaries
    nsalstr = nQstr + nTHstr + nHstr + nMHstr + nMHVstr + nCSstr &
            + nNHstr + nNHVstr + nNTHstr + nNTHVstr
    allocate(sal_str(nsalstr))
    sal_str%inc = 1 !Initialize
    isal = 0
    
!--- River BC -------------------------------------------
    do i=1,nQstr        
      call sal_bnd_chk(Q_str(i)%bidfile,Q_str(i)%bidpath,&
         Q_str(i)%ncells,Q_str(i)%cells,Q_str(i)%faces,isal)  
    enddo

!--- Tidal/Harmonic BC ---------------------------------
    do i=1,nTHstr
      call sal_bnd_chk(TH_str(i)%bidfile,TH_str(i)%bidpath,&
         TH_str(i)%ncells,TH_str(i)%cells,TH_str(i)%faces,isal)
    enddo
    
!--- Single Water Level BC --------------------------
    do i=1,nHstr  
      call sal_bnd_chk(H_str(i)%bidfile,H_str(i)%bidpath,&
         H_str(i)%ncells,H_str(i)%cells,H_str(i)%faces,isal)
    enddo

!--- Multiple Water Level BC ------------------------------
    do i=1,nMHstr
      call sal_bnd_chk(MH_str(i)%bidfile,MH_str(i)%bidpath,&
         MH_str(i)%ncells,MH_str(i)%cells,MH_str(i)%faces,isal)
    enddo
    
!--- Multiple Water Level and Velocity BC -------------------------
    do i=1,nMHVstr
      call sal_bnd_chk(MHV_str(i)%bidfile,MHV_str(i)%bidpath,&
         MHV_str(i)%ncells,MHV_str(i)%cells,MHV_str(i)%faces,isal)
    enddo 

!--- Cross-shore BC (Type CS=6) --------------------------
    do i=1,nCSstr 
      call sal_bnd_chk(CS_str(i)%bidfile,CS_str(i)%bidpath,&
         CS_str(i)%ncells,CS_str(i)%cells,CS_str(i)%faces,isal)  
    enddo
    
!--- Nested Water Level BC (Type 7=NH) ----------------------------
    do i=1,nNHstr    
      call sal_bnd_chk(NH_str(i)%bidfile,NH_str(i)%bidpath,&
         NH_str(i)%ncells,NH_str(i)%cells,NH_str(i)%faces,isal)  
    enddo        
        
!--- Nested Water Level and Velocity BC (Type 8=NHV) ---------------
    do i=1,nNHVstr    
      call sal_bnd_chk(NHV_str(i)%bidfile,NHV_str(i)%bidpath,&
         NHV_str(i)%ncells,NHV_str(i)%cells,NHV_str(i)%faces,isal)  
    enddo 
    
!--- Nested Tidal Database WSE BC (Type 9=NTH) -----------------------
    do i=1,nNTHstr    
      call sal_bnd_chk(NTH_str(i)%bidfile,NTH_str(i)%bidpath,&
         NTH_str(i)%ncells,NTH_str(i)%cells,NTH_str(i)%faces,isal)  
    enddo 
    
!--- Nested Tidal Database WSE and Velocity BC (Type 10=NTHV) ----------  
    do i=1,nNTHVstr    
      call sal_bnd_chk(NTHV_str(i)%bidfile,NTHV_str(i)%bidpath,&
         NTHV_str(i)%ncells,NTHV_str(i)%cells,NTHV_str(i)%faces,isal)  
    enddo 
    
    !Check
    if(isal/=nsalstr)then
      call diag_print_error('Problem determining the salinity boundary conditions')
    endif
    
    return
    endsubroutine sal_bc_init
    
!************************************************************************
    subroutine sal_bnd_chk(bidfile,bidpath,nbndcells,icells,kfaces,isal)
! Checkes the flow boundary condition file and paths 
! for salinity information
! written by Alex Sanchez, USACE-CHL
!************************************************************************
#include "CMS_cpp.h"
    use size_def
    use prec_def, only: ikind
    use geo_def, only: cell2cell,ncface
    use comvarbl, only: mpfile
    use in_lib, only: read_xys
    use diag_lib
    use sal_def
    use bnd_def
#ifdef XMDF_IO
    use xmdf
#endif
    implicit none
    !Input/Output
    character(len=*),intent(in) :: bidfile,bidpath
    integer, intent(in) :: nbndcells
    integer, intent(in) :: icells(nbndcells),kfaces(nbndcells)
    integer, intent(inout) :: isal
    !Internal Variables
    integer :: nstimes,ierr,PID,GID
    character(len=100) :: apath,afile,aext,abnd
    character(len=200) :: asalfile
    real(ikind),pointer :: stimes(:),sval(:)
    
    call fileext(bidfile,aext)
    select case (aext)
    case('h5')
#ifdef XMDF_IO
      call XF_OPEN_FILE(bidfile,READONLY,PID,ierr)
      call XF_OPEN_GROUP(PID,bidpath,GID,ierr)
      call XF_GET_PROPERTY_NUMBER(GID,'Sal_Times',nstimes,ierr)
      if(ierr<0)then
        call diag_print_error('Salinity must be specified at multi-wse boundary conditions')
      endif
      isal = isal + 1
      sal_str(isal)%ntimes = nstimes
      allocate(sal_str(isal)%val(nstimes))
      allocate(sal_str(isal)%timesal(nstimes))
      call XF_READ_PROPERTY_FLOAT(GID,'Sal_Times',nstimes,sal_str(isal)%timesal(1),ierr)
      call XF_READ_PROPERTY_FLOAT(GID,'Salinity',nstimes,sal_str(isal)%val(1),ierr)
      call XF_CLOSE_GROUP(GID,ierr)
      call XF_CLOSE_FILE(PID,ierr)
#endif
    case('bid')
      isal = isal + 1                           !First time through this is 0.  Needs to be incremented to read from correct file.
      call fileparts (bidfile,apath,afile,aext)
      write(abnd,'(I0)') isal
      asalfile = trim(apath) // trim(afile) // '_sal_' // trim(abnd) // '.xys'
      call read_xys(asalfile,nstimes,stimes,sval)
      sal_str(isal)%ntimes = nstimes
      allocate(sal_str(isal)%val(nstimes))
      allocate(sal_str(isal)%timesal(nstimes))
      sal_str(isal)%timesal = stimes
      sal_str(isal)%val = sval
        
    end select
    
    !Cell id's and faces
    sal_str(isal)%ncells = nbndcells
    allocate(sal_str(isal)%cells(nbndcells))
    allocate(sal_str(isal)%faces(nbndcells))
    sal_str(isal)%cells = icells
    sal_str(isal)%faces = kfaces
    
    return
    endsubroutine sal_bnd_chk
    
!**************************************************************
    subroutine bndsaleval
! Evaluate salinity at boundary cell strings
! Alex Sanchez, USACE-CHL
!**************************************************************    
    use size_def
    use flow_def
    use comvarbl, only: timehrs,ramp
    use sal_def
    use const_def, only: eps
    use prec_def
    implicit none
    integer :: i,ntimes,inc  
    real(ikind) :: fac
    
    do i=1,nsalstr
      if(sal_str(i)%ntimes==0) cycle !Boundary value is constant
      !find out where we are in the time/value arrays
      inc = sal_str(i)%inc
      ntimes = sal_str(i)%ntimes
      do while((timehrs+eps)>=sal_str(i)%timesal(inc+1) .and. inc<ntimes)
        inc = inc + 1
        sal_str(i)%inc = inc
      enddo
      fac = (timehrs-sal_str(i)%timesal(inc))/ &
            (sal_str(i)%timesal(inc+1)-sal_str(i)%timesal(inc))
      sal_str(i)%salbnd = (1.0-fac)*sal_str(i)%val(inc) + &
                fac*sal_str(i)%val(inc+1)
    enddo
    
    return
    endsubroutine bndsaleval
    
!***********************************************************************
    subroutine bound_sal
! Applies boundary condition to salinity transport equation
! Weiming Wu, NCCHE;  Alex Sanchez, USACE-CHL
!***********************************************************************
    use size_def
    use geo_def
    use flow_def
    use sal_def
    use comvarbl, only: ramp
    use prec_def
    implicit none
    integer :: i,j,k,isal,nck,idrym

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
!        apsal=0.0
        susal0(i)=0.0
        do k=1,ncface(i)
          nck=cell2cell(k,i)
          if(idrym>=1)then
            acoef(k,i)=iwet(nck)
          else
            acoef(k,i)=dsxy(k,i)
          endif
          if(nck>ncells) acoef(k,i)=0.0
          !apsal=apsal+acoef(k,i)
          !susal0(i)=susal0(i)+acoef(k,i)*sal(nck)
        enddo
        !susal0(i)=susal0(i)-apsal*sal(i)
        sp(i)=0.0
      endif
    enddo
 !$OMP END PARALLEL DO
 
!--- Salinity boundary ---------------------
    do isal=1,nsalstr
      do j=1,sal_str(isal)%ncells
        i=sal_str(isal)%cells(j)
        !if(i==0) cycle !Land cell   
         k=sal_str(isal)%faces(j)
        !do k=1,ncface(i)
          nck=cell2cell(k,i)
          !if(nck>ncells)then
            if(flux(k,i)<0.0)then !Inflow
               sal(nck)=ramp*sal_str(isal)%salbnd+(1-ramp)*sal_str(isal)%salbnd0(j)
               susal0(i)=susal0(i)+acoef(k,i)*sal(nck)
               sp(i)=sp(i)-acoef(k,i)
            else  !Outflow
              sal(nck)=sal(i)
            endif
            acoef(k,i)=0.0
          !endif
        !enddo
      enddo
    enddo   

!--- Structures --------------------
    call struct_sal
    
    return
    endsubroutine bound_sal
    