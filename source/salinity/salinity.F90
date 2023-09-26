!=================================================================
! CMS Salinity Transport routines
!
! Contains the following:
!   sal_default - Sets default variables for salinity transport
!   sal_cards   - Reads the salinity transport cards from the 
!                 control file
!   sal_init    - Initializes the salinity transport variables
!   sal_bc_init - Initializes the salinity boundary conditions
!   sal_laplace - Solves a lapace equation using a scatter set 
!                 as boundary conditions to obtain the initial 
!                 condition for salinity transport
!   coeffinitsal - Calculates the coefficients for salinity 
!                 initial conditions which is a Laplace equation
!   sal_print   - Prints the salinity transport settings to the
!                 screen and diagnositic file
!   sal_block   - Reads a salinity boundary condition block
!   salstr_resize - Resizes the salinity transport BC variable 
!   sal_imp     - Solves the salinity transport equation
!   coeffsourcesink_sal - Assembles the coefficient matrix and 
!                   R.H.S. for salinity transport
!   bound_sal   - Applies BC's to salinity transport equation
!   bndsaleval  - Evaluate salinity at boundary cell strings
!   boundinit_sal - Applies boundary condition to salinity 
!                 transport equation for initial condition
!
! written by Alex Sanchez, USACE-CH
!            Weiming Wu, NCCHE
!        and Honghai, Li, USACE-CHL
!=================================================================
    
!**************************************************************
    subroutine sal_default
! Sets default variables for salinity transport
! written by Alex Sanchez, USACE-CHL
!**************************************************************    
    use sal_def
    implicit none
    
    saltrans = .false.    
    nsalstr = 0    
    itermaxsal = 30
    tolsal = 1.0e-5
    salic = 35.0 !ppt
    icsal = 1    !1-Constant, 2-Dataset, 3-Scatter set
    schmidtsal = 1.0 !Schmidt number for salinity
    
    return
    end subroutine sal_default

!**************************************************************
    subroutine sal_cards(cardname,foundcard)
! Reads the salinity cards.
! written by Alex Sanchez, USACE0CHL
!**************************************************************    
    use geo_def, only: grdfile
    use comvarbl, only: flowpath
    use flow_def, only: watersalt
    use out_def, only: outlist
    use sal_def
    implicit none
    
    integer :: i,ierr
    character(len=32) :: cardname,cdum
    character*120, allocatable :: TEMPPATH(:)
    real :: rdum
    logical :: foundcard
    
    foundcard = .true.
    select case(cardname)
    case('WATER_SALINITY')
      call card_scalar(77,'ppt','ppt',salic,ierr)
      watersalt = salic
      
    case('SALINITY_SCHMIDT_NUMBER')  
      backspace(77)
      read(77,*) cardname, schmidtsal    
      
    case('SALT_OUT_TIMES_LIST')
      backspace(77)
      read(77,*) cardname, outlist(6)%ilist  
          
    case('CALC_SALINITY')
      call card_boolean(77,saltrans,ierr)     
      
    case('SALINITY_MAX_ITERATIONS','SALINITY_MAX_ITER')
      backspace(77)
      read(77,*) cardname, itermaxsal      
      
    case('SALINITY_IC_CONSTANT','SALINITY_IC')
      call card_scalar(77,'ppt','ppt',salic,ierr)  !Constant initial concentration
      salinitconst = .true.
      icsal = 1  
      
    case('SALINITY_IC_DATASET') !Spatially variable dataset
      call card_dataset(77,grdfile,flowpath,salfile,salpath,1)
      salinitconst = .false.  
      icsal = 2
        
    case('SALINITY_IC_ASCII','SALINITY_IC_SCATTERSET')
      backspace(77)
      read(77,*) cardname, salobsfile
      salinitconst = .false.
      icsal = 3      
      
    case('SALINITY_CELLSTRING')  
      nsalstr = nsalstr + 1
      if(nsalstr>1)then 
        allocate(temppath(nsalstr-1))
        do i=1,nsalstr-1
          temppath(i) = saldrvpath(i)
        enddo
        deallocate(saldrvpath)
        allocate(saldrvpath(nsalstr))
        do i=1,nsalstr-1        
          saldrvpath(i) = temppath(i)
        enddo
        deallocate(TEMPPATH)
      else
        allocate(saldrvpath(1))  
      endif
      backspace(77)
      read(77,*) cardname, cdum, saldrvpath(nsalstr)  
      
    case('SALINITY_BOUNDARY_BEGIN', 'SALT_BOUNDARY_BEGIN')
      call sal_bnd_block
      
    CASE ('SALINITY_CALC_INTERVAL')
      backspace(77)
      READ (77,*) CARDNAME, rdum        !Not really used for Implicit - ignore the value 05/04/2016 meb
        
    case default
      foundcard = .false.
          
    end select
    
    return
    end subroutine sal_cards
    
!*******************************************    
    subroutine sal_bnd_block
!*******************************************
    use sal_def
    use comvarbl, only: mpfile,flowpath
    implicit none
    integer :: K,ierr
    character(len=32) :: cardname
    logical :: foundcard
    
    call salstr_resize
    
    foundcard = .true.
d1: do k=1,10
      foundcard = .true.
      read(77,*,iostat=ierr) cardname
      if(ierr/=0) exit
      if(cardname(1:1)=='!' .or. cardname(1:1)=='#') cycle
      select case(cardname)  
      case('SALINITY_CURVE')
        call card_dataset(77,mpfile,flowpath,sal_str(nsalstr)%salfile,sal_str(nsalstr)%salpath,1)
        
      case('SALINITY_CONSTANT')
        backspace(77)
        read(77,*) cardname,sal_str(nsalstr)%salbnd
        
      case('SALINITY_BOUNDARY_END','SALT_BOUNDARY_END','END')
          exit d1
            
      case default
        foundcard = .false.
        
      end select
    enddo d1
      
    return
    end subroutine sal_bnd_block
    
!**********************************************************    
    subroutine sal_block(isaltype,salfile,salpath,salbnd)
!**********************************************************
    use comvarbl, only: mpfile,flowpath
    use prec_def
    implicit none
    !Input/Output
    integer,intent(inout) :: isaltype
    character(len=*),intent(inout) :: salfile,salpath
    real(ikind),intent(inout) :: salbnd
    !Internal Variables
    integer :: k,ierr
    character(len=32) :: cardname
    logical :: foundcard
    
    foundcard = .true.
d1: do k=1,10
      foundcard = .true.
      read(77,*,iostat=ierr) cardname
      if(ierr/=0) exit
      if(cardname(1:1)=='!' .or. cardname(1:1)=='#') cycle
      select case(cardname)  
      case('SALINITY_CURVE')
        call card_dataset(77,mpfile,flowpath,salfile,salpath,1)
        isaltype = 2
        
      case('SALINITY_CONSTANT')
        backspace(77)
        read(77,*) cardname,salbnd
        isaltype = 1
        
      case('SALINITY_END','SALT_END','SAL_END','END')
          exit d1
            
      case default
        foundcard = .false.
        
      end select
    enddo d1
      
    return
    end subroutine sal_block
    
!************************************************************    
    subroutine salstr_resize
! Resizes the salinity transport boundary condition variable    
!************************************************************
    use sal_def
    use prec_def
    implicit none
    integer :: i
    type(sal_driver), allocatable :: sal_temp(:)    
    
    nsalstr = nsalstr + 1
    if(nsalstr==1)then  
      allocate(sal_str(nsalstr))
    else
      allocate(sal_temp(nsalstr-1))
      do i=1,nsalstr-1
        sal_temp(i) = sal_str(i)
      enddo
      deallocate(sal_str)
      allocate(sal_str(nsalstr))
      do i=1,nsalstr-1        
        sal_str(i) = sal_temp(i)
      enddo
      deallocate(sal_temp)
    endif
    
    !Initialize and set default values
    sal_str(nsalstr)%bidfile = ''
    sal_str(nsalstr)%bidpath = ''
    sal_str(nsalstr)%ncells = 0
    sal_str(nsalstr)%ntimes = 0
    sal_str(nsalstr)%inc = 1
  
    return
    end subroutine salstr_resize
    
!****************************************************************************    
    subroutine sal_print()
! Prints the salinity transport settings to the screen and diagnositic file
! written by Alex Sanchez, USACE-CHL
!****************************************************************************
    use sal_def
    use diag_def, only: dgunit,dgfile
    implicit none
    integer :: i,iunit(2)

787 format(' ',A,T40,A)
153 format(' ',A,T40,F0.3)
236 format(' ',A,T40,I0)
    
    iunit = (/6,dgunit/)
    
    open(dgunit,file=dgfile,access='append') 
    do i=1,2
      write(iunit(i),*)
      if(.not.saltrans)then
        write(iunit(i),787)   'Salinity Transport:','OFF'   
      else
        write(iunit(i),787)   'Salinity Transport:','ON'    
        if(icsal==1)then
          write(iunit(i),153) '  Initial Salinity:',salic
        elseif(icsal==2)then
          write(iunit(i),787) '  Initial Condition File:',trim(salfile)
          write(iunit(i),787) '  Initial Condition Dataset:',trim(salpath)
        else
          write(iunit(i),787) '  Initial Condition Scatter:',trim(salobsfile) 
        endif
        write(iunit(i),236)   '  Salinity Max Iterations:',itermaxsal
      endif
    enddo
    close(dgunit)
    
    return
    end subroutine sal_print

!**************************************************************
    subroutine sal_init
! Solves the salinity transport equation
! Alex Sanchez, USACE-CHL; Weiming Wu, NCCHE
!**************************************************************        
#include "CMS_cpp.h"
    use size_def
    use geo_def, only: grdfile
    use comvarbl, only:  ntsch
    use sal_def
#ifdef XMDF_IO
    use in_xmdf_lib, only: readscalh5
#endif   
    use in_lib, only: readscalTxt
    
    implicit none
    character(len=10):: aext
    integer :: ierr,i,j,isal

    !------ Boundary Conditions -------------------------
    call sal_bc_init
    
    !------ Initialize Condition ------------------------
    allocate(sal(ncellsD),sal1(ncellsD))   
    sal=0.0; sal1=0.0
    allocate(dsalx(ncellsD),dsaly(ncellsD))
    dsalx=0.0; dsaly=0.0
    
    !----- Static Source term -----------------------------------------
    allocate(susal0(ncellsD))
    susal0=0.0
    
    !---- Residuals ---------------
    allocate(rssal(ncellsD))
    rssal = 0.0
    
    !----- Initial Condition -----------------------------------
    if(icsal==1)then !Constant value
      sal = salic
    elseif(icsal==2)then !Specified dataset
      call fileext(trim(salfile),aext)      
      select case (aext)
      case('h5')
#ifdef XMDF_IO      
      call readscalh5(salfile,salpath,sal,ierr)
#endif
      case('txt')
        call readscalTxt(salfile,sal,ierr)
      end select
      
      if(ierr/=0)then
        write(*,*) 'ERROR: Check path for salinity initial concentration dataset'
        write(*,*) 'Press <enter> key to continue.'
        read(*,*)
        stop
      endif
    else              !Laplace interpolation of scatter set        
     call sal_laplace
    endif
    
    !--- Previous Time Step Salinities -------------
    sal1 = sal
    if(ntsch==2)then
      allocate(sal2(ncellsD))
      sal2=sal1
    endif
    
    !--- Boundary Initial Values ------------------------
    do isal=1,nsalstr
      allocate(sal_str(isal)%salbnd0(sal_str(isal)%ncells))
      do j=1,sal_str(isal)%ncells
        i=sal_str(isal)%cells(j)
        sal_str(isal)%salbnd0(j)=sal(i)
      enddo  
    enddo
    
    return
    end subroutine sal_init

!**************************************************************
    subroutine sal_laplace
! Solves a lapace equation using a scatter set as boundary conditions
! to obtain the initial condition for salinity
! Weiming Wu, NCCHE; Alex Sanchez, USACE-CHL
!**************************************************************    
    use size_def
    use geo_def, only: xOrigin,yOrigin,azimuth_fl,x,y,dx,dy,icol,irow
    use comvarbl, only: ctime,stime,rmom,ntsch
    use flow_def, only: acoef,ap,sp,su
    use sal_def
    use const_def, only: small,deg2rad
    use diag_lib
    use prec_def   
    implicit none
    integer :: i,j,k,itersal
    real(ikind) :: tempxyz(10000,3),azimuth_r,cosang,sinang
    real(ikind) :: x_global,y_global,sumdis,sumsal,factobs,distobs2
    logical :: ok    
    
    !Solve laplace equation
    !Read scatter set
    inquire(file=salobsfile,exist=ok)
    if(.not.ok)then
      call diag_print_error('Salinity initial condition scatterset not found: ', &
        salobsfile,'  Check CMS-Flow card: SALINITY_IC_ASCII')
    endif
    open(65,file=salobsfile)
    read(65,*) !Skip header file
    nsalobs = 0
    do i=1,10000
      read(65,*,end=342) (tempxyz(i,j),j=1,3)
      nsalobs = nsalobs + 1
    enddo
342 close(65)

    numsalbnd = nsalstr + nsalobs !Obervation cells are like boundaries
    allocate(salobspts(nsalobs,4))
    do i=1,nsalobs
      salobspts(i,1:3) = tempxyz(i,1:3)
    enddo
      
    !Convert global to local coordinates
    azimuth_r = azimuth_fl*deg2rad
    cosang = cos(azimuth_r)
    sinang = sin(azimuth_r)
    do i=1,nsalobs
      x_global = salobspts(i,1)
      y_global = salobspts(i,2)
      salobspts(i,1) =  (x_global-xOrigin)*cosang+(y_global-yOrigin)*sinang !x-local
      salobspts(i,2) = -(x_global-xOrigin)*sinang+(y_global-yOrigin)*cosang !y-local
    enddo
      
    !Find cells containing points
    do j=1,nsalobs
      salobspts(j,4) = -1
      do i=1,ncells        
        if(abs(salobspts(j,1)-x(i))<=0.5*dx(i) .and. &
           abs(salobspts(j,2)-y(i))<=0.5*dy(i))then
          salobspts(j,4) = i
          write(*,*) i,icol(i),irow(i)
          exit
        endif
      enddo
      if(salobspts(j,4)<0)then
        call diag_print_error('Found salinity initial condition measurement point outside of computational grid')
      endif
    enddo
      
    !Initialize using average of observation points  
!    salic = 0.0
!    do i=1,nsalobs
!      salic = salic + salobspts(i,3)
!    enddo
!    do isal=1,nsalstr
!      salic = salic + sal_str(isal)%val(1)
!    enddo
!    salic = salic/real(nsalobs+nsalstr)   
!!    if(nsalobs>=1) then
!!      salic = salic/real(nsalobs)   
!!    else
!!      do isal=1,nsalstr
!!        salic = salic + sal_str(isal)%val(1)
!!      enddo
!!      salic = salic/real(nsalstr)   
!!    endif
!    sal = salic 
    
    do i=1,ncellsD
       sumdis=0.0
       sumsal=0.0
       do k=1,nsalobs
          distobs2=(x(i)-salobspts(k,1))**2+(y(i)-salobspts(k,2))**2+small      
          factobs=1.0/distobs2
          sumdis=sumdis+factobs
          sumsal=sumsal+factobs*salobspts(k,3)
       enddo
       sal(i)=sumsal/sumdis    !Inverse-distance-squared interpolation   !Wu
    enddo

552 format(5x,I10,4x,E13.4,1x)    
    
    ctime = stime       
    call bndsaleval        
    write(*,*) 'Salinity Initial Condition' 
    write(*,*) '          Iteration    Residual'
    do itersal=1,itermaxsal*3
      call coeffinitsal                  !for Laplace equation  
      call boundinit_sal
      call solve(acoef,su,sp,rssal,sal,5)  
      if(mod(itersal,5)==0) write(*,552) itersal,rmom(5)    
      if(rmom(5)<=tolsal) exit
    enddo 
    
    !Copy to ghost cells at forcing BC's   
    call bndcopy2ghost(sal)
    
    return
    end subroutine sal_laplace

!**************************************************************
    subroutine sal_imp
! Solves the salinity transport equation
! Weiming Wu, NCCHE; Alex Sanchez, USACE-CHL
!**************************************************************    
#include "CMS_cpp.h"    
    use comp_lib
    use comvarbl, only: ndsch,rmom,skewcor,timehrs
    use der_def, only: nder,nlim,gow
    use der_lib, only: der_grad_eval    
    use diag_def
    use diag_lib
    use flow_def
    use sal_def
    use size_def
!! added 06/22/2016    
    use out_def
#ifdef XMDF_IO    
    use out_lib, only: writescalh5
#endif
!!
    
    implicit none  
    integer :: itersal, nn, i
    logical, save :: IC_written = .false.
    character(len=200) :: apath
    real(ikind) :: tempsal(ncellsD)    !Added 6/21/2016 meb

    !if (.not.IC_written) then
    !  nn = len_trim(simlabel)
    !  apath = simlabel(1:nn)//'/'
    !  call writescalh5(outlist(6)%afile,apath,'Salinity_IC',sal,'ppt',timehrs,1)
    !  IC_written = .true.
    !endif

    !Evaluate salinity boundary conditions
    call bndsaleval         
    !Compute matrix coefficients and source/sink terms
    select case(ndsch) 
    case(2); call coeffsourcesink_sal(hybridcoef)
    case(3); call coeffsourcesink_sal(powerlawcoef)
    case(4); call coeffsourcesink_sal(exponentialcoef)
    case default; call coeffsourcesink_sal(upwindcoef)
    end select
    !Apply boundary conditions
    call bound_sal  
    !Check matrix and source terms
    if(debug_mode) call check_variables(5)
    
    !call diag_print_message(' Salinity: itersal  ResidualSal')
    call diag_print_message(' Salinity: iteration     residual')
    rmom(5)=100.0
    itersal=0
    tempsal = sal    !Added 6/21/2016 meb
    do itersal=1,itermaxsal
      !Reset source term without deferred corrections
      su=susal0
      !Compute deferred corrections
      if(ncellsimple==ncells)then !No gradients required
        select case(ndsch) !Anti-diffusion corrections
        case(5); call defcorhlpa(sal,su)
        case(6); call defcorgamma(gammadefcor,sal,su)
        case(7); call defcorgamma(cubistadefcor,sal,su)
        case(8); call defcorgamma(alvsmartdefcor,sal,su)
        case(9); call defcorgamma(hoabdefcor,sal,su)
        end select
      else 
        call der_grad_eval(gow,nlim,sal,dsalx,dsaly)
        select case(ndsch) !Anti-diffusion corrections
        case(5); call defcorhlpagrad(sal,dsalx,dsaly,su)
        case(6); call defcorgammagrad(gammadefcor,sal,dsalx,dsaly,su)
        case(7); call defcorgammagrad(cubistadefcor,sal,dsalx,dsaly,su)
        case(8); call defcorgammagrad(alvsmartdefcor,sal,dsalx,dsaly,su)
        case(9); call defcorgammagrad(hoabdefcor,sal,dsalx,dsaly,su)
        case default; if(skewcor) call defcorparagrad(dsalx,dsaly,su) !Skewness correction 
        end select
      endif
      !Solve matrix for salinity
      call solve(acoef,su,sp,rssal,sal,5)

      !Screen Output
      if(itersal==1 .or. mod(itersal,10)==0)then
        write(*,555) itersal,rmom(5)
      endif  
      !Check for convergence
      if((rmom(5)<tolsal) .and. itersal>5)then
        exit  
      endif
    enddo
    
!$OMP PARALLEL DO PRIVATE (i)                !Added 6/21/2016 meb
    do i=1,ncellsD
      if(iwet(i)==0) sal(i) = tempsal(i)
    enddo
!$OMP END PARALLEL DO
    
    if(mod(itersal,10)/=0) write(*,555) itersal,rmom(5)
    open(dgunit,file=dgfile,access='append')
    write(dgunit,555) itersal,rmom(5)
    close(dgunit)                

!    call bndcopy2ghost(sal) !Copy sal to dummy/ghost cells

555 format(4x,I10,1x,E13.4,1x)
    
    call sal_step_stat
    
    return
    end subroutine sal_imp
    
!**************************************************************
    subroutine coeffinitsal
! Calculates the coefficients for salinity initial conditions
! which is a Laplace equation
! written by Alex Sanchez, USACE-CHL
!**************************************************************    
    use size_def
    use geo_def, only: ncface,dsxy
    use flow_def, only: acoef,su,sp
    use sal_def
    implicit none
    integer :: i,k
    
    do i=1,ncells
       do k=1,ncface(i)
         acoef(k,i)=dsxy(k,i)
       enddo  
       su(i) = 0.0
       sp(i) = 0.0
    enddo        
    
    return
    end subroutine coeffinitsal

!***********************************************************************
    subroutine boundinit_sal
! Applies boundary condition to salinity transport equation for initial condition
! Alex Sanchez, USACE-CHL; Weiming Wu, NCCHE
!***********************************************************************
    use size_def
    use geo_def
    use flow_def
    use struct_def
    use comvarbl
    use sal_def
    implicit none
    integer :: i,j,k,isal

!   Dry nodes
    do i=1,ncells
      if(iwet(i)==0) then !Dry
        sal(i)=0.0
        acoef(:,i)=0.0
        sp(i)=-1.0
        su(i)=0.0                     
      endif
    enddo

!   Wet cells next to dry cells      
    do i=1,ncells
      if(iwet(i)==1) then !Wet
        do k=1,ncface(i)
          if(iwet(cell2cell(k,i))==0) then
            acoef(k,i)=0.0
          endif
        enddo             
      endif
    enddo
  
 !  Salinity boundary
    do isal=1,nsalstr
      do j=1,sal_str(isal)%ncells
        i=sal_str(isal)%cells(j)
        if(i==0) cycle !Land cell
        acoef(:,i)=0.0
        sp(i)=-1.0
        su(i)=sal_str(isal)%salbnd                     
      enddo
    enddo   

!   Observation Stations
    do isal=1,nsalobs
      i=salobspts(isal,4)
      acoef(:,i)=0.0
      sp(i)=-1.0
      su(i)=salobspts(isal,3)
    enddo
    
    return
    end subroutine boundinit_sal
    
!*************************************************************************
    subroutine coeffsourcesink_sal(schmcoef)
! Assembles the coefficient matrix and R.H.S. for salinity transport
! written by Alex Sanchez, USACE-CHL
!*************************************************************************
    use size_def
    use geo_def, only: areap,ncface,cell2cell,dsxy,idirface
    use flow_def, only: h1,h2,hk,iwet,acoef,sp,flux,visk,viskfl
    use comvarbl, only: dtime,ntsch,ctsch1,ctsch2
    use sal_def
    use comp_lib
    use prec_def
    implicit none
    integer :: i,k
    real(ikind):: ddk,val,fac1,fac2,dtimeinv,schmidtsalinv
    
    interface
      function schmcoef(dk,fk)
        use prec_def
        implicit none
        real(ikind),intent(in) :: dk,fk
        real(ikind) :: schmcoef
      end function
    endinterface    

    dtimeinv=1.0/dtime
    schmidtsalinv=1.0/schmidtsal
    if(ntsch==1)then
!$OMP PARALLEL DO PRIVATE(i,k,ddk,fac1)
      do i=1,ncells
        fac1 = iwet(i)*areap(i)*h1(i)*dtimeinv   !MEB modified 06/21/2016, so that only wet cells are included in salinity calculations
        sp(i)=-fac1
        susal0(i)=fac1*sal1(i)
        do k=1,ncface(i)          
          ddk=viskfl(k,i)*hk(k,i)*dsxy(k,i)*schmidtsalinv
          acoef(k,i)=schmcoef(ddk,flux(k,i))
        enddo
      enddo   
!$OMP END PARALLEL DO
    else
!$OMP PARALLEL DO PRIVATE(i,k,val,fac1,fac2,ddk)    
      do i=1,ncells
        val=iwet(i)*areap(i)*dtimeinv            !MEB modified 06/21/2016, so that only wet cells are included in salinity calculations
        fac1=val*ctsch1*h1(i)
        fac2=val*ctsch2*h2(i)
        sp(i)=-fac1+fac2
        susal0(i)=fac1*sal1(i)-fac2*sal2(i)
        do k=1,ncface(i)          
          ddk=viskfl(k,i)*hk(k,i)*dsxy(k,i)*schmidtsalinv
          acoef(k,i)=schmcoef(ddk,flux(k,i))
        enddo 
      enddo
!$OMP END PARALLEL DO
    endif
    
    return
    end subroutine coeffsourcesink_sal

    
!******************************************************************
    subroutine sal_step_stat
! Calculates the salinity transport statistics for each time step.
! written by Alex Sanchez, USACE-CHL
!******************************************************************   
    use size_def
    use geo_def, only: mapid
    use flow_def, only: iwet
    use sal_def
    use diag_def, only: dgfile,dgunit
    use prec_def
    implicit none
    integer :: i,idsalmax,idsalmax2,idsalmin,idsalmin2,nwet
    real(ikind) :: salmax,salmax2,salmin,salmin2,salavg
    
    !Initialize
    salmax=-1.0e6; salmax2=-1.0e6; idsalmax=1; idsalmax2=1
    salmin= 1.0e6; salmin2= 1.0e6; idsalmin=1; idsalmin2=1
    salavg=0.0; nwet = 0
    
!!$OMP PARALLEL PRIVATE(idsalmax2,salmax2,idsalmin2,salmin2)
!!$OMP DO PRIVATE(i) REDUCTION(+:salavg) REDUCTION(+:nwet)
    do i=1,ncells
      if(iwet(i)==1)then
        nwet=nwet+1
        salavg=salavg+sal(i)        
        if(sal(i)>salmax2)then
          salmax2=sal(i)
          idsalmax2=i
        endif
        if(sal(i)<salmin2)then
          salmin2=sal(i)
          idsalmin2=i
        endif
      endif
    enddo
!!$OMP END DO
!!$OMP CRITICAL
    if(salmax2>salmax)then
      salmax=salmax2
      idsalmax=idsalmax2
    endif
    if(salmin2<salmin)then
      salmin=salmin2
      idsalmin=idsalmin2
    endif
!!$OMP END CRITICAL
!!$OMP END PARALLEL

    nwet = max(nwet,1)
    salavg = salavg/real(nwet,kind=ikind)
    !write(*,626) salavg,idsalmax,salmax,idsalmin,salmin    
    idsalmax = mapid(idsalmax)
    idsalmin = mapid(idsalmin)
    
626 format('   salavg=',F6.3,', salmax(',I6,')=',F6.3,', salmin(',I6,')=',F6.3)
    !if(abs(Ctmax)>1.0e5 .or. abs(dzbxtr)>1.0e3) return    
    open(dgunit,file=dgfile,access='append')
    write(*,626)      salavg,idsalmax,salmax,idsalmin,salmin
    write(dgunit,626) salavg,idsalmax,salmax,idsalmin,salmin
    close(dgunit)      
    
    return
    end subroutine sal_step_stat    
    