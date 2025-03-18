!=================================================================
! CMS Heat Transfer routines
!
! Contains the following:
!   heat_default - Sets default variables for heat transfer
!   heat_cards   - Reads the heat transfer cards from the 
!                 control file
!   heat_init    - Initializes the heat transfer variables
!   heat_bc_init - Initializes the heat boundary conditions
!   heat_laplace - Solves a lapace equation using a scatter set 
!                 as boundary conditions to obtain the initial 
!                 condition for heat transfer
!   coeffinitheat - Calculates the coefficients for heat 
!                 initial conditions which is a Laplace equation
!   heat_print   - Prints the heat transfer settings to the
!                 screen and diagnositic file
!   heat_block   - Reads a heat boundary condition block
!   heatstr_resize - Resizes the heat transfer BC variable 
!   heat_imp     - Solves the heat transfer equation
!   coeffsourcesink_heat - Assembles the coefficient matrix and 
!                   R.H.S. for heat transfer
!   bound_heat   - Applies BC's to heat transfer equation
!   bndheateval  - Evaluate heat at boundary cell strings
!   boundinit_heat - Applies boundary condition to heat 
!                 transfer equation for initial condition
!
! written by Weiming Wu, Clarkson University
!        and Honghai, Li, USACE-CHL
!=================================================================
    
!**************************************************************
    subroutine heat_default
! Sets default variables for heat transfer

!**************************************************************    
    use heat_def
    implicit none
    
    heattrans = .false.    
    nheatstr = 0    
    itermaxheat = 30
    tolheat = 1.0e-5
    heatic = 10.0 !deg C
    icheat = 1    !1-Constant, 2-Dataset, 3-Scatter set
    schmidtheat = 1.0 !Schmidt number for heat
    heatpath = "PROPERTIES/Model Params/TemperatureParameters"
    
    return
    end subroutine heat_default

!**************************************************************
    subroutine heat_cards(cardname,foundcard)
! Reads the heat cards.
! meb 05/03/2016 - I don't like the term HEAT.  Using TEMPERATURE 
!   for INPUT/OUTPUT cards as much as possible
!**************************************************************    
    use geo_def, only: grdfile
    use comvarbl, only: flowpath
    use flow_def, only: watertemp
    use out_def, only: outlist
    use heat_def
    implicit none
    
    integer :: i,ierr
    character(len=32) :: cardname,cdum
    character*120, allocatable :: TEMPPATH(:)
    real :: rdum
    logical :: foundcard
    
    interface
      subroutine card_dataset(inunit,defaultfile,defaultpath,datafile,datapath,ndim,isboundary)	  
        integer,intent(in) :: inunit
        character(len=*),intent(in) :: defaultfile,defaultpath
        character(len=*),intent(inout) :: datafile,datapath
        integer, intent(in) :: ndim
        logical, intent(in), optional :: isboundary
      end subroutine
    end interface
        
    foundcard = .true.
    select case(cardname)
    case('WATER_HEAT')                              !Not implemented
      call card_scalar(77,'c','c',heatic,ierr)
      watertemp = heatic
      
    case('HEAT_SCHMIDT_NUMBER')                     !Not implemented
      backspace(77)
      read(77,*) cardname, schmidtheat    
      
    case('HEAT_OUT_TIMES_LIST')                     !Not implemented
      backspace(77)
      read(77,*) cardname, outlist(6)%ilist  
          
    case('CALC_HEAT', 'CALC_TEMPERATURE')
      !meb 05/03/2016
      !CMS uses constant temperature by default.  Added "icheat" and "heatinitconst" 
      !variables here in lieu of using the HEAT_IC_CONSTANT card below
      call card_boolean(77,heattrans,ierr)
      heatinitconst = .true.
      icheat = 1
      
    case('HEAT_MAX_ITERATIONS','HEAT_MAX_ITER','TEMPERATURE_MAX_ITER')  !Not implemented
      backspace(77)
      read(77,*) cardname, itermaxheat      
      
    case('HEAT_IC_CONSTANT','HEAT_IC')              !Not implemented
      call card_scalar(77,'c','c',heatic,ierr)      !Constant initial temperature
      heatinitconst = .true.
      icheat = 1  
      
    case('HEAT_IC_DATASET', 'TEMPERATURE_IC_DATASET') !Spatially variable dataset
      call card_dataset(77,grdfile,flowpath,heatfile,heatpath,1)
      heatinitconst = .false.  
      icheat = 2
        
    case('HEAT_IC_ASCII','HEAT_IC_SCATTERSET')      !Not implemented
      backspace(77)
      read(77,*) cardname, heatobsfile
      heatinitconst = .false.
      icheat = 3      
      
    case('HEAT_CELLSTRING', 'TEMPERATURE_CELLSTRING')  !Not implemented
      nheatstr = nheatstr + 1
      if(nheatstr>1)then 
        allocate(temppath(nheatstr-1))
        do i=1,nheatstr-1
          temppath(i) = heatdrvpath(i)
        enddo
        deallocate(heatdrvpath)
        allocate(heatdrvpath(nheatstr))
        do i=1,nheatstr-1        
          heatdrvpath(i) = temppath(i)
        enddo
        deallocate(TEMPPATH)
      else
        allocate(heatdrvpath(1))  
      endif
      backspace(77)
      read(77,*) cardname, cdum, heatdrvpath(nheatstr)  
      
    case('HEAT_BOUNDARY_BEGIN','TEMPERATURE_BEGIN')     !Not implemented
      call heat_bnd_block
      
    CASE ('TEMPERATURE_CALC_INTERVAL')                  !Only used for Explicit
      backspace(77)
      READ (77,*) CARDNAME, rdum        !Not really used for Implicit - ignore the value 05/04/2016 meb
        
    case default
      foundcard = .false.
          
    end select
    
    return
    end subroutine heat_cards
    
!*******************************************    
    subroutine heat_bnd_block
!*******************************************
    use heat_def
    use comvarbl, only: mpfile,flowpath
    implicit none
    integer :: K,ierr
    integer :: ii,ntimes
    character(len=32) :: cardname
    logical :: foundcard
    
    interface
      subroutine card_dataset(inunit,defaultfile,defaultpath,datafile,datapath,ndim,isboundary)	  
        integer,intent(in) :: inunit
        character(len=*),intent(in) :: defaultfile,defaultpath
        character(len=*),intent(inout) :: datafile,datapath
        integer, intent(in) :: ndim
        logical, intent(in), optional :: isboundary
      end subroutine
    end interface
        
!    call heatstr_resize
    
    foundcard = .true.
d1: do k=1,10
      foundcard = .true.
      read(77,*,iostat=ierr) cardname
      if(ierr/=0) exit
      if(cardname(1:1)=='!' .or. cardname(1:1)=='#') cycle
      select case(cardname)  
      case('HEAT_CURVE', 'TEMPERATURE_CURVE')
        call card_dataset(77,mpfile,flowpath,heat_str(nheatstr)%heatfile,heat_str(nheatstr)%heatpath,1)
        
      case('HEAT_FILE_CURVE','TEMPERATURE_FILE_CURVE')
         call heatstr_resize
         read(51,*) ntimes
         heat_str(nheatstr)%ntimes = ntimes
         allocate (heat_str(nheatstr)%timeheat(ntimes),heat_str(nheatstr)%val(ntimes))
         do ii=1,ntimes
            read(51,*) heat_str(nheatstr)%timeheat(ii),heat_str(nheatstr)%val(ii)
         enddo
         
      case('HEAT_CONSTANT')
        backspace(77)
        read(77,*) cardname,heat_str(nheatstr)%heatbnd
        
      case('HEAT_BOUNDARY_END','END')
          exit d1
            
      case default
        foundcard = .false.
        
      end select
    enddo d1
      
    return
    end subroutine heat_bnd_block
    
!************************************************************    
    subroutine heatstr_resize
! Resizes the heat transfer boundary condition variable    
!************************************************************
    use heat_def
    use prec_def
    implicit none
    integer :: i
    type(heat_driver), allocatable :: heat_temp(:)    
    
    nheatstr = nheatstr + 1
    if(nheatstr==1)then  
      allocate(heat_str(nheatstr))
    else
      allocate(heat_temp(nheatstr-1))
      do i=1,nheatstr-1
        heat_temp(i) = heat_str(i)
      enddo
      deallocate(heat_str)
      allocate(heat_str(nheatstr))
      do i=1,nheatstr-1        
        heat_str(i) = heat_temp(i)
      enddo
      deallocate(heat_temp)
    endif
    
    !Initialize and set default values
    heat_str(nheatstr)%bidfile = ''
    heat_str(nheatstr)%bidpath = ''
    heat_str(nheatstr)%ncells = 0
    heat_str(nheatstr)%ntimes = 0
    heat_str(nheatstr)%inc = 1
  
    return
    end subroutine heatstr_resize
    
!****************************************************************************    
    subroutine heat_print()
! Prints the heat transfer settings to the screen and diagnositic file
! written by Alex Sanchez, USACE-CHL
!****************************************************************************
    use heat_def
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
      if(.not.heattrans)then
        write(iunit(i),787)   'Temperature Transfer:','OFF'   
      else
        write(iunit(i),787)   'Temperature Transfer:','ON'    
        if(icheat==1)then
          write(iunit(i),153) '  Initial Water Temperature:',heatic
        elseif(icheat==2)then
          write(iunit(i),787) '  Initial Condition File:',trim(heatfile)
          write(iunit(i),787) '  Initial Condition Dataset:',trim(heatpath)
        else
          write(iunit(i),787) '  Initial Condition Scatter:',trim(heatobsfile) 
        endif
        write(iunit(i),236)   '  Temperature Max Iterations:',itermaxheat
      endif
    enddo
    close(dgunit)
    
    return
    end subroutine heat_print

!**************************************************************
    subroutine heat_init
! Solves the heat transfer equation
! 
!**************************************************************        
#include "CMS_cpp.h"
    use size_def
    use geo_def, only: grdfile
    use comvarbl, only:  ntsch
    use heat_def
#ifdef XMDF_IO
    use in_xmdf_lib, only: readscalh5
#endif   
    use in_lib, only: readscalTxt
    
    implicit none
    character(len=10) :: aext
    integer :: ierr,i,j,iheat

    !------ Boundary Conditions -------------------------
    call heat_bc_init
    
    !------ Initialize Condition ------------------------
    allocate(heat(ncellsD),heat1(ncellsD))   
    allocate(dheatx(ncellsD),dheaty(ncellsD))
    heat=0.0; heat1=0.0
    dheatx=0.0; dheaty=0.0
    
    !----- Static Source term -----------------------------------------
    allocate(suheat0(ncellsD))
    suheat0=0.0
    
    !---- Residuals ---------------
    allocate(rsheat(ncellsD))
    rsheat = 0.0
    
    !----- Initial Condition -----------------------------------
    if(icheat==1)then !Constant value
      heat = heatic
    elseif(icheat==2)then !Specified dataset
      call fileext(trim(heatfile),aext)      
      select case (aext)
      case('h5')
#ifdef XMDF_IO      
      call readscalh5(heatfile,heatpath,heat,ierr)
#endif
      case('txt')
        call readscalTxt(heatfile,heat,ierr)
      end select
      
      if(ierr/=0)then
        write(*,*) 'ERROR: Check path for temperature initial concentration dataset'
        write(*,*) 'Press <enter> key to continue.'
        read(*,*)
        stop
      endif
    else              !Laplace interpolation of scatter set        
     call heat_laplace
    endif
    
    !--- Previous Time Step Temperatures -------------
    heat1 = heat
    if(ntsch==2)then
      allocate(heat2(ncellsD))
      heat2=heat1
    endif
    
    !--- Boundary Initial Values ------------------------
    do iheat=1,nheatstr
      allocate(heat_str(iheat)%heatbnd0(heat_str(iheat)%ncells))
      do j=1,heat_str(iheat)%ncells
        i=heat_str(iheat)%cells(j)
        heat_str(iheat)%heatbnd0(j)=heat(i)
      enddo  
    enddo

    !--- Read solar radiation data ----------------------
#ifdef XMDF_IO      
    call read_heat_solar_xmdf
#else    
    call read_heat_solar      
#endif

    return
    end subroutine heat_init

!**************************************************************
    subroutine heat_laplace
! Solves a lapace equation using a scatter set as boundary conditions
! to obtain the initial condition for heat
! Weiming Wu, NCCHE; Alex Sanchez, USACE-CHL
!**************************************************************    
    use size_def
    use geo_def, only: xOrigin,yOrigin,azimuth_fl,x,y,dx,dy,icol,irow
    use comvarbl, only: ctime,stime,rmom,ntsch
    use flow_def, only: acoef,ap,sp,su
    use heat_def
    use const_def, only: small,deg2rad
    use diag_lib
    use prec_def   
    implicit none
    integer :: i,j,k,iterheat
    real(ikind) :: tempxyz(10000,3),azimuth_r,cosang,sinang
    real(ikind) :: x_global,y_global,sumdis,sumheat,factobs,distobs2
    logical :: ok    
    
   !Solve laplace equation
    !Read scatter set
    inquire(file=heatobsfile,exist=ok)
    if(.not.ok)then
      call diag_print_error('Temperature initial condition scatterset not found: ', &
        heatobsfile,'  Check CMS-Flow card: TEMPERATURE_IC_ASCII')
    endif
    open(65,file=heatobsfile)
    read(65,*) !Skip header file
    nheatobs = 0
    do i=1,10000
      read(65,*,end=342) (tempxyz(i,j),j=1,3)
      nheatobs = nheatobs + 1
    enddo
342 close(65)

    numheatbnd = nheatstr + nheatobs !Obervation cells are like boundaries
    allocate(heatobspts(nheatobs,4))
    do i=1,nheatobs
      heatobspts(i,1:3) = tempxyz(i,1:3)
    enddo
      
    !Convert global to local coordinates
    azimuth_r = azimuth_fl*deg2rad
    cosang = cos(azimuth_r)
    sinang = sin(azimuth_r)
    do i=1,nheatobs
      x_global = heatobspts(i,1)
      y_global = heatobspts(i,2)
      heatobspts(i,1) =  (x_global-xOrigin)*cosang+(y_global-yOrigin)*sinang !x-local
      heatobspts(i,2) = -(x_global-xOrigin)*sinang+(y_global-yOrigin)*cosang !y-local
    enddo
      
    !Find cells containing points
    do j=1,nheatobs
      heatobspts(j,4) = -1
      do i=1,ncells        
        if(abs(heatobspts(j,1)-x(i))<=0.5*dx(i) .and. &
           abs(heatobspts(j,2)-y(i))<=0.5*dy(i))then
          heatobspts(j,4) = i
          write(*,*) i,icol(i),irow(i)
          exit
        endif
      enddo
      if(heatobspts(j,4)<0)then
        call diag_print_error('Found temperature initial condition measurement point outside of computational grid')
      endif
    enddo
      
    !Initialize using average of observation points  
!    heatic = 0.0
!    do i=1,nheatobs
!      heatic = heatic + heatobspts(i,3)
!    enddo
!    do iheat=1,nheatstr
!      heatic = heatic + heat_str(iheat)%val(1)
!    enddo
!    heatic = heatic/real(nheatobs+nheatstr)   
!!    if(nheatobs>=1) then
!!      heatic = heatic/real(nheatobs)   
!!    else
!!      do iheat=1,nheatstr
!!        heatic = heatic + heat_str(iheat)%val(1)
!!      enddo
!!      heatic = heatic/real(nheatstr)   
!!    endif
!    heat = heatic 
    
    do i=1,ncellsD
       sumdis=0.0
       sumheat=0.0
       do k=1,nheatobs
          distobs2=(x(i)-heatobspts(k,1))**2+(y(i)-heatobspts(k,2))**2+small      
          factobs=1.0/distobs2
          sumdis=sumdis+factobs
          sumheat=sumheat+factobs*heatobspts(k,3)
       enddo
       heat(i)=sumheat/sumdis    !Inverse-distance-squared interpolation   !Wu
    enddo

552 format(5x,I10,4x,E13.4,1x)    
    
    ctime = stime       
    call bndheateval        
    write(*,*) 'Temperature Initial Condition' 
    write(*,*) '          iteration  Residual'
    do iterheat=1,itermaxheat*3
      call coeffinitheat                  !for Laplace equation  
      call boundinit_heat
      call solve(acoef,su,sp,rsheat,heat,5)  
      if(mod(iterheat,5)==0) write(*,552) iterheat,rmom(5)    
      if(rmom(5)<=tolheat) exit
    enddo 
    
    !Copy to ghost cells at forcing BC's   
    call bndcopy2ghost(heat)
    
    return
    end subroutine heat_laplace

!**************************************************************
    subroutine heat_imp
! Solves the heat transfer equation
! Weiming Wu, NCCHE; Alex Sanchez, USACE-CHL
!**************************************************************    
    use comp_lib
    use comvarbl, only: ndsch,rmom,skewcor
    use der_def, only: nder,nlim,gow
    use der_lib, only: der_grad_eval    
    use diag_def
    use diag_lib
    use flow_def
    use heat_def
    use size_def
    !use struct_def
    implicit none  
    integer :: iterheat, i
    real(ikind) :: tempheat(ncellsD)

    !Evaluate heat boundary conditions
    call bndheateval         

    call heatdatainterpol  ! Determine solar radiation at current time step
    
    !Compute matrix coefficients and source/sink terms
    select case(ndsch) 
    case(2);      call coeffsourcesink_heat(hybridcoef)
    case(3);      call coeffsourcesink_heat(powerlawcoef)
    case(4);      call coeffsourcesink_heat(exponentialcoef)
    case default; call coeffsourcesink_heat(upwindcoef)
    end select

    call heatfluxsrc
    
    !Apply boundary conditions
    call bound_heat  
    !Check matrix and source terms
    if(debug_mode) call check_variables(5)
    
    !call diag_print_message(' Heat: iterheat  ResidualHeat')
    call diag_print_message(' Temperature: iteration     residual')
    rmom(5)=100.0
    iterheat=0
    tempheat = heat
    do iterheat=1,itermaxheat
      !Reset source term without deferred corrections
      su=suheat0
      !Compute deferred corrections
      if(ncellsimple==ncells)then !No gradients required
        select case(ndsch) !Anti-diffusion corrections
        case(5); call defcorhlpa(heat,su)
        case(6); call defcorgamma(gammadefcor,heat,su)
        case(7); call defcorgamma(cubistadefcor,heat,su)
        case(8); call defcorgamma(alvsmartdefcor,heat,su)
        case(9); call defcorgamma(hoabdefcor,heat,su)
        end select
      else 
        call der_grad_eval(gow,nlim,heat,dheatx,dheaty)
        select case(ndsch) !Anti-diffusion corrections
        case(5); call defcorhlpagrad(heat,dheatx,dheaty,su)
        case(6); call defcorgammagrad(gammadefcor,heat,dheatx,dheaty,su)
        case(7); call defcorgammagrad(cubistadefcor,heat,dheatx,dheaty,su)
        case(8); call defcorgammagrad(alvsmartdefcor,heat,dheatx,dheaty,su)
        case(9); call defcorgammagrad(hoabdefcor,heat,dheatx,dheaty,su)
        case default; if(skewcor) call defcorparagrad(dheatx,dheaty,su) !Skewness correction 
        end select
      endif
      !Solve matrix for heat
      call solve(acoef,su,sp,rsheat,heat,5)
      !Screen Output
      if(iterheat==1 .or. mod(iterheat,10)==0)then
        write(*,555) iterheat,rmom(5)
      endif  
      !Check for convergence
      if((rmom(5)<tolheat) .and. iterheat>5)then
        exit  
      endif
    enddo
    
!$OMP PARALLEL DO PRIVATE (i)
    do i=1,ncellsD
      if(iwet(i)==0) heat(i) = tempheat(i)
    enddo
!$OMP END PARALLEL DO

    
    if(mod(iterheat,10)/=0) write(*,555) iterheat,rmom(5)
    open(dgunit,file=dgfile,access='append')
    write(dgunit,555) iterheat,rmom(5)
    close(dgunit)                

!    call bndcopy2ghost(heat) !Copy heat to dummy/ghost cells

555 format(4x,I10,1x,E13.4,1x)
    
    call heat_step_stat
    
    return
    end subroutine heat_imp
    
!**************************************************************
    subroutine coeffinitheat
! Calculates the coefficients for heat initial conditions
! which is a Laplace equation
! written by Alex Sanchez, USACE-CHL
!**************************************************************    
    use size_def
    use geo_def, only: ncface,dsxy
    use flow_def, only: acoef,su,sp
    use heat_def
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
    end subroutine coeffinitheat

!***********************************************************************
    subroutine boundinit_heat
! Applies boundary condition to heat transfer equation for initial condition
! Alex Sanchez, USACE-CHL; Weiming Wu, NCCHE
!***********************************************************************
    use size_def
    use geo_def
    use flow_def
    use struct_def
    use comvarbl
    use heat_def
    implicit none
    integer :: i,j,k,iht

!   Dry nodes
    do i=1,ncells
      if(iwet(i)==0) then !Dry
        heat(i)=0.0
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
  
 !  Heat boundary
    do iht=1,nheatstr
      do j=1,heat_str(iht)%ncells
        i=heat_str(iht)%cells(j)
        if(i==0) cycle !Land cell
        acoef(:,i)=0.0
        sp(i)=-1.0
        su(i)=heat_str(iht)%heatbnd                     
      enddo
    enddo   

!   Observation Stations
    do iht=1,nheatobs
      i=heatobspts(iht,4)
      acoef(:,i)=0.0
      sp(i)=-1.0
      su(i)=heatobspts(iht,3)
    enddo
    
    return
    end subroutine boundinit_heat
    
!*************************************************************************
    subroutine coeffsourcesink_heat(schmcoef)
! Assembles the coefficient matrix and R.H.S. for heat transfer
! written by Alex Sanchez, USACE-CHL
!*************************************************************************
    use size_def
    use geo_def, only: areap,ncface,cell2cell,dsxy,idirface
    use flow_def, only: h1,h2,hk,iwet,acoef,sp,flux,visk,viskfl
    use comvarbl, only: dtime,ntsch,ctsch1,ctsch2
    use heat_def
    use comp_lib
    use prec_def
    implicit none
    integer :: i,k
    real(ikind):: ddk,val,fac1,fac2,dtimeinv,schmidtheatinv
    
    interface
      function schmcoef(dk,fk)
        use prec_def
        implicit none
        real(ikind),intent(in) :: dk,fk
        real(ikind) :: schmcoef
      end function
    endinterface    

    dtimeinv=1.0/dtime
    schmidtheatinv=1.0/schmidtheat   !Warning:  The 'schmidtheat' variable has not been initialized.   MEB - 062016
    if(ntsch==1)then
!$OMP PARALLEL DO PRIVATE(i,k,ddk)
      do i=1,ncells
        fac1=areap(i)*h1(i)*dtimeinv
        sp(i)=-fac1
        suheat0(i)=fac1*heat1(i)
        do k=1,ncface(i)          
          ddk=viskfl(k,i)*hk(k,i)*dsxy(k,i)*schmidtheatinv
          acoef(k,i)=schmcoef(ddk,flux(k,i))
        enddo
      enddo   
!$OMP END PARALLEL DO
    else
!$OMP PARALLEL DO PRIVATE(i,k,val,fac1,fac2,ddk)    
      do i=1,ncells
        val=areap(i)*dtimeinv
        fac1=val*ctsch1*h1(i)
        fac2=val*ctsch2*h2(i)
        sp(i)=-fac1+fac2
        suheat0(i)=fac1*heat1(i)-fac2*heat2(i)
        do k=1,ncface(i)          
          ddk=viskfl(k,i)*hk(k,i)*dsxy(k,i)*schmidtheatinv
          acoef(k,i)=schmcoef(ddk,flux(k,i))
        enddo 
      enddo
!$OMP END PARALLEL DO
    endif
    
    return
    end subroutine coeffsourcesink_heat
    
!******************************************************************
    subroutine heat_step_stat
! Calculates the heat transfer statistics for each time step.
! written by Alex Sanchez, USACE-CHL
!******************************************************************   
    use size_def
    use geo_def, only: mapid
    use flow_def, only: iwet
    use heat_def
    use diag_def, only: dgfile,dgunit
    use prec_def
    implicit none
    integer :: i,idheatmax,idheatmax2,idheatmin,idheatmin2,nwet
    real(ikind) :: heatmax,heatmax2,heatmin,heatmin2,heatavg
    
    !Initialize
    heatmax=-1.0e6; heatmax2=-1.0e6; idheatmax=1; idheatmax2=1
    heatmin= 1.0e6; heatmin2= 1.0e6; idheatmin=1; idheatmin2=1
    heatavg=0.0; nwet = 0
    
!!$OMP PARALLEL PRIVATE(idheatmax2,heatmax2,idheatmin2,heatmin2)
!!$OMP DO PRIVATE(i) REDUCTION(+:heatavg) REDUCTION(+:nwet)
    do i=1,ncells
      if(iwet(i)==1)then
        nwet=nwet+1
        heatavg=heatavg+heat(i)        
        if(heat(i)>heatmax2)then
          heatmax2=heat(i)
          idheatmax2=i
        endif
        if(heat(i)<heatmin2)then
          heatmin2=heat(i)
          idheatmin2=i
        endif
      endif
    enddo
!!$OMP END DO
!!$OMP CRITICAL
    if(heatmax2>heatmax)then
      heatmax=heatmax2
      idheatmax=idheatmax2
    endif
    if(heatmin2<heatmin)then
      heatmin=heatmin2
      idheatmin=idheatmin2
    endif
!!$OMP END CRITICAL
!!$OMP END PARALLEL

    nwet = max(nwet,1)
    heatavg = heatavg/real(nwet,kind=ikind)
    !write(*,626) heatavg,idheatmax,heatmax,idheatmin,heatmin    
    idheatmax = mapid(idheatmax)
    idheatmin = mapid(idheatmin)
    
626 format('   Avg Temp=',F0.3,', Max Temp(',I6,')=',F0.3,', Min Temp(',I6,')=',F0.3)
    !if(abs(Ctmax)>1.0e5 .or. abs(dzbxtr)>1.0e3) return    
    open(dgunit,file=dgfile,access='append')
    write(*,626)      heatavg,idheatmax,heatmax,idheatmin,heatmin
    write(dgunit,626) heatavg,idheatmax,heatmax,idheatmin,heatmin
    close(dgunit)      
    
    return
    end subroutine heat_step_stat    
    