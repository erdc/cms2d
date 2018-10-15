!==================================================================
! CMS Hot-start routines
!
! Contains the following:
!   hot_default - Sets the hotstart module default settings
!   hot_cards   - Reads the hotstart cards from the control file
!   hot_init    - Initializes the hotstart module
!   hot_print   - Prints the hotstart settings to the diagnostic
!                 file and the screen
!
! written by Weiming, Wu, NCCHE,
!            Alex Sanchez, USACE-CHL
!            Mitchel Brown, USACE-CHL
!==================================================================    

!******************************************************************    
    subroutine hot_default
! Sets hot start default valus
! writtin by Alex Sanchez, USACE-CHL    
!******************************************************************
    use hot_def
    use comvarbl, only: flowpath
    implicit none   
    integer :: npath
    
    coldstart =  .true.
    hot_timehr = .false.
    hot_recur =  .false.
    hot_out   =  .false.
    setconc2eq = .false.
    hotdt = 0.0    !hrs
    hottime = 0.0  !hrs     
    timeout = -1.0 !hrs
    ahotfile = 'Auto_Hot_Start.h5'
!!    hotpath = 'Datasets/'//simlabel
    hotpath = 'Datasets/'               !'Dataset'  - had to change because the path wasn't written correctly
    icfile  = 'Initial_Cond_File.h5'
    icpath  = 'Datasets/'               !'Dataset'  - had to change because the path wasn't written correctly
    ictime  = -999.0  !hrs
    npath   = len_trim(flowpath)
    hotfile = flowpath(1:npath) // 'Hot_Start.h5'      
    icpres = .false.
    icwse  = .false.
    icvel  = .false.
    icwet  = .false.
    icflux = .false.
    
    return
    endsubroutine hot_default
    
!*************************************************************   
    subroutine hot_cards(cardname,foundcard)
! Reads hot start data from Model Parameters file
! written by Alex Sanchez, USACE-CHL
!*************************************************************     
    use comvarbl, only: flowpath,tjulday0
    use time_lib, only: calendar2julian,julianday2calendarmonthday
    use hot_def
    use prec_def
    implicit none
    integer :: iyrhot,imohot,idayhot,ihrhot,iminhot,isechot,ierr
    real(ikind) :: tjuldayhot
    character(len=37) :: cardname    
    character :: apath*200,aname*200,aext*10,atemp*200
    logical :: foundcard
    
    foundcard = .true.
    selectcase(cardname)
      !----- Initial Condition (Input) ----------------------------------------
      case('INITIAL_STARTUP_FILE','INITIAL_CONDITION_FILE')
        call card_dataset(77,icfile,flowpath,icfile,icpath)
        atemp = icfile
        call uppercase(atemp) 
        if(atemp(1:7)=='DEFAULT' .or. atemp(1:4)=='NONE')then
          icfile = 'NONE'
          coldstart = .true.
        else
          coldstart = .false.       
          call fileparts(icfile,apath,aname,aext)
          if(len_trim(apath)==0)then
            icfile = trim(flowpath) // icfile
          endif
        endif
        
      case('INITIAL_CONDITION_TIME')
        call card_scalar(77,'hrs','hrs',ictime,ierr)
        coldstart = .false.
        
      case('INITIAL_CONDITION_DATE_TIME')
        call card_datetime(88,iyrhot,imohot,idayhot,ihrhot,iminhot,isechot) !YYYY-MM-DD HH:MM:SS UTC
        call calendar2julian(iyrhot,imohot,idayhot,ihrhot,iminhot,isechot,tjuldayhot)
        ictime = (tjuldayhot - tjulday0)*24 !hours
        
      !--- Hot Start (Output) -------------------------------------------  
      case('HOT_START_OUTPUT_FILE','HOT_START_FILE')
        backspace(77)
        read(77,*) cardname, hotfile   
        call fileparts(hotfile,apath,aname,aext)
    	if(len_trim(apath)==0)then
	      hotfile = trim(flowpath) // hotfile
	    endif       
        
      case('HOT_START_TIME')
        call card_scalar(77,'hrs','hrs',hottime,ierr)
        hot_timehr = .true.        
        hot_out = .true.  

      case('AUTO_HOT_START_INTERVAL','HOT_START_INTERVAL','HOT_START_OUTPUT_INTERVAL')
        backspace(77)
        read(77,*) cardname, hotdt  	
        hot_recur = .true.	  
        hot_out = .true.   
        
      case default
        foundcard = .false.           
        
    endselect
    
    return
    endsubroutine hot_cards

!******************************************************************    
    subroutine hot_read()
! Reads hot start information
! written by Alex Sanchez, USACE-CHL    
! modified by Mitch Brown - 05/18/2012
!******************************************************************
#include "CMS_cpp.h"  
    use flow_def, only: eta,u,v
    use hot_def
    use comvarbl, only: timehrs,stime,ctime
    use diag_lib
    implicit none
    character(len=10) :: aext
    logical :: ok
    
    !Check if hot_start file exists
    inquire(file=icfile,exist=ok)
    if(.not.ok)then
      call diag_print_error('Unable to find Initial Condition File: ',icfile)
    endif
    
    call fileext(icfile,aext)
    if(aext(1:2)=='h5')then
#ifdef XMDF_IO
      call hot_read_xmdf
#else
      call diag_print_error('Cannot read hot start from *.h5 file without XMDF libraries')
#endif
    else
      call hot_read_dat !Not implemented yet
    endif   

    return
    endsubroutine hot_read
    
!********************************************************************
    subroutine hot_read_dat()
! Reads initial condition file from an SMS ASCII dataset (*.dat) file 
! written by Alex Sanchez, USACE-CHL    
!********************************************************************    
    use size_def
    use flow_def, only: eta,u,v,h,p,iwet,grav
    use hot_def
    use in_def
    use in_lib, only: read_dat
    use comvarbl, only: timehrs
    use interp_lib, only: interp_scal_node2cell
    use unitconv_lib, only: unitconv_scal
    use diag_lib
    implicit none
    integer :: nscal,nvec
    type(scaldattype), pointer :: scaldat(:)
    type(vecdattype),  pointer :: vecdat(:)
    real(4) :: etemp(ncellsfull),utemp(ncellsfull),vtemp(ncellsfull),wtemp(ncellsfull)
    real(ikind) :: wet(ncellsD)
    integer :: i
    
    call read_dat(icfile,nscal,scaldat,nvec,vecdat)
    
    !--- Water Levels -----
    do i=1,nscal
      selectcase(scaldat(i)%name)
      case('Water_Elevation','Water_Surface_Elevation','Water_Level','wse','WSE')
        if(scaldat(i)%nd/=ncellsfull)then
          call diag_print_error('Invalid initial conditions size',&
            '  Size of dataset does not match grid.')
        endif
        !Water level
        etemp(:) = scaldat(i)%val(:,scaldat(i)%nt)        
        if(ncellpoly>0)then
          call interp_scal_node2cell(etemp,eta) !Interpolate node to cell centers
        else
          call map_scal_full2active(etemp,eta) !Convert from full to active grid 
        endif
        !Wet/dry
        if(scaldat(i)%istatus==1)then
          wtemp(:) = scaldat(i)%stat(:,scaldat(i)%nt)  
          if(ncellpoly>0)then
            call interp_scal_node2cell(wtemp,wet) !Interpolate node to cell centers
          else
            call map_scal_full2active(wtemp,wet) !Convert from full to active grid 
          endif
          iwet = wet
        endif
        p = eta*grav
        timehrs = scaldat(i)%time(scaldat(i)%nt)  !Start time
        call unitconv_scal(scaldat(i)%time_units,'hrs',timehrs)
        exit
      endselect
    enddo
    
    !--- Current Velocities -----
    do i=1,nvec
      selectcase(vecdat(i)%name)
      case('Current_Velocity','Water_Velocity','Depth-Averaged_Current_Velocity','Depth-Averaged_Velocity')
        if(vecdat(i)%nd/=ncellsfull)then
          call diag_print_error('Invalid initial conditions size',&
            '  Size of dataset does not match grid')
        endif
        !Velocities
        utemp(:) = vecdat(i)%val(:,1,vecdat(i)%nt)
        vtemp(:) = vecdat(i)%val(:,2,vecdat(i)%nt)        
        if(ncellpoly>0)then
          call interp_scal_node2cell(utemp,u) !Interpolate node to cell centers
          call interp_scal_node2cell(vtemp,v) !Interpolate node to cell centers
        else
          call map_scal_full2active(utemp,u) !Convert from full to active grid 
          call map_scal_full2active(vtemp,v) !Convert from full to active grid 
        endif
        !Wet/dry
        if(vecdat(i)%istatus==1)then
          wtemp(:) = vecdat(i)%stat(:,vecdat(i)%nt)
          if(ncellpoly>0)then
            call interp_scal_node2cell(wtemp,wet) !Interpolate node to cell centers
          else
            call map_scal_full2active(wtemp,wet) !Convert from full to active grid 
          endif
          iwet = wet          
        endif
        u = u*iwet
        v = v*iwet
        timehrs = vecdat(i)%time(vecdat(i)%nt)     !Start time
        call unitconv_scal(vecdat(i)%time_units,'hrs',timehrs)
        exit
      endselect
    enddo
    
    return
    endsubroutine hot_read_dat

!******************************************************************    
    subroutine hot_read_xmdf()
! Reads hot start information
! written by Alex Sanchez, USACE-CHL    
! modified by Mitch Brown - 05/18/2012
!******************************************************************
#include "CMS_cpp.h"
#ifdef XMDF_IO
    use size_def
    use geo_def, only: cell2cell,zb,zb0
    use flow_def, only: u,v,u1,v1,uv,p,p1,eta,h,iwet,flux,grav,gravinv
    use sed_def, only: sedtrans,nsed,Ctk,Ctk1,Ctkstar,icapac,nlay,bedlay,&
        pbk,pbk1,db,db1,d50,d90,diam,diamlim,logdiamlim,dbmax,dmconst,zb1
    use sal_def, only: saltrans,sal
    use heat_def, only: heattrans, heat
    use hot_def
    use in_lib, only: readscallasth5,readveclasth5
    use out_def, only: outlist,simlabel
    use xmdf
    !use ifport, only: system
    use DFLIB, only: systemqq
    use diag_lib
    use prec_def
    implicit none
    integer :: res !Added MEB - 05/18/2012 for hot start fixes
    integer :: i,j,k,ks,fid,gid,nn,ierr,ntimes,npath
    real(ikind) :: temphr,temp(ncellsD)
    real(8), allocatable :: timesd(:)
    real(ikind), allocatable :: d35(:)
    !real(ikind) :: var(ncellsD)
    character(len=200) :: apath,thepath,thename,tfile !Added tfile MEB - 05/18/2012 for hot start fixes
    character(len=5) :: apbk,alay
    character(len=10) :: theext 
    logical :: ok
    character(len=30) :: wse_names(5)
    character(len=300) :: msg
    
    tfile=trim(icfile)   
    call fileparts(icfile,thepath,thename,theext)
    
!****************** MEB - 05/18/2012
! XMDF error keeps hot_start.h5 from being deleted after reading

    ! Check if temp.h5 file exists, if so, delete it
    tfile=trim(thepath)//'temp.h5'
    inquire(file=tfile, exist=ok)
    if(ok)then
      open(100,file=tfile,iostat=ierr)
      close(100,status='delete',iostat=ierr)
    endif

    ! Make a copy of hot_start.h5 to temp.h5 and read from that.
    !res = system('copy '//trim(icfile)//' '//trim(tfile)//' > file')
    res = systemqq('copy '//trim(icfile)//' '//trim(tfile)//' > file')
    inquire(file='file', exist=ok)
    if(ok)then            !remove scratch file
      open(100,file='file',iostat=ierr)
      close(100,status='delete',iostat=ierr)
    endif
!******************     

    !Write to Diag file
    call diag_print_message(' ','*** Reading Hot Start File ***')
    
    npath = len_trim(icpath)
    ntimes = 0      !initialize variable
    temphr = -1.0

747 format(1x,A,A50)
        
!---- Water Elevation --------------------------------
    wse_names(1) = 'Water_Elevation'
    wse_names(2) = 'Water_Surface_Elevation'
    wse_names(3) = 'Water_Level'
    wse_names(4) = 'eta'
    wse_names(5) = 'wse'
    do i=1,5
      apath = icpath(1:npath) // trim(wse_names(i))
      call readscallasth5(tfile,apath,ntimes,eta,reftimehot,ictime,ierr) !use tfile instead of icfile
      if(ierr>=0) exit
    enddo
    if(ierr>=0)then
      icwse = .true.  
      write(msg,'(A,A)') 'Read water level: ',trim(apath)
      call diag_print_message(msg)
      p = eta*grav
      ictime = temphr !Start time in hours 
    endif
    if(ntimes>1)then
      write(msg,'(A,A)') 'Found multiple time steps in: ',trim(icfile)
      call diag_print_message(msg,'  Using last time step')
    endif    
    
!---- Water Pressure ----------------------------------
    apath = icpath(1:npath) //'Water_Pressure'  
    call readscallasth5(tfile,apath,ntimes,p,reftimehot,temphr,ierr)
    if(ierr>=0)then
      icpres = .true.
      write(msg,'(A,A)') 'Read water pressure: ',trim(apath)
      call diag_print_message(msg,'   Water level calculated from water pressure')
      if(.not.icwse)then
        eta = p*gravinv
      endif
      ictime = temphr !Start time in hours      
    endif
    
!---- Current Velocity -------------------------
    apath = icpath(1:npath) //'Current_Velocity'  
    call readveclasth5(tfile,apath,ntimes,u,v,reftimehot,temphr,ierr)
    if(ierr<0)then
      apath = icpath(1:npath) // 'u'  
      call readscallasth5(tfile,apath,ntimes,u,reftimehot,temphr,ierr)
      if(ierr==0)then
        apath = icpath(1:npath) // 'v'  
        call readscallasth5(tfile,apath,ntimes,v,reftimehot,temphr,ierr)
      endif
    endif
    if(ierr<0)then
      icvel = .false.  
      call diag_print_warning('Unable to find current velocities ',&
        '    Use one of the following paths and names: ',&
        '      "Datasets\Current_Velocity" ',&
        '      "Datasets\u" ',&
        '      "Datasets\v" ',&  
        '   Setting current velocities to zero')
      u = 0.0      
      v = 0.0
    else
      icvel = .true.   
      write(msg,'(A,A)') 'Read current velocity: ',trim(apath)
      call diag_print_message(msg)
      ictime = temphr !Start time in hours  
    endif
    
!---- Wet/dry --------------------------------------------------------    
    apath = icpath(1:npath) // 'Wet'  
    call readscallasth5(tfile,apath,ntimes,temp,reftimehot,temphr,ierr)
    if(ierr>=0)then
      icwet = .true.
      write(msg,*) '  Read wet/dry states: ',apath  
      call diag_print_message(msg)  
      iwet = temp
    endif
    
!--- Fluxes ----------------------------------------------------------
    do k=1,nmaxfaces
      if(k<=9)then
         write(apath,'(A,A4,I1)') icpath(1:npath),'Flux',k
      else
         write(apath,'(A,A4,I1)') icpath(1:npath),'Flux',k
      endif
      call readscallasth5(tfile,apath,ntimes,temp,reftimehot,temphr,ierr)      
      if(ierr<0)then
        icflux = .false.
        exit
      else
        icflux = .true.  
        write(msg,*) '  Read cell-face fluxes: ',trim(apath)
        call diag_print_message(msg)  
        do i=1,ncellsD
          flux(k,i) = temp(i)
        enddo        
      endif
    enddo    
      
!---- Sediment Transport ------------------------
    if(sedtrans)then !~~~~~~~ Single Grain size ~~~~~~~~~~~~~~~
      !Water depths
      apath = icpath(1:npath) // 'Depth'  
      call readscallasth5(tfile,apath,ntimes,zb,reftimehot,temphr,ierr)    
      if(ierr<0)then
        call diag_print_warning('Unable to find depths ',&
          '   Using grid depths as depths')
      else
        write(msg,'(A,A)') 'Read initial water depths: ',trim(apath)
        call diag_print_message(msg)
        zb = -zb
        zb1 = zb        
      endif
      
!      !Initial Water depths
!      apath = string(1:nn) // 'Initial_Depth'  
!      call readscalsteph5(icfile,apath,1,temphr,zb0,ierr)    
!      if(ierr<0)then
!        open(dgunit,file=dgfile,access='append') 
!        write(*,*)      'WARNING: Unable to find initial depths'
!        write(dgunit,*) 'WARNING: Unable to find initial depths'
!        write(*,*)      '   Using grid depths as initial depths'
!        write(dgunit,*) '   Using grid depths as initial depths'    
!        close(dgunit)
!        zb0 = zb
!      else
!        open(dgunit,file=dgfile,access='append') 
!        write(*,747)       '  Read initial water depths: ',apath
!        write(dgunit,747)  '  Read initial water depths: ',apath
!        close(dgunit)
!        zb0 = -zb0
!      endif
      
      !Sediment concentrations
      if(nsed==1)then !Single grain size
        apath = icpath(1:npath) //'Concentration'  
        call readscallasth5(tfile,apath,ntimes,Ctk(:,1),reftimehot,temphr,ierr)
        if(ierr<0)then
          call diag_print_warning('Unable to find initial sediment concentrations ',&
            '   Setting initial sediment concentrations to equilibrium concentration')
          setconc2eq = .true.      
        else
          write(msg,'(A,A)') 'Read sediment concentrations: ',trim(apath)
          call diag_print_message(msg)
        endif  
        
      else  !~~~~~~~~~ Multiple grain sizes ~~~~~~~~~~~~~~~~~~
        
62  format('_',I2.2)
71  format('(',I1,')')
72  format('(',I2,')')               
        
        !Layer thickness
loop1: do j=1,nlay
          if(j<=9)then
            write(alay,71) j
          else
            write(alay,72) j
          endif
          apath = icpath(1:npath) // 'Thickness' // alay
          call readscallasth5(tfile,apath,ntimes,db(:,j),reftimehot,temphr,ierr)
          if(ierr<0)then 
            apath = icpath(1:npath) // 'Thickness ' // alay
            call readscallasth5(tfile,apath,ntimes,db(:,j),reftimehot,temphr,ierr)
            if(ierr<0)then 
              apath = icpath(1:npath) // 'Thickness_' // alay
              call readscallasth5(tfile,apath,ntimes,db(:,j),reftimehot,temphr,ierr)
            endif            
          endif
          if(ierr<0)then  
            call diag_print_warning('Unable to find initial bed layer thickness',&
               '  Setting thickness of layers to maximum thickness except for mixing layer thickness')
            db = dbmax
            db(:,2) = db(:,1)-dmconst
            db(:,1) = dmconst   
            exit
          else
            write(msg,'(A,A)') 'Read bed layer thickness: ',trim(apath)
            call diag_print_message(msg)
          endif
        enddo loop1         
        
        !---- Fractions ----------
        ok = .true.
loopj:  do j=1,nlay
          if(j<=9)then
            write(alay,71) j
          else
            write(alay,72) j
          endif  
          do ks=1,nsed
            write(apbk,62) ks
            apath = icpath(1:npath) //'Fraction' // trim(apbk) // alay 
            call readscallasth5(tfile,apath,ntimes,pbk(:,ks,j),reftimehot,temphr,ierr)            
            if(ierr<0)then
              call diag_print_warning('Problem reading bed composition',&
               '  Setting the bed composition to the default values')  
              ok = .false.
              exit loopj
            else
              write(msg,'(A,A)') 'Read bed layer thickness: ',trim(apath)
              call diag_print_message(msg)  
            endif
          enddo !ks
        enddo loopj !j         
        
        !Look for D35,D50,D90   
        if(.not.ok)then
          ok = .true.  
          apath = icpath(1:npath) // 'D35'
          if(.not.allocated(d35)) allocate(d35(ncellsD))
          call readscallasth5(tfile,apath,ntimes,d35,reftimehot,temphr,ierr)          
          if(ierr==0)then
            apath = icpath(1:npath) // 'D90'
            call readscallasth5(tfile,apath,ntimes,d90,reftimehot,temphr,ierr)              
            if(ierr<0)then
              call diag_print_warning('Unable to find initial D90 dataset')
            else
              write(msg,'(A,A)') 'Read D90 dataset: ',trim(apath)
              call diag_print_message(msg)
            endif
            apath = icpath(1:npath) // 'D50'
            call readscallasth5(tfile,apath,ntimes,d50,reftimehot,temphr,ierr)  
            call diag_print_message('  Read initial D35, D50 and D90 datasets ',&
              '    Calculating bed composition based on D35, D50, and D90 datasets')
            diam = diam*1000.0      
            diamlim = diamlim*1000.0
            call bed_d35d50d90(nsed,diam,diamlim,d35,d50,d90,pbk(:,:,1)) !Calculates grain size distribution from D35,D50,D90 
            do j=2,nlay
              pbk(:,:,j) = pbk(:,:,1)  
            enddo            
            deallocate(d35)    
            d50 = d50/1000.0    
            d90 = d90/1000.0
            diam = diam/1000.0
            diamlim = diamlim/1000.0
          else
            call diag_print_warning('Unable to find initial D35 dataset')
            apath = icpath(1:npath) // 'D50'
            call readscallasth5(tfile,apath,ntimes,d50,reftimehot,temphr,ierr)                
            if(ierr<0)then
              call diag_print_warning('Unable to find initial D50 dataset')
            else
              call diag_print_message('Read D50 dataset ',&
                '   Calculating bed composition based on D50 dataset and standard deviation')
              diam = diam*1000.0
              diamlim = diamlim*1000.0      
              call bed_d50sigma(nsed,diam,diamlim,d50,bedlay(1)%geostddev,pbk(:,:,1)) !Calculates grain size fractions based on D50 and standard deviation
              do j=1,nlay
                pbk(:,:,j) = pbk(:,:,1)  
              enddo
              d50 = D50/1000.0    
              d90 = D90/1000.0
              diam = diam/1000.0
              diamlim = diamlim/1000.0
            endif
          endif     
        endif
        
        !----- Concentrations -------------------
        apath = icpath(1:npath) // 'Concentration'  
        call readscallasth5(tfile,apath,ntimes,Ctk(:,1),reftimehot,temphr,ierr)
        if(ierr<0)then
!          open(dgunit,file=dgfile,access='append') 
!          write(*,*)      'WARNING: Unable to find initial sediment concentrations'
!          write(dgunit,*) 'WARNING: Unable to find initial sediment concentrations'
!          write(*,*)      '   Setting initial sediment concentrations to equilibrium concentration'
!          write(dgunit,*) '   Setting initial sediment concentrations to equilibrium concentration'
!          close(dgunit)     
          do ks=1,nsed
            write(apbk,62) ks
            apath = icpath(1:npath) // 'Concentration_' // apbk
            call readscallasth5(tfile,apath,ntimes,Ctk(:,ks),reftimehot,temphr,ierr)
            if(ierr<0)then
              call diag_print_warning('Unable to find initial sediment concentrations',&
                '   Setting initial sediment concentrations to equilibrium concentration')
              setconc2eq = .true. 
              exit
            endif 
          enddo
        else
          write(msg,*) 'Read sediment concentration: ',trim(apath)
          call diag_print_message(msg)
          do ks=1,nsed
            Ctkstar(:,ks) = Ctk(:,1)*pbk(:,ks,1)          
          enddo  
          Ctk = Ctkstar
        endif   
        
        logdiamlim=log(diamlim)
        !Calculate necessary sediment percentiles for calculations
        !Others are calculated only as needed for output
        call sedpercentile(50,d50)
        call sedpercentile(90,d90)
        
        Ctk1=Ctk
        pbk1=pbk(:,:,1)
        db1=db
      endif
    endif
    
!---- Salinity Transport ------------------------
    if(saltrans)then
      apath = icpath(1:npath) // 'Salinity'
      call readscallasth5(tfile,apath,ntimes,sal,reftimehot,temphr,ierr)
      if(ierr<0)then
        call diag_print_warning('Unable to find initial salinity concentrations',&
           '   Setting salinity concentrations to zero')
        sal = 0.0
      else
        write(msg,*) 'Read salinity concentrations: ',trim(apath)
        call diag_print_message(msg)
      endif
    endif
    
!---- Heat Transfer ------------------------
    if(heattrans)then
      apath = icpath(1:npath) // 'Temperature'
      call readscallasth5(tfile,apath,ntimes,sal,reftimehot,temphr,ierr)
      if(ierr<0)then
        call diag_print_warning('Unable to find initial temperatures',&
           '   Setting temperatures to 15 deg')
        heat = 15.0
      else
        write(msg,*) 'Read temperatures: ',trim(apath)
        call diag_print_message(msg)
      endif
    endif
    
!--- Get Last Output Time ----------------------------------
755 format(1x,A,F12.3,A)  
    timeout = -1.0
    call XF_OPEN_FILE(trim(outlist(1)%afile),READONLY,fid,ierr)             
    if(ierr<0)then
      call diag_print_message(' Solution file is missing',&
        '  Will create new solution file')
    else
      nn=len_trim(simlabel)
      apath = simlabel(1:nn)//'/Water_Elevation'
      call XF_OPEN_GROUP(fid,trim(apath),gid,ierr)  
      if(ierr<0)then
        apath = simlabel(1:nn)//'/Current_Velocity'
        call XF_OPEN_GROUP(fid,trim(apath),gid,ierr)  
      endif
      call XF_GET_DATASET_NUM_TIMES(gid,ntimes,ierr)
      allocate(timesd(ntimes))
      call XF_GET_DATASET_TIMES(gid,ntimes,timesd,ierr)
      timeout = timesd(ntimes) !Hours
      deallocate(timesd)
!    call XF_CHANGE_SCALAR_VALUES_TIMESTEP_FLOAT(FileId, TimestepIndex, NumValsToEdit, &
!                                                  Indices, NewValues, ierr)    
      call XF_CLOSE_GROUP(gid,ierr)  !Close dataset    
      call XF_CLOSE_FILE(fid,ierr)   !Close XMDF file        
    endif
    
    !if(.not.fluxreadin)then       
    !  call diag_print_message('Determining fluxes from water levels and current velocities')
    !  call interp_scal_cell2face(p,1,pk,dpx,dpy)
    !  hk = max(hmin,pk*gravinv-zbk)
    !  where(hk<hmin) hk = hmin
    !  call flow_convflux
    !endif
    
!--- Delete Temporary file ---------------------
    call delete_file(tfile)

#endif    
    return
    endsubroutine hot_read_xmdf
    
!******************************************************************    
    subroutine hot_init()
! Initializes any missing variables icluding Ctk and pbk
! and takes care of boundaries
! written by Alex Sanchez, USACE-CHL
!
! only called when hot start is to be read   !MEB  02/24/2016
!******************************************************************    
#include "CMS_cpp.h"
    use size_def
    use geo_def, only: cell2cell,zb
    use flow_def, only: u,v,u1,v1,uv,p,p1,eta,h,h1,hmin,iwet,iwet1,ponding,grav,gravinv
    use comvarbl, only: ctime,stime,timehrs,timesecs
    use sed_def, only: nsed,icapac,Ctk,Ctk1,Ctkstar,CtstarP,pbk,&
       singlesize,variableD50,d50,diam,mhe,iHidExpForm,varsigma
    use stat_def, only: calc_stats,tstat
    use bnd_def
    use hot_def
    use out_def, only: simlabel
    use sal_def
    use heat_def
    use diag_lib
    implicit none
    integer :: i,j,k,ibnd,nck,iwse,isal,iheat
    character(len=100) :: msg2,msg3
            
    !if(.not.icwse .and. icpres)then !Only pressure specified
    !  eta = p*gravinv
    !elseif(icwse .and. .not.icpres)then !Only water level specified
    !  p = eta*grav      
    !elseif(.not.icwse .and. .not.icpres)then
    !  call  diag_print_warning('Unable to find initial water level/pressure',&
    !    '    Use one of the following paths and names: ',&
    !    '      "Datasets\Water_Level" ',&
    !    '      "Datasets\Water_Elevation" ',&
    !    '      "Datasets\Water_Surface_Elevation"',&
    !    !'      "Datasets\Water_Pressure"',&
    !    '    Setting water level to zero')
    !  eta = 0.0   
    !  p = 0.0      
    !endif   
    !  
    !if(.not.icvel)then
    !  call diag_print_warning('Unable to find current velocities ',&
    !    '    Use one of the following paths and names: ',&
    !    '      "Datasets\Current_Velocity" ',&
    !    '      "Datasets\u" ',&
    !    '      "Datasets\v" ',&  
    !    '   Setting current velocities to zero')
    !  u = 0.0      
    !  v = 0.0
    !else
    !  write(msg,11) 'Read current velocity: ',trim(apath)
    !  call diag_print_message(msg)
    !  ictime = temphr !Start time in hours  
    !endif  
    
    !--- Hot Start Time ---------------------------------------------
755 format(1x,A,F10.3,A)    
    !Overwrites time in icfile. The time in icfile may be incorrectly set to zero by SMS
    if(ictime>-1.0e-20)then !if specified
      timehrs = ictime !Initial condition time
      !call diag_print_message('Initial Condition Time Specified...')
    else 
      if(timeout>-1.0e-20)then !if specified
        timehrs = timeout !Used time in solution file
        !call diag_print_message('Using last output time as hot start time...')
      else  
        timehrs = 0.0  !no initial condition time estimate available so set to default
        !call diag_print_warning('Initial Condition Time set to zero...')        
      endif  
    endif
    stime = timehrs*3600.0
    ctime = stime
    timesecs = dble(ctime)
    hstarttime = timehrs
!!    nspinup=max(nspinup,5) !********
    
    write(msg2,755) '  Hot Start Time:   ',timehrs,' hours'
    write(msg3,755) '  Last Output Time: ',timeout,' hours'
    call diag_print_message(' ',msg2,msg3)   
    
    !--- Simulation Statistics ----
    if(calc_stats)then
      if(ctime>=tstat(1) .and. ctime<=tstat(2))then         
#ifdef XMDF_IO
        call stat_read !Save for later use
#else
        call diag_print_error('Cannot read simulation statistics from *.h5 file without XMDF libraries')
#endif
      endif   
    endif
    
    !!!Wetting and drying
    !call flow_wetdry(0)
    !!If all cells are dry, then assume closed domain and turn on ponding
    !if(sum(1-iwet(1:ncells))==ncells)then 
    !  ponding=.true.
    !  call flow_wetdry(0)
    !endif 
    !iwet1=iwet    
    !eta = iwet*p*gravinv-999.0*(1-iwet)
    
    !call der_grad_eval(goa,0,zb,dzbx,dzby) !Bed-slope
    !call interp_scal_cell2face(zb,0,zbk,dzbx,dzby)
    !call flow_grad_interp
      
    !---- Apply wet/dry ------------
    !Dry nodes     
    do i=1,ncellsD
      h(i) = p(i)*gravinv - zb(i)
      if(h(i)<hmin+1.0e-5)then    
        h(i) = hmin
        iwet(i) = 0  
        u(i) = 0.0
        v(i) = 0.0  
        eta(i) = -999.0
        if(.not.icpres)then
          p(i) = (h(i)+zb(i)-1.0e-4)*grav
        endif
      else
        iwet(i) = 1
        eta(i)=p(i)*gravinv
      endif
      uv(i)=sqrt(u(i)*u(i)+v(i)*v(i))
      u1(i)=u(i)
      v1(i)=v(i)      
      h1(i)=h(i)
      p1(i)=p(i)  
      iwet1(i)=iwet(i)
    enddo
    
    !All Forcing Boundaries
    do ibnd=1,nbndstr
      do j=1,bnd_str(ibnd)%ncells
        i=bnd_str(ibnd)%cells(j)
        k=bnd_str(ibnd)%faces(j)
        nck=cell2cell(k,i)
        u(nck)=u(i)
        v(nck)=v(i)
        p(nck)=p(i)
        eta(nck)=eta(i)
        zb(nck)=zb(i)
        h(i)=max(hmin,p(i)*gravinv-zb(i))
        h(nck)=h(i)        
      enddo
    enddo

    call flow_update
    
    !-- Set concentration to equilibrium concentrations ----
    if(setconc2eq)then
      !Hiding and Exposure
      selectcase(iHidExpForm) !varsigma(i,k)
      case(1); call HidExpEgiazaroff    !Egiazaroff (1965)
      case(2); call HidExpParker        !Parker et al. (1982) and others
      case(3); call HidExpWu            !Wu et al. (2000)
      case(4); call HidExpAshidaMichiue !Ashida and Michiue 1980
      case(5); call HidExpHayashi       !Hayashi et al. 1980
      endselect
      
  	  !Transport Capacity
      selectcase(icapac)  
      case(1); call sedcapac_lundcirp !Lund-CIRP          
      case(2); call sedcapac_vanrijn  !Van Rijn           
      case(3); call sedcapac_watanabe !Watanabe  
      case(4); call sedcapac_soulsby  !Soulsby (1997)
      case(5); call wucapac           !Wu et al. 2000 (under testing)
      endselect   
!!      CtstarP = 0.0
      !Concentrations
      Ctkstar = CtstarP*pbk(:,:,1)
      Ctk = Ctkstar    
      Ctk1 = Ctk
    endif
    
!--- Save initial water levels and current velocities -----
    !Tidal/Harmonic boundary
    do iwse=1,nTHstr
      do j=1,TH_STR(iwse)%ncells
        i=TH_STR(iwse)%cells(j)
        TH_STR(iwse)%wsebnd0(j)=eta(i)
      enddo  
    enddo
    
    !Single Water Level BC
    do iwse=1,nHstr
      do j=1,H_STR(iwse)%ncells
        i=H_STR(iwse)%cells(j)
        H_STR(iwse)%wsebnd0(j)=eta(i)
      enddo
    enddo
    
    !Multiple Water Level BC
    do iwse=1,nMHstr
      do j=1,MH_STR(iwse)%ncells
        i=MH_STR(iwse)%cells(j)
        MH_STR(iwse)%wsebnd0(j)=eta(i)
      enddo
    enddo
    
    !Multiple Water Level and Velocity BC
    do iwse=1,nMHVstr
      do j=1,MHV_STR(iwse)%ncells
        i=MHV_STR(iwse)%cells(j)
        MHV_STR(iwse)%wsebnd0(j)=eta(i)
        MHV_STR(iwse)%ubnd0(j)=u(i)
        MHV_STR(iwse)%vbnd0(j)=v(i)
      enddo
    enddo
    
    !Nested Water Level BC
    do iwse=1,nNHstr
      do j=1,NH_STR(iwse)%ncells
        i=NH_STR(iwse)%cells(j)
        NH_STR(iwse)%wsebnd0(j)=eta(i)
      enddo
    enddo
    
    !Nested Water Level and Velocity BC
    do iwse=1,nNHVstr
      do j=1,NHV_STR(iwse)%ncells
        i=NHV_STR(iwse)%cells(j)
        NHV_STR(iwse)%wsebnd0(j)=eta(i)
        NHV_STR(iwse)%ubnd0(j)=u(i)
        NHV_STR(iwse)%vbnd0(j)=v(i)
      enddo
    enddo
    
    !Tidal Database Water Level BC
    do iwse=1,nNTHstr
      do j=1,NTH_STR(iwse)%ncells
        i=NTH_STR(iwse)%cells(j)
        NTH_STR(iwse)%wsebnd0(j)=eta(i)
      enddo
    enddo
    
    !Tidal Database Water Level and Velocity BC
    do iwse=1,nNTHVstr
      do j=1,NTHV_STR(iwse)%ncells
        i=NTHV_STR(iwse)%cells(j)
        NTHV_STR(iwse)%wsebnd0(j)=eta(i)
        NTHV_STR(iwse)%ubnd0(j)=u(i)
        NTHV_STR(iwse)%vbnd0(j)=v(i)
      enddo
    enddo  
    
!--- Salinity Boundary Condition ----------------------------------------------------  
    do isal=1,nsalstr
      do j=1,sal_str(isal)%ncells
        i=sal_str(isal)%cells(j)
        sal_str(isal)%salbnd0(j) = sal(i)
      enddo
    enddo     
    
!--- Temperature Boundary Condition ----------------------------------------------------  
    do iheat=1,nheatstr
      do j=1,heat_str(iheat)%ncells
        i=heat_str(iheat)%cells(j)
        heat_str(iheat)%heatbnd0(j) = heat(i)
      enddo
    enddo     

    call diag_print_message('*** Hot Start Initialization Complete ***')   
    
    return
    endsubroutine hot_init
    
!******************************************************************
    subroutine hot_write
! Writes hot start h5 file for restarting CMS
! written by Alex Sanchez, USACE-CHL
!******************************************************************
#include "CMS_cpp.h"
    use hot_def, only: hottime,hot_timehr,&
        hotdt,hot_recur,icfile,hotpath
    use comvarbl, only: timehrs
    use diag_def
    use diag_lib
    use prec_def
    implicit none
    real(ikind) :: eps
    logical :: ihot
    
    eps = 1.e-10    
    ihot = .false. 
    if(hot_timehr .and. abs(timehrs-hottime)<=eps) ihot = .true.      
    if(hot_recur .and. mod(timehrs,hotdt)<=eps) ihot = .true.   
    if(.not.ihot) return              
        
    write(msg,'(A,F13.4,A)') 'Writing Hotstart File at: ',timehrs,' hrs'
    call diag_print_message(' ',msg)
    
#ifdef XMDF_IO
    call hot_write_xmdf !XMDF hot file
#else
    !call hot_write_bin !Binary hot file, not implemented yet
    call diag_print_error('Cannot write out hot start file without XMDF libraries')    
#endif

    return
    endsubroutine hot_write

!******************************************************************
    subroutine hot_write_xmdf()
! Writes hot start h5 file for restarting CMS
! written by Alex Sanchez, USACE-CHL
!******************************************************************
#include "CMS_cpp.h"
#ifdef XMDF_IO
    use size_def
    use hot_def, only: hotfile,hottime,hot_timehr,&
        hotdt,hot_recur,icfile,hotpath
    use geo_def, only: zb
    use flow_def, only: u,v,u1,v1,p,p1,eta,iwet,iwet1,flux,flux1,grav
    use comvarbl, only: ctime,timehrs,flowpath,ntsch
    use sed_def, only: sedtrans,nsed,Ctk,nlay,pbk,db
    use sal_def, only: saltrans,sal
    use heat_def, only: heattrans, heat
    use stat_def, only: flowstats,sedstats,salstats,heatstats
    use out_def, only: simlabel
    use out_lib, only: writescalh5,writevech5
    use prec_def
!!    use ifport  !Only for intel compiler
    use xmdf
    implicit none
    integer :: j,ks
#ifdef DEV_MODE
    integer :: i,k
#endif
    real(ikind) :: var(ncellsD)
    character(len=100) :: aname
    character(len=5) :: apbk,alay
    
    !delete previous hotstart file so there is only one record 
    open(100,file=hotfile)
    close(100,status='DELETE')  
  
    !--- Hydrodynamics ----
    call writescalh5(hotfile,hotpath,'Water_Pressure',p,'m^2/s^2',timehrs,1)   
    call writescalh5(hotfile,hotpath,'Water_Elevation',eta,'m',timehrs,1)
    call writevech5(hotfile,hotpath,'Current_Velocity',u,v,'m/s',timehrs,1)
     
#ifdef DEV_MODE    
    !Wet/dry
    var = iwet
    call writescalh5(hotfile,hotpath,'Wet',var,'',timehrs,1)    
    
    !Fluxes
    do k=1,nmaxfaces
      do i=1,ncellsD
        var(i) = flux(k,i)
      enddo
      if(k<=9)then
         write(aname,'(A,I1)') 'Flux',k
      else
         write(aname,'(A,I2)') 'Flux',k
      endif 
      call writescalh5(hotfile,hotpath,aname,var,'m^3/s',timehrs,1)
    enddo
#endif

!    !Second-order temporal scheme
!    if(ntsch==2)then
!      call writescalh5(hotfile,hotpath,'Water_PressureS',p1,'m^2/s^2',timehrs,1)   
!!      call writescalh5(hotfile,apath,'Water_Elevation',eta,'m',timehrs,1)
!      call writevech5(hotfile,hotpath,'Current_VelocityS',u1,v1,'m/s',timehrs,1)
!      var = iwet1
!      call writescalh5(hotfile,hotpath,'WetS',var,'',timehrs,1)    
!      do k=1,nmaxfaces
!        var = flux1(:,k)
!        if(k<=9)then
!           write(aname,'(A,I1)') 'FluxS',k
!        else
!           write(aname,'(A,I1)') 'FluxS',k
!        endif 
!        call writescalh5(hotfile,hotpath,aname,var,'m^3/s',timehrs,1)
!      enddo  
!    endif
    
    !--- Sediment Transport ----
    if(sedtrans)then      
      !Water depths
      var = -zb
      call writescalh5(hotfile,hotpath,'Depth',var,'m',timehrs,1)  
!      !Initial Water depths
!      var = -zb0
!      call writescalh5(hotfile,apath,'Initial_Depth',var,'m',timehrs,1)  
      !Sediment concentrations
      if(nsed==1)then !Single grain size
        aname = 'Concentration'
        call writescalh5(hotfile,hotpath,aname,Ctk(:,1),'kg/m^3',timehrs,1)  
      else              !Multiple grain size

62  format('_',I2.2)
71  format(1x,'(',I1,')')
72  format(1x,'(',I2,')')        
        !Concentrations
        do ks=1,nsed
            write(apbk,62) ks
            aname = 'Concentration'// apbk
          call writescalh5(hotfile,hotpath,aname,Ctk(:,ks),'kg/m^3',timehrs,1)  
        enddo
        
        !Bed material composition and layer thickness
        do j=1,nlay
          if(j<=9)then
            write(alay,71) j
          else
            write(alay,72) j
          endif  
          aname = 'Thickness' // alay !Bug Fix, changed apath to aname
          call writescalh5(hotfile,hotpath,aname,db(:,j),'m',timehrs,1)  
          do ks=1,nsed
            write(apbk,62) ks
            aname = 'Fraction' // trim(apbk) // alay
            call writescalh5(hotfile,hotpath,aname,pbk(:,ks,j),'none',timehrs,1)          
          enddo !j
        enddo !ks      
      endif
    endif
    
    !--- Salinity Transport ----
    if(saltrans)then
      call writescalh5(hotfile,hotpath,'Salinity',sal,'ppt',timehrs,1)  
    endif
    
    !--- Heat Transfer ----
    if(heattrans)then
      call writescalh5(hotfile,hotpath,'Heat',sal,'ppt',timehrs,1)  
    endif    
#endif

    return
    endsubroutine hot_write_xmdf

!!******************************************************************
!    subroutine hot_write_bins
!!******************************************************************
!    
!    return
!    endsubroutine hot_write_bin
    
!***********************************************************************    
    subroutine hot_print()
! Prints the Hot Start settings to the screen and diagnostic file
! written by Alex Sanchez, USACE-CHL
!***********************************************************************    
    use hot_def
    use diag_def, only: dgunit,dgfile
    implicit none
    integer :: i,iunit(2)
    
645 format(' ',A,F10.2,A)
745 format(' ',A,F6.2,A)
888 format(' ',A)
222 format(' ',A,A)
    
    iunit = (/6, dgunit/)
    
    open(dgunit,file=dgfile,access='append') 
    do i=1,2	
	  write(iunit(i),*)
	  if(coldstart)then
	    !write(iunit(i),'(A)') 'Hot Start:                      OFF'
        write(iunit(i),888) 'Start Mode:                     COLD'
	  else
	    !write(iunit(i),'(A)') 'Hot Start:                      ON'
        write(iunit(i),888) 'Start Mode:                     HOT'
	  endif
  	  if(hot_out)then
        write(iunit(i),888)   'Hot Start Output: '
	    write(iunit(i),222) '  File:                         ',trim(hotfile)
	    if(hot_timehr)then 
	      write(iunit(i),645)   '  Time:                      ',hottime,' hrs'
	    endif
	    if(hot_recur)then
	      write(iunit(i),745)   '  Recurring Interval:          ',hotdt,' hrs'
	    endif
	  endif
	enddo
	close(dgunit)

    return
    endsubroutine hot_print
    