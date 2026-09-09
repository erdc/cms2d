!****************************************************************************************
    program CMS2D
! Coastal Modeling System (CMS)
!
! See 'CMS Terms and Conditions.txt' and 'LICENSE.md'
! 
! Instructions on how to compile in readme.txt
! Code changes recorded in logsheet.txt
!
! written by     
!   Weiming Wu   NCCHE     - hydrodynamics + sediment transport
!   Alex Sanchez USACE-CHL - sediment transport + hydrodynamics
!   Lihwa Lin    USACE-CHL - wave transformation
!   Honghai Li   USACE-CHL - auxiliary subroutines
!   Mitch Brown  USACE-CHL - auxiliary subroutines
! notes
!   noptset == 1 - CMS Wave only            nfsch == 0 - Implicit solver
!   noptset == 2 - CMS Flow only            nfsch == 1 - Explicit solver
!   noptset == 3 - CMS Flow/Wave Steering
!****************************************************************************************
#include "CMS_cpp.h"
#ifdef UNIT_TEST
    use CMS_test
#endif
    use cms_def,    only: noptset, dtsteer
    use hot_def,    only: coldstart
    use geo_def,    only: idmap,zb,x
    use sed_def,    only: db,d50,nlay,d90,pbk,nsed
    use prec_def,   only: ikind
    use diag_def,   only: dgunit, dgfile, msg, msg2
    use diag_lib,   only: diag_print_error
    use comvarbl,   only: ctime,Version,Revision,release,developmental,rdate,nfsch,machine,major_version,minor_version,bugfix
    use size_def,   only: ncellsD
    use dredge_def, only: dredging
    implicit none
    
    integer k,j,ID,i,jlay,iper,loc,lstr
    character(len=20) :: astr,first,second
    real depthT
    real(ikind), allocatable :: dper(:)
    
    interface
      character(len=20) function Int2Str(k)
        integer, intent(in) :: k
      end function
    end interface

    !NOTE: Change variables below to update CMS header information
    version  = 5.4           ! CMS version         !For interim version
    revision = 8             ! Revision number
    bugfix   = 0             ! Bugfix number
    rdate    = '09/09/2026'

    !Manipulate to get major and minor versions - MEB  09/15/2020
    call split_real_to_integers (version, 2, major_version, minor_version)  !Convert version to two integer portions before and after the decimal considering 2 digits of precision.
  
#ifdef _WIN32
    machine='Windows'
#elif defined (__linux) 
    machine='Linux' 
#else
    machine='Unknown'
#endif
    
#ifdef DEV_MODE
    release = .false.
#else
    release = .true.
#endif

developmental = .false.      !Change this to .false. for truly RELEASE code   meb  05/11/2022
    
#ifdef UNIT_TEST
    call CMS_test_run
    stop
#else
    call steering_default 

    call print_header     !screen and debug file header 
    call get_com_arg      !get command line arguments

    if(noptset==1) then   !CMS-Wave only
!Figure out how to call wave_only_print
      call sim_start_print  !start timer here
      dtsteer=3.0
      ctime=0.0
      coldstart=.true.
      open(dgunit,file=dgfile,access='append') !Note: Needs to be open for CMS-Wave
      call CMS_Wave_inline  !added to be able to run inline wave model for checking results.
      close(dgunit)
      call sim_end_print
    else !if(noptset==2 .or. noptset==3) then
      call SOLUTION_SCHEME_OPTION()                               !STACK:
      if(nfsch < 0 .or. nfsch > 1) then
        msg  = 'Cannot determine flow mode (EXP or IMP)'
        msg2 = 'Solution scheme option nfsch = '//Int2Str(nfsch)
        call diag_print_error(msg, msg2)                          !STACK:
      endif

      IF(nfsch.eq.1) then   !CALL EXPLICIT
        call CMS_FLOW_EXP_GRIDTYPE   
      ELSE                  !CALL IMPLICIT
        call CMS_Flow
      ENDIF
    endif
#endif

    end program CMS2D
    
!*************************************
    subroutine get_com_arg
! Gets the command line arguments
! or runs an interactive input
!*************************************  
    use geo_def,   only: grdfile,telfile
    use cms_def,   only: cmsflow, cmswave, wavsimfile, wavepath, wavename, noptset, dtsteer, noptwse, noptvel, noptzb
    use comvarbl,  only: casename, flowpath, mpfile, ctlfile
    use diag_def,  only: dgfile, msg, interactive
    use diag_lib,  only: diag_print_error
    use hot_def,   only: coldstart
    use out_def,   only: write_sup,write_tecplot,write_ascii_input
    use steer_def, only: auto_steer
    
    implicit none
    integer :: i,k,narg,nlenwav,nlenflow,ncase,ierr,iloc
    character :: cardname*37,aext*10, answer,laext*10   !added variable to hold lowercase version of the extension for the case statement.
    character(len=200) :: astr,apath,aname
    logical :: ok
    
    interface
      function toLower (astr)
        character(len=*),intent(in) :: astr
        character(len=len(astr)) :: toLower
      end function

      function findCard(aFile,aCard,aValue)
        character(len=*),intent(in)    :: aFile
        character(len=*),intent(in)    :: aCard
        character(len=100),intent(out) :: aValue
        logical :: findCard
      end function
    end interface    
      
    narg = command_argument_count()
    
    !Parse for --non-interactive flag
    call getarg(narg,astr)
    iloc = index(astr,'--non-interactive')
    if (iloc > 0) then
      print*,'Found --non-interactive flag, skipping read from keyboard on Errors'
      interactive = .false.
      narg = narg - 1
    endif
    
    do i=0,min(narg,2)
      if(i==0 .and. narg==0)then      !CMS was called with no arguments
        write(*,*) ' '
        write(*,*) 'Enter CMS-Flow Card File, CMS-Wave Sim File, or "Tools" and Press <RETURN>'
        write(*,*) ' '
        read(*,*) astr
      elseif(i==0 .and. narg>0)then
        cycle  
      else
        call getarg(i,astr)
      endif

      call fileparts(astr,apath,aname,aext)           
      astr = toLower(trim(astr))              !moved below previous line to retain the exact filename - meb 05/15/2020
      laext = toLower(trim(aext))             !needed to compare the lower-case version of the extension but retain the original case - meb 05/21/2020
      if (astr == 'inline' .or. astr == 'tools') laext=astr
      select case(laext)
        case('cmcards') !Flow model
          ctlfile = trim(aname) // '.' // trim(aext)
          flowpath = apath
          casename = aname
          inquire(file=ctlfile,exist=ok)
          if(.not.ok)then
            write(msg,*) trim(ctlfile),' does not exist'
            call diag_print_error(msg)
          endif    
          cmsflow = .true.
          !Search for Steering Cards
          open(77,file=ctlfile)
          do
            read(77,*,iostat=ierr) cardname    
            if(ierr/=0) exit
            call steering_cards(cardname)
          enddo
          close(77)

        case('sim') !Wave model
          WavSimFile = trim(aname) // '.' // trim(aext)  !'.sim'
          Wavepath=apath
          wavename=aname
          inquire(file=WavSimFile,exist=ok)
          if(.not.ok) then
            write(msg,*) trim(astr),' does not exist'
            call diag_print_error(msg)
          endif    
          cmswave = .true.          
          noptset = 3 
          
        case('tools')
          call CMS_tools_dialog

        case('flp')  !Old format input files       
          casename = aname !(1:ind-1)   
          flowpath = apath
          inquire(file=astr,exist=ok)
          if(.not.ok) then
            write(msg,*) trim(astr),' does not exist'
            call diag_print_error(msg)
          endif    
          cmsflow = .true.
          coldstart = .true.  
          write_sup = .true.
          write_tecplot = .true.
          write_ascii_input = .false.  ! Not overwrite input
        case default
          write(msg,*) 'File not found: ',trim(astr)
          call diag_print_error(msg)
      end select      
    enddo

!CMS was called with no arguments but user entered Flow parameter filename, also ask for Wave info
    if(cmsflow .and. .not.cmswave .and. narg==0)then
      write(*,*) ' '
      write(*,*) 'Type name of CMS-Wave Sim File'
      write(*,*) 'or type 0 for none and Press <RETURN>'
      read(*,*) astr
      if(astr(1:1)/='0' .and. astr(1:4)/='none')then
        call fileparts(astr,apath,aname,aext)          
        inquire(file=astr,exist=ok)
        if(.not.ok) then
          write(msg,*) trim(astr),' does not exist'
          call diag_print_error(msg)
        endif
        WavSimFile = trim(aname) // '.sim'
        Wavepath=apath
        wavename=aname
        cmswave = .true.          
        noptset = 3           
      endif
!CMS was called with no arguments but user entered Wave parameter filename, also ask for Flow info
    elseif(.not.cmsflow .and. cmswave .and. narg==0)then
      write(*,*) ' '
      write(*,*) 'Type name of CMS-Flow Card File'
      write(*,*) '  or type 0 for none and Press <RETURN>'
      write(*,*) ' '
      read(*,*) astr         
      if(astr(1:1)/='0' .and. astr(1:4)/='none')then
        inquire(file=astr,exist=ok)
        if(.not.ok) then
          write(msg,*) trim(astr),' does not exist'
          call diag_print_error(msg)
        endif
        ctlfile = trim(aname) // '.cmcards'
        flowpath = apath
        casename = aname  
        cmsflow = .true.
        !Search for Steering Cards
        open(77,file=astr)
        do
          read(77,*,iostat=ierr) cardname
          if(ierr/=0) exit
          call steering_cards(cardname)              
        enddo
        close(77)         
      endif
    endif

    if(cmswave .and. .not.cmsflow .and. narg<=1)then
      noptset = 1  !CMS-Wave only
      casename = wavename
    elseif(.not.cmswave .and. cmsflow .and. narg<=1)then  
      noptset = 2  !CMS-Flow only
    elseif(cmswave .and. cmsflow .and. narg>=2)then  
      noptset = 3  !CMS-Flow and CMS-Wave
!    elseif(cmswave .and. cmsflow .and. narg==1)then  
!      noptset = 4  !CMS-Flow and wave input
    endif    
    
    !Get steering interval if specified
    if(narg>=3)then !Alex, bug fix, changed == to >=
      call getarg(3,astr)
      read(astr,*) dtsteer          
      dtsteer = dtsteer*3600.0  !Convert from hours to seconds 
    endif
    
    if(noptset==3 .and. dtsteer<0.0)then
      if(narg==0)then
        write(*,*) 'Type the steering interval'
        write(*,*) 'in hours and Press <RETURN>'
        write(*,*) ' '
        read(*,*) dtsteer
        if (dtsteer > 0) then
            dtsteer = dtsteer*3600.0  !Convert from hours to seconds 
        elseif (dtsteer == 0) then
            dtsteer = 10800.0
        else
            auto_steer = .true.
            dtsteer = 10800.0 !Initially just set to 3-hour default.
        endif
      else
        dtsteer = 10800.0 !3 hours ********************
      endif
    endif

    !Wave water level
    if(narg>=4)then
      call getarg(4,astr)
      read(astr,*) noptwse       
    endif
    
    !Wave current velocity
    if(narg>=5)then
      call getarg(5,astr)
      read(astr,*) noptvel       
    endif
    
    !Wave bed elevation
    if(narg>=6)then
      call getarg(6,astr)
      read(astr,*) noptzb       
    endif 
    
    if(noptset==3)then
      if(narg==0)then
        if(WavSimFile(1:1) == ' ') then
          write(*,*) ' '
          write(*,*) 'Select the method for estimating the wave '
          write(*,*) 'water levels from below and Press <RETURN>' 
          write(*,*) '  0 - wse(wave_time,wave_grid)=0.0'
          write(*,*) '  1 - wse(wave_time,wave_grid)=wse(flow_time,flow_grid)'
          write(*,*) '  2 - wse(wave_time,wave_grid)=tide(wave_time,flow_grid)'
          write(*,*) '  3 - wse(wave_time,wave_grid)=wse(flow_time,flow_grid) '
          write(*,*) '        +tide(wave_time)-tide(flow_time)'
          read(*,*) noptwse
          noptwse = max(min(noptwse,3),0)
          write(*,*) ' '
          write(*,*) 'Select the method for estimating the wave '
          write(*,*) 'current velocities from below and Press <RETURN>' 
          write(*,*) '  0 - vel(wave_time,wave_grid)=0.0'
          write(*,*) '  1 - vel(wave_time,wave_grid)=vel(flow_time,flow_grid)'
          read(*,*) noptvel
          noptvel = max(min(noptvel,1),0)
          write(*,*) ' '
          write(*,*) 'Select the method for estimating the wave '
          write(*,*) 'bed elevations from below and Press <RETURN>' 
          write(*,*) '  0 - zb(wave_grid)=zb(wave_grid)'
          write(*,*) '  1 - zb(wave_time,wave_grid)=zb(flow_time,flow_grid)'
          read(*,*) noptzb
          noptzb = max(min(noptzb,1),0)
        else
          write(*,*) ' '
          write(*,*) '--Using steering card values or wave steering parameter defaults'
        endif
      endif
    endif
    
    !If wave path is empty than use path for flow    
    nlenwav = len_trim(wavepath)
    nlenflow = len_trim(flowpath)
    if(nlenwav==0 .and. nlenflow>0)then
      wavepath = flowpath !*******************  
      nlenwav = len_trim(wavepath)
    endif
    
    ncase = len_trim(casename)  

    if(cmsflow .and. nlenflow>0)then
      dgfile  = flowpath(1:nlenflow) // dgfile    
    elseif(cmswave .and. nlenwav>0 .and. .not.cmsflow)then !if only waves and there is a path, put diagnostic file there.
      dgfile  = wavepath(1:nlenwav) // dgfile
    endif
    
    aext=''
    if(noptset>=2)then !Only needed if running flow
      ctlfile  = flowpath(1:nlenflow) // ctlfile     
      mpfile   = flowpath(1:nlenflow) // casename(1:ncase) // '_mp.h5'   

      if (findCard(ctlfile,'GRID_FILE',grdfile)) then      !Get the value for GRID_FILE from the ctlfile
        call removequotes(grdfile)                           !If result is in quotes, remove the quotes so it can properly call next subroutine
        call fileext(grdfile,aext)                           !Return the extension of the gridfile
      endif
      if (aext /= 'cart') telfile  = flowpath(1:nlenflow) // casename(1:ncase) // '.tel'         !If gridfile extension is 'cart', leave TELFILE = ' '
    endif
    
  end subroutine get_com_arg
  