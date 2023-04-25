!****************************************************************************************
    program CMS2D
! Coastal Modeling System (CMS)
!
! CMS License in 'Coastal Modeling System Terms and Conditions.txt'
! Instructions on how to compile in readme.txt
! Code changes recorded in logsheet.txt
!
! written by     
!   Weiming Wu   NCCHE     - hydrodynamics + sediment transport
!   Alex Sanchez USACE-CHL - sediment transport + hydrodynamics
!   Mitch Brown  USACE-CHL - auxiliary subroutines
!   Lihwa Lin    USACE-CHL - wave transformation
! notes
!   noptset==1 - CMS Wave only
!   noptset==2 - CMS Flow only
!   noptset==3 - CMS Flow/Wave Steering
!****************************************************************************************
#include "CMS_cpp.h"
#ifdef UNIT_TEST
    use CMS_test
#endif
    use cms_def
    use diag_def
    use comvarbl,   only: ctime,Version,Revision,release,developmental,rdate,nfsch,machine,major_version, minor_version, bugfix
    use hot_def,    only: coldstart
    use geo_def,    only: idmap,zb,x
    use sed_def,    only: db,d50,nlay,d90,pbk,nsed
    use size_def,   only: ncellsD
    use dredge_def, only: dredging


    implicit none
    integer k,j,ID,i,jlay,iper,loc,lstr
    character(len=20) :: astr,first,second
    real depthT
    real(ikind), allocatable :: dper(:)

    !Code version - moved here for easier modification when new descriptions are added
    !NOTE: Change variables Below to update header information
    version  = 5.3           ! CMS version         !For interim version
    revision = 4             ! Revision number
    bugfix   = 0             ! Bugfix number
    rdate    = '04/10/2023'

    !Manipulate to get major and minor versions - MEB  09/15/20
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
      call sim_start_print  !start timer here

      !Might check to see if dtsteer already set
      dtsteer=3.0
      ctime=0.0
      coldstart=.true.

      open(dgunit,file=dgfile,access='append') !Note: Needs to be open for CMS-Wave

      ! Lihwa fixed issues requiring two separate wave models (inline and stand-alone).  Only calling 'inline' option now.
      call CMS_Wave_inline  !added to be able to run inline wave model for checking results.
!      if (inlinewave) then
!        call CMS_Wave_inline  !added to be able to run inline wave model for checking results.
!      else
!        call CMS_Wave !(noptset,nsteer,dtsteer,ctime,coldstart)
!      endif
      
      close(dgunit)

      call sim_end_print
      
    else !if(noptset==2 .or. noptset==3) then
      !if(n2Dor3D==2) then
    !STACK:  OPEN & READ FLOW CARD FILE TO SEE IF IT EXPLICIT OR IMPLICIT
      !nfsch = 0, implicit     !STACK:
      !nfsch = 1, explicit     !STACK:        
      !nfsch = 0  !STACK:  implicit is defualt, nfsch updated from card file in SOLUTION_SCHEME_OPTION
      call SOLUTION_SCHEME_OPTION()       !STACK:
      if(nfsch.lt.0.or.nfsch.gt.1) then
        write(*,*)'error in determnining flow mode (EXP or IMP)'  !STACK:
        write(*,*)'solution scheme option nfsch = ',nfsch  !STACK:
      endif

      IF(nfsch.eq.1) then   !CALL EXPLICIT
        call CMS_FLOW_EXP_GRIDTYPE   
      ELSE                  !CALL IMPLICIT
        call CMS_Flow
      ENDIF
         !call CMS_Flow_Shock
      !else
         !call CMS3D_Flow
      !endif
    endif
#endif

    end program CMS2D
    
    
!*************************************
    subroutine get_com_arg
! Gets the command line arguments
! or runs an interactive input
!*************************************  
    use geo_def, only: grdfile,telfile
    use cms_def
    use comvarbl
    use diag_def
    use hot_def, only: coldstart
    use out_def, only: write_sup,write_tecplot,write_ascii_input
    use steer_def, only: auto_steer
    implicit none
    integer :: i,k,narg,nlenwav,nlenflow,ncase,ierr
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
            write(*,*) 'ERROR: ',trim(ctlfile),' does not exist'
            write(*,*) 'Press any key to continue.'
            read(*,*)
            stop
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
            write(*,*) 'ERROR: ',trim(astr),' does not exist'
            write(*,*) 'Press any key to continue.'
            read(*,*)
            stop
          endif    
          cmswave = .true.          
          noptset = 3 
          
        case('inline')
          inlinewave = .true.
          narg = narg -1
          
        case('tools')
          call CMS_tools_dialog

        case('flp')  !Old format input files       
          casename = aname !(1:ind-1)   
          flowpath = apath
          inquire(file=astr,exist=ok)
          if(.not.ok) then
            write(*,*) 'ERROR: ',trim(astr),' does not exist'
            write(*,*) 'Press any key to continue.'
            read(*,*)
            stop
          endif    
          cmsflow = .true.
          coldstart = .true.  
          write_sup = .true.
          write_tecplot = .true.
          write_ascii_input = .false.  ! Not overwrite input
        case default
          write(*,*) 'File not found: ',trim(astr)
          write(*,*) 'Press any key to continue.'
          read(*,*)
          stop
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
          write(*,*) 'ERROR: ',trim(astr),' does not exist'
          write(*,*) 'Press any key to continue.'
          read(*,*)
          stop
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
          write(*,*) 'ERROR: ',trim(astr),' does not exist'
          write(*,*) 'Press any key to continue.'
          read(*,*)
            stop
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
    
    return
    end subroutine get_com_arg
        
!************************************************************************    
    subroutine print_header
!************************************************************************        
    use comvarbl, only: version,revision,release,developmental,rdate,machine,major_version,minor_version,bugfix
    use diag_def
    
    implicit none
    integer      :: iunit(2),i
    character*22 :: string

7009  format(' **********************************************************')
7011  format('              U.S. Army Corps of Engineers                 ')
7012  format('            Coastal Inlets Research Program                ')
7013  format('                Coastal Modeling System                    ')
7014  format('       CMS2D, Version ',I0,'.',I0,'.',I0,1X,A,1X,A)
8014  format('       CMS2D, Version ',I0,'.',I0,'.',I0,'.',I0,1X,A,1X,A)
7114  format('          This version is for testing purposes only!       ')
7015  format(' Coupled Hydrodynamic, Wave, and Sediment Transport Model  ')
7016  format('               Last updated - ',A10)
7017  format('       For the latest version of CMS please visit          ')
7018  format('        https://cirpwiki.info/wiki/CMS_Releases            ')
      
9001  format('      By using this software the user has agreed to the    ')
9002  format('      terms and conditions of CMS license agreement.       ') 
9003  format('      A copy of the license can be obtained from the       ')
9004  format('      website shown above.                                 ')
     
8001  format('      This software is distributed on an "AS IS" basis     ')
8002  format('      WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND,        ')
8003  format('      either express or implied.                           ')
          
    !Declare file names    
    dgfile = 'CMS_DIAG.txt' !Diagnostic file is always in flow path
    dgunit = 9
    
    open(dgunit,file=dgfile,STATUS='unknown')   
    iunit = (/6,dgunit/)
    do i=1,2
      write(iunit(i),*)
      write(iunit(i),7009)
      write(iunit(i),7011)
      write(iunit(i),7012)
      write(iunit(i),7013)
      if(developmental) then    !DEVELOPMENTAL - this overrides the 'release' setting.
        string='DEVELOPMENTAL for'
        write(iunit(i),7114)
      elseif(.not.release)then  !BETA
        string='BETA for'
        write(iunit(i),7114)
      else                      !RELEASE
        string='RELEASE for'
      endif

      !Adding logic to show information for a bug fix and print it.      MEB  09/15/20
      if (bugfix == 0) then
        write(iunit(i),7014) major_version,minor_version,revision,trim(string),trim(machine)
      else
        write(iunit(i),8014) major_version,minor_version,revision,bugfix,trim(string),trim(machine)
      endif
 
      write(iunit(i),7016) rdate !Last revision date
      write(iunit(i),7017)
      write(iunit(i),7018)
      write(iunit(i),*)
      write(iunit(i),9001)
      write(iunit(i),9002)
      write(iunit(i),9003)
      write(iunit(i),9004)
      !write(iunit(i),8001)
      !write(iunit(i),8002)
      !write(iunit(i),8003)
      write(iunit(i),7009)
      write(iunit(i),*)
    enddo
    close(dgunit)        
    
    return
    end subroutine print_header

!*************************************************************
    subroutine sim_start_print
! Calculates and prints the elapsed cpu time
! by Mitch Brown and modified by Alex Sanchez, USACE-ERDC-CHL
!*************************************************************
#include "CMS_cpp.h"
    use cms_def
    use diag_def, only: dgunit,dgfile
    use comvarbl, only: ctime,ctime1
    use time_lib 
    
#ifdef PROFILE
    use watch_lib
#endif   
    use prec_def
    implicit none
    integer :: ita(8),i,iunit(2)
    character(len=200) :: datetimestr

645 format(' *************************')
684 format('    START OF SIMULATION   ')
            
740 format(I4,'-',I2.2,'-',I2.2,1x,I2.2,':',I2.2,':',I2.2)
750 format(2x,A)
    
    timebegin = time_cpu() !sec
    call date_and_time(values=ita)
    timestart = time_jul(ita) !sec
    timenow = timestart
    ctime1 = ctime
    
    call time_cal2str(datetimestr,ita(1),ita(2),ita(3),ita(5),ita(6),ita(7))
    
    open(dgunit,file=dgfile,access='append') 
    iunit = (/6,dgunit/)
    do i=1,2
      write(iunit(i),*)
      write(iunit(i),645)
      write(iunit(i),684)
      write(iunit(i),750) trim(datetimestr)
      write(iunit(i),645)  
    enddo
    close(dgunit)
    
#ifdef PROFILE
    call watch_start('Simulation')
#endif

    return
    end subroutine sim_start_print
        
!************************************************************************
    subroutine sim_end_print
! Calculates and prints the elapsed cpu time
!
! by Mitch Brown and modified by Alex Sanchez, USACE-ERDC-CHL
!************************************************************************
#include "CMS_cpp.h"
    use cms_def
    use comvarbl, only: stimet
    use diag_def, only: dgunit,dgfile
    use time_lib
#ifdef PROFILE
    use watch_lib
#endif    
    use prec_def
    implicit none
    integer :: ita(8),i,iunit(2)
    real*8 :: time_dur,speed
    character(len=200) :: clocktimestr,cputimestr,datetimestr

970 format(' *********************************')
971 format('         END OF SIMULATION        ')
740 format(2x,I4,'-',I2.2,'-',I2.2,1x,I2.2,':',I2.2,':',I2.2)    
840 format('  - Clock time: ',I0,' min, ',F7.4,' s')
841 format('  - Clock time: ',I0,' hrs, ',I0,' min, ',F7.4,' s')     
940 format('  - Computational Speed: ',F0.2)
720 format('  - CPU time:   ',I0,' min, ',F7.4,' s')
721 format('  - CPU time:   ',I0,' hrs, ',I0,' min, ',F7.4,' s')   
    
850 format('  - Clock time: ',A)    
730 format('  - CPU time:   ',A)
750 format(2x,A)

#ifdef PROFILE
    call watch_stop('Simulation')
#endif
    
    !Wall clock time
    call date_and_time(values=ita)
    
    time_dur = time_jul(ita) - timestart !sec
    call time_sec2str(time_dur,clocktimestr)

    call time_cal2str(datetimestr,ita(1),ita(2),ita(3),ita(5),ita(6),ita(7))
    
    !Speed 
    speed = dble(stimet)/max(time_dur,0.00001) 
    
    !CPU time
    time_dur = time_cpu() - timebegin
    call time_sec2str(time_dur,cputimestr)
    
    open(dgunit,file=dgfile,access='append')    
    iunit = (/6,dgunit/)     
    do i=1,2
      write(iunit(i),*)    
      write(iunit(i),970)
      write(iunit(i),971)
      write(iunit(i),750) trim(datetimestr)
      write(iunit(i),850) trim(clocktimestr)
      if(cmsflow)then
        write(iunit(i),940) speed        
      endif
      write(iunit(i),730) trim(cputimestr) 
      write(iunit(i),970)
    enddo
    close(dgunit)
    
#ifdef PROFILE
    call watch_output
    call watch_destroy
#endif   
    
    return
    end subroutine sim_end_print

!********************************************************************************
    subroutine cms_print
! Prints the general CMS settings to the screen and diagnostic file
! written by Alex Sanchez, USACE-CHL
!********************************************************************************    
    use comvarbl,  only: flowpath,ctlfile, input_ver, advfile, read_adv, SMS_ver
    use cms_def,   only: noptset,noptwse,noptvel,noptzb,wavsimfile,wavepath,  &
       dtsteer,radpath,wavpath,perpath,dirpath,disspath,                      &
       noptxtrpfl,xtrpdistfl,noptxtrpwav,xtrpdistwav
    use diag_def
    use geo_def,   only: wgrdfile
    use tool_def,  only: vstrlz    
    
    implicit none
    integer :: iunit(2),i, first, second
    character :: aname*200,apath*200,aext*10,astring*200,dstring*200,adate*8,atime*10,azone*5

341 format(' ',A,T40,A,A)     !Added for vstrlz function results
342 format(' ',A,T40,I0,'.',I0,A)  
887 format(' ',A,T40,A)
764 format(' ',A,T40,F0.3,A)  
    
    call DATE_AND_TIME (date=adate,time=atime,zone=azone)
    dstring=adate(5:6)//'/'//adate(7:8)//'/'//adate(1:4)
    dstring=trim(dstring)//' '//atime(1:2)//':'//atime(3:4)//' '//azone

    open(dgunit,file=dgfile,access='append') 
    iunit = (/6,dgunit/)
    
    do i=1,2
      write(iunit(i),*) 
      write(iunit(i),887) 'Actual Start Date/Time: ',trim(dstring)
      if(noptset>=2)then   
        if (flowpath /= '') write(iunit(i),887)    'CMS-Flow Path:',trim(flowpath)
        call fileparts(ctlfile,apath,aname,aext)
        astring=trim(aname) // '.' // aext
        write(iunit(i),887)  'CMS-Flow Card File:',trim(astring)
        
        call split_real_to_integers(input_ver,2,first,second)
        write(iunit(i),342)  'CMS Input Version:',first,second    
        
        call split_real_to_integers(SMS_ver,2,first,second)
        if (SMS_ver == -1) then
          write(iunit(i),887)  'SMS Version used:','Unknown (13.0 or previous)'
        else  
          write(iunit(i),342)  'SMS Version used:',first,second
        endif
      endif
      if (read_adv) then
        call fileparts(advfile,apath,aname,aext)
        astring=trim(aname) // '.' // aext
        write(iunit(i),887)    'Advanced Card File used:',trim(astring)
      endif
      if(noptset==1 .or. noptset==3)then
        write(iunit(i),887)    'CMS-Wave Path:',trim(wavepath)
        call fileparts(WavSimFile,apath,aname,aext)
        astring=trim(aname) // '.' // aext
        write(iunit(i),887)    'CMS-Wave Sim File:',trim(astring)   
      endif  
      if(noptset==3)then  
        write(iunit(i),887)    'Steering Mode:','ON'
        write(iunit(i),341)    'Steering Interval:',trim(vstrlz(dtsteer/3600.0,'(F0.3)')),' hrs'
        write(iunit(i),887)    'Wave-to-Flow Coupling:'
        write(iunit(i),887)    '  Temporal Interpolation:','LINEAR'
        !write(iunit(i),764)      '  Extrapolation Distance:   ',xtrpdistwav,' m'         
        write(iunit(i),*)
        write(iunit(i),887)    'Wave-to-Flow Coupling:'
        write(iunit(i),887)    '  Temporal Extrapolation: '
        write(iunit(i),887)    '    Water Level: '
        select case(noptwse)
        case(0)
          write(iunit(i),887)  '      wse(wave_time,wave_grid) = 0.0'     
        case(1)
          write(iunit(i),887)  '      wse(wave_time,wave_grid) = wse(flow_time,flow_grid)'    
        case(2)
          write(iunit(i),887)  '      wse(wave_time,wave_grid) = tide(wave_time,flow_grid)'
        case(3)
          write(iunit(i),887)  '      wse(wave_time,wave_grid) = wse(flow_time,flow_grid) '
          write(iunit(i),887)  '             + tide(wave_time) - tide(flow_time)'
        end select 
        write(iunit(i),887)    '    Current Velocities:'
        select case(noptvel)
        case(0)
          write(iunit(i),887)  '      vel(wave_time,wave_grid) = 0.0'     
        case(1)
          write(iunit(i),887)  '      vel(wave_time,wave_grid) = vel(flow_time,flow_grid)'
        end select
        write(iunit(i),887)    '    Bed Elevation: '
        select case(noptzb)
        case(0)
          write(iunit(i),887)  '      zb(wave_grid) = zb(wave_grid)'    
        case(1)
          write(iunit(i),887)  '      zb(wave_time,wave_grid) = zb(flow_time,flow_grid)'     
        case(2)
          write(iunit(i),887)  '      zb(wave_time,wave_grid) = zb(start_time,wave_grid) '    
          write(iunit(i),887)  '           + zb(flow_time,flow_grid) - zb(start_time,flow_grid)'    
        end select
        !write(iunit(i),764)    '  Extrapolation Distance:   ',xtrpdistfl,' m'        
      elseif(noptset==4)then    
        write(iunit(i),*) ' '
        write(iunit(i),887)    'Forcing with Single Wave Condition:',trim(wgrdfile)
        write(iunit(i),887)    '  Wave Height Dataset:'     ,'- '//trim(wavpath)
        write(iunit(i),887)    '  Wave Period Dataset:'     ,'- '//trim(perpath)      
        write(iunit(i),887)    '  Wave Direction Dataset:'  ,'- '//trim(dirpath)      
        write(iunit(i),887)    '  Wave Dissipation Dataset:','- '//trim(disspath)      
        write(iunit(i),887)    '  Radiation Stress Dataset:','- '//trim(radpath)
      endif      
   !   call fileparts(mpfile,apath,aname,aext)
      !astring=trim(aname) // '.' // aext
      !write(iunit(i),787)             'Parameters File:             ',trim(astring)
    enddo
    close(dgunit)
    
    return
    end subroutine cms_print
    
!********************************************************************************
    subroutine wave_only_print (simfile,iprpp,icur,ibreak,irs,kout,ibnd,  &
          iwet,ibf,iark,iarkr,akap,bf,ark,arkr,iwvbk,nonln,igrav,irunup,  &
          imud,iwnd,isolv,ixmdf,iproc,iview,iroll)
! Prints the known CMS-Wave parameters to the screen and diagnostic file
! written by Mitchell Brown, USACE-CHL
!********************************************************************************    
    use diag_def
    use tool_def, only: vstrlz
    
    implicit none
    integer, intent(in)      :: iprpp,icur,ibreak,irs,ibnd,iwet,ibf,iark,iarkr,iwvbk
    integer, intent(in)      :: igrav,irunup,imud,iwnd,isolv,ixmdf,iproc,iview,iroll,kout
    integer, intent(in)      :: nonln
    real, intent(in)         :: akap,bf,ark,arkr
    character(len=*), intent(in) :: simfile
    
    integer :: iunit(2),i
    logical :: isOpen
    character :: aname*200,apath*200,aext*10,astring*200

342 format(' ',A,T40,F0.2,A)   !also 890
354 format(' ',A,T40,A,A)     !Added for vstrlz function results
887 format(' ',A,T40,A)        !also 888
889 format(' ',A,T40,I0,A)     !also 891
    
    inquire(unit=dgunit,opened=isOpen)
    if(isOpen) close(dgunit)
    open(dgunit,file=dgfile,access='append') 
    iunit = (/6,dgunit/)
    call fileparts(simfile,apath,aname,aext)
    astring=trim(aname) // '.' // aext
    do i=1,2
      write(iunit(i),887)     'CMS-Wave Sim File:',trim(astring)
      write(iunit(i),887)     '  Mode: Waves Only'
      select case (iview)
        case(0)
          write(iunit(i),887) "  Half-Plane Spectral Forcing"
        case(1)
          write(iunit(i),887) "  Full-Plane Spectral Forcing (ignore 'wave.spc')"
        case(2)
          write(iunit(i),887) "  Full-Plane Spectral Forcing (read 'wave.spc')"
      end select
      select case (iprpp)
        case(0)
          write(iunit(i),887) 'Wave Propagation:','Waves and Wind'
        case(1)
          write(iunit(i),887) 'Wave Propagation:','Waves (neglect wind input)'
        case(-1)
          write(iunit(i),887) 'Wave Propagation:','Waves and Wind (Fast Mode)'
      end select

      select case (icur)
        case(0)
          write(iunit(i),887) 'Wave modified by current:','OFF'
        case(1)
          write(iunit(i),887) 'Wave modified by current:','ON, multiple currents read sequentially'
        case(2)
          write(iunit(i),887) 'Wave modified by current:','ON, first current set only'
      end select

      select case (ibreak)
        case(0)
          write(iunit(i),887) 'Breaking/Dissipation output:','OFF'
        case(1)
          write(iunit(i),887) 'Breaking/Dissipation output:','ON, Breaking Indices written'
        case(2)
          write(iunit(i),887) 'Breaking/Dissipation output:','ON, Dissipation Fluxes written'
      end select

      select case (irs)
        case(0)
          write(iunit(i),887) 'Radiation Stress Output:','OFF'
        case(1)
          write(iunit(i),887) 'Radiation Stress Output:','ON'
        case(2)
          write(iunit(i),887) 'Radiation Stress and Wave Setup:','Maximum water level written'
      end select

      select case (kout)
        case(0)
          write(iunit(i),887) 'Spectral and Parameter Output:','OFF'
        case default
          write(iunit(i),887) 'Spectral and Parameter Output:','ON'
          write(iunit(i),889) '  Output for ',kout,' cells'
      end select

      select case (ibnd)
        case(0)
          write(iunit(i),887) 'Nested Grid:','OFF'
        case(1)
          write(iunit(i),887) 'Nested Grid:','ON (Linear Interpolation)'
        case(2)
          write(iunit(i),887) 'Nested Grid:','ON (Morphic Interpolation)'
      end select

      select case (iwet)
        case(0)
          write(iunit(i),887) 'Wetting/Drying:','ON'
        case(1)
          write(iunit(i),887) 'Wetting/Drying:','OFF'
      end select

      select case (ibf)
        case(0)
          write(iunit(i),887) 'Bottom Friction:','OFF'
        case(1)
          write(iunit(i),887) 'Bottom Friction:','ON, Constant Darcy-Weisbach'
          write(iunit(i),354) '  Value: ',trim(vstrlz(bf,'(f0.2)'))
        case(2)
          write(iunit(i),887) 'Bottom Friction:','ON, Variable Darcy-Weisbach'
        case(3)
          write(iunit(i),887) "Bottom Friction:","ON, Constant Manning's"
          write(iunit(i),354) '  Value: ',trim(vstrlz(bf,'(f0.2)'))
        case(4)
          write(iunit(i),887) "Bottom Friction:","ON, Variable Manning's"
      end select

      select case (iark)
        case(0)
          write(iunit(i),887) 'Forward Reflection:','OFF'
        case(1)
          write(iunit(i),887) 'Forward Reflection:','ON, Constant'
          write(iunit(i),354) '  Value: ',trim(vstrlz(ark,'(f0.2)'))
        case(2)
          write(iunit(i),887) 'Forward Reflection:','ON, Variable'
      end select

      select case (iarkr)
        case(0)
          write(iunit(i),887) 'Backward Reflection:','OFF'
        case(1)
          write(iunit(i),887) 'Backward Reflection:','ON, Constant'
          write(iunit(i),354) '  Value: ',trim(vstrlz(arkr,'(f0.2)'))
        case(2)
          write(iunit(i),887) 'Backward Reflection:','ON, Variable'
      end select
      
      select case (nonln)
        case(0)
          write(iunit(i),887) 'Nonlinear Wave-Wave Interaction:','OFF'
        case(1)
          write(iunit(i),887) 'Nonlinear Wave-Wave Interaction:','ON'
      end select
      
      select case (igrav)
        case(0)
          write(iunit(i),887) 'Infragravity Waves:','OFF'
        case(1)
          write(iunit(i),887) 'Infragravity Waves:','ON'
      end select
      
      select case (irunup)
        case(0)
          write(iunit(i),887) 'Wave Runup Calculations:','OFF'
        case(1)
          write(iunit(i),887) 'Wave Runup Calculations:','ON (relative to absolute datum)'
        case(2)
          write(iunit(i),887) 'Wave Runup Calculations:','ON (relative to updated MWL)'
      end select
      
      select case (imud)
        case(0)
          write(iunit(i),887) "Muddy Bottom Calculations:","ON (read from 'mud.dat')"
        case(1)
          write(iunit(i),887) 'Muddy Bottom Calculations:','OFF'
      end select
      
      select case (iwnd)
        case(0)
          write(iunit(i),887) "Wind Forcing:","ON (if 'wind.dat' exists)"
        case(1)
          write(iunit(i),887) "Wind Forcing:","OFF"
        case(2)
          write(iunit(i),887) "Wind Forcing:","ON (dismiss incident wave inflation under stronger wind forcing)"
      end select
      
      select case (isolv)
        case(0)
          write(iunit(i),887) "Matrix Solver:","GSR"
        case(1)
          write(iunit(i),887) "Matrix Solver:","ADI (no parallelization)"
      end select
      
      select case (ixmdf)
        case(0)
          write(iunit(i),887) "Output:","ASCII"
        case(1)
          write(iunit(i),887) "Output:","XMDF"
        case(2)
          write(iunit(i),887) "Input/Output:","XMDF"
        case(-1)
          write(iunit(i),887) "Output:","ASCII"
      end select
      
      select case (iproc)
        case(0)
          write(iunit(i),887)   "Number of Processors:","1"
        case(1)
          write(iunit(i),887)   "Number of Processors:","1"
        case default
          if (isolv.eq.0) then
            write(iunit(i),889) "Number of Processors:",iproc
            write(iunit(i),887) "  Note: Processors should approximately equal Total Rows/300"
          else
            write(iunit(i),887) "Number of Processors:","1 (no parallelization with ADI)"
          endif
      end select
      
      select case (iwvbk)
        case(0)
          write(iunit(i),887) 'Wave Breaking Formula:','Extended Goda'
        case(1)
          write(iunit(i),887) 'Wave Breaking Formula:','Extended Miche'
        case(2)
          write(iunit(i),887) 'Wave Breaking Formula:','Battjes and Janssen'
        case(3)
          write(iunit(i),887) 'Wave Breaking Formula:','Chawla and Kirby'
        case(4)
          write(iunit(i),887) 'Wave Breaking Formula:','Battjes and Janssen (2007)'
        case(5)
          write(iunit(i),887) 'Wave Breaking Formula:','Miche (original)'
        case(6)
          write(iunit(i),887) 'Wave Breaking Formula:','Lifting Breaking'
      end select

      write(iunit(i),889) 'Wave Roller effect:',iroll,' (min. 0, max. 4)'
      write(iunit(i),354) 'Diffraction Intensity Factor:',trim(vstrlz(akap,'(f0.2)')),' (min. 0, max. 4)'
      write(iunit(i),887) ''
    enddo

    write(dgunit,*) '*** Starting CMS-Wave Run ***'
    close(dgunit)
    
    return
    end subroutine wave_only_print
     
