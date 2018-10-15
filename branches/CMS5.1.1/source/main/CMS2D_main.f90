!****************************************************************************************
	program CMS2D
! Coastal Modeling System (CMS)
!
! CMS License in 'Coastal Modeling System Terms and Conditions.txt'
! Instructions on how to compile in readme.txt
! Code changes recorded in logsheet.txt
!
! written by     
!   Weimings Wu, NCCHE - hydrodynamics + sediment transport
!   Alex Sanchez, USACE-CHL - sediment transport + hydrodynamics
!   Lihwa Lin, USACE-CHL - wave transformation
!   Mitch Brown, USACE-CHL - auxiliary subroutines
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
    use comvarbl,   only: ctime,Version,Revision,release,rdate,nfsch
    use hot_def,    only: coldstart
    use geo_def,    only: idmap,zb,x
    use sed_def,    only: db,d50,nlay,d90,pbk,nsed
    use size_def,   only: ncellsD
    
#ifdef DREDGE    
    use dredge_def, only: dredging
#endif

    implicit none
    integer k,j,ID,i,jlay,iper
    real depthT
    real(ikind), allocatable :: dper(:)

    !Code version - moved here for easier modification when new descriptions are added
    !NOTE: Change variables Below to update header information
    Version  = 5.1   !CMS version
    Revision = 1     !Revision number
#ifdef DEV_MODE
    release  = .false.
#else
    release = .true.
#endif
    rdate    = '10/04/2018'
    
    !n2Dor3D=2    !=2 for 2D; =3 for 3D    
    
#ifdef UNIT_TEST
    call CMS_test_run
    stop
#else
    call steering_default 
    call get_com_arg      !get command line arguments
    call print_header     !screen and debug file header  
    
	if(noptset==1) then   !CMS-Wave only, use Stand alone code   MEB  10/15/2018
	  call sim_start_print  !start timer here
!Might check to see if dtsteer already set
      dtsteer=3.0
!
      ctime=0.0
      coldstart=.true.

      open(dgunit,file=dgfile,access='append') !Note: Needs to be open for CMS-Wave
      call CMS_Wave !(noptset,nsteer,dtsteer,ctime,coldstart)       
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

#ifdef DREDGE
    !special output for drege-multiple grains size CHETN'
    if (dredging) then
      allocate(dper(ncellsD))
      iper = 50
      open(unit=4044,file='beddist.txt')
      write(4044,*)'X Z d50 F1 F2 F3'
      do k=2517,2539
        id = idmap(k)        
        depthT= -6.8 + 20*(zb(ID)+6.8)
        j=1
        jlay = j
        call sedpercentile_bedlayer(iper,jlay,dper)        
        write(4044,"(6e15.7)")X(ID),depthT,1000*dper(ID),(pbk(id,i,j),i=1,nsed) 
        j=1
        jlay = j
        call sedpercentile_bedlayer(iper,jlay,dper)          
        depthT = depthT - 20*db(id,j)/2.0               
        write(4044,"(6e15.7)")X(ID),depthT,1000*dper(ID),(pbk(id,i,j),i=1,nsed)        
        do j=2,nlay
          jlay = j
          call sedpercentile_bedlayer(iper,jlay,dper)              
          depthT = depthT - 20*db(id,j-1)/2.0 - 20*db(id,j)/2.0
          write(4044,"(6e15.7)")X(ID),depthT,1000*dper(ID),(pbk(id,i,j),i=1,nsed)
        enddo
        j=nlay
        jlay = j
        call sedpercentile_bedlayer(iper,jlay,dper)          
        depthT = depthT - 20*db(id,j)/2.0               
        write(4044,"(6e15.7)")X(ID),depthT,1000*dper(ID),(pbk(id,i,j),i=1,nsed)          
      enddo
    endif !dredging
#endif
    
    endprogram CMS2D
    
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
    implicit none
    integer :: i,k,narg,nlenwav,nlenflow,ncase,ierr
    character :: cardname*37,aext*10, answer
    character(len=200) :: astr,apath,aname
    logical :: ok
      
!!    narg = iargc() !Not supported by PGI
    narg = command_argument_count()
!    if(narg==2) then  
!       write(*,*) 'WARNING: No steering interval specified'
!       write(*,*) '         Using default value of ',dtsteer/3600.0, ' hours'
!       write(*,*)        
!    endif      
    
    do i=0,min(narg,2)
      if(i==0 .and. narg==0)then
!CMS was called with no arguments
        write(*,*) ' '
        write(*,*) 'Type name of CMS-Flow Card File or '
        write(*,*) '  CMS-Wave Sim File and Press <RETURN>'
        write(*,*) ' '
        read(*,*) astr
      elseif(i==0 .and. narg>0)then
        cycle  
      else
        call getarg(i,astr)
      endif
      call fileparts(astr,apath,aname,aext)                    
      selectcase(aext)
        case('cmcards') !Flow model
          ctlfile = trim(aname) // '.cmcards'
          flowpath = apath
          casename = aname
          inquire(file=astr,exist=ok)
	      if(.not.ok)then
            write(*,*) 'ERROR: ',trim(astr),' does not exist'
            write(*,*) 'Press any key to continue.'
	        read(*,*)
	        stop
	      endif	
          cmsflow = .true.
          !Search for Steering Cards
          open(77,file=astr)
          do k=1,1000
	        read(77,*,iostat=ierr) cardname	
            if(ierr/=0) exit
	        call steering_cards(cardname)
          enddo
          close(77)

        case('sim') !Wave model
          WavSimFile = trim(aname) // '.sim'
          Wavepath=apath
          wavename=aname
          inquire(file=astr,exist=ok)
  	      if(.not.ok) then
            write(*,*) 'ERROR: ',trim(astr),' does not exist'
            write(*,*) 'Press any key to continue.'
	        read(*,*)
	        stop
	      endif	
          cmswave = .true.          
          noptset = 3 

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
      endselect      
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
        do k=1,1000
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
        dtsteer = dtsteer*3600.0  !Convert from hours to seconds 
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
    
    !Declare file names    
    dgfile = 'CMS_DIAG.txt' !Diagnostic file is always in flow path
    dgunit = 9
    ncase = len_trim(casename)  
    !dgfile  = casename(1:ncase) // '_diag.txt'
    if(cmsflow .and. nlenflow>0)then
      dgfile  = flowpath(1:nlenflow) // dgfile    
    elseif(cmswave .and. nlenwav>0 .and. .not.cmsflow)then !if only waves and there is a path, put diagnostic file there.
      dgfile  = wavepath(1:nlenwav) // dgfile
    endif
    
    if(noptset>=2)then !Only needed if running flow
      ctlfile  = flowpath(1:nlenflow) // ctlfile     
      grdfile  = flowpath(1:nlenflow) // casename(1:ncase) // '_grid.h5'
      mpfile   = flowpath(1:nlenflow) // casename(1:ncase) // '_mp.h5'    
      telfile  = flowpath(1:nlenflow) // casename(1:ncase) // '.tel'      
    endif
    
    return
    endsubroutine get_com_arg
		
!************************************************************************    
    subroutine print_header
!************************************************************************        
    use comvarbl, only: version,revision,release,rdate
    use diag_def
    
    implicit none
    integer :: iunit(2),i

7009  format(' **********************************************************')
7011  format('              U.S. Army Corps of Engineers                 ')
7012  format('            Coastal Inlets Research Program                ')
7013  format('                Coastal Modeling System                    ')
7014  format('           CMS2D, Version ',F5.2,'.',I2.2,1X,A)
7114  format('          This version is for testing purposes only!       ')
7015  format(' Coupled Hydrodynamic, Wave, and Sediment Transport Model  ')
7016  format('               Last updated - ',A10)
7017  format('       For the latest version of CMS please visit          ')
7018  format('          http://cirp.usace.army.mil/products/             ')

9001  format('      By using this software the user has agreed to the    ')
9002  format('      terms and conditions of CMS license agreement.       ') 
9003  format('      A copy of the license can be obtained from the       ')
9004  format('      website shown above.                                 ')
     
8001  format('      This software is distributed on an "AS IS" basis     ')
8002  format('      WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND,        ')
8003  format('      either express or implied.                           ')
          
    open(dgunit,file=dgfile,STATUS='unknown')   
    iunit = (/6,dgunit/)
    do i=1,2
      write(iunit(i),*)
      write(iunit(i),7009)
      write(iunit(i),7011)
      write(iunit(i),7012)
      write(iunit(i),7013)
      if(.not.release)then                            !BETA
        write(iunit(i),7014) version,revision,'BETA'
        write(iunit(i),7114)
      else                                               !RELEASE
        write(iunit(i),7014) version,revision,'RELEASE'
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
    endsubroutine print_header

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

645 format('*************************')
684 format('  START OF SIMULATION ')
740 format(2x,I4,'-',I2.2,'-',I2.2,1x,I2.2,':',I2.2,':',I2.2)
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
    endsubroutine sim_start_print
		
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

970 format('*********************************')
971 format('  END OF SIMULATION ')
740 format(2x,I4,'-',I2.2,'-',I2.2,1x,I2.2,':',I2.2,':',I2.2)    
840 format('  - Clock time: ',I0,' min, ',F7.4,' s')
841 format('  - Clock time: ',I0,' hrs, ',I0,' min, ',F7.4,' s')     
940 format('  - Computational Speed: ',F8.2)
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
    speed = dble(stimet)/time_dur
    
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
    endsubroutine sim_end_print

 !********************************************************************************
    subroutine cms_print
! Prints the general CMS settings to the screen and diagnostic file
! written by Alex Sanchez, USACE-CHL
!********************************************************************************    
    use comvarbl, only: flowpath,ctlfile, input_ver
    use cms_def, only: noptset,noptwse,noptvel,noptzb,wavsimfile,wavepath,&
       dtsteer,radpath,wavpath,perpath,dirpath,disspath,&
       noptxtrpfl,xtrpdistfl,noptxtrpwav,xtrpdistwav
    use diag_def
    
    implicit none
    integer :: iunit(2),i
    character :: aname*200,apath*200,aext*10,astring*200

342 format(' ',A,F5.2,A)    
887 format(' ',A,1x,A)
764 format(' ',A,F8.3,A)
888 format(' ',A)    
    
    open(dgunit,file=dgfile,access='append') 
    iunit = (/6,dgunit/)
    do i=1,2
      if(noptset>=2)then   
        write(iunit(i),*)      
        write(iunit(i),887)  'CMS-Flow Path:                 ',trim(flowpath)
        call fileparts(ctlfile,apath,aname,aext)
	    astring=trim(aname) // '.' // aext
        write(iunit(i),887)  'CMS-Flow Card File:            ',trim(astring)
        write(iunit(i),342)  'Input Version:                 ',input_ver     
      endif
      if(noptset==1 .or. noptset==3)then
        write(iunit(i),*)
        write(iunit(i),887)  'CMS-Wave Path:                 ',trim(wavepath)
        call fileparts(WavSimFile,apath,aname,aext)
	    astring=trim(aname) // '.' // aext
        write(iunit(i),887)  'CMS-Wave Sim File:             ',trim(astring)   
      endif  
      if(noptset==3)then  
        write(iunit(i),*)
        write(iunit(i),888)    'Steering Mode:                  ON'
        write(iunit(i),764)      'Steering Interval:           ',dtsteer/3600.0,' hrs'
        write(iunit(i),888)    'Wave-to-Flow Coupling:'
        write(iunit(i),888)    '  Temporal Interpolation:       LINEAR'
        !write(iunit(i),764)      '  Extrapolation Distance:   ',xtrpdistwav,' m'         
        write(iunit(i),*)
        write(iunit(i),888)    'Wave-to-Flow Coupling:'
        write(iunit(i),888)    '  Temporal Extrapolation: '
        write(iunit(i),888)    '    Water Level: '
        selectcase(noptwse)
        case(0)
          write(iunit(i),888)  '      wse(wave_time,wave_grid) = 0.0'     
        case(1)
          write(iunit(i),888)  '      wse(wave_time,wave_grid) = wse(flow_time,flow_grid)'    
        case(2)
          write(iunit(i),888)  '      wse(wave_time,wave_grid) = tide(wave_time,flow_grid)'
        case(3)
          write(iunit(i),888)  '      wse(wave_time,wave_grid) = wse(flow_time,flow_grid) '
          write(iunit(i),888)  '             + tide(wave_time) - tide(flow_time)'
        endselect 
        write(iunit(i),888)    '    Current Velocities:'
        selectcase(noptvel)
        case(0)
          write(iunit(i),888)  '      vel(wave_time,wave_grid) = 0.0'     
        case(1)
          write(iunit(i),888)  '      vel(wave_time,wave_grid) = vel(flow_time,flow_grid)'
        endselect
        write(iunit(i),888)    '    Bed Elevation: '
        selectcase(noptzb)
        case(0)
          write(iunit(i),888)  '      zb(wave_grid) = zb(wave_grid)'    
        case(1)
          write(iunit(i),888)  '      zb(wave_time,wave_grid) = zb(flow_time,flow_grid)'     
        case(2)
          write(iunit(i),888)  '      zb(wave_time,wave_grid) = zb(start_time,wave_grid) '    
          write(iunit(i),888)  '           + zb(flow_time,flow_grid) - zb(start_time,flow_grid)'    
        endselect
        !write(iunit(i),764)  '  Extrapolation Distance:   ',xtrpdistfl,' m'        
      elseif(noptset==4)then    
        write(iunit(i),887)  'Wave Height Dataset:           ',trim(wavpath)
        write(iunit(i),887)  'Wave Period Dataset:           ',trim(perpath)      
        write(iunit(i),887)  'Wave Direction Dataset:        ',trim(dirpath)      
        write(iunit(i),887)  'Wave Dissipation Dataset:      ',trim(disspath)      
        write(iunit(i),887)  'Radiation Stress Dataset:      ',trim(radpath)
      endif      
   !   call fileparts(mpfile,apath,aname,aext)
	  !astring=trim(aname) // '.' // aext
	  !write(iunit(i),787)             'Parameters File:             ',trim(astring)
    enddo
    close(dgunit)
    
    return
    endsubroutine cms_print
    
     