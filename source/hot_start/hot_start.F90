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
!            Mitchell Brown, USACE-CHL
!==================================================================    

!******************************************************************    
    subroutine hot_default
! Sets hot start default valus
! written by Alex Sanchez, USACE-CHL    
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
    
    autohotname = 'AutoHotStart'
    hotname     = 'SingleHotStart'
    autohotfile = trim(AutoHotName)//autohot_inc//'.h5'    !Adding '_1' to the automatic hot start file name for XMDF, to prepare to alternate writing to the file.
    hotfile     = trim(HotName)//'.h5'
    autohotpath = 'Datasets/'
    hotpath     = 'Datasets/'    

    icfile = 'Initial_Cond_File.h5'
    icpath = 'Datasets/'    
    ictime = -999.0  !hrs
    npath  = len_trim(flowpath)
    icpres = .false.
    icwse  = .false.
    icvel  = .false.
    icwet  = .false.
    icflux = .false.
    
    add_duration_HS = .false.  !MEB  04/29/2022
    
    return
    end subroutine hot_default
    
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
    
    integer     :: iyrhot,imohot,idayhot,ihrhot,iminhot,isechot,ierr
    real(ikind) :: tjuldayhot
    character   :: apath*200,aname*200,aext*10,atemp*200
    logical     :: foundcard
    character(len=37) :: cardname    
    
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
      !----- Initial Condition (Input) ----------------------------------------
      !If another is added to the list, modify the appropriate line in 'input.F90' for subroutine 'card_dataset'
      case('INITIAL_STARTUP_FILE','INITIAL_CONDITION_FILE','HOT_START_SIMULATION') 
        call card_dataset(77,icfile,flowpath,icfile,icpath,1)
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
        
      case('EXTEND_DURATION_RUN_FOR_HOTSTART')
        call card_boolean(77,add_duration_HS,ierr)
        
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
        !backspace(77)
        !read(77,*) cardname, hotdt      
        call card_scalar(77,'hrs','hrs',hotdt,ierr)
        hot_recur = .true.      
        hot_out = .true.   
        
      case default
        foundcard = .false.           
        
    end select
    
    return
    end subroutine hot_cards

!******************************************************************    
    subroutine hot_read()
! Reads hot start information
! written by Alex Sanchez, USACE-CHL    
! modified by Mitch Brown - 05/18/2012
!******************************************************************
#include "CMS_cpp.h"  
    use hot_def
    use diag_lib
    use flow_def, only: eta,u,v
    use comvarbl, only: timehrs,stime,ctime
    implicit none

    character(len=10) :: aext
    logical :: ok
    
    !Check if hot_start file exists
    inquire(file=icfile,exist=ok)
    if(.not.ok)then
      call diag_print_error('Unable to find Initial Condition File: ',icfile)
    endif
    
    call fileext(icfile,aext)
    select case(aext)
    case('h5')
#ifdef XMDF_IO
      call hot_read_xmdf
#else
      call diag_print_error('Cannot read hot start from *.h5 file without XMDF libraries')
#endif
    
    case('sup')
      call hot_read_sup(icfile) !Implementation in progress
      
    case default
      call diag_print_error('Unknown file type with extension, '//trim(aext))
    end select

    return
    end subroutine hot_read
    
!********************************************************************
    subroutine hot_read_sup(supfile)
! Reads initial condition file from an SMS ASCII dataset (*.dat) file 
! written by Alex Sanchez, USACE-CHL    
! completed by Mitchell Brown, 03/20/2018
!********************************************************************    
    use diag_lib,   only: diag_print_error, diag_print_message
    use diag_def,   only: debug_mode
    use flow_def,   only: eta,u,v,h,p,iwet,grav,gravinv,flux
    use sed_def,    only: sedtrans,zb1,nsed,ctk,db,pbk
    use sal_def,    only: saltrans,sal
    use heat_def,   only: heattrans,heat
    use geo_def,    only: zb
    use in_def,     only: scaldattype,vecdattype
    use in_lib,     only: read_dat
    use hot_def,    only: icpres,icwse,icvel,icwet,icflux,ictime
    use size_def,   only: ncellsfull,ncellpoly,ncellsD,ncells
    use prec_def,   only: ikind
    use interp_lib, only: interp_scal_node2cell
    use comvarbl,   only: timehrs    
    use unitconv_lib, only: unitconv_scal    
    
    implicit none
    character(len=*), intent(in) :: supfile
    character(len=200) :: cardname, datfile,astring
    character :: ictimeunits*10
    integer :: kunit,ierr,istart,iend,imid,ival,ks,j, i,k,nd
    logical :: foundfile, founddataset
    logical :: icsingle = .false., icmulti = .false.
    integer :: nscal,nvec
    type(scaldattype), pointer :: scaldat(:)
    type(vecdattype),  pointer :: vecdat(:)
    real(4) :: etemp(ncellsD),utemp(ncellsD),vtemp(ncellsD),wtemp(ncellsD), umax, vmax
    real(ikind) :: temp(ncellsD)
    
    kunit = 500
    ierr = 0
    inquire(FILE=trim(supfile),exist=foundfile)
    if (foundfile) then
      open(kunit,file=trim(supfile),status='OLD')
    else
      call diag_print_error('Could not open file: ',trim(supfile))
      ierr=-1
    endif
    
    write(*,*) 'Starting Hot Start'

    if(debug_mode)then   !Set number of cells to read from files.  Same as in write routines.
      nd = ncellsD
    else
      nd = ncells
    endif

    do while(ierr==0)
      read(kunit,*,iostat=ierr) cardname
      if(ierr==-1) exit !End of File
      
      if(cardname(1:4).eq.'DATA') then
        backspace(kunit)
        read(kunit,*) cardname, datfile
        istart=index(trim(datfile),'_',.true.)+1  !from the back end, look for '_' and get the index of that location
        iend  =index(trim(datfile),'.',.true.)-1  !from the back end, look for '.' and get the index of that location
        astring=datfile(istart:iend)
        
        select case(astring(1:4))
        case ('eta ')   !Water Elevation
          call read_dat (datfile,nscal,scaldat,nvec,vecdat)
          if(scaldat(1)%nd /= nd)then
            call diag_print_error('Invalid initial conditions size','  Size of dataset does not match grid.')
          endif
          if(nscal.ne.1) call diag_print_error('Problem reading Water_Elevation data from '//trim(datfile))
          
          icwse = .true.
          call diag_print_message ('   Read Initial Water Level:            '//trim(datfile))
          
          do i=1,nd
            etemp(i) = scaldat(1)%val(i,scaldat(1)%nt)
          enddo
          do i=nd+1,ncellsd
            etemp(i) = 0.0
          enddo

          eta = etemp
   !Reading input for "ncellsD" (Active) cells already.  Removing this section.  meb  02/21/2019
          !if(ncellpoly>0)then
          !  call interp_scal_node2cell(etemp,eta) !Interpolate node to cell centers
          !else
          !  call map_scal_full2active(etemp,eta) !Convert from full to active grid 
          !endif
          
          !write(*,*) 'IC - Before P assignment, p(1) = ',p(1)
       !section added 07/06/2018  meb
          if (.not.icpres) then 
            p = eta*grav           !don't overwrite the pressures if already read in.  
          endif
          !write(*,*) 'IC - After P assignment,  p(1) = ',p(1)
          if(ictime < 0) then
            ictime = scaldat(1)%time(scaldat(1)%nt)  !Modify start time to match the initial conditions file.
            call unitconv_scal(scaldat(1)%time_units,'hrs',ictime)          
          elseif (ictime .ne. scaldat(1)%time(scaldat(1)%nt)) then
            call diag_print_error('Differing times in Hot Start files','Cannot continue')
          endif
       !Deallocate after
          deallocate(scaldat)
           
        case ('vel ')   !Current Velocity
          call read_dat (datfile,nscal,scaldat,nvec,vecdat)
          if(vecdat(1)%nd /= nd)then
            call diag_print_error('Invalid initial conditions size','  Size of dataset does not match grid.')
          endif
          if(nvec.ne.1) call diag_print_error('Problem reading Current_Velocity data from '//trim(datfile))
          
          icvel = .true.
          call diag_print_message ('   Read Initial Current Velocities:     '//trim(datfile))

          !Velocities
          do i=1,nd                                      
            utemp(i) = vecdat(1)%val(i,1,vecdat(1)%nt)
            vtemp(i) = vecdat(1)%val(i,2,vecdat(1)%nt)        
          enddo
          do i=nd+1,ncellsd   !If nd==ncellsD, do nothing more
            !etemp(i) = 0.0   !removed 03/25/21  Copied from scalar, but not correctly modified for vector  MEB
            utemp(i) = 0.0    !added
            vtemp(i) = 0.0    !added
          enddo
          
          u=utemp
          v=vtemp

          !umax = maxval(u)   !added for testing purposes  03/25/2021  MEB
          !vmax = maxval(v)
          
          !if(ncellpoly>0)then
          !  call interp_scal_node2cell(utemp,u) !Interpolate node to cell centers
          !  call interp_scal_node2cell(vtemp,v) !Interpolate node to cell centers
          !else
          !  call map_scal_full2active(utemp,u) !Convert from full to active grid 
          !  call map_scal_full2active(vtemp,v) !Convert from full to active grid 
          !endif

          if(ictime < 0) then
            ictime = vecdat(1)%time(vecdat(1)%nt)  !Modify start time to match the initial conditions file.
            call unitconv_scal(vecdat(1)%time_units,'hrs',ictime)          
          elseif (ictime .ne. vecdat(1)%time(vecdat(1)%nt)) then
            call diag_print_error('Differing times in Hot Start files','Cannot continue')
          endif
          deallocate(vecdat)

        case ('p   ')   !Water Pressure
          call read_dat (datfile,nscal,scaldat,nvec,vecdat)
          if(.not. icwse) then     !Don't overwrite eta, if it was already read in.
            if(scaldat(1)%nd /= nd)then
              call diag_print_error('Invalid initial conditions size','  Size of dataset does not match grid.')
            endif
            if(nscal.ne.1) call diag_print_error('Problem reading Water_Pressure data from '//trim(datfile))
          
            icpres = .true.
            call diag_print_message ('   Read Initial Water Pressure :        '//trim(datfile))
            do i=1,nd
              etemp(i) = scaldat(1)%val(i,scaldat(1)%nt)
            enddo
            do i=nd+1,ncellsd
              etemp(i) = 0.0
            enddo

            p=etemp
            !if(ncellpoly>0)then
            !  call interp_scal_node2cell(etemp,p) !Interpolate node to cell centers
            !else
            !  call map_scal_full2active(etemp,p) !Convert from full to active grid 
            !endif
            
            eta = p*gravinv
          endif  

          if(ictime < 0) then
            ictime = scaldat(1)%time(scaldat(1)%nt)  !Modify start time to match the initial conditions file.
            call unitconv_scal(scaldat(1)%time_units,'hrs',ictime)          
          elseif (ictime .ne. scaldat(1)%time(scaldat(1)%nt)) then
            call diag_print_error('Differing times in Hot Start files','Cannot continue')
          endif
          deallocate(scaldat)
            
        case ('wet ')   !Wet/Dry
          call read_dat (datfile,nscal,scaldat,nvec,vecdat)
          if(scaldat(1)%nd /= nd)then
            call diag_print_error('Invalid initial conditions size','  Size of dataset does not match grid.')
          endif
          if(nscal.ne.1) call diag_print_error('Problem reading Wet/Dry states from '//trim(datfile))
          
          icwet = .true.
          call diag_print_message ('   Read wet/dry states:                 '//trim(datfile))

          do i=1,nd
            etemp(i) = scaldat(1)%val(i,scaldat(1)%nt)
          enddo
          do i=nd+1,ncellsd
            etemp(i) = 0.0
          enddo

          temp = etemp
          !if(ncellpoly>0)then
          !  call interp_scal_node2cell(etemp,temp) !Interpolate node to cell centers
          !else
          !  call map_scal_full2active(etemp,temp) !Convert from full to active grid 
          !endif
          
          do i=1,nd
            iwet(i) = int(temp(i))
          enddo
            
          deallocate(scaldat)
          
        case ('Flux')   !Fluxes                           - Multiple possible
          call read_dat (datfile,nscal,scaldat,nvec,vecdat)
          if(scaldat(1)%nd /= nd)then
            call diag_print_error('Invalid initial conditions size','  Size of dataset does not match grid.')
          endif
          if(nscal.ne.1) call diag_print_error('Problem reading cell-face fluxes from '//trim(datfile))
          
          if (.not. icflux) call diag_print_message ('   Read cell-face fluxes:               '//trim(datfile))
          icflux = .true. !Switch to true so this only prints once
          
          do i=1,nd
            etemp(i) = scaldat(1)%val(i,scaldat(1)%nt)
          enddo
          do i=nd+1,ncellsd
            etemp(i) = 0.0
          enddo
          
          temp = etemp
          !if(ncellpoly>0)then
          !  call interp_scal_node2cell(etemp,temp) !Interpolate node to cell centers
          !else
          !  call map_scal_full2active(etemp,temp) !Convert from full to active grid 
          !endif
          
          !Should be more than one face.  Determine which one this is
          istart = index(datfile,'Flux')+4
          iend   = index(datfile,'.',.true.)-1
          astring=datfile(istart:iend)
          read(astring,'(I1)') ival

          !Save values to the appropriate dimension of the array
          flux(ival,:)=temp          
                
        case ('dept')   !Water Depth (at hotstart time)
          if(sedtrans) then
            call read_dat (datfile,nscal,scaldat,nvec,vecdat)
            if(scaldat(1)%nd /= nd)then
              call diag_print_error('Invalid initial conditions size','  Size of dataset does not match grid.')
            endif
            if(nscal.ne.1) call diag_print_error('Problem reading Depth from '//trim(datfile))
          
            icwse = .true.
            call diag_print_message ('   Read Initial Water Depths:           '//trim(datfile))
            
            do i=1,nd
              etemp(i) = scaldat(1)%val(i,scaldat(1)%nt)
            enddo
            do i=nd+1,ncellsd
              etemp(i) = 0.0
            enddo
            
            zb=etemp
            !if(ncellpoly>0)then
            !  call interp_scal_node2cell(etemp,zb) !Interpolate node to cell centers
            !else
            !  call map_scal_full2active(etemp,zb) !Convert from full to active grid 
            !endif
            
            deallocate(scaldat)
            zb = -zb
            zb1 = zb
          endif  
              
        case ('conc')   !Sediment Concentration           - Multiple possible
          if(sedtrans) then
            call read_dat (datfile,nscal,scaldat,nvec,vecdat)
            if(scaldat(1)%nd /= nd)then
              call diag_print_error('Invalid initial conditions size','  Size of dataset does not match grid.')
            endif
            if(nscal.ne.1) call diag_print_error('Problem reading Initial Sediment Concentrations from '//trim(datfile))
            
            call diag_print_message ('   Read Initial Sediment Concentration: '//trim(datfile))

            do i=1,nd
              etemp(i) = scaldat(1)%val(i,scaldat(1)%nt)
            enddo
            do i=nd+1,ncellsd
              etemp(i) = 0.0
            enddo
            
            temp = etemp
            !if(ncellpoly>0)then
            !  call interp_scal_node2cell(etemp,temp) !Interpolate node to cell centers
            !else
            !  call map_scal_full2active(etemp,temp) !Convert from full to active grid 
            !endif
            deallocate(scaldat)
          
            if(nsed==1 .and. .not.icsingle) then  !Single Grain Size
              Ctk(:,1)=temp
              icsingle=.true.
            else              !Multiple Grain Sizes
              !Determine which concentration file this is
              istart = index(datfile,'conc')+4
              iend   = index(datfile,'.',.true.)-1
              astring=datfile(istart:iend)
              read(astring,'(I1)') ks

              !Save values to the appropriate dimension of the array
              Ctk(:,ks)=temp          
              icmulti=.true.
            endif
          endif  

        case ('thic')   !Layer thickness                  - Multiple possible
          if(sedtrans) then
            call read_dat (datfile,nscal,scaldat,nvec,vecdat)
            if(scaldat(1)%nd /= nd)then
              call diag_print_error('Invalid initial conditions size','  Size of dataset does not match grid.')
            endif
            if(nscal.ne.1) call diag_print_error('Problem reading Initial Layer Thicknesses from '//trim(datfile))
            
            call diag_print_message ('   Read Initial Layer Thicknesses:      '//trim(datfile))

            do i=1,nd
              etemp(i) = scaldat(1)%val(i,scaldat(1)%nt)
            enddo
            do i=nd+1,ncellsd
              etemp(i) = 0.0
            enddo
            
            temp = etemp
            !if(ncellpoly>0)then
            !  call interp_scal_node2cell(etemp,temp) !Interpolate node to cell centers
            !else
            !  call map_scal_full2active(etemp,temp) !Convert from full to active grid 
            !endif
            deallocate(scaldat)
          
            !Multiple Grain Sizes
            !Determine which thickness file this is
            istart = index(datfile,'thick')+5
            iend   = index(datfile,'.',.true.)-1
            astring=datfile(istart:iend)
            ival = len_trim(astring)
            select case(ival)
            case (1)
              read(astring,'(I1)') j
            case (2)
              read(astring,'(I2)') j
            end select

            !Save values to the appropriate dimension of the array
            db(:,j)=temp          
          endif  
            
        case ('frac')   !Fraction of grain size per layer - Multiple possible
          if(sedtrans) then
            call read_dat (datfile,nscal,scaldat,nvec,vecdat)
            if(scaldat(1)%nd /= nd)then
              call diag_print_error('Invalid initial conditions size','  Size of dataset does not match grid.')
            endif
            if(nscal.ne.1) call diag_print_error('Problem reading Bed Composition from '//trim(datfile))
            
            call diag_print_message ('   Read Initial Bed Composition:        '//trim(datfile))

            do i=1,nd
              etemp(i) = scaldat(1)%val(i,scaldat(1)%nt)
            enddo
            do i=nd+1,ncellsd
              etemp(i) = 0.0
            enddo
            
            temp = etemp
            !if(ncellpoly>0)then
            !  call interp_scal_node2cell(etemp,temp) !Interpolate node to cell centers
            !else
            !  call map_scal_full2active(etemp,temp) !Convert from full to active grid 
            !endif
            deallocate(scaldat)
          
            !Multiple Grain Sizes
            !Determine which bed composition file this is
            istart = index(datfile,'frac')+4
            imid   = index(datfile,'-',.true.)-1
            iend   = index(datfile,'.',.true.)-1
            astring=datfile(istart:imid)
            ival=len_trim(astring)
            select case (ival)
            case (1)
              read(astring,'(I1)') j   !Layer number
            case (2)
              read(astring,'(I2)') j   !Layer number
            end select

            astring=datfile(imid+2:iend)
            ival=len_trim(astring)
            select case (ival)
            case (1)
              read(astring,'(I1)') ks  !Size class number
            case (2)
              read(astring,'(I2)') ks  !Size class number
            end select
          
            !Save values to the appropriate dimension of the array
            pbk(:,ks,j)=temp
          endif  

        case ('sal ')   !Salinity Concentration
          if(saltrans) then
            call read_dat (datfile,nscal,scaldat,nvec,vecdat)
            if(scaldat(1)%nd /= nd)then
              call diag_print_error('Invalid initial conditions size','  Size of dataset does not match grid.')
            endif
            if(nscal.ne.1) call diag_print_error('Problem reading Salinity Concentration from '//trim(datfile))
          
            icwse = .true.
            call diag_print_message ('   Read Initial Salinity Concentration: '//trim(datfile))

            do i=1,nd
              etemp(i) = scaldat(1)%val(i,scaldat(1)%nt)
            enddo
            do i=nd+1,ncellsd
              etemp(i) = 0.0
            enddo
            
            sal = etemp
            !if(ncellpoly>0)then
            !  call interp_scal_node2cell(etemp,sal) !Interpolate node to cell centers
            !else
            !  call map_scal_full2active(etemp,sal) !Convert from full to active grid 
            !endif
            deallocate(scaldat)
          endif  

        case ('heat')   !Temperature
          if(heattrans) then
            call read_dat (datfile,nscal,scaldat,nvec,vecdat)
            if(scaldat(1)%nd /= nd)then
              call diag_print_error('Invalid initial conditions size','  Size of dataset does not match grid.')
            endif
            if(nscal.ne.1) call diag_print_error('Problem reading Temperature from '//trim(datfile))
          
            icwse = .true.
            call diag_print_message ('   Read Initial Temperatures:           '//trim(datfile))

            do i=1,nd
              etemp(i) = scaldat(1)%val(i,scaldat(1)%nt)
            enddo
            do i=nd+1,ncellsd
              etemp(i) = 0.0
            enddo
            
            heat = etemp
            !if(ncellpoly>0)then
            !  call interp_scal_node2cell(etemp,heat) !Interpolate node to cell centers
            !else
            !  call map_scal_full2active(etemp,heat) !Convert from full to active grid 
            !endif
            deallocate(scaldat)
          endif  

        end select
      endif
     
    enddo 
    close(kunit)
    
    return
    end subroutine hot_read_sup
    
!********************************************************************
    subroutine hot_read_dat()
! Reads initial condition file from an SMS ASCII dataset (*.dat) file 
! written by Alex Sanchez, USACE-CHL    
! completed by Mitchell Brown, 03/20/2018
!********************************************************************    
    use size_def
    use diag_lib
    use hot_def
    use in_def
    use flow_def,     only: eta,u,v,h,p,iwet,grav
    use in_lib,       only: read_dat
    use comvarbl,     only: timehrs
    use interp_lib,   only: interp_scal_node2cell
    use unitconv_lib, only: unitconv_scal
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
      select case(scaldat(i)%name)
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
      end select
    enddo
    
    !--- Current Velocities -----
    do i=1,nvec
      select case(vecdat(i)%name)
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
      end select
    enddo
    
    return
    end subroutine hot_read_dat

!******************************************************************    
    subroutine hot_read_xmdf()
! Reads hot start information
! written by Alex Sanchez, USACE-CHL    
! modified by Mitch Brown - 05/18/2012
!******************************************************************
#include "CMS_cpp.h"
#ifdef XMDF_IO
    use size_def
    use hot_def
    use xmdf
    use diag_lib
    use prec_def
    use geo_def,  only: cell2cell,zb,zb0
    use flow_def, only: u,v,u1,v1,uv,p,p1,eta,h,iwet,flux,grav,gravinv
    use sed_def,  only: sedtrans,nsed,Ctk,Ctk1,Ctkstar,icapac,nlay,bedlay,pbk,pbk1
    use sed_def,  only: db,db1,d50,d90,diam,diamlim,logdiamlim,dbmax,dmconst,zb1
    use sal_def,  only: saltrans,sal
    use heat_def, only: heattrans, heat
    use out_def,  only: outlist,simlabel
    use DFLIB,    only: systemqq
    use in_xmdf_lib, only: readscallasth5,readveclasth5
    
    implicit none
    integer :: res !Added MEB - 05/18/2012 for hot start fixes
    integer :: i,j,k,ks,fid,gid,nn,ierr,ntimes,npath,fileCopied
    real(ikind) :: temphr,temp(ncellsD)
    real(8), allocatable :: timesd(:)
    real(ikind), allocatable :: d35(:)
    !real(ikind) :: var(ncellsD)
    character(len=200) :: apath,thepath,thename,tfile !Added tfile MEB - 05/18/2012 for hot start fixes
    character(len=5) :: apbk,alay
    character(len=10) :: theext 
    logical :: ok
    character(len=30) :: wse_names(5)
    character(len=300) :: msg,astring
    
    tfile=trim(icfile)   
    call fileparts(icfile,thepath,thename,theext)
    
!****************** MEB - 05/18/2012
! XMDF error keeps ICFILE, AutoHotStart.h5 or SingleHotStart.h5, from being deleted after reading

    ! Check if temp.h5 file exists, if so, delete it
    tfile=trim(thepath)//'temp.h5'
    inquire(file=tfile, exist=ok)
    if(ok)then
      open(100,file=tfile,iostat=ierr)
      close(100,status='delete',iostat=ierr)
    endif

    ! Make a copy of ICFILE to temp.h5 and read from that. 
!This may need to change for linux MEB 10/16/2018
    !res = system('copy '//trim(icfile)//' '//trim(tfile)//' > file')
    astring = 'copy "' // trim(icfile) // '" "' // trim(tfile) // '" > file'
    res = systemqq(trim(astring))

    inquire(file='file', exist=ok)
    !Checking for existence of 'file' is not good enough.  
    !Read 'file' and look for the text 'copied' on Windows
    !Modify for appropriate message on Linux
    if(ok)then            !remove scratch file
      open(100,file='file',iostat=ierr)
      read(100,'(A300)') astring
      fileCopied=index(astring,'opied')
      if (filecopied .gt. 0) then 
        close(100,status='delete',iostat=ierr)
      else
        call diag_print_error(' ','*** Error copying IC file to temporary file ***','See "file" for message')
      endif
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
      !ictime = temphr !Start time in hours 
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
         write(apath,'(A,A4,I2)') icpath(1:npath),'Flux',k    !If over 9, then there are two digits to the Flux - meb 03/20/2018
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
            apath = icpath(1:npath) //'Fraction' // trim(apbk) // trim(alay)
            call readscallasth5(tfile,apath,ntimes,pbk(:,ks,j),reftimehot,temphr,ierr)            
            if(ierr<0)then                                                        !Fix for missing space or underscore  MEB  03/15/2021
              apath = icpath(1:npath) // 'Fraction' // trim(apbk) // ' ' // trim(alay)
              call readscallasth5(tfile,apath,ntimes,db(:,j),reftimehot,temphr,ierr)
              if(ierr<0)then 
                apath = icpath(1:npath) // 'Fraction' // trim(apbk) // '_' // trim(alay)
                call readscallasth5(tfile,apath,ntimes,db(:,j),reftimehot,temphr,ierr)
              endif            
            endif
            if(ierr<0)then
              call diag_print_warning('Problem reading bed composition',&
               '  Setting the bed composition to the default values')  
              ok = .false.
              exit loopj
            else
              write(msg,'(A,A)') 'Read bed composition: ',trim(apath)   !Copy/paste error.  This still printed "bed layer thickness" instead of "bed composition"  MEB  03/15/2021
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
            apath = icpath(1:npath) // 'Concentration' // apbk                         !The underscore is already added with the '62' format  MEB  03/15/2021
            call readscallasth5(tfile,apath,ntimes,Ctk(:,ks),reftimehot,temphr,ierr)
            if(ierr<0)then
              call diag_print_warning('Unable to find initial sediment concentrations',&
                '   Setting initial sediment concentrations to equilibrium concentration')
              setconc2eq = .true. 
              exit
            else
              write(msg,'(A,A)') 'Read sediment concentration: ',trim(apath)
              call diag_print_message(msg)
            endif 
          enddo
        else
          write(msg,'(A,A)') 'Read sediment concentration: ',trim(apath)
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
      call readscallasth5(tfile,apath,ntimes,heat,reftimehot,temphr,ierr)  !Fixed 03/20/2018 - changed 'sal' to 'heat'
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
    
!Test to see if file exists in current working directory
    call XF_OPEN_FILE(trim(outlist(1)%afile),READONLY,fid,ierr)             
    
    if (ierr<0)then        !No solution file exists in the current, so use the time from the initial condition.
      call diag_print_message('','Hot Start: Previous solution file(s) not found. New solution file(s) will be created.')
      timeout = ictime
    else                   !Otherwise a solution file was found, obtain the last time from that file.
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
    end subroutine hot_read_xmdf
    
!******************************************************************    
    subroutine hot_init()
! Initializes any missing variables including Ctk and pbk
! and takes care of boundaries
! written by Alex Sanchez, USACE-CHL
!
! only called when hot start is to be read   !MEB  02/24/2016
!******************************************************************    
#include "CMS_cpp.h"
    use size_def
    use bnd_def
    use hot_def
    use sal_def
    use heat_def
    use diag_lib
    use geo_def,  only: cell2cell,zb
    use flow_def, only: u,v,u1,v1,uv,p,p1,eta,h,h1,hmin,iwet,iwet1,ponding,grav,gravinv
    use comvarbl, only: ctime,stime,timehrs,timesecs
    use sed_def,  only: nsed,icapac,Ctk,Ctk1,Ctkstar,CtstarP,pbk,singlesize
    use sed_def,  only: variableD50,d50,diam,mhe,iHidExpForm,varsigma
    use stat_def, only: calc_stats,tstat
    use out_def,  only: simlabel
    implicit none
    
    integer :: i,j,k,ibnd,nck,iwse,isal,iheat
    character(len=100) :: msg2,msg3
            
      
    !--- Hot Start Time ---------------------------------------------
755 format(1x,A,F10.3,A)    
    !Overwrites time in icfile. The time in icfile may be incorrectly set to zero by SMS
    if(ictime>-1.0e-20)then !if specified
      timehrs = ictime !Initial condition time
    else 
      if(timeout>-1.0e-20)then !if specified
        timehrs = timeout !Used time in solution file
      else  
        timehrs = 0.0  !no initial condition time estimate available so set to default
      endif  
    endif
    stime = timehrs*3600.0     !stime is elapsed time in seconds
    ctime = stime
    timesecs = dble(ctime)
    hstarttime = timehrs
    
    write(msg2,755) '  Hot Start Time:   ',timehrs,' hours'
    call diag_print_message(' ',msg2)   
    
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
    
    !---- Apply wet/dry ------------
    !Dry nodes     
    do i=1,ncellsD
      !write(*,*) 'Hot Init - p(i) ',p(i),i
      h(i) = p(i)*gravinv - zb(i)
      if(h(i)<hmin+1.0e-5)then    !if dry
        h(i) = hmin
        iwet(i) = 0  
        u(i) = 0.0
        v(i) = 0.0  
        eta(i) = -999.0
        if(.not.icpres)then
          p(i) = (h(i)+zb(i)-1.0e-4)*grav
        endif
      else                        !if wet
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
      select case(iHidExpForm) !varsigma(i,k)
      case(1); call HidExpEgiazaroff    !Egiazaroff (1965)
      case(2); call HidExpParker        !Parker et al. (1982) and others
      case(3); call HidExpWu            !Wu et al. (2000)
      case(4); call HidExpAshidaMichiue !Ashida and Michiue 1980
      case(5); call HidExpHayashi       !Hayashi et al. 1980
      end select
      
        !Transport Capacity
      select case(icapac)  
      case(1); call sedcapac_lundcirp !Lund-CIRP          
      case(2); call sedcapac_vanrijn  !Van Rijn           
      case(3); call sedcapac_watanabe !Watanabe  
      case(4); call sedcapac_soulsby  !Soulsby (1997)
      case(5); call wucapac           !Wu et al. 2000 (under testing)
      case(6); call sedcapac_c2shore   !C2SHORE (bdj)
      end select   
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
    end subroutine hot_init
    
!******************************************************************
    subroutine hot_write
! Writes hot start h5 file for restarting CMS
! written by Alex Sanchez, USACE-CHL
!******************************************************************
#include "CMS_cpp.h"
    use diag_def
    use diag_lib
    use prec_def
    use hot_def,  only: hottime,hot_timehr,hot_recur,hotdt,hotname,autohotname,hotfile,hotpath,autohotfile,autohotpath,autohot_inc
    use comvarbl, only: timehrs,casename
    use out_def,  only: write_sup
#ifdef _WIN32    
    use IFPORT
#endif
    
    implicit none
    real(ikind) :: eps
    logical :: iSingleHot,iAutoHot,found,created
    character(len=200) :: hotdirpath, aname
    
    eps = 1.e-10    
    iSingleHot = .false. ; iAutoHot = .false.
    if(hot_timehr .and. abs(timehrs-hottime)<=eps) iSingleHot = .true.   
    if(hotdt.gt.0.0) then 
      if(hot_recur  .and. mod(timehrs*60,hotdt*60)<=eps)   iAutoHot   = .true.   !multiplied by 60 so it would output correctly with smaller increments. MEB 08/19/2024
    endif  
    if(.not.iSingleHot .and. .not.iAutoHot) return              
    hotdirpath='ASCII_HotStart'

    if(write_sup)then
      if (iSingleHot) then
        aname = trim(hotdirpath)//'/'//hotname
        write(msg,'(A,F13.4,A)') 'Writing Single Hotstart File at:    ',timehrs,' hrs'
        call diag_print_message(' ',msg)
        call write_hotstart_sup(hotdirpath,aname)  !Super File *.sup and individual ASCII Dataset *.dat files
      endif
      if (iAutoHot)   then
        aname = trim(hotdirpath)//'/'//autohotname
        write(msg,'(A,F13.4,A)') 'Writing Recurring Hotstart File at: ',timehrs,' hrs'
        call diag_print_message(' ',msg)
        call write_hotstart_sup(hotdirpath,aname)  !Super File *.sup and individual ASCII Dataset *.dat files
      endif
    else
#ifdef XMDF_IO
      if (iSingleHot) then
        write(msg,'(A,F13.4,A)') 'Writing Single Hotstart File at:    ',timehrs,' hrs'
        call diag_print_message(' ',msg)
        call hot_write_xmdf (hotfile,hotpath)         !XMDF single hot file
      endif
      if (iAutoHot)   then
        write(msg,'(A,F13.4,A)') 'Writing '//autohotfile(1:14)//' File at: ',timehrs,' hrs'
        call diag_print_message(' ',msg)
        call hot_write_xmdf (autohotfile,autohotpath) !XMDF recurring hot file
      endif
#endif
    endif

    !added to write out two files - MEB 08/19/2024
    if (autohot_inc == '_1') then              !First one just written, renumber to second.
      autohot_inc = '_2'
      autohotfile = autohotfile(1:13)//'2.h5'
    else                                       !Second one just written, renumber back to first.
      autohot_inc = '_1'
      autohotfile = autohotfile(1:13)//'1.h5'
    endif
        
    return
    end subroutine hot_write

!***********************************************************************    
    subroutine hotstart_file_init
! Initializes the ASCII hotstart files
! written by Mitchell Brown, USACE-CHL - 03/16/18
! update 06/07/19 - If not hot start, you will not even get into this routine.  
!***********************************************************************
    use out_lib,  only: write_xy_file
    use comvarbl, only: casename
    use diag_def, only: debug_mode
    use diag_lib, only: diag_print_error
    use hot_def,  only: autohotfile, hotfile, autohotname, hotname, hot_out, hot_timehr, hot_recur, coldstart,icfile,autohot_inc
#ifdef _WIN32
    use ifport,   only: makedirqq
#endif
    implicit none
    
    character(len=200) :: hotdirpath, aname, filesup, filexy, tfile, astring
    logical            :: iSingleHot, iAutoHot, found, created, res, ok
    integer            :: nunit,ierr,filecopied

    iSingleHot = .false. ; iAutoHot = .false.
    if(hot_timehr) iSingleHot = .true.        !Changed from 'hot_out'.  'hot_timehr' governs Single time hot start
    if(hot_recur)  iAutoHot   = .true.   
    
    !Save all these files to a subdirectory named "ASCII_HotStart"
    hotdirpath='ASCII_HotStart'

    !Change from XMDF filename 
    autohotfile = trim(AutoHotName)//'.sup'   !keep outputting only 1 ASCII recurring hotstart file for now.
    hotfile     = trim(HotName)//'.sup'
    tfile       = trim(hotdirpath)//'/InitialCondition.sup'
#ifdef _WIN32    
    inquire(directory=trim(hotdirpath), exist=found)
    if(.not.found) then
      created=MakeDirQQ(trim(hotdirpath))
      if(.not.created)then
        call diag_print_error('Failed to create subdirectory- '//trim(hotdirpath))
      endif
    endif
#else
    inquire(file=trim(hotdirpath), exist=found)
    if (.not.found) then 
      call system('mkdir '//(trim(hotdirpath)))
    endif
#endif
    !If this is a hotstart, write a copy of the file to a new name (long way because of the difficulty of '\' vs '/' in file paths)
    if(.not.coldstart) then
      open(100,file=trim(icfile))
      open(101,file=trim(tfile))
      do
        read(100,'(A)',iostat=ierr) astring
        if(ierr<0) exit
        write(101,'(A)') trim(astring)
      enddo
      close(100)
      close(101)
      !Also reset the name of the ICFILE to the new name for reading later
      icfile=trim(tfile)
    else
      inquire(file=trim(tfile),exist=found)
      if (found) then
        open(100,file=trim(tfile))
        close(100,status='delete')
      endif
    endif

101 format('SUPER')
102 format('SCAT2D  "',A,'"')       

    if (iSingleHot) then
      aname = trim(hotdirpath)//'/'//hotname
      call write_xy_file(aname,casename)         !XY coordinate file *.xy  
      filesup = trim(aname) // '.sup'    
      filexy  = trim(aname) // '.xy'
      nunit = 46
      open(nunit,file=filesup)
      write(nunit,101)
      write(nunit,102) trim(filexy)    
      close(nunit)
    endif  
    if (iAutoHot)   then
      aname = trim(hotdirpath)//'/'//autohotname
      call write_xy_file(aname,casename)         !XY coordinate file *.xy  
      filesup = trim(aname) // '.sup'    
      filexy  = trim(aname) // '.xy'
      nunit = 46
      open(nunit,file=filesup)
      write(nunit,101)
      write(nunit,102) trim(filexy)    
      close(nunit)
    endif

    return
    end subroutine hotstart_file_init

!******************************************************************
    subroutine hot_write_xmdf(outfile,outpath)
! Writes hot start h5 file for restarting CMS
! written by Alex Sanchez, USACE-CHL
!******************************************************************
#include "CMS_cpp.h"
#ifdef XMDF_IO
    use size_def
    use prec_def
    use xmdf
    use hot_def,  only: hotfile,hottime,hot_timehr,hotdt,hot_recur,icfile,hotpath,autohot_inc
    use geo_def,  only: zb
    use flow_def, only: u,v,u1,v1,p,p1,eta,iwet,iwet1,flux,flux1,grav
    use comvarbl, only: ctime,timehrs,flowpath,ntsch
    use sed_def,  only: sedtrans,nsed,Ctk,nlay,pbk,db
    use sal_def,  only: saltrans,sal
    use heat_def, only: heattrans, heat
    use stat_def, only: flowstats,sedstats,salstats,heatstats
    use out_def,  only: simlabel
    use out_lib,  only: writescalh5,writevech5
    implicit none
    
    character(len=*), intent(in) :: outfile, outpath
    integer :: j,ks
#ifdef DEV_MODE
    integer :: i,k
#endif
    real(ikind) :: var(ncellsD)
    character(len=100) :: aname
    character(len=5) :: apbk,alay
    
    !delete previous hotstart file so there is only one record 
    open(100,file=outfile)
    close(100,status='DELETE')  
  
    !--- Hydrodynamics ----
    call writescalh5(outfile,outpath,'Water_Pressure',p,'m^2/s^2',timehrs,1)   
    call writescalh5(outfile,outpath,'Water_Elevation',eta,'m',timehrs,1)
    call writevech5(outfile,outpath,'Current_Velocity',u,v,'m/s',timehrs,1)
     
#ifdef DEV_MODE    
    !Wet/dry
    var = iwet
    call writescalh5(outfile,outpath,'Wet',var,'',timehrs,1)    
    
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
      call writescalh5(outfile,outpath,aname,var,'m^3/s',timehrs,1)
    enddo
#endif

62  format('_',I2.2)
71  format(1x,'(',I1,')')
72  format(1x,'(',I2,')')
    
    !--- Sediment Transport ----
    if(sedtrans)then      
      !Water depths
      var = -zb
      call writescalh5(outfile,outpath,'Depth',var,'m',timehrs,1)  
!      !Initial Water depths
      !Sediment concentrations
      if(nsed==1)then !Single grain size
        aname = 'Concentration'
        call writescalh5(outfile,outpath,aname,Ctk(:,1),'kg/m^3',timehrs,1)  
      else            !Multiple grain size
        !Concentrations
        do ks=1,nsed
          write(apbk,62) ks
          aname = 'Concentration'// apbk
          call writescalh5(outfile,outpath,aname,Ctk(:,ks),'kg/m^3',timehrs,1)  
        enddo
        
        !Bed material composition and layer thickness
        do j=1,nlay
          if(j<=9)then
            write(alay,71) j
          else
            write(alay,72) j
          endif  
          aname = 'Thickness' // alay !Bug Fix, changed apath to aname
          call writescalh5(outfile,outpath,aname,db(:,j),'m',timehrs,1)  
          do ks=1,nsed
            write(apbk,62) ks
            aname = 'Fraction' // trim(apbk) // alay 
            call writescalh5(outfile,outpath,aname,pbk(:,ks,j),'none',timehrs,1)          
          enddo !j
        enddo !ks      
      endif
    endif
    
    !--- Salinity Transport ----
    if(saltrans)then
      call writescalh5(outfile,outpath,'Salinity',sal,'ppt',timehrs,1)  
    endif
    
    !--- Heat Transfer ----
    if(heattrans)then
      call writescalh5(outfile,outpath,'Heat',heat,'ppt',timehrs,1)  
    endif    
#endif

    return
    end subroutine hot_write_xmdf

!***********************************************************************    
    subroutine write_hotstart_sup (hotdirpath, aname)
! writes the SMS Super ASCII HotStart Files
! written by Mitchell Brown, USACE-CHL - 03/16/18
! - starting with individual hotstart files and will eventually write all into one.
!***********************************************************************
    use size_def
    use prec_def
    use hot_def,  only: hotfile,hottime,hot_timehr,hotdt,hot_recur,icfile,hotpath
    use geo_def,  only: zb
    use flow_def, only: u,v,u1,v1,p,p1,eta,iwet,iwet1,flux,flux1,grav
    use comvarbl, only: ctime,timehrs,flowpath,ntsch,casename
    use sed_def,  only: sedtrans,nsed,Ctk,nlay,pbk,db
    use sal_def,  only: saltrans,sal
    use heat_def, only: heattrans, heat
    use stat_def, only: flowstats,sedstats,salstats,heatstats
    use out_def,  only: simlabel
    use out_lib,  only: write_scal_dat_file,write_vec_dat_file
    use diag_lib, only: diag_print_error
    implicit none

    integer            :: i,j,k,ks,nunit,ierr
    real(ikind)        :: var(ncellsD)
    character(len=5)   :: apbk,alay    
    character(len=100) :: lname,sname,filesup,filexy,acmd
    character(len=200) :: apath,anarg,aline
    character(len=*), intent(in) :: aname,hotdirpath
    
101 format('SUPER')
102 format('SCAT2D  "',A,'"')   
    
    !delete previous hotstart DAT files, so there is only one record in each
    !!NOTE: when 'HotStart' is found in the filename, the previous file will now be deleted by the 
    !!'write_scal_dat_file' and 'write_vel_dat_file' routines.
    !!MEB - 03/16/2018
    
    !--- Hydrodynamics ----
    call write_scal_dat_file(aname,'Water_Pressure','p',p)   
    call write_scal_dat_file(aname,'Water_Elevation','eta',eta)
    call write_vec_dat_file (aname,'Current_Velocity','vel',u,v)

    !Wet/dry
    var = iwet
    call write_scal_dat_file(aname,'Wet','wet',var)
    
#ifdef DEV_MODE
    !Fluxes
    do k=1,nmaxfaces
      do i=1,ncellsD
        var(i) = flux(k,i)
      enddo
      if(k<=9)then
        write(lname,('A,I1)') 'Flux',k
      else
        write(lname,('A,I2)') 'Flux',k
      endif
      call write_scal_dat_file(aname,lname,lname,var)
#endif

62  format('_',I2.2)
71  format(1x,'(',I1,')')
72  format(1x,'(',I2,')')        

    !--- Sediment Transport ----
    if(sedtrans)then      
      !Water depths
      var = -zb
      call write_scal_dat_file(aname,'Depth','depth',var)

      !Sediment concentrations
      if(nsed==1)then !Single grain size
        call write_scal_dat_file(aname,'Concentration','conc',Ctk(:,1))
      else              !Multiple grain size
        !Concentrations
        do ks=1,nsed
          write(apbk,62) ks
          lname = 'Concentration'// apbk
          sname = 'conc'// apbk
          call write_scal_dat_file(aname,lname,sname,Ctk(:,ks))  
        enddo
        
        !Bed material composition and layer thickness
        do j=1,nlay
          if(j<=9)then
            write(alay,71) j
          else
            write(alay,72) j
          endif  
          lname = 'Thickness' // alay !Bug Fix, changed apath to aname
          sname = 'thick' // alay
          call write_scal_dat_file(aname,lname,sname,db(:,j))  
          do ks=1,nsed
            write(apbk,62) ks
            lname = 'Fraction' // trim(apbk) // alay
            sname = 'frac' // trim(apbk) // alay  
            call write_scal_dat_file(aname,lname,sname,pbk(:,ks,j))
          enddo !j
        enddo !ks      
      endif
    endif
        
    !--- Salinity Transport ----
    if(saltrans)then
      call write_scal_dat_file(aname,'Salinity','sal',sal)
    endif
    
    !--- Heat Transfer ----
    if(heattrans)then
      call write_scal_dat_file(aname,'Heat','heat',heat)
    endif    

    return
    end subroutine write_hotstart_sup
   
!***********************************************************************    
    subroutine hot_print()
! Prints the Hot Start settings to the screen and diagnostic file
! written by Alex Sanchez, USACE-CHL
!***********************************************************************    
    use hot_def
    use diag_def, only: dgunit,dgfile
    use tool_def, only: vstrlz
    implicit none
    
    integer :: i,iunit(2)
    
645 format(' ',A,T40,F0.2,A)
222 format(' ',A,T40,A,A)
223 format(' ',A)
    
    iunit = (/6, dgunit/)
    
    open(dgunit,file=dgfile,access='append') 
    do i=1,2    
      write(iunit(i),*)

      if(coldstart)then
        write(iunit(i),222)   'Start Mode:','COLD START'
      else
        write(iunit(i),222)   'Start Mode:','HOT START'
        write(iunit(i),222)   '  Initial Conditions File:',trim(icfile)
      endif

      if(hot_out)then
        if(hot_timehr)then 
          write(iunit(i),222) 'Single Hot Start Output: '
          write(iunit(i),222) '  File:',trim(hotfile)
          write(iunit(i),222) '  Time:',trim(vstrlz(hottime,'(f0.2)')),' hrs'
        endif
        if(hot_recur)then
          write(iunit(i),222) 'Recurring Hot Start Output: '
          write(iunit(i),222) '  File:',trim(autohotfile)
          write(iunit(i),222) '  Recurring Interval:',trim(vstrlz(hotdt,'(f0.2)')),' hrs'
        endif
        if(add_duration_HS)then
          write(iunit(i),222) 'Extend DURATION_RUN with IC Time:','ON'
          write(iunit(i),223) '  - Total time of simulation will be extended by the Time value found in the IC file.' 
        endif
      endif
    enddo
    close(dgunit)

    return
    end subroutine hot_print
    