!***********************************************************************
    subroutine prestart_EXP_shared
! by Weiming Wu, NCCHE, Oct. 2008
! modified by Alex Sanchez, USACE-CHL
!***********************************************************************
#include "CMS_cpp.h"
    use size_def, only: ncells
    use geo_def,  only: telfile, igridtype, grdfile, zb, dzbx, dzby, zbk
    use flow_def, only: iwet, iwet1, eta, p, gravinv, ponding, flowvolbal, h0, h, hmin
    use comvarbl, only: ctlfile, input_ver, advfile, read_adv, ctime, nthr
    use cms_def,  only: ncards, cardlist
    use sed_def,  only: sedtrans
    use sal_def,  only: saltrans
    use heat_def, only: heattrans
    use hot_def,  only: coldstart, hot_out, hot_timehr, hotfile, hot_recur, autohotfile, icflux, timeout
    use der_def,  only: goa
    use out_def,  only: write_sup, obs_cell, save_point
    use der_lib,  only: der_grad_eval
    use interp_lib, only: interp_scal_cell2face
    use dredge_def, only: dredging
    use diag_lib, only: diag_print_warning, diag_print_message, diag_print_error
#ifdef DEV_MODE
    use q3d_def
    use veg_def
#endif
#ifdef PROFILE
    use watch_lib
#endif

    implicit none
    integer :: k,ierr,icount,ilen,i
    character :: aLetter
    character(len=37) :: cardname,aext
    character(len=200) :: aname,apath,astring
    logical :: foundcard,foundfile
    
    !Get CMS Version number for setting defaults
    open(77,file=ctlfile,status='unknown')
    do
      read(77,*,iostat=ierr) cardname
      if(ierr/=0) exit
      if(cardname(1:11)=='CMS_VERSION')then !Get CMS Version
        backspace(77)
        read(77,*) cardname, input_ver      !MEB changed 'ver' to 'input_ver' 6/26/2016
        exit
      endif
    enddo
    close(77)
        
!--- Set defaults -------------------------
    write(*,*) 'Setting Module Defaults'
    call geo_default    !Geospatial
    call flow_default   !Hydrodynamics
    call solv_default   !Sparse matrix solvers
    call fric_default   !Friction
    call bnd_default    !Boundary conditions
    call struct_default !Structures
    call hot_default    !Hot start
    call sed_default    !Sediment transport 
    call sal_default    !Salinity
    call heat_default   !Temperature
    call met_default    !Meteorological
    call rol_default    !Roller
    call out_default    !Output
    call stat_default   !Statistics
    call der_default    !Derivatives
    call diag_default
    call wave_default
#ifdef DEV_MODE
    call q3d_default    !Quasi-3D 
    call veg_default    !Vegetation
#endif

    call dredge_default !dredge module

#ifdef PROFILE
    call watch_default  !Watches
#endif

!--- Read Input and Advanced Card Files -----------------------------   
    do icount=1,2
      if (icount == 1) then                              !First time through process from Input '*.cmcards' file
        call fileparts(ctlfile,apath,aname,aext) 
        astring = trim(aname) // '.' // trim(aext)
        write(*,*) 'Reading CMS-Flow Card File: ',trim(astring)
        open(77,file=ctlfile,status='unknown')
      else                                               !Second time through process from Advanced 'advanced.cmcards' file.  Note: this will overwrite previously read cards if they exist in the file.  
        inquire(file=advfile,exist=foundfile)
        if(.not.foundfile) exit  !no advanced cards to use
        call fileparts(advfile,apath,aname,aext) 
        astring = trim(aname) // '.' // trim(aext)
        write(*,*) 'Reading CMS-Flow Advanced Card File: ',trim(astring)
        open(77,file=advfile,status='unknown')
        read_adv = .true.
      endif
      do
        read(77,*,iostat=ierr) cardname
        if(ierr/=0) exit
        if(cardname(1:14)=='END_PARAMETERS') exit
        if(cardname(1:1)=='!' .or. cardname(1:1)=='#' .or. cardname(1:1)=='*') cycle
        call ignr_cards(cardname,foundcard);          if(foundcard) cycle
        call geo_cards(cardname,foundcard,.true.);    if(foundcard) cycle
        call flow_cards(cardname,foundcard,.true.);   if(foundcard) cycle
        call bnd_cards(cardname,foundcard);           if(foundcard) cycle
        call struct_cards(cardname,foundcard);        if(foundcard) cycle
        call fric_cards(cardname,foundcard);          if(foundcard) cycle
        call hot_cards(cardname,foundcard);           if(foundcard) cycle
        call sed_cards(cardname,foundcard);           if(foundcard) cycle
        call sal_cards(cardname,foundcard);           if(foundcard) cycle
        call heat_cards(cardname,foundcard);          if(foundcard) cycle
        call met_cards(cardname,foundcard);           if(foundcard) cycle
        call rol_cards(cardname,foundcard);           if(foundcard) cycle
        call out_cards(cardname,foundcard);           if(foundcard) cycle
        call diag_cards(cardname,foundcard);          if(foundcard) cycle
        call wave_cards(cardname,foundcard);          if(foundcard) cycle
#ifdef DEV_MODE
        call q3d_cards(cardname,foundcard);    if(foundcard) cycle
        call veg_cards(cardname,foundcard);    if(foundcard) cycle
#endif
        call dredge_cards(cardname,foundcard); if(foundcard) cycle
        call der_cards(cardname,foundcard);    if(foundcard) cycle
        call stat_cards(cardname,foundcard);   if(foundcard) cycle
        call balancecheck_cards(cardname,foundcard); if(foundcard) cycle
        if(.not.foundcard)then
          call addUnknownCard2List(77)
        endif
      enddo
      close(77)
    enddo  
    
    if(nCards > 0) then
      call diag_print_warning('Advanced cards found:')
      do i=1,nCards
        call diag_print_message('  '//trim(cardList(i)%cardname))
      enddo
!     write(*,*) 'Continue? (Y,n)'
!     read(*,*) aLetter
!     if (aLetter == 'n' .or. aLetter == 'N') then
!       call diag_print_message('Correct issue(s) and restart')
!       stop
!     endif
    endif

    if(coldstart .and. hot_out)then
      if(hot_timehr) then
        inquire(file=hotfile,exist=foundfile)
        if(foundfile)then
          open(100,file=hotfile)
          close(100,status='delete')
        endif
      endif  
      if(hot_recur) then
        inquire(file=autohotfile,exist=foundfile)
        if(foundfile)then
          open(100,file=autohotfile)
          close(100,status='delete')
        endif
      endif
    endif  

!--- Read Grid File and setup geometric variables ------------------------------
    write(*,*) 'Reading Grid'
    inquire(file=telfile,exist=foundfile)
    if(foundfile) igridtype = 1
    selectcase(igridtype)
    case(0)
      call fileext(grdfile,aext) 
      if(aext(1:2)=='h5')then
#ifdef XMDF_IO
      call read_grid_xmdf !Read non-telescoping grid from XMDF file
#else
      call diag_print_error('Cannot read XMDF grid without XMDF libraries')
#endif
      elseif(aext(1:4)=='cart')then
        call read_grid_cart
      endif
    case(1); call read_tel       !Read telescoping grid from ASCII file        
    !case(2); call read_2dm       !Read unstructured 2d mesh  
    !case(3); call read_curv     !Not implimented yet
    case default
      call diag_print_error('Unsupported Input Grid Type')
    endselect

!--- Initialize -------------------------------
    write(*,*) 'Initializing Modules'
    call geo_init                    !Geospatial variables
    call bnd_init                    !Boundary conditions
    call struct_init                 !Structures
    !call solv_init                   !Sparse matrix iterative solvers
    call flow_init                   !Hydrodynamics
    call interp_init                 !Linear Interpolation
    call der_init                    !Derivatives
    call fric_init                   !Bottom and wall friction
    if(write_sup .and. (hot_out .or. hot_recur)) call hotstart_file_init  !Initialize ASCII HotStart files - modified to check if hot-start files are to be written.
    call met_init                    !Meteorologic
    if(sedtrans) call sed_init       !Sediment transport
    if(saltrans) call sal_init       !Salinity    
    if(heattrans) call heat_init     !Temperature    
    call out_init                    !Global Output         
    if(obs_cell) call obs_cell_init  !Observation cells 
    if(save_point) call save_point_init  !Save Point cells, Mitch 5/8/2012
    call stat_init                   !Simulation Statistics
#ifdef DEV_MODE
    !if(q3d) call q3d_init            !Quasi-3D
    if(veg) call veg_init            !Vegetation
#endif    
    if(dredging) call dredge_init    !Dredge operations
#ifdef PROFILE
    call watch_init
#endif
    
!--- Print to screen and Diagnostic File --------------------------------
    call cms_print         !CMS Parameters
    call diag_print        !Diagnositics
    call geo_print         !Geospatial
    call flow_print        !Hydrodynamics
    call bnd_print         !Boundaries
    call struct_print      !Structures
    call fric_print        !Friction
    call met_print         !Wind and atmospheric pressure
    call sed_print         !Sediment transport
    call sal_print         !Salinity transport
    call heat_print        !Heat transfer
    call out_print         !Output
    call stat_print        !Simulation statistics
#ifdef DEV_MODE
    if(q3d) call q3d_print !Quasi-3D
    if(veg) call veg_print !Vegetation
#endif

    call wave_print                                   !moved out of DEV_MODE conditional call - Mitch 06/23/2017

    if(dredging) call dredge_print !dredge module   
    call rol_print         !Surface roller    
    call hot_print         !Hot start
#ifdef PROFILE
    call watch_print
#endif

    write(*,*) '**********************************************************'  !blank line after all diagnostic printing
    write(*,*) ' '

!--- Read Hot Start --------------------------------------------
    if(.not.coldstart)then    !Hot Start
      call hot_read           !Overwrites initialized variables      
      ctime = timeout*3600.0  !Added 02/13/2019 meb  !Must know ctime for waves to update 
      !call hot_init           !Initializes any missing variables such as Ctk, pbk                
    endif
    
!--- Waves ---------------------------    
    call wave_init         !Waves, Note: needs to be after hot_read

!$  if(nthr>=1) call omp_set_num_threads(nthr)

    write(*,*) 'NTHR =', nthr
    
!--- Initialize Hot Start -----------------------------------------
    if(.not.coldstart)then !Hot Start    
      call hot_init        !Initializes any missing variables such as Ctk, pbk                
    endif
    
!--- Wetting and drying and Secondary Variables --------------
    call flow_wetdry(0)
    !If all cells are dry, then assume closed domain and turn on ponding
    if(sum(1-iwet(1:ncells))==ncells) then 
      ponding=.true.
      call flow_wetdry(0)
    endif 
    iwet1 = iwet    
    eta = iwet*p*gravinv-999.0*(1-iwet)
    if(flowvolbal) then
      h0 = h
      where(h0<hmin) 
        h0 = hmin
      endwhere  
    endif
    call der_grad_eval(goa,0,zb,dzbx,dzby) !Bed-slope
    call interp_scal_cell2face(zb,0,zbk,dzbx,dzby)
    call flow_grad_interp
    !if(coldstart)then !Hot Start
    if(coldstart .or. .not.icflux) then !Hot Start
      call flow_convflux
    endif
    call fric_eval  
    call flow_eddyvis
#ifdef DEV_MODE
    if(q3d) call q3d_flow
#endif

    return
    endsubroutine prestart_EXP_shared

