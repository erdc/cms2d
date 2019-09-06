!===========================================================================
!  CMS Bed and Wall Friction routines
! 
! Contains:
!   fric_default - Sets the default values to the friction related variables
!   fric_cards   - Reads wind data from Model Parameters file
!   fric_init    - Friction variable initialization
!   fric_print   - Prints the friction variables to the diagnostic 
!                  file and the screen
!   fric_bed     - Calculates bed friction term using one of five
!                  wave-current interaction models
!   fric_eval    - Updates the bed friction variables
!
! written by Alex Sanchez, USACE-CHL
!===========================================================================

!*********************************************************************
    subroutine fric_default()
! Sets the default values to the friction related variables
!
! written by Alex Sanchez, USACE-CHL    
!*********************************************************************
    use fric_def
    implicit none

!Bottom friction dataset
    mbedfric = 0    !1-User specified cfrict,2-Manning's equation,3-Roughness height,0-None
    constbotfric = .false. !Constant bottom friction input
    fricscale = 1.0 !Scaling factor for input roughness
    fricfile = ''
    fricpath = ''

!Wave bottom boundry layer streaming
    bbl_stream = .false.
    
!Bed slope friction factor
    fricbedslope = .false. !Include bed slope factor in bed friction

!Wave-current mean bottom friction
    mwavcurint = 1   !Wave-current interaction model
    cfricwav = 0.65  !Used in wave bottom friction term, xbeach uses 0.67, 0.65 obtained form comparing against other models

!Wall friction
    wallfric = .true.
    !nwallfric = 1   !0-None, 1-Use local Mannings, 2-Constant
    !iwallfric(1) = 1
    !iwallfric(2) = 0
    !wallman = 0.01  !Wall Manning's N coefficient
    wallfac = 1.0
    
!Dynamic roughness
    roughscalegrain = 1.0
    roughscaleripple = 1.0
    roughscalemegaripple = 1.0
    roughscaledune = 1.0
    roughscaletrans = 1.0
    mripple = 1    !Method for ripples
    mripplewave = 1    !Method for wave-generated ripples, 0-none, 1-van Rijn (1984), 2-Soulsby and Whitehouse (2005)
    mripplecurrent = 1 !Method for current-generated rippples, 0-none, 1-Soulsby (1997), 2-Soulsby and Whitehouse (2005), 3-Raudkivi (1998)
    mdune = 1          !Method for dunes (from currents), 0-none, 1-van Rijn
    mroughtranscur= 1  !Method for transport roughness
    biodegradhalflife = 50.0*3600.0 !Biodegradation half-life [sec]
    
    return
    endsubroutine fric_default

!*************************************************************   
    subroutine fric_cards(cardname,foundcard)
! Reads wind data from Model Parameters file
! written by Alex Sanchez, USACE-CHL
!*************************************************************
    use comvarbl, only: flowpath
    use geo_def, only: grdfile
    use fric_def
    use prec_def
    use diag_lib
    implicit none
    integer :: i,ierr
    character(len=37) :: cardname,cdum
    logical :: foundcard
    
    foundcard = .true.
    selectcase(cardname)    
    !=== Wall friction ======================
    case('USE_WALL_FRICTION_TERMS','USE_WALL_FRICTION')
      call card_boolean(77,wallfric,ierr)
      !Note: wallfric=.true. is a partial slip bc
      ! and wallfric=.false. is a free slip bc
      
    case('WALL_FRICTION_FACTOR')
      backspace(77)
      read(77,*) cardname, wallfac
          
   ! case('WALL_MANNINGS_N_CONSTANT')
      !backspace(77)
      !read(77,*) cardname, wallman
      !nwallfric = 2  
          
    !=== Bottom Friction =====================  
    case('BOTTOM_ROUGHNESS_METHOD')
      backspace(77)
      read(77,*) cardname,cdum
      selectcase(cdum)
      case('DYNAMIC','AUTOMATIC')
        mbedfric = 0
      case('COEFFICIENT','FRICTION_COEFFICIENT')
        mbedfric = 1
      case('MANNINGS')
        mbedfric = 2
      case('ROUGHNESS','ROUGHNESS_HEIGHT','HEIGHT')
        mbedfric = 3
      case('LINEAR_FRICTION_COEFFICIENT')
        mbedfric = 4
      case default
        call diag_print_error('Invalid bottom roughness method')
      endselect
    
    case('RIPPLE_METHOD')     !Changed from RIPPPLE_METHOD which looked like a typo.
      backspace(77)
      read(77,*) cardname,cdum
      selectcase(cdum)
      case('SOULSBY')
        mripple = 1
      case('RAUDKIVI')
        mripple = 2
      case('VAN_RIJN')
        mripple = 3  
      case('SOULSBY_WHITEHOUSE','SOULSBY-WHITEHOUSE')
        mripple = 4
      case default
        call diag_print_error('Invalid wave ripple method')
      endselect
      
    case('TRANSPORT_ROUGHNESS_METHOD')
      backspace(77)
      read(77,*) cardname,cdum
      selectcase(cdum)
      case('SOULSBY')
        mripple = 1
      case('RAUDKIVI')
        mripple = 2
      case('SOULSBY_WHITEHOUSE','SOULSBY-WHITEHOUSE')
        mripple = 2
      case default
        call diag_print_error('Invalid wave ripple method')
      endselect
         
         
    !case('WAVE_RIPPPLE_METHOD','RIPPLE_WAVE_METHOD')  
    !  backspace(77)
    !  read(77,*) cardname,cdum
    !  selectcase(cdum)
    !  case('VAN_RIJN','VANRIJN','VAN-RIJN')
    !    mripplewave = 1
    !  case('SOULSBY_WHITEHOUSE','SOULSBY-WHITEHOUSE')
    !    mripplewave = 2
    !  case default
    !    call diag_print_error('Invalid wave ripple method')
    !  endselect
    
    !case('CURRENT_RIPPPLE_METHOD','RIPPLE_CURRENT_METHOD')  
    !  backspace(77)
    !  read(77,*) cardname,cdum
    !  selectcase(cdum)
    !  case('SOULSBY')
    !    mripplecurrent = 1
    !  case('SOULSBY_WHITEHOUSE','SOULSBY-WHITEHOUSE')
    !    mripplecurrent = 2
    !  case default
    !    call diag_print_error('Invalid current ripple method')
    !  endselect
    
    case('RIPPLE_ROUGHNESS_SCALE_FACTOR','RIPPLE_ROUGHNESS_SCALING_FACTOR',&
        'RIPPLE_ROUGH_SCALE_FACTOR')  
      backspace(77)
      read(77,*) cardname,roughscaleripple
      
    case('DUNE_ROUGHNESS_SCALE_FACTOR','DUEN_ROUGHNESS_SCALING_FACTOR',&
         'DUNE_ROUGH_SCALE_FACTOR')  
      backspace(77)
      read(77,*) cardname,roughscaledune
      
    case('GRAIN_ROUGHNESS_SCALE_FACTOR','GRAIN_ROUGHNESS_SCALING_FACTOR',&
         'GRAIN_ROUGH_SCALE_FACTOR')  
      backspace(77)
      read(77,*) cardname,roughscalegrain
    
    case('TRANSPORT_ROUGHNESS_SCALE_FACTOR','TRANSPORT_ROUGHNESS_SCALING_FACTOR')  
      backspace(77)
      read(77,*) cardname,roughscaletrans
      
    case('MANNINGS_N_DATASET','MANNING_N_DATASET','MANNINGS_DATASET','MANNNING_DATASET')
      call card_dataset(77,grdfile,flowpath,fricfile,fricpath,1)
      mbedfric = 2     
      constbotfric = .false.
    
    case('MANNINGS_N_CONSTANT','MANNINGS_N_VALUE')  
      backspace(77)
      read(77,*) cardname,cbotfric    
      mbedfric = 2     
      constbotfric = .true.
          
    case('BOTTOM_FRICTION_COEF_DATASET','BOTTOM_FRICTION_DATASET','FRICTION_COEFFICIENT_DATASET')
      call card_dataset(77,grdfile,flowpath,fricfile,fricpath,1)
      mbedfric = 1 
      constbotfric = .false.
        
    case('BOTTOM_FRICTION_COEF_CONSTANT','BOTTOM_FRICTION_CONSTANT','FRICTION_COEFFICIENT_CONSTANT',&
       'BOTTOM_FRICTION_COEF_VALUE','BOTTOM_FRICTION_VALUE','FRICTION_COEFFICIENT_VALUE')
      backspace(77)
      read(77,*) cardname, cbotfric    
      mbedfric = 1 
      constbotfric = .true.          
        
    case('LINEAR_BOTTOM_FRICTION_COEF_CONSTANT','LINEAR_FRICTION_COEF_CONSTANT')
      backspace(77)
      read(77,*) cardname, cbotfric            
      mbedfric = 4  
      constbotfric = .true.         
          
    case('ROUGHNESS_HEIGHT_DATASET','ROUGHNESS_DATASET')
      call card_dataset(77,grdfile,flowpath,fricfile,fricpath,1)
      mbedfric = 3    
      constbotfric = .false.
        
    case('ROUGHNESS_HEIGHT_CONSTANT','ROUGHNESS_CONSTANT','ROUGHNESS_HEIGHT_VALUE','ROUGHNESS_VALUE')
      call card_scalar(77,'m','m',cbotfric,ierr)
      mbedfric = 3    
      constbotfric = .true.            
         
    case('BED_FRICTION_COEF')    !Wu, OLD
      backspace(77)
      read(77,*) cardname, mbedfric   !1-C_b, 2-Manning's n, 3-roughness height
          
    case('BOTTOM_FRICTION_INPUT') !Alex, OLD
      backspace(77)
      read(77,*) cardname, cdum
      do i=1,size(abedfric)
        if(cdum==abedfric(i))then
          mbedfric = i
          exit
        endif
      enddo                   
    
    !=== Wave bottom boundary layer streaming ======
    case('FLOW_BBL_STREAMING')
      call card_boolean(77,bbl_stream,ierr)    
      
    !==== Wave bottom friction ===========================        
    case('WAVE_BOTTOM_FRICTION_COEFFICIENT','WAVE_BOTTOM_FRICTION_COEF')
      backspace(77)
      read(77,*) cardname, cfricwav         
            
    case('WAVE-CURRENT_MEAN_STRESS','WAVE_CURRENT_MEAN_STRESS')
      backspace(77)
      read(77,*) cardname, cdum 
      do i=1,size(awavcur)
        if(cdum==awavcur(i))then
          mwavcurint = i
          exit
        endif
      enddo
        
    !=== Bed slope friction factor ============================
    case('BED_SLOPE_FRICTION_FACTOR','BEDSLOPE_FRICTION_FACTOR')
      call card_boolean(77,fricbedslope,ierr)       
        
    case default
      foundcard = .false.  
                          
    endselect
    
    return
    endsubroutine fric_cards

!***********************************************************************
    subroutine fric_init()
! Friction variable initialization
! written by Alex Sanchez, USACE
!***********************************************************************
#include "CMS_cpp.h"
    use size_def, only: ncells,ncellsD,nmaxfaces
    use geo_def, only: ncface,cell2cell,dzbx,dzby
    use flow_def, only: h,grav
    use comvarbl, only: nfsch
    use const_def, only: small
    use fric_def
    use fric_lib, only: fric_wavcur_init
    use sed_def, only: sedtrans,d50,d90,singlesize,singleD50
    use diag_lib
#ifdef XMDF_IO
    use in_xmdf_lib, only: readscalh5
#endif   
    use in_lib, only: readscalTxt
    implicit none
    integer :: i,k,nck,ierr
    character(len=10) :: aext
    
    !Allocate and initialize
    allocate(coefman(ncellsD),cfrict(ncellsD),z0(ncellsD))    
    allocate(uelwc(ncellsD))  
    allocate(bsxy(ncellsD),bsvel(ncellsD))  
    coefman=0.0; cfrict=0.0; z0=1.0e-5
    uelwc=0.0
    bsxy=0.0; bsvel=0.0
    if(nfsch==0)then  !For implicit scheme
      allocate(cbcfuwcap(ncellsD)) 
      cbcfuwcap=0.0
    else
      allocate(cbcfuwc(ncellsD)) 
      cbcfuwc=0.0  
    endif
    
    !Bottom roughness
    !if(len_trim(fricpath)==0 .and. .not.constbotfric)then !None specified
    !  call diag_print_warning('Invalid bottom friction specification',&
    !    '  Using a constant Manning''s N equal to 0.025')
    !  constbotfric = .true.
    !  cbotfric = 0.025
    !  mbedfric = 2
    !endif
    call fileext(trim(fricfile),aext)
    
    
    if(constbotfric)then !Constant value
      selectcase(mbedfric)
      case(1); cfrict = fricscale*cbotfric
      case(2); coefman = fricscale*cbotfric
      case(3); z0 = fricscale*cbotfric/30.0 !Convert from physical roughness height to hydraulic roughness length *********
      endselect  
    else                 !User-defined Dataset
      select case (aext)    !Added case to accept ASCII input of grid-specific datasets - meb 022118

      case('h5')
#ifdef XMDF_IO
      selectcase(mbedfric)
      case(1) !Friction coefficient
        call readscalh5(fricfile,fricpath,cfrict,ierr)
        if(ierr<0) call rough_read_error_msg(fricfile,fricpath)  
        !where(cfrict<1.0e-6) cfrict=1.0e-6
        cfrict = fricscale*cfrict
      case(2) !Manning's
        call readscalh5(fricfile,fricpath,coefman,ierr)
        if(ierr<0) call rough_read_error_msg(fricfile,fricpath)  
        !where(coefman<1.0e-6) coefman=1.0e-6
        coefman = fricscale*coefman
      case(3) !Roughness height
        call readscalh5(fricfile,fricpath,z0,ierr)
        if(ierr<0) call rough_read_error_msg(fricfile,fricpath)  
        where(z0<1.0e-5) z0=1.0e-5
        z0 = z0/30.0 !Convert from physical roughness height to hydraulic roughness length *********
        z0 = fricscale*z0
      endselect      
#endif
      case('txt') 
        selectcase(mbedfric)
        case(1) !Friction coefficient
          call readscalTxt(fricfile,cfrict,ierr)
          if(ierr<0) call rough_read_error_msg2(fricfile)  
          !where(cfrict<1.0e-6) cfrict=1.0e-6
          cfrict = fricscale*cfrict
        case(2) !Manning's
          call readscalTxt(fricfile,coefman,ierr)
          if(ierr<0) call rough_read_error_msg2(fricfile)  
          !where(coefman<1.0e-6) coefman=1.0e-6
          coefman = fricscale*coefman
        case(3) !Roughness height
          call readscalTxt(fricfile,z0,ierr)
          if(ierr<0) call rough_read_error_msg2(fricfile)  
          where(z0<1.0e-5) z0=1.0e-5
          z0 = z0/30.0 !Convert from physical roughness height to hydraulic roughness length *********
          z0 = fricscale*z0
        endselect      
      endselect  
    endif

    !Bed slope friction factor
    allocate(cmb(ncellsD))
    cmb = 1.0  !Default is 1.0. if(fricbedslope) then it is calculated based on the bed slope
    
    !call fric_eval

    !Copy bottom friction to dummy cells
    do i=1,ncells+1,ncellsD
      do k=1,ncface(i)
        nck=cell2cell(k,i)
        if(nck<=ncells)then
          cfrict(i)=cfrict(nck)
          coefman(i)=coefman(nck)
          z0(i)=z0(nck)
          exit
        endif
      enddo
    enddo 
    
    !Mean wave-current bottom friction
    if(mwavcurint>1) call fric_wavcur_init !Initialize variables
    
    !Dynamic bottom friction
    if(mbedfric==0)then
      allocate(riplen(ncellsD),riphgt(ncellsD),ripdir(ncellsD))
      riplen = 0.0; riphgt = 0.0; ripdir = 0.0
      allocate(dunelen(ncellsD),dunehgt(ncellsD))
      dunelen = 0.0; dunehgt = 0.0
      allocate(megariplen(ncellsD),megariphgt(ncellsD))
      megariplen = 0.0; megariphgt = 0.0
      allocate(rksr(ncellsD),rksmr(ncellsD),rksd(ncellsD),rksg(ncellsD),rkstc(ncellsD))
      rksr = 0.0
      rksmr = 0.0
      rksd = 0.0
      rksg = 0.0
      rkstc = 0.0
      if(.not.sedtrans)then
        allocate(d50(ncellsD),d90(ncellsD))
        if(singlesize)then
          d50 = singleD50/1000.0
        else            
          d50 = 0.2/1000.0
        endif
        d90 = d50*1.5
      endif
    endif
    
    return
    endsubroutine fric_init
    
!**************************************************    
    subroutine rough_read_error_msg(afile,apath)
!**************************************************
    use diag_lib
    implicit none
    character(len=200),intent(in) :: afile
    character(len=200),intent(in) :: apath
    character(len=200) :: msg2,msg3
    
    write(msg2,*) '       File: ',trim(afile)
    write(msg3,*) '       Path: ',trim(apath)    
    call diag_print_error('Could not find roughness dataset ',msg2,msg3,&
      '  Check input files and restart')
    
    stop         
    endsubroutine

!**************************************************    
    subroutine rough_read_error_msg2(afile)
!**************************************************
    use diag_lib
    implicit none
    character(len=200),intent(in) :: afile
    character(len=200) :: msg2,msg3
    
    write(msg2,*) '       File: ',trim(afile)
    call diag_print_error('Could not find roughness dataset ',msg2,&
      '  Check input files and restart')
    
    stop         
    endsubroutine

!**************************************************
    subroutine fric_print()
! Prints the friction related setup parameters to
! the screen and the diagnositic file
!
! written by Alex Sanchez, USACE-CHL    
!**************************************************   
    use fric_def
    use diag_def, only: dgunit,dgfile
    use cms_def,  only: noptset
    use tool_def, only: vstrlz
    implicit none
    integer :: i,iunit(2)
    character(len=200) :: apath,aname,astring
    character(len=10) :: aext
    
342 format(' ',A,T40,F0.2,A)
353 format(' ',A,T40,F0.3,A)
354  format(' ',A,T40,A,A)    !Added for vstrlz function results
345 format(' ',A,T40,F0.5)
787 format(' ',A,T40,A)       
    
    call fileparts(fricfile,apath,aname,aext)
    astring=trim(aname) // '.' // aext
    iunit = (/6, dgunit/)
    do i=1,2
      if(i==2) open(dgunit,file=dgfile,access='append')   
      write(iunit(i),*)
      write(iunit(i),787)       '  Bottom and Wall Friction'
      selectcase(mbedfric)
        case(0)
          write(iunit(i),787)   '    Dynamic bottom friction:','ON'
          if(mripplewave==1)then
            write(iunit(i),787) '    Wave-ripple method:','VAN_RIJN'
          elseif(mripplewave==2)then
            write(iunit(i),787) '    Wave-ripple method:','SOULSBY_WHITEHOUSE'
          endif
          if(mripplecurrent==1)then
            write(iunit(i),787) '    Current-ripple method:','SOULSBY'
          elseif(mripplecurrent==2)then
            write(iunit(i),787) '    Current-ripple method:','SOULSBY_WHITEHOUSE'
          endif
          write(iunit(i),787)   '    Roughness scaling factors: '
          write(iunit(i),354)   '      Ripple:',trim(vstrlz(roughscaleripple,'(F0.5)'))
          write(iunit(i),354)   '      Dune:',trim(vstrlz(roughscaledune,'(F0.5)'))
          write(iunit(i),354)   '      Grain:',trim(vstrlz(roughscalegrain,'(F0.5)'))
          write(iunit(i),354)   '      Transport:',trim(vstrlz(roughscaletrans,'(F0.5)'))
        case(1)
          if(constbotfric)then
            write(iunit(i),354) '    Constant Friction Coeff:',trim(vstrlz(cbotfric,'(F0.3)'))
          else
            write(iunit(i),787) '    Friction Coeff File:',trim(astring)
            write(iunit(i),787) '    Friction Coeff Path:',trim(fricpath)
          endif  
        case(2)
          if(constbotfric)then
            write(iunit(i),354) '    Constant Mannings N:',trim(vstrlz(cbotfric,'(F0.3)'))
          else
            write(iunit(i),787) '    Mannings N File:',trim(astring)
            write(iunit(i),787) '    Mannings N Path:',trim(fricpath)
          endif  
        case(3)  
          if(constbotfric)then
            write(iunit(i),354) '    Constant Roughness Height:',trim(vstrlz(cbotfric*1000,'(F0.3)')),' mm'          
          else
            write(iunit(i),787) '    Roughness Height File:',trim(astring)
            write(iunit(i),787) '    Roughness Height Path:',trim(fricpath)
          endif
        case(4)  
          if(constbotfric)then
            write(iunit(i),354) '    Constant Linear Friction Coeff: ',trim(vstrlz(cbotfric,'(F0.3)'))
          else
            write(iunit(i),787) 'ERROR: Linear bottom friction must be spatially constant'          
            stop
          endif    
      endselect    
      if(wallfric)then
        write(iunit(i),787)     '    Wall Friction:','ON'
        write(iunit(i),354)     '    Wall Friction Factor:',trim(vstrlz(wallfac,'(F0.3)'))
      else
        write(iunit(i),787)     '    Wall Friction:','OFF'
      endif    
      if(fricbedslope)then
        write(iunit(i),787)     '    Bed-slope Friction Factor:','ON'
      else
        write(iunit(i),787)     '    Bed-slope Friction Factor:','OFF'
      endif    
      if(noptset>=3)then
          write(iunit(i),787)   '    Wave-Current Mean Bottom Shear Stress Model:',awavcur(mwavcurint)
        if(mwavcurint==1)then
          write(iunit(i),354)   '    Wave Bottom Friction Coefficient:',trim(vstrlz(cfricwav,'(F0.3)'))
        endif
      endif 
      if(bbl_stream)then
        write(iunit(i),787)     '    Bottom boundary layer streaming:','ON'
      endif
    enddo    
    close(dgunit)            

    return        
    endsubroutine fric_print
    
!***********************************************************************
    function fric_bed(hi,Cd,z0i,ui,vi,usi,vsi,Uws,Uwr,Twr,Dw) result(cbu)
! Calculates bed friction term using one of five
! wave-current interaction models
! written by Alex Sanchez, USACE-CHL
!***********************************************************************
    use flow_def, only: rhow,viscos
    use fric_def, only: mwavcurint,cfricwav
    use fric_lib
    use const_def, only: sqrttwo
    use prec_def
    implicit none
    real(ikind),intent(in) :: hi,Cd,z0i,ui,vi,usi,vsi,Uws,Uwr,Twr,Dw
    real(ikind) :: Uc2,Uc,Uwc,taum,cbu
    
    selectcase(mwavcurint) 
    case(1) !Quad
      Uc2=(ui-usi)**2+(vi-vsi)**2 !Eulerian velocity squared  
      Uwc=sqrt(Uc2+cfricwav*Uws*Uws)      
    case(2) !Soulsby's 1995 Data2 method
      call fric_wavcurmean_data2(rhow,Cd,z0i,ui,vi,usi,vsi,Uwr,Twr,Uc,Uwc,taum)
    case(3) !Soulsby's 1995 Data13 method        
      call fric_wavcurmean_data13(rhow,Cd,z0i,ui,vi,usi,vsi,Uwr,Twr,Dw,Uc,Uwc,taum)
    case(4) !Huynh-Thanh and Temperville 1991
      call fric_wavcurmean_HT91(rhow,hi,z0i,ui,vi,usi,vsi,Uwr,Twr,Dw,Uc,Uwc,taum)
    case(5) !Fredsoe 1984
      call fric_wavcurmean_F84(rhow,Cd,z0i,ui,vi,usi,vsi,Uwr,Twr,Dw,Uc,Uwc,taum)        
    case(6)
      call fric_wavcurmean_DSK88(rhow,Cd,z0i,ui,vi,usi,vsi,Uwr,Twr,Dw,Uc,Uwc,taum)
    case default !(7)
      call fric_wavcurmean_GM79(rhow,Cd,z0i,ui,vi,usi,vsi,Uwr,Twr,Dw,Uc,Uwc,taum)
    endselect   
    cbu=Cd*Uwc
    
    return
    endfunction fric_bed
    
!***********************************************************************
    subroutine fric_eval()
! Calculates current mangitud and bed shear variables
! Updates coefman,z0,cfrict,uv,uelwc,bsxy,bsvel
! last modified April 2, 2010
! by Weiming Wu, Jan. 2010
!***********************************************************************
#include "CMS_cpp.h"
    use size_def
    use geo_def, only: dzbx,dzby,areap !!,zb
    use flow_def, only: h,u,v,us,vs,uv,rhow,viscos,grav !!,p,hmin,iwet
    use fric_def
    use fric_lib
    use sed_def, only: rhosed,calcmorph,d50,d90
    use wave_flowgrid_def, only: Worb,Worbrep,Wper,Wang,Wunitx,Wunity
    use const_def, only: small,sqrttwo !!,gravinv
    use cms_def, only: noptset
    use comvarbl, only: dtime
    implicit none
    integer :: i
    real(ikind), parameter :: hmax = 10.0 ![m]
    !real(ikind) :: hold,fac
    logical :: isnankind

    selectcase(mbedfric) ! user specified  C_b
    case(1)
!$OMP PARALLEL DO PRIVATE(i)            
      do i=1,ncells
        coefman(i) = fric_conv_drag2man(grav,cfrict(i),h(i),hmax)
        z0(i) = fric_conv_drag2length(cfrict(i),h(i))
      enddo
!$OMP END PARALLEL DO

    case(2) !Manning's equation  
!$OMP PARALLEL DO PRIVATE(i)       
      do i=1,ncells
        cfrict(i) = fric_conv_man2drag(grav,coefman(i),h(i),hmax)          
        z0(i) = fric_conv_drag2length(cfrict(i),h(i)) 
#ifdef DIAG_MODE
        if(cfrict(i)>0.5 .or. isnankind(cfrict(i)))then
          continue
        endif
#endif        
      enddo
!$OMP END PARALLEL DO
       
    case(3) !Roughness height
!$OMP PARALLEL DO PRIVATE(i)             
      do i=1,ncells
        cfrict(i) = fric_conv_length2drag(z0(i),h(i))
        coefman(i) = fric_conv_drag2man(grav,cfrict(i),h(i),hmax)
      enddo
!$OMP END PARALLEL DO

    case(4) !Linear bottom friction
!$OMP PARALLEL DO PRIVATE(i)             
      do i=1,ncells
        cfrict(i) = cbotfric*h(i)/max(uelwc(i),1e-15)
        coefman(i) = fric_conv_drag2man(grav,cfrict(i),h(i),h(i))
        z0(i) = fric_conv_drag2length(cfrict(i),h(i)) 
      enddo
!$OMP END PARALLEL DO

    case(0) !Dynamic bottom friction
!$OMP PARALLEL DO PRIVATE(i) 
      do i=1,ncells
        cfrict(i) = fric_conv_length2drag(z0(i),h(i))
        coefman(i) = fric_conv_drag2man(grav,cfrict(i),h(i),hmax)
      enddo
!$OMP END PARALLEL DO

    endselect
      
    !Update bottom slope friction coefficient
    if(fricbedslope .and. calcmorph)then !Only recalculate if morphology change is on
!$OMP PARALLEL DO PRIVATE(i)    
        do i=1,ncells
        cmb(i)=sqrt(1.0+dzbx(i)*dzbx(i)+dzby(i)*dzby(i))
      enddo
!$OMP END PARALLEL DO
    endif

!!#ifdef DEV_MODE
!!    !Increase friction to avoid negative water depths
!!    do i=1,ncells
!!      if(iwet(i)==1)then
!!        hold=p(i)*gravinv-zb(i)
!!        if(hold<hmin)then
!!          fac = 1.0+exp(0.005*(hmin-hold)/hmin)
!!          cfrict(i) = cfrict(i)*fac
!!        endif
!!      endif
!!    enddo
!!#endif
    
!Note: The methods below except 1 were designed for regular waves but
! but are being implimented here for random waves. 
! The bottom orbital velocity is specified here based on the significant
! wave height and not the root-mean-squared wave like like Delft-3D
    if(noptset>=3)then
      selectcase(mwavcurint) 
      case(1) !Wu 2009
!$OMP PARALLEL DO PRIVATE(i)                 
        do i=1,ncells
          call fric_wavcurmean_quad(rhow,cfrict(i),u(i),v(i),us(i),vs(i),&  
               Worb(i),uv(i),uelwc(i),bsxy(i),bsvel(i))
        enddo         
!$OMP END PARALLEL DO
             
      case(2) !Soulsby's 1995 Data2 method
!$OMP PARALLEL DO PRIVATE(i)
        do i=1,ncells
          call fric_wavcurmean_data2(rhow,cfrict(i),z0(i),u(i),v(i),&
               us(i),vs(i),Worbrep(i),Wper(i),uv(i),uelwc(i),bsxy(i))
          bsvel(i)=sqrt(bsxy(i)/rhow+small)   
        enddo      
!$OMP END PARALLEL DO
                                  
      case(3) !Soulsby's 1995 Data13 method
!$OMP PARALLEL DO PRIVATE(i)              
        do i=1,ncells
          call fric_wavcurmean_data13(rhow,cfrict(i),z0(i),u(i),v(i),&
               us(i),vs(i),Worbrep(i),Wper(i),Wang(i),uv(i),uelwc(i),bsxy(i))
          bsvel(i)=sqrt(bsxy(i)/rhow+small)   
        enddo
!$OMP END PARALLEL DO
             
      case(4) !Huynh-Thanh and Temperville 1991
!$OMP PARALLEL DO PRIVATE(i)              
        do i=1,ncells
          call fric_wavcurmean_HT91(rhow,h(i),z0(i),u(i),v(i),us(i),vs(i),&
               Worbrep(i),Wper(i),Wang(i),uv(i),uelwc(i),bsxy(i))
          bsvel(i)=sqrt(bsxy(i)/rhow+small)   
        enddo
!$OMP END PARALLEL DO    
       
      case(5) !Fredsoe 1984
!$OMP PARALLEL DO PRIVATE(i)               
        do i=1,ncells
          call fric_wavcurmean_F84(rhow,cfrict(i),z0(i),u(i),v(i),us(i),vs(i),&
               Worbrep(i),Wper(i),Wang(i),uv(i),uelwc(i),bsxy(i))
          bsvel(i)=sqrt(bsxy(i)/rhow+small)   
        enddo
!$OMP END PARALLEL DO        

      case(6) !Davies et al. 1988
!$OMP PARALLEL DO PRIVATE(i)               
        do i=1,ncells
          call fric_wavcurmean_DSK88(rhow,cfrict(i),z0(i),u(i),v(i),us(i),vs(i),&
               Worbrep(i),Wper(i),Wang(i),uv(i),uelwc(i),bsxy(i))
          bsvel(i)=sqrt(bsxy(i)/rhow+small) 
        enddo
!$OMP END PARALLEL DO

      case(7) !Grant and Madsen 1979
!$OMP PARALLEL DO PRIVATE(i)               
        do i=1,ncells
          call fric_wavcurmean_GM79(rhow,cfrict(i),z0(i),u(i),v(i),us(i),vs(i),&
               Worbrep(i),Wper(i),Wang(i),uv(i),uelwc(i),bsxy(i))
          bsvel(i)=sqrt(bsxy(i)/rhow+small)   
        enddo
!$OMP END PARALLEL DO   

      endselect
         
    else    !No waves   
!$OMP PARALLEL DO PRIVATE(i)    
      do i=1,ncells
        call fric_curmean_quad(rhow,cfrict(i),u(i),v(i),uv(i),uelwc(i),bsxy(i),bsvel(i))                              
      enddo
!$OMP END PARALLEL DO
    endif
    
    !!cfrict=0.0
    !!uelwc=0.0
    !!bsxy=0.0
    !!bsvel=0.0
    
    return
    endsubroutine fric_eval
    
!*******************************************************************************
    subroutine fric_rough_eval
!*******************************************************************************
    use size_def
    use geo_def, only: dzbx,dzby,areap !!,zb
    use flow_def, only: h,u,v,us,vs,uv,rhow,viscos,grav !!,p,hmin,iwet
    use fric_def
    use fric_lib
    use sed_def, only: rhosed,calcmorph,d50,d90
    use wave_flowgrid_def, only: Worb,Worbrep,Wper,Wang,Wunitx,Wunity
    use const_def, only: small,sqrttwo !!,gravinv
    use cms_def, only: noptset
    use comvarbl, only: dtime
    implicit none
    integer :: i
    real(ikind) :: rksc,Uw,Tw,Dw,dx
    
    Uw = 0.0
    Tw = 8.0
    Dw = 0.0
    do i=1,ncells
      if(noptset>=3)then
        Uw = worbrep(i)
        Tw = wper(i)
        Dw = wang(i)
      endif
      !rhos,d50,d90,rhow,vsk,h,Uc,U,V,Uw,Tw,Wx,Wy
      dx = sqrt(areap(i))
      call fric_roughness(rhosed,d50(i),d90(i),&    !Sediment
        rhow,viscos,h(i),uv(i),u(i),v(i),&  !Hydro
        Uw,Tw,Dw,&                         !Waves
        riphgt(i),riplen(i),ripdir(i),megariphgt(i),megariplen(i),dunehgt(i),dunelen(i),& !Bedforms
        rksr(i),rksmr(i),rksd(i),rksg(i),rkstc(i),&  !Roughnesses
        biodegradhalflife,dx,dtime)  !Biodegradation, spatial resolution, and time step
      rksc = rksr(i) + rksmr(i) + rksd(i) + rksg(i) + rkstc(i) !Total roughhness height [m]
      z0(i) = rksc/30.0
    enddo
    
    return
    endsubroutine fric_rough_eval
        
!*******************************************************************************
    subroutine fric_roughness( &
        rhos,d50,d90,rhow,&  !Sediment
        vsk,h,Uc,U,V,&  !Hydro
        Uw,Tw,Dw,&   !Waves
        rhgt,rlen,rdir,mrhgt,mrlen,dhgt,dlen, & !Bedforms
        rksr,rksmr,rksd,rksg,rkstc,& !Roughnesses
        Tbhalf,dx,dt)
! Computes the total bed roughness height (felt by currents)
!
! Description:
!   The roughness felt by currents includes wave and current-generated
!   rippples, dunes, skin friction, and current-related sediment transport
!
! Author: Alex Sanchez, USACE-CHL
!*******************************************************************************
    use prec_def
    use const_def, only: pi
    use flow_def, only: grav
    use fric_def, only: roughscaleripple,&
       roughscaledune,roughscaletrans,roughscalegrain,&
       mripple,mdune,mroughtranscur,roughscalemegaripple
    use fric_lib
    use sed_def, only: d2dstar
    use sed_lib, only: shields_soulsby
    implicit none
    !Input
    real(ikind),intent(in) :: Tbhalf,dx,dt
    real(ikind),intent(in) :: rhos,d50,d90,rhow,vsk,h,Uc,U,V,Uw,Tw,Dw
    !Input/Output
    real(ikind),intent(inout) :: rhgt,rlen,rdir,mrhgt,mrlen,dhgt,dlen
    !Output
    real(ikind),intent(out) :: rksr,rksmr,rksd,rksg,rkstc
    !Internal
    real(ikind) :: taup,Cdg,thetap,thetacr
    real(ikind) :: s,Ts,z0g,dstar50,r,fac
    
    s = rhos/rhow
    
    !Grain roughness
    rksg = roughscalegrain*fric_rough_grain(d90) ![m]
    
    !Skin shear stress and shields parameter
    z0g = rksg/30.0
    Cdg = fric_conv_length2drag(z0g,h)
    taup = rhow*Cdg*Uc*Uc
    thetap = taup/((rhos-rhow)*grav*d50)
    
    !Critical shields for d50
    dstar50 = d50*((s-1.0)*grav)**0.3333333/vsk**0.66666667
    thetacr = shields_soulsby(dstar50)
    
    !Dune roughness
    Ts = thetap/thetacr - 1.0
    call fric_dune_vanrijn(h,d50,Ts,dhgt,dlen)
    rksd = roughscaledune*fric_rough_dune(dhgt,dlen)
    !Make sure dunes are a subgrid feature
    r = dlen/dx
    fac = 0.5+0.5*cos(pi*min(r,1.0))
    rksd = rksd*fac
    
    !Megaripples
    call fric_megaripple(h,Uc,Uw,mrhgt,mrlen)
    !mrhgt = min(0.02*h,0.1)
    !mrlen = min(1*h,5.0)
    rksmr = roughscalemegaripple*fric_rough_ripple(mrhgt,mrlen)
    
    !Notes: 
    ! 1. Should be separated into current-only and current and waves cases
    ! 2. Should keep track of ripple orientation as well
    
    !Current ripple roughness
    if(mripple==1)then !Soulsby (1997)
      call fric_ripple_current_soulsby(d50,rlen,rhgt)
    elseif(mripple==2)then !Raudkivi (1998)
      call fric_ripple_current_raudkivi(d50,rhgt,rlen)
    elseif(mripple==3)then !Van Rijn (2007)
      !call fric_ripple_wave_vanrijn(s,d50,Uw,Tw,riphgtw,riplenw) !Includes both wave and current related ripples  
    elseif(mripple==4)then !Soulsby and Whitehouse (2005)  
      call fric_ripple_soulsbywhitehouse(grav,s,d50,dstar50,thetacr,&
         h,Uc,Uw,Tw,Dw,Tbhalf,rhgt,rlen,rdir,dt)
    !else
    endif  
    rksr = roughscaleripple*fric_rough_ripple(rhgt,rlen) ![m]
    
    !Transport roughness (related to currents)
    if(mroughtranscur==1)then
      rkstc = roughscaletrans*fric_rough_trans_nielson(d50,thetap,thetacr) ![m]
    elseif(mroughtranscur==2)then
      rkstc = roughscaletrans*fric_rough_trans_wiberg(d50,Ts) ![m]
    endif
    
    !Total roughness
    !rksc = rksrc + rksrwc + rksd + rksg + rkstc !Total roughhness height [m]
    !rksc = sqrt(rksrc**2 + rksrw**2 + rksd**2 + rksg**2 + rkstc**2) !Total roughhness height [m]
    !rksc = max(rksc,1.0e-4)
        
    return
    endsubroutine fric_roughness
    