!********************************************************************************    
      Subroutine prestart_EXP
!********************************************************************************    
      use EXP_Global_def,    only: fac_uw, fac_dw, advect, advcoef, mixing, time, waves, sedtransexp, watanabe, lundcirp, adeq
      use EXP_Global_def,    only: isedform, dtsed, dtmorph, dstar, wsfall, thetacr, thetac, taucr, ustcr, a0, dt
      USE EXP_transport_def, only: slpfac, xks, rhowdiv8, a0divrhowgrav, tcr, twopi, rhowdiv2
      use flow_def,  only: rhow, grav, viscos
      use cms_def,   only: noptset
      use prec_def,  only: ikind
      use const_def, only: pi,deg2Rad
      use comvarbl,  only: ctlfile,rampdur,iday,tmax,ramp,timehrs,ndsch
      use out_def,   only: save_point,obs_cell,write_capac,write_fracsusp,write_conc
      use sal_def,   only: saltrans,sal
      use stat_def,  only: flowstats,sedstats,salstats
      use sed_def,   only: rhosed,sedclass,variableD50,sedtrans,poros,dcoeff
      use hot_def,   only: coldstart
      use diag_def,  only: dgunit, dgfile
          
      implicit none
      real(ikind) s1,singleD50

      fac_UW = ADVCOEF   !Explicit Solver Parameter: upwind extrapolation
      fac_DW = 1.0_ikind

      if(ndsch .eq. 0) then
        ADVECT = .false.
      else
        ADVECT = .true.
      endif

      Mixing = .true.

      if(coldstart) then
        time = 0.0_ikind
        timeHRS = time
        tmax = Tmax
      endif
      
      if(noptset == 3) waves = .true.

      !sediment transport options
      sedtransEXP = .false.     
      watanabe = .false.
      lundcirp= .false.
      adeq = .false.
      if(isedform .ge. 1 .and. isedform .le. 3) sedtransEXP = .true.
      
      if(sedtransEXP)then
        write(*,*)'isedform = ',isedform
        write(*,*)'dtsed = ',dtsed
        write(*,*)'dtmorph = ',dtmorph   
        poros = 1.0/(1.0-poros)  !this is done only for excplicit code (implciit sets solid = 1-poros)
        if(.not. variableD50) singleD50 = sedclass(1)%diam
        SLPFAC = DCOEFF  !/tan(30*3.14159/180)  ! used in bedload slope terms
        !SLPFACC = -0.99*SLPFAC
      
        S1 = (RHOSED-RHOW)/RHOW                             !SPECIFIC GRAVITY
        DSTAR = singleD50*(S1*Grav/viscos**2.0_ikind)**(1.0_ikind/3.0_ikind)  !DIMENSIONLESS GRAIN SIZE
        if(wsfall.le.0.0) then                              !FALL SPEED (SOULSBY)
          WSFALL = viscos/singleD50*((10.36_ikind**2.0_ikind+ 1.049D0*DSTAR**3.0_ikind)**0.5_ikind-10.36_ikind)
        endif
        !CRITICAL SHIELDS PARAMETER FOR INCIPIENT MOTION (SOULSBY AND WHITEHOUSE)
        THETACR = 0.30D0/(1.0_ikind+1.2_ikind*DSTAR) +0.055D0*(1.0_ikind-EXP(-0.02_ikind*DSTAR))
        thetac = thetacr
        TAUCR = THETACR*(RHOSED-RHOW)*Grav*singleD50        !CRITICAL SHEAR STRESS
        USTCR = SQRT(TAUCR/RHOW)                           !CRITICAL SHEAR VELOCITY
      
        select case (isedform)  !choose the method and set consants
        case(1) 
          watanabe = .true.
          write_fracsusp = .false.
          write_conc = .false.
          write_capac = .false.          
          xks = 2.5*singled50
          rhowdiv8 = rhow/8.0_ikind
          A0divrhowgrav = A0/(rhow*grav)
          tcr = (rhosed-rhow)*grav*singleD50*thetac
          twopi = 2.0*pi
          rhowdiv2 = rhow/2.0_ikind
        case(2)
          lundcirp = .true.
          write_fracsusp = .false.
          write_conc = .false.
          write_capac = .false. 
        case(3)
          adeq = .true.
          write_capac = .false.          
        end select

        dtsed = max(dtsed,dt)
        dtmorph = max(dtmorph,dt)
    
        !Chris is using the elapsed times but I would prefer
        !to use time steps which are multiples of each other
        !Round dtsed to lowest multiple of dt   
        dtsed = floor(dtsed/dt)*dt           !Alex, Oct. 9, 2008
        dtmorph = floor(dtmorph/dtsed)*dtsed !Alex, Oct. 9, 2008  
      endif !sedtrans
      
      return
      end subroutine

    !*************************************************************   
    subroutine balancecheck_cards(cardname,foundcard)
    ! Reads mass/mom balance toggle and interval form advanced card file
    ! modified by Chris Reed starting with similar routines created by Alex Sanchez
    !*************************************************************
    use BalanceCheck_def, only: balcheck_int, balcheck
    
    implicit none
    integer :: ierr
    character(len=37) :: cardname, cdum
    logical :: foundcard
    
    foundcard = .true.
    select case (cardname)      
    case('BalanceCheck')
      backspace(77)
      read(77,*) cardname, BalCheck_int 
      balcheck = .true.
    case default
      foundcard = .false.
    end select
    
    return
    end subroutine balancecheck_cards