!============================================================================
! Implicit Flow Model routines
!
! Contains:
!   flow_imp - Implicit flow solver main routine
!   flow_srcsnk_static - Static components of the hydro source and 
!        sink terms are calculated once and stored for efficiency 
!   coeffsourcesink_uv - Calculates amtrix coefficients, sources, 
!       and sinks for u and v
!   coefsourcesink_pp - Calculates the matrix coefficients, sources, 
!       and sinks for the pressure-correction equation 
!   pvfc - Correct velocity, mass flux and pressure
!   flow_grad_interp - Calculates the flow gradients and interpolated values
!   flow_recalculate - Reinitializes time step variables and 
!       reduces the time step
!   velface_eta - Updates the cell center u and v velocities using face fluxes  
!       and calcualtes the water surface elevation from the pressure
!   check_conv - Checks the implicit hydrodynamic model solver convergence and 
!    outputs the variable with the convergence state
!   flow_pred - Predictor step for implicit temporal solution scheme
!   check_momentum - Checks the momentum equation matrix and R.H.S.
!   check_variables - Checks the matrix and R.H.S.
!
! written by Weiming Wu, NCCHE
!   and Alejandro Sanchez, USACE-CHL
!============================================================================

!***********************************************************************
    recursive subroutine flow_imp
! Implicit Flow solver
! written by Weiming Wu, Clarkson University
!            Alejandro Sanchez, USACE-CHL
!***********************************************************************
#include "CMS_cpp.h"
    use flow_def
    use comvarbl
    use cms_def,  only: cmswave
    use comp_lib, only: zerocoef, hybridcoef, powerlawcoef, exponentialcoef, upwindcoef, gammadefcor, cubistadefcor, alvsmartdefcor, hoabdefcor
    use diag_lib, only: diag_print_error, diag_print_message
    use diag_def, only: msg, msg2, debug_mode, dgunit, dgfile
    use size_def, only: ncells, ncellsimple
    use flow_lib, only: number_wet_cells
    use met_def,  only: presconst, windconst, windvar, presvar, windsta, pressta, iwndlagr, rain_evap
    use solv_def, only: iconv
    
#ifdef DEV_MODE
    use q3d_def
#endif

#ifdef PROFILE
    use watch_lib
#endif

    implicit none
    integer :: i

!3   format('ntime=',i8,',  dt=',F8.3,', time=',1pe12.5)
3   format('ntime=',i8,',  dt=',F8.3,', time=',1pe12.5,', Active Cells=',i0)
4   format(1x,i8,3(3x,1pe12.4))
5   format('Forcing ramp % = ',F6.2)
    
#ifdef PROFILE
    call watch_start('flow_imp prep')
#endif

    !Hydro Boundary conditions
    call bnd_eval
    
    !Flow Wetting and drying
#ifdef PROFILE
    call watch_start('flow_wetdry')
#endif    
    call flow_wetdry(0)
#ifdef PROFILE
    call watch_stop('flow_wetdry')
#endif 
    
    !Update derivative coefficients
    call der_update

    !Wave Wetting and Drying
    if(cmswave) call wave_wetdry
     
    !Wind and Atmospheric Pressure
    if(windconst) call windcurve_eval    !Spatially constant wind field
    if(presconst) call prescurve_eval    !Spatially constant atmospheric pressure    
    if(windvar .or. presvar) call windpresfield_eval !Spatially variable wind and pressure fields   
    if(windsta .or. pressta) call metsta_eval   !Spatially variable wind and pressure stations
   
    !Print to screen
!    write(msg,3) ntime,dtime,ctime

    write(msg,3) ntime,dtime,ctime,number_wet_cells()
    if (ramp.lt.1) then
      write(msg2,5) ramp*100    
      call diag_print_message(' ',msg2,msg,' Flow: iter    p_res          U_res          V_res')
    else
      call diag_print_message(' ',msg,' Flow: iter    p_res          U_res          V_res')
    endif
    
#ifdef PROFILE
    call watch_start('flow_srcsnk_static')
#endif    
    call flow_srcsnk_static !ssu0,ssv0,sspp0,spuv0,sppp0
#ifdef PROFILE
    call watch_stop('flow_srcsnk_static')
    call watch_stop('flow_imp prep')
    call watch_start('flow_imp outer')
#endif

    do niter=1,maxit !outer iteration loop 
      call bndriver_vel       !Adjust velocity at river boundaries
      if(iwndlagr==1) call wind_lagrangian_update !Lagrangian wind reference frame

!--- U and V momentum Equations --------------------------------------
#ifdef PROFILE
      call watch_start('flow_imp outer uv')
      call watch_start('flow_imp outer uv assem')
#endif
      select case(ndsch) !coefficients, sources, and sinks
      case(0); call coeffsourcesink_uv(zerocoef)
      case(2); call coeffsourcesink_uv(hybridcoef)              !Chris' code has 1
      case(3); call coeffsourcesink_uv(powerlawcoef)          !Chris' code has 2
      case(4); call coeffsourcesink_uv(exponentialcoef)          !Chris' code has 3
      case default; call coeffsourcesink_uv(upwindcoef)
      end select
      !Deferred anti-diffusion corrections
      if(ncellsimple<ncells)then !Not a simple Cartesian grid
        select case(ndsch) 
        case(0); !Do nothing. No deferred correction for this case    
        case(5); call defcorhlpagradvec(u,v,dux,duy,dvx,dvy,su,sv)                    !Chris' code has 4
        !case(5); call defcorhlpagradvecnew(u,v,dux,duy,dvx,dvy,su,sv)
        case(6); call defcorgammagradvec(gammadefcor,u,v,dux,duy,dvx,dvy,su,sv)        !Chris' code has 5
        case(7); call defcorgammagradvec(cubistadefcor,u,v,dux,duy,dvx,dvy,su,sv)      !Chris' code has 6
        case(8); call defcorgammagradvec(alvsmartdefcor,u,v,dux,duy,dvx,dvy,su,sv)  !Chris' code has 7 
        case(9); call defcorgammagradvec(hoabdefcor,u,v,dux,duy,dvx,dvy,su,sv)        !Chris' code has 8
        case default; call defcorparagradvec(dux,duy,dvx,dvy,su,sv) !Only deferred skewness corrections
        end select
      else !Simple Cartesian grid
        select case(ndsch)
        case(0); !Do nothing. No deferred correction for this case    
        case(5); call defcorhlpavec(u,v,su,sv)
        case(6); call defcorgammavec(gammadefcor,u,v,su,sv)
        case(7); call defcorgammavec(cubistadefcor,u,v,su,sv)
        case(8); call defcorgammavec(alvsmartdefcor,u,v,su,sv)
        case(9); call defcorgammavec(hoabdefcor,u,v,su,sv)
        end select
      endif
      !if(debug_mode) call check_momentum !Check Variables
      call bound_uv !Boundary Conditions
      if(debug_mode) call check_momentum !Check Variables
#ifdef PROFILE
      call watch_stop('flow_imp outer uv assem')
      call watch_start('flow_imp outer uv solve')
#endif
      call solve(acoef,su,spu,rsu,u,2) !Solve for U
      if(iconv==0)then
        if(dtvar)then
          call flow_recalculate
          call flow_imp
          return
        else
          if(mod(niter,nprint)/=0)then
            write(msg,4) niter,(rmom(i),i=1,3)  
            call diag_print_message(msg)
          endif  
          call flow_step_stat(1)  
          call write_debug  
          call diag_print_error('Model is divergent',&
            '  Please check model setup and restart run')  
        endif
      endif
      call solve(acoef,sv,spv,rsv,v,3) !Solve for V
      if(iconv==0)then        
        if(dtvar)then
          call flow_recalculate
          call flow_imp
          return
        else
          if(mod(niter,nprint)/=0)then
            write(msg,4) niter,(rmom(i),i=1,3)    
            call diag_print_message(msg)
          endif
          call flow_step_stat(1)  
          call write_debug  
          call diag_print_error('Model is divergent',&
            '  Please check model setup and restart run')
        endif
      endif

#ifdef PROFILE
      call watch_stop('flow_imp outer uv')
      call watch_stop('flow_imp outer uv solve')
      call watch_start('flow_imp outer pp')
      call watch_start('flow_imp outer pp assem')
#endif
      
!--- Pressure-correction equation ------------------------
      call coefsourcesink_pp
      call bound_pp
      if(debug_mode) call check_variables(1)
#ifdef PROFILE
      call watch_stop('flow_imp outer pp assem')
      call watch_start('flow_imp outer pp solve')
#endif
      call solve(acoef,su,sp,rsp,pp,1)
      if(iconv==0)then
        if(dtvar)then
          call flow_recalculate
          call flow_imp
          return
        else
          if(mod(niter,nprint)/=0)then
            write(msg,4) niter,(rmom(i),i=1,3)    
            call diag_print_message(msg)
          endif 
          call flow_step_stat(1)  
          call write_debug  
          call diag_print_error('Model is divergent',&
            '  Please check model setup and restart run')
        endif  
      endif
#ifdef PROFILE
      call watch_stop('flow_imp outer pp solve')
#endif      
      call pvfc
#ifdef PROFILE
      call watch_stop('flow_imp outer pp')
      call watch_start('flow_imp outer var')
      call watch_start('flow_grad_interp')
#endif
      call flow_grad_interp !dpx,dpy,dux,duy,hk,pk
#ifdef PROFILE
      call watch_stop('flow_grad_interp')
#endif
      call bndriver_vel       !Adjust velocity at river boundaries
      call bndflux
#ifdef PROFILE
      call watch_start('fric_eval')
#endif   
      call fric_eval
#ifdef PROFILE
      call watch_stop('fric_eval')
      call watch_start('flow_eddyvis')
#endif
      call flow_eddyvis
#ifdef PROFILE
      call watch_stop('flow_eddyvis')
#endif       
#ifdef DEV_MODE
      if(q3d) call q3d_flow
#endif      
      if(debug_mode) call flow_step_stat(0)
#ifdef PROFILE
      call watch_stop('flow_imp outer var')
#endif
      
      if(niter==2) write(6,4) niter-1,(rmom(i),i=1,3)
      if(mod(niter,nprint)==0) write(6,4) niter,(rmom(i),i=1,3)

      call check_conv
      if(iconv==0)then !Divergent
        if(mod(niter,nprint)/=0) write(6,4) niter,(rmom(i),i=1,3)
        open(dgunit,file=dgfile,access='append') 
        write(dgunit,4) niter,(rmom(i),i=1,3)
        close(dgunit)
        call flow_step_stat(1)
        call write_debug
        if(dtvar)then
          call flow_recalculate
          call flow_imp
          return
        else
          call diag_print_error('Model is divergent',&
            '  Please check model setup and restart run')  
        endif
      elseif(iconv==2)then !Converged
        if(mod(niter,nprint)/=0 .and. .not.debug_mode)then
          write(6,4) niter,(rmom(i),i=1,3)
        endif
        open(dgunit,file=dgfile,access='append') 
        write(dgunit,4) niter,(rmom(i),i=1,3)
        close(dgunit)
        exit  
      elseif(niter==maxit)then
        if(mod(niter,nprint)/=0) write(6,4) niter,(rmom(i),i=1,3)
        open(dgunit,file=dgfile,access='append') 
        write(dgunit,4) niter,(rmom(i),i=1,3)
        close(dgunit)
      endif
      !!if(niter==1) call write_debug
    enddo !outer iteration loop

#ifdef PROFILE
    call watch_stop('flow_imp outer')
    call watch_start('flow_imp adjust')
#endif

    if(rain_evap) call flow_rainevap !******************************
    
    call velface_eta       !Determine u and v using face fluxes and water surface elevation    
    call fric_eval         !Finalize current magnitude and bed shear stresses
    call fluxthrustructure !Flow through structures    
    call flow_step_stat(1) !update statistics and print to screen
    if(flowvolbal) call flow_balance

#ifdef PROFILE
    call watch_stop('flow_imp adjust')
#endif   
    
    return
    end subroutine flow_imp
    

!*******************************************************
    subroutine flow_srcsnk_static
! Static components of the hydro source and sink terms are
! calculated once and stored for efficiency
!*******************************************************    
    use geo_def
    use met_def
    use flow_def
    use size_def,  only: ncells
    use cms_def,   only: noptset
    use comvarbl,  only: dtime,ntsch,ctsch1,ctsch2,ramp
    use const_def, only: pi
    use prec_def,  only: ikind
    use wave_flowgrid_def, only: wavestrx,wavestry

    implicit none
    integer :: i
    real(ikind):: fac1,fac2,val,gravdtimeinv,dtimeinv
    
    gravdtimeinv = 1.0/(grav*dtime)
    dtimeinv = 1.0/dtime
!$OMP PARALLEL    
!--- Static Sources and Sinks ---------------
    if(ntsch==1)then
!$OMP DO PRIVATE(i,fac1)
      do i=1,ncells
        fac1=areap(i)*h1(i)*dtimeinv
        spuv0(i)=-fac1  !Same for u and v
        ssu0(i)=fac1*u1(i)
        ssv0(i)=fac1*v1(i)        
        sppp0(i)=-areap(i)*gravdtimeinv !must be non-positive
        sspp0(i)=-p1(i)*sppp0(i)
      enddo   
!$OMP END DO
    else
!$OMP DO PRIVATE(i,val,fac1,fac2)    
      do i=1,ncells
        val=areap(i)*dtimeinv
        fac1=val*ctsch1*h1(i)
        fac2=val*ctsch2*h2(i)
        spuv0(i)=-fac1+fac2  !Same for u and v
        ssu0(i)=fac1*u1(i)-fac2*u2(i)
        ssv0(i)=fac1*v1(i)-fac2*v2(i)
        sppp0(i)=-areap(i)*gravdtimeinv !must be non-positive
        sspp0(i)=(-ctsch1*p1(i)+ctsch2*p2(i))*sppp0(i)
      enddo
!$OMP END DO
    endif
    
!Note: Preciptation and Evaporation moved to flow_rainevap
!    !Precipitation and Evaporation
!    if(rain_evap)then
!      val = ramp*(rain-evap)/3600.0 !Note conversion from m/hr to m/s
!!$OMP DO PRIVATE(i)
!      do i=1,ncells
!        sspp0(i)=sspp0(i)+val*areap(i)*iwet(i) !source/sink
!      enddo
!!$OMP END DO        
!    endif
    
    !Wave forcing
    if(noptset>=3)then
!$OMP DO PRIVATE(i)       
      do i=1,ncells
        ssu0(i)=ssu0(i)+wavestrx(i)*areap(i)
        ssv0(i)=ssv0(i)+wavestry(i)*areap(i)
      enddo    
!$OMP END DO 
    endif

    !Wind (Eulerian Reference Frame)
    if(iwndlagr==0)then
      if(windconst)then
!$OMP DO PRIVATE(i)    
        do i=1,ncells    
          ssu0(i)=ssu0(i)+tauwx*areap(i)
          ssv0(i)=ssv0(i)+tauwy*areap(i)
        enddo
!$OMP END DO  
      elseif(windvar .or. windsta)then    
!$OMP DO PRIVATE(i)    
        do i=1,ncells    
          ssu0(i)=ssu0(i)+tauwindx(i)*areap(i)
          ssv0(i)=ssv0(i)+tauwindy(i)*areap(i)
        enddo
!$OMP END DO        
      endif
    endif
!$OMP END PARALLEL

    return        
    end subroutine flow_srcsnk_static
    
!***********************************************************************
    subroutine coeffsourcesink_uv(schmcoef)
! Calculates matrix coefficients, sources, and sinks for u and v
! written by Alex Sanchez, USACE-CHL
!***********************************************************************
#include "CMS_cpp.h"
    use comvarbl
    use geo_def
    use prec_def, only: ikind
    use size_def, only: ncells
    use cms_def,  only: noptset
    use der_def,  only: gow
    use der_lib,  only: dx2d,dy2d
    use flow_def, only: flux,u,v,us,vs,h,hk,visk,fc,waveflux,iwet,dpx,dpy,acoef,su,sv,sp,ssu0,ssv0,spuv0,apuareap,sumu,uv,rhow
    use fric_def, only: cbcfuwcap,cfrict,uelwc,cmb,z0,bbl_stream
    use fric_lib, only: fric_normapprough,fric_streaming_stress
    use met_def,  only: iwndlagr,windconst,windvar,presvar,windsta,cdWndareap,wndx,wndy,uwind,vwind,pressatmdx,pressatmdy,pressta
    use wave_flowgrid_def, only: wunitx,wunity,worbrep,wlen,wper
#ifdef DEV_MODE
    use q3d_def
    use veg_def
#endif

    implicit none
    integer :: i,k
    real(ikind) :: ddk,densitinv,za,asum,apu,val,dvarx,dvary
    
    interface
      function schmcoef(dk,fk)
        use prec_def
        implicit none
        real(ikind),intent(in) :: dk,fk
        real(ikind) :: schmcoef
      end function
    endinterface
    
!$OMP PARALLEL
!$OMP DO PRIVATE(i,k,ddk)
!!!$OMP DO PRIVATE(i,k,ddk,apu,asum)
    do i=1,ncells
      su(i)=ssu0(i)+( fc(i)*v(i)-dpx(i))*h(i)*areap(i)
      sv(i)=ssv0(i)+(-fc(i)*u(i)-dpy(i))*h(i)*areap(i)
      !!su(i)=ssu0(i)-( fc(i)*v(i)-dpx(i))*zb(i)*areap(i) !Linearized, for analytical cases
      !!sv(i)=ssv0(i)-(-fc(i)*u(i)-dpy(i))*zb(i)*areap(i) !Linearized, for analytical cases
      cbcfuwcap(i)=cmb(i)*cfrict(i)*uelwc(i)*areap(i)
      sp(i)=spuv0(i)-cbcfuwcap(i)
      !Note: by setting visk and flux to zero at dry faces, acoef is also equal to zero
      do k=1,ncface(i)
        ddk=visk(k,i)*hk(k,i)*dsxy(k,i)
        acoef(k,i)=schmcoef(ddk,flux(k,i))
      enddo
      !!!Save variables for momentum interpolation and pressure-correction equation
      !!!Notes: These are the same for u and v. Wind sink term is ignored for the sake of speed
      !!asum=sum(acoef(1:ncface(i),i))
      !!apu=iwet(i)*relax(2)/max((asum-sp(i)),1.0e-6) !=1/ap(i)
      !!apuareap(i)=apu*areap(i) !Pressure-correction and momentum interpolation
      !!sumu(i)=asum*apu        !Pressure-correction equation
    enddo 
!$OMP END DO 
    
    !Wind (Lagrangian Reference Frame)
    if(iwndlagr==1)then
      if(windconst)then
!$OMP DO PRIVATE(i)  
        do i=1,ncells
          su(i)=su(i)+cdWndareap(i)*(wndx+us(i)) !source
          sv(i)=sv(i)+cdWndareap(i)*(wndy+vs(i)) !source
          sp(i)=sp(i)-cdWndareap(i)      !sink, must be non-positive
        enddo
!$OMP END DO
      elseif(windvar .or. windsta)then
!$OMP DO PRIVATE(i)  
        do i=1,ncells
          su(i)=su(i)+cdWndareap(i)*(uwind(i)+us(i)) !source
          sv(i)=sv(i)+cdWndareap(i)*(vwind(i)+us(i)) !source
          sp(i)=sp(i)-cdWndareap(i)          !sink, must be non-positive
        enddo
!$OMP END DO      
      endif
    endif   
    
    !Wave velocity forcing
    if(noptset>=3 .and. waveflux)then 
!$OMP DO PRIVATE(i)  
      do i=1,ncells       
        su(i)=su(i)+cbcfuwcap(i)*us(i)
        sv(i)=sv(i)+cbcfuwcap(i)*vs(i)
      enddo
!$OMP END DO
    endif
            
    !Atmospheric pressure gradients
    if(presvar .or. pressta)then
      densitinv=1.0/rhow
!$OMP DO PRIVATE(i)  
      do i=1,ncells       
        su(i)=su(i)-pressatmdx(i)*h(i)*areap(i)*densitinv
        sv(i)=sv(i)-pressatmdy(i)*h(i)*areap(i)*densitinv
      enddo
!$OMP END DO  
    endif
    
#ifdef DEV_MODE
    !Bottom streaming
    if(noptset>=3 .and. waveflux .and. bbl_stream)then 
!$OMP DO PRIVATE(i)  
      do i=1,ncells
        !za = h(i)*fric_normapprough(uelwc(i),uv(i),cfrict(i)) !Apparent roughness
        !za = 0.00533*max(worbrep(i),1.0e-5)**2.25  !Apparent roughness
        !za = 1.6667e-05 !=0.2/1000/12
        za = z0(i)
        val = fric_streaming_stress(rhow,za,worbrep(i),wper(i),wlen(i))*areap(i)/rhow
        su(i)=su(i)+val*wunitx(i)
        sv(i)=sv(i)+val*wunity(i)
      enddo
!$OMP END DO
    endif

    !3D Dispersion Terms
    if(q3d .and. q3d_to_flow)then
!$OMP DO PRIVATE(i,dvarx,dvary)
      do i=1,ncells
        call dx2d(gow,i,f3dxx,dvarx)
        call dy2d(gow,i,f3dxy,dvary)
        f3du(i) = -(dvarx+dvary)
        su(i)=su(i)-(dvarx+dvary)*areap(i)*ramp
        call dx2d(gow,i,f3dxy,dvarx)
        call dy2d(gow,i,f3dyy,dvary)
        f3dv(i) = -(dvarx+dvary)
        sv(i)=sv(i)-(dvarx+dvary)*areap(i)*ramp
      enddo
!$OMP END DO      
    endif

    !Vegetation
    if(veg)then
!$OMP SINGLE        
      val = alphac*0.5*Cdv*Dv
!$OMP END SINGLE      
!$OMP DO PRIVATE(i)       
      do i=1,ncells
        sp(i)=sp(i)-val*nv(i)*(min(hv,h(i))**2)/h(i)*uelwc(i)*areap(i)
      enddo  
!$OMP END DO           
    endif
#endif
    
!!$OMP DO PRIVATE(i,apu,asum)
!    do i=1,ncells
!      !Save variables for momentum interpolation and pressure-correction equation
!      !Notes: These are the same for u and v. Wind and vegetation included in sp term
!      asum=sum(acoef(1:ncface(i),i))
!      apu=iwet(i)*relax(2)/max((asum-sp(i)),1.0e-6) !=1/ap(i)
!      apuareap(i)=apu*areap(i) !Pressure-correction and momentum interpolation
!      sumu(i)=asum*apu        !Pressure-correction equation
!    enddo 
!!$OMP END DO     

!$OMP END PARALLEL    

!#ifdef DEV_MODE
!    !Cross-diffusion terms
!    call flow_crossdiff
!#endif    

    return
    end subroutine coeffsourcesink_uv

!***********************************************************************
    subroutine coefsourcesink_pp()
! Calculates the matrix coefficients, sources, and sinks for the 
! pressure-correction equation 
! written by Weiming Wu, NCCHE, Oct. 2008
! Modified by Alex Sanchez, USACE, 2013
!***********************************************************************
#include "CMS_cpp.h"
    use size_def,   only: ncells, ncellsd, ncelljoint, ncellpoly, nmaxfaces
    use geo_def
    use flow_def
    use der_lib,    only: der_grad_eval,der_gradvec_eval
    use comvarbl,   only: dtime,niter,maxit,relax,ntsch,ctsch,ctsch1,ctsch2,skewcor,rmom
    use cms_def,    only: noptset
    use const_def,  only: small
    use prec_def,   only: ikind
    use interp_def, only: fintp
    
#ifdef DIAG_MODE
    use diag_def
    use diag_lib
#endif

    implicit none
    integer :: i,j,k,ii,nck,jcn
    real(ikind) :: fk,fp,apuvolpf,sumuf,asum,apu
    real(ikind), parameter :: facskew = 0.2 !**********
    !real(ikind) :: wck,wgk
    real(ikind) :: gravdtime,fcorr(nmaxfaces,ncellsD)

    gravdtime = grav*dtime
    
!$OMP PARALLEL

!$OMP DO PRIVATE(i,apu,asum)
    do i=1,ncells
      !Save variables for momentum interpolation and pressure-correction equation
      !Notes: These are the same for u and v. Wind and vegetation included in sp term
      asum=sum(acoef(1:ncface(i),i))
      apu=iwet(i)*relax(2)/max((asum-sp(i)),1.0e-6) !=1/ap(i)
      apuareap(i)=apu*areap(i) !Pressure-correction and momentum interpolation
      sumu(i)=asum*apu        !Pressure-correction equation  
      Hu(i)=iwet(i)*(u(i)+apuareap(i)*dpx(i)*h(i))
        Hv(i)=iwet(i)*(v(i)+apuareap(i)*dpy(i)*h(i))
      !!Hu(i)=iwet(i)*(h(i)*u(i)+apuareap(i)*dpx(i)*h(i)*h(i))
        !!Hv(i)=iwet(i)*(h(i)*v(i)+apuareap(i)*dpy(i)*h(i)*h(i))
    enddo
!$OMP END DO
    
!--- Coefficient Matrix -----------------------
!Variables are interpolated using first order method
!This approach is much more stable and since the pressure correction
!goes to zero at convergence, there is no effect on the results.
!$OMP DO PRIVATE(i,j,k,nck,fk,fp,apuvolpf,sumuf)
      do i=1,ncells
        do j=1,nxyface(i) !No repeat cell faces
          k=kxyface(j,i)    
          nck=cell2cell(k,i) !Forward connectivity
          if(iwet(i)*iwet(nck)==0)then !Closed Boundaries
            acoef(k,i)=0.0
          elseif(nck>ncells)then !Open Boundaries
            acoef(k,i)=hk(k,i)*hk(k,i)*ds(k,i)*apuareap(i)/(1.0-sumu(i))/dn(k,i)
          else
            fk=fintp(k,i); fp=1.0-fk
            apuvolpf=fk*apuareap(nck)+fp*apuareap(i)
            sumuf=fk*sumu(nck)+fp*sumu(i)
            acoef(k,i)=hk(k,i)*hk(k,i)*ds(k,i)*apuvolpf/(1.0-sumuf)/dn(k,i)
            !!acoef(k,i)=hk(k,i)*hk(k,i)*ds(k,i)*(fk*apuareap(nck)/(1.0-sumu(nck))&
            !!          +fp*apuareap(i)/(1.0-sumu(i)))/dn(k,i) !gives almost identical results
          endif
          acoef(llec2llec(k,i),nck)=acoef(k,i) !Matrix is symmetric
        enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL

!--- Flux (Momentum Interpolation) ----------------------------------------
    if(skewcor)then !Inlude second order skewness corrections (linear reconstructions)
    !if(.false.)then !Inlude second order skewness corrections (linear reconstructions)
      call flow_imp_flux_skew     !Modified Rhie and Chow type interpolation with skewness corrections
      !call flow_imp_flux_skew_rc !Original Rhie and Chow type interpolation with skewness corrections
    else !No skewness corrections (no gradients are needed for linear reconstructions)
      call flow_imp_flux_noskew     !Modified Rhie and Chow type interpolation without skewness corrections
      !call flow_imp_flux_noskew_rc !Original Rhie and Chow type interpolation without skewness corrections
    endif

    !Check flux boundaries
#ifdef DIAG_MODE
    do i=1,ncells
      do k=1,ncface(i)
        nck=cell2cell(k,i)  
        if(iwet(i)*iwet(nck)==0 .and. abs(flux(k,i))>1.0e-10)then
          write(msg2,*) ' Cell: ',i  
          write(msg3,*) ' Neighbor: ',nck
          call diag_print_error('Non-zero flux between wet/dry cells at ',msg2,msg3)
        endif
      enddo    
    enddo
#endif
    
!**** ONLY FOR TESTING **************
!    if(noptset>=3)then
!!$OMP DO PRIVATE(i,ii,k,nck,nck1,q2k)
!      do ii=1,ncelljoint
!         i=idcelljoint(ii)  
!         do k=1,ncface(i)-1
!           nck=cell2cell(k,i)
!           nck1=cell2cell(i,k+1)
!           if(idirface(k,i)==idirface(i,k+1) .and. iwet(nck)==1 .and. iwet(nck1)==1)then !Same direction
!             q2k=(flux(k,i)/hk(k,i)+flux(i,k+1)/hk(i,k+1))/2.0
!             flux(k,i)=q2k*hk(k,i)
!             flux(i,k+1)=q2k*hk(i,k+1)
!             flux(llec2llec(k,i),nck)=flux(k,i)
!             flux(nck1,llec2llec(k+1,i))=flux(i,k+1)
!           endif
!         enddo
!      enddo
!!$OMP END DO
!    endif
!**** ONLY FOR TESTING **************

    call bndflux

!--- Source/sink terms -------------------------------------------------    
!$OMP PARALLEL    
    if(ntsch==1)then
!$OMP DO PRIVATE(i)      
      do i=1,ncells
        sp(i)=sppp0(i) !must be non-positive
        su(i)=sspp0(i)-sum(flux(1:ncface(i),i))+p(i)*sp(i)
      enddo      
!$OMP END DO
    else
!$OMP DO PRIVATE(i)      
      do i=1,ncells
        sp(i)=sppp0(i) !must be non-positive
        su(i)=sspp0(i)-sum(flux(1:ncface(i),i))+ctsch*p(i)*sp(i) !Other temporal terms are in sspp0
      enddo      
!$OMP END DO   
    endif

!Deferred skewness correction for pressure correction equation 
!Put into source term. The pressure correction pp here is from
!the previous time step.
    if(skewcor)then
!$OMP DO PRIVATE(i,ii)
      do ii=1,ncelljoint
        i=idcelljoint(ii)
        fcorr(1:ncface(i),i)=0.0
      enddo
!$OMP END DO
!$OMP DO PRIVATE(i)   
      do i=1,ncellpoly
        fcorr(1:ncface(i),i)=0.0
      enddo
!$OMP END DO
!$OMP DO PRIVATE(i,j,ii,k,nck,jcn)
      do ii=1,ncelljoint
         i=idcelljoint(ii)
         if(iwet(i)==0) cycle
         do j=1,nxface(i)
           k=kxface(j,i) 
           nck=cell2cell(k,i)
           if(iwet(nck)==0) cycle
           jcn=llec2llec(k,i)  !Backwards connectivity
           if(abs(rpy(k,i))>1.0e-6)then !coarse to fine
             fcorr(k,i)=-acoef(k,i)*rpy(k,i)*dppy(i)
             fcorr(jcn,nck)=-fcorr(k,i)
           elseif(abs(rpy(jcn,nck))>1.0e-6)then !fine to coarse
             fcorr(k,i)=acoef(k,i)*rpy(jcn,nck)*dppy(nck) !Note: acoef(k,i)=acoef(jcn,nck)
             fcorr(jcn,nck)=-fcorr(k,i)
           endif
         enddo
         do j=1,nyface(i)
           k=kyface(j,i)
           nck=cell2cell(k,i)
           if(iwet(nck)==0) cycle
           jcn=llec2llec(k,i)  !Backwards connectivity
           if(abs(rpx(k,i))>1.0e-6)then !coarse to fine  
             fcorr(k,i)=-acoef(k,i)*rpx(k,i)*dppx(i)
             fcorr(jcn,nck)=-fcorr(k,i)
           elseif(abs(rpx(jcn,nck))>1.0e-6)then !fine to coarse
             fcorr(k,i)=acoef(k,i)*rpx(jcn,nck)*dppx(nck)
             fcorr(jcn,nck)=-fcorr(k,i)
           endif
         enddo
      enddo
!$OMP END DO
!$OMP DO PRIVATE(i,j,k,nck,jcn)    
      do i=1,ncellpoly
         if(iwet(i)==0) cycle
         do j=1,nxyface(i)
           k=kxyface(j,i)  
           nck=cell2cell(k,i)  !Foward connectivity
           if(iwet(nck)==0) cycle
           jcn=llec2llec(k,i)  !Backwards connectivity
           fcorr(k,i)=acoef(k,i)*(rpx(jcn,nck)*dppx(nck)+rpy(jcn,nck)*dppy(nck) &
                      -rpx(k,i)*dppx(i)-rpy(k,i)*dppy(i)) !*Note: acoef(k,i)=acoef(jcn,nck)
           fcorr(jcn,nck)=-fcorr(k,i)
         enddo
      enddo
!$OMP END DO
!$OMP DO PRIVATE(i,ii)
      do ii=1,ncelljoint
        i=idcelljoint(ii)
        su(i)=su(i)+facskew*sum(fcorr(1:ncface(i),i)) !Note: the 0.2 factor is justified because pp will be reduced after each iteration
      enddo
!$OMP END DO
!$OMP DO PRIVATE(i)   
      do i=1,ncellpoly
        su(i)=su(i)+facskew*sum(fcorr(1:ncface(i),i)) !Note: the 0.2 factor is justified because pp will be reduced after each iteration
      enddo
!$OMP END DO
    endif
    
!$OMP END PARALLEL
    
    return
    end subroutine coefsourcesink_pp
    
!***********************************************************************
    subroutine pvfc()
! Correct velocity, mass flux and pressure
! written by Weiming Wu, NCCHE, Oct. 2008
! Modified by Alex Sanchez, USACE-CHL
!***********************************************************************
#include "CMS_cpp.h"
    use size_def
    use geo_def
    use flow_def
    use bnd_def
    use comvarbl
    use wave_flowgrid_def
    use der_def
    use der_lib, only: der_grad_eval,slopelim
    use const_def, only: small
#ifdef DIAG_MODE
    use diag_def
#endif    
    use diag_lib
    use prec_def
    implicit none
    integer :: i,ii,j,k,nck,jcn,iriv,icsh
    real(ikind) :: val !,snode
    
!!--- Wet/dry Boundary ------------------------------
!!Assign pp to dry cells based on average of wet neighboring cells
!!$OMP DO PRIVATE(i,k,nck,jcn,val,snode) 
!    do i=1,ncells
!      if(iwet(i)==1) cycle
!      val=0.0; snode=0.0
!      do k=1,ncface(i)
!        nck=cell2cell(k,i)
!        if(nck>0 .and. nck<=ncells)then
!          if(iwet(nck)==1)then
!            !val=val+pp(nck) !Zero-gradient
!            jcn=llec2llec(k,i)
!            val=val+pp(nck)+2.0*(rx(jcn,nck)*dppx(nck)+ry(jcn,nck)*dppy(nck)) !Linear extrapolation
!            snode=snode+1.0
!          endif
!        endif
!      enddo
!      if(snode>0.9) pp(i)=val/snode !Average of neighboring wet cells
!    enddo
!!$OMP END DO

!!---Ghost Cells -------------
!    do i=ncells+1,ncellsD
!      if(iwet(i)==1) cycle
!      val=0.0; snode=0.0
!      do k=1,ncface(i)
!        nck=cell2cell(k,i)
!        if(nck>0 .and. nck<=ncells)then
!          if(iwet(nck)==1)then
!            val=val+pp(nck)
!            snode=snode+1.0
!          endif
!        endif
!      enddo
!      if(snode>0.9) pp(i)=val/snode !Average of neighboring wet cells
!    enddo    

!--- Wall BC ---------------------------------------------------------
    do j=1,W_str%ncells
      i=W_str%cells(j)
      k=W_str%faces(j)
      pp(cell2cell(k,i))=pp(i)
      !pp(cell2cell(k,i))=pp(i)+2.0*(rx(k,i)*dppx(i)+ry(k,i)*dppy(i)) !Linear extrapolation
    enddo

!--- River Boundary -------------------------------
    do iriv=1,nQstr
      do j=1,Q_STR(iriv)%ncells
        i=Q_STR(iriv)%cells(j)
        k=Q_STR(iriv)%faces(j)
        pp(cell2cell(k,i))=pp(i)
        !pp(cell2cell(k,i))=pp(i)+2.0*(rx(k,i)*dppx(i)+ry(k,i)*dppy(i)) !Linear extrapolation
      enddo
    enddo
      
!--- Cross-shore boundary -------------------------
    do icsh=1,nCSstr
      do j=1,CS_STR(icsh)%ncells
        i=CS_STR(icsh)%cells(j)
        k=CS_STR(icsh)%faces(j)
        pp(cell2cell(k,i))=pp(i)
      enddo
    enddo

    call der_grad_eval(gow,nlim,pp,dppx,dppy) !1-Zero-gradient at wet/dry faces
    !call interp2facegrad(pp,dppx,dppy,ppk)
!!$OMP PARALLEL DO PRIVATE(i)
!    do i=1,ncells
!      dppx(i)=sum(wcx(i,1:ncx(i))*pp(icx(i,1:ncx(i))))
!      dppy(i)=sum(wcy(i,1:ncy(i))*pp(icy(i,1:ncy(i))))
!    enddo 
!!$OMP END PARALLEL DO

!$OMP PARALLEL
!$OMP DO PRIVATE(i,j,k,nck,val)
    do i=1,ncells
      if(iwet(i)==0) cycle  
      val=h(i)*apuareap(i)/(1.0-sumu(i)+small) !Same for u and v
      u(i)=u(i)-val*dppx(i) !Velocity correction
      v(i)=v(i)-val*dppy(i)
      do j=1,nxyface(i)
        k=kxyface(j,i)
        nck=cell2cell(k,i)
        if(iwet(nck)==0) cycle
        flux(k,i)=flux(k,i)-acoef(k,i)*(pp(nck)-pp(i))  !Flux correction 
        flux(llec2llec(k,i),nck)=-flux(k,i)    
      enddo
    enddo
!$OMP END DO

!Skewness correction using linear reconstruction at joint cells
    if(skewcor)then
!$OMP DO PRIVATE(i,ii,j,k,nck,jcn)   
      do ii=1,ncelljoint
        i=idcelljoint(ii)
        if(iwet(i)==0) cycle
        do j=1,nxface(i) !East and West sides
          k=kxface(j,i)
          nck=cell2cell(k,i) !Forward connectivity
          !if(iwet(nck)==0) cycle
          if(nck>ncells .or. iwet(nck)==0) cycle
          jcn=llec2llec(k,i) !Backwards connectivity
          if(abs(rpy(k,i))>1.0e-6)then !Coarse to fine
            flux(k,i)=flux(k,i)+acoef(k,i)*rpy(k,i)*dppy(i)         
            flux(jcn,nck)=-flux(k,i)      
          elseif(abs(rpy(jcn,nck))>1.0e-6)then !Fine to coarse
            flux(k,i)=flux(k,i)-acoef(k,i)*rpy(jcn,nck)*dppy(nck)
            flux(jcn,nck)=-flux(k,i)  
          endif
        enddo
        do j=1,nyface(i) !North and South faces
          k=kyface(j,i)
          nck=cell2cell(k,i) !Forward connectivity
          !if(iwet(nck)==0) cycle
          if(nck>ncells .or. iwet(nck)==0) cycle
          jcn=llec2llec(k,i) !Backwards connectivity
          if(abs(rpx(k,i))>1.0e-6)then !Coarse to fine
            flux(k,i)=flux(k,i)+acoef(k,i)*rpx(k,i)*dppx(i)         
            flux(jcn,nck)=-flux(k,i) 
          elseif(abs(rpx(jcn,nck))>1.0e-6)then !Fine to coarse
            flux(k,i)=flux(k,i)-acoef(k,i)*rpx(jcn,nck)*dppx(nck)
            flux(jcn,nck)=-flux(k,i)    
          endif
        enddo         
      enddo
!$OMP END DO
!$OMP DO PRIVATE(i,j,k,nck,jcn)   
      do i=1,ncellpoly
        if(iwet(i)==0) cycle
        do j=1,nxyface(i) !East and West sides
          k=kxyface(j,i)
          nck=cell2cell(k,i) !Forward connectivity
          !if(iwet(nck)==0) cycle
          if(nck>ncells .or. iwet(nck)==0) cycle
          jcn=llec2llec(k,i) !Backward connectivity
          flux(k,i)=flux(k,i)+acoef(k,i)*(rpx(k,i)*dppx(i)+rpy(k,i)*dppy(i) &
             -rpx(jcn,nck)*dppx(nck)-rpy(jcn,nck)*dppy(nck))
          flux(jcn,nck)=-flux(k,i)
        enddo
      enddo
!$OMP END DO
    endif

!$OMP DO PRIVATE(i)
    do i=1,ncellsD
      p(i)=p(i)+relax(8)*pp(i) !SIMPLEC usually doesn't need relaxation (relax(8)=1.0)
      h(i)=max(hmin,p(i)*gravinv-zb(i)) !Update water depth
      pp(i)=facpp*pp(i) !Initialize for next iteration 
    enddo
!$OMP END DO   
!$OMP END PARALLEL
    
!--- Wall BC ---------------------------------------------------------
    do j=1,W_str%ncells
      i=W_str%cells(j)
      k=W_str%faces(j)
      nck=cell2cell(k,i)
      p(nck)=p(i)
      !p(nck)=p(i)+2.0*(rx(k,i)*dpx(i)+ry(k,i)*dpy(i)) !Linear extrapolation    !Chris' code had this
      !p(nck)=2.0*p(i)-p(cell2cell(llec2llec(k,i),i))  !linear extrapolation
      !h(nck)=max(hmin,p(nck)*gravinv-zb(nck))
    enddo

!---- River boundary ---------------------------
    do iriv=1,nQstr
      do j=1,Q_STR(iriv)%ncells
        i=Q_STR(iriv)%cells(j)
        k=Q_STR(iriv)%faces(j)
        nck=cell2cell(k,i)
        p(nck)=p(i)
        !p(nck)=p(i)+2.0*(rx(k,i)*dpx(i)+ry(k,i)*dpy(i)) !Linear extrapolation    !Chris' code had this
        !p(nck)=2.0*p(i)-p(cell2cell(llec2llec(k,i),i))  !linear extrapolation
        h(nck)=max(hmin,p(nck)*gravinv-zb(nck))
      enddo
    enddo
      
!---- Cross-shore boundary --------------    
    do iriv=1,nCSstr
      do j=1,CS_STR(iriv)%ncells
        i=CS_STR(iriv)%cells(j)
        k=CS_STR(iriv)%faces(j)
        nck=cell2cell(k,i)
        if(flux(k,i)<=0.0)then !Inflow
          p(nck)=p(i)
          !p(nck)=p(i)+2.0*(rx(k,i)*dpx(i)+ry(k,i)*dpy(i)) !Linear extrapolation
          !p(nck)=2.0*p(i)-p(cell2cell(llec2llec(k,i),i))  !linear extrapolation
          h(nck)=max(hmin,p(nck)*gravinv-zb(nck))
        endif
      enddo
    enddo

    !do i=1,ncellsD
    !  h(i)=max(hmin,p(i)*gravinv-zb(i)) !Update water depth
    !enddo
 
    !call flow_wetdry(1)  !Added by Wu
    !call der_update
    !
    !do i=1,ncells          !Added by Wu 
    !  flux(1:ncface(i),i)=iwet(i)*iwet(cell2cell(1:ncface(i),i))*flux(1:ncface(i),i)
    !enddo
    
!!--- Treat wet/dry boundary ------ 
!!!$OMP PARALLEL DO PRIVATE(i,k,nck,val,snode)
!!    do i=1,ncells
!!      if(iwet(i)==1) cycle
!!      val=0.0; snode=0.0
!!      do k=1,ncface(i)
!!        nck=cell2cell(k,i)
!!        if(iwet(nck)==1)then
!!          val=val+p(nck)
!!          snode=snode+1.0
!!        endif
!!      enddo
!!      if(snode>0.9) p(i)=val/max(snode,1.0) !Average of neighboring wet cells
!!    enddo
!!!$OMP END PARALLEL DO

    !Check fluxes and coefficients
#ifdef DIAG_MODE
    do i=1,ncells
      do k=1,ncface(i)
        nck=cell2cell(k,i)
        if(iwet(i)*iwet(nck)==0)then
          if(abs(flux(k,i))>1.0e-15)then
            write(msg2,*) '  Cell: ',mapid(i)
            write(msg3,*) '  Neighbor: ',mapid(nck)
            write(msg4,*) '  flux(k,i)=',flux(k,i)  
            call diag_print_error('Non-zero flux between wet and dry cells after pvfc',&
              msg2,msg3,msg4)
          endif  
          if(nck>ncells .and. abs(acoef(k,i))>1.0e-15)then
            write(msg2,*) '  Cell: ',mapid(i)
            write(msg3,*) '  Neighbor: ',mapid(nck)
            write(msg4,*) '  acoef(k,i)=',acoef(k,i)  
            call diag_print_error('Invalid coefficient between wet and dry cells after pvfc',&
              msg2,msg3,msg4)            
          endif
          !if(nck<=ncells .and. abs(acoef(k,i)-1.0)>1.0e-15)then
          !  write(msg2,*) '  Cell: ',mapid(i)
          !  write(msg3,*) '  Neighbor: ',mapid(nck)
          !  write(msg4,*) '  acoef(k,i)=',acoef(k,i)  
          !  call diag_print_error('Invalid coefficient between wet and dry cells after pvfc',&
          !    msg2,msg3,msg4)
          !endif
        endif
      enddo
    enddo  
#endif

    return
    end subroutine pvfc
    
!*******************************************************
    subroutine flow_grad_interp()
! Calculates the flow gradients and interpolated values
! written by Alex Sanchez, USACE-CHL
!*******************************************************
    use size_def, only: ncells,nnodes
    use geo_def, only: ncface,zbk,areap,dsx,dsy !,cell2cell,llec2llec,nxyface,kxyface
    use flow_def
    use der_def
    use der_lib
    use interp_lib, only: interp_scal_cell2face,interp_scal_cell2node,interp_scal_node2face
    implicit none
    integer :: i,k
    !integer :: j,nck,jcn
    !real*4 :: pn(nnodes)
    
!Notes: 
! In the past the water depth at the cell faces was calculated
! by interpolating using the cell-center water depths and gradients.
! The new approach is to interpolate the water level at the cell
! faces and then calculate the cell-face water depths using
! the cell-face bed elevation.
! Because the bed elevation only changes once every hydro time step,
! the second approach is cheaper.
! It is also more stable and preferable for unstructured meshes where
! the cell-face depth can be obtained from the face nodes.
! In the case of the pressure, it is best to always use 
! the Green-Gauss method because it is conservative.
! The face values are calculated using previous iteration gradients. 
! This approach is also very cheap in that it requires much less
! calculations and no temporary variables.
    
    !--- Current velocity ----
    call der_gradvec_eval(gow,nlim,u,v,dux,duy,dvx,dvy)

    !--- Depth and pressure ------------
    !Linear interpolation and extrapolation of cell-face values (most accurate and inexpensive)
    if(nlim>0) call slopelim(p,nlim,dpx,dpy) !Apply slope limiters here
    call interp_scal_cell2face(p,1,pk,dpx,dpy) !New approach, Note: pk is bounded by p
    !!call interp_scal_cell2face(p,-1,pk,dpx,dpy) !New approach, Note: pk is bounded by p
    !!call interp_scal_cell2face(p,0,pk,dpx,dpy) !New approach, Note: pk is bounded by p
    
!$OMP PARALLEL DO PRIVATE(i,k)
    do i=1,ncells
      dpx(i)=0.0; dpy(i)=0.0 
      do k=1,ncface(i)
        dpx(i) = dpx(i)+dsx(k,i)*pk(k,i)
        dpy(i) = dpy(i)+dsy(k,i)*pk(k,i)
        hk(k,i) = max(hmin,pk(k,i)*gravinv-zbk(k,i))
      enddo
      dpx(i)=dpx(i)/areap(i)
      dpy(i)=dpy(i)/areap(i)
    enddo
!$OMP END PARALLEL DO

    return
    end subroutine flow_grad_interp
    
!***********************************************************************
    subroutine flow_recalculate()
! Reinitializes time step variables and reduces the time step
! by Weiming Wu, NCCHE, Oct. 2008
! modified by Alex Sanchez, USACE-CHL
!***********************************************************************
    use size_def
    use geo_def, only: zb,zbk,ncface,areap,dsx,dsy
    use flow_def
    use comvarbl
    use sed_def
    use sal_def
    use heat_def
    use der_def
    use der_lib, only: der_gradvec_eval
    use interp_lib, only: interp_scal_cell2face
    use diag_lib
    implicit none
    integer :: i,k

    ntsch = 1
    mtime = 1
    jtime = jtime+1
    timesecs = max(dble(stime),timesecs-deltime)
    deltime = dtimebeg/2**jtime !Avoids precision errors
    timesecs = timesecs + deltime
    dtime = real(deltime,kind=ikind)
    ctime = real(timesecs,kind=ikind)
    timehrs = ctime/3600.0 !Bug fix, Alex May 31, 2010, added line
    if(dtime<=dtimebeg/2**8)then
      call write_debug  
      call diag_print_error('Model is divergent',&
        '  Please check model setup and restart run')
    endif
    write(*,*)
    write(*,*) 'Time Step Reduced. Recalculating ...'
    
!--- Hydro --------------   
!$OMP PARALLEL
!$OMP DO PRIVATE(i)   
    do i=1,ncellsD
      u(i)=u1(i)
      v(i)=v1(i)
      p(i)=p1(i)
      h(i)=h1(i)
      pp(i)=0.0
      dppx(i)=0.0
      dppy(i)=0.0
      iwet(i)=iwet1(i)
      flux(:,i)=flux1(:,i)
    enddo
!$OMP END DO    
    
!--- Sediment -------------------
    if(sedtrans)then
!$OMP DO PRIVATE(i)    
      do i=1,ncellsD
        zb(i)=zb1(i)
        Ctk(i,:)=Ctk1(i,:)
        if(ibt>0) btk(i,:)=btk1(i,:)
        if(.not.singlesize)then
          pbk(i,:,1)=pbk1(i,:) !Mixing layer
          db(i,:)=db1(i,:)     !Bed layer thickness
        endif
      enddo        
!$OMP END DO      
    endif
     
!--- Salinity --------------------
    if(saltrans)then
!$OMP DO PRIVATE(i)
      do i=1,ncellsD  
        sal(i)=sal1(i)
      enddo  
!$OMP END DO
    endif

!--- Temperature --------------------
    if(heattrans)then
!$OMP DO PRIVATE(i)
      do i=1,ncellsD  
        heat(i)=heat1(i)
      enddo  
!$OMP END DO
    endif
!$OMP END PARALLEL

    call flow_wetdry(0)
    call bndflux
    
    !call flow_grad_interp
    call der_gradvec_eval(gow,nlim,u,v,dux,duy,dvx,dvy) !nder-1 Cell-based Green-Gauss (1st order for skewed cells)
    call interp_scal_cell2face(p,1,pk) !No second-order corrections
!$OMP PARALLEL DO PRIVATE(i,k)   
    do i=1,ncells
      dpx(i)=0.0; dpy(i)=0.0 
      do k=1,ncface(i)
        dpx(i) = dpx(i)+dsx(k,i)*pk(k,i)
        dpy(i) = dpy(i)+dsy(k,i)*pk(k,i)
        hk(k,i) = max(hmin,pk(k,i)*gravinv-zbk(k,i))
      enddo
      dpx(i)=dpx(i)/areap(i)
      dpy(i)=dpy(i)/areap(i)
    enddo
!$OMP END PARALLEL DO
        
    call fric_eval
    call flow_eddyvis
    
    return
    end subroutine flow_recalculate

!***********************************************************************
    subroutine velface_eta()
! Updates the cell center u and v velocities using face fluxes and 
! calcualtes the water surface elevation from the pressure
! written by Weiming Wu, NCCHE
! modified by Alex Sanchez, USACE-CHL
!***********************************************************************
#include "CMS_cpp.h"
    use comvarbl
    use cms_def
    use bnd_def
    use diag_def
    use diag_lib
    use flow_def
    use geo_def
    use met_def, only: rain,evap
    use prec_def
    use size_def
    use struct_def
    use wave_flowgrid_def
    implicit none
    integer :: i,j,k,iwse,iriv,ibnd
    !integer :: nck
    real(ikind) :: qfluxchk,gravdtime,hold,dhre
    real(ikind) :: uf,vf,umin,vmin,umax,vmax,vel,ctschinv
    !real(ikind) :: fx,fy,sx,sy
    !real(ikind) :: uc,vc,uvnorm,upar,vpar,unorm,vnorm,fac
    !real(ikind) :: icsh,im

!The loop below applies a limiting procedure to the cell-centered 
!current velocities. The velocities are limited to be within the bounds 
!of the cell-faces. This is physically correct and improves the quality 
!of the solution when there are stong source terms (waves). This is
!believed to be caused by the corrector step in the SIMPLEC algorithm.

    if(ncellpoly>0)then
!$OMP PARALLEL DO PRIVATE(i,k,uf,vf,umin,vmin,umax,vmax,vel)  
      do i=1,ncells
        if(h(i)>hmin .and. iwet(i)==1)then
        !if(iwet(i)==1)then
          eta(i)=p(i)*gravinv 
        else
          h(i)=max(hmin,h(i))
          u(i)=0.0; v(i)=0.0
          us(i)=0.0; vs(i)=0.0
          uv(i)=0.0
          eta(i)=-999.0
        endif
        viskfl(1:ncface(i),i)=visk(1:ncface(i),i)
      enddo
!$OMP END PARALLEL DO        
    else    
!$OMP PARALLEL DO PRIVATE(i,k,uf,vf,umin,vmin,umax,vmax) !,vel,fx,fy,sx,sy  
      do i=1,ncells
        if(h(i)>hmin .and. iwet(i)==1)then
        !if(iwet(i)==1)then
          umin= 1.0e6; umax=-1.0e6
          vmin= 1.0e6; vmax=-1.0e6
          do k=1,ncface(i)
            vel=flux(k,i)/(hk(k,i)*ds(k,i))
            uf=fnx(k,i)*vel
            vf=fny(k,i)*vel
            umin=min(umin,uf)
            umax=max(umax,uf)
            vmin=min(vmin,vf)
            vmax=max(vmax,vf)
          enddo
          !if(u(i)>umax .or. u(i)<umin .or. v(i)>vmax .or. v(i)<vmin)then
          !  continue  
          !endif
          u(i)=min(max(u(i),umin),umax)
          v(i)=min(max(v(i),vmin),vmax)
          !if(u(i)>umax .or. u(i)<umin)then
          !  !fx=0.0
          !  !do k=1,ncface(i)
          !  !  fx=fx+fnx(k,i)*flux(k,i)/hk(k,i)
          !  !enddo
          !  !u(i)=0.5*fx/dy(i)
          !  fx=0.0; sy=0.0
          !  do k=1,ncface(i)
          !    fx=fx+fnx(k,i)*flux(k,i)/hk(k,i)
          !    sy=sy+abs(dsy(k,i))
          !  enddo
          !  u(i)=fx/sy
          !endif
          !if(v(i)>vmax .or. v(i)<vmin)then
          !  !fy=0.0
          !  !do k=1,ncface(i)
          !  !  fy=fy+fny(k,i)*flux(k,i)/hk(k,i)
          !  !enddo
          !  !v(i)=0.5*fy/dx(i)
          !  fy=0.0; sx=0.0
          !  do k=1,ncface(i)
          !    fy=fy+fny(k,i)*flux(k,i)/hk(k,i)
          !    sx=sx+abs(dsx(k,i))
          !  enddo
          !  v(i)=fy/sx
          !endif
          eta(i)=p(i)*gravinv
          uv(i)=sqrt((u(i)-us(i))**2+(v(i)-vs(i))**2)
        else
          h(i)=max(hmin,h(i))
          u(i)=0.0; v(i)=0.0
          us(i)=0.0; vs(i)=0.0
          uv(i)=0.0
          eta(i)=-999.0
        endif
        viskfl(1:ncface(i),i)=visk(1:ncface(i),i)
      enddo
!$OMP END PARALLEL DO
    endif
        
    !if(ncellpoly>0 .and. noptset>=3)then
    !  do i=1,ncells
    !    if(waveibr(i)*iwet(i)==1 .and. &
    !      minval(iwet(cell2cell(1:ncface(i),i)))==0)then 
    !      umin= 1.0e6; umax=-1.0e6
    !      vmin= 1.0e6; vmax=-1.0e6
    !      do k=1,ncface(i)
    !        vel=flux(k,i)/hk(k,i)/ds(k,i)  
    !        if(abs(fnx(k,i))>1.0e-5)then
    !          uf=vel/fnx(k,i)
    !          umin=min(umin,uf)
    !          umax=max(umax,uf)
    !        endif
    !        if(abs(fny(k,i))>1.0e-5)then
    !          vf=vel/fny(k,i)
    !          vmin=min(vmin,vf)
    !          vmax=max(vmax,vf)
    !        endif
    !      enddo  
    !      u(i)=min(max(u(i),umin),umax)
    !      v(i)=min(max(v(i),vmin),vmax)     
    !    endif
    !  enddo
    !endif
    
    !if(ncellpoly>0 .and. noptset>=3)then
    !  do i=1,ncells
    !    if(waveibr(i)*iwet(i)>0.3)then !In surf zone
    !      uc = u(i)-us(i) !Current velocity
    !      vc = v(i)-vs(i)
    !      do k=1,ncface(i)
    !        nck=cell2cell(k,i)
    !        if(iwet(nck)==0 .or. nck>ncells)then !Next to dry cell
    !          uvnorm = uc*fnx(k,i)+vc*fny(k,i)  !Wall normal current velocity magnitude
    !          unorm = uvnorm*fnx(k,i)
    !          vnorm = uvnorm*fny(k,i)
    !          upar = uc-unorm  !Wall parallel velocity x-component
    !          vpar = vc-vnorm  !Wall parallel velocity y-component
    !          fac = abs(uvnorm)/sqrt(upar*upar+vpar*vpar)
    !          fac = min(fac,0.1)
    !          uc = upar+unorm*fac
    !          vc = vpar+vnorm*fac
    !          u(i) = uc+us(i)
    !          v(i) = vc+vs(i)
    !        endif
    !      enddo
    !    endif
    !  enddo
    !endif

7897 format(A,I8,A,F9.4,A,F9.4,A,F9.4)
    
!#ifdef DIAG_MODE
!    !k=0
!    !do i=1,ncells
!    !  hold=p(i)*gravinv-zb(i)
!    !  if(iwet(i)==1 .and. hold<0.0)then
!    !    write(msg2,7897), '  Cell: ',mapid(i),', hold=',hold,', h(i)=',h(i),', zb(i)=',zb(i)  
!    !    if(k==0)then
!    !      k = 1  
!    !      call diag_print_warning('Negative water depth',msg2)          
!    !    else
!    !      call diag_print_message(msg2)   
!    !    endif
!    !  endif
!    !enddo
!    k=0
!    do i=1,ncells
!      hold=p(i)*gravinv-zb(i)
!      if(iwet(i)==1 .and. hold<hmin)then
!        if(ncellpoly>0)then
!          write(msg2,7897) '  Cell: ',i,', hold=',hold,', h(i)=',h(i),', zb(i)=',zb(i)
!        else
!          write(msg2,7897) '  Cell: ',mapid(i),', hold=',hold,', h(i)=',h(i),', zb(i)=',zb(i)  
!        endif
!        if(k==0)then
!          k = 1
!          call diag_print_warning('Low wet water depth',msg2)
!        else
!          call diag_print_message(msg2)   
!        endif
!      endif
!    enddo
!#endif
    
#ifdef DEV_MODE
!Eliminate any mass conservation due to numerical precision
!by applying the continuity equation using the final fluxes
    if(volcor)then
      k=0
      dhre = ramp*(rain-evap)*dtime/3600.0
      gravdtime = grav*dtime
      ctschinv = 1.0/ctsch
      open(dgunit,file=dgfile,access='append') 
!$OMP PARALLEL DO PRIVATE(i,hold)          
di:   do i=1,ncells
        !Skip dry cells
        if(iwet(i)==0) cycle di
        !Skip boundary cells
        do ibnd=1,nbndstr
          do j=1,bnd_str(ibnd)%ncells
            if(i==bnd_str(ibnd)%cells(j)) cycle di
          enddo
        enddo
        !Compute water depth based on continuity equation
        hold = h(i)
        if(ntsch==1)then
          h(i) = h1(i) + dhre - dtime*sum(flux(1:ncface(i),i))/areap(i)   
        else
          h(i) = ctschinv*(ctsch1*h1(i) - ctsch2*h2(i) + dhre - dtime*sum(flux(1:ncface(i),i))/areap(i))
        endif
        !Check depth
        if(h(i)<0.0)then
          write(msg2,7897), '  Cell: ',mapid(i),', hold=',hold,', h(i)=',h(i),', zb(i)=',zb(i)
          if(k==0)then
            k = 1  
            call diag_print_warning('Negative water depth after volume correction',msg2)
          else
            call diag_print_message(msg2)
          endif  
          h(i) = hmin
        elseif(h(i)<hmin)then
          h(i) = hmin   
          !iwet(i) = 0
          !do k=1,ncface(i)
          !  flux(k,i) = 0.0
          !  flux(llec2llec(k,i),cell2cell(k,i)) = 0.0
          !enddo
          !u(i) = 0.0; v(i) = 0.0
          !eta(i) = -999.0
          !p(i) = grav*(h(i)+zb(i))
          !h(i) = hmin
        !endif
        elseif(abs(hold/h(i)-1.0)>0.001)then
          write(msg2,7897), '  Cell: ',mapid(i),', hold=',hold,', h(i)=',h(i),', zb(i)=',zb(i)
          if(k==0)then
            k = 1  
            call diag_print_warning('Large volume correction',msg2)
          else
            call diag_print_message(msg2)
          endif
          h(i) = hold  
        else              
          eta(i) = h(i)+zb(i)
          p(i) = eta(i)*grav
          u(i) = u(i)*hold/h(i)
          v(i) = v(i)*hold/h(i)
          uv(i)=sqrt((u(i)-us(i))**2+(v(i)-vs(i))**2)
        endif
      enddo di
!$OMP END PARALLEL DO
      close(dgunit)
    endif
#endif
    
!--- River Boundary Check -----------------
    do iriv=1,nQstr
      qfluxchk=0.0
      do j=1,Q_str(iriv)%ncells
        i=Q_str(iriv)%cells(j)
        k=Q_str(iriv)%faces(j)
        qfluxchk = qfluxchk - flux(k,i)
      enddo
      if (Q_str(iriv)%qflux .ne. 0.0) then                           !Avoid divide by zero  MEB  04/08/2021
        qfluxchk=abs((qfluxchk-Q_str(iriv)%qflux)/Q_str(iriv)%qflux) !Normalized error
        if(qfluxchk>0.001)then
          write(msg2,*) '  Flux Boundary:   ',Q_str(iriv)%idnum
          write(msg3,*) '  Specified Flux:  ',Q_str(iriv)%qflux,' m^3/s'
          write(msg4,*) '  Calculated Flux: ',qfluxchk,' m^3/s' 
          call diag_print_warning('Problem distributing flow at flux boundary ',msg2,msg3,msg4)
        endif
      endif
    enddo

!--- Tidal BC -----------------------------------------------------
    !do iwse=1,nTHstr
    !  call bndvelxtrap(TH_STR(iwse)%ncells,TH_STR(iwse)%cells,TH_STR(iwse)%faces)
    !enddo  

!--- Single Water Level BC -------------------------------      
    !do iwse=1,nHstr
    !  call bndvelxtrap(H_STR(iwse)%ncells,H_STR(iwse)%cells,H_STR(iwse)%faces)
    !enddo

!--- Multiple Water Level BC -------------------------------      
    !do iwse=1,nMHstr
    !  call bndvelxtrap(MH_STR(iwse)%ncells,MH_STR(iwse)%cells,MH_STR(iwse)%faces)
    !enddo
    
!!--- Multiple Water Level and Velocity BC -------------------------------      
!    do iwse=1,nMHstr
!      call bndvelxtrap(MHV_STR(iwse)%ncells,MHV_STR(iwse)%cells,MHV_STR(iwse)%faces)
!    enddo    
    
!--- Cross-shore boundary -------------------------------------------
    !do icsh=1,nCSstr
    !  do im=1,CS_str(icsh)%ncells  
    !    i=CS_str(icsh)%cells(im)
    !    k=CS_str(icsh)%faces(im) 
    !    nck=cell2cell(k,i)
    !    if(mod(idirface(k,i),2)==0) then
    !      if(uxsh(im)*real(idirface(k,i)-3)<=0.0) then
    !        u(nck)=u(i)*h(i)/h(nck)*iwet(i)
    !        v(nck)=v(i)*h(i)/h(nck)*iwet(i)
    !      endif
    !    else
    !      if(vxsh(im)*real(idirface(k,i)-2)<=0.0) then
    !        u(nck)=u(i)*h(i)/h(nck)*iwet(i)
    !        v(nck)=v(i)*h(i)/h(nck)*iwet(i)
    !      endif
    !    endif
    !  enddo
    !enddo

    call struct_velbnd

    return
    end subroutine velface_eta

!******************************************************************      
    subroutine check_conv()
! Checks the implicit hydrodynamic model solver convergence and 
! outputs the variable iconv with one convergence states:
!   iconv = 0 -> Divergent. Reduce time step and recalculate
!   iconv = 1 -> Continue without changing time step
!   iconv = 2 -> Converged. Exit outer loop
!   iconv = 3 -> Continue but with reduced time step
!
! written by Alex Sanchez, USACE-CHL
!******************************************************************    
    use comvarbl
    use solv_def, only: iconv
    use diag_def
    use diag_lib
    use prec_def
    implicit none
    integer :: nhalf !iconv=0 - Divergent, iconv=1 - Continue, iconv=2 - Converged
    real(ikind) :: rmomuv1,rmomuv,rmommax
    !logical :: isnankind
    
455 format(i3,1x,i5,6(3x,1pe12.5))

    iconv = 1
    if(niter<=2)then   
      rmom0(1)=max(rmom(1),1e-15)
      rmom0(2)=max(rmom(2),1e-15)
      rmom0(3)=max(rmom(3),1e-15)
      rmom1=rmom0
    endif  
    rmommax=max(rmom(1),rmom(2),rmom(3))
    
    !if(rmommax>1.0 .or. abs(uxtrm)>velmax .or. abs(vxtrm)>velmax .or.abs(pxtrm)>presmax)then
    !if(rmommax>1.0)then
    if(rmommax>5.0)then
      iconv = 0 !Divergent
      return
    endif
    
!    return  !**********************************
    
    if(niter<=5) return
    
    if(debug_mode)then
      if(rmommax>0.1 .or. abs(uxtrm)>velmax .or. abs(vxtrm)>velmax .or.abs(pxtrm)>presmax)then
        iconv = 0 !Divergent
        return
      endif
    endif
        
    nhalf=maxit/2
        
    !Ratio
    rp=rmom(1)/rmom0(1)
    ruv=max(rmom(2)/rmom0(2),rmom(3)/rmom0(3))
!    if(niter>=nhalf .and. rp<=1.0e-5 .and. ruv<=1.0e-4 .and. rmom(1)<1.0e-5)then
!      iconv=2 !Converged
!      if(debug_mode) write(*,*) 1,niter,rp,ruv
!      return
!    else
    if((rp>=rmomratiop*10.0 .or. ruv>=rmomratiouv*10.0) .and. rmom(1)>rmomtargetp)then
      iconv = 0 !Divergent
      if(debug_mode)then
        write(msg,455) 2,niter,rp,ruv,rmom(1)
        call diag_print_message(msg)
      endif
      return
    endif
!!    if(niter==nhalf .and. (rp>=1.0 .or. ruv>=2.0) .and. rmom(1)>1.0e-4)then
!!      iconv=0 !Divergent
!!      if(debug_mode) write(*,455) 3,niter,rp,ruv,rmom(1)
!!      return      
    if(niter==maxit .and. (rp>=rmomratiop .or. ruv>=rmomratiouv) .and. rmom(1)>rmommaxp)then
      iconv = 0 !Divergent
      if(debug_mode)then
        write(msg,455) 4,niter,rp,ruv,rmom(1)
        call diag_print_message(msg) 
      endif
      return 
    endif    
    if(niter==maxit .and. (rp>=rmomratiop/5.0 .or. ruv>=rmomratiouv/5.0) .and. rmom(1)>rmomtargetp)then
      if(dtvar)then
        iconv = 3 !Reduce time step
      else  
        iconv = 0 !Divergent
      endif  
      if(debug_mode)then
        write(msg,455) 12,niter,rp,ruv,rmom(1)
        call diag_print_message(msg)
      endif  
      return
    endif
    
    !Absolute Error
    rmomuv=max(rmom(2),rmom(3))
    if(rmom(1)<=rmomminp .and. rmomuv<=rmomminuv)then
      iconv = 2 !Converged
      if(debug_mode)then
        write(msg,455) 5,niter,rmom(1),rmomuv
        call diag_print_message(msg)  
      endif  
      return
    endif 
!    elseif(rmom(1)>=1.0e-4 .or. rmomuv>=1.0e-4)then
!      ii=0 !Divergent
!      if(debug_mode) write(*,*) 6,niter,rmom(1),rmomuv
!      return
!!    if(niter>=nhalf .and. (rmom(1)>1.0e-2 .or. rmomuv>=1.0e-1))then
!!      iconv=0 !Divergent
!!      if(debug_mode) write(*,455) 7,niter,rmom(1),rmomuv
!!      return
!!    endif
    if(niter==maxit .and. (rmom(1)>rmommaxp .or. rmomuv>rmommaxuv))then
      iconv = 0 !Divergent
      if(debug_mode)then
        write(msg,455) 8,niter,rmom(1),rmomuv 
        call diag_print_message(msg)      
      endif  
      return
    endif
    if(niter==maxit .and. (rmom(1)>rmomtargetp .or. rmomuv>rmomtargetuv))then
      if(dtvar)then  
        iconv = 3 !Reduce time step
      else  
        iconv = 0 !Divergent
      endif  
      if(debug_mode)then
        write(msg,455) 13,niter,rp,ruv,rmom(1),rmomuv
        call diag_print_message(msg)                
      endif  
      return
    endif
    
    if(niter<=nhalf .or. mod(niter,2)/=0) return
    
    !Absolute Error Change
    rchp=rmom1(1)-rmom(1)
    rmomuv1=max(rmom1(2),rmom1(3))
    rchuv=rmomuv1-rmomuv
    if(rchp<0.0 .or. rchuv<0.0) return
    !if(rchp<=rmomabschgp .and. rchuv<=rmomabschguv)then
    if(rchp<=rmomminp .and. rchuv<=rmomminuv)then
      iconv = 2 !Converged
      if(debug_mode)then
        write(msg,455) 9,niter,rchp,rchuv,rmom(1),rmomuv
        call diag_print_message(msg)     
      endif  
      return
    endif
    
    !!Relative error change    
    !rchp=rchp/rmom1(1)
    !rchuv=rchuv/max(rmom1(2),rmom1(3))
    !if(rchp<=0.05 .and. rchuv<=0.01 .and. rmom(1)<=1.0e-4)then
    !  iconv = 2 !Converged
    !  if(debug_mode)then
    !    open(dgunit,file=dgfile,access='append')  
    !    write(*,455)      10,niter,rchp,rchuv,rmom(1),rmomuv
    !    write(dgunit,455) 10,niter,rchp,rchuv,rmom(1),rmomuv
    !    close(dgunit)
    !  endif  
    !  return
    !endif  
!!    if(rchp>=2.0 .or. rchuv>=5.0)then  
!!      ii=0 !Divergent
!!      if(debug_mode) write(*,455) 11,niter,rchp,rchuv,rmom(1),rmomuv
!!      return
!!    endif    
    
    rmom1=rmom
    
    return
    end subroutine check_conv
    
!**************************************************************************
    subroutine flow_pred
! Predictor step for implicit temporal solution scheme
! Does a simple prediction of the next time step using using an explicit
! formulation which is used as the initial guess for the implicit scheme
! rather than the final solution of the previous time step.
! Still under testing....
! written by Alex Sanchez, USACE-CHL
!**************************************************************************
#include "CMS_cpp.h"
#ifdef DEV_MODE
    use size_def
    use geo_def, only: idirface,ncface,cell2cell,areap,dx,dy,ds,&
        llec2llec,nxyface,kxyface,fnx,fny,dsxy,zb
    use flow_def
    use fric_def, only: cbcfuwcap,wallfric,z0,wallfac,cfrict,uelwc
    use fric_lib, only: wall_coef
    use wave_flowgrid_def, only: worb,worbrep,wper,wang  
    use interp_def, only: fintp
    use comvarbl
    use comp_lib
    use diag_def
    use diag_lib
    use met_def, only: iwndlagr,windconst,windvar,presvar,&
         wndx,wndy,uwind,vwind,cdWndareap,pressatmdx,pressatmdy,&
         tauwx,tauwy,tauwindx,tauwindy
    use wave_flowgrid_def, only: wavestrx,wavestry    
    use cms_def, only: noptset
    use prec_def
    implicit none
    integer :: i,j,k,nck,ierr
    real(ikind) :: sumanu1,sumanv1,sumacoef,delta,vparl,z0wall
    real(ikind) :: rmsh,rmsu,rmsv,forcex,forcey,gammawall
    real(ikind) :: ww(3),psi,cbuwc,uni,unk,val,ddk,fric_bed

    !Calculate RMSE of predicted and calculated values
    rmsh=0.0; rmsu=0.0; rmsv=0.0
    if(pred_corr)then      
      do i=1,ncells
        rmsh=rmsh+(hpred(i)-h(i))**2
        rmsu=rmsu+(upred(i)-u(i))**2
        rmsv=rmsv+(vpred(i)-v(i))**2
      enddo
    else
      do i=1,ncells
        rmsh=rmsh+(h1(i)-h(i))**2
        rmsu=rmsu+(u1(i)-u(i))**2
        rmsv=rmsv+(v1(i)-v(i))**2
      enddo
    endif
    rmsh=sqrt(rmsh/real(ncells,kind=ikind))
    rmsu=sqrt(rmsu/real(ncells,kind=ikind))
    rmsv=sqrt(rmsv/real(ncells,kind=ikind))  
      
787 format(3(3x,1pe12.4))     
    write(msg,787,iostat=ierr) rmsh,rmsu,rmsv
    call diag_print_message('RMSE of predicted values: ',msg)
    
    !Make new prediction
    if(.true.)then
    psi=0.5
    ww(1)=1.0-psi/2.0
    ww(2)=1.0-psi
    ww(3)=psi/2.0
!$OMP PARALLEL    
!$OMP DO PRIVATE(i,cbuwc,val)      
    do i=1,ncells
      !hpred(i)=h(i)-dtime*(h(i)*dux(i)+u(i)*dhx(i)+h(i)*dvy(i)+v(i)*dhy(i))
      hpred(i)=h(i)-iwet(i)*dtime/areap(i)*sum(flux(1:ncface(i),i)) !Water depth prediction
      !hpred(i)=1.0/ww(1)*(h(i)*ww(2)+h1(i)*ww(3) &
      !  -iwet(i)*dtime/areap(i)*sum(flux(1:ncface(i),i))) !Water depth prediction
      !if(h(i)<1.0e-6)then
      !  hpred(i)=1.0e-6
      !endif
      if(hpred(i)<hmin)then
        hpred(i)=hmin  
        upred(i)=0.0
        vpred(i)=0.0
      else
        !call d2xy2d(i,dux,duy,d2uix2,d2uiy2)
        !call d2xy2d(i,dvx,dvy,d2vix2,d2viy2)
        !upred(i)=u(i)+dtime*(-u(i)*dux(i)-v(i)*duy(i)+forcex(i)/h(i)+fc(i)*v(i)-dpx(i)+vis(i)*(d2uix2+d2uiy2)) !Note: Mixing terms
        !vpred(i)=v(i)+dtime*(-u(i)*dvx(i)-v(i)*dvy(i)+forcey(i)/h(i)-fc(i)*u(i)-dpy(i)+vis(i)*(d2vix2+d2viy2))
        upred(i)=u(i)+dtime*(-u(i)*dux(i)-v(i)*duy(i)+forcex(i)/h(i)+fc(i)*v(i)-dpx(i)) !Note: Mixing terms
        vpred(i)=v(i)+dtime*(-u(i)*dvx(i)-v(i)*dvy(i)+forcey(i)/h(i)-fc(i)*u(i)-dpy(i))
        !upred(i)=upred(i)/(1.0+cfrict(i)*uelwc(i)/h(i)) !Bottom friction treated semi-implicitly for stability
        !vpred(i)=vpred(i)/(1.0+cfrict(i)*uelwc(i)/h(i)) !Bottom friction treated semi-implicitly for stability
        !upred(i)=1.0/ww(1)*(u(i)*ww(2)+u1(i)*ww(3) &
        !    -dtime*(u(i)*dux(i)+v(i)*duy(i)-forcex(i)/h(i)-fc(i)*v(i)+dpx(i)))        
        !vpred(i)=1.0/ww(1)*(v(i)*ww(2)+v1(i)*ww(3) &
        !    -dtime*(u(i)*dvx(i)+v(i)*dvy(i)-forcey(i)/h(i)+fc(i)*u(i)+dpy(i)))
        if(noptset>=3)then
          cbuwc=fric_bed(h(i),cfrict(i),z0(i),vpred(i),vpred(i),&
                us(i),vs(i),worb(i),worbrep(i),wper(i),wang(i)) 
        else
          cbuwc=cfrict(i)*sqrt(upred(i)*upred(i)+vpred(i)*vpred(i))
        endif
        val=1.0+dtime*cbuwc/h(i)
        upred(i)=upred(i)/val !Bottom friction treated semi-implicitly for stability
        vpred(i)=vpred(i)/val !Bottom friction treated semi-implicitly for stability
      endif
    enddo      
!$OMP END DO   
!$OMP END PARALLEL
    else !Conservative form
!$OMP PARALLEL
!$OMP DO PRIVATE(i)      
    do i=1,ncells
      su(i)=-cbcfuwcap(i)*u(i)+( fc(i)*v(i)-dpx(i))*h(i)*areap(i)
      sv(i)=-cbcfuwcap(i)*v(i)+(-fc(i)*u(i)-dpy(i))*h(i)*areap(i)
      sp(i)=-cbcfuwcap(i)
      do k=1,ncface(i)
        ddk=visk(k,i)*hk(k,i)*dsxy(k,i)
        acoef(k,i)=hybridcoef(ddk,flux(k,i))
      enddo
    enddo  
!$OMP END DO   

    !Wind (Lagrangian Reference Frame)
    if(iwndlagr==1)then
      if(windconst)then
!$OMP DO PRIVATE(i)  
        do i=1,ncells
          su(i)=su(i)+cdWndareap(i)*(wndx-u(i)+us(i)) !source/sink
          sv(i)=sv(i)+cdWndareap(i)*(wndy-v(i)+vs(i)) !source/sink
          sp(i)=sp(i)-cdWndareap(i)      !sink, must be non-positive
        enddo
!$OMP END DO        
      elseif(windvar)then
!$OMP DO PRIVATE(i)  
        do i=1,ncells
          su(i)=su(i)+cdWndareap(i)*(uwind(i)-u(i)+us(i)) !source/sink
          sv(i)=sv(i)+cdWndareap(i)*(vwind(i)-v(i)+vs(i)) !source/sink
          sp(i)=sp(i)-cdWndareap(i)      !sink, must be non-positive
        enddo
!$OMP END DO
      endif
    else !Eulerian Referenc Frame
      if(windconst)then
!$OMP DO PRIVATE(i)    
        do i=1,ncells    
          su(i)=su(i)+tauwx*areap(i)
          sv(i)=sv(i)+tauwy*areap(i)
        enddo
!$OMP END DO  
      elseif(windvar)then    
!$OMP DO PRIVATE(i)    
        do i=1,ncells    
          su(i)=su(i)+tauwindx(i)*areap(i)
          sv(i)=sv(i)+tauwindy(i)*areap(i)
        enddo
!$OMP END DO        
      endif        
    endif   
   
    !Wave forcing
    if(noptset>=3)then
!$OMP DO PRIVATE(i)       
      do i=1,ncells
        su(i)=su(i)+wavestrx(i)*areap(i)
        sv(i)=sv(i)+wavestry(i)*areap(i)
      enddo    
!$OMP END DO 
      !Wave velocity forcing
      if(waveflux)then
!$OMP DO PRIVATE(i)  
        do i=1,ncells       
          su(i)=su(i)+cbcfuwcap(i)*us(i) 
          sv(i)=sv(i)+cbcfuwcap(i)*vs(i) 
        enddo
!$OMP END DO
      endif
    endif  
    
    !Atmospheric pressure gradients
    if(presvar)then
!$OMP DO PRIVATE(i)  
      do i=1,ncells       
        su(i)=su(i)-pressatmdx(i)*h(i)*areap(i)/rhow
        sv(i)=sv(i)-pressatmdy(i)*h(i)*areap(i)/rhow
      enddo
!$OMP END DO  
    endif
!$OMP END PARALLEL

    !Deferred Corrections
    select case(ndsch) !Anti-diffusion corrections
    case(5)                                                      !Chris 4
      call defcorhlpagrad(u,dux,duy,su)
      call defcorhlpagrad(v,dvx,dvy,sv)
    case(6)                                                      !Chris 5
      call defcorgammagrad(gammadefcor,u,dux,duy,su)
      call defcorgammagrad(gammadefcor,v,dvx,dvy,sv)
    case(7)                                                      !Chris 6
      call defcorgammagrad(cubistadefcor,u,dux,duy,su)
      call defcorgammagrad(cubistadefcor,v,dvx,dvy,sv)
    case(8)                                                      !Chris 7
      call defcorgammagrad(alvsmartdefcor,u,dux,duy,su)
      call defcorgammagrad(alvsmartdefcor,v,dvx,dvy,sv)  
    case(9)                                                      !Chris 8
      call defcorgammagrad(hoabdefcor,u,dux,duy,su)
      call defcorgammagrad(hoabdefcor,v,dvx,dvy,sv)
    end select  
    if(skewcor)then !Skewness corrections
      call defcorparagradvec(dux,duy,dvx,dvy,su,sv)
    endif
    
!--- Dry nodes --------------------------------------------------------
    if(wallfric)then
!$OMP PARALLEL DO PRIVATE(i,k,nck,delta,vparl,z0wall,gammawall)            
      do i=1,ncells
        if(iwet(i)==1)then
           spu(i)=sp(i)
           spv(i)=sp(i) 
           do k=1,ncface(i)
             nck=cell2cell(k,i)
             if(iwet(nck)==0)then    !side is dry            
               !vparl=sqrt(u(i)*u(i)+v(i)*v(i))    
               z0wall=z0(i)*wallfac
               if(idirface(k,i)==1.or.idirface(k,i)==3)then   !north/south face
                 vparl=abs(u(i))                                  
                 delta=0.5*dy(i)
                 !call wall_gamma(viscos,delta,z0wall,vparl,gammawall) !Valid for smooth to rough flow
                 gammawall = wall_coef(delta,z0wall)*vparl
                 spu(i)=spu(i)-gammawall*ds(k,i)*h(i)
               else                 !east/west face
                 vparl=abs(v(i))
                 delta=0.5*dx(i)
                 !call wall_gamma(viscos,delta,z0wall,vparl,gammawall) !Valid for smooth to rough flow
                 gammawall = wall_coef(delta,z0wall)*vparl
                 spv(i)=spv(i)-gammawall*ds(k,i)*h(i)
               endif
               acoef(k,i)=0.0
             endif
           enddo
        endif
      enddo
!$OMP END PARALLEL DO
    else !No wall friction
!$OMP PARALLEL DO PRIVATE(i)
      do i=1,ncells
        spu(i)=sp(i)
        spv(i)=sp(i)
      enddo     
!$OMP END PARALLEL DO      
    endif    

!$OMP PARALLEL DO PRIVATE(i,k,sumanu1,sumanv1,sumacoef)  
    do i=1,ncells
      !Water level
      hpred(i)=h(i)-iwet(i)*dtime/areap(i)*sum(flux(1:ncface(i),i)) !Water depth prediction
      !Current velocities
      ap(i)=hpred(i)*areap(i)/dtime  
      sumanu1=0.0; sumanv1=0.0; sumacoef=0.0
      do k=1,ncface(i)
        nck=cell2cell(k,i)  
        sumanu1=sumanu1+acoef(k,i)*u(nck)
        sumanv1=sumanv1+acoef(k,i)*v(nck)
        sumacoef=sumacoef+acoef(k,i)
      enddo
      upred(i)=iwet(i)*(sumanu1+(ap(i)-sumacoef)*u(i)+su(i))/(ap(i)-spu(i))
      vpred(i)=iwet(i)*(sumanv1+(ap(i)-sumacoef)*v(i)+sv(i))/(ap(i)-spv(i))
    enddo
!$OMP END PARALLEL DO
    endif

!!$OMP PARALLEL DO PRIVATE(i)
!    do i=1,ncells  
!      h(i)=hpred(i)
!      eta(i)=h(i)+zb(i)
!      p(i)=grav*eta(i)
!      u(i)=upred(i)
!      v(i)=vpred(i)
!    enddo      
!!$OMP END PARALLEL DO
!
!!$OMP PARALLELDO PRIVATE(i,j,k,nck,uni,unk)      
!    do i=1,ncells      
!      do j=1,nxyface(i)
!         k=kxyface(j,i)
!         nck=cell2cell(k,i)
!         hk(k,i)=fintp(k,i)*h(nck)+(1.0-fintp(k,i))*h(i)
!         uni=fnx(k,i)*u(i)+fny(k,i)*v(i)
!         unk=fnx(k,i)*u(nck)+fny(k,i)*v(nck)
!         flux(k,i)=iwet(i)*iwet(nck)*ds(k,i)*hk(k,i)*(fintp(k,i)*unk+(1.0-fintp(k,i))*uni) !Outward flux
!         flux(llec2llec(k,i),nck)=-flux(k,i)         
!       enddo
!    enddo
!!$OMP END PARALLEL DO

#endif
    return
    end subroutine flow_pred

!***************************************************
    subroutine check_momentum
!***************************************************
    use size_def
    use geo_def, only: mapid,ncface,cell2cell
    use flow_def
    use comvarbl, only: flowpath
    use solv_def, only: iconv
    use diag_lib
    implicit none
    integer :: i,k,idunit
    logical :: isnankind
    
    idunit = 999
    iconv = 1
    do i=1,ncells
      if(isnankind(su(i)))then
        call diag_print_warning('Source term su(i) = NaN','for U velocity component')
        call diag_print_var(i,2)
        iconv=0
      endif
      if(isnankind(sv(i)))then
        call diag_print_warning('Source term sv(i) = NaN','for V velocity component')
        call diag_print_var(i,3)
        iconv=0
      endif
      if(isnankind(spu(i)))then
        call diag_print_warning('Source term spu(i) = NaN','for U velocity component')
        call diag_print_var(i,2)
        iconv=0
      endif
      if(isnankind(spv(i)))then
        call diag_print_warning('Source term spv(i) = NaN','for V velocity component')
        call diag_print_var(i,3)
        iconv=0
      endif
      do k=1,ncface(i)
        if(isnankind(acoef(k,i)))then
          call diag_print_warning('Coefficient term acoef(k,i) = NaN','for velocity fields')  
          call diag_print_var(i,2)
          iconv=0
        endif
        if(iwet(i)*iwet(cell2cell(k,i))==0)then
          if(abs(flux(k,i))>1.0e-15)then
            call diag_print_warning('Non-zero flux between at dry boundary')
            call diag_print_var(i,2)
            iconv=0
          endif
          if(abs(acoef(k,i))>=1.0e-15)then
            call diag_print_warning('Non-zero coefficient at dry boundary for velocity fields')
            call diag_print_var(i,2)
            iconv=0
          endif
        endif
      enddo !k
    enddo !i    
    if(iconv==0)then       
      call write_debug
    endif
    
    return
    end subroutine check_momentum    
    
!***************************************************
    subroutine check_variables(neq)
!***************************************************
    use size_def
    use geo_def, only: mapid,ncface
    use flow_def
    use comvarbl, only: flowpath
    use solv_def, only: iconv
    use diag_lib
    implicit none
    integer :: i,k,idunit,neq
    logical :: isnankind,debug_write
    
    idunit = 999
    iconv = 1
    debug_write = .false.
    do i=1,ncells
      if(isnankind(su(i)))then
        write(*,*) 'su(i)=NaN'
        call diag_print_var(i,neq)
        iconv=0
      endif
      if(isnankind(sp(i)))then
        write(*,*) 'sp(i)=NaN'
        call diag_print_var(i,neq)
        iconv=0
      endif
      do k=1,ncface(i)
        if(isnankind(acoef(k,i)))then
          write(*,*) 'acoef(k,i)=NaN'
          call diag_print_var(i,neq)
          iconv=0
        endif
      enddo !k
      if(iconv==0 .and. .not.debug_write)then       
        call write_debug       
        !write(*,*) 'Press <RETURN> to continue'
        !read(*,*)    
        !stop
        debug_write = .true.
      endif
    enddo !i    

    return
    end subroutine check_variables
    
!****************************************************    
    subroutine flow_crossdiff
! Calculates the cross-diffusion terms and adds them
! to the source term
! written by Alex Sanchez, USACE-CHL
!****************************************************    
    use prec_def
    use size_def, only: ncells,ncellsD
    use geo_def, only: areap
    use flow_def, only: iwet,vis,h,dux,duy,dvx,dvy,su,sv
    use der_def
    implicit none
    integer :: i,k,nck
    real(ikind) :: vhux(ncellsD),vhvx(ncellsD)
    real(ikind) :: vhuy(ncellsD),vhvy(ncellsD)
    real(ikind) :: dvhuxx,dvhvxy,dvhuyx,dvhvyy
    real(ikind) :: val,wck

!$OMP PARALLEL
!$OMP DO PRIVATE(i,val)
    !Temporary arrays
    do i=1,ncells
      val = vis(i)*h(i)
      vhux(i) = val*dux(i)
      vhvx(i) = val*dvx(i)
      vhuy(i) = val*duy(i)
      vhvy(i) = val*dvy(i)
    enddo
!$OMP END DO 
!$OMP DO PRIVATE(i,k,nck,wck,dvhuxx,dvhvxy,dvhuyx,dvhvyy)
    do i=1,ncells
      dvhuxx=0.0; dvhvxy=0.0
      dvhuyx=0.0; dvhvyy=0.0
      do k=1,gow%ncx(i)
        nck = gow%icx(k,i)
        wck = gow%wcx(k,i)
        if(nck>ncells .or. iwet(nck)==0)then
          dvhuxx = dvhuxx + wck*vhux(nck)
          dvhuyx = dvhuyx + wck*vhuy(nck)
        else
          dvhuxx = dvhuxx + wck*vhux(i)
          dvhuyx = dvhuyx + wck*vhuy(i)
        endif
      enddo !k
      do k=1,gow%ncy(i)
        nck = gow%icy(k,i)
        wck = gow%wcy(k,i)
        if(nck>ncells .or. iwet(nck)==0)then
          dvhvxy = dvhvxy + wck*vhvx(nck)
          dvhvyy = dvhvyy + wck*vhuy(nck)
        else
          dvhvxy = dvhvxy + wck*vhvx(i)
          dvhvyy = dvhvyy + wck*vhuy(i)  
        endif
      enddo !k
      su(i) = su(i) + (dvhuxx+dvhvxy)*areap(i)
      sv(i) = sv(i) + (dvhvyy+dvhuyx)*areap(i)
    enddo
!$OMP END DO
!$OMP END PARALLEL

    return
    end subroutine flow_crossdiff
    
!******************************************************************************    
    subroutine flow_volres(volerr)
! Checks the local volume balance at each cell and outputs the error as 
! a water volume.
! written by Alex Sanchez, USACE-CHL
!******************************************************************************
    use comvarbl, only: ntsch,ctsch,ctsch1,ctsch2,dtime,ramp
    !!use const_def, only: gravinv  
    use flow_def, only: h,h1,h2,flux
    use geo_def, only: areap,ncface 
    use met_def, only: rain_evap,rain,evap
    use prec_def
    use size_def, only: ncells,ncellsD
    implicit none
    integer :: i
    real(ikind) :: volerr(ncellsD),deta
    
!--- Local Volume Balance ------------------    
    if(ntsch==1)then
!$OMP PARALLEL DO PRIVATE(i)         
      do i=1,ncells
        volerr(i) = dtime*sum(flux(1:ncface(i),i)) &
          + areap(i)*(h(i)-h1(i))
      enddo
!$OMP END PARALLEL DO         
    else
!$OMP PARALLEL DO PRIVATE(i)        
      do i=1,ncells
        volerr(i) = dtime*sum(flux(1:ncface(i),i)) &
           + areap(i)*(ctsch*h(i)-ctsch1*h1(i)+ctsch2*h2(i))
      enddo
!$OMP END PARALLEL DO      
    endif
    
    if(rain_evap)then
      deta = ramp*(rain-evap)*dtime/3600.0  
!$OMP PARALLEL DO PRIVATE(i)        
      do i=1,ncells
        volerr(i) = volerr(i) - deta*areap(i)
      enddo
!$OMP END PARALLEL DO      
    endif
    
    return
    end subroutine flow_volres
    
!******************************************************************************    
    subroutine flow_balance()
! Checks the global water volume balance 
! written by Alex Sanchez, USACE-CHL
!******************************************************************************
    use comvarbl, only: ntsch,ctsch,ctsch1,ctsch2,dtime,ramp
    use flow_def, only: h,h0,h1,h2,flux,u,v,iwet,iwet1,gravinv,&
       volH2Ocumstg,volH2Ocumbnd,volH2Ocumrain,volH2Ocumevap,volH2Ocumnorm,volH2Ocumflux
    use geo_def, only: areap,ncface,cell2cell,fnx,fny,ds
    use met_def, only: rain_evap,rain,evap
    use prec_def
    use size_def, only: ncells,ncellsD
    use bnd_def, only: nbndstr,bnd_str
    use diag_def
    use diag_lib
    implicit none
    integer :: i,j,k,ibnd,nck
    real(ikind) :: dhrain,dhevap,volH2O,vel
    real(ikind) :: volH2Ostg,volH2Obnd,volH2Oflux,volH2Orain,volH2Oevap
    real(ikind) :: volH2Oerr,volH2Ocumerr,volH2Obal,volH2Ocumbal,volH2Onet,volH2Ocumnet
    
    !Initialize variables
    volH2Ostg = 0.0     !Current time step global volume balance for storage term
    volH2Obnd = 0.0     !Current time step global volume net efflux
    volH2Orain = 0.0    !Current time step global volume balance due to rain
    volH2Oevap = 0.0    !Current time step global volume balance due to evaporation
    volH2Oflux = 0.0
    !volH2Ocumstg = 0.0
    
    if(ntsch==1)then   
      do i=1,ncells
        !if(iwet(i)==0) cycle
        !if(iwet(i)*iwet1(i)==0) cycle
        volH2Ostg = volH2Ostg + areap(i)*(h(i)-h1(i))
        !volH2Ocumstg = volH2Ocumstg + areap(i)*(h(i)-h0(i))
        volH2Oflux = volH2Oflux + dtime*sum(flux(1:ncface(i),i))
      enddo       
    else    
      do i=1,ncells
        !if(iwet(i)==0) cycle
        !if(iwet(i)*iwet1(i)==0) cycle
        volH2Ostg = volH2Ostg + areap(i)*(ctsch*h(i)-ctsch1*h1(i)+ctsch2*h2(i))
        !volH2Ocumstg = volH2Ocumstg + areap(i)*(h(i)-h0(i))
        volH2Oflux = volH2Oflux + dtime*sum(flux(1:ncface(i),i))
      enddo  
    endif
    volH2Ocumstg = volH2Ocumstg + volH2Ostg
    volH2Ocumflux = volH2Ocumflux + volH2Oflux
    
    !All Boundaries
    do ibnd=1,nbndstr
      do j=1,bnd_str(ibnd)%ncells
        i = bnd_str(ibnd)%cells(j)
        !if(iwet(i)==0) cycle
        k = bnd_str(ibnd)%faces(j)
        nck = cell2cell(k,i)
        vel=(fnx(k,i)*u(i)+fny(k,i)*v(i))
        !volH2O = dtime*ds(k,i)*h(i)*vel
        volH2O = dtime*flux(k,i)
        volH2Obnd = volH2Obnd + volH2O
      enddo
    enddo  
    volH2Ocumbnd = volH2Ocumbnd + volH2Obnd
    
    if(rain_evap)then
      dhrain = ramp*rain*dtime/3600.0  
      dhevap = -ramp*evap*dtime/3600.0
      do i=1,ncells
        !if(iwet(i)==0) cycle
        volH2Orain = volH2Orain + dhrain*areap(i)
        volH2Oevap = volH2Oevap + dhevap*areap(i)        
      enddo
      volH2Ocumrain = volH2Ocumrain + volH2Orain
      volH2Ocumevap = volH2Ocumevap + volH2Oevap
    endif
        
    !Continuity equation -> storage + flux + rain - evaporation = 0
    volH2Obal = volH2Ostg + volH2Obnd - volH2Orain + volH2Oevap  !Current time step
    volH2Onet = volH2Ostg + volH2Oflux - volH2Orain + volH2Oevap  !Current time step
    volH2Ocumbal = volH2Ocumstg + volH2Ocumbnd - volH2Ocumrain + volH2Ocumevap  !Cumulative
    volH2Ocumnet = volH2Ocumstg + volH2Ocumflux - volH2Ocumrain + volH2Ocumevap  !Cumulative
    
    !Percentage error
    !volH2O = 0.0
    !do i=1,ncells
    !  volH2O = h(i)*areap(i)
    !enddo  
    volH2O = volH2Orain - volH2Oevap + abs(volH2Ostg) + abs(volH2Obnd) !Normalizing volume (magnitude)    
    volH2Oerr = 100.0*volH2Obal/max(volH2O,1.0e-10)
    !volH2O = volH2Ocumrain - volH2Ocumevap + abs(volH2Ocumstg) + abs(volH2Ocumbnd) !Normalizing volume (magnitude)
    volH2Ocumnorm = max(volH2Ocumnorm,abs(volH2Ocumstg))
    volH2O = volH2Ocumnorm
    volH2Ocumerr = 100.0*volH2Ocumbal/max(volH2O,1.0e-10)
    
    !Print results               
751 format(' Water Balance Volumes, cu m')
752 format('              Storage, Net Efflux, Evaporation, Precipitation, Balance    ')
753 format(' Current   ',1x,1pe10.3,2x,1pe10.3,2x,1pe10.3,4x,1pe10.3,2x,1pe10.3)
754 format(' Cumulative',1x,1pe10.3,2x,1pe10.3,2x,1pe10.3,4x,1pe10.3,2x,1pe10.3)   
844 format(' Current Water Percent Mass Error: ',1x,1pe10.3)    
855 format(' Cumulative Water Percent Mass Error: ',1x,1pe10.3)    
    
    write(msg2,751)
    write(msg2,752)
    write(msg3,753) volH2Ostg,volH2Obnd,volH2Oevap,volH2Orain,volH2Obal
    write(msg3,753) volH2Ostg,volH2Oflux,volH2Oevap,volH2Orain,volH2Onet
    write(msg4,754) volH2Ocumstg,volH2Ocumbnd,volH2Ocumevap,volH2Ocumrain,volH2Ocumbal
    write(msg5,754) volH2Ocumstg,volH2Ocumflux,volH2Ocumevap,volH2Ocumrain,volH2Ocumnet
    write(msg6,844) volH2Oerr   
    write(msg7,855) volH2Ocumerr  
    call diag_print_message(' ',msg2,msg3,msg4,msg5,msg6,msg7) 
    
    return
    end subroutine flow_balance    
        
!****************************************************    
    subroutine flow_rainevap()
! Adds the precipitation (rain) and evaporation
! to the water surface elevation and corrects the
! current velocities to conserve fluxes
! This is the simplist approach for dealing with the
! surface water fluxes. An alternative approach is
! to put source/sink terms in the continuity equation
! written by Alex Sanchez, USACE-CHL 
! modified for rainfall by Chris Reed, 06/28/2016
!****************************************************
    use comvarbl, only: dtime,ramp
    use geo_def, only: zb
    use flow_def, only: eta,p,u,v,h,iwet,hmin,grav,gravinv
    use met_def, only: rain,evap
    use prec_def
    use size_def, only: ncellsD
    use met_def, only: rain_time,rf_unit         !06/28/2016
    use sal_def, only: saltrans,sal,sal1         !06/28/2016
    implicit none
    integer :: i
    real(ikind) :: dpre,dpre_evap,hold
    
    rain_time = rain_time + dtime
    if(rain_time .ge. 3600) then
        rain_time = 0.0
        read(rf_unit,*)rain,evap
    endif

    dpre_evap = grav*ramp*(rain-evap)*dtime/3600.0
    dpre      = grav*ramp*rain*dtime/3600.0

    if(saltrans) then    
!$OMP PARALLEL DO PRIVATE(i,hold)
      do i=1,ncellsD
        if(iwet(i)==1) then
          p(i)=p(i)+dpre_evap
        else
          p(i)=p(i)+dpre
        endif
        eta(i)=p(i)*gravinv
        hold=h(i)
        h(i)=max(hmin,p(i)*gravinv-zb(i))
        u(i)=u(i)*hold/h(i)
        v(i)=v(i)*hold/h(i)
        sal1(i) = sal1(i)*hold/h(i)
      enddo
!$OMP END PARALLEL DO
    else
!$OMP PARALLEL DO PRIVATE(i,hold)
      do i=1,ncellsD
        if(iwet(i)==1) then
          p(i)=p(i)+dpre_evap
        else
          p(i)=p(i)+dpre
        endif
        eta(i)=p(i)*gravinv
        hold=h(i)
        h(i)=max(hmin,p(i)*gravinv-zb(i))
        u(i)=u(i)*hold/h(i)
        v(i)=v(i)*hold/h(i)
      enddo
!$OMP END PARALLEL DO
    endif    
    
    return
    end subroutine flow_rainevap