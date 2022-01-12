!===============================================================================
! CMS Sediment Transport Total-load Correction Factor Routines
!  
!  Contains the following:
!     btklogexp  - Calculates the total-load correction factor
!                  based on a logorithmic current velocity and an 
!                  exponential sediment concentration vertical profiles
!     btklogrouse - Calculates the total load correction factor
!                  based on a logarithmic current velocity and 
!                  Rouse sediment concentration vertical profiles.
!
! written by Alex Sanchez, USACE-CHL
!===============================================================================

!***************************************************************************   
    subroutine btklogexp()
! Calculates the total-load correction factor
! based on a logorithmic current velocity and an 
! exponential concentration vertical profiles
! written by Alex Sanchez, USACE-CHL
!***************************************************************************    
    use size_def    
    use flow_def, only: vis,h,uv,iwet,rhow,grav
    use fric_def, only: z0,coefman,cfrict,uelwc,bsxy
    use sed_def, only: nsed,btk,wsfall,schmidt,taucr,diam,rsk,&
        d50,dstar,pbk,CtstarP,bsk,ubk,epsvk,bsk,specgrav
    use sed_lib, only: bedvel_vanrijn_wu
    use beta_lib, only: bslogexp_table !,bslogexp
    use fric_lib, only: fric_normapprough
    use prec_def
    implicit none
    integer :: i,ks
    real(ikind) :: zap,phip,taup,cnp,bbk
  
!$OMP PARALLEL DO PRIVATE(i,KS,zap,phip,bbk,taup,cnp)
    do i=1,ncells
      if(iwet(i)==0)then
        btk(i,:)=1.0
        bsk(i,:)=1.0
        ubk(i,:)=0.0
      else
        zap=fric_normapprough(uelwc(i),uv(i),cfrict(i)) !Normalized apparent roughness
        !zap=z0(i)/h(i)  !Only for no waves
        cnp=0.05*d50(i)**0.16666667                !Manning corresponding to grain roughness
        if (coefman(i) .eq. 0) then
          taup=bsxy(i)*min(cnp/1.0e-6,1.0)**1.5      !If coefman(i) is zero, use very small number instead   MEB 01/12/2022
        else
          taup=bsxy(i)*min(cnp/coefman(i),1.0)**1.5  !Grain shear stress, min is used in case cn is less than cnp
        endif          
        do ks=1,nsed
          !phip=wsfall(ks)*h(i)/(vis(i)/schmidt) 
          phip=wsfall(ks)*h(i)/max(epsvk(i,ks),1.0e-6) !More correct
          !!bsk(i,ks)=bslogexp(zap,phip)               !Suspended load correction factor, analytical
          bsk(i,ks)=bslogexp_table(zap,phip)           !Suspended load correction factor, table
          ubk(i,ks)=bedvel_vanrijn_wu(grav,taup,taucr(ks),specgrav,diam(ks))   !Bed load velocity, m/s
          bbk=min(ubk(i,ks)/max(uv(i),1.0e-6),1.0)      !physical restriction, bed load cannot move faster than current velocity
          bbk=max(bbk,0.2)    
          btk(i,ks)=1.0/(rsk(i,ks)/bsk(i,ks)+(1.0-rsk(i,ks))/bbk) !Total load correction factor
          btk(i,ks)=max(btk(i,ks),0.1)
          btk(i,ks)=min(btk(i,ks),1.0)
        enddo !ks
      endif  
    enddo !i
!$OMP END PARALLEL DO

    return
    endsubroutine btklogexp

!***************************************************************************   
    subroutine btklogrouse
! Calculates the total load correction factor
! based on a logarithmic current velocity and 
! Rouse sediment concentration vertical profiles.
! written by Alex Sanchez, USACE-CHL
!***************************************************************************    
    use size_def
    use flow_def, only: vis,h,uv,iwet,rhow,grav
    use fric_def, only: z0,coefman,cfrict,uelwc,bsxy,bsvel
    use sed_def, only: nsed,btk,wsfall,taucr,diam,rsk,d50,dstar,&
        CtstarP,pbk,epsvk,bsk,ubk,specgrav
    use sed_lib, only: bedvel_vanrijn_wu
    use beta_lib, only: bslogrouse_table
    use fric_lib, only: fric_normapprough
    use prec_def
    implicit none
    integer :: i,ks
    real(ikind) :: zap,r,taup,cnp,Qc,bbk
   
!$OMP PARALLEL DO PRIVATE(i,ks,zap,r,taup,cnp,Qc,bbk)   
    do i=1,ncells
      if(iwet(i)==0)then
        btk(i,:)=1.0
        bsk(i,:)=1.0
        ubk(i,:)=0.0
      else
        zap=fric_normapprough(uelwc(i),uv(i),cfrict(i)) !Normalized apparent roughness
        !zap=z0(i)/h(i)  !Only for no waves
        cnp=0.05*d50(i)**0.16666667            !Manning corresponding to grain roughness 0.05=1/20
        if (coefman(i) .eq. 0) then
          taup=bsxy(i)*min(cnp/1.0e-6,1.0)**1.5      !If coefman(i) is zero, use very small number instead   MEB 01/12/2022
        else
          taup=bsxy(i)*min(cnp/coefman(i),1.0)**1.5  !Grain shear stress, min is used in case cn is less than cnp
        endif          
        Qc=uv(i)*h(i)
        do ks=1,nsed
          r=wsfall(ks)/(0.4*bsvel(i))  !Rouse number
          bsk(i,ks)=bslogrouse_table(zap,r)   !Suspended load correction factor, table    
          !!qbedk=pbk(i,ks,1)*(1.0-rsk(i,ks))*CtstarP(i,ks)*Qc   !Bed load transport, kg/m/s
          !!ubk(i,ks)=bedvel(taup,taucr(ks),diam(ks),dstar(ks),qbedk)   !Bed load velocity, m/s
          ubk(i,ks)=bedvel_vanrijn_wu(grav,taup,taucr(ks),specgrav,diam(ks))   !Bed load velocity, m/s
          ubk(i,ks)=max(min(ubk(i,ks),uv(i)),1e-20) !physical restriction, bed load cannot move faster than current velocity
          bbk=min(ubk(i,ks)/max(uv(i),1.0e-6),1.0)      !physical restriction, bed load cannot move faster than current velocity
          bbk=max(bbk,0.2)    
          btk(i,ks)=1.0/(rsk(i,ks)/bsk(i,ks)+(1.0-rsk(i,ks))/bbk) !Total load correction factor
          btk(i,ks)=max(btk(i,ks),0.1)
          btk(i,ks)=min(btk(i,ks),1.0)
        enddo !ks
      endif  
    enddo !i
!$OMP END PARALLEL DO

    return
    endsubroutine btklogrouse
