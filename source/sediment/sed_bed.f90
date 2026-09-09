!===============================================================================
! CMS Bed Routines
!
! Contains the following:
!   Bed change, sorting and gradation
!     bedchange     - Calculates the bed change terms for single-size sediment transport
!     bedchangesort - Calculates the total and fractional bed change, mixing layer 
!                     thickness and sorting for nonuniform sediment transport
!     bedgrad       - Calculates bed material gradation for sub-mixing layers
!     mixing_layer_min - Estimates the mixing layer thichness based 
!                        on the ripple dimensions and grain size
!   Bed Composition
!     d50sigma2diamlim - Calculates the diameter limits based on the d50
!                        dataset and geometric standard deviation
!     bed_d50sigma - Calculates sediment fractions for whole grid based on
!                    d50 dataset and geometric standard deviation sigma
!                    assuming a log-normal distribution.
!     bed_d16d50d84 - Calculates the bed material composition based on the 
!                     d16, d50, d84 datasets assuming a log-normal distribution.
!     bed_d35d50d90 - Calculates the bed material composition based on the 
!                     d35, d50, d90 datasets assuming a log-normal distribution.
!     sedpercentile - Calculates sediment grain size percentiles 
!                     from a fractional distribution  
!
!  Other
!    dper_read_error_msg
!
! written by Alex Sanchez, USACE-CHL
!===============================================================================

!**************************************************************************
    subroutine bedchange
! Calculates the bed change terms for single-size sediment transport
! written by Alex Sanchez, USACE-ERDC-CHL  
! Revised by Weiming Wu, Clarkson University
!**************************************************************************
    use size_def, only: ncells,ncellsD
    use comvarbl, only: dtime, timehrs
    use sed_def, only: scalemorph,solid,rhosed,alphat,&
        wsfall,Ctk,Ctkstar,Sb,dzb,cohesivesed,Etkstar,wsfallcohsed,alphacoh,rhobed
     use prec_def
     implicit none
     integer :: i
     real(ikind) :: fac
    
    if(cohesivesed)then
!$OMP PARALLEL DO PRIVATE(i)
      do i=1,ncells
        dzb(i)=scalemorph*dtime/rhobed(i,1)   &
              *((alphacoh(i)*wsfallcohsed(i)*Ctk(i,1)-Etkstar(i,1))+Sb(i,1))   !Note idry is already included in alphat
      enddo  
!$OMP END PARALLEL DO
    else
      fac = scalemorph*dtime/solid/rhosed     !for mass balance, scalemorph=1.0     
!$OMP PARALLEL DO PRIVATE(i)
      do i=1,ncells
        dzb(i)=fac*(alphat(i)*wsfall(1)*(Ctk(i,1)-Ctkstar(i,1))+Sb(i,1))   !Note iwet is already included in alphat
      enddo  
!$OMP END PARALLEL DO
    endif
    
    return
    end subroutine bedchange
    
!**************************************************************************
    subroutine check_hardbottom
! Checks the hardbottom and corrects the concentration potential 
! if the hardbottom depth is exceeded
!
! written by Alex Sanchez, USACE-ERDC-CHL  
!**************************************************************************
    use sed_def, only: nhard,idhard,hardzb,CtstarP,zb1,dzb,EtstarP
    use geo_def, only: zb
    implicit none
    integer :: ih,i
   
    do ih=1,nhard                 
       i=idhard(ih)
       if((zb1(i)+dzb(i))<hardzb(i))then
         CtstarP(i,:) = 0.7*CtstarP(i,:)
         EtstarP(i,:) = 0.7*EtstarP(i,:)   !Wu
        !Move hardbottom, small mass balance error
!         zb(i)=hardzb(i)
!         dzi=zb(i)-zb1(i)
!         do ks=1,nsed
!           dzbk(i,ks) = dzi*dzbk(i,ks)/dzb(i)
!         enddo
!         dzb(i)=dzi         
       endif    
    enddo !ih
   
    return
    end subroutine check_hardbottom
    
!*************************************************************************************
    subroutine bedchangesort()
! Calculates the total and fractional bed change, mixing layer thickness and sorting
! for nonuniform sediment transport
! written by Alex Sanchez, USACE-ERDC-CHL  
!*************************************************************************************    
#include "CMS_cpp.h"
    use size_def, only: ncells
    use flow_def, only: iwet
    use comvarbl, only: dtime
    use diag_def
    use diag_lib
    use sed_def
    use prec_def
    implicit none
    
    integer :: i
    real(ikind) :: fterm(nsed),fterm2(nsed),fterm3(nsed)
    real(ikind) :: znum(nsed),zdem(nsed),Ctkp(nsed),aws(nsed)
    real(ikind) :: pbksum,pbktemp(nsed),rxb,rxm,dbdz,fmdt
#ifdef DIAG_MODE
    logical :: isnankind
#endif

    !Relaxation factors
     rxb = 0.5
    rxm = 1.0-rxb 

    fmdt = scalemorph*dtime !Apply morphologic scaling factor here
    
!$OMP PARALLEL DO PRIVATE(i,dbdz,Ctkp,aws,fterm,fterm2,fterm3,znum,zdem,pbktemp,pbksum) REDUCTION(MAX:errpbk)   
    do i=1,ncells              
!!       if(iwet(i)==0 .or. zb(i)<=hardzb(i)+1.0e-6) then
       if(iwet(i)==0)then
         dzb(i) = 0.0
        dzbk(i,:) = 0.0
         cycle               
       endif
      !Bed change
      dbdz = dzb(i) - db(i,1) + db1(i,1) !Rise or fall of lower bound of mixing layer
      if(dbdz>=0.0)then        
        pbkstar(i,:) = pbk1(i,:) !from last time step
      else
        pbkstar(i,:) = pbk(i,:,2) !note: pbk(i,ks,2) does not change until bedgrad subroutine   
      endif
      !Ctkp = Ctkstar(i,:)/max(pbk(i,:,1),1.0e-20) !kg/m^3    !HLI 01/13/2017
      Ctkp = CtstarP(i,:)                                       !added by Wu, 2026-3-18
      aws = alphat(i)*wsfall(:)
      fterm = rhosed*solid*db(i,1) + fmdt*aws*Ctkp !kg/m^2
      fterm2 = db1(i,1)*pbk1(i,:) + (db(i,1)-db1(i,1))*pbkstar(i,:)
      fterm2 = Ctk(i,:)*db(i,1) - Ctkp*fterm2
      fterm2 = fmdt*(aws*fterm2 + Sb(i,:)*db(i,1)) !kg/m^2 **************
      fterm3 = fmdt*aws*Ctkp*pbkstar(i,:) !kg/m^2
      znum = fterm2/fterm 
      zdem = fterm3/fterm 
      dzb(i) = sum(znum)/(1.0-sum(zdem)) !Total bed change      
      dzbk(i,:) = (fterm2 + dzb(i)*fterm3)/fterm !Fractional bed change
#ifdef DIAG_MODE
      if(abs(dzb(i))>1.0e-6)then
        !dzbk(i,:) = dzbk(i,:)*dzb(i)/sum(dzbk(i,:))
        if(abs(sum(dzbk(i,:))-dzb(i))/dzb(i)>0.1 .and. itersed>10)then
          write(msg2,*) 'dzb(i) = ',dzb(i)
          write(msg3,*) 'dzbk(i,:) = ',dzbk(i,:)
          write(msg4,*) 'sum(znum) = ',sum(znum)
          write(msg5,*) 'sum(zdem) = ',sum(zdem)
          call diag_print_error('Problem calculating bed change',msg2,msg3,msg4,msg5)
          call diag_print_var(i,4)
        endif
      endif
#endif
      !Sorting
      pbktemp = pbk(i,:,1) !Save old values for convergence test
      dbdz = dzb(i) - db(i,1) + db1(i,1) !Update rise or fall of lower bound of mixing laye, note dzb from previous iteration    
      pbk(i,:,1) = (dzbk(i,:) + db1(i,1)*pbk1(i,:) - dbdz*pbkstar(i,:))/db(i,1)            
      pbk(i,:,1) = pbk(i,:,1)*rxb + pbktemp*rxm !Relaxation      
      !where(pbk(i,:,1)<1.0e-20) pbk(i,:,1)=1.0e-20     !Opened by Wu, 2026-3-18 
      pbksum = sum(max(pbk(i,:,1),1.0e-20)) !sum fraction    !HLI 01/13/2017
#ifdef DIAG_MODE
      if(abs(pbksum-1.0)>0.1 .or. isnankind(pbksum))then
        call diag_print_error('Problem calculating bed composition')
        call diag_print_var(i,4)
      endif
#endif
      pbk(i,:,1) = pbk(i,:,1)/pbksum !normalize
      errpbk = max(errpbk,maxval(abs(pbk(i,:,1)-pbktemp))) !check error
      !!if(errpbk>0.5)then
      !!  continue
      !!endif
    enddo !i 
!$OMP END PARALLEL DO

    return
    end subroutine bedchangesort
    
!**************************************************************************
    subroutine bedchangesort_cohesivesed
    ! Calculates the total and fractional bed change, mixing layer thickness and sorting
    ! for nonuniform sediment transport
    ! written by Weiming Wu, Clarkson University
    !**************************************************************************    
    use size_def, only: ncells
    use flow_def, only: iwet
    use comvarbl, only: dtime
    use diag_def
    use diag_lib
    use sed_def
    use prec_def
    implicit none
    integer :: i,ih,ks,js,itg,im,jstar,iwr,irm
    real(ikind) :: fac,dbdm,dzi,zbi,val
    real(ikind) :: fterm(nsed),fterm2(nsed),fterm3(nsed)
    !real(ikind) :: znum(nsed),zdem(nsed),Ctkp(nsed),aws(nsed)   !Closed by Wu
    real(ikind) :: znum(nsed),zdem(nsed),Etkp(nsed),aws(nsed)     !Wu
    real(ikind) :: pbksum,pbktemp(nsed),rxb,rxm,fmdt
 	
    rxb = 0.5
    rxm = 1.0-rxb 

    fmdt = scalemorph*dtime !Apply morphologic scaling factor here
    
!$OMP PARALLEL DO PRIVATE(i,dbdm,Etkp,aws,fterm,fterm2,fterm3,znum,zdem,pbktemp,pbksum) REDUCTION(MAX:errpbk)   
    do i=1,ncells  	  	  
       !!   if(idry(i).eq.0 .or. zb(i).le.hardzb(i)+1.0e-6) then
       if(iwet(i).eq.0)then
          dm(i) = 0.0; dmk(i,:) = 0.0
          dzb(i) = 0.0; dzbk(i,:) = 0.0 
          cycle    	   	
       endif
 	  
       !...Mixing Layer 
       dbms(i,1) = db(i,1)*rhobed(i,1)
       
       !...Bed change 	     	  	      
       dbdm = dm(i) - dbms(i,1) + dbms1(i,1) !Rise or fall of lower bound of mixing layer       
       if(dbdm.ge.0.0)then        
          pbkstar(i,:) = pbk1(i,:) !from last time step
       else
          pbkstar(i,:) = pbk(i,:,2) !note: pbk(i,ks,2) does not change until bedgrad subroutine   
       endif
       !Etkp=Etkstar(i,:)/max(pbk(i,:,1),1.0e-20) !kg/m^3
       Etkp = EtstarP(i,:)
       aws = alphat(i)*wsfall(:)
       if(cohesivesed) aws(1) = alphacoh(i)*wsfallcohsed(i)      
       fterm = dbms(i,1) + fmdt*Etkp !kg/m^2
       fterm2 = dbms1(i,1)*pbk1(i,:) + (dbms(i,1)-dbms1(i,1))*pbkstar(i,:)
       fterm2 = aws*Ctk(i,:)*dbms(i,1) - Etkp*fterm2
       fterm2 = fmdt*(fterm2 + Sb(i,:)*dbms(i,1)) !kg/m^2 **************
       fterm3 = fmdt*Etkp*pbkstar(i,:) !kg/m^2
       znum = fterm2/fterm         
       zdem = fterm3/fterm 	         
       dm(i) = sum(znum)/(1.0-sum(zdem)) !Total bed change      
       dmk(i,:) = (fterm2 + dm(i)*fterm3)/fterm !Fractional bed change      
 
       !...Sorting      
       pbktemp = pbk(i,:,1) !Save old values for convergence test 
       dbdm = dm(i) - dbms(i,1) + dbms1(i,1)
       pbk(i,:,1) = (dmk(i,:) + dbms1(i,1)*pbk1(i,:) - dbdm*pbkstar(i,:))/dbms(i,1)            
       pbk(i,:,1) = pbk(i,:,1)*rxb + pbktemp*rxm !Relaxation      
       !where(pbk(i,:,1).lt.1.0e-20) pbk(i,:,1)=1.0e-20 
       pbksum=sum(max(pbk(i,:,1),1.0e-20)) !sum fraction         
       pbk(i,:,1)=max(pbk(i,:,1),1.0e-20)/pbksum     
       !pbk(i,:,1)=pbk(i,:,1)/pbksum     
       errpbk=max(errpbk,maxval(abs(pbk(i,:,1)-pbktemp))) !check error
       
       dzb(i) = dm(i)/rhobed(i,1)   
       dzbk(i,:) = dmk(i,:)/rhobed(i,1) 
     
    enddo !i 
!$OMP END PARALLEL DO

    return
    end subroutine bedchangesort_cohesivesed
    
!**************************************************************************
    subroutine bedgrad()
! Calculates bed material gradation for sub-mixing layers
! written by Alex Sanchez, USACE-ERDC-CHL  
!**************************************************************************
    use size_def   
    use flow_def
    use const_def, only: small
    use sed_def, only: db, db1, dzb, pbk, pbkstar, dbmax, dbmin, nlay
    use prec_def
    implicit none
    integer :: i
    real(ikind) :: dbdz

!$OMP PARALLEL    
!$OMP DO PRIVATE(i,dbdz)                        
    do i=1,ncells
      if(iwet(i)==0)then
        db(i,:) = db1(i,:)
        cycle !db(i,js), and pbk(i,ks,js=2,nlay) stay the same
      endif
      dbdz = dzb(i) - db(i,1) + db1(i,1) !Rise or fall of lower bound of mixing layer
      db(i,2) = db1(i,2) + dbdz  !Calculate second layer thickness  **************
      pbk(i,:,2) = (db1(i,2)*pbk(i,:,2) + dbdz*pbkstar(i,:))/db(i,2)
      if(db(i,2)<=dbmax .and. db(i,2)>=dbmin)then !Second layer ok
        db(i,3:nlay) = db1(i,3:nlay) 
        !pbk(i,:,3:nlay) = pbk(i,:,3:nlay) stays the same
      elseif(db(i,2)>dbmax)then !Second layer too thick, split into two              
        !First merge last two layers   
        db(i,nlay) = db1(i,nlay) + db1(i,nlay-1) 
        pbk(i,:,nlay) = (db1(i,nlay)*pbk(i,:,nlay) + & 
          db1(i,nlay-1)*pbk(i,:,nlay-1))/db(i,nlay)
        !Move index of layers 4 to nlay-1 
        db(i,4:nlay-1) = db1(i,3:nlay-2)
        pbk(i,:,4:nlay-1) = pbk(i,:,3:nlay-2)
        !Second and third layer compositions equal
        pbk(i,:,3) = pbk(i,:,2)
        if(dbdz>dbmin)then !dbdz is large enough make it the second layer
          db(i,2) = dbdz
          db(i,3) = db1(i,2)
        else !split second layer into equal parts
          db(i,2) = 0.5*db(i,2)   !Note that db(i,2) = db1(i,2) + dbdz was set above
          db(i,3) = db(i,2)
        endif                                         
      else !Second layer too thin
        !Merge second with third layer
        pbk(i,:,2) = (db(i,2)*pbk(i,:,2) + &
          db1(i,3)*pbk(i,:,3))/(db(i,2) + db1(i,3))
        db(i,2) = db(i,2) + db1(i,3)
        
        !Method 2
        if(db1(i,nlay)>=dbmax)then !split last layer
          db(i,3:nlay-2) = db1(i,4:nlay-1)
          pbk(i,:,3:nlay-2) = pbk(i,:,4:nlay-1)
          db(i,nlay-1) = 0.5*db1(i,nlay)
          db(i,nlay) = db(i,nlay-1)
          pbk(i,:,nlay-1) = pbk(i,:,nlay)
          !pbk(i,:,nlay) = pbk(i,:,nlay) stays the same
        else !add bottom layer
          db(i,3:nlay-1) = db1(i,4:nlay)
          pbk(i,:,3:nlay-1) = pbk(i,:,4:nlay)
          db(i,nlay) = db1(i,nlay)
          pbk(i,:,nlay) = pbk(i,:,nlay)
        endif
                
        !!--- Method 3 ------
        !pbk(i,:,3:nlay-1) = pbk(i,:,4:nlay)
        !db(i,3:nlay-1) = db1(i,4:nlay)
        !!Add layer at bottom
        !pbk(i,:,nlay) = pbk(i,:,nlay)
        !db(i,nlay) = db1(i,nlay)
        
      endif
    enddo !i
!$OMP END DO
!$OMP END PARALLEL
    
    !Check total depth of layers
!    if(hardbottom)then
!      do ih=1,nhard
!        i=idhard(ih)
!        if(iwet(i)==0) cycle
!        zbotlay=zb(i)-db(i,1) !Bottom elevation of mixing layer
!        do js=2,nlay-1
!          zbotlay=zbotlay-db(i,js)
!          if(zbotlay<=hardzb(i))then
!            db(i,js)=db(i,js)-hardzb(i)+zbotlay     
!            db(i,js)=max(small,db(i,js))
!          endif
!        enddo !js
!      enddo !ih
!    endif
    
    return
    end subroutine bedgrad
    
!**************************************************************************
    subroutine bedgrad_consolid()
    ! Calculates bed material gradation for sub-mixing layers
    ! written by  Weiming Wu, Clarkson University
    !**************************************************************************
    use size_def   
    use flow_def
    use const_def, only: small, great
    use comvarbl, only: dtime
    use sed_def, only: db, db1, dzb, pbk, pbk1, pbkstar, dbmax, dbmin, nlay, scalemorph
    use sed_def, only: dm,dmk,dbms,dbms1,rhobed,rhobedcoh,tconsolid,tconsolid1  !,tconsolid0(:)
    use sed_def, only: arhobed,prhobed,rhobedcoh0,tdexp,rhobedcohtdexp,rhobedcoh1yr,betarhobed,rhobednoncoh
    use prec_def
    implicit none
    integer :: i,j,klay
    real(ikind) :: dbdm,tconsolidstar,rhobedcohold,fmdt
   
    fmdt = scalemorph*dtime !Apply morphologic scaling factor here
    
!$OMP PARALLEL    
!$OMP DO PRIVATE(i,dbdm,tconsolidstar)                        
    do i=1,ncells
      if(iwet(i)==0)then
        dbms(i,:)=dbms1(i,:)
        tconsolid(i,:)=tconsolid1(i,:)+fmdt
        cycle !db(i,js), and pbk(i,ks,js=2,nlay) stay the same
      endif
      dbdm = dm(i) - dbms(i,1) + dbms1(i,1) !Rise or fall of lower bound of mixing layer     
      if(dbdm.ge.0.0)then              
        tconsolidstar = tconsolid1(i,1)
      else
        tconsolidstar = tconsolid1(i,2)   
      endif
      dbms(i,2) = dbms1(i,2) + dbdm  !Calculate second layer thickness  **************  
      if(dbms(i,2).le.dbmax*rhobed(i,2) .and. dbms(i,2).ge.dbmin*rhobed(i,2))then !Second layer ok                              
        if(dbdm.ge.0.0)then        
          tconsolid(i,2) = (dbms1(i,2)*pbk(i,1,2)*tconsolid1(i,2)+dbdm*pbkstar(i,1)*tconsolidstar)  &
                          /(dbms1(i,2)*pbk(i,1,2)+dbdm*pbkstar(i,1)+small)
        endif
        tconsolid(i,2) = tconsolid1(i,2) + fmdt
        pbk(i,:,2) = (dbms1(i,2)*pbk(i,:,2)+dbdm*pbkstar(i,:))/dbms(i,2)
        dbms(i,3:nlay) = dbms1(i,3:nlay)     
        tconsolid(i,3:nlay) = tconsolid1(i,3:nlay) + fmdt
      elseif(dbms(i,2).gt.dbmax*rhobed(i,2))then !Second layer too thick, split into two              
          !First merge last two layers   
        dbms(i,nlay) = dbms1(i,nlay) + dbms1(i,nlay-1)     
        tconsolid(i,nlay) = (dbms1(i,nlay)*pbk(i,1,nlay)*tconsolid1(i,nlay)   &
                            +dbms1(i,nlay-1)*pbk(i,1,nlay-1)*tconsolid1(i,nlay-1))  &
                          /(dbms1(i,nlay)*pbk(i,1,nlay)+dbms1(i,nlay-1)*pbk(i,1,nlay-1)+small)  &
                          +fmdt
        pbk(i,:,nlay) = (dbms1(i,nlay)*pbk(i,:,nlay)+dbms1(i,nlay-1)*pbk(i,:,nlay-1))/dbms(i,nlay)
          !Move index of layers 4 to nlay-1        
        dbms(i,4:nlay-1) = dbms1(i,3:nlay-2)
        pbk(i,:,4:nlay-1) = pbk(i,:,3:nlay-2)     
        tconsolid(i,4:nlay-1) = tconsolid1(i,3:nlay-2) + fmdt
           !Calculate second and third layers
        if(dbdm.gt.dbmin*rhobed(i,2))then !dbdz is large enough make it the second layer
          dbms(i,2) = dbdm
          dbms(i,3) = dbms1(i,2)
          pbk(i,:,3) = pbk(i,:,2)
          !!pbk(i,:,2)=pbk1(i,:)
          pbk(i,:,2) = pbkstar(i,:)
          tconsolid(i,3) = tconsolid1(i,2) + fmdt
          !!tconsolid(i,2)=tconsolid1(i,1)+fmdt
          tconsolid(i,2) = tconsolidstar + fmdt
        else !dbdz is small so second layer split into equal parts
          if(dbdm.ge.0.0)then        
            tconsolid(i,2) = (dbms1(i,2)*pbk(i,1,2)*tconsolid1(i,2)+dbdm*pbkstar(i,1)*tconsolidstar)  &
                          /(dbms1(i,2)*pbk(i,1,2)+dbdm*pbkstar(i,1)+small) + fmdt
          endif
          pbk(i,:,2) = (dbms1(i,2)*pbk(i,:,2)+dbdm*pbkstar(i,:))/dbms(i,2)
          dbms(i,3) = 0.5*dbms(i,2)   !Note that db(i,2)=db1(i,2)+dbdz was set above
          dbms(i,2) = dbms(i,3)
          tconsolid(i,3) = tconsolid1(i,2)
          pbk(i,:,3) = pbk(i,:,2)

          !!dbms(i,2)=0.5*dbms(i,2)   !Note that db(i,2)=db1(i,2)-dbdz was set above
          !!dbms(i,3)=dbms(i,2)
          !!!Note that for db(i,2).gt.dbmax, dbdz.le.0.0  must be true so pbkstar=pbk1(i,:)
          !!tconsolid(i,3)=tconsolid(i,2)+fmdt
          !!tconsolid(i,2)= ((dbms(i,2)-dbdm)*pbk(i,1,2)*tconsolid(i,2)+dbdm*pbk1(i,1)*tconsolid(i,1))  &
          !!                 /((dbms(i,2)-dbdm)*pbk(i,1,2)+dbdm*pbk1(i,1)+small) &
          !!               +fmdt
          !!pbk(i,:,3)=pbk(i,:,2)
          !!pbk(i,:,2)=((dbms(i,2)-dbdm)*pbk(i,:,2)+dbdm*pbk1(i,:))/dbms(i,2)   
        endif                                         
      else !Second layer too thin, merge with third layer
        if(dbdm.ge.0.0)then        
          tconsolid(i,2) = (dbms1(i,2)*pbk(i,1,2)*tconsolid1(i,2)+dbdm*pbkstar(i,1)*tconsolidstar)  &
                          /(dbms1(i,2)*pbk(i,1,2)+dbdm*pbkstar(i,1)+small)
        endif
        pbk(i,:,2) = (dbms1(i,2)*pbk(i,:,2)+dbdm*pbkstar(i,:))/dbms(i,2)
        tconsolid(i,2)= (dbms(i,2)*pbk(i,1,2)*tconsolid1(i,2)+dbms1(i,3)*pbk(i,1,3)*tconsolid1(i,3))  &
                         /(dbms(i,2)*pbk(i,1,2)+dbms1(i,3)*pbk(i,1,3)+small) &
                       +fmdt
        pbk(i,:,2)=(dbms(i,2)*pbk(i,:,2)+dbms1(i,3)*pbk(i,:,3))/(dbms(i,2)+dbms1(i,3))
        dbms(i,2)=dbms(i,2)+dbms1(i,3)  
        !!pbk(i,:,2)=((dbms1(i,2)-dbdz)*pbk(i,:,2)+dbms1(i,3)*pbk(i,:,3))/dbms(i,2)
        if(dbms1(i,nlay)>=dbmax*rhobed(i,nlay))then !split last layer
          dbms(i,3:nlay-2) = dbms1(i,4:nlay-1)
          pbk(i,:,3:nlay-2) = pbk(i,:,4:nlay-1)
          tconsolid(i,3:nlay-2)=tconsolid1(i,4:nlay-1)+fmdt
          dbms(i,nlay-1) = 0.5*dbms1(i,nlay)
          dbms(i,nlay) = dbms(i,nlay-1)
          pbk(i,:,nlay-1) = pbk(i,:,nlay)
          !pbk(i,:,nlay) = pbk(i,:,nlay) stays the same
          tconsolid(i,nlay-1)=tconsolid1(i,nlay)+fmdt
          tconsolid(i,nlay)=tconsolid1(i,nlay)+fmdt
        else !add bottom layer
          dbms(i,3:nlay-1) = dbms1(i,4:nlay)
          pbk(i,:,3:nlay-1) = pbk(i,:,4:nlay)
          tconsolid(i,3:nlay-1)=tconsolid1(i,4:nlay)+fmdt
          dbms(i,nlay) = dbms1(i,nlay)
          pbk(i,:,nlay) = pbk(i,:,nlay)
          tconsolid(i,nlay)=tconsolid1(i,nlay)+fmdt
          !tconsolid(i,nlay) = great 
        endif    
      endif
           
      if(dbdm.le.0.0) then
        if(dmk(i,1).ge.0.0) then
          tconsolid(i,1)= (dbms1(i,1)*pbk1(i,1)*tconsolid1(i,1)-dbdm*pbkstar(i,1)*tconsolidstar)  &
                         /(dbms(i,1)*pbk(i,1,1)-dbdm*pbkstar(i,1)+dmk(i,1)+small) + fmdt
        else
          tconsolid(i,1)= (max(0.0,dbms1(i,1)*pbk1(i,1)+dmk(i,1))*tconsolid1(i,1)  &
                           -dbdm*pbkstar(i,1)*tconsolidstar)  &
                         /(max(0.0,dbms(i,1)*pbk(i,1,1)+dmk(i,1))-dbdm*pbkstar(i,1)+small) + fmdt
        endif    
      else  
        if(dmk(i,1).ge.0.0) then
          tconsolid(i,1)= (dbms1(i,1)*pbk1(i,1)*tconsolid1(i,1))  &
                         /(dbms(i,1)*pbk(i,1,1)+dmk(i,1)+small) + fmdt
        else
          tconsolid(i,1)= tconsolid1(i,1) + fmdt
        endif    
      endif     
      !tconsolid(i,1)= (max(0.0,dbms1(i,1)*pbk1(i,1)+min(0.0,dmk(i,1)))*tconsolid(i,1)  &
      !                     - min(0.0,dbdm)*pbkstar(i,1)*tconsolidstar)  &
      !               /(max(0.0,dbms(i,1)*pbk(i,1,1)+min(0.0,dmk(i,1)))     &
      !                     - min(0.0,dbdm)*pbkstar(i,1) + max(0.0,dmk(i,1))+small) + fmdt     
    enddo !i
!$OMP END DO
!$OMP END PARALLEL     

    do i=1,ncells
    do j=2,nlay
      pbk(i,:,j) = max(pbk(i,:,j),small)/sum(max(pbk(i,:,j),small))  
    enddo
    enddo
    
    do i=1,ncells
      do klay=1,nlay
        if(tconsolid(i,klay).le.tdexp) then  !less than tdexp, use exponential form
          rhobedcoh(i,klay)=rhobedcoh0*(1.0-arhobed*exp(-prhobed*tconsolid(i,klay)))/(1.0-arhobed)
        elseif((tconsolid(i,klay).gt.tdexp).and.(tconsolid(i,klay).lt.31536000.0)) then  !less than 1 yr
          rhobedcoh(i,klay)=rhobedcohtdexp+(rhobedcoh1yr-rhobedcohtdexp)  & 
                                    *(tconsolid(i,klay)-tdexp)/(31536000.0-tdexp)        
        else
          rhobedcoh(i,klay)=rhobedcoh1yr+betarhobed*alog10(tconsolid(i,klay)/31536000.0)
        endif
        rhobed(i,klay)=1.0/(pbk(i,1,klay)/rhobedcoh(i,klay)+(1.0-pbk(i,1,klay))/rhobednoncoh)
      enddo
    enddo

    do i=1,ncells
      !dzb(i)=dm(i)/rhobed(i,1)   ! Near-bed exchange
      do klay=1,nlay
        if(tconsolid(i,klay).le.tdexp) then  !less than tdexp, use exponential form
          rhobedcohold=rhobedcoh0*(1.0-arhobed*exp(-prhobed*(tconsolid(i,klay)-fmdt)))/(1.0-arhobed)
        elseif((tconsolid(i,klay).gt.tdexp).and.(tconsolid(i,klay).lt.31536000.0)) then  !less than 1 yr
          rhobedcohold=rhobedcohtdexp+(rhobedcoh1yr-rhobedcohtdexp)  &
                                    *(tconsolid(i,klay)-fmdt-tdexp)/(31536000.0-tdexp)        
        else
          rhobedcohold=rhobedcoh1yr+betarhobed*alog10((tconsolid(i,klay)-fmdt)/31536000.0)
        endif
        dzb(i)=dzb(i)+dbms(i,klay)*pbk(i,1,klay)*(1.0/rhobedcoh(i,klay)-1.0/rhobedcohold)   !Consolidation
      enddo
    enddo

    return
    end subroutine bedgrad_consolid
    
!*******************************************************************************
    subroutine mixing_layer
! Estimate mixing layer thickness based on the ripple dimensions and grain size
! modified from Lund-CIRP subroutine written by Magnus Larson 
! written by Alex Sanchez, USACE-ERDC-CHL 
!*******************************************************************************       
    use size_def
    use sed_def, only: s1grav,d50,dzb,&
       dmconst,mixlayconst,db1min,dbmax,db
    use geo_def, only: zb
    use wave_flowgrid_def, only: Worb,Wper
    use cms_def, only: noptset
    use const_def, only: small,twopi
    !use comvarbl
    use prec_def
    implicit none
    integer :: i
    real(ikind) :: riph,pw,ripwh,Aw
    
    if(mixlayconst)then
      !$omp parallel do private(i)
      do i=1,ncells
        db(i,1)=min(max(dmconst,2.0*d50(i)),dbmax) !Note: dmconst>=db1min
      enddo  
      !$omp end parallel do
    else 
      !Note that D50 use here is from previous time step
      if(noptset>=3)then
        if(.not.allocated(Worb)) call wave_flgrid_init
        
        !$omp parallel do private(i,riph,pw,ripwh,Aw)  
        do i=1,ncells
          riph=142.8571429*d50(i)             !Current ripple height, 142.8571429=1000.0/7.0   
          pw=Worb(i)*Worb(i)/(s1grav*d50(i)+small)  !Wave mobility parameter
          if(pw<=10.0)then
            Aw=Worb(i)*Wper(i)/twopi          !amplitude of near-bed wave exursion
            ripwh=0.220*Aw                    !Wave ripple height
          elseif(pw<250.0-1.0e-6)then !to avoid underflow
            Aw=Worb(i)*Wper(i)/twopi   
            ripwh=2.8e-13*Aw*(250.0-pw)**5    !Wave ripple height
          else
            ripwh=0.0                         !Wave ripple height
          endif                
          riph=max(ripwh,riph)
          db(i,1)=min(max(0.5*riph,2.0*d50(i),db1min),dbmax)
        enddo
        !$omp end parallel do
      else
        !$omp parallel do private(i,riph)   
        do i=1,ncells
          riph=142.8571429*d50(i)            !Current ripple height, 142.8571429=1000.0/7.0
          db(i,1)=min(max(0.5*riph,2.0*d50(i),db1min),dbmax)
          !db(i,1)=min(max(71.43*d50(i),db1min),dbmax) !Note: 71.43=0.5*riph=142.8571429=1000.0/7.0/2.0
        enddo
        !$omp end parallel do
      endif  
    endif    
    
    return
    end subroutine mixing_layer    

!********************************************************
    subroutine d50sigma2diamlim(nd,d50lay,stdlay,dlim)
! Calculates the diameter limits based
! on the d50 dataset and geometric standard deviation
!********************************************************
    use size_def
    use math_lib, only: logninv
    use prec_def
    implicit none
    !Input/Output
    integer,intent(in) :: nd !Number of sediment size classes (diameters)    
    real(ikind),intent(in) :: d50lay(ncellsD) !Layer median grain size for all cells
    real(ikind),intent(in) :: stdlay      !Layer geometric standard deviation
    real(ikind),intent(out) :: dlim(nd+1) !Characteristic diameters and bounds
    !Internal variables
    integer :: i,ks
    real(ikind) :: d50low,d50high,mu,d1,d2
    real(ikind) :: logd1,val,sig,prob
    
!The upper and lower limits of the histogram
!are determined with the maximum and mininum values
!values of the d50 dataset 
    d50low = 1000.0
    d50high = -1000.0
    do i=1,ncells     
      d50low = min(d50low,d50lay(i))
      d50high = max(d50high,d50lay(i))
    enddo
    mu = log(d50low)
    sig = log(stdlay)
    prob = 0.001
    d1 = logninv(prob,mu,sig)
    d1 = max(d1,0.05) !Limit for non-cohesive sediments ******************
    mu = log(d50high)
    prob = 0.98
    d2 = logninv(prob,mu,sig)
    
!Size class limits using log-scale    
    logd1 = log(d1)
    val = (log(d2)-logd1)/float(nd-1)
    do ks=1,nd+1
      dlim(ks) = exp(logd1+(ks-1.5)*val)
    enddo            

    return
    end subroutine d50sigma2diamlim

!***************************************************************************
    subroutine bed_d50sigma(nd,d,dlim,d50lay,stdlay,pbklay)
! Calculates sediment fractions for whole grid based on
! the d50 dataset and geometric standard deviation sigma
! assuming a log-normal distribution.
!
! Input:
!  nd - # of sediment size classes
!  d - Sediment size class diameters [m]
!  dlim - Sediment size class limits [m]
!  d50lay - Bed layer 50th percentile diameters for a layer
!  stdlay - Geometric standard deviation for whole layer
!
! Output:
!   pbklay - Fractional composition for a bed layer
!
! Subprograms: lognpdf
!
! Revision History:
!   11/27/2012 Alex Sanchez
!      Changed to work with any bed layer instead of just the surface layer
!***************************************************************************
    use size_def, only: ncells,ncellsD
    use math_lib, only: lognpdf
    use prec_def
    implicit none
    !Input/Output
    integer,intent(in) :: nd  !Number of sediment size classes (diameters)
    real(ikind),intent(in) :: d(nd),dlim(nd+1) !Characteristic diameters and diameter bounds
    real(ikind),intent(in) :: d50lay(ncellsD)       !50th percentile diameter for a layer
    real(ikind),intent(in) :: stdlay                !Geometric standard deviation
    real(ikind),intent(out) :: pbklay(ncellsD,nd) !Fractional composition for a bed layer
    !Internal Variables
    integer :: i,k
    real(ikind) :: mu,sig  
    
    sig = log(stdlay)
    do i=1,ncells
      mu = log(d50lay(i))
      do k=1,nd
        pbklay(i,k) = (dlim(k+1)-dlim(k))*lognpdf(d(k),mu,sig)
!        pbklay(i,k) = max(pbklay(i,k),1.0e-20)              !HLI 01/13/2017
      enddo      
      pbklay(i,:) = pbklay(i,:)/sum(max(pbklay(i,:),1.0e-20)) !Check conservations     !HLI 01/13/2017
    enddo            

    return
    end subroutine bed_d50sigma

!**************************************************************************
    subroutine bed_d16d50d84(nd,d,dlim,d16lay,d50lay,d84lay,pbklay)
! Calculates the bed material composition based 
! on the d16, d50, d84 datasets assuming a log-normal distribution.
!
! Input:
!  d16lay,d50lay,d84lay - Bed layer percentile diameters
!
! Output:
!   pbklay - Fractional composition for a bed layer
!
! Subprograms: lognpdf
!
! Revision History:
!   11/27/2012 Alex Sanchez
!      Changed to work with any bed layer instead of just the surface layer
!**************************************************************************
    use size_def, only: ncells,ncellsD
    use math_lib, only: lognpdf
    use prec_def
    implicit none
    !Input/Output
    integer,intent(in) :: nd  !Number of sediment size classes (diameters)
    real(ikind),intent(in) :: d(nd),dlim(nd+1) !Characteristic diameters and diameter bounds
    real(ikind),intent(in) :: d16lay(ncellsD),d50lay(ncellsD),d84lay(ncellsD) !Bed layer percentile diameters
    real(ikind),intent(out) :: pbklay(ncellsD,nd)  !Fractional composition for a bed layer
    !Internal variables
    integer :: i,k
    real(ikind) :: mu,sig

    do i=1,ncells
      sig = log(sqrt(d84lay(i)/d16lay(i)))
      mu = log(d50lay(i))
      do k=1,nd        
        pbklay(i,k) = (dlim(k+1)-dlim(k))*lognpdf(d(k),mu,sig)
!        pbklay(i,k) = max(pbklay(i,k),1.0e-20)              !HLI 01/13/2017
      enddo
      pbklay(i,:) = pbklay(i,:)/sum(max(pbklay(i,:),1.0e-20)) !Check conservation       !HLI 01/13/2017
    enddo

    return
    end subroutine bed_d16d50d84
    
!***************************************************************************  
    subroutine bed_d35d50d90(nd,d,dlim,d35lay,d50lay,d90lay,pbklay)
! Calculates the bed material composition based 
! on the d35, d50, d90 datasets assuming a log-normal distribution.
!
! Input:
!  d35lay,d50lay,d90lay - Bed layer percentile diameters
!
! Output:
!   pbklay(ncellsD,nsed) - Fractional composition for a bed layer
!
! Subprograms: lognpdf
!
! Revision History:
!   11/27/2012 Alex Sanchez
!      Changed to work with any bed layer instead of just the surface layer
!**************************************************************************
    use size_def, only: ncells,ncellsD
    use math_lib, only: lognpdf
    use prec_def
    implicit none
    !Input/Output
    integer,intent(in) :: nd  !Number of sediment size classes (diameters)
    real(ikind),intent(in) :: d(nd),dlim(nd+1) !Characteristic diameters and diameter bounds
    real(ikind),intent(in) :: d35lay(ncellsD),d50lay(ncellsD),d90lay(ncellsD) !Bed layer percentile diameters
    real(ikind),intent(out) :: pbklay(ncellsD,nd) !Fractional composition for a bed layer
    !Internal variables
    integer :: i,k
    real(ikind) :: mu,sig
       
    do i=1,ncells
      sig = log((d90lay(i)/d35lay(i))**0.61)   
      mu = log(d50lay(i))      
      do k=1,nd
        pbklay(i,k) = (dlim(k+1)-dlim(k))*lognpdf(d(k),mu,sig)        
!        pbklay(i,k) = max(pbklay(i,k),1.0e-20)              !HLI 01/13/2017
      enddo
      pbklay(i,:) = pbklay(i,:)/sum(max(pbklay(i,:),1.0e-20)) !Check conservation   !HLI 01/13/2017   
    enddo

    return
    end subroutine bed_d35d50d90
    
!************************************************************************
    subroutine bed_perdiam(j)
! Calculates the bed layer j composition pbk(1:ncells,1:nsed,j)
! using percentile diameter datasets
!
! written by Alex Sanchez, USACE-CHL
!************************************************************************
#include "CMS_cpp.h"
    use size_def, only: ncellsD,ncells
    use sed_def
    use diag_def
    use diag_lib
#ifdef XMDF_IO 
    use in_xmdf_lib, only: readscalh5
#endif   
    use in_lib, only: readscalTxt
    use prec_def
    use geo_def, only: mapid                 !MEB  01/19/2022
    
    implicit none
    
    !Input/Output
    integer,intent(in) :: j  !Bed layer id
    
    !Internal variables
    integer :: i,ii,ks,ierr,istart
    integer :: idpr(nperdiam)
    real(ikind) :: plim(nsed+1),sumpbk,fac
    real(ikind) :: per(nperdiam),dper(nperdiam)
    real(ikind) :: dtemp(ncellsD)
    
    integer :: nfrac, frac_array(ncellsD)    !MEB  01/19/2022
    character(len=200) :: file,path
    character(len=10)  :: aext
    
    interface
      function Int2Str (k)
        integer,intent(in) :: k
        character(len=20) :: Int2Str
      end function
    end interface    
    
    !Read percentile diameter datasets
    !and determine number of input datasets
    nperinp = 0    
    do i=1,nperdiam
      if(bedlay(j)%perdiam(i)%inp)then
        file = bedlay(j)%perdiam(i)%file  
        path = bedlay(j)%perdiam(i)%path

        call fileext(trim(file),aext)      
        select case (aext)
        case('h5')
#ifdef XMDF_IO     
        call readscalh5(file,path,dtemp,ierr) 
#endif       
        case('txt')
          call readscalTxt(file,dtemp,ierr)
        end select
        
        if(ierr<0)then
          bedlay(j)%perdiam(i)%inp = .false.    
          call dper_read_error_msg(file,path)
        else
          allocate(bedlay(j)%perdiam(i)%dper(ncellsD))
          bedlay(j)%perdiam(i)%dper = dtemp !Bed layer percentile diameter [mm]
          nperinp = nperinp + 1
          idpr(nperinp) = i
        endif
      endif
    enddo
    
    !Calculate fractional bed composition
    per=0; dper=0.0 !Initialize
    nfrac = 0; frac_array = 0
    do i=1,ncells
      !Collect percentiles and diameters
      do ii=1,nperinp
        per(ii) = real(iper(idpr(ii)),kind=ikind)
        dper(ii) = bedlay(j)%perdiam(idpr(ii))%dper(i)
      enddo
      !Linearly interpolate percentiles at dlim  
      plim = 0.0
      istart = 1
      do ks=1,nsed+1        
        if(diamlim(ks)<dper(1))then !Below range
          plim(ks) = 0.0
          cycle
        endif
        if(diamlim(ks)>dper(nperinp))then !Above range
          plim(ks:nsed+1) = 100.0
          exit
        endif
        do ii=istart,nperinp-1
          if(diamlim(ks)>dper(ii) .and. diamlim(ks)<=dper(ii+1))then
            fac = (diamlim(ks)-dper(ii))/(dper(ii+1)-dper(ii))
            fac = max(min(fac,1.0),0.0)
            plim(ks) = (1.0-fac)*per(ii) + fac*per(ii+1)
            istart = ii
            exit
          endif
        enddo
      enddo
      !Differenciate percentiles to get fractions
      do ks=1,nsed
        pbk(i,ks,j) = (plim(ks+1)-plim(ks))/100.0
!        pbk(i,ks,1) = max(pbk(i,ks,1),1.0e-10) !To avoid divide by zero  !HLI 01/13/2017
      enddo
      sumpbk=sum(max(pbk(i,:,j),1.0e-20))                                  !HLI 01/13/2017
      if(sumpbk<0.8)then
        nfrac = nfrac + 1
        frac_array(nfrac) = i
      endif
      pbk(i,:,j) = pbk(i,:,j)/sumpbk
    enddo

    !Added this section to consolidate all of the messages for warnings of this type     MEB  01/19/2022
99  format(15(i0,x))
100 format('Sum of fractions is less than 0.8 for ',i0,' cells.')   
101 format('Cell IDs written to file: "fraction_warning.txt"')
    if(nfrac .gt. 0) then
      write(msg2,100) nfrac
      write(msg3,101) 
      call diag_print_warning(msg2,msg3,'')
      open(200,file='fraction_warning.txt',status='unknown') 
      write(200,99) (mapid(frac_array(i)),i=1,nfrac)  !changed to the ID of the cell as in SMS.  
      close(200)
    endif
    
    return
    end subroutine bed_perdiam    
   
!**************************************************************************
    subroutine sedpercentile(iper,dper)
! Calculates sediment grain size percentiles from a fractional distribution    
! written by  Weiming Wu, NCCHE; Alex Sanchez, USACE-ERDC-CHL  
!**************************************************************************
    use size_def
    use sed_def, only: nsed,pbk,diam,diamlim,logdiamlim
    use prec_def
    implicit none
    !Input/Output
    integer,intent(in) :: iper
    real(ikind),intent(out) :: dper(ncellsD)
    !Internal variables
    integer :: i,ks,kks,ks0,ks2
    real(ikind) :: pbcum(nsed+1),fac,per  
      
    per = real(iper,kind=ikind)/100.0
    
!$OMP PARALLEL DO PRIVATE(i,ks,kks,ks0,ks2,pbcum,fac)           
    do i=1,ncells
      !call sed_dper(nsed,diam,diamlim,logdiamlim,pbk(i,:,1),per,dper(i))    
      dper(i) = diamlim(1) 
      pbcum(1) = 0.0
      do ks=1,nsed
        pbcum(ks+1) = pbcum(ks) + pbk(i,ks,1)
      enddo
      if(pbcum(nsed+1)<=per)then
        dper(i) = diamlim(nsed+1)
        cycle
      endif
      do ks=1,nsed
        if(pbcum(ks+1)>=per .and. pbcum(ks)<per)then
          !fac = (pbcum(ks+1)-per)/(pbcum(ks+1)-pbcum(ks))
          !dper(i) = exp((1.0-fac)*logdiamlim(ks+1)+fac*logdiamlim(ks))
          !Transverse duplicate diameters
          do kks=ks,1,-1
            if(abs(diam(kks)-diam(ks))<1.0e-6)then
              ks0 = kks
            endif
          enddo
          do kks=ks,nsed
            if(abs(diam(kks)-diam(ks))<1.0e-6)then
              ks2 = kks+1
            endif
          enddo
          exit
        endif
      enddo !ks
      fac = (pbcum(ks2)-per)/(pbcum(ks2)-pbcum(ks0))
      dper(i) = exp((1.0-fac)*logdiamlim(ks2)+fac*logdiamlim(ks0))
    enddo
!$OMP END PARALLEL DO

    return
    end subroutine sedpercentile
    
!*******************************************************************************
    subroutine dper_read_error_msg(afile,apath)
!*******************************************************************************
    use diag_def, only: dgfile,dgunit, msg, msg2, msg3, msg4
    use diag_lib, only: diag_print_error
    implicit none
    
    integer :: i,iunit(2)
    character(len=200) :: afile,apath    
    
    write(msg,*)  'Could not find percentile dataset '
    write(msg2,*) '       File: ',trim(afile)
    write(msg3,*) '       Path: ',trim(apath)
    write(msg4,*) 'Check input files and restart'
    call diag_print_error(msg, msg2, msg3, msg4)
    
    end subroutine dper_read_error_msg
    
!*******************************************************************************
    subroutine getFormatType(filename, newFormat, ierr)
!*******************************************************************************
! This routine looks at the 10th character of the first dataset returned to see
! if it is an 'S' which matches the new format: 'Fraction_S1_L1' instead of the
! old: 'Fraction_01 (1)' and returns True or False
!*******************************************************************************
    use XMDF
    use const_def, only: READONLY
    
    character(len=*), intent(in) :: filename
    logical, intent(out)         :: newFormat
    integer, intent(out)         :: ierr
    integer(XID) :: fid, gid, nfid, ndid
    integer      :: nSDsets, sMaxPathLength
    character(len=200),allocatable :: paths(:)
    character(len=200) :: sPath(1)
    character(len=1) :: aChar
    
    newFormat = .false.
    call XF_OPEN_FILE(trim(filename), READONLY, fid, ierr)
    call XF_OPEN_GROUP(fid, 'Datasets', gid, ierr)
    call XF_GET_SCALAR_DATASETS_INFO (gid, nSDsets, sMaxPathLength, ierr)
    allocate(paths(nSDsets))
    paths = ' '; sPath = ' '; aChar = ' '
    call XF_GET_SCALAR_DATASET_PATHS (gid, nSDsets, sMaxPathLength, sPath, ierr)
    if (ierr .le. 0) return
    
    aChar = sPath(1)(10:10)
    if(aChar == 'S') newFormat = .true.
 
    deallocate(paths)
    
    return
    end subroutine getFormatType
    
!*******************************************************************************
    subroutine read_pbk(j,pbkfilelay,pbkpathlay)
! Reads the bed material composision file
! written by Alex Sanchez, USACE-ERDC-CHL 
!******************************************************************************* 
#include "CMS_cpp.h"
    use size_def
    use sed_def
#ifdef XMDF_IO     
    use in_xmdf_lib, only: readscalh5
#endif   
    use in_lib, only: readscalTxt
    implicit none
    
    integer :: i,j,ks,ierr
    character(len=200) :: apath, altpath
    character(len=5)   :: apbk,alay,altpbk,altlay
    character(len=10)  :: aext
    character(len=*)   :: pbkfilelay,pbkpathlay
    logical            :: isNewFormat
    
    !Determine if using the new dataset naming or old.
    call getFormatType(pbkfilelay, isNewFormat, ierr)
    
    if(isNewFormat) then                   !New Fractional dataset naming
      if(j<=9)then
        write(alay,"('_L',I1)") j
      else
        write(alay,"('_L',I2)") j
      endif
    else
      if(j<=9)then
        write(alay,"(1x,'(',I1,')')") j
      else
        write(alay,"(1x,'(',I2,')')") j
      endif
    endif   
    
    do ks=1,nsed
      if(isNewFormat) then                 !New Fractional dataset naming         
        if (nsed <= 9) then 
          write(apbk,"('_S',I1)") ks                                   
        else
          write(apbk,"('_S',I2)") ks
        endif
      else
        write(apbk,"('_',I2.2)") ks
      endif
      apath = trim(pbkpathlay) // trim(apbk) // alay

      call fileext(trim(pbkfilelay),aext)      
      select case (aext)
      case('h5')
#ifdef XMDF_IO         
        call readscalh5(pbkfilelay,apath,pbk(:,ks,j),ierr)                   
        
        if(ierr < 0) then
          apath = trim(pbkpathlay) // trim(apbk)                             
          call readscalh5(pbkfilelay,apath,pbk(:,ks,j),ierr)
        endif
#endif
      case('txt')
        call readscalTxt(pbkfilelay,pbk(:,ks,j),ierr)    
        if(ierr/=0)then
          apath = trim(pbkpathlay) // trim(apbk)
          call readscalTxt(pbkfilelay,pbk(:,ks,j),ierr)
        endif
      end select
        
    enddo
    
    !Make sure fractions sum 1.0
    do i=1,ncells
!      do ks=1,nsed
!        pbk(i,ks,j) = max(pbk(i,ks,j),1.0e-20)
!      enddo
      pbk(i,:,j) = pbk(i,:,j)/sum(max(pbk(i,:,j),1.0e-20))     !HLI 01/13/2017
    enddo
    
    return
    end subroutine read_pbk    
    
