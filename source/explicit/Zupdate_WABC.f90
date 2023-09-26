!*************************************************************
      subroutine update_WABC()
!*************************************************************
      use EXP_Global_def,    only: dt, drydep
      USE EXP_bndcond_def,   only: wabc, swabc, mwabc, twabc
      use wave_flowgrid_def, only: wavestrx, wavestry
      use bnd_def,   only: nhstr, nthstr, nmhstr, nmhvstr, h_str, th_str, mh_str, mhv_str
      use geo_def,   only: dx, dy, zb
      use met_def,   only: windvar, windconst, tauwindx, tauwindy
      use flow_def,  only: eta, grav
      use prec_def,  only: ikind
      use const_def, only: pi 
      !USE SYNOPTIC_VARS    !Alex, Sep 23, 2009
      
      implicit none
      !local vars  
      integer i,j,id1,id2,loc,id 
      real(ikind) deltax,hplus,hminExp,hgt,tauwindt,radval
      real(ikind) depth1,depth2,volq,dep_ave,detadx  
      
      !calculate advective terms
      if(WABC .or. windconst .or. windvar) then
        if(nHstr .gt. 0) then
          !update momentum equation
          do i=1,nHstr
            do j=1,H_str(i)%ncells - 1
              id1 = H_str(i)%cells(j)
              id2 = H_str(i)%cells(j+1)
              if(id1.eq.id2) cycle
              DELTAX = SWABC(i)%Del(j)
              HPLUS = -zb(id1)+eta(id1) !+ SWABC(i)%eta(j) !swabc is added to 
              hminExp = -zb(id2)+eta(id2) !+ SWABC(i)%eta(j+1) !eta in wse bc routine
              DETADX = ((HPLUS**2.D0 - hminExp**2.D0) - (HPLUS+hminExp)*(-zb(Id1)+zb(id2)))/DELTAX
              HGT = (HPLUS+hminExp)/2.D0
              loc = swabc(I)%rad(j)
              id = SWABC(i)%cell(j)
              if(loc.eq.1) then
              tauwindT = tauwindX(id)*SWABC(i)%sgn(j)
              RADVAL = SWABC(i)%SGN(J)*wavestrx(id)              
              else
              tauwindT = tauwindY(id)*SWABC(i)%sgn(j)
              RADVAL = SWABC(i)%SGN(J)*wavestry(id)          
              endif

              SWABC(i)%QN(j)=(SWABC(i)%Q(j)+DT*(-GRAV*DETADX/2.D0+RADval  + 0.0*tauwindT) )/(1+0.02*DT*abs(SWABC(i)%Q(j))/HGT**2 )
    
              !check wetting and drying
              DEPTH1 = -zb(id1)+ETA(id1) 
              DEPTH2 = -zb(id2)+ETA(id2)
              IF(DEPTH1.LE.DRYDEP.and.DEPTH2.LE.DRYDEP)THEN
                SWABC(i)%QN(j) = 0.0
              ELSEIF(DEPTH1.LE.DRYDEP.and.DEPTH2.GT.DRYDEP)THEN
                IF(SWABC(i)%QN(j) .LT. 0.d0) THEN      
                  SWABC(i)%QN(j) = 0.D0
                ENDIF
              ELSEIF(DEPTH1.GT.DRYDEP.and.DEPTH2.LE.DRYDEP)THEN
                IF(SWABC(i)%QN(j) .GT. 0.D0) THEN
                  SWABC(i)%QN(j) = 0.D0
                ENDIF
              ENDIF

              !update wse
              VOLQ = SWABC(i)%QN(j)*DT*SWABC(i)%wdth(j)
              SWABC(i)%Q(j) = SWABC(i)%QN(j)
              SWABC(i)%eta(j) = SWABC(i)%eta(j)+VOLQ/(dx(id1)*dy(id1))
              SWABC(i)%eta(j+1)=SWABC(i)%eta(j+1)-VOLQ/(dx(id2)*dy(id2))
            enddo
          enddo

        endif ! H_single

        if(nTHstr .gt. 0) then
          !update momentum equation
          i=1
          do j=1,nTHstr
              id1 = TH_str(i)%cells(j)
              id2 = TH_str(i)%cells(j+1)
            DELTAX = TWABC(i)%Del(j)
            HPLUS = -zb(id1)+eta(id1) !+ TWABC(i)%eta(j) !TWABC is added to 
            hminExp = -zb(id2)+eta(id2) !+ TWABC(i)%eta(j+1) !eta in wse bc routine
            DETADX = ((HPLUS**2.D0 - hminExp**2.D0) - (HPLUS+hminExp)*(-zb(Id1)+zb(id2)))/DELTAX
            HGT = (HPLUS+hminExp)/2.D0
            loc = TWABC(I)%rad(j)
            id = TWABC(i)%cell(j)

              if(loc.eq.1) then
              tauwindT = tauwindX(id)*TWABC(i)%sgn(j)
               RADVAL = TWABC(i)%SGN(J)*wavestrx(id)           
              else
              tauwindT = tauwindY(id)*TWABC(i)%sgn(j)
              RADVAL = TWABC(i)%SGN(J)*wavestry(id)           
              endif            
            TWABC(i)%QN(j) = (TWABC(i)%Q(j)+DT*(-GRAV*DETADX/2.D0+RADval + 0.0*tauwindT) )/(1+0.02*DT*abs(TWABC(i)%Q(j))/HGT**2 )

            !check wetting and drying
            DEPTH1 = -zb(id1)+ETA(id1) 
            DEPTH2 = -zb(id2)+ETA(id2)
            IF(DEPTH1.LE.DRYDEP.and.DEPTH2.LE.DRYDEP)THEN
              TWABC(i)%QN(j) = 0.0
            ELSEIF(DEPTH1.LE.DRYDEP.and.DEPTH2.GT.DRYDEP)THEN
              IF(TWABC(i)%QN(j) .LT. 0.d0) THEN      
                TWABC(i)%QN(j) = 0.D0
              ENDIF
            ELSEIF(DEPTH1.GT.DRYDEP.and.DEPTH2.LE.DRYDEP)THEN
              IF(TWABC(i)%QN(j) .GT. 0.D0) THEN
                TWABC(i)%QN(j) = 0.D0
              ENDIF
            ENDIF

            !update wse
            VOLQ = TWABC(i)%QN(j)*DT*TWABC(i)%wdth(j)
            TWABC(i)%Q(j) = TWABC(i)%QN(j)
            TWABC(i)%eta(j) = TWABC(i)%eta(j)+VOLQ/(dx(id1)*dy(id1))
            TWABC(i)%eta(j+1) = TWABC(i)%eta(j+1)-VOLQ/(dx(id2)*dy(id2))
          enddo

        endif ! H_tide

        if(nMHstr .gt. 0) then

          !update momentum equation
          do i=1,nMHVstr
            do j=1,MH_str(i)%ncells - 1
              id1 = MH_str(i)%cells(j)
              id2 = MH_str(i)%cells(j+1)
              DELTAX = MWABC(i)%Del(j)
              Dep_Ave = (-zb(id1) -zb(id2)+eta(id1)+eta(id2) - mwabc(i)%eta(j) - mwabc(i)%eta(j+1))/2.0
              HPLUS = dep_ave + mwabc(i)%eta(j) !-zb(id1)+eta(id1) !+ MWABC(i)%eta(j) !swabc is added to 
              hminExp = dep_ave + mwabc(i)%eta(j+1) !-zb(id2)+eta(id2) !+ MWABC(i)%eta(j+1) !eta in wse bc routine
              DETADX = ((HPLUS**2.D0 - hminExp**2.D0) - 0.0*(HPLUS+hminExp)*(-zb(Id1)+zb(id2)) )/DELTAX
            HGT = (HPLUS+hminExp)/2.D0
            loc = MWABC(I)%rad(j)
              id = MWABC(i)%cell(j)

              if(loc.eq.1) then
                tauwindT = tauwindX(id)*MWABC(i)%sgn(j)
                RADVAL = MWABC(i)%SGN(J)*wavestrx(id)             
              else
                tauwindT = tauwindY(id)*MWABC(i)%sgn(j)
                RADVAL = MWABC(i)%SGN(J)*wavestry(id)            
              endif              
              MWABC(i)%QN(j)=(MWABC(i)%Q(j)+DT*(-GRAV*DETADX/2.D0 + 0.0*RADval + 0.0*tauwindT) )/(1+0.02*DT*abs(MWABC(i)%Q(j))/HGT**2 )

              !check wetting and drying
              DEPTH1 = -zb(id1)+ETA(id1) 
              DEPTH2 = -zb(id2)+ETA(id2)
              IF(DEPTH1.LE.DRYDEP.and.DEPTH2.LE.DRYDEP)THEN
                MWABC(i)%QN(j) = 0.0
              ELSEIF(DEPTH1.LE.DRYDEP.and.DEPTH2.GT.DRYDEP)THEN
                IF(MWABC(i)%QN(j) .LT. 0.d0) THEN      
                  MWABC(i)%QN(j) = 0.D0
                ENDIF
              ELSEIF(DEPTH1.GT.DRYDEP.and.DEPTH2.LE.DRYDEP)THEN
                IF(MWABC(i)%QN(j) .GT. 0.D0) THEN
                  MWABC(i)%QN(j) = 0.D0
                ENDIF
              ENDIF

              !update wse
              VOLQ = MWABC(i)%QN(j)*DT*MWABC(i)%wdth(j)
              MWABC(i)%Q(j) = MWABC(i)%QN(j)
              MWABC(i)%eta(j) = MWABC(i)%eta(j)+VOLQ/(dx(id1)*dy(id1))
              MWABC(i)%eta(j+1)=MWABC(i)%eta(j+1)-VOLQ/(dx(id2)*dy(id2))
            enddo
          enddo

        endif ! H_multi
        
        if(nMHVstr .gt. 0) then

          !update momentum equation
          do i=1,nMHVstr
            do j=1,MH_str(i)%ncells - 1
              id1 = MHV_str(i)%cells(j)
              id2 = MHV_str(i)%cells(j+1)
              DELTAX = MWABC(i)%Del(j)
              Dep_Ave = (-zb(id1) + (-zb(id2)) +eta(id1)+eta(id2) - mwabc(i)%eta(j) - mwabc(i)%eta(j+1))/2.0
              HPLUS = dep_ave + mwabc(i)%eta(j) !-zb(id1)+eta(id1) !+ MWABC(i)%eta(j) !swabc is added to 
              hminExp = dep_ave + mwabc(i)%eta(j+1) !-zb(id2)+eta(id2) !+ MWABC(i)%eta(j+1) !eta in wse bc routine
              DETADX = ((HPLUS**2.D0 - hminExp**2.D0) - 0.0*(HPLUS+hminExp)*(-zb(Id1)+zb(id2)) )/DELTAX
            HGT = (HPLUS+hminExp)/2.D0
            loc = MWABC(I)%rad(j)
              id = MWABC(i)%cell(j)

              if(loc.eq.1) then
                tauwindT = tauwindX(id)*MWABC(i)%sgn(j)
                RADVAL = MWABC(i)%SGN(J)*wavestrx(id)             
              else
                tauwindT = tauwindY(id)*MWABC(i)%sgn(j)
                RADVAL = MWABC(i)%SGN(J)*wavestry(id)            
              endif              
              MWABC(i)%QN(j)=(MWABC(i)%Q(j)+DT*(-GRAV*DETADX/2.D0 + 0.0*RADval + 0.0*tauwindT) )/(1+0.02*DT*abs(MWABC(i)%Q(j))/HGT**2 )

              !check wetting and drying
              DEPTH1 = -zb(id1)+ETA(id1) 
              DEPTH2 = -zb(id2)+ETA(id2)
              IF(DEPTH1.LE.DRYDEP.and.DEPTH2.LE.DRYDEP)THEN
                MWABC(i)%QN(j) = 0.0
              ELSEIF(DEPTH1.LE.DRYDEP.and.DEPTH2.GT.DRYDEP)THEN
                IF(MWABC(i)%QN(j) .LT. 0.d0) THEN      
                  MWABC(i)%QN(j) = 0.D0
                ENDIF
              ELSEIF(DEPTH1.GT.DRYDEP.and.DEPTH2.LE.DRYDEP)THEN
                IF(MWABC(i)%QN(j) .GT. 0.D0) THEN
                  MWABC(i)%QN(j) = 0.D0
                ENDIF
              ENDIF

              !update wse
              VOLQ = MWABC(i)%QN(j)*DT*MWABC(i)%wdth(j)
              MWABC(i)%Q(j) = MWABC(i)%QN(j)
              MWABC(i)%eta(j) = MWABC(i)%eta(j)+VOLQ/(dx(id1)*dy(id1))
              MWABC(i)%eta(j+1)=MWABC(i)%eta(j+1)-VOLQ/(dx(id2)*dy(id2))
            enddo
          enddo
        endif ! H_multi        
      endif  ! Radstr

      return  
      end subroutine
