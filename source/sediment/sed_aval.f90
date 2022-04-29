!********************************************************************
    subroutine avalanche
! Calculates the avalanching using a relaxation approach
! written by Alex Sanchez, USACE-CHL;  Weiming Wu, NCCHE

!********************************************************************
    use size_def
    use geo_def
    use sed_def
    use diag_lib
    use prec_def
    implicit none        
    integer :: i,ii,j,k,ks,na,ncw,ncs,nck,iter,nat,iddzbavmax,iddzbavmaxt
    real(ikind) :: slope,a_repose2,val
    real(ikind) :: zbav(ncellsD),dzbavmax,dzbavmaxt
    real(ikind) :: dzbav(nmaxfaces,ncellsD) !Bed changes are saved at faces in order to avoid race conditions
    character(len=100) :: msg
    logical :: hdr
    
656 format(12x,I8,3x,I9)
234 format(3x,'dzbav(',I6,') = ',F9.5)
    
    !************************************************************************************************
    !Note: The sign of the slope should be changed to be more clear (slope>0 should be upslope)
    !and a corresponding change in the dzbav calculation
    !************************************************************************************************
    
    hdr = .false.
    a_repose2 = a_repose + 1.0e-4 !Used to avoid cycling due to precision error

!$OMP PARALLEL DO PRIVATE(i)
    do i=1,ncells
      zbav(i) = zb(i)  
      dzbav(1:ncface(i),i) = 0.0
    enddo
!$OMP END PARALLEL DO    
    
    do iter=1,nmaxaval
      !Note: For some reason the reduction OMP clause is not working for this case. 
      !Instead the reduction is done manually with nat and the critical section below. 
      na = 0; nat = 0
!$OMP PARALLEL FIRSTPRIVATE(nat)
      !--- Initialize avalanche bed changes --------------
      !Note: saved at cell faces in order to avoid a race condition with OMP
      !--- Calculate avalanche bed changes ------------------
!$OMP DO PRIVATE(i,ii,ncw,ncs,slope,val)
      do ii=1,ncellsimple      
        i = idcellsimple(ii)
        !West simple face
        ncw = cell2cell(4,i)
        if(ncw<=ncells)then  
          slope = (zbav(ncw)-zbav(i))/dc(4,i)
          if(slope<-a_repose2 .and. zbav(i)>hardzb(i))then !Downslope from i to ncw, slope<0
            val = 0.5*(slope+a_repose) !0.5 from mass balance  ! Wu: dc(i==>ncw)/(dx(i)+dx(ncw))=0.5 for simple mesh
            dzbav(4,i) = dx(ncw)*val
            dzbav(2,ncw) = -dx(i)*val
            nat = nat + 1
          elseif(slope>a_repose2 .and. zbav(ncw)>hardzb(ncw))then !Upslope, slope>0
            val = 0.5*(slope-a_repose) !0.5 from mass balance  ! Wu: dc(i==>ncw)/(dx(i)+dx(ncw))=0.5 for simple mesh
            dzbav(4,i) = dx(ncw)*val
            dzbav(2,ncw) = -dx(i)*val
            nat = nat + 1
          endif       
        endif   
        !South simple face
        ncs = cell2cell(3,i)
        if(ncs<=ncells)then  
          slope = (zbav(ncs)-zbav(i))/dc(3,i)
          if(slope<-a_repose2 .and. zbav(i)>hardzb(i))then !Downslope from i to ncs, slope<0
            val = 0.5*(slope+a_repose) !0.5 from mass balance  ! Wu: dc(i==>ncw)/(dx(i)+dx(ncw))=0.5 for simple mesh
            dzbav(3,i) = dy(ncs)*val
            dzbav(1,ncs) = -dy(i)*val
            nat = nat + 1
          elseif(slope>a_repose2 .and. zbav(ncs)>hardzb(ncs))then !Upslope, slope>0
            val = 0.5*(slope-a_repose) !0.5 from mass balance  ! Wu: dc(i==>ncw)/(dx(i)+dx(ncw))=0.5 for simple mesh
            dzbav(3,i) = dy(ncs)*val
            dzbav(1,ncs) = -dy(i)*val
            nat = nat + 1
          endif
        endif 
      enddo  !ncellsimple
!$OMP END DO
!$OMP DO PRIVATE(i,ii,j,k,nck,slope,val)
      do ii=1,ncelljoint
        i = idcelljoint(ii)
        do j=1,nxyface(i) !No repeat cell faces
          k = kxyface(j,i)
          nck = cell2cell(k,i)
          if(nck<=ncells)then  
            slope = (zbav(nck)-zbav(i))/dc(k,i)
            if(slope<-a_repose2 .and. zbav(i)>hardzb(i))then !Downslope from i to nck, slope<0
              val = dc(k,i)*(slope+a_repose)/(areap(i)+areap(nck))  !Wu
              dzbav(k,i) = val*areap(nck)
              dzbav(llec2llec(k,i),nck) = -val*areap(i)
              nat = nat + 1
            elseif(slope>a_repose2 .and. zbav(nck)>hardzb(nck))then !Upslope, slope>0
              val = dc(k,i)*(slope-a_repose)/(areap(i)+areap(nck))   !Wu
              dzbav(k,i) = val*areap(nck)
              dzbav(llec2llec(k,i),nck) = -val*areap(i)
              nat = nat + 1
            endif
          endif
        enddo !j-no-repeat face
      enddo  !ncelljoint
!$OMP END DO
!$OMP DO PRIVATE(i,ii,j,k,nck,slope,val)
      do i=1,ncellpoly
        do j=1,nxyface(i) !No repeat cell faces
          k = kxyface(j,i)
          nck = cell2cell(k,i)
          if(nck<=ncells)then  
            slope = (zbav(nck)-zbav(i))/dc(k,i)
            if(slope<-a_repose2 .and. zbav(i)>hardzb(i))then !Downslope from i to nck, slope<0
              val = dc(k,i)*(slope+a_repose)/(areap(i)+areap(nck))  !Wu
              dzbav(k,i) = val*areap(nck)
              dzbav(llec2llec(k,i),nck) = -val*areap(i)
              nat = nat + 1
            elseif(slope>a_repose2 .and. zbav(nck)>hardzb(nck))then !Upslope, slope>0
              val = dc(k,i)*(slope-a_repose)/(areap(i)+areap(nck))   !Wu
              dzbav(k,i) = val*areap(nck)
              dzbav(llec2llec(k,i),nck) = -val*areap(i)
              nat = nat + 1
            endif
          endif
        enddo !j-no-repeat face
      enddo !ncellpoly
!$OMP END DO
!$OMP CRITICAL
      na = na + nat
!$OMP END CRITICAL
!$OMP END PARALLEL
      if(na==0) exit
    !--- Update bed elevation -----
!$OMP PARALLEL DO PRIVATE(i)
      do i=1,ncells
        zbav(i) = zbav(i) + relax_aval*sum(dzbav(1:ncface(i),i))
        dzbav(1:ncface(i),i) = 0.0
      enddo
!$OMP END PARALLEL DO
      if(.not.hdr)then
        write(msg,656) iter,na
        call diag_print_message(' Avalanching:  Iteration    Number',msg)
        hdr = .true.
      endif
      if(mod(iter,50)==0)then
        write(msg,656) iter,na
        call diag_print_message(msg)
      endif
    enddo !iter

    if(na>0)then
      if(mod(iter,50)/=0 .and. iter<nmaxaval)then
        write(msg,656) iter,na
        call diag_print_message(msg)
      endif
      iddzbavmax = 0
      iddzbavmaxt = 0
      dzbavmaxt = 0.0
      dzbavmax = 0.0
!!$OMP PARALLEL
!!$OMP DO PRIVATE(i,val,iddzbavmaxt,dzbavmaxt)
      do i=1,ncells
        val = zb(i) - zbav(i)
        dzb(i) = dzb(i) + val
        zb(i) = zbav(i)
        if(abs(val)>abs(dzbavmaxt))then
          dzbavmaxt = val
          iddzbavmaxt = i
        endif
      enddo
!!$OMP END DO
      !Distribute avalanching according for fractional bed change
      if(nsed>1)then
!!$OMP DO PRIVATE(i,ks,val)
        do i=1,ncells
          val = zb(i) - zbav(i)  
          do ks=1,nsed
            dzbk(i,ks) = dzbk(i,ks) + pbk(i,ks,1)*val 
          enddo
        enddo
!!$OMP END DO
      endif
!!$OMP CRITICAL
    if(abs(dzbavmaxt)>abs(dzbavmax))then
      dzbavmax = dzbavmaxt
      iddzbavmax = iddzbavmaxt
    endif
!!$OMP END CRITICAL
!!$OMP END PARALLEL
      if(allocated(mapid))then
        iddzbavmax = mapid(iddzbavmax)
      endif
      write(msg,234) iddzbavmax,dzbavmax
      call diag_print_message(msg)
    endif
    
    return
    end subroutine avalanche