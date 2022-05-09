!=================================================================
! CMS Dredging Routines
!
! Contains the following:
!   dredge_op - checks and modifies status of dredging operations activity
!   trigger_it - checks to see if dredging should be triggered
!   dredge_it - moves sediment for any active dredging operation
!   place_it - places any dredged material
! written by Chris Reed
!=================================================================
    
!**************************************************************
    subroutine dredge_op()
! checks and modifies status of dredging operations activity
! written by Chris Reed
!**************************************************************    
    use dredge_def
    use geo_def, only: zb,dx,dy  
    use sed_def, only: hardbottom
    implicit none
    integer :: i,k,NCnt,ID_t
    real(ikind) :: vol_tot,diff,bed_t
    
    DO k=1,ndredge_operations 
      !This is done here to calculate volume in the drege source areas (differnt thatn remaining dredge colume for method 1)
      Ncnt = dredge_operations(k)%NumDredgeAreaCells      
      vol_tot = 0.0
      do i=1,Ncnt
        ID_t=dredge_operations(k)%DredgeAreaCells(i)
        bed_t= -zb(ID_t)
        Diff = dredge_operations(k)%Dredge_Depth(i) -  bed_t
        vol_tot = vol_tot + max(diff,0.0)*dx(ID_t)*dy(ID_t)
      enddo 
      DredgeTS_Vars(k,8) = vol_tot  !for printing         
      DredgeTS_Vars(k,7) = vol_tot  !for printing  - shows volume in templete prior to triggering dredging
    
      !check dredging status and adjust it if required
      !if dredging is active, see if it can be stopped if all cells are below dredge_depth
      if(dredge_operations(k)%active) then  !check to see if it can be turned off
        dredge_operations(k)%active = .false.
        do i=1,dredge_operations(k)%NumDredgeAreaCells
          if(zb(dredge_operations(k)%DredgeAreaCells(i)) > -dredge_operations(k)%Dredge_Depth(i)) dredge_operations(k)%active = .true.
        enddo
      endif
     
      !if dredging is not active, then see if needs to become active
      if(.not. dredge_operations(k)%active) call trigger_it(k) !check to see if it needs to be activated  
      if(dredge_operations(k)%active .and. dredge_operations(k)%Trigger_approach == 4) call cancel_it(k)  !only check for stopping if option 4 is used

      IF(dredge_operations(k)%active) THEN  
        !if dredging is active, calculate removal volume (cell by cell) and then put it in the placement area
        call dredge_it(k)
        call Placement_allocation(k)
      ENDIF
     
      !cumulative outputs updated here just in case dredging has ceased and active = .false.
      DredgeTS_Vars(k,6) = dredge_operations(k)%placement_excess 
      DredgeTS_Vars(k,4) = dredge_operations(k)%total_dredged_vol      
      DredgeTS_Vars(k,5) = dredge_operations(k)%total_placed_vol 
     
      if (write_dredge_diag) call WriteDredgeTS(k)
     
      if(hardbottom) call dredge_hardbottom()
     
    ENDDO

    return
    end subroutine dredge_op

!**************************************************************
    subroutine trigger_it(k)
! checks and modifies status of dredging operations activity
! written by Chris Reed
!**************************************************************    
    use dredge_def
    use geo_def, only: zb,dx,dy  
    use comvarbl, only: timehrs
    implicit none
    integer i,k
    real(ikind) sum_area,sum_area_exceed,area,sum_vol,ratio
    
    select case(dredge_operations(k)%trigger_approach)        
    case(1)  !any cell above trigger depth
      do i=1,dredge_operations(k)%NumDredgeAreaCells
        if(zb(dredge_operations(k)%DredgeAreaCells(i)) > -dredge_operations(k)%Trigger_Depth) dredge_operations(k)%active = .true.
      enddo
         
    case(2)  !volume above dredge depth is greater than specified volume
      sum_vol = 0.0   
      do i=1,dredge_operations(k)%NumDredgeAreaCells
        area = dx(dredge_operations(k)%DredgeAreaCells(i))*dy(dredge_operations(k)%DredgeAreaCells(i))
        sum_vol = sum_vol + area*max(dredge_operations(k)%Dredge_Depth(i) + zb(dredge_operations(k)%DredgeAreaCells(i)),0.0)
      enddo
      if(sum_vol >= dredge_operations(k)%Trigger_Vol) dredge_operations(k)%active = .true.  
        
    case(3)  !specified percentage of area above trigger depth        
      sum_area_exceed = 0.0
      sum_area = 0.0
      do i=1,dredge_operations(k)%NumDredgeAreaCells
        area = dx(dredge_operations(k)%DredgeAreaCells(i))*dy(dredge_operations(k)%DredgeAreaCells(i))
        sum_area = sum_area + area
        if(zb(dredge_operations(k)%DredgeAreaCells(i)) > -dredge_operations(k)%Trigger_Depth) sum_area_exceed = sum_area_exceed + area
      enddo   
      ratio = 100*sum_area_exceed/(sum_area + 1.e-20)
      if(ratio >=  dredge_operations(k)%Trigger_Percentage ) dredge_operations(k)%active = .true.
          
    case(4)  !within a specified dredging time interval
      dredge_operations(k)%active = .false.    
      do i=1,dredge_operations(k)%num_trigger_intervals
        if (timehrs >= dredge_operations(k)%Trigger_Start(i) .and. timehrs < dredge_operations(k)%trigger_finish(i)) then
          dredge_operations(k)%active = .true.
          !write(*,*)k,i,dredge_operations(k)%Trigger_Start(i),dredge_operations(k)%trigger_finish(i),dredge_operations(k)%active 
        endif
      enddo
            
    case default
      write(*,*)'Approach specified for triggering dredging is not recognized'
      write(*,*)'Approach specified is: ',dredge_operations(k)%trigger_approach
        
    end select
        
    return
    end subroutine trigger_it    

!**************************************************************
    subroutine cancel_it(k)
! checks and modifies status of dredging operations activity
! written by Chris Reed
!**************************************************************    
    use dredge_def
    use geo_def, only: zb,dx,dy  
    use comvarbl, only: timehrs
    implicit none
    integer i,k

    
   !operation is currently set to .true., so see if it should be stopped
   !because specified active interval is finsihed
   if(dredge_operations(k)%num_trigger_intervals > 1) then
     do i=1,dredge_operations(k)%num_trigger_intervals-1
       if(timehrs >= dredge_operations(k)%Trigger_Finish(i) .and. timehrs < dredge_operations(k)%trigger_start(i+1)) then
         dredge_operations(k)%active = .false.
       endif
     enddo
   endif
   if(timehrs >= dredge_operations(k)%Trigger_Finish(dredge_operations(k)%num_trigger_intervals)) dredge_operations(k)%active = .false.            
        
   return
   end subroutine cancel_it       
    
    
!**************************************************************
   subroutine dredge_it(k)
! moves sediment for any active dredging operation
! written by Chris Reed
!**************************************************************    
    use dredge_def
    use geo_def, only: zb,dx,dy
    use sed_def, only: singlesize    
    implicit none  
    integer k,ncnt,i,ival,cell_inc
    real(ikind) Diff,vol_tot,vol_limit,area,tmp_vol,tot_vol_dredged  !,dep_inc
    real(ikind) vol_cell,vol_available,dep
    real*8 sumvol   
    !real(ikind), allocatable :: bed(:),del(:),bed1(:)
    real*8, allocatable :: bed(:),del(:),bed1(:)    
    integer, allocatable :: cells(:)
    
    !copy values into local arrays for easier coding
    Ncnt = dredge_operations(k)%NumDredgeAreaCells    
    allocate(cells(Ncnt),bed(Ncnt),Del(Ncnt),bed1(Ncnt))   

    do i=1,Ncnt
      cells(i) = dredge_operations(k)%DredgeAreaCells(i) 
      bed(i) = -zb(cells(i)) 
      bed_change(i) = 0.0
    enddo     
    
    ! write(*,*)'------------------------------------------------------'
    ! write(*,*)'ncnt = ',ncnt
    ! write(777,*)'vol_tot = ',vol_tot
    
    !remove volume from cells based on selected approach
    select case(dredge_operations(k)%dredge_approach)
        
    case(1) !remove from shallowest cells first
      !calculate potential dredge volume
      vol_tot = 0.0
      do i=1,Ncnt
        Diff = dredge_operations(k)%Dredge_Depth(i) -  bed(i)
        vol_tot = vol_tot + max(diff,0.0)*dx(cells(i))*dy(cells(i))
      enddo 
      DredgeTS_Vars(k,7) = vol_tot  !for printing
    
      !limit volume by less of dredge rate and potential volume
      vol_limit = dredge_operations(k)%rate*dredge_time_lapse
      vol_tot = min(vol_limit,vol_tot)
      !dredge_operations(k)%placement_vol = vol_tot
      DredgeTS_Vars(k,9) = vol_tot     !for printing       
        
      call sort2DP(bed,cells,ncnt)
      bed1=bed       

      !calculate depth intervals between successive cells
      do i=1,Ncnt-1
        del(i) = bed(i+1)-bed(i)
      enddo
      del(Ncnt) = 0.0
      !find the cells that will be dredged for this time interval
      area=0
      tmp_vol = 0.0
      ival = 0
      do i=1,ncnt
        area = area + dx(cells(i))*dy(cells(i))
        diff = area*del(i)
        tmp_vol = tmp_vol + diff
        !write(777,"(i5,)")i,tmp_vol,vol_tot,del(i),bed(i+1),bed(i)
        if(tmp_vol > vol_tot) then  !found cells
          ival = i
          dep = bed(i) + (vol_tot - (tmp_vol-diff))/area
          !delta = (vol_tot - (tmp_vol-diff))/area
          !write(777,*)'dep = ',dep
          go to 100    !This should probably be an EXIT
        endif
      enddo
100   continue 
      if(ival == 0) then  !depth is set as remaining volume below last cell
        ival = ncnt
        dep = bed(ncnt) + (vol_tot - tmp_vol)/area
        !delta = (vol_tot - tmp_vol)/area
      endif
      !write(777,*)'ival= ',ival
      tot_vol_dredged = 0
      do i=1,ival
        area = dx(cells(i))*dy(cells(i))
        !dep_inc = dep-bed(i)
        bed_change(i) = dep-bed(i)
        !tot_vol_dredged = tot_vol_dredged + dep_inc*area
        tot_vol_dredged = tot_vol_dredged + bed_change(i)*area        
        !bed(i) = bed(i) + dep_inc
        bed(i) = bed(i) + bed_change(i)
        !zb(cells(i)) = zb(cells(i)) - dep_inc
        zb(cells(i)) = zb(cells(i)) - bed_change(i)      
        !write(777,"(i5,4e15.7)")i,bed(i),zb(cells(i)),tot_vol_dredged        
      enddo
        
      !dredge_operations(k)%placement_vol = tot_vol_dredged
      !for printing
      sumvol=0.0
      do i=1,ncnt
        sumvol=sumvol + (bed1(i)-bed(i))*dx(cells(i))*dy(cells(i))
      enddo        
      DredgeTS_Vars(k,1) = -sumvol  !tot_vol_dredged       
      dredge_operations(k)%placement_vol= -sumvol
      dredge_operations(k)%total_dredged_vol = dredge_operations(k)%total_dredged_vol - sumvol
       
    case(2) ! remove from cell closest to starting cell first, cells ranked by proximity in dredge_init()
      vol_tot = 0.0
      
      cell_inc = dredge_operations(k)%cell_inc
      
      do i=dredge_operations(k)%cell_inc,Ncnt
        Diff = dredge_operations(k)%Dredge_Depth(i) -  bed(i)
        vol_tot = vol_tot + max(diff,0.0)*dx(cells(i))*dy(cells(i))
      enddo 
      DredgeTS_Vars(k,7) = vol_tot  !for printing
    
      !limit volume by less of dredge rate and potential volume
      vol_limit = dredge_operations(k)%rate*dredge_time_lapse
      vol_tot = min(vol_limit,vol_tot)
      !dredge_operations(k)%placement_vol = vol_tot
      DredgeTS_Vars(k,9) = vol_tot     !for printing      
      
      bed1=bed
      vol_available = vol_tot
      do i=dredge_operations(k)%cell_inc,ncnt  
        if(vol_available > 0.0) then  !dredge cell 
          if(bed(i) < dredge_operations(k)%Dredge_Depth(i)) then ! this cell still needs dredging
            vol_cell = dx(cells(i))*dy(cells(i))*(dredge_operations(k)%Dredge_Depth(i) -  bed(i))
            if(vol_cell > vol_available) then
              !zb(cells(i)) = zb(cells(i)) - vol_available/(dx(cells(i))*dy(cells(i)))
              !bed(i) = bed(i) + vol_available/(dx(cells(i))*dy(cells(i)))
              bed_change(i) = vol_available/(dx(cells(i))*dy(cells(i)))
              zb(cells(i)) = zb(cells(i)) - bed_change(i)
              bed(i) = bed(i) + bed_change(i)              
              vol_available = 0.0
              cell_inc = i               
            else
              !zb(cells(i)) = zb(cells(i)) - vol_cell/(dx(cells(i))*dy(cells(i)))
              !bed(i) = bed(i) +  vol_cell/(dx(cells(i))*dy(cells(i)))  
              bed_change(i) = vol_available/(dx(cells(i))*dy(cells(i))) 
              zb(cells(i)) = zb(cells(i)) - bed_change(i)
              bed(i) = bed(i) +  bed_change(i)               
              vol_available = vol_available - vol_cell
              cell_inc = i+1
            endif
          endif
        else                                 
          go to 200    !This should probably be an EXIT statement rather than a GOTO.   
        endif  
      enddo
200   continue        
      dredge_operations(k)%cell_inc  = cell_inc      
        
      !see if dregding is complete
      if(dredge_operations(k)%cell_inc > ncnt) then
        dredge_operations(k)%active = .false.
        dredge_operations(k)%cell_inc = 1  
      endif
      if(dredge_operations(k)%cell_inc == ncnt) then
        if(bed(ncnt) >= dredge_operations(k)%Dredge_Depth(ncnt)) then
          dredge_operations(k)%active = .false.
          dredge_operations(k)%cell_inc = 1   
        endif
      endif
        
      !dredge_operations(k)%placement_vol = vol_tot - vol_available
        
      !for printing
      sumvol=0.0
      do i=1,ncnt
        sumvol=sumvol + (bed1(i)-bed(i))*dx(cells(i))*dy(cells(i))
      enddo        
      DredgeTS_Vars(k,1) = -sumvol  !tot_vol_dredged         
      !DredgeTS_Vars(k,1) = vol_tot - vol_available       
      dredge_operations(k)%placement_vol = -sumvol   
      dredge_operations(k)%total_dredged_vol = dredge_operations(k)%total_dredged_vol - sumvol                
        
    case default
      write(*,*)'Approach specified for dredging source area is not recognized'
      write(*,*)'Approach specified is: ',dredge_operations(k)%dredge_approach
        
    end select   
       
    if(.not. singlesize) call dredge_bed_sort(ncnt,cells)   
    
    
    return
    end subroutine dredge_it
    
!**************************************************************
    subroutine Placement_allocation(k)
! checks and modifies status of dredging operations activity
! written by Chris Reed
!**************************************************************    
    use dredge_def
    use geo_def, only: zb,dx,dy   
    implicit none
    integer k,j,m,num_remain
    real(ikind) placement_vol,excess,amount_placed,sum_percentage
    real(ikind), allocatable :: placement_Vol_byarea(:),PAallocation(:)    

     excess=0.0
     
     !do this so that for no placement areas all dredge materail goes to 'excess'
     if(dredge_operations(k)%NumPlacementAreas == 0) then
         excess = dredge_operations(k)%placement_vol
         DredgeTS_Vars(k,3) = excess
         !dredge_operations(k)%placement_excess = dredge_operations(k)%placement_excess + excess
     endif
     
     IF(dredge_operations(k)%NumPlacementAreas == 1) THEN   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       j=1 
       placement_vol = dredge_operations(k)%placement_vol
       call place_it(k,j,placement_vol,amount_placed,excess)  !placed dredged volume for kth operation and jth placement area
       DredgeTS_Vars(k,2) = amount_placed
       DredgeTS_Vars(k,3) = excess
       dredge_operations(k)%placement_excess = dredge_operations(k)%placement_excess + excess
       dredge_operations(k)%total_placed_vol = dredge_operations(k)%total_placed_vol + amount_placed     
     ELSE
       if(dredge_operations(k)%PAmethod == 1) then  ! allocate by order!!!!!!!!!!!!!!!!!!!!!!!
         placement_vol = dredge_operations(k)%placement_vol 
         DredgeTS_Vars(k,2) = 0.0
         do j=1,dredge_operations(k)%NumPlacementAreas
           call place_it(k,j,placement_vol,amount_placed,excess)  !placed dredged volume for kth operation and jth placement area
           placement_vol = excess
           DredgeTS_Vars(k,2) =  DredgeTS_Vars(k,2) + amount_placed
           dredge_operations(k)%placed_vol_by_area(j) = amount_placed
           dredge_operations(k)%total_placed_by_area(j) = dredge_operations(k)%total_placed_by_area(j) + amount_placed
           dredge_operations(k)%total_placed_vol = dredge_operations(k)%total_placed_vol + amount_placed    
         enddo 
         !DredgeTS_Vars(k,2) = amount_placed  !dredge_operations(k)%placement_vol - excess
         DredgeTS_Vars(k,3) = excess
         dredge_operations(k)%placement_excess = dredge_operations(k)%placement_excess + excess      
       elseif(dredge_operations(k)%PAmethod .eq. 2) then  ! allocate by percentage!!!!!!!!!!!!!!!!!!!!!!!
         allocate(placement_Vol_byarea(dredge_operations(k)%NumPlacementAreas))
         allocate(PAallocation(dredge_operations(k)%NumPlacementAreas))
         do j=1,dredge_operations(k)%NumPlacementAreas
           placement_Vol_byarea(j) = dredge_operations(k)%placement_vol*dredge_operations(k)%DredgePlacementAllocation(j)/100.0  
         enddo
         do j=1,dredge_operations(k)%NumPlacementAreas
           PAallocation(j) = dredge_operations(k)%DredgePlacementAllocation(j)
         enddo
     
         DredgeTS_Vars(k,2) =  0.0
         do j=1,dredge_operations(k)%NumPlacementAreas
           placement_vol =  placement_Vol_byarea(j) + excess*PAallocation(j)/100.0   
           call place_it(k,j,placement_vol,amount_placed,excess)  !placed dredged volume for kth operation and jth placement area
           DredgeTS_Vars(k,2) =  DredgeTS_Vars(k,2) + amount_placed
           dredge_operations(k)%placed_vol_by_area(j) = amount_placed
           dredge_operations(k)%total_placed_by_area(j) = dredge_operations(k)%total_placed_by_area(j) + amount_placed     
           dredge_operations(k)%total_placed_vol = dredge_operations(k)%total_placed_vol + amount_placed
           if(excess .gt. 0.0 .and. j .lt. dredge_operations(k)%NumPlacementAreas) then  !re-adjust %s based on subsequent placement areas
             sum_percentage = 0.0
             num_remain = 0
             do m=j+1,dredge_operations(k)%NumPlacementAreas
               sum_percentage = sum_percentage + PAallocation(m)
               num_remain = num_remain + 1
             enddo  
             if(sum_percentage .gt. 0.0) then !reallocate
               do m=j+1,dredge_operations(k)%NumPlacementAreas
                 !write(*,*)'reall1: ',m,PAallocation(m),excess              
                 PAallocation(m) = 100*PAallocation(m)/sum_percentage
                 !write(*,*)'reall2: ',m,PAallocation(m)
               enddo 
             else
               do m=j+1,dredge_operations(k)%NumPlacementAreas
                 PAallocation(m) = 100.0/num_remain
               enddo               
             endif
          endif
        enddo
        !DredgeTS_Vars(k,2) = amount_placed  !dredge_operations(k)%placement_vol - excess
        DredgeTS_Vars(k,3) = excess
        dredge_operations(k)%placement_excess = dredge_operations(k)%placement_excess + excess
      endif
    ENDIF
     
    return
    end subroutine Placement_allocation        
        
    
!**************************************************************
   subroutine place_it(k,j,placement_vol,amount_placed,excess)
! moves sediment for any active dredging operation
! written by Chris Reed
!**************************************************************    
    use dredge_def
    use geo_def, only: zb,dx,dy
    use sed_def, only: singlesize     
    implicit none  
    integer k,ncnt,i,j,cell_inc,m
    real(ikind) dump_depth,area,placement_vol,capacity,diff,excess,amount_placed,placement_vol_remaining
    real(ikind), allocatable :: bed(:),bed1(:)
    integer, allocatable :: cells(:)
    logical, allocatable :: flag(:)
    
    amount_placed = 0.0
    excess=0.0
    bed_change = 0.0
    
     !copy values into local arrays for easier coding   
    Ncnt = dredge_operations(k)%NumPlacementAreaCells(j)    
    allocate(cells(Ncnt),bed(Ncnt),bed1(Ncnt),flag(Ncnt))   
      do i=1,Ncnt
      cells(i) = dredge_operations(k)%PlacementAreaCells(j,i) 
      bed(i) = -zb(cells(i)) 
      enddo     
    
    !remove volume from cells based on selected approach
    select case(dredge_operations(k)%placement_approach(j))
        
    case(1) !place evenly accross bed, but do not exceed water depth limit
        
      call sort2(bed,cells,ncnt) !this sorting facilitates looking for depth limit exceedances
      bed1 = bed
      !write(*,*)'in placement'
      !write(*,*)k,j,dredge_operations(k)%placement_vol        
      !write(*,*)k,j,dredge_operations(k)%placement_area(j)

      !calculate capacity of area
      capacity = 0.0
      do i=1,ncnt
        diff = bed(i) - dredge_operations(k)%placement_limit(j,i)
        capacity = capacity + max(diff,0.0)*dx(cells(i))*dy(cells(i))
      enddo 
        
      !write(77,*)'capacity = ',capacity
        
      if(capacity .le. 0) then
        excess = placement_vol
        !dredge_operations(k)%placement_excess = dredge_operations(k)%placement_excess + dredge_operations(k)%placement_vol            
        !DredgeTS_Vars(k,2) = 0.0
        !DredgeTS_Vars(k,3) = placement_vol
        !DredgeTS_Vars(k,6) = dredge_operations(k)%placement_excess  
        !write(77,*)'zero cap, returning'
        return
      endif
        
      !there is capacity but not enough, then calculate what can be placed
      if(capacity > placement_vol) then  !capacity limited
        !excess = placement_vol  - capacity
        !dredge_operations(k)%placement_vol  = capacity
        !dredge_operations(k)%placement_excess = dredge_operations(k)%placement_excess + excess 
        !write(77,*)'cap is limiting palcement: ',excess
        !else  ! entire placement_vol can be placed in area
        !excess = 0.0
        capacity = placement_vol
        !write(77,*)'room for all:',dredge_operations(k)%placement_vol,capacity
      endif
        
      !work form shallowest cell to deepest cell to place material
      do i=1,ncnt
        !based on remaining cells calculate the thickness
        area = 0.0
        do m=i,ncnt
          area =area + dx(cells(m))*dy(cells(m))
        enddo
        dump_depth = capacity/area
        !write(77,*)'dumpdepth = ',dump_depth,capacity,area
        !see if depth is OK for ith cell
        if(bed(i) - dump_depth < dredge_operations(k)%placement_limit(j,i)) then  !need to truncate cells placement
          diff = bed(i) - dredge_operations(k)%placement_limit(j,i)
          capacity = capacity - max(diff,0.0)*dx(cells(i))*dy(cells(i))
          bed(i) = min(bed(i),dredge_operations(k)%placement_limit(j,i))
        else !all remaining cells can handle dump_depth, so place it
          do m=i,ncnt
            bed(m) = bed(m) - dump_depth    
          enddo
          goto 100
        endif
      enddo
        
100   continue        

      !put modified bed back into zb array
      do i=1,ncnt
        bed_change(i)  = -zb(cells(i)) - bed(i)
        zb(cells(i)) = -bed(i)
      enddo
        
      !for printing
      amount_placed=0.0
      do i=1,ncnt
        amount_placed=amount_placed + (bed1(i)-bed(i))*dx(cells(i))*dy(cells(i))
      enddo   
        
      !if(k == 1 .and. j == 1 .and. amount_placed < 30) then
      !  write(*,*)'amount_placed = ',amount_placed,placement_vol, dredge_operations(k)%DredgePlacementAllocation(j)
      !  stop
      !endif
        
      excess = placement_vol - amount_placed
      !DredgeTS_Vars(k,2) = sumvol
      !DredgeTS_Vars(k,3) = excess
      !DredgeTS_Vars(k,6) = dredge_operations(k)%placement_excess
        
    case(2) !place at cell closest to starting cell first, up to thikness, and do not exceed limits
            !cells ranked by proximity in dredge_init()
            
      bed1 = bed
      placement_vol_remaining = placement_vol
      !placement_vol = dredge_operations(k)%placement_vol
      cell_inc = dredge_operations(k)%placement_cell_inc(j)                   !added by Chris Reed
      do i=dredge_operations(k)%placement_cell_inc(j),ncnt  
        if(placement_vol_remaining > 0.0) then  !place material into cell 
          if(bed(i) > dredge_operations(k)%placement_limit(j,i)) then ! this cell still has capacity
            capacity = dx(cells(i))*dy(cells(i))*(bed(i) - dredge_operations(k)%placement_limit(j,i))
            if(capacity > placement_vol_remaining) then
              !zb(cells(i)) = zb(cells(i)) + placement_vol_remaining/(dx(cells(i))*dy(cells(i)))
              !bed(i) = bed(i) - placement_vol_remaining/(dx(cells(i))*dy(cells(i)))
              bed_change(i) = placement_vol_remaining/(dx(cells(i))*dy(cells(i)))
              zb(cells(i)) = zb(cells(i)) + bed_change(i)
              bed(i) = bed(i) - bed_change(i)              
              placement_vol_remaining = 0.0
              cell_inc = i               
            else
              !zb(cells(i)) = zb(cells(i)) + capacity/(dx(cells(i))*dy(cells(i)))
              !bed(i) = bed(i) - capacity/(dx(cells(i))*dy(cells(i)))
              bed_change(i) = capacity/(dx(cells(i))*dy(cells(i)))
              zb(cells(i)) = zb(cells(i)) + bed_change(i)
              bed(i) = bed(i) - bed_change(i)              
              placement_vol_remaining = placement_vol_remaining - capacity
              cell_inc = i+1
            endif
          endif
        else
          go to 200         
        endif  
      enddo
      !if we end up here then we cannot place all of the material and we need to put it into excess
      !dredge_operations(k)%placement_excess = dredge_operations(k)%placement_excess + placement_vol         
200   continue        
      dredge_operations(k)%placement_cell_inc(j)  = cell_inc 
      !excess = placement_vol
        
      !for printing
      amount_placed=0.0
      do i=1,ncnt
        amount_placed=amount_placed + (bed1(i)-bed(i))*dx(cells(i))*dy(cells(i))
      enddo   
      excess = placement_vol - amount_placed              
      !DredgeTS_Vars(k,2) = sumvol
      !DredgeTS_Vars(k,3) = placement_vol_init - sumvol
      !dredge_operations(k)%placement_excess = dredge_operations(k)%placement_excess + dredge_operations(k)%placement_vol - sumvol       
      !DredgeTS_Vars(k,6) = dredge_operations(k)%placement_excess
        
    case default
      write(*,*)'Approach specified for dredge material placement is not recognized'
      write(*,*)'Approach specified is: ',dredge_operations(k)%placement_approach
        
    end select

    if(.not. singlesize) call place_bed_sort(ncnt,cells)    

    return
    end subroutine place_it    
    
!**************************************************************
  subroutine dredge_implicit_update
! updates imlicit varaibles effected by depth changes due
! to dredge operations
! written by Chris Reed
!**************************************************************    
    use size_def
    use geo_def, only: idmap,zb,zbk,dzbx,dzby
    use flow_def
    use sed_def
    use comp_lib
    use comvarbl
    use der_def, only: nder,nlim,goa
    use der_lib, only: der_grad_eval
    use interp_lib, only: interp_scal_cell2face
    use prec_def 
    implicit none
    integer i
    real(ikind):: val
    
    !=== Bed-slopes ====--------------===========
    call der_grad_eval(goa,0,zb,dzbx,dzby) !Bed-slope

    !=== Update bed elevation at cell faces =====
    call interp_scal_cell2face(zb,0,zbk,dzbx,dzby)

    !=== Correct concentrations for depth changes =======
 !$OMP PARALLEL DO PRIVATE(i,val)
    do i=1,ncells
      val=h(i)                        !save old total water depth
      h(i)=max(hmin,p(i)/grav-zb(i))  !New total water depth
      !if(h(i)>2*hmin .and. val>2*hmin .and. abs(dzb(i))<0.25*hmin)then
      if(h(i)>2*hmin .and. val>2*hmin)then
        Ctk(i,:) = val*Ctk(i,:)/h(i)
      endif
    enddo
 !$OMP END PARALLEL DO

 !$OMP PARALLEL DO PRIVATE(i)
    do i=1,ncells
      dzb(i)=0.0
      dzbk(i,:)=0.0
    enddo
 !$OMP END PARALLEL DO


    !=== Apply boundary condition to zb ====
    call bndzb

    return
    end subroutine dredge_implicit_update
    
!**************************************************************
   subroutine WriteDredgeTS(k)
! moves sediment for any active dredging operation
! written by Chris Reed
!**************************************************************    
    use dredge_def
    use comvarbl, only: timehrs
    use geo_def, only: zb,zb0
    implicit none 
    integer :: k,j,m,NumPlacementAreas

51  FORMAT (f15.4,8(',',e15.7))
102 FORMAT (f15.4,1000(',',e15.7))

    NumPlacementAreas = dredge_operations(k)%NumPlacementAreas

    !MEB 04/24/2017 Switched from using the FORMAT backslash character to concatenate because of issues with GNU compiler.
    if(dredge_operations(k)%numplacementareas == 1) then
      write(DredgeUnit(k),102) timehrs,(DredgeTS_Vars(k,j),j=1,8)
    else
      write(DredgeUnit(k),102) timehrs,(DredgeTS_Vars(k,j),j=1,8), (dredge_operations(k)%placed_vol_by_area(m), dredge_operations(k)%total_placed_by_area(m), m=1,NumPlacementAreas)
    endif

    !Reset for next interval
    DredgeTS_Vars(k,:) = 0.0
    
    return
    end subroutine WriteDredgeTS
    
!**************************************************************
   subroutine Dredge_hardbottom
! checks to see if dredging has occured through hardbottom and adjusts
! the haarbottom depth to the new dredge depth
! written by Chris Reed
!**************************************************************    
    use prec_def
    use sed_def, only: poros,hardbed,nhard,idhard
    use geo_def, only: zb
    use sed_Def, only: hardbottom
    implicit none
    integer:: i,id
    real(ikind):: tdepth
    
    do i=1,nhard 
      id = idhard(i)
      tdepth = -zb(id)
      if(tdepth>hardbed(i)) then
        hardbed(i) = tdepth
      endif
    enddo
    
    return
    end subroutine Dredge_hardbottom
    
!*********************************************************************
! sorts x, y follows, x is real, y is integer
!*********************************************************************
 SUBROUTINE SORT2(X,Y,N) 
!     SUBROUTINE TO SORT TWO VECTORS  
!     SORTING IS DONE ACCORDING TO VALUE OF X
!     uses Shell algorithm
!     
      use prec_def 
      integer N
      real(ikind) X(N),TEMPX  
      INTEGER  H,Y(N),TEMPY,I,J
!
!   first find increment
      H=1
 10   H=3*H+1
      IF(H.LE.N) GO TO 10
!
!   main loop
!
 20   H=H/3     
!
!   sort pass at increment H
      DO 45 I=H+1,N
        TEMPX=X(I)
        TEMPY=Y(I)
        J=I
!       find proper place and make a hole there
 30     IF( X(J-H) > TEMPX ) THEN
          X(J)=X(J-H)
          Y(J)=Y(J-H)
          J=J-H
          IF(J>H) GO TO 30
        ENDIF
!       insert new element
        X(J)=TEMPX
        Y(J)=TEMPY
 45   CONTINUE  
!
      IF(H>1) GO TO 20      
      RETURN
    END   
    
 !*********************************************************************
! sorts x, y follows, x is real, y is integer
!*********************************************************************
 SUBROUTINE SORT2DP(X,Y,N) 
!     SUBROUTINE TO SORT TWO VECTORS  
!     SORTING IS DONE ACCORDING TO VALUE OF X
!     uses Shell algorithm
!       
      Integer N
      REAL*8  X(N)
      real TEMPX
      INTEGER  H,Y(N),TEMPY,I,J
!
!   first find increment
      H=1
 10   H=3*H+1
      IF(H.LE.N) GO TO 10
!
!   main loop
!
 20   H=H/3     
!
!   sort pass at increment H
      DO 45 I=H+1,N
        TEMPX=X(I)
        TEMPY=Y(I)
        J=I
!       find proper place and make a hole there
 30     IF( X(J-H) > TEMPX ) THEN
          X(J)=X(J-H)
          Y(J)=Y(J-H)
          J=J-H
          IF(J>H) GO TO 30
        ENDIF
!       insert new element
        X(J)=TEMPX
        Y(J)=TEMPY
 45   CONTINUE  
!
      IF(H>1) GO TO 20      
      RETURN
      END        
    
!****************************************************    
    subroutine dredge_eval()
!****************************************************
    use dredge_def
    use comvarbl, only: dtime
    implicit none
    
    dredge_time_lapse = dredge_time_lapse + dtime
    if(dredge_time_lapse >= dredge_interval) then
      call dredge_Op
      dredge_time_lapse = 0
    endif
    
    return
    end subroutine dredge_eval
    
!****************************************************    
    subroutine dredge_bed_sort(ncnt,cells)
!****************************************************    
    use size_def, only: ncells
    use sed_def
     use prec_def
    use dredge_def, only: bed_change,dredge_mat_gradation
    implicit none

    integer ncnt,lay_max,ival
    integer cells(ncnt) 
    integer i,id,k
    real depth,thickness,sum1(nsed),sum2
    real db_temp(nlay) 
      
    sum1 = 0
    sum2 = 0

    DO id=1,ncnt
      i = cells(id) 
      db_temp(:) = db(i,:)
      if(bed_change(id) .gt. db(i,1) ) then !change is larger than active layer thickness
        !add fractional content contribution for each layer up to bed_change depth
        depth = db(i,1)
        sum1(:) = db(i,1)*pbk(i,:,1)
        lay_max = nlay
        do k=2,nlay
          depth = depth + db(i,k)
          if(bed_change(id) .gt. depth) then
            sum1(:) = sum1(:) + db(i,k)*pbk(i,:,k)
          else
            thickness = bed_change(id) - (depth - db(i,k))
            sum1(:) = sum1(:) + thickness*pbk(i,:,k)
            lay_max = k
            go to 10
          endif
        enddo
10      continue          

        !reset gradation
        if(lay_max .eq. nlay) then  !dredged past lower layer so now all layers have last layer sorting
          do k=1,nlay-1
            pbk(i,:,k) = pbk(i,:,nlay)
          enddo
        else  !partial layer becomes active and lower layers shifted
          db(i,1) = thickness  
          pbk(i,:,1) = pbk(i,:,lay_max)
          do k=2,nlay
            ival = min(nlay,lay_max+k-1)
            db(i,k)  = db(i,ival)
            pbk(i,:,k) = pbk(i,:,ival)          
          enddo
        endif
        
        !adjust layer thickness if needed
        if(db(i,2)<=dbmax .and. db(i,2)>=dbmin)then !Second layer ok
          !db(i,3:nlay) = db1(i,3:nlay) 
        elseif(db(i,2)>dbmax)then !Second layer too thick, split into two              
          !First merge last two layers   
          pbk(i,:,nlay) = (db(i,nlay)*pbk(i,:,nlay) + & 
          db(i,nlay-1)*pbk(i,:,nlay-1))/(db(i,nlay) + db(i,nlay-1) )
          db(i,nlay) = db(i,nlay) + db(i,nlay-1) 
          !Move index of layers 4 to nlay-1 
          db(i,4:nlay-1) = db(i,3:nlay-2)
          pbk(i,:,4:nlay-1) = pbk(i,:,3:nlay-2)
          !Second and third layer compositions equal
          pbk(i,:,3) = pbk(i,:,2)                                        
        else !Second layer too thin
          !Merge second with third layer
          pbk(i,:,2) = (db(i,2)*pbk(i,:,2) + &
          db(i,3)*pbk(i,:,3))/(db(i,2) + db(i,3))
          db(i,2) = db(i,2) + db(i,3)    
          !add layer to bottom to make up for merged 2/3 layer
          if(db(i,nlay)>=dbmax)then !split last layer
            db(i,3:nlay-2) = db1(i,4:nlay-1)
            pbk(i,:,3:nlay-2) = pbk(i,:,4:nlay-1)
            db(i,nlay-1) = 0.5*db(i,nlay)
            db(i,nlay) = db(i,nlay-1)
            pbk(i,:,nlay-1) = pbk(i,:,nlay)
            !pbk(i,:,nlay) = pbk(i,:,nlay) stays the same
          else !add bottom layer
            db(i,3:nlay-1) = db(i,4:nlay)
            pbk(i,:,3:nlay-1) = pbk(i,:,4:nlay)
            db(i,nlay) = db(i,nlay)
            pbk(i,:,nlay) = pbk(i,:,nlay)
          endif
        endif    
      else !change is less than activer layer thickness
        !set dredge material gradation to that of active layer 
        sum1(:) = sum1(:) + bed_change(id)*pbk(i,:,1)        
        db(i,1) = db(i,1) - bed_change(id)
     
        !adjust layer thickness if needed     
        ! no need to adjust layers
      endif
     
      !write(*,"(2i4,8e10.3)")id,i,bed_change(id),db(i,1),pbk(i,:,1)
         
    ENDDO
     
    sum2 = sum(sum1(:))
    !write(*,"(2i4,5e10.3)")id,i,sum2,sum1(:)
    
    if (sum2 .ne. 0.0) then 
      dredge_mat_gradation(:) = sum1(:)/sum2  !normalized  
    endif
    !write(*,*)'dregmatgrad: ',dredge_mat_gradation(:)

    return
    end subroutine dredge_bed_sort
    
!****************************************************    
    subroutine place_bed_sort(ncnt,cells)
!****************************************************    
    use size_def, only: ncells
    use sed_def
     use prec_def
    use dredge_def, only: bed_change,dredge_mat_gradation      
    implicit none

    integer ncnt,i,id
    integer cells(ncnt),dbinit(ncells)
    real sum1(nsed),sum2  
      
    DO id=1,ncnt
      i = cells(id)
      dbinit(id) = db(i,1)
     
      !put placed material into upper (active) layer
      sum1(:) = db(i,1)*pbk(i,:,1) + bed_change(id)*dredge_mat_gradation(:)    
      sum2 = sum(sum1(:))
      pbk(i,:,1) = sum1(:)/sum2
      db(i,1) = db(i,1) + bed_change(id)
     
      !adjust layer thickness if needed  
      if(db(i,1)>dbmax)then !first layer too thick, split into two              
        !First merge last two layers   
        pbk(i,:,nlay) = (db(i,nlay)*pbk(i,:,nlay) + & 
          db(i,nlay-1)*pbk(i,:,nlay-1))/(db(i,nlay) + db(i,nlay-1) )
        db(i,nlay) = db(i,nlay) + db(i,nlay-1)                
        !Move index of layers 4 through nlay-1 
        db(i,3:nlay-1) = db1(i,2:nlay-2)
        pbk(i,:,3:nlay-1) = pbk(i,:,2:nlay-2)
        !Second and third layer compositions equal
        pbk(i,:,2) = pbk(i,:,1)                                           
      endif     
      
      if(db(i,2)>dbmax)then !Second layer too thick, split into two              
        !First merge last two layers   
        pbk(i,:,nlay) = (db(i,nlay)*pbk(i,:,nlay) + & 
          db(i,nlay-1)*pbk(i,:,nlay-1))/(db(i,nlay) + db(i,nlay-1) )
        db(i,nlay) = db(i,nlay) + db(i,nlay-1)                
        !Move index of layers 4 through nlay-1 
        db(i,4:nlay-1) = db1(i,3:nlay-2)
        pbk(i,:,4:nlay-1) = pbk(i,:,3:nlay-2)
        !Second and third layer compositions equal
        pbk(i,:,3) = pbk(i,:,2)                                           
      endif   
      
      !write(*,"(2i4,6e10.3)")id,i,bed_change(id),dbinit(id),db(i,1),pbk(i,:,1)
    ENDDO
     
    return
    end subroutine place_bed_sort   
    
    
!**************************************************************************
    subroutine sedpercentile_bedlayer(iper,klay,dper)
! Calculates sediment grain size percentiles from a fractional distribution    
! written by  Weiming Wu, NCCHE; Alex Sanchez, USACE-ERDC-CHL  
!**************************************************************************
    use size_def
    use sed_def, only: nsed,pbk,diam,diamlim,logdiamlim
    use prec_def
    implicit none
    !Input/Output
    integer,intent(in) :: iper,klay
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
        pbcum(ks+1) = pbcum(ks) + pbk(i,ks,klay)
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
    end subroutine sedpercentile_bedlayer
