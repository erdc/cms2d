!********************************************************************************
!           routine to set the active array for explicit solver
!           active(*,1) = false   no u-momentum calculation
!           active(*,2) = false   no v-momentum calculation
!           active(*,3) = false   no mass balance calculation 
!
!********************************************************************************    
    SUBROUTINE Initialize_ActivityArray_tel
    use EXP_Global_def
    USE EXP_bndcond_def 
    use geo_def, only: dx,dy,cell2cell
    use size_def     
    use flow_def
    use comvarbl
    use out_def, only: goutfile
    use bnd_def
    use sal_def
      
    IMPLICIT NONE
    !LOCAL VARIABLES
    INTEGER I,J,JJ,ID_CELL
      
    
    allocate (active(ncellsD,3))
    do i=1,ncells
      !active(i,1) = .true.
      !active(i,2) = .true.
      active(i,3) = .true.
    enddo    
    do i=ncells+1,ncellsD
      !active(i,1) = .false.
      !active(i,2) = .false.
      active(i,3) = .false.
    enddo     
      
    !set activity for wall cells on south and west ends of grid
    !do i=1,ncells
    !  if(cell2cell(4,i).gt.ncells) active(i,1) = .false.
    !  if(cell2cell(3,i).gt.ncells) active(i,2) = .false.
    !enddo      

    if(nHstr .gt. 0) then
      !set activity based on wse bc cells
      do i=1,nHstr
        do j=1,H_Str(I)%NCELLS
          id_cell = H_Str(I)%Cells(j)
          !no h equation for WSE cells
          active(id_cell,3)=.false.
          !no u calc for outside (west) edge of WSE cell
          !if(cell2cell(4,id_cell).gt.NCells) active(id_cell,1)=.false.
          !no v calc for outside (south) edge of WSE cell 
          !if(cell2cell(3,id_cell).gt.NCells) active(id_cell,2)=.false.
          !jj=max(j-1,1)
          !if(cell2cell(3,id_cell) .eq. H_Str(i)%cells(jj) ) active(id_cell,2)=.false. 
          !jj=min(j+1,H_Str(I)%ncells) 
          !if(cell2cell(3,id_cell) .eq. H_Str(i)%cells(jj))  active(id_cell,2)=.false. 
          !jj=max(j-1,1) 
          !if(cell2cell(4,id_cell) .eq. H_Str(i)%cells(jj))  active(id_cell,1)=.false. 
          !jj=min(j+1,H_Str(I)%ncells) 
          !if(cell2cell(4id_cell) .eq. H_Str(i)%cells(jj))  active(id_cell,1)=.false. 
        enddo
      enddo
    Endif !(H_single)	

    if(nTHstr .gt. 0) then
      i=1
      !set activity based on wse bc cells
      do j=1,TH_str(i)%ncells
        id_cell = TH_Str(I)%Cells(j)
        !no h equation for WSE cells
        active(id_cell,3)=.false.
        !no u calc for outside (west) edge of WSE cell
        !if(cell2cell(4,id_cell).gt.NCells) active(id_cell,1)=.false.
        !!no v calc for outside (south) edge of WSE cell 
        !if(cell2cell(3,id_cell).gt.NCells) active(id_cell,2)=.false.
        !jj=max(j-1,1)
        !if(cell2cell(3,id_cell) .eq. TH_Str(I)%Cells(jj)) active(id_cell,2)=.false.
        !jj=min(j+1,TH_str(i)%ncells)
        !if(cell2cell(3,id_cell) .eq. TH_Str(I)%Cells(jj)) active(id_cell,2)=.false.
        !jj=max(j-1,1)
        !if(cell2cell(4,id_cell) .eq. TH_Str(I)%Cells(jj)) active(id_cell,1)=.false.
        !jj=min(j+1,TH_str(i)%ncells)
        !if(cell2cell(4,id_cell) .eq. TH_Str(I)%Cells(jj)) active(id_cell,1)=.false.
      enddo
    Endif !(H_tides)	
		
    if(nMHstr .gt. 0) then
      !set activity based on wse bc cells
      do i=1,nMHstr
        do j=1,MH_Str(I)%NCELLS
          id_cell = MH_Str(I)%Cells(j)
          !no h equation for WSE cells
          active(id_cell,3)=.false.
          !no u calc for outside (west) edge of WSE cell
          !if(cell2cell(4,id_cell).gt.NCells) active(id_cell,1)=.false.
          !no v calc for outside (south) edge of WSE cell 
          !if(cell2cell(3,id_cell).gt.NCells) active(id_cell,2)=.false.
          !jj=max(j-1,1)
          !if(cell2cell(3,id_cell) .eq. MH_Str(i)%cells(jj)) active(id_cell,2)=.false.
          !jj=min(j+1,MH_Str(I)%ncells)
          !if(cell2cell(3,id_cell) .eq. MH_Str(i)%cells(jj)) active(id_cell,2)=.false.
          !jj=max(j-1,1)
          !if(cell2cell(4,id_cell) .eq. MH_Str(i)%cells(jj)) active(id_cell,1)=.false.
          !jj=min(j+1,MH_Str(I)%ncells)
          !if(cell2cell(4,id_cell) .eq. MH_Str(i)%cells(jj)) active(id_cell,1)=.false.
        enddo
      enddo
    Endif !(H_multi)
	
    if(nMHVstr .gt. 0) then
      !set activity based on wse bc cells
      do i=1,nMHVstr
        do j=1,MHV_str(I)%NCELLS
          id_cell = MHV_str(I)%Cells(j)
          !no h equation for WSE cells
          active(id_cell,3)=.false.
          !no u calc for outside (west) edge of WSE cell
          !if(cell2cell(4,id_cell).gt.NCells) active(id_cell,1)=.false.
          !no v calc for outside (south) edge of WSE cell 
          !if(cell2cell(3,id_cell).gt.NCells) active(id_cell,2)=.false.
          !jj=max(j-1,1)
          !if(cell2cell(3,id_cell) .eq. MHV_str(i)%cells(jj)) active(id_cell,2)=.false.
          !jj=min(j+1,MHV_str(I)%ncells)
          !if(cell2cell(3,id_cell) .eq. MHV_str(i)%cells(jj)) active(id_cell,2)=.false.
          !jj=max(j-1,1)
          !if(cell2cell(4,id_cell) .eq. MHV_str(i)%cells(jj)) active(id_cell,1)=.false.
          !jj=min(j+1,MHV_str(I)%ncells)
          !if(cell2cell(4,id_cell) .eq. MHV_str(i)%cells(jj)) active(id_cell,1)=.false.
        enddo
      enddo
    endif
      
    RETURN
    END SUBROUTINE 
      
