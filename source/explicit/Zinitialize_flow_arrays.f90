!********************************************************************************    
  subroutine initialize_flow_arrays
!********************************************************************************    
#include "CMS_cpp.h"    
    use EXP_Global_def, only: linktodummies, qx, qxn, qy, qyn, etan, ue, ve, advectx, advecty, cdx, cdy, rhoprim, num_fg_m_cells
    use EXP_Global_def, only: advect, mixing, fuu, fuv, gvv, gvu, drydep, num_fg_a_cells, ncn, nce, ncs, ncw, fg_a_cells, fg_m_cells
    use size_def, only: ncells, ncellsd
    use flow_def, only: eta,iwet
    use met_def, only: tauwindx,tauwindy,pressatm
    use wave_flowgrid_def, only: wavestrx,wavestry
    use geo_def, only: zb,cell2cell 

    implicit none 
    !local variables
    integer i,icnt,ii,jj,ncn2,itag
    real totdepth
 
    !create the linktodummies array and populate
    !this array is used to quickly copy grid interior values to dummy cells
    icnt = 0.0
    do i=1,ncells
      if(cell2cell(1,i) .gt. ncells) icnt = icnt+1
      if(cell2cell(2,i) .gt. ncells) icnt = icnt+1    
      if(cell2cell(3,i) .gt. ncells) icnt = icnt+1
      if(cell2cell(4,i) .gt. ncells) icnt = icnt+1
    enddo
    allocate(linktodummies(icnt))
    icnt = 0.0
    do i=1,ncells
      if(cell2cell(1,i) .gt. ncells)then
        icnt = icnt+1
        linktodummies(icnt) = i
      endif
      if(cell2cell(2,i) .gt. ncells)then
        icnt = icnt+1
        linktodummies(icnt) = i
      endif   
      if(cell2cell(3,i) .gt. ncells)then
        icnt = icnt+1
        linktodummies(icnt) = i
      endif
      if(cell2cell(4,i) .gt. ncells)then
        icnt = icnt+1
        linktodummies(icnt) = i
      endif
    enddo
   
    allocate (qx(ncellsD),qxn(ncellsD),qy(ncellsD),qyn(ncellsD))
    allocate (etan(ncellsD),uE(ncellsD),vE(ncellsD))
    allocate (ADVECTX(NcellsD),ADVECTY(NcellsD)) 
    allocate (cdx(ncellsD),cdy(ncellsD))
    allocate(RHOPrim(ncellsD))  

    DO I=1,NCELLSD
      QX(I) = 0.0
      QY(I) = 0.0
      QXN(I) = 0.0
      QYN(I) = 0.0
      ETA(I) = 0.0
      ETAN(I) = 0.0
      UE(I) = 0.0
      VE(I) = 0.0
      CDX(i) = 0.005
      CDY(i) = 0.005
      ADVECTX(I) = 0.0
      ADVECTY(I) = 0.0
      iWET(I) = 1
      RHOPrim(I) = 1.0          !CR - 01/21/2009
    ENDDO    
      
    eta = 0.0  !9.47
    etan = eta

    if(ADVECT .or. mixing) then
      allocate (Fuu(ncellsD),Fuv(NcellsD),Gvv(NcellsD),Gvu(NcellsD))
      do i=1,ncellsD
        Fuu(i) = 0.0
        Fuv(i) = 0.0
        Gvv(i) = 0.0
        Gvu(i) = 0.0
      enddo
    endif
      
    !modify inital eta array for wetting and drying
    do i=1,ncells
      totdepth = eta(i)-zb(i)
      if(totdepth.le.(0.9*drydep)) then                                     
        eta(i) = (0.9d0*drydep) + zb(i)
        etan(i) = eta(i)
      endif
    enddo
    !copy changes to dummy cells
    ii=0
    do i=ncells+1,ncellsD
      ii=ii+1
      jj = linktodummies(ii)
      eta(i) = eta(jj)
      etan(i) = etan(jj)
      zb(i) = zb(jj)       
    enddo    
        
    !MAKE LIST OF DUMMY CELLS THAT NEED FUV AND GVU CALCULATED
    !FOR ADVECTIVE TERMS    
    num_FG_A_cells = 0
    if(ADVECT) then  !see if list FG cells is required
      do i=1,ncells
        ncn = cell2cell(1,i)
        if(ncn.gt.ncells) then
          ncw = cell2cell(4,i)
          if(ncw.le.ncells) then
            ncn2 = cell2cell(1,ncw)
            if(ncn2.le.ncells) then
              nce = cell2cell(2,ncn2)
              if(nce.gt.ncells) then
                num_FG_A_cells = num_FG_A_cells + 1
              endif
            endif 
          endif 
        endif
      enddo  
    endif 
#ifdef DEBUG    
    write(DGUNIT,*)'number of special advective FG cells is: ',num_FG_A_cells
#endif
    if(num_FG_A_cells .gt. 0) then          !make list of dummy cells that need Fuv and Guv calcs
      allocate(FG_A_cells(num_FG_A_cells,4))
      num_FG_A_cells = 0
      do i=1,ncells
        ncn = cell2cell(1,i)
        if(ncn.gt.ncells) then
          ncw = cell2cell(4,i)
          if(ncw.le.ncells) then
            ncn2 = cell2cell(1,ncw)
            if(ncn2.le.ncells) then
              nce = cell2cell(3,ncn2)
              if(nce.gt.ncells) then
                num_FG_A_cells = num_FG_A_cells + 1
                FG_A_cells(num_FG_A_cells,1) = ncn
                FG_A_cells(num_FG_A_cells,2) = nce
                FG_A_cells(num_FG_A_cells,3) = i
                FG_A_cells(num_FG_A_cells,4) = ncn2
              endif
            endif 
          endif 
        endif
      enddo  
    endif     

    !MAKE LIST OF REGULAR CELLS THAT NEED FUV AND GVU SET TO ZERO - REQUIRED FOR
    !SYMMETRY AND TO PREVENT MOMENTUM DIFFUSION ACCROSS BOUNDARIES
    num_FG_M_cells = 0
    do i=1,ncells
      itag = 0
      ncs = cell2cell(3,i)
      if(ncs.gt.ncells) itag = 1
      if(cell2cell(4,ncs) .gt. ncells) itag = 1
      ncw = cell2cell(4,i)
      if(ncw.gt.ncells) itag = 1
      if(cell2cell(3,ncw) .gt. ncells) itag = 1
      if(itag.eq.1) num_FG_M_cells = num_FG_M_cells + 1
    enddo
#ifdef DEBUG
    write(DGUNIT,*)'number of special mixing FG cells is: ',num_FG_M_cells
#endif
    allocate(Fg_m_cells(num_FG_M_cells))
    if(num_FG_M_cells.gt.0) then
      num_FG_M_cells = 0
      do i=1,ncells
        itag = 0
        ncs = cell2cell(3,i)
        if(ncs.gt.ncells) itag = 1
        if(cell2cell(4,ncs) .gt. ncells) itag = 1
        ncw = cell2cell(4,i)
        if(ncw.gt.ncells) itag = 1
        if(cell2cell(3,ncw) .gt. ncells) itag = 1
        if(itag.eq.1) then
          num_FG_M_cells = num_FG_M_cells + 1
          FG_M_cells(num_Fg_M_cells) = i
        endif
      enddo    
    endif
     
    !These arrays are from the IMPLICIT code, and may or may not be dimensioned,
    !however they are needed for the explicit flowcode even if wind or waves are 
    !not simulated
    if ( .not. allocated(pressatm)) then
      allocate(pressatm(ncellsD))
      pressatm = 0.0
    endif
    if ( .not. allocated(wavestrx)) then
      allocate(wavestrx(ncellsD))
      wavestrx = 0.0
    endif
    if ( .not. allocated(wavestry)) then
      allocate(wavestry(ncellsD))
      wavestry = 0.0
    endif      
    if ( .not. allocated(tauwindx)) then
      allocate(tauwindx(ncellsD))
      tauwindx = 0.0
    endif
    if ( .not. allocated(tauwindy)) then
      allocate(tauwindy(ncellsD))
      tauwindy = 0.0
    endif   

  return
  end subroutine        