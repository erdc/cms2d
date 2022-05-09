    subroutine initialize_flow_arrays_tel
#include "CMS_cpp.h"    
    use EXP_Global_def
    USE EXP_transport_def 
    USE EXP_bndcond_def
    use size_def    
    use flow_def, only: eta,iwet
    use met_def, only: tauwindx,tauwindy   
    use wave_flowgrid_def, only: wavestrx,wavestry
    use geo_def, only: zb,cell2cell,idirface,dx,dy,mapid,icol,irow
    use met_Def, only: pressatm
    use exp_telescoping
    use sal_def, only: saltrans,sal        
    use bnd_def
      
    implicit none 
    !local variables
    integer i,icnt,ii,jj,ncn2,itag,k,j,id1,id2,IDopt1,IDopt2,L
    integer id3,id4,IT1,IT2,IB1,IB2,IR1,IR2,IL1,IL2
    integer  n_cnt,e_cnt,s_cnt,w_cnt,nfaces
    real totdepth
    logical id1dummy,id2dummy,id3dummy,id4dummy
    logical reg1,reg2,reg3,reg4,reg5,reg6
 
    !create the linktodummies array and populate
    !this array is used to quickly copy grid interior values to dummy cells
    icnt = 0.0
    do i=1,ncells
      if(cell2cell(1,i) .gt. ncells) icnt = icnt+1
      if(cell2cell(2,i) .gt. ncells) icnt = icnt+1    
      if(cell2cell(3,i) .gt. ncells) icnt = icnt+1
      if(cell2cell(4,i) .gt. ncells) icnt = icnt+1
      if(cell2cell(5,i) .gt. ncells) icnt = icnt+1
      if(cell2cell(6,i) .gt. ncells) icnt = icnt+1    
    enddo
    num_linktodummies = icnt
    allocate(linktodummiesTel(2,icnt))
    icnt = 0.0
    do i=1,ncells
      if(cell2cell(1,i) .gt. ncells)then
        icnt = icnt+1
        linktodummiesTel(1,icnt) = i
        linktodummiesTel(2,icnt) = cell2cell(1,i)     
      endif
      if(cell2cell(2,i) .gt. ncells)then
        icnt = icnt+1
        linktodummiesTel(1,icnt) = i
        linktodummiesTel(2,icnt) = cell2cell(2,i)       
      endif   
      if(cell2cell(3,i) .gt. ncells)then
        icnt = icnt+1
        linktodummiesTel(1,icnt) = i
        linktodummiesTel(2,icnt) = cell2cell(3,i)       
      endif
      if(cell2cell(4,i) .gt. ncells)then
        icnt = icnt+1
        linktodummiesTel(1,icnt) = i
        linktodummiesTel(2,icnt) = cell2cell(4,i)        
      endif
      if(cell2cell(5,i) .gt. ncells)then
        icnt = icnt+1
        linktodummiesTel(1,icnt) = i
        linktodummiesTel(2,icnt) = cell2cell(5,i)        
      endif
      if(cell2cell(6,i) .gt. ncells)then
        icnt = icnt+1
        linktodummiesTel(1,icnt) = i
        linktodummiesTel(2,icnt) = cell2cell(6,i)        
      endif    
    enddo
    
    ! copy dx,dy to dummy cells
    do i=1,num_linktodummies
      id1 = linktodummiesTel(1,i)
      id2 = linktodummiesTel(2,i)
      dx(id2) = dx(id1)
      dy(id2) = dy(id1)    
    enddo        

    !put cell mapping info into cellmap array
    !this array combines implicit code cell2cell and idirface arrays
    allocate(cellmap(8,ncells))
    do i=1,ncells 
      k=0
      n_cnt = 0
      e_cnt = 0
      s_cnt = 0
      w_cnt = 0
     
      do j=1,nmaxfaces
        select case(idirface(j,i))
        case(1)
          n_cnt=N_cnt+1
        case(2)
          e_cnt=e_cnt+1
        case(3)
          s_cnt=s_cnt+1
        case(4)
          w_cnt=w_cnt+1
        end select
      enddo
    
      if(n_cnt .eq. 1) then
        k=k+1
        cellmap(1,i) = cell2cell(k,i)
        cellmap(2,i) = 0   
      elseif(n_cnt .eq. 2) then
        k=k+1
        cellmap(1,i) = cell2cell(k,i) 
        k=k+1
        cellmap(2,i) = cell2cell(k,i)   
      endif
    
      if(e_cnt .eq. 1) then
        k=k+1
        cellmap(3,i) = cell2cell(k,i)
        cellmap(4,i) = 0   
      elseif(e_cnt .eq. 2) then
        k=k+1
        cellmap(3,i) = cell2cell(k,i) 
        k=k+1
        cellmap(4,i) = cell2cell(k,i)   
      endif    
    
      if(s_cnt .eq. 1) then
        k=k+1
        cellmap(5,i) = cell2cell(k,i)
        cellmap(6,i) = 0   
      elseif(s_cnt .eq. 2) then
        k=k+1
        cellmap(5,i) = cell2cell(k,i) 
        k=k+1
        cellmap(6,i) = cell2cell(k,i)   
      endif   
    
      if(w_cnt .eq. 1) then
        k=k+1
        cellmap(7,i) = cell2cell(k,i)
        cellmap(8,i) = 0   
      elseif(w_cnt .eq. 2) then
        k=k+1
        cellmap(7,i) = cell2cell(k,i) 
        k=k+1
        cellmap(8,i) = cell2cell(k,i)   
      endif    
    enddo
    
    !when idntify faces below, the cell faces will be identified
    !for use int he mawss continuity equations solution
    allocate(cellfaces(8,ncells))
    cellfaces=0
      
    !determine the face array for x 
    nfaces =0
    do i = 1,ncells
      if(cellmap(7,i) .gt. 0) nfaces = nfaces + 1
      if(cellmap(8,i) .gt. 0) nfaces = nfaces + 1
      !check for dummy cell
      if(cellmap(3,i) .gt. ncells) nfaces = nfaces + 1
    enddo
    
    !allocate (xFace(0:nfaces))
    allocate (xface_CadvF(8,0:nfaces),xface_cells(4,0:nfaces),xface_advF(4,0:nfaces))
    allocate (xface_q(0:nfaces),xface_qn(0:nfaces),xface_Length(0:nfaces),xface_gcoef(4,0:nfaces),xface_vel(0:nfaces),xface_advdif_I(2,0:nfaces),xface_advdif_C(0:nfaces))
    allocate (xface_wall(0:nfaces),xface_grad1(0:nfaces),xface_basic_orientation(0:nfaces))   
    allocate(xface_wet(0:nfaces))
    allocate(xface_flux(0:nfaces))
    
    nfaces =0
    xFace_q(0)= 0    !this and the next line facilitate mass conservation equation solution
    xFace_length(0)=0 
    xface_advdif_I(1,0) = 0  
    xface_advdif_I(2,0) = 0      
    xface_vel(0)=0
    xface_advdif_C(0) = 0    
    xface_CadvF=0      !MEB added 4/19/2016
    xface_wet = .true.
    xface_flux = 0.0
    do i = 1,ncells
      if(cellmap(7,i) .gt. 0) then
        nfaces = nfaces + 1
        xFace_cells(2,nfaces) = cellmap(7,i)
        xFace_cells(1,nfaces) = i
        xface_wall(nfaces) = .false.
        xface_Length(nfaces) = min(dy(i),dy(cellmap(7,i)))
        cellfaces(7,i) = nfaces     
        if(cellmap(7,i).gt.ncells) then
          xface_wall(nfaces) = .true.
        else  
          if(cellmap(3,cellmap(7,i)) .eq. i) cellfaces(3,cellmap(7,i)) = nfaces  
          if(cellmap(4,cellmap(7,i)) .eq. i) cellfaces(4,cellmap(7,i)) = nfaces      
        endif
      endif
      if(cellmap(8,i) .gt. 0) then
        nfaces = nfaces + 1
        xFace_cells(2,nfaces) = cellmap(8,i)
        xFace_cells(1,nfaces) = i
        xface_wall(nfaces) = .false.
        xface_Length(nfaces) = min(dy(i),dy(cellmap(8,i)))  
        cellfaces(8,i) = nfaces     
        if(cellmap(8,i).gt.ncells) then
          xface_wall(nfaces) = .true.   
        else
          cellfaces(3,cellmap(8,i)) = nfaces     
        endif
      endif
      !check for dummy cell to right and if so, add face for dummy cell
      if(cellmap(3,i) .gt. ncells) then
        nfaces = nfaces + 1
        xFace_cells(2,nfaces) = i
        xFace_cells(1,nfaces) = cellmap(3,i)
        xface_wall(nfaces) = .true.
        cellfaces(3,i) = nfaces  
        xface_Length(nfaces) = min(dy(i),dy(cellmap(3,i)))          
      endif
    enddo    
    numxFaces=nfaces
    
    !determine the face array for y 
    nfaces =0
    do i = 1,ncells
      if(cellmap(5,i) .gt. 0) nfaces = nfaces + 1
      if(cellmap(6,i) .gt. 0) nfaces = nfaces + 1
      !check for dummy cell
      if(cellmap(1,i) .gt. ncells) nfaces = nfaces + 1
    enddo
    
    !allocate (yFace(0:nfaces))  
    allocate (yface_CadvF(8,0:nfaces),yface_cells(4,0:nfaces),yface_advF(4,0:nfaces))
    allocate (yface_q(0:nfaces),yface_qn(0:nfaces),yface_Length(0:nfaces),yface_gcoef(4,0:nfaces),yface_vel(0:nfaces),yface_advdif_I(2,0:nfaces),yface_advdif_C(0:nfaces))
    allocate (yface_wall(0:nfaces),yface_grad1(0:nfaces),yface_basic_orientation(0:nfaces))     
    allocate(yface_wet(0:nfaces))
    allocate(yface_flux(0:nfaces))
    nfaces =0
    yFace_q(0) = 0  !this and the next line facilitate mass conservation equation solution
    yFace_length(0)=0
    yface_vel(0) =0
    yface_advdif_I(1,0) = 0
    yface_advdif_I(2,0) = 0   
    yface_advdif_C(0) = 0
    yface_CadvF=0      !MEB added 4/19/2016    
    yface_wet = .true.
    yface_flux = 0.0
    do i = 1,ncells
      if(cellmap(5,i) .gt. 0) then
        nfaces = nfaces + 1
        yFace_cells(2,nfaces) = cellmap(5,i)
        yFace_cells(1,nfaces) = i
        yface_wall(nfaces) = .false.
        yface_Length(nfaces) = min(dx(i),dx(cellmap(5,i)))  
        cellfaces(5,i) = nfaces     
        if(cellmap(5,i).gt.ncells) then
          yface_wall(nfaces) = .true.
        else
          if(cellmap(1,cellmap(5,i)) .eq. i) cellfaces(1,cellmap(5,i)) = nfaces  
          if(cellmap(2,cellmap(5,i)) .eq. i) cellfaces(2,cellmap(5,i)) = nfaces            
        endif
      endif
      if(cellmap(6,i) .gt. 0) then
        nfaces = nfaces + 1
        yFace_cells(2,nfaces) = cellmap(6,i)
        yFace_cells(1,nfaces) = i
        yface_wall(nfaces) = .false.
        yface_Length(nfaces) = min(dx(i),dx(cellmap(6,i)))  
        cellfaces(6,i) = nfaces     
        if(cellmap(6,i).gt.ncells) then
          yface_wall(nfaces) = .true.   
        else
          cellfaces(1,cellmap(6,i)) = nfaces   
        endif
      endif
      !check for dummy cell
      if(cellmap(1,i) .gt. ncells) then
        nfaces = nfaces + 1
        yFace_cells(2,nfaces) = i
        yFace_cells(1,nfaces) = cellmap(1,i)
        yface_wall(nfaces) = .true.
        yface_Length(nfaces) = min(dx(i),dx(cellmap(1,i)))      
        cellfaces(1,i) = nfaces      
      endif
    enddo     
    numyFaces=nfaces   
    
    !allocate coreolis velocity arrays
    allocate(QXc(numyfaces),QYc(numxFaces))
    
    !get template for inline u advection terms  !!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
    do I = 1,numXfaces
      xface_advF(:,i) = 0
      if(.not. xface_wall(i)) then
        !identify upstream face(s)           
        ID1=xface_cells(1,i)
        if(cellfaces(4,id1).gt.0) then
          xface_advF(1,i) = cellfaces(3,id1)
          xface_advF(2,i) = cellfaces(4,id1)   
        else
          xface_advF(1,i) = cellfaces(3,id1)  
          xface_advF(2,i) = 0          
        endif
        !identify downstream face(s)
        ID2=xface_cells(2,i)
        if(cellfaces(8,id2).gt.0) then
          xface_advF(3,i) = cellfaces(7,id2)
          xface_advF(4,i) = cellfaces(8,id2)   
        else
          xface_advF(3,i) = cellfaces(7,id2)  
          xface_advF(4,i) = 0           
        endif       
      elseif( xface_cells(1,i) .le. ncells) then  
        !identify upstream face(s)           
        ID1=xface_cells(1,i)
        if(cellfaces(4,id1).gt.0) then
          xface_advF(1,i) = cellfaces(3,id1)
          xface_advF(2,i) = cellfaces(4,id1)   
        else
          xface_advF(1,i) = cellfaces(3,id1)  
          xface_advF(2,i) = 0         
        endif
      endif !not a wall
    enddo
    
    !get template for inline v advection terms !!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do I = 1,numYfaces
      if(.not. yface_wall(i)) then
        !identify upstream face(s)   
        ID1=yface_cells(1,i)
        if(cellfaces(2,id1).gt.0) then
          yface_advF(1,i) = cellfaces(1,id1)
          yface_advF(2,i) = cellfaces(2,id1)   
        else
          yface_advF(1,i) = cellfaces(1,id1)  
          yface_advF(2,i) = 0           
        endif
        !identify downstream face(s)
        ID2=yface_cells(2,i)          
        if(cellfaces(6,id2).gt.0) then
          yface_advF(3,i) = cellfaces(5,id2)
          yface_advF(4,i) = cellfaces(6,id2)   
        else
          yface_advF(3,i) = cellfaces(5,id2)  
          yface_advF(4,i) = 0         
        endif 
      elseif (yface_cells(1,i) .le. ncells) then
        !identify upstream face(s)           
        ID1=yface_cells(1,i)
        if(cellfaces(2,id1).gt.0) then
          yface_advF(1,i) = cellfaces(1,id1)
          yface_advF(2,i) = cellfaces(2,id1)   
        else
          yface_advF(1,i) = cellfaces(1,id1)  
          yface_advF(2,i) = 0         
        endif         
      endif
    enddo
     
     
    !open(unit=670,file='Xadf.txt')
    !do i=1,numXfaces
    !  write(670,"(5i12)") i,(xface_advF(k,i),K=1,4)
    !enddo
    !open(unit=680,file='Yadf.txt')     
    !do i=1,numYfaces
    !  write(680,"(5i12)") i,(yface_advF(k,i),K=1,4)
    !enddo  
     
     
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    !set up indexes for duv/dy (x-momentum equation)
    do i=1,numXfaces
      ID1 = xface_cells(1,i)
      ID2 = xface_cells(2,i)
      IF( .not.  xface_wall(i)) then  
        IF(dx(ID1) .gt. 1.1*dx(ID2) ) THEN  !ID1 > ID2
          IF(cellfaces(3,ID2) .eq. cellfaces(8,ID1)) THEN  !it is a upper face
            !horizontal faces    
            xface_CadvF(1,i) = cellfaces(1,ID1)
            xface_CadvF(7,i) = cellfaces(1,ID1)
            if(cellfaces(2,ID2) .gt. 0) then
              xface_CadvF(2,i) = cellfaces(2,ID2)
              xface_CadvF(8,i) = cellfaces(2,ID2)      
            else
              xface_CadvF(2,i) = cellfaces(1,ID2) 
              xface_CadvF(8,i) = cellfaces(1,ID2)        
            endif            
            !vertical face at top
            IT2 = cellmap(1,ID2)
            if(cellmap(2,ID2) .gt. 0) IT2 = cellmap(2,ID2)
            if(it2.gt.ncells) then
              xface_CadvF(5,i) = i
              xface_CadvF(3,i) = i                
            else
              if(cellfaces(4,IT2) .gt. 0) then
                xface_CadvF(5,i) = cellfaces(4,IT2)
                xface_CadvF(3,i) = cellfaces(4,IT2) 
              else
                xface_CadvF(5,i) = cellfaces(3,IT2)   
                xface_CadvF(3,i) = cellfaces(3,IT2)              
              endif
            endif
            !vertical face at bottom   
            xface_CadvF(4,i) = cellfaces(7,ID1)   
            xface_CadvF(6,i) = cellfaces(7,ID1)              
          ELSE !it is a lower face
            !horizontal faces    
            xface_CadvF(1,i) = cellfaces(1,ID1)
            xface_CadvF(7,i) = cellfaces(5,ID1)
            if(cellfaces(2,ID2) .gt. 0) then
              xface_CadvF(2,i) = cellfaces(2,ID2)
              xface_CadvF(8,i) = cellfaces(2,ID2)      
            else
              xface_CadvF(2,i) = cellfaces(1,ID2) 
              xface_CadvF(8,i) = cellfaces(1,ID2)        
            endif 
            !vertical face at top         
            xface_CadvF(5,i) = cellfaces(8,ID1)
            xface_CadvF(3,i) = cellfaces(8,ID1)          
            !vertical face at bottom          
            IB2 = cellmap(5,ID2)
            if(ib2.gt.ncells) then
              xface_CadvF(6,i) = i
              xface_CadvF(4,i) = i                
            else          
              xface_CadvF(6,i) = cellfaces(3,IB2)   
              xface_CadvF(4,i) = cellfaces(3,IB2) 
           endif
         ENDIF !end of upper/lower face 
       ELSEIF(dx(ID2) .gt. 1.1*dx(ID1) ) THEN  !ID2 > ID1
         IF(cellfaces(3,ID2) .eq. cellfaces(7,ID1)) THEN  !upper face
           !horizontal faces    
           xface_CadvF(1,i) = cellfaces(1,ID1)
           xface_CadvF(7,i) = cellfaces(1,ID1)
          if(cellfaces(2,ID2) .gt. 0) then
            xface_CadvF(2,i) = cellfaces(2,ID2)
            xface_CadvF(8,i) = cellfaces(2,ID2)      
          else
            xface_CadvF(2,i) = cellfaces(1,ID2) 
            xface_CadvF(8,i) = cellfaces(1,ID2)        
          endif            
          !vertical face at top
          IT1 = cellmap(1,ID1)
          if(it1.gt.ncells) then
            xface_CadvF(5,i) = i
            xface_CadvF(3,i) = i                
          else         
            xface_CadvF(5,i) = cellfaces(7,IT1)   
            xface_CadvF(3,i) = cellfaces(7,IT1)
          endif
          !vertical face at bottom   
          xface_CadvF(4,i) = cellfaces(4,ID2)   
          xface_CadvF(6,i) = cellfaces(4,ID2)              
        ELSE !lower face
          !horizontal faces    
          xface_CadvF(1,i) = cellfaces(1,ID1)
          xface_CadvF(7,i) = cellfaces(1,ID1)
          if(cellfaces(2,ID2) .gt. 0) then
            xface_CadvF(2,i) = cellfaces(2,ID2)
          else
            xface_CadvF(2,i) = cellfaces(1,ID2)              
          endif
          xface_CadvF(8,i) = cellfaces(5,ID2)           
          !vertical face at top         
          xface_CadvF(5,i) = cellfaces(3,ID2)
          xface_CadvF(3,i) = cellfaces(3,ID2)          
          !vertical face at bottom          
          IB2 = cellmap(5,ID2)
          if(ib2.gt.ncells) then
            xface_CadvF(6,i) = i
            xface_CadvF(4,i) = i                
          else          
            xface_CadvF(6,i) = cellfaces(3,IB2)   
            xface_CadvF(4,i) = cellfaces(3,IB2) 
          endif
        ENDIF !end of upper/lower face           
      ELSE  !they are the same size
        !horizontal faces    
        xface_CadvF(1,i) = cellfaces(1,ID1)
        xface_CadvF(7,i) = cellfaces(1,ID1)
        if(cellfaces(2,ID2) .gt. 0) then
          xface_CadvF(2,i) = cellfaces(2,ID2)
          xface_CadvF(8,i) = cellfaces(2,ID2)      
        else
          xface_CadvF(2,i) = cellfaces(1,ID2) 
          xface_CadvF(8,i) = cellfaces(1,ID2)        
        endif 
        !vertical face at top
        IT1 = cellmap(1,ID1)
        IT2 = cellmap(1,ID2)      
        if(cellmap(2,ID2) .gt. 0) IT2 = cellmap(2,ID2)
          if(it1.gt.ncells.or.it2.gt.ncells) then
            xface_CadvF(5,i) = i   
            xface_CadvF(3,i) = i      
          else         
            if(IT1 .eq. IT2) then
              xface_CadvF(5,i) = cellfaces(7,IT2)   
              xface_CadvF(3,i) = cellfaces(3,IT1)
              if(cellfaces(4,IT1) .gt. 0) xface_CadvF(3,i) = cellfaces(4,IT1)
            else
              xface_CadvF(5,i) = cellfaces(7,IT1)   
              xface_CadvF(3,i) = cellfaces(7,IT1)                
            endif
          endif !if t1,it2 > ncells
          !vertical face at bottom
          IB2 = cellmap(5,ID2)
          IB1 = cellmap(5,ID1)
          if(ib1.gt.ncells.or.ib2.gt.ncells) then
            xface_CadvF(6,i) = i   
            xface_CadvF(4,i) = i      
          else
            if(cellmap(6,ID1) .gt. 0) IB1 = cellmap(6,ID1)  
            if(IB1.eq.IB2) then   
              xface_CadvF(4,i) = cellfaces(3,IB1)   
              xface_CadvF(6,i) = cellfaces(7,IB1)
              if(cellfaces(8,IB1) .gt. 0) xface_CadvF(6,i) = cellfaces(8,IB1)           
            else
              xface_CadvF(6,i) = cellfaces(3,IB2)   
              xface_CadvF(4,i) = cellfaces(3,IB2)            
            endif
          endif  !if t1,it2 > ncells
        ENDIF ! end of ID1 <=> ID2
      ENDIF !if face a wall
    enddo  !all x faces

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    !set up indexes for duv/dx (y-momentum equation)
    do i=1,numYfaces
      ID1 = yface_cells(1,i)
      ID2 = yface_cells(2,i)
      IF( .not.  yface_wall(i)) then  
        IF(dy(ID1) .gt. 1.1*dy(ID2) ) THEN  !ID1 > ID2
          IF(cellfaces(1,ID2) .eq. cellfaces(5,ID1)) THEN  !it is a right face
            !vertical faces    
            yface_CadvF(1,i) = cellfaces(3,ID2)
            yface_CadvF(7,i) = cellfaces(3,ID2)
            if(cellfaces(4,ID1) .gt. 0) then
              yface_CadvF(2,i) = cellfaces(4,ID1)
              yface_CadvF(8,i) = cellfaces(4,ID1)      
            else
              yface_CadvF(2,i) = cellfaces(3,ID1) 
              yface_CadvF(8,i) = cellfaces(3,ID1)        
            endif            
            !horizontal face at right
            IR2 = cellmap(3,ID2)
            if(ir2.gt.ncells) then
              yface_CadvF(5,i) = i
              yface_CadvF(3,i) = i  
            else         
              yface_CadvF(3,i) = cellfaces(1,IR2)   
              yface_CadvF(5,i) = cellfaces(1,IR2)   
            endif
            !horizontal face at left   
            yface_CadvF(4,i) = cellfaces(6,ID1)
            yface_CadvF(6,i) = cellfaces(6,ID1)             
          ELSE !it is a left face
            !vertical faces    
            yface_CadvF(1,i) = cellfaces(3,ID2)
            yface_CadvF(7,i) = cellfaces(3,ID2)
            if(cellfaces(4,ID1) .gt. 0) then
              yface_CadvF(2,i) = cellfaces(4,ID1)
              yface_CadvF(8,i) = cellfaces(7,ID1)      
            else
              yface_CadvF(2,i) = cellfaces(3,ID1) 
              yface_CadvF(8,i) = cellfaces(7,ID1)        
            endif 
            !horizontal face at right       
            yface_CadvF(5,i) = cellfaces(5,ID1)
            yface_CadvF(3,i) = cellfaces(5,ID1)          
            !horizontal face at left         
            IL1 = cellmap(7,ID1)
            if(il1.gt.ncells) then
              yface_CadvF(6,i) = i
              yface_CadvF(4,i) = i                
            else         
              yface_CadvF(6,i) = cellfaces(5,IL1)   
              yface_CadvF(4,i) = cellfaces(5,IL1) 
            endif
          ENDIF !end of right/left face 
        ELSEIF(dx(ID2) .gt. 1.1*dx(ID1) ) THEN  !ID2 > ID1
          IF(cellfaces(2,ID2) .eq. cellfaces(5,ID1)) THEN  !right face
            !verical faces    
            yface_CadvF(1,i) = cellfaces(3,ID2)
            yface_CadvF(7,i) = cellfaces(3,ID2)
            if(cellfaces(4,ID1) .gt. 0) then
              yface_CadvF(2,i) = cellfaces(4,ID1)
              yface_CadvF(8,i) = cellfaces(4,ID1)      
            else
              yface_CadvF(2,i) = cellfaces(3,ID1) 
              yface_CadvF(8,i) = cellfaces(3,ID1)        
            endif            
            !horizontal face at right 
            IR2 = cellmap(3,ID2)
            if(ir2.gt.ncells) then
              yface_CadvF(5,i) = i
              yface_CadvF(3,i) = i                
            else          
              yface_CadvF(5,i) = cellfaces(1,IR2)   
              yface_CadvF(3,i) = cellfaces(1,IR2) 
            endif
            !horizontal face at left   
            IL1 = cellmap(7,ID1)
            if(il1.gt.ncells) then
              yface_CadvF(6,i) = i
              yface_CadvF(4,i) = i                
            else          
              yface_CadvF(4,i) = cellfaces(5,IL1)   
              yface_CadvF(6,i) = cellfaces(5,IL1)  
            endif
          ELSE !left face
            !vertical faces 
            if(cellfaces(8,ID2) .gt. 0) then              
              yface_CadvF(1,i) = cellfaces(3,ID2)
              yface_CadvF(7,i) = cellfaces(8,ID2)
            else
              yface_CadvF(1,i) = cellfaces(3,ID2)
              yface_CadvF(7,i) = cellfaces(7,ID2)              
            endif
            if(cellfaces(4,ID1) .gt. 0) then
              yface_CadvF(2,i) = cellfaces(4,ID1)
              yface_CadvF(8,i) = cellfaces(4,ID1)              
            else
              yface_CadvF(2,i) = cellfaces(3,ID1) 
              yface_CadvF(8,i) = cellfaces(3,ID1)            
            endif
            !horizontal face at right       
            yface_CadvF(5,i) = cellfaces(2,ID2)
            yface_CadvF(3,i) = cellfaces(2,ID2)          
            !horizontal face at left         
            IL1 = cellmap(7,ID1)
            if(il1.gt.ncells) then
              yface_CadvF(6,i) = i
              yface_CadvF(4,i) = i                
            else         
              yface_CadvF(6,i) = cellfaces(5,IL1)   
              yface_CadvF(4,i) = cellfaces(5,IL1)  
            endif
          ENDIF !end of right/left face           
        ELSE  !they are the same size
          !vertical faces    
          yface_CadvF(1,i) = cellfaces(3,ID2)
          yface_CadvF(7,i) = cellfaces(3,ID2)
          if(cellfaces(4,ID1) .gt. 0) then
            yface_CadvF(2,i) = cellfaces(4,ID1)
            yface_CadvF(8,i) = cellfaces(4,ID1)      
          else
            yface_CadvF(2,i) = cellfaces(3,ID1) 
            yface_CadvF(8,i) = cellfaces(3,ID1)        
          endif 
          !horizontal face at right
          IR2 = cellmap(3,ID2)
          IR1 = cellmap(3,ID1)
          if(cellmap(4,ID1) .gt. 0) IR1 = cellmap(4,ID1)
          if(ir1.gt.ncells.or.ir2.gt.ncells) then
            yface_CadvF(5,i) = i   
            yface_CadvF(3,i) = i    
          else
            if(IR1 .eq. IR2) then
              yface_CadvF(5,i) = cellfaces(1,IR1)   
              yface_CadvF(3,i) = cellfaces(5,IR1)
              if(cellfaces(6,IR1) .gt. 0) yface_CadvF(3,i) = cellfaces(6,IR1)
            else
              yface_CadvF(5,i) = cellfaces(1,IR2)   
              yface_CadvF(3,i) = cellfaces(1,IR2)                
            endif
          endif !ir1,ir2>ncells
          !horizontal face at left 
          IL1 = cellmap(7,ID1)
          IL2 = cellmap(7,ID2)
          if(cellmap(8,ID2) .gt. 0) IL2 = cellmap(8,ID2)   
          if(il1.gt.ncells.or.il2.gt.ncells) then
            yface_CadvF(6,i) = i   
            yface_CadvF(4,i) = i    
          else
            if(IL1.eq.IL2) then   
              yface_CadvF(4,i) = cellfaces(5,IL2)   
              yface_CadvF(6,i) = cellfaces(1,IL1)
              if(cellfaces(2,IL1) .gt. 0) yface_CadvF(6,i) = cellfaces(2,IL1)           
            else
              yface_CadvF(6,i) = cellfaces(5,IL1)   
              yface_CadvF(4,i) = cellfaces(5,Il1)            
            endif
          endif !end il1,il2>ncells
        ENDIF ! end of ID1 <=> ID2
      ENDIF !if face a wall 
    enddo  !all y faces   
    
    
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! !set up indexes for duv/dx terms 
    !do i=1,numYfaces
    ! IF( .not.  yface_wall(i)) then  
    !    !use classification to determine the DadvF array values
    !
    !     
    !  !vertical faces   
    !  yface_CadvF(1,i) = cellfaces(3,yface_cells(2,i))
    !  if(cellfaces(4,yface_cells(1,i)) .gt. 0) then
    !  yface_CadvF(2,i) = cellfaces(4,yface_cells(1,i))      
    !  else
    !  yface_CadvF(2,i) = cellfaces(3,yface_cells(1,i)) 
    !  endif
    !  !yface_CadvF(3,i) = cellfaces(7,yface_cells(1,i))
    !  !if(cellfaces(8,yface_cells(2,i)) .gt. 0) then
    !  !yface_CadvF(4,i) = cellfaces(8,yface_cells(2,i))      
    !  !else
    !  !yface_CadvF(4,i) = cellfaces(7,yface_cells(2,i)) 
    !  !endif
    !  
    !  !horizontal face at right  
    !  ID1 = cellmap(3,yface_cells(2,i))
    !  !if(i .eq. 12643)write(*,*)'id1',yface_cells(2,i),cellmap(3,yface_cells(2,i))
    !  if(cellmap(4,yface_cells(1,i)) .gt. 0 ) then
    !  ID2 = cellmap(4,yface_cells(1,i))
    !  !if(i .eq. 12643)write(*,*)'id2a',yface_cells(1,i),cellmap(4,yface_cells(1,i))
    !  else
    !  ID2 = cellmap(3,yface_cells(1,i)) 
    !  !if(i .eq. 12643)write(*,*)'id2b',yface_cells(1,i),cellmap(3,yface_cells(1,i))
    !  endif
    !  IF(ID1 .EQ. ID2) THEN
    !  yface_CadvF(5,i) = cellfaces(1,id1)
    !  if(cellfaces(6,id1) .eq. 0) then
    !  yface_CadvF(3,i) = cellfaces(5,id1)
    !  else
    !  yface_CadvF(3,i) = cellfaces(6,id1)    
    !  endif    
    !  ELSE
    !  ID1dummy = .true.
    !  ID2dummy = .true.     
    !  if(ID1 .gt. ncells) ID1dummy = .false.
    !  if(ID2 .gt. ncells) ID2dummy = .false. 
    !  if(id1dummy)  yface_CadvF(5,i) = cellfaces(1,id1)
    !  if(id1dummy .and.  id2dummy) then
    !    if(dy(yface_cells(2,i)) .gt. 1.1*dy(yface_cells(1,i)) ) yface_CadvF(5,i) = cellfaces(5,id2)
    !    if(dy(yface_cells(1,i)) .gt. 1.1*dy(yface_cells(2,i)) ) yface_CadvF(5,i) = cellfaces(1,id1)
    !  endif
    !  if(.not. id1dummy .and. .not. id2dummy) yface_CadvF(5,i) = i  !use same value to prevent diffusion (advection will be zero since at boundary)
    !  !if(i .eq. 12643)  write(*,*)i,id1dummy,id2dummy
    !  if(.not. id1dummy .and. id2dummy) then
    !     ! if(i .eq. 12643)write(*,*)i,id2
    !    if(cellfaces(6,id2) .gt. 0 ) then
    !      yface_CadvF(5,i) = cellfaces(6,id2)  
    !    else
    !      yface_CadvF(5,i) = cellfaces(5,id2)  
    !    endif
    !  endif
    !  yface_CadvF(3,i) = yface_CadvF(5,i)
    !  ENDIF
    ! 
    !  !horizontal face at left
    !  ID4 = cellmap(7,yface_cells(1,i))
    !  if(cellmap(8,yface_cells(2,i)) .gt. 0 ) then
    !  ID3 = cellmap(8,yface_cells(2,i))
    !  else
    !  ID3 = cellmap(7,yface_cells(2,i))   
    !  endif
    !  IF(ID3 .EQ. ID4)THEN
    !  yface_CadvF(6,i) = cellfaces(5,id4)
    !  if(cellfaces(2,id3) .eq. 0) then
    !  yface_CadvF(4,i) = cellfaces(1,id3)
    !  else
    !  yface_CadvF(4,i) = cellfaces(2,id3)    
    !  endif       
    !  ELSE
    !  ID3dummy = .true.
    !  ID4dummy = .true.     
    !  if(ID3 .gt. ncells) ID3dummy = .false.
    !  if(ID4 .gt. ncells) ID4dummy = .false.          
    !  if(id4dummy)  yface_CadvF(6,i) = cellfaces(5,id4)
    !  if(.not. id3dummy .and. .not. id4dummy) yface_CadvF(6,i) = i  !use same value to prevent diffusion (advectin will be zero since at boundary)
    !  if(.not. id4dummy .and. id3dummy) then
    !    if(cellfaces(2,ID3) .gt. 0 ) then
    !      yface_CadvF(6,i) = cellfaces(2,id3)  
    !    else
    !      yface_CadvF(6,i) = cellfaces(1,id3)  
    !    endif
    !  endif 
    !  yface_CadvF(4,i) = yface_CadvF(6,i)      
    !  ENDIF
    !  
    ! endif  !if wall
    !enddo  !all y faces   
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !list of wall-designated faces that need inline advective fluxes
    numspecX = 0
    numspecY = 0
    if(nHstr .gt. 0) then
      do i = 1,nHstr  !for each cell string
        do j=1,H_str(i)%NCells    !for each cell in string
          ID1 = H_str(i)%Cells(j)
          if(cellmap(7,id1) .gt. ncells) numspecX = numSpecX + 1
          if(cellmap(5,id1) .gt. ncells) numspecY = numSpecY + 1          
        enddo
      enddo ! end of each cell string
    endif  !H_single  
       
    allocate(specX(numspecX),specY(numspecY))    
    numspecX = 0
    numspecY = 0
    if(nHstr .gt. 0) then
      do i = 1,nHstr  !for each cell string
        do j=1,H_str(i)%NCells    !for each cell in string
          ID1 = H_str(i)%Cells(j)
          if(cellmap(7,id1) .gt. ncells)then
            numspecX = numSpecX + 1  
            specX(numspecX) = cellfaces(7,id1)
          endif
          if(cellmap(5,id1) .gt. ncells) then
            numspecY = numSpecY + 1            
            specY(numSpecY) = cellfaces(5,id1)              
          endif
        enddo
      enddo ! end of each cell string
    endif  !H_single   
    
    do ii=1,numspecX
      i=specX(ii)
      !identify upstream face(s)           
      ID1=xface_cells(1,i)
      xface_advF(1,i) = cellfaces(3,id1)  
      xface_advF(2,i) = 0           
    enddo 
    do ii=1,numspecY
      i=specY(ii)
      !identify upstream face(s)           
      ID1=yface_cells(1,i)
      yface_advF(1,i) = cellfaces(1,id1)  
      yface_advF(2,i) = 0           
    enddo        
        
#ifdef DEBUG
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !DEBUG-CHECk OUTUT
    open(unit=3050,file='cellmap.csv')
    do i=1,ncells
        write(3050,"(8(i7,','),i7)")i,(cellmap(j,i),j=1,8)
    enddo
    close(3050)
    open(unit=3050,file='cellfaces.csv')
    do i=1,ncells
        write(3050,"(8(i7,','),i7)")i,(cellfaces(j,i),j=1,8)
    enddo 
    close(3050)
    open(unit=3050,file='xfaces.csv')
    do i=1,numxfaces
        write(3050,"(3(i20,','),i20)")i,(xface_cells(j,i),j=1,2)
    enddo 
    close(3050)    
    open(unit=3050,file='yfaces.csv')
    do i=1,numyfaces
        write(3050,"(3(i20,','),i20)")i,(yface_cells(j,i),j=1,2)
    enddo 
    close(3050)       
#endif    

    do i=1,numxfaces
      xface_q(i)=0
      xface_qn(i)=0
      xface_vel(i)=0
    enddo
    do i=1,numyfaces
      yface_q(i)=0
      yface_qn(i)=0
      yface_vel(i)=0
    enddo
    
    allocate (etan(ncellsD))
    allocate (ADVECTX(NcellsD),ADVECTY(NcellsD)) 
    allocate (cdx(ncellsD),cdy(ncellsD))
    allocate(RHOPrim(ncellsD))  
    DO I=1,NCELLSD
      ETA(I) = 0.0
      ETAN(I) = 0.0
      CDX(i) = 0.005
      CDY(i) = 0.005
      ADVECTX(I) = 0.0
      ADVECTY(I) = 0.0
      iWET(I) = 1
      RHOPrim(I) = 1.0          !CR - 01/21/2009
    ENDDO    

    if(ADVECT .or. mixing .or. saltrans) then
      allocate (Fuu(0:numxfaces),Fuv(0:numxfaces),Gvv(0:numyfaces),Gvu(0:numyfaces))
      Fuu = 0.0
      Fuv = 0.0
      Gvv = 0.0
      Gvu = 0.0
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
    do i=1,num_linktodummies
      id1 = linktodummiesTel(1,i)
      id2 = linktodummiesTel(2,i)
      eta(id2) = eta(id1)
      etan(id2) = etan(id1)
      zb(id2) = zb(id1)       
    enddo    

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
    
    !Make lists of regular and telescoping boudnary cells and faces
    !These regular cell and face IDs are processed using more efficient algorithms
    !since they have simpler connectivity
    
    !for cells
    numREGcells = 0
    numTBcells = 0
    do i=1,ncells
      if(cellmap(2,i) .eq. 0 .and. cellmap(4,i) .eq. 0 .and. &
        cellmap(6,i) .eq. 0 .and. cellmap(8,i) .eq. 0) then
        numREGcells = numREGcells + 1
      else
        numTBcells = numTBcells + 1
      endif
    enddo
    allocate (REGcells(numREgcells),TBcells(numTBcells))
    numREGcells = 0
    numTBcells = 0      
    do i=1,ncells
      if(cellmap(2,i) .eq. 0 .and. cellmap(4,i) .eq. 0 .and. &
        cellmap(6,i) .eq. 0 .and. cellmap(8,i) .eq. 0) then
        numREGcells = numREGcells + 1
        REGcells(numREGcells) = i
      else
        numTBcells = numTBcells + 1
        TBcells(numTBcells) = i        
      endif
    enddo  
    
    !for Xfaces
    numREGXfaces = 0
    numTBXfaces = 0
    do i=1,numXfaces
      reg1 = .false.
      reg2 = .false.
      reg3 = .false.
      reg4 = .false.
      reg5 = .false.
      reg6 = .false.
      if(xface_advF(2,i) .eq. 0 .and. xface_advF(4,i) .eq. 0) reg1 = .true.
      if(xface_CadvF(1,i) .eq. xface_CadvF(7,i)) reg2 = .true.
      if(xface_CadvF(2,i) .eq. xface_CadvF(8,i)) reg3 = .true.      
      if(xface_CadvF(3,i) .eq. xface_CadvF(5,i)) reg4 = .true.
      if(xface_CadvF(4,i) .eq. xface_CadvF(6,i)) reg5 = .true.  
      ID1 = xface_cells(1,i)
      ID2 = xface_cells(2,i)
      if(dx(ID1) .gt. 0.9*dx(ID2) .and. dx(ID1) .lt. 1.1*dx(ID2) ) reg6 = .true.   
      if(reg1 .and. reg2 .and. reg3 .and. reg4 .and. reg5 .and. reg6) then
        !if(xface_cells(1,i) .le. ncells .and. xface_cells(2,i) .le. ncells ) 
        numREGXfaces = numREGXfaces + 1
      else
        !if(xface_cells(1,i) .le. ncells .and. xface_cells(2,i) .le. ncells )
        numTBXfaces = numTBXfaces + 1
      endif
    enddo
    allocate (REGXfaces(numREGXfaces),TBXfaces(numTBXfaces))
    numREGXfaces = 0
    numTBXfaces = 0      
    do i=1,numXfaces
      reg1 = .false.
      reg2 = .false.
      reg3 = .false.
      reg4 = .false.
      reg5 = .false.
      reg6 = .false.
      if(xface_advF(2,i) .eq. 0 .and. xface_advF(4,i) .eq. 0) reg1 = .true.
      if(xface_CadvF(1,i) .eq. xface_CadvF(7,i)) reg2 = .true.
      if(xface_CadvF(2,i) .eq. xface_CadvF(8,i)) reg3 = .true.      
      if(xface_CadvF(3,i) .eq. xface_CadvF(5,i)) reg4 = .true.
      if(xface_CadvF(4,i) .eq. xface_CadvF(6,i)) reg5 = .true.  
      ID1 = xface_cells(1,i)
      ID2 = xface_cells(2,i)
      if(dx(ID1) .gt. 0.9*dx(ID2) .and. dx(ID1) .lt. 1.1*dx(ID2) ) reg6 = .true.   
      if(reg1  .and. reg2 .and. reg3 .and. reg4 .and. reg5 .and. reg6) then
        !if(xface_cells(1,i) .le. ncells .and. xface_cells(2,i) .le. ncells ) then  
        numREGXfaces = numREGXfaces + 1
        REGXfaces(numREGXfaces) = i
        !endif
      else
        !if(xface_cells(1,i) .le. ncells .and. xface_cells(2,i) .le. ncells ) then 
        numTBXfaces = numTBXfaces + 1
        TBXfaces(numTBXfaces) = i    
        !endif
      endif
    enddo      
    
    !for Yfaces
    numREGYfaces = 0
    numTBYfaces = 0
    do i=1,numYfaces
      reg1 = .false.
      reg2 = .false.
      reg3 = .false.
      reg4 = .false.
      reg5 = .false.
      reg6 = .false.
      if(yface_advF(2,i) .eq. 0 .and.yface_advF(4,i) .eq. 0) reg1 = .true.
      if(yface_CadvF(1,i) .eq. yface_CadvF(7,i)) reg2 = .true.
      if(yface_CadvF(2,i) .eq. yface_CadvF(8,i)) reg3 = .true.      
      if(yface_CadvF(3,i) .eq. yface_CadvF(5,i)) reg4 = .true.
      if(yface_CadvF(4,i) .eq. yface_CadvF(6,i)) reg5 = .true.      
      ID1 = yface_cells(1,i)
      ID2 = yface_cells(2,i)
      if(dy(ID1) .gt. 0.9*dy(ID2) .and. dy(ID1) .lt. 1.1*dy(ID2) ) reg6 = .true.   
      if(reg1  .and. reg2 .and. reg3 .and. reg4 .and. reg5 .and. reg6) then
        numREGYfaces = numREGYfaces + 1
      else
        numTBYfaces  = numTBYfaces + 1
      endif
    enddo
    allocate (REGYfaces(numREGYfaces),TBYfaces(numTBYfaces))
    numREGYfaces = 0
    numTBYfaces = 0      
    do i=1,numYfaces
      reg1 = .false.
      reg2 = .false.
      reg3 = .false.
      reg4 = .false.
      reg5 = .false.
      reg6 = .false.
      if(yface_advF(2,i) .eq. 0 .and.yface_advF(4,i) .eq. 0) reg1 = .true.
      if(yface_CadvF(1,i) .eq. yface_CadvF(7,i)) reg2 = .true.
      if(yface_CadvF(2,i) .eq. yface_CadvF(8,i)) reg3 = .true.      
      if(yface_CadvF(3,i) .eq. yface_CadvF(5,i)) reg4 = .true.
      if(yface_CadvF(4,i) .eq. yface_CadvF(6,i)) reg5 = .true.      
      ID1 = yface_cells(1,i)
      ID2 = yface_cells(2,i)
      if(dy(ID1) .gt. 0.9*dy(ID2) .and. dy(ID1) .lt. 1.1*dy(ID2) ) reg6 = .true.   
      if(reg1  .and. reg2 .and. reg3 .and. reg4 .and. reg5 .and. reg6) then
        numREGYfaces = numREGYfaces + 1
        REGYfaces(numREGYfaces) = i
      else
        numTBYfaces  = numTBYfaces + 1
        TBYfaces(numTBYfaces)  = i        
      endif
    enddo     

#ifdef DEBUG      
    write(*,*)'numregcells = ',numREgcells,numTBcells,ncells
    write(*,*)'numregXfaces = ',numREgXfaces,numTBXfaces     
    write(*,*)'numregYfaces = ',numREgYfaces,numTBYfaces    
      
    open(unit=3050,file='TYPE_yfaces.csv')
    do i=1,numREGyfaces
      write(3050,*)REGYfaces(i),", reg"
    enddo 
    do i=1,numTByfaces
      write(3050,*)TBYfaces(i),", tb"
    enddo     
    close(3050)   
#endif

    !for all TB faces
    allocate(yface_side(numTByfaces),xface_side(numTBxfaces))
    yface_side = 1
    xface_side = 1
    do ii=1,numTByfaces
      i=TByfaces(ii)
      if(dy(yface_cells(2,i)) .gt. 1.1*dy(yface_cells(1,i)) )then !determine if "i" face is right or left face  left = 1, right = 2
        ID1 = yface_advF(3,i)
        if(yface_advF(1,ID1) .eq. i) then
          yface_side(ii) = 1 
        else
          yface_side(ii) = 2
        endif
      endif 
    enddo  
    do ii=1,numTBxfaces
      i=TBxfaces(ii)
      if(dx(xface_cells(2,i)) .gt. 1.1*dx(xface_cells(1,i)) )then !determine if "i" face is right or left face  top = 1, bottom = 2
        ID1 = xface_advF(3,i)
        if(xface_advF(1,ID1) .eq. i) then
           xface_side(ii) = 1 
        else
           xface_side(ii) = 2
        endif
      endif       
    enddo  
    
    call EXP_gradient_indexes_tel()
    return
    end subroutine        
    
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
      subroutine EXP_gradient_indexes_tel()
      use size_def, only: ncells    
      use geo_def, only:  dy 
      use exp_telescoping
      
      implicit none
      
      integer :: i
      
    !MEB - 05/10/2016
    !This should work with Square Telescoping Cells, but not for rectangular telescoping cells.  

    
    !set up coefficients for gradient term in momentum equations for x faces
    do i = 1,numxfaces   

    if(.not. xface_wall(i)) then  
      if(DY(xface_cells(2,i)) .gt. 1.1*DY(xface_cells(1,i))) then !then grid1 = .false., one left cell, two right cells
        xface_grad1(i) = .false.
        xface_basic_orientation(i) = .true.  !setting HMIN
        if(cellmap(3,xface_cells(2,i)) .eq. xface_cells(1,i)) then  ! then right cell is upper
          if(cellmap(1,xface_cells(2,i)) .gt. ncells) then  !3rd cell is dummy cell, use simple formula
            xface_gcoef(1,i) = 1.0
            xface_gcoef(2,i) = 1.0
            xface_gcoef(3,i) =  0.0
            xface_gcoef(4,i) =  0.0 
            xface_cells(3,i) = cellmap(1,xface_cells(2,i))  !this does not matter since xface_gcoef(3) = 0.0
            xface_cells(4,i) = cellmap(1,xface_cells(2,i))  !this does not matter since xface_gcoef(4) = 0.0              
          elseif(cellmap(2,xface_cells(2,i)) .gt. 0) then !3rd cell is small  (not considering second potential bottom cell does not matter)
            xface_gcoef(1,i) =  1.
            xface_gcoef(2,i) =  0.75
            xface_gcoef(3,i) =  0.125
            xface_gcoef(4,i) =  0.125   
            xface_cells(3,i) = cellmap(2,xface_cells(2,i))
            xface_cells(4,i) = cellmap(1,xface_cells(2,i))             
          else !3rd cell is big
            xface_gcoef(1,i) =  1.0
            xface_gcoef(2,i) =  0.75
            xface_gcoef(3,i) =  0.25
            xface_gcoef(4,i) =  0.0 
            xface_cells(3,i) = cellmap(1,xface_cells(2,i))
            xface_cells(4,i) = cellmap(1,xface_cells(2,i)) !this does not matter since xface_gcoef(4) = 0.0             
          endif         
        elseif(cellmap(4,xface_cells(2,i)) .eq. xface_cells(1,i)) then !then right cell is lower
          if(cellmap(5,xface_cells(2,i)) .gt. ncells) then  !3rd cell is dummy cell, use simple formula
            xface_gcoef(1,i) =  1.0
            xface_gcoef(2,i) =  1.0
            xface_gcoef(3,i) =  0.0
            xface_gcoef(4,i) =  0.0  
            xface_cells(3,i) = cellmap(5,xface_cells(2,i))   !dummy location since goef(3,i) is zero
            xface_cells(4,i) = cellmap(5,xface_cells(2,i)) !this does not matter since xface_gcoef(4) = 0.0                     
          elseif(cellmap(6,xface_cells(2,i)) .gt. 0) then !3rd cell is small  (not considering second potential bottom cell does not matter)
            xface_gcoef(1,i) =  1.0
            xface_gcoef(2,i) =  0.75
            xface_gcoef(3,i) =  0.125
            xface_gcoef(4,i) =  0.125  
            xface_cells(3,i) = cellmap(5,xface_cells(2,i))
            xface_cells(4,i) = cellmap(6,xface_cells(2,i))                   
          else !3rd cell is big
            xface_gcoef(1,i) =  1.0
            xface_gcoef(2,i) =  0.75
            xface_gcoef(3,i) =  0.25
            xface_gcoef(4,i) =  0.0   
            xface_cells(3,i) = cellmap(5,xface_cells(2,i))
            xface_cells(4,i) = cellmap(5,xface_cells(2,i)) !this does not matter since xface_gcoef(4) = 0.0                     
          endif
        else
          write(*,*)'problem with grid indexing, xface i = ,',i
        endif
         
      elseif(DY(xface_cells(2,i)) .lt. 0.9*DY(xface_cells(1,i))) then  !the grid1 = .false., two left cells, one right cell
        xface_grad1(i) = .false.
        xface_basic_orientation(i) = .false.  !setting HPLUS
        if(cellmap(8,xface_cells(1,i)) .eq. xface_cells(2,i)) then  ! then left cell is upper
          if(cellmap(1,xface_cells(1,i)) .gt. ncells) then  !3rd cell is dummy cell, use simple formula
            xface_gcoef(2,i) =  1.0
            xface_gcoef(1,i) =  1.0
            xface_gcoef(3,i) =  0.0
            xface_gcoef(4,i) = 0.0
            xface_cells(3,i) = cellmap(1,xface_cells(1,i))
            xface_cells(4,i) = cellmap(1,xface_cells(1,i))  !this does not matter since xface_gcoef(4) = 0.0             
          elseif(cellmap(2,xface_cells(1,i)) .gt. 0) then !3rd cell is small  (not considering second potential upper cell does not matter)
            xface_gcoef(2,i) =  1.0
            xface_gcoef(1,i) =  0.75
            xface_gcoef(3,i) =  0.125
            xface_gcoef(4,i) =  0.125
            xface_cells(3,i) = cellmap(1,xface_cells(1,i))
            xface_cells(4,i) = cellmap(2,xface_cells(1,i))              
          else !3rd cell is big
            xface_gcoef(2,i) =  1.0
            xface_gcoef(1,i) =  0.75
            xface_gcoef(3,i) =  0.25
            xface_gcoef(4,i) =  0.0  
            xface_cells(3,i) = cellmap(1,xface_cells(1,i))
            xface_cells(4,i) = cellmap(1,xface_cells(1,i))             
          endif         
        elseif(cellmap(7,xface_cells(1,i)) .eq. xface_cells(2,i)) then !then left cell is lower
          if(cellmap(5,xface_cells(1,i)) .gt. ncells) then  !3rd cell is dummy cell, use simple formula
            xface_gcoef(2,i) =  1.0
            xface_gcoef(1,i) =  1.0
            xface_gcoef(3,i) =  0.0
            xface_gcoef(4,i) =  0.0
            xface_cells(3,i) = cellmap(5,xface_cells(1,i))
            xface_cells(4,i) = cellmap(5,xface_cells(1,i))  !this does not matter since xface_gcoef(4) = 0.0 
          elseif(cellmap(6,xface_cells(1,i)) .gt. 0) then !3rd cell is small  (not considering second potential bottom cell does not matter)
            xface_gcoef(2,i) =  1.0
            xface_gcoef(1,i) =  0.75
            xface_gcoef(3,i) =  0.125
            xface_gcoef(4,i) =  0.125
            xface_cells(3,i) = cellmap(5,xface_cells(1,i))
            xface_cells(4,i) = cellmap(6,xface_cells(1,i))  
          else !3rd cell is big
            xface_gcoef(2,i) =  1.0
            xface_gcoef(1,i) =  0.75
            xface_gcoef(3,i) =  0.25
            xface_gcoef(4,i) =  0.0  
            xface_cells(3,i) = cellmap(5,xface_cells(1,i))
            xface_cells(4,i) = cellmap(5,xface_cells(1,i))
          endif
        else
          write(*,*)'problem with grid indexing, xface i = ,',i
        endif       
      else
        xface_grad1(i) = .true. 
        xface_gcoef(1,i) = 1.0
        xface_gcoef(2,i) = 1.0
        xface_gcoef(3,i) = 0.0
        xface_gcoef(4,i) = 0.0 
        xface_cells(3,i) = xface_cells(2,i)  !this does not matter since xface_gcoef(3) = 0.0
        xface_cells(4,i) = xface_cells(1,i)  !this does not matter since xface_gcoef(3) = 0.0             
      endif
    endif
    enddo

    !set up coefficients for gradient term in momentum equations for y faces    

    do i = 1,numyfaces

    yface_cells(3,i) = 0        
    if(.not. yface_wall(i)) then   
    yface_basic_orientation(i) = .true.   !setting HMIN
       
       if(DY(yface_cells(2,i)) .gt. 1.1*DY(yface_cells(1,i))) then !then grid1 = .false., one south cell, two north cells
       yface_grad1(i) = .false.
       
         if(cellmap(1,yface_cells(2,i)) .eq. yface_cells(1,i)) then  ! then north cell is left
            if(cellmap(7,yface_cells(2,i)) .gt. ncells) then  !3rd cell is dummy cell, use simple formula
              yface_gcoef(1,i) =  1.0
              yface_gcoef(2,i) =  1.0
              yface_gcoef(3,i) =  0.0
              yface_gcoef(4,i) =  0.0
              yface_cells(3,i) = cellmap(7,yface_cells(2,i))
              yface_cells(4,i) = cellmap(7,yface_cells(2,i))              
            elseif(cellmap(8,yface_cells(2,i)) .gt. 0) then !3rd cell is small  (not considering second potential bottom cell does not matter)
              yface_gcoef(1,i) =  1.0
              yface_gcoef(2,i) = 0.75
              yface_gcoef(3,i) = 0.125
              yface_gcoef(4,i) = 0.125 
              yface_cells(3,i) = cellmap(8,yface_cells(2,i))
              yface_cells(4,i) = cellmap(7,yface_cells(2,i))                 
            else !3rd cell is big
              yface_gcoef(1,i) =  1.0
              yface_gcoef(2,i) =  0.75
              yface_gcoef(3,i) =  0.25
              yface_gcoef(4,i) =  0.0
              yface_cells(3,i) = cellmap(7,yface_cells(2,i))
              yface_cells(4,i) = cellmap(7,yface_cells(2,i))                 
            endif         
         elseif(cellmap(2,yface_cells(2,i)) .eq. yface_cells(1,i)) then !then north cell is right
            if(cellmap(3,yface_cells(2,i)) .gt. ncells) then  !3rd cell is dummy cell, use simple formula
              yface_gcoef(1,i) =  1.0
              yface_gcoef(2,i) =  1.0
              yface_gcoef(3,i) =  0.0
              yface_gcoef(4,i) =  0.0
              yface_cells(3,i) = cellmap(3,yface_cells(2,i))
              yface_cells(4,i) = cellmap(3,yface_cells(2,i))             
            elseif(cellmap(4,yface_cells(2,i)) .gt. 0) then !3rd cell is small  (not considering second potential bottom cell does not matter)
              yface_gcoef(1,i) =  1.0
              yface_gcoef(2,i) =  0.75
              yface_gcoef(3,i) =  0.125
              yface_gcoef(4,i) =  0.125
              yface_cells(3,i) = cellmap(3,yface_cells(2,i))
              yface_cells(4,i) = cellmap(4,yface_cells(2,i))                
            else !3rd cell is big
              yface_gcoef(1,i) =  1.0
              yface_gcoef(2,i) =  0.75
              yface_gcoef(3,i) =  0.25
              yface_gcoef(4,i) =  0.0 
              yface_cells(3,i) = cellmap(3,yface_cells(2,i))
              yface_cells(4,i) = cellmap(3,yface_cells(2,i))                 
            endif
         else
           write(*,*)'probelm with grid indexing, yface i = ,',i
         endif
       elseif(DY(yface_cells(2,i)) .lt. 0.9*DY(yface_cells(1,i))) then !the grid1 = .false., two south cells, one north cell
         yface_grad1(i) = .false.
         yface_basic_orientation(i) = .false.    
       
         if(i.eq.569) write(*,*)'A ', yface_basic_orientation(i)
         if(cellmap(6,yface_cells(1,i)) .eq. yface_cells(2,i)) then  ! then south cell is left
           if(cellmap(7,yface_cells(1,i)) .gt. ncells) then  !3rd cell is dummy cell, use simple formula
             yface_gcoef(2,i) =  1.0
             yface_gcoef(1,i) =  1.0
             yface_gcoef(3,i) =  0.0
             yface_gcoef(4,i) =  0.0
             yface_cells(3,i) = cellmap(7,yface_cells(1,i))
             yface_cells(4,i) = cellmap(7,yface_cells(1,i))   
             if(i.eq.569) write(*,*)'B '
          elseif(cellmap(8,yface_cells(1,i)) .gt. 0) then !3rd cell is small  (not considering second potential upper cell does not matter)
             yface_gcoef(2,i) =  1.0
             yface_gcoef(1,i) =  0.75
             yface_gcoef(3,i) =  0.125
             yface_gcoef(4,i) =  0.125
             yface_cells(3,i) = cellmap(7,yface_cells(1,i))
             yface_cells(4,i) = cellmap(8,yface_cells(1,i))  
             if(i.eq.569) write(*,*)'C '              
           else !3rd cell is big             
             yface_gcoef(2,i) =  1.0
             yface_gcoef(1,i) =  0.75
             yface_gcoef(3,i) =  0.25
             yface_gcoef(4,i) =  0.0
             yface_cells(3,i) = cellmap(7,yface_cells(1,i))
             yface_cells(4,i) = cellmap(7,yface_cells(1,i))  
             if(i.eq.569) write(*,*)'AA ', yface_cells(3,i), yface_cells(4,i)               
           endif         
         elseif(cellmap(5,yface_cells(1,i)) .eq. yface_cells(2,i)) then !then south cell is right
           if(i.eq.569) write(*,*)'D '  
           if(cellmap(3,yface_cells(1,i)) .gt. ncells) then  !3rd cell is dummy cell, use simple formula
             yface_gcoef(2,i) =  1.0
             yface_gcoef(1,i) =  1.0
             yface_gcoef(3,i) =  0.0
             yface_gcoef(4,i) =  0.0
             yface_cells(3,i) = cellmap(3,yface_cells(1,i))
             yface_cells(4,i) = cellmap(3,yface_cells(1,i))             
           elseif(cellmap(4,yface_cells(1,i)) .gt. 0) then !3rd cell is small  (not considering second potential bottom cell does not matter)
             yface_gcoef(2,i) =  1.0
             yface_gcoef(1,i) =  0.75
             yface_gcoef(3,i) =  0.125
             yface_gcoef(4,i) =  0.125
             yface_cells(3,i) = cellmap(4,yface_cells(1,i))
             yface_cells(4,i) = cellmap(3,yface_cells(1,i))              
           else !3rd cell is big
             yface_gcoef(2,i) =  1.0
             yface_gcoef(1,i) =  0.75
             yface_gcoef(3,i) =  0.25
             yface_gcoef(4,i) =  0.0
             yface_cells(3,i) = cellmap(3,yface_cells(1,i))
             yface_cells(4,i) = cellmap(3,yface_cells(1,i))             
           endif
         else
           write(*,*)'probelm with grid indexing, yface i = ,',i
         endif       
       else
         yface_grad1(i) = .true. 
         yface_gcoef(1,i) = 1
         yface_gcoef(2,i) = 1
         yface_gcoef(3,i) = 0.0
         yface_gcoef(4,i) = 0.0
         yface_cells(3,i) = cellmap(3,yface_cells(1,i))
         yface_cells(4,i) = cellmap(3,yface_cells(1,i))  
       endif
    endif

    enddo    

    return
    end subroutine EXP_gradient_indexes_tel