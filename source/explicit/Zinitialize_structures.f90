!********************************************************************************    
      subroutine initialize_structures_CWR
!********************************************************************************    
      use EXP_Global_def,     only: etan, drydep, adeq
      USE EXP_transport_def,  only: salt, adss, cohes, cohesive
      use EXP_Structures_def, only: structures, srm_on, hasstruct, haspor, hashgt, hashc, hasa, hasb, srmu, srm, srmv, cul_on, cul
      use sal_def,  only: saltrans
      use flow_def, only: eta
      use sed_def,  only: nhard, idhard, hardbed
      use size_def, only: ncells, ncellsd
      use geo_def, only: dx,dy,idmap,zb,cell2cell
      
      implicit none
      integer i,j,id,id2,ncellsU,ncellsV,kp,ij,id1
      real totdepth

      if(structures) then
        if(SRM_on) then
          ALLOCATE (HASSTRUCT(ncellsD),HASPOR(ncellsD),HASHGT(ncellsD),HASHC(ncellsD),HASA(ncellsD),HASB(ncellsD))
          HASSTRUCT = .false.
          !get CMS cell ID for each structure cell
          do j=1,SRM%ncells
            ID = SRM%CELLS(j)
            ID2 = idmap(ID)
            SRM%CELLS(j) = ID2
            SRM%IVAL(j) = ID  !this is used to output during steering
            HASSTRUCT(ID2) = .true.
            HASPOR(id2) = SRM%por(j)
            HASHGT(id2) = SRM%hgt(j)
            HASHC(id2) = SRM%HC(j)
            HASA(id2) = SRM%A(j)
            HASB(id2) = SRM%B(j)            
          enddo

          !identify cells that have U vel modified
          !first see how many
          ncellsU = 0
          do j=1,SRM%ncells
            ID = SRM%cells(j)
            ncellsU = ncellsU + 1
            ID2 = cell2cell(2,id)
            if(.not. hasstruct(id2)) ncellsU  = ncellsU + 1
          enddo
          !then allocate and process
          SRMU%ncells = ncellsU
          allocate (SRMU%cells(ncellsU),SRMU%por(ncellsU),SRMU%hgt(ncellsU), &
                    SRMU%L(ncellsU),SRMU%hc(ncellsU),SRMU%a(ncellsU),SRMU%b(ncellsU))
          ncellsU=0
          do j=1,SRM%ncells
            ID = SRM%cells(j)
            ncellsU = ncellsU + 1
            SRMU%cells(ncellsU) = ID
            ID2 = cell2cell(4,id)
            if(hasstruct(id2)) then
              SRMU%L(ncellsU) = 0.5*(dx(id)+dx(id2))
              SRMU%por(ncellsU) = 2./(1./haspor(id)+1./haspor(id2))
              SRMU%HC(ncellsU) = 2./(1./hasHC(id)+1./hasHC(id2))           
              SRMU%hgt(ncellsU) = 0.5*(hashgt(id)+hashgt(id)) 
              SRMU%a(ncellsU) = hasa(id)
              SRMU%b(ncellsU) = hasb(id)         
            else
              SRMU%L(ncellsU) = 0.5*dx(id)
              SRMU%por(ncellsU) = haspor(id)
              SRMU%HC(ncellsU) = hasHC(id)
              SRMU%hgt(ncellsU) = hashgt(id)
              SRMU%a(ncellsU) = hasa(id)
              SRMU%b(ncellsU) = hasb(id)           
            endif
            ID2 = cell2cell(2,id)
            if(.not. hasstruct(id2)) then
              ncellsU  = ncellsU + 1 
              SRMU%cells(ncellsU) = ID2
              SRMU%L(ncellsU) = 0.5*dx(id)
              SRMU%por(ncellsU) = haspor(id)
              SRMU%HC(ncellsU) = hasHC(id)          
              SRMU%hgt(ncellsU) = hashgt(id)
              SRMU%a(ncellsU) = hasa(id)
              SRMU%b(ncellsU) = hasb(id)          
            endif
          enddo     
          write(*,*)'Structures - U face mod = ',ncellsU,SRMU%ncells
      
          !identify cells that have V vel modified
          !first see how many
          ncellsV = 0
          do j=1,SRM%ncells
            ID = SRM%cells(j)
            ncellsV = ncellsV + 1
            ID2 = cell2cell(1,id)
            if(.not. hasstruct(id2)) ncellsV  = ncellsV + 1
          enddo
          !then allocate and process
          srmV%ncells = ncellsV
          allocate (srmV%cells(ncellsV),srmV%por(ncellsV),srmV%hgt(ncellsV),  &
                    srmV%L(ncellsV),srmV%hc(ncellsV),srmV%a(ncellsV),srmV%b(ncellsV))
          ncellsV = 0
          do j=1,SRM%ncells
            ID = SRM%cells(j)
            ncellsV = ncellsV + 1
            srmV%cells(ncellsV) = ID
            ID2 = cell2cell(3,id)
            if(hasstruct(id2)) then
              srmV%L(ncellsV) = 0.5*(dy(id)+dy(id2))
              srmV%por(ncellsV) = 2./(1./haspor(id)+1./haspor(id2))
              srmV%HC(ncellsV) = 2./(1./hasHC(id)+1./hasHC(id2))           
              srmV%hgt(ncellsV) = 0.5*(hashgt(id)+hashgt(id)) 
              srmV%a(ncellsV) = hasa(id)
              srmV%b(ncellsV) = hasb(id)                    
            else
              srmV%L(ncellsV) = 0.5*dy(id)
              srmV%por(ncellsV) = haspor(id)
              srmV%HC(ncellsV) = hasHC(id)         
              srmV%hgt(ncellsV) = hashgt(id)
              srmV%a(ncellsV) = hasa(id)
              srmV%b(ncellsV) = hasb(id)            
            endif
            ID2 = cell2cell(1,id)
            if(.not. hasstruct(id2)) then
              ncellsV  = ncellsV + 1 
              srmV%cells(ncellsV) = ID2
              srmV%L(ncellsV) = 0.5*dy(id)
              srmV%por(ncellsV) = haspor(id)
              srmV%HC(ncellsV) = hasHC(id)         
              srmV%hgt(ncellsV) = hashgt(id)
              srmV%a(ncellsV) = hasa(id)
              srmV%b(ncellsV) = hasb(id)           
            endif        
          enddo 
          write(*,*)"SRM - u faces:"
          write(*,*)(SRMU%cells(kp),kp=1,ncellsU)
          write(*,*)"SRM - v faces:"
          write(*,*)(srmV%cells(kp),kp=1,ncellsV)      
          write(*,*)'Structures - V face mod = ',ncellsV,srmV%ncells     
        
          DEALLOCATE (HASSTRUCT,HASPOR,HASHGT,HASHC,HASA,HASB)
      
          !modify cell bottom elevation for RM structure cells
          do j=1,SRM%ncells
            ID = SRM%cells(j) 
            SRM%HGT(J) = -zb(ID)
          enddo
          do j=1,SRM%ncells
            ID = SRM%cells(j) 
            zb(ID) = -SRM%base
          enddo
          !if RM cell is hardbottom, then modify hardbottom to be new depth
          do j=1,SRM%ncells
            ID = SRM%cells(j) 
            do ij=1,nhard
              id2 = idhard(ij)
              if(id.eq.id2) hardbed(ij) = -zb(id)
            enddo      
          enddo      
          !modify initial water elevation (was set to top of RM previously)
          do j=1,SRM%ncells
            ID = SRM%cells(j) 
            eta(ID) = 0.0
            etan(id) = 0.0
            totdepth = eta(id)+ (-zb(id))                     !Added parens around unary operator - 4/20/17 MEB
            if(totdepth.le.0.9*drydep) then    !CWRSALFIX
              eta(id) = 0.9d0*drydep +zb(id)
              etan(id) = eta(id)
            endif
          enddo      
          IF(SALTrans) THEN
            do j=1,SRM%ncells
              ID = SRM%cells(j) 
              salt(id)%vol = (etan(id)-zb(id))*dx(id)*dy(id)
            enddo 
          ENDIF        
          IF(ADEQ) THEN
            do j=1,SRM%ncells
              ID = SRM%cells(j) 
              adss(id)%vol = (etan(id)-zb(id))*dx(id)*dy(id)
            enddo 
          ENDIF
          IF(COHESIVE) THEN
            do j=1,SRM%ncells
              ID = SRM%cells(j) 
              cohes(id)%vol = (etan(id)-zb(id))*dx(id)*dy(id)
            enddo 
          ENDIF            
        endif !srm_on
      
        if(CUL_on) then
          allocate (cul%flow(cul%num)) 
          if(saltrans) allocate(cul%saltF(cul%num))      
          do j=1,cul%num
            ID = cul%cells1(j)
            ID2 = idmap(ID)
            cul%cells1(j) = ID2
            ID = cul%cells2(j)
            ID2 = idmap(ID)
            cul%cells2(j) = ID2
            cul%flow(j) = 0.0
          enddo  
          !find number of culvert cells with unique ids
          ALLOCATE (HASSTRUCT(ncells))
          HASSTRUCT = .false.     
          do j=1,cul%num
            id1 = cul%cells1(j)
            id2 = cul%cells2(j)
            HASSTRUCT(id1) = .true.
            HASSTRUCT(id2) = .true.
          enddo  
          cul%numui = 0    
          do i=1,ncells
            if(hasstruct(i))  cul%numui = cul%numui + 1
          enddo
          allocate (cul%cellUI(cul%numUI))
          cul%numui = 0
          do i=1,ncells
            if(hasstruct(i)) then
              cul%numui = cul%numui + 1      
              cul%cellUI(cul%numui) = i
            endif
          enddo
          deallocate (hasstruct)
      
          do i=1,cul%numui
            write(*,*)i,cul%cellUI(i)
          enddo
        endif !cul%on
      endif !structures      
      
      return
      end subroutine