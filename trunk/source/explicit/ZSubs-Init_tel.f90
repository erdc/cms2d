!********************************************************************************
!   initalizes arrays for explicit algorithm flow and salinity boundary conditions
!   determines if cell string is for U or V flow
!   determines sign to assing to flow input so positive flow is into grid
!   redstributes flow input so it is uniform over string (imapct only wehn cells
!       are not of uniform width
!   re-assigns salintiy concentratino BC to dummy cells for for advection into 
!       grid
!*******************************************************************************
      SUBROUTINE Initialize_BC_tel
      use EXP_Global_def
      USE EXP_bndcond_def 
      use geo_def, only: dx,dy,cell2cell
      use size_def     
      use flow_def
      use comvarbl
      use out_def, only: goutfile
      use bnd_def
      use sal_def
      use EXP_TELESCOPING
      USE diag_def, only: dgunit,dgfile
      
      IMPLICIT NONE
      !LOCAL VARIABLES
      INTEGER i,j,id_t
      INTEGER ICNT_E,ICNT_N,ICNT_W,ICNT_S,ICNTT
      REAL(ikind) sum
      LOGICAL isOpen


      if(nQstr .gt. 0) then
                
        !determine orientation of cell string to determine to which face the flow is applied
        allocate (QstringEXP(nQstr))
        do i=1,nQstr
          icnt_E = 0
          icnt_N = 0
          icnt_W = 0
          icnt_S = 0
          QstringEXP(I)%SGN= 1
          do j=1,Q_Str(I)%NCELLS
            ID_T = Q_Str(I)%cells(j)
            if(cellmap(3,ID_T) .gt. NCELLS) icnt_E = icnt_E + 1
            if(cellmap(7,ID_T) .gt. NCELLS) icnt_W = icnt_W + 1
            if(cellmap(1,ID_T) .gt. NCELLS) icnt_N = icnt_N + 1
            if(cellmap(5,ID_T) .gt. NCELLS) icnt_S = icnt_S + 1
          enddo
          icntT = max(icnt_N,icnt_S,icnt_E,icnt_W)

          !reassign boundary conditions on east and north to adjacent dummy cells
          !depending on the results of the string orientation
          !also determine if Q applied to v or u face (via vface value true or false)
          if(icntT.eq.icnt_N) then
            QstringEXP(I)%Vface = .true.
            QstringEXP(I)%sgn = -1
            !put the Q value on a m2/s basis - divide .qcurv by total lenght of cell string            
            sum=0
            do j=1,Q_Str(I)%NCELLS
              ID_T = Q_Str(I)%cells(j)
              Q_str(I)%cells(j) = cellfaces(1,ID_T)
              sum = sum + dx(ID_T)
            enddo  
            do j=1,Q_str(i)%Ntimes
            Q_str(i)%qcurv(j)	=   Q_str(i)%qcurv(j)/sum
            enddo            
          endif
          if(icntT.eq.icnt_E) then
            QstringEXP(I)%sgn = -1        
            QstringEXP(I)%vface = .false.
             !put the Q value on a m2/s basis - divide .qcurv by total lenght of cell string           
            sum = 0.0
            do j=1,Q_Str(I)%NCELLS
              ID_T = Q_Str(I)%cells(j)
              Q_str(I)%cells(j) = cellfaces(3,ID_T)
              sum = sum + dy(ID_T)
            enddo
            do j=1,Q_str(i)%Ntimes
            Q_str(i)%qcurv(j)	=   Q_str(i)%qcurv(j)/sum
            enddo
          endif
          ! no need to shift cells for W and S cases, since these cell faces are the 
          !proper ones for assigning flows for cellstrings on west or south edges         
          if(icntT.eq.icnt_S) then
              QstringEXP(I)%vface = .true.
           !put the Q value on a m2/s basis - divide .qcurv by total lenght of cell string     
            sum=0.0
            do j=1,Q_Str(I)%NCELLS
              ID_T = Q_Str(I)%cells(j)
              Q_str(I)%cells(j) = cellfaces(5,ID_T) 
              sum = sum + dx(ID_T)
            enddo
            do j=1,Q_str(i)%Ntimes
            Q_str(i)%qcurv(j)	=   Q_str(i)%qcurv(j)/sum
            enddo              
          endif
          if(icntT.eq.icnt_W) then
              QstringEXP(I)%vface = .false.
            !put the Q value on a m2/s basis - divide .qcurv by total lenght of cell string               
            sum = 0.0
            do j=1,Q_Str(I)%NCELLS
              ID_T = Q_Str(I)%cells(j)
              Q_str(I)%cells(j) = cellfaces(7,ID_T)             
              sum = sum + dy(ID_T)
            enddo
            do j=1,Q_str(i)%Ntimes
            Q_str(i)%qcurv(j)	=   Q_str(i)%qcurv(j)/sum
            enddo
        endif
          
         !reassign boundary condition cells to dummy cells based on cell string location         
         !if(saltrans .and. nsalstr .gt. 0) then
         ! if(icntT.eq.icnt_N) then
         !   do j=1,sal_str(I)%NCELLS
         !     ID_T = sal_str(I)%cells(j)
         !     sal_str(I)%cells(j) = cellmap(1,ID_T)
         !   enddo
         ! endif
         ! if(icntT.eq.icnt_E) then
         !   do j=1,sal_str(I)%NCELLS
         !     ID_T = sal_str(I)%cells(j)
         !     sal_str(I)%cells(j) = cellmap(3,ID_T)
         !   enddo
         ! endif
         ! if(icntT.eq.icnt_S) then
         !   do j=1,sal_str(I)%NCELLS
         !     ID_T = sal_str(I)%cells(j)
         !     sal_str(I)%cells(j) = cellmap(4,ID_T)
         !   enddo
         ! endif
         ! if(icntT.eq.icnt_W) then
         !   do j=1,sal_str(I)%NCELLS
         !     ID_T = sal_str(I)%cells(j)
         !     sal_str(I)%cells(j) = cellmap(7,ID_T)
         !   enddo
         ! endif  
         !endif       
        
        enddo !end of NQdriver

        inquire(unit=DGUNIT,opened=isOpen)
        if(.not. isOpen) then
          open(dgunit,file=dgfile,access='append') 
        endif
        write(DGUNIT,*)'finished initializing flow time series bcs'
        close (DGUNIT)
        
        Endif !(Q_single)

      RETURN
      END SUBROUTINE !INITIALIZE_BC
      
    !********************************************************************************
      SUBROUTINE INITIALIZE_WABC_tel
#include "CMS_cpp.h"      
      use EXP_Global_def 
      USE EXP_bndcond_def  
      use size_def
      use flow_def   
      use cms_def, only: radpath,noptset
      use comvarbl, only: timehrs
      use sed_def
      use met_def, only: windconst,windvar
      use bnd_def
      use geo_def, only: dx,dy
      use EXP_TELESCOPING, only: cellmap
      
      IMPLICIT NONE
      !local variables
      INTEGER i,j,NDIM,ID1,ID2
      character str*3
     
      !if flow and waves are simulated, then WABC is turned on
      if(noptset==3) WABC = .true.
      write(*,*)'WABC = ',wabc
      write(*,*)'noptset = ',noptset
      
      !the WABC also applies for wind stresses, both spatially constant and variable
      if(WABC .or. windvar .or. windconst) then           

        if(nHstr .gt. 0) then    !initialize WABC
          ALLOCATE (SWABC(nHstr))
          do i=1,nHstr
            Ndim = H_Str(i)%ncells
            ALLOCATE (SWABC(i)%Q(Ndim),SWABC(i)%ETA(Ndim),               &
              SWABC(i)%del(Ndim),SWABC(i)%SGN(Ndim),SWABC(i)%CELL(Ndim), &
              SWABC(i)%RAD(Ndim),SWABC(i)%WDTH(Ndim),SWABC(i)%QN(Ndim))

            WRITE(STR,'(I0)')I
            SWABC(i)%name = 'SWABC_'//trim(STR)	
            do j=1,H_Str(i)%ncells - 1
              id1 = H_Str(i)%cells(j)
              id2 = H_Str(i)%cells(j+1)
              if(id1.eq.id2) cycle
              if(id2.eq.cellmap(5,id1) .or. id2.eq.cellmap(6,id1)) then  ! cell string directed to the south
                SWABC(i)%sgn(j) = 1
                SWABC(i)%rad(j) = 2
                SWABC(i)%cell(j) = id1
                SWABC(i)%Q(j) = 0
                SWABC(i)%QN(j) = 0
                SWABC(i)%eta(j) = 0
                SWABC(i)%del(j) = (dy(id1)+dy(id2))/2.
                SWABC(i)%wdth(j) = min(dx(id1),dx(id2))
              endif
              if(id2.eq.cellmap(1,id1) .or. id2.eq.cellmap(2,id1) ) then  ! cell string directed to the north
                SWABC(i)%sgn(j) = -1
                SWABC(i)%rad(j) = 2
                SWABC(i)%cell(j) = id2
                SWABC(i)%Q(j) = 0
                SWABC(i)%QN(j) = 0
                SWABC(i)%eta(j) = 0
                SWABC(i)%del(j) = (dy(id1)+dy(id2))/2.
                SWABC(i)%wdth(j) = min(dx(id1),dx(id2))
              endif	
               if(id2.eq.cellmap(7,id1) .or. id2.eq.cellmap(8,id1) ) then   ! cell string directed to the west
                SWABC(i)%sgn(j) = 1
                SWABC(i)%rad(j) = 1
                SWABC(i)%cell(j) = id1
                SWABC(i)%Q(j) = 0
                SWABC(i)%QN(j) = 0
                SWABC(i)%eta(j) = 0
                SWABC(i)%del(j) = (dx(id1)+dx(id2))/2.
                SWABC(i)%wdth(j) = min(dy(id1),dy(id2))
              endif
               if(id2.eq.cellmap(3,id1) .or. id2.eq.cellmap(4,id1) ) then   ! cell string directed to the east
                SWABC(i)%sgn(j) = -1                   
                SWABC(i)%rad(j) = 1
                SWABC(i)%cell(j) = id2
                SWABC(i)%Q(j) = 0
                SWABC(i)%QN(j) = 0
                SWABC(i)%eta(j) = 0
                SWABC(i)%del(j) = (dx(id1)+dx(id2))/2.
                SWABC(i)%wdth(j) = min(dy(id1),dy(id2))
              endif
            enddo
            j=H_Str(i)%Ncells
            SWABC(i)%sgn(j) = 0
            SWABC(i)%rad(j) = 0
            SWABC(i)%cell(j) = 0
            SWABC(i)%Q(j) = 0
            SWABC(i)%QN(j) = 0
            SWABC(i)%eta(j) = 0
            SWABC(i)%del(j) = 0
            SWABC(i)%wdth(j) = 0
            
             do j=1,H_Str(i)%ncells - 1
              id1 = H_Str(i)%cells(j)
              id2 = H_Str(i)%cells(j+1)
              !write(*,*)'---------------------------------------'
              !write(*,*)'j,id1,id2 ',j,id1,id2
              !write(*,*)'sgn ',  SWABC(i)%sgn(j) 
              !write(*,*)'rad ',  SWABC(i)%rad(j) 
              !write(*,*)'cell ',  SWABC(i)%cell(j) 
              !write(*,*)'del',  SWABC(i)%del(j) 
              !write(*,*)'width ',  SWABC(i)%wdth(j) 
             enddo  

          enddo !ndriver 
          
          
        endif !H_single

        if(nTHstr .gt. 0) then
          ALLOCATE (TWABC(1))
          i=1
          Ndim = TH_Str(i)%ncells
          ALLOCATE (TWABC(i)%Q(Ndim),TWABC(i)%ETA(Ndim),                &
            TWABC(i)%del(Ndim),TWABC(i)%SGN(Ndim),TWABC(i)%CELL(Ndim),  &
            TWABC(i)%RAD(Ndim),TWABC(i)%WDTH(Ndim),TWABC(i)%QN(Ndim))
          WRITE(STR,'(I0)')I
          TWABC(i)%name = 'TWABC_'//trim(STR)
          do j=1,TH_Str(i)%ncells - 1
            id1 = TH_Str(I)%Cells(j)
            id2 = TH_Str(I)%Cells(j+1)
            if(id2.eq.cellmap(5,id1) .or. id2.eq.cellmap(6,id1)) then
              TWABC(i)%sgn(j) = 1
              TWABC(i)%rad(j) = 2
              TWABC(i)%cell(j) = id1
              TWABC(i)%Q(j) = 0
              TWABC(i)%QN(j) = 0
              TWABC(i)%eta(j) = 0
              TWABC(i)%del(j) = (dy(id1)+dy(id2))/2.
              TWABC(i)%wdth(j) = min(dx(id1),dx(id2))
            endif
            if(id2.eq.cellmap(1,id1) .or. id2.eq.cellmap(2,id1)) then
              TWABC(i)%sgn(j) = -1
              TWABC(i)%rad(j) = 2
              TWABC(i)%cell(j) = id2
              TWABC(i)%Q(j) = 0
              TWABC(i)%QN(j) = 0
              TWABC(i)%eta(j) = 0
              TWABC(i)%del(j) = (dy(id1)+dy(id2))/2.
              TWABC(i)%wdth(j) = min(dx(id1),dx(id2))
            endif	
            if(id2.eq.cellmap(7,id1) .or. id2.eq.cellmap(8,id1)) then
              TWABC(i)%sgn(j) = 1
              TWABC(i)%rad(j) = 1
              TWABC(i)%cell(j) = id1
              TWABC(i)%Q(j) = 0
              TWABC(i)%QN(j) = 0
              TWABC(i)%eta(j) = 0
              TWABC(i)%del(j) = (dx(id1)+dx(id2))/2.
              TWABC(i)%wdth(j) = min(dy(id1),dy(id2))
            endif
            if(id2.eq.cellmap(3,id1) .or. id2.eq.cellmap(4,id1)) then
              TWABC(i)%sgn(j) = -1
              TWABC(i)%rad(j) = 1
              TWABC(i)%cell(j) = id2
              TWABC(i)%Q(j) = 0
              TWABC(i)%QN(j) = 0
              TWABC(i)%eta(j) = 0
              TWABC(i)%del(j) = (dx(id1)+dx(id2))/2.
              TWABC(i)%wdth(j) = min(dy(id1),dy(id2))
            endif
          Enddo
          j=TH_Str(i)%ncells
          TWABC(i)%sgn(j) = 0
          TWABC(i)%rad(j) = 0
          TWABC(i)%cell(j) = 0
          TWABC(i)%Q(j) = 0
          TWABC(i)%QN(j) = 0
          TWABC(i)%eta(j) = 0
          TWABC(i)%del(j) = 0
          TWABC(i)%wdth(j) = 0
        endif !H_tide

        if(nMHstr .ge. 0) then
          ALLOCATE (MWABC(nMHstr))
          do i=1,nMHstr
            Ndim = MH_Str(i)%ncells	
            ALLOCATE (MWABC(i)%Q(Ndim),MWABC(i)%ETA(Ndim),               &
              MWABC(i)%del(Ndim),MWABC(i)%SGN(Ndim),MWABC(i)%CELL(Ndim), &
              MWABC(i)%RAD(Ndim),MWABC(i)%WDTH(Ndim),MWABC(i)%QN(Ndim))
            WRITE(STR,'(I0)')I
            MWABC(i)%name = 'MWABC_'//trim(STR)
            do j=1,MH_Str(i)%ncells - 1
              id1 = MH_Str(i)%cells(j)
              id2 = MH_Str(i)%cells(j+1)
              if(id2.eq.cellmap(5,id1) .or. id2.eq.cellmap(6,id1)) then
                MWABC(i)%sgn(j) = 1
                MWABC(i)%rad(j) = 2
                MWABC(i)%cell(j) = id1
                MWABC(i)%Q(j) = 0
                MWABC(i)%QN(J) = 0
                MWABC(i)%eta(j) = 0
                MWABC(i)%del(j) = (dy(id1)+dy(id2))/2.
                MWABC(i)%wdth(j) = min(dx(id1),dx(id2))
              endif
              if(id2.eq.cellmap(1,id1) .or. id2.eq.cellmap(2,id1)) then
                MWABC(i)%sgn(j) = -1
                MWABC(i)%rad(j) = 2
                MWABC(i)%cell(j) = id2
                MWABC(i)%Q(j) = 0
                MWABC(i)%QN(J) = 0
                MWABC(i)%eta(j) = 0
                MWABC(i)%del(j) = (dy(id1)+dy(id2))/2.
                MWABC(i)%wdth(j) = min(dx(id1),dx(id2))
              endif	
              if(id2.eq.cellmap(7,id1) .or. id2.eq.cellmap(8,id1)) then
                MWABC(i)%sgn(j) = 1
                MWABC(i)%rad(j) = 1
                MWABC(i)%cell(j) = id1
                MWABC(i)%Q(j) = 0
                MWABC(i)%QN(J) = 0
                MWABC(i)%eta(j) = 0
                MWABC(i)%del(j) = (dx(id1)+dx(id2))/2.
                MWABC(i)%wdth(j) = min(dy(id1),dy(id2))
              endif
              if(id2.eq.cellmap(3,id1) .or. id2.eq.cellmap(4,id1)) then
                MWABC(i)%sgn(j) = -1
                MWABC(i)%rad(j) = 1
                MWABC(i)%cell(j) = id2
                MWABC(i)%Q(j) = 0
                MWABC(i)%QN(J) = 0
                MWABC(i)%eta(j) = 0
                MWABC(i)%del(j) = (dx(id1)+dx(id2))/2.
                MWABC(i)%wdth(j) = min(dy(id1),dy(id2))
              endif
            enddo
            j=MH_Str(i)%Ncells
            MWABC(i)%sgn(j) = 0
            MWABC(i)%rad(j) = 0
            MWABC(i)%cell(j) = 0
            MWABC(i)%Q(j) = 0
            MWABC(i)%QN(J) = 0
            MWABC(i)%eta(j) = 0
            MWABC(i)%del(j) = 0
            MWABC(i)%wdth(j) = 0
          enddo !nmdriver 
        endif !H_multi

#ifdef DEBUG
        write(DGUNIT,*)'finished initalizing radiation stress'
#endif
      endif  ! radstr or wndsim or wind_spatial		 CR
      
      RETURN
      END SUBROUTINE INITIALIZE_WABC_tel
      
!********************************************************************************
      SUBROUTINE INITIALIZE_TRANSPORT_Tel
      use EXP_Global_def 
      USE EXP_bndcond_def
      USE EXP_transport_def
      use size_def    
      use flow_def
      use comvarbl
      use bnd_def
      use sal_def  
      use sed_def, only: rhosed,hardbed  
      use hot_def  
      use geo_def, only: dx,dy,idmap,zb,cell2cell
      use sed_def, only: sedtrans,ct,qtx,qty,rs
      use EXP_TELESCOPING
       
      IMPLICIT NONE
      INTEGER i,j,K,id_t
      real fac, dmax_E,dmax_D,dmin_E,dmin_D
    
    !IMPLICIT SOLVER ARRAY SAL IS USED FOR Salinity Concentration
    ! the explicit variables for diffusion and time averaged flux are retained
      if(saltrans) then
        allocate(SALT(ncellsD))
        salt%diffC=0
        salt%qx = 0
        salt%qy = 0
        do i=1,ncellsD
          salt(i)%vol = (eta(i)-zb(i))*dx(i)*dy(i)
        enddo
        if(.not. allocated(xTransQ) ) allocate(xTransQ(numxfaces))
        if(.not. allocated(yTransQ) ) allocate(yTransQ(numyfaces))       
        xTransQ = 0
        yTransQ = 0
      Endif      
      
    !sedtrans is implicit code switch for NET
    !the explicit NET is built on the ADEQ advection/diffsuion algorithm
    !and uses some of those varaible for processing - they are allocated
    !and initialized here
        if(sedtrans) then
          allocate(ADSS(ncellsD))
          ADSS%conc=0
          ADSS%concN=0
          ADSS%diffC=0
          ADSS%eros=0
          ADSS%depo=0
          ADSS%qx = 0
          ADSS%qy = 0
          do i=1,ncellsD
          ADSS(i)%vol = (eta(i)-zb(i))*dx(i)*dy(i)
          enddo
        if(.not. allocated(xSedTransQ) ) allocate(xSedTransQ(numxfaces))
        if(.not. allocated(ySedTransQ) ) allocate(ySedTransQ(numyfaces))  
        xSedTransQ = 0
        ySedTransQ = 0
          allocate(acoef(nmaxfaces,ncellsD))          !ncellsD,nmaxfaces switched 10/27/2015
        Endif

        
! the explicit code also supports AD (with Lund-CIRP) and TL (with Lund-CIRP and WATANABE)
! the required varaibles are allocated and initialized here
! some implicit variables are also allocated, since the explicit values are mapped to implicit varaibles for 
! writing to output files

      if(sedtransEXP) then  
          if(isedform == 3) then   !AD scheme
          allocate(qsx(NCellsD),qsy(NcellsD),bed(ncellsD))
          tsed_elapse = 0.0
          tmorph_elapse = 0.0
          qsx=0
          qsy=0
          bed=0
          if( .not. allocated(adss)) allocate(ADSS(ncellsD))
          ADSS%conc=0
          ADSS%concN=0
          ADSS%diffC=0
          ADSS%eros=0
          ADSS%depo=0
          ADSS%qx = 0
          ADSS%qy = 0
          do i=1,ncellsD
            ADSS(i)%vol = (eta(i)-zb(i))*dx(i)*dy(i)
          enddo
          if( .not. allocated(ct)) allocate (ct(ncellsD))
          if( .not. allocated(qtx)) allocate (qtx(ncellsD))
          if( .not. allocated(qty)) allocate (qty(ncellsD))         
          if( .not. allocated(rs)) allocate (rs(ncellsD)) 
          endif
          if(isedform == 1 .or. isedform == 2) then   !TL (Watanabe or Lund-CIRP)
          allocate(qsx(NCellsD),qsy(NcellsD),bed(ncellsD))
          tsed_elapse = 0.0
          tmorph_elapse = 0.0
          qsx=0
          qsy=0
          bed=0
          if( .not. allocated(qtx)) allocate (qtx(ncellsD))
          if( .not. allocated(qty)) allocate (qty(ncellsD))         
          endif          
      endif      

      
      ! cexplicit code cohesive algorithm initialization here
        if(cohesive) then
          cohesive_read = .false.
          !if(elev_read) cohesive_read = .true.
          allocate(COHES(ncellsD))
          COHES%conc=0
          COHES%concN=0
          COHES%diffC=0
          COHES%eros=0
          COHES%depo=0
          COHES%qx = 0
          COHES%qy = 0
          COHES%tbmax = 0
          do i=1,ncellsD
            COHES(i)%vol = (eta(i)-zb(i))*dx(i)*dy(i)
          enddo
          csfils(1) = 'cohesive.dat'
          csfils(2) = 'cohesive.dat'
          
          if(CHparms%Tcrit_variable) then
          !set E and D according to depth
          do i=1,ncellsD
          if(-zb(i).ge.CHparms%Depth(1)) then
          COHES(i)%Tcrit_E = CHparms%Tcrit_Ev(1)
          COHES(i)%Tcrit_D = CHparms%Tcrit_Dv(1)
          endif        
          do k=1,CHparms%numDEPTHS-1
          if(-zb(i).le.CHparms%Depth(k).and. -zb(i)  .gt.CHparms%Depth(k+1))then
          fac = (-zb(i)-CHparms%Depth(k+1))/(CHparms%Depth(k)-CHparms%Depth(k+1))
          COHES(i)%Tcrit_E = CHparms%Tcrit_Ev(k+1)+fac*(CHparms%Tcrit_Ev(k)-CHparms%Tcrit_Ev(k+1))
          COHES(i)%Tcrit_D = CHparms%Tcrit_Dv(k+1)+fac*(CHparms%Tcrit_Dv(k)-CHparms%Tcrit_Dv(k+1)) 
          endif
          enddo
          if(-zb(i).le.CHparms%Depth(CHparms%numDEPTHs)) then
          COHES(i)%Tcrit_E = CHparms%Tcrit_Ev(CHparms%numDEPTHS)
          COHES(i)%Tcrit_D = CHparms%Tcrit_Dv(CHparms%numDEPTHS)
          endif         
          enddo         
          else
          !set all value so same input value
          COHES%Tcrit_E = CHparms%Tcrit_E
          COHES%Tcrit_D = CHparms%Tcrit_D
          endif
          
          dmin_E = 9999
          dmin_D = 9999
          dmax_E = 0
          dmax_D = 0
          do i=1,ncellsD
          dmax_E = MAX(cohes(i)%TCriT_E,dmax_E)
          DMIN_e = MIN(COHES(I)%Tcrit_E,dmin_E)
           dmax_D = MAX(cohes(i)%TCriT_D,dmax_D)
          DMIN_D = MIN(COHES(I)%Tcrit_D,dmin_D) 
          enddo
          write(*,*)dmax_E,dmin_E
          write(*,*)dmax_D,dmin_D                  
          
        Endif 
       
        if(nQstr .eq. 0) cohes_flow_bc = .false.  !no need for cohes bcs if no flow bcs

        if(cohes_flow_bc) then
          maxunit=maxunit+1
          cohes_unit = maxunit
          open(unit=cohes_unit,file=cohes_flow_bc_file,status='old')

          read(cohes_unit,*)!header
          allocate (cohes_bc(2,nQstr),cohes_bc_time(2))
          read(cohes_unit,*)cohes_bc_time(1),(cohes_bc(1,j),j=1,nQstr)
          read(cohes_unit,*)cohes_bc_time(2),(cohes_bc(2,j),j=1,nQstr)
          do i=1,2
            do j=1,nQstr
              cohes_bc(i,j) = cohes_bc(i,j)/(rhosed*1000)  !convert from mg/L to volume concentration
            enddo
          enddo
          !determine cells that require cohes bcs (piigy back of of flow bcs)
          !first allocate cell stoarge array - need max length of flow bc strings
          ndummy=0
          do i=1,nQstr
            ndummy = max(ndummy,Q_Str(i)%ncells)
          enddo
          allocate(cohes_bc_cells(nQstr,ndummy))
          !now assign cell values
          do i=1,nQstr
            if(QstringEXP(I)%Vface .eqv. .true.) then !then on north or south face  	
              do j=1,Q_Str(I)%NCELLS
                ID_T = Q_Str(I)%cells(j)
                if(ID_T.gt.ncells) cohes_bc_cells(i,j)= ID_T  !on the north face
                if(ID_T.le.ncells) cohes_bc_cells(i,j)= cell2cell(3,ID_T) !on the south face
             enddo
           endif
           if(QstringEXP(I)%Vface .eqv. .false.) then !then on east or west face  	
             do j=1,Q_Str(I)%NCELLS
               ID_T = Q_Str(I)%cells(j)
               if(ID_T.gt.ncells) cohes_bc_cells(i,j)= ID_T  !on the east face
               if(ID_T.le.ncells) cohes_bc_cells(i,j)= cell2cell(4,ID_T) !on the west face
              enddo
            endif	 
          enddo      
        endif !end cohes_flow_bc
      
      RETURN
      END SUBROUTINE !INITIALIZE TRANSPORT

