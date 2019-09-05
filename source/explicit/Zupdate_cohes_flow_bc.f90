
      subroutine update_cohes_flow_bc()
    use EXP_Global_def 
      USE EXP_bndcond_def
      use bnd_def
      USE EXP_transport_def 
      use flow_def
    use comvarbl, only: timehrs      
      use sed_def, only: rhosed 
      use geo_def, only: dx,dy,cell2cell
          
      implicit none  
      real(ikind) fac,value
      integer i,j,ii,id,ido   

      do while (timeHRS.ge.cohes_bc_time(2))
        cohes_bc_time(1) = cohes_bc_time(2)
        do i=1,nQstr
          cohes_bc(1,i) = cohes_bc(2,i)
        enddo
        read(cohes_unit,*)cohes_bc_time(2),(cohes_bc(2,j),j=1,nQstr)
        do j=1,nQstr
          cohes_bc(2,j) = cohes_bc(2,j)/(rhosed*1000)  !convert from mg/L to volume concentration
        enddo     
      enddo

      fac=(timeHRS-cohes_bc_time(1))/(cohes_bc_time(2)-cohes_bc_time(1))
      do i=1,nQstr
        value = cohes_bc(1,i) + fac*(cohes_bc(2,i)-cohes_bc(1,i))
        do j=1,Q_str(i)%ncells
          id = cohes_bc_cells(i,j)
          cohes(id)%conc = value
        enddo
      enddo
        
      !calculate FUU and GUU array for use in conc transport calcs (doen next)
      !some of these may be re-calculated in conc trans routines - but not all of them
      do j = 1,nQstr  !for each cell string
        if(QstringEXP(j)%vface) then
          IDO = Q_str(j)%NCells
          do ii=1,IDO    !calculate GVV
            i=Q_str(j)%cells(ii)    
            if(COHES(i)%qy.gt.0) then
              ncs = cell2cell(3,i)
              Gvv(i) = COHES(i)%qy*COHES(ncs)%Conc*dx(i)
            else
              Gvv(i) = COHES(i)%qy*COHES(i)%Conc*dx(i)
            endif    
          enddo
        else
          IDO = Q_str(j)%NCells
          do ii=1,IDO    !calculate FUU
            i=Q_str(j)%cells(ii)
            if(COHES(i)%qx.gt.0) then
              ncw = cell2cell(4,i)
              Fuu(i) = COHES(i)%qx*COHES(ncw)%Conc*dy(i)
            else
              Fuu(i) = COHES(i)%qx*COHES(i)%Conc*dy(i)
            endif
          enddo
        endif
      enddo ! end of NQdriver   
      
      end subroutine
