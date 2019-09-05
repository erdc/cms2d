      subroutine CMS_FLOW_EXP_GridType()

      !Program selects routines based on gridtype       
  
    use geo_def, only: igridtype
      implicit none
      
!      write(*,*)'igridtype = ',igridtype
      if(igridtype .eq. 0) then
!           Write(*,*)'      using explicit cartesian grid solver'         
          call CMS_FLOW_EXP_C
      elseif(igridtype .eq. 1) then
!          Write(*,*)'      using explicit telescoping grid solver'
          call CmS_FLOW_EXP_T
      elseif(igridtype .gt. 1) then
          write(*,*)'Grid Type not supported by the explicit solver'
          write(*,*)'Program terminated'
          write(*,*)         
      endif

      
      END SUBROUTINE 
 
