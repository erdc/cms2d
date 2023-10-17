!***********************************************************************
      subroutine ST_hardbottom_COHES()
!***********************************************************************
      use EXP_Global_def,    only: hardbottom
      USE EXP_transport_def, only: bed, cohes
      use sed_def, only: poros,hardbed,scalemorph,nhard,idhard
      use geo_def, only: zb 
           
      implicit none
      integer i,id
      real tdepth 

      !check for hardbottom conditions and adjust qs's as required
      if(hardbottom) then
        do i=1,nhard
          id = idhard(i)
          tdepth = -zb(id)- bed(id)*scalemorph*POROS
          if(tdepth.gt.hardbed(i)) then
            COHES(id)%eros = 0
          endif
        enddo
      endif    

      return
      end subroutine
