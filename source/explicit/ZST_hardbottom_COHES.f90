      subroutine ST_hardbottom_COHES()
    use EXP_Global_def 
      USE EXP_bndcond_def
      USE EXP_transport_def    
      use sed_def, only: poros,hardbed,scalemorph,nhard,idhard
      use flow_def
      use comvarbl 
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
      
      end subroutine
