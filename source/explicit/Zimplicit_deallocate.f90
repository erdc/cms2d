!***********************************************************************
    subroutine Zimplicit_var_deallocate()
!***********************************************************************  
    !use size_def
    use geo_def,  only: nxface, nyface, kxface, kyface, dc, ds
    use flow_def, only: iwet1, iextrap, p, p1, h1, dhx, dhy, u1, v1, pp, dppx
    use flow_def, only: sumu, su, sp, sppp0, hu, u2, v2, p2, h2, huk, hvk, viskfl

    !now deallocate implicit arrays that are not needed
    if (allocated(nxface)) deallocate(nxface)
    if (allocated(nyface)) deallocate(nyface)
    if (allocated(kxface)) deallocate(kxface)
    if (allocated(kyface)) deallocate(kyface)
    !if (allocated(nxyface)) deallocate(nxyface)      !This is needed MEB 07/11/23
    !if (allocated(kxyface)) deallocate(kxyface)      !This is needed MEB 07/11/23
    !if (allocated(ncface)) deallocate(ncface)  
    !if (allocated(loconect)) deallocate(loconect)
    !if (allocated(idirface)) deallocate(idirface)
    !if (allocated(idjface)) deallocate(idjface)          

    !State
    if (allocated(iwet1))    deallocate(iwet1)
    if (allocated(iextrap))  deallocate(iextrap)
    !if (allocated(iwet))    deallocate(iwet)
    !if (allocated(icorner)) deallocate(icorner)            
    
    !Interpolation
    !if (allocated(fintp))   deallocate(fintp)   
    
    !Geometry
    if (allocated(dc))       deallocate(dc)
    if (allocated(ds))       deallocate(ds)
    !if (allocated(dsxy))     deallocate(dsxy)        !This is needed MEB 07/11/23
    !if (allocated(volp))     deallocate(volp)   
    !if (allocated(dnorm))    deallocate(dnorm)
    !if (allocated(dpara))    deallocate(dpara)  
    
    !Shared implicit and explicit flow variables
    if (allocated(p))        deallocate(p)
    if (allocated(p1))       deallocate(p1)
    if (allocated(h1))       deallocate(h1)
    if (allocated(dhx))      deallocate(dhx)
    if (allocated(dhy))      deallocate(dhy)          
    if (allocated(u1))       deallocate(u1)
    if (allocated(v1))       deallocate(v1)
    !if (allocated(zb))       deallocate(zb)          !This is needed MEB 07/11/23
    !if (allocated(cdflux))   deallocate(cdflux)
    !if (allocated(uv))       deallocate(uv)    
    !if (allocated(piks))     deallocate(piks)
    !if (allocated(dpx))      deallocate(dpx)
    !if (allocated(dpy))      deallocate(dpy)  
    !if (allocated(hks))      deallocate(hks)
    !if (allocated(uiks))     deallocate(uiks)
    !if (allocated(dux))      deallocate(dux)
    !if (allocated(duy))      deallocate(duy)    
    !if (allocated(vik))      deallocate(vik)
    !if (allocated(dvx))      deallocate(dvx)
    !if (allocated(dvy))      deallocate(dvy)                
    
    !Implicit solver    
    if (allocated(pp))       deallocate(pp)
    if (allocated(dppx))     deallocate(dppx)
    if (allocated(sumu))     deallocate(sumu)      
    if (allocated(su))       deallocate(su)
    if (allocated(sp))       deallocate(sp) 
    if (allocated(sppp0))    deallocate(sppp0)
    if (allocated(Hu))       deallocate(Hu)
    if (allocated(u2))       deallocate(u2)
    if (allocated(v2))       deallocate(v2)
    if (allocated(p2))       deallocate(p2)
    if (allocated(h2))       deallocate(h2)
    !if (allocated(acoef))    deallocate(acoef)       !This is needed MEB 07/11/23
    !if (allocated(ppk))      deallocate(ppk)
    !if (allocated(apuvolp))  deallocate(apuvolp)   !Not allocatable variable
    !if (allocated(suu0))     deallocate(suu0)
    !if (allocated(suv0))     deallocate(suv0)
    !if (allocated(supp0))    deallocate(supp0) 
    !if (allocated(spu0))     deallocate(spu0)
    !if (allocated(spv0))     deallocate(spv0)
    !if (allocated(Hv))       deallocate(Hv)  
    !if (allocated(rsp))      deallocate(rsp)
    !if (allocated(rsu))      deallocate(rsu)
    !if (allocated(rsv))      deallocate(rsv)       
    
    !Total Fluxes for u and v equations
    if (allocated(Huk))      deallocate(Huk)
    if (allocated(Hvk))      deallocate(Hvk)        
    
    !diffusion
    if (allocated(viskfl))   deallocate(viskfl)
    !if (allocated(vis))      deallocate(vis)
    !if (allocated(visk))     deallocate(visk)

    return
    end subroutine Zimplicit_var_deallocate
