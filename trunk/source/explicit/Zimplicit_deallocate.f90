!***********************************************************************
    subroutine Zimplicit_var_deallocate()
!
!***********************************************************************  
    use size_def
    use geo_def
    use flow_def
    use struct_def
    use met_def
    use solv_def
    use comvarbl
	use cms_def
    use flow_wavegrid_def
	use wave_flowgrid_def
    use bnd_def
    use sed_def
    use sal_def
    use fric_def
    use stat_def
    use veg_def
    use hot_def
    use der_def
    use out_def
    use q3d_def
    !use der_lib, only: der_grad_eval
    !use intp_lib, only: interp_coef_cell2face
    use const_def   
!    use FLOW_WAVE_INTERP, only: wavex,wavey     

    !now deallocate implicit arrays that are not needed
    !if (allocated(ncface)) deallocate(ncface)  
    !if (allocated(loconect)) deallocate(loconect)
    !if (allocated(idirface)) deallocate(idirface)
    !if (allocated(idjface)) deallocate(idjface)          
    if (allocated(nxface)) deallocate(nxface)
    if (allocated(nyface)) deallocate(nyface)
    if (allocated(nxyface)) deallocate(nxyface)
    if (allocated(kxface)) deallocate(kxface)
    if (allocated(kyface)) deallocate(kyface)
    if (allocated(kxyface)) deallocate(kxyface)    
    !State
    !if (allocated(iwet))    deallocate(iwet)
    if (allocated(iwet1))    deallocate(iwet1)
    if (allocated(iextrap))  deallocate(iextrap)
    !if (allocated(icorner)) deallocate(icorner)            
    
    !Interpolation
    !if (allocated(fintp))   deallocate(fintp)   
    
    !Geometry
    !if (allocated(volp))    deallocate(volp)   
    if (allocated(dc))       deallocate(dc)
    if (allocated(ds))       deallocate(ds)
    if (allocated(dsxy))     deallocate(dsxy) 
    !if (allocated(dnorm))   deallocate(dnorm)
    !if (allocated(dpara))   deallocate(dpara)  
    
    !Shared implicit and explicit flow variables
    !if (allocated(cdflux))  deallocate(cdflux)
    if (allocated(zb))       deallocate(zb)
    !if (allocated(uv))      deallocate(uv)    
    if (allocated(p))        deallocate(p)
    if (allocated(p1))       deallocate(p1)
    !if (allocated(piks))    deallocate(piks)
    !if (allocated(dpx))     deallocate(dpx)
    !if (allocated(dpy))     deallocate(dpy)  
	if (allocated(h1))       deallocate(h1)
	!if (allocated(hks))     deallocate(hks)
	if (allocated(dhx))      deallocate(dhx)
	if (allocated(dhy))      deallocate(dhy)          
    if (allocated(u1))       deallocate(u1)
    !if (allocated(uiks))    deallocate(uiks)
    !if (allocated(dux))     deallocate(dux)
    !if (allocated(duy))     deallocate(duy)	
    if (allocated(v1))       deallocate(v1)
    !if (allocated(vik))     deallocate(vik)
    !if (allocated(dvx))     deallocate(dvx)
    !if (allocated(dvy))     deallocate(dvy)	            
    
	!Implicit solver    
	if (allocated(acoef))    deallocate(acoef)
    if (allocated(pp))       deallocate(pp)
    !if (allocated(ppk))     deallocate(ppk)
    if (allocated(dppx))     deallocate(dppx)
    !if (allocated(apuvolp)) deallocate(apuvolp)   !Not allocatable variable
    if (allocated(sumu))     deallocate(sumu)      
    if (allocated(su))       deallocate(su)
    if (allocated(sp))       deallocate(sp) 
    !if (allocated(suu0))    deallocate(suu0)
    !if (allocated(suv0))    deallocate(suv0)
    !if (allocated(supp0))   deallocate(supp0) 
    !if (allocated(spu0))    deallocate(spu0)
    !if (allocated(spv0))    deallocate(spv0)
    if (allocated(sppp0))    deallocate(sppp0)
    if (allocated(Hu))       deallocate(Hu)
    !if (allocated(Hv))      deallocate(Hv)  
    !if (allocated(rsp))     deallocate(rsp)
    !if (allocated(rsu))     deallocate(rsu)
    !if (allocated(rsv))     deallocate(rsv)       
	if (allocated(u2))       deallocate(u2)
	if (allocated(v2))       deallocate(v2)
	if (allocated(p2))       deallocate(p2)
	if (allocated(h2))       deallocate(h2)
    
    !Total Fluxes for u and v equations
	if (allocated(Huk))      deallocate(Huk)
	if (allocated(Hvk))      deallocate(Hvk)        
    
    !diffusion
	!if (allocated(vis))     deallocate(vis)
	if (allocated(viskfl))   deallocate(viskfl)
	!if (allocated(visk))    deallocate(visk)

    return
    end subroutine Zimplicit_var_deallocate
