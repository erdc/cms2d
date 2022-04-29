!===================================================================
module flow_wavegrid_def
! Flow variables on the wave grid.
! Used to pass information from the flow model to the wave model
!===================================================================
    use prec_def
    implicit none
    save

    real(ikind),allocatable:: xwavfl(:,:),ywavfl(:,:) !Wave grid on flow coordinate system
    real(ikind),allocatable:: depwave(:,:)            !Still water depth on wave grid [m], positive is wet
    real(ikind),allocatable:: depwave0(:,:)           !Initial still water depth on wave grid
    real(ikind),allocatable:: etawave(:,:)            !Mean water level on wave grid [m]
    real(ikind),allocatable:: hwave(:,:)              !Total water depth on wave grid [m],hwave=depwave+etawave
    real(ikind),allocatable:: uwave(:,:),vwave(:,:)   !Current velocities on the wave grid [m/s], with respect to local coordinate system    
    real(ikind),allocatable:: coefintp_flwav(:,:,:)   !Interpolation coefficients from flow to wave grid
    integer,    allocatable:: iiflcell(:,:,:)         !Interpolation indeces from flow to wave grid        
    character(len=200) :: flwav_intpcoef_file
    
end module flow_wavegrid_def 

!===================================================================
module wave_wavegrid_def
! Wave variables on the wave grid. 
! Used to pass information from the wave model to the flow model 
!===================================================================
    use prec_def
    implicit none
    save
    
    integer :: idatewave
    integer :: nwavei,nwavej,nwaveij                  !Wave grid size [-]
    real(ikind) :: xwav0,ywav0                        !Wave grid origin [m]
    real(ikind) :: azimuth_wav                        !Wave grid azimuth [deg]
    real(ikind),allocatable:: xwave(:),ywave(:)       !Wave coordinates on the wave grid (local coordinate sytem)
    real(ikind),allocatable:: dxwav(:),dywav(:)       !Wave grid resolution [m]
    real(ikind),allocatable:: wxrs1(:,:),wyrs1(:,:)   !Wave radiation stress gradients [m^2/s^2]
    real(ikind),allocatable:: wheight(:,:)            !Wave height [m]
    real(ikind),allocatable:: wperiod(:,:)            !Wave peak period, [sec]
!Lihwa removed - test   12/18/2020
    !real(ikind),allocatable:: wibr(:,:)               !Wave breaking index, 0-no breaking, 1-Breaking
    real(ikind),allocatable:: wdiss(:,:)              !Wave dissipation [N/m/s]
    real(ikind),allocatable:: wcos(:,:),wsin(:,:)     !Wave unit vectors [-]
    
end module wave_wavegrid_def

!===================================================================
module wave_flowgrid_def
!===================================================================    
    use prec_def
    implicit none
    save
    
    integer :: IDate0
    
    real(ikind):: tswave1,tswave2,tide1,tide2                        !Times [hrs] 
    real(ikind),allocatable:: xflwav(:),yflwav(:)                    !Flow grid on wave coordinate system [m]
    real(ikind),allocatable:: Whgt(:),Whgt1(:),Whgt2(:)              !Significant Wave height [m]
    real(ikind),allocatable:: Wper(:),Wper1(:),Wper2(:)              !Peak Wave period [sec]
    real(ikind),allocatable:: wavediss(:),wavediss1(:),wavediss2(:)  !Wave dissipation [N/m/s]
!Lihwa removed - test   12/18/2020
    !real(ikind),allocatable:: waveibr(:),waveibr1(:),waveibr2(:)     !Breaking index (0 or 1)
    real(ikind),allocatable:: wavestrx(:),wavestrx1(:),wavestrx2(:)  !Wave radiation stress gradient in x-direction [m^2/s^2]
    real(ikind),allocatable:: wavestry(:),wavestry1(:),wavestry2(:)  !Wave radiation stress gradient in y-direction [m^2/s^2]
    real(ikind),allocatable:: Wunitx(:),Wunitx1(:),Wunitx2(:)        !Wave unit vector in x-direction [-]
    real(ikind),allocatable:: Wunity(:),Wunity1(:),Wunity2(:)        !Wave unit vector in x-direction [-]
    real(ikind),allocatable:: Wang(:)                                !Wave angle, going to  [rad]
    real(ikind),allocatable:: Wlen(:)                                !Wavelength correponding to peak wave period [m]
    real(ikind),allocatable:: Worb(:)                                !Significant bottom orbital velocity amplitude [m]
    real(ikind),allocatable:: Worbrep(:)                             !Root-mean-squared error bottom orbital velocity amplitude [m]
    real(ikind),allocatable:: wetsteer(:)                            !Used for steering
    real(ikind),allocatable:: Ssr(:)                                 !Surface roller energy x 2 (Ssr=2*Essr) [N/m]
    real(ikind),allocatable:: coefintp_wavfl(:,:)                    !Wave-to-flow Interpolation coefficients
    real(ikind),allocatable:: ueff(:),veff(:)                        !Effective current velocities used in wave model [m/s]
    integer,allocatable:: ijwavcell(:,:) 
    character(len=200) :: wavfl_intpcoef_file
    
    !For testing
    real(ikind),allocatable:: wavstrx(:),wavstry(:) !**************************************
    
    !Constant wave parameters for testing and idealized cases
    logical :: constant_waves
    real(ikind) :: waveheight, waveperiod, wavedir 
    
end module wave_flowgrid_def 