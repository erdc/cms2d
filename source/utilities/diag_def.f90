!===================================================================================
module diag_def
!===================================================================================
    implicit none
    
    !Diagnostic file and variables
    integer :: dgunit
    character(len=200) :: dgfile  !Diagnostic file    
    logical :: debug_mode  !Turns on or off the debug mode    
    integer :: debug_level !Placeholder. Not implemented yet
    
    character(len=200) :: msg,msg2,msg3,msg4,msg5,msg6,msg7,msg8,msg9,msg10,msg11,msg12,msg13
    
end module diag_def    
    