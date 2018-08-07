
module tools
implicit none

contains

!
function thing(x)
implicit none

integer, intent(in) :: x
integer :: thing         
        thing=2*x
end function

!
subroutine doubleit(x)
implicit none 
integer, intent(inout) :: x
        x=2*x
end subroutine doubleit

!
end module tools


