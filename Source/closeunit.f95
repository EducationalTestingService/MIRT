!   See if unit to be closed.
subroutine closeunit(unit)

use unitdef
implicit none
integer,intent(in)::unit
unitcount(unit)=unitcount(unit)-1
if(unitcount(unit)==0)close(unit)
end subroutine closeunit

