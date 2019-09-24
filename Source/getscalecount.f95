!
!	Obtain counts numberscales of scale scores per weighted sum.
subroutine getscalecount(numberscales)
implicit none
integer,intent(out)::numberscales
!	io is the flag for input errors.
integer::io
!	numberscales is the number of scales.
namelist/scalecounts/numberscales
!	Default is 0.
numberscales=0
read(*,nml=scalecounts,iostat=io)
if(io/=0) stop 'Unable to read number of scale scores desired.'
!	Force counts to be at least 0.
numberscales=max(0,numberscales)
return
end subroutine getscalecount	