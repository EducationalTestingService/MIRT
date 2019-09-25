!	write a label and a floating point number x in csv format to unit.
subroutine writelabelfl(label,unit,x)
implicit none
character(len=*),intent(in)::label
integer,intent(in)::unit
real(kind=8),intent(in)::x
!	buff is used to process nv.
!	io is a flag for a writing error.
character(len=25)::buff
integer::io
write(buff,'(g25.16e3)')x
write(unit,'(3a)',iostat=io)trim(adjustl(label)),",",trim(adjustl(buff))
if(io/=0)stop "Output failure"
return
end subroutine writelabelfl
