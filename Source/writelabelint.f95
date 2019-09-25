!	write a label and an integer nvin csv format to unit.
subroutine writelabelint(label,nv,unit)
implicit none
character(len=*),intent(in)::label
integer,intent(in)::nv,unit
!	buff is used to process nv.
!	io is a flag for a writing error.
character(len=12)::buff
integer::io
write(buff,'(i12)')nv
write(unit,'(3a)',iostat=io)trim(adjustl(label)),",",trim(adjustl(buff))
if(io/=0)stop "Output failure"
return
end subroutine writelabelint
