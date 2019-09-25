!	title is written on unitno.
subroutine outputtitle(title,unitno)
implicit none
character(len=80),intent(in)::title
integer,intent(in)::unitno
!	io error flag
integer::io
write(unitno,'(a)',iostat=io) trim(adjustl(title)) 
if(io/=0)stop "Output failure"
end subroutine outputtitle
