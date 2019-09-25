!	Print gradients.
!	id is the observation id.
!	paramname is the array of parameter names.
!	unitgrad is the output unit.
!	useid indicates if individual identifications are present in id.
!	gradsdes is the array of gradients.	
subroutine outputgrad(id,paramname,unitgrad,useid,gradsdes)
implicit none
character(len=*),intent(in)::id(:)
character(len=*),intent(in)::paramname(:)
integer::unitgrad
logical,intent(in)::useid
real(kind=8),intent(in)::gradsdes(:,:)

!	Buffers for writing.
character(len=4)::buffer
character(len=32)::bufferid
character(len=25),allocatable::buffm(:)
character(len=7)::writefmt
!	col is column number.
!	io is error flag.
!	obs is observation number.

integer::col,io,obs
col=size(paramname)


allocate(buffm(col),stat=io)
if(io/=0) stop "Allocation failure for arrays for gradient printing."
write(unitgrad,'(a)') "Gradients"
write(buffer,'(i4)') 1+2*col
writefmt='('//trim(adjustl(buffer))//'a)'
write(unit=unitgrad,fmt=writefmt,iostat=io) "ID",(",",""//trim(adjustl(paramname(col))),col=1,size(paramname))
if(io/=0) stop "Writing of gradients failed."
	
do obs=1,size(gradsdes,2)
	do col=1,size(paramname)
		write(buffm(col),'(g25.16e3)')gradsdes(col,obs)
	end do
	if(useid)then
		bufferid=id(obs)
	else
		write(bufferid,'(i12)') obs
	end if
	write(unit=unitgrad,fmt=writefmt,iostat=io) trim(adjustl(bufferid)),(",",trim(adjustl(buffm(col))),col=1,size(paramname))
	
	if(io/=0) stop "Writing of EAPs failed."	
end do
return
end subroutine outputgrad
