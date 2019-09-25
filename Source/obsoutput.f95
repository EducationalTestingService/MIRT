!	print observed scale score.
!	comment is used for identification.
!	factorname gives names of elements of the vector.
!	id is individual id.
!	unitobs unitobs is the unit for output.
!	obsid indicates if ids are put out.
!   presence indicates what is missing.
!	observed gives observations.

subroutine obsoutput(comment,factorname,id,unitobs,obsid,presence,observed)
implicit none
character(len=*),intent(in)::comment
character(len=*),intent(in)::factorname(:)
character(len=*),intent(in),optional::id(:)
integer,intent(in)::unitobs
logical,intent(in)::obsid,presence(:)
real(kind=8),intent(in)::observed(:,:)
!	Buffers for writing.
character(len=4)::buffer
character(len=32)::bufferid
character(len=25),allocatable::buffm(:)
!	format
character(len=7)::writefmt
!	col is column number.
!	io is error flag.
!	obs is observation number.
!	row is row number.	
integer::col,io,obs,row
col=size(factorname)


allocate(buffm(col+1),stat=io)
if(io/=0) stop "Allocation failure for arrays for output of observations."
write(unitobs,'(a)') comment
write(buffer,'(i4)') 2*(2*col)
writefmt='('//trim(adjustl(buffer))//'a)'
write(unit=unitobs,fmt=writefmt,iostat=io) "ID",",","Presence",",",(trim(adjustl(factorname(row))),",",row=1,size(factorname))
if(io/=0) stop "Writing of conditional means and covariances failed."
	
do obs=1,size(presence)
    write(buffm(1),'(L1)')presence(obs)
	do row=1,size(factorname)
		write(buffm(row+1),'(g25.16e3)')observed(row,obs)
	end do
	if(obsid)then
		bufferid=id(obs)
	else
		write(bufferid,'(i12)') obs
	end if
	write(unit=unitobs,fmt=writefmt,iostat=io) trim(adjustl(bufferid)),",",(trim(adjustl(buffm(row))),",",row=1,size(factorname)+1)
	
	if(io/=0) stop "Writing of observed scores failed."
end do
return



end subroutine obsoutput
