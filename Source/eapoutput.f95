!	print eap data for latent vector.
!	comment is used for identification.
!	factorname gives names of elements of the latent vector.
!	id is individual id.
!	unit uniteap is the unit for output.
!	eapid indicates if ids are put out.
!	covtheta are conditional covariance matrices.
!	meantheta are eaps. 

subroutine eapoutput(comment,factorname,id,uniteap,eapid,covtheta,meantheta)
implicit none
character(len=*),intent(in)::comment
character(len=32),intent(in)::factorname(:)
character(len=32),intent(in),optional::id(:)
integer,intent(in)::uniteap
logical,intent(in)::eapid
real(kind=8),intent(in)::covtheta(:,:,:),meantheta(:,:)

!	Buffers for writing.
character(len=4)::buffer
character(len=32)::bufferid
character(len=25),allocatable::buffm(:),buffc(:,:)
!	format
character(len=7)::writefmt
!	col is column number.
!	io is error flag.
!	obs is observation number.
!	row is row number.	
integer::col,io,obs,row
col=size(factorname)


allocate(buffm(col),buffc(col,col),stat=io)
if(io/=0) stop "Allocation failure for arrays for output of conditional means and covariances."
write(uniteap,'(2a)') "Conditional expectations and covariances ",comment
write(buffer,'(i4)') 1+2*col*(col+1)
writefmt='('//trim(adjustl(buffer))//'a)'
write(unit=uniteap,fmt=writefmt,iostat=io) "ID",(",","Mean_"//trim(adjustl(factorname(row))),row=1,size(factorname)),&
	((",","Cov_"//trim(adjustl(factorname(row)))//"_"//trim(adjustl(factorname(col))),&
	col=1,size(factorname)),row=1,size(factorname))
if(io/=0) stop "Writing of conditional means and covariances failed."
	
do obs=1,size(meantheta,2)
	do row=1,size(factorname)
		write(buffm(row),'(g25.16e3)')meantheta(row,obs)
		do col=1,size(factorname)
			write(buffc(row,col),'(g25.16e3)')covtheta(row,col,obs)
		end do
	end do
	if(eapid)then
		bufferid=id(obs)
	else
		write(bufferid,'(i12)') obs
	end if
	write(unit=uniteap,fmt=writefmt,iostat=io) trim(adjustl(bufferid)),(",",trim(adjustl(buffm(row))),row=1,size(factorname)),&
			((",",trim(adjustl(buffc(row,col))),col=1,size(factorname)),row=1,size(factorname))
	
	if(io/=0) stop "Writing of conditional means and covariances failed."	
end do
return



end subroutine eapoutput
