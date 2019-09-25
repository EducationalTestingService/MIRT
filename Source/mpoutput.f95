!	print maximum posterior likelihood data for latent vector.
!	factorname gives names of elements of the latent vector.
!	id is individual id.
!	unitmp is the unit for output.
!	eapid indicates if ids are put out.
!	alpha are locations of maxima.
!	cholnhess is the array of modified Cholesky decompositions.
subroutine mpoutput(factorname,id,unitmp,eapid,alpha,cholnhess)
implicit none
interface
!	Compute inverse for n by n matrix with modified Cholesky
!	decomposition tri.
	function invert(tri)
		implicit none
		real(kind=8),intent(in)::tri(:,:)
		real(kind=8)::invert(size(tri,1),size(tri,1))
	end function invert
end interface
character(len=32),intent(in)::factorname(:)
character(len=32),intent(in),optional::id(:)
integer,intent(in)::unitmp
logical,intent(in)::eapid
real(kind=8),intent(in)::alpha(:,:),cholnhess(:,:,:)

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
!	Inverse location is inv.
real(kind=8)::inv(size(alpha,1),size(alpha,1))

col=size(factorname)


allocate(buffm(col),buffc(col,col),stat=io)
if(io/=0) stop "Allocation failure for arrays for maximum posterior printing."
write(unitmp,'(a)') "Maximum posterior likelihood estimate and inverse information matrix"
write(buffer,'(i4)') 1+2*col*(col+1)
writefmt='('//trim(adjustl(buffer))//'a)'
write(unit=unitmp,fmt=writefmt,iostat=io) "ID",(",","Estimate_"//trim(adjustl(factorname(row))),row=1,size(factorname)),&
	((",","Inverse_"//trim(adjustl(factorname(row)))//"_"//trim(adjustl(factorname(col))),&
	col=1,size(factorname)),row=1,size(factorname))
if(io/=0) stop "Writing of maximum posterior likelihood estimates failed."
	
do obs=1,size(alpha,2)
	inv=invert(cholnhess(:,:,obs))
	do row=1,size(factorname)
		write(buffm(row),'(g25.16e3)')alpha(row,obs)
		do col=1,size(factorname)
			write(buffc(row,col),'(g25.16e3)')inv(row,col)
		end do
	end do
	if(eapid)then
		bufferid=id(obs)
	else
		write(bufferid,'(i12)') obs
	end if
	write(unit=unitmp,fmt=writefmt,iostat=io) trim(adjustl(bufferid)),(",",trim(adjustl(buffm(row))),row=1,size(factorname)),&
			((",",trim(adjustl(buffc(row,col))),col=1,size(factorname)),row=1,size(factorname))
	
	if(io/=0) stop "Writing of maximum posterior likelihood estimates failed."	
end do
return



end subroutine mpoutput
