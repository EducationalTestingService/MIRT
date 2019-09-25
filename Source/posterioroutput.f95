!	print posterior distribution data for latent vector.
!	factorname gives names of elements of the latent vector.
!	id is ids.
!	Output unit is unitpost.
!	eapid indicates if ids are put out.
!	postdensity are weights.
!	theta are adaptive quadrature points.  
subroutine posterioroutput(factorname,id,unitpost,eapid,postdensity,theta)
implicit none
character(len=32),intent(in),optional::id(:)
character(len=32),intent(in)::factorname(:)

integer,intent(in)::unitpost

logical,intent(in)::eapid
real(kind=8),intent(in)::postdensity(:,:),theta(:,:,:)
!	buffers for writing.
character(len=12)::buff
character(len=32)::buff1
character(len=25)::buffer(size(postdensity,1)*(size(factorname)+1))
character(len=17)::writefmt
!	counter for buffer.
!	Dimension number.
!	io error flag.
!	Observation number.
!	Number of quadrature point.
integer::counter,dimno,io,obs,quad
write(unitpost,'(a)') "Posterior distributions"
!	Number quadrature points
write(buff,'(i12)') 2*size(postdensity,1)*(size(factorname)+1)
writefmt='(a,'//trim(adjustl(buff))//'a)'
do quad=1,size(postdensity,1)
	write(buffer(quad),'(i4)') quad
end do
write(unit=unitpost,fmt=writefmt,iostat=io) "ID",&
	((",",trim(adjustl(factorname(dimno)))//"_"//trim(adjustl(buffer(quad))),&
			dimno=1,size(factorname)),&
			",","Weight_"//trim(adjustl(buffer(quad))),&
			quad=1,size(postdensity,1))
if(io/=0) stop "Failure to write posterior information."			
do obs=1,size(postdensity,2)
	if(eapid)then
		buff1=id(obs)
	else
		write(buff,'(i12)') obs
		buff1=trim(adjustl(buff))
	end if
	counter=1
	do quad=1,size(postdensity,1)
		do dimno=1,size(factorname)
			write(buffer(counter),'(g25.16)') theta(dimno,quad,obs)
			counter=counter+1
		end do
		write(buffer(counter),'(g25.16)') postdensity(quad,obs)
		counter=counter+1
	end do
	write(unit=unitpost,fmt=writefmt,iostat=io) buff1,&
		(",",trim(adjustl(buffer(counter))),counter=1,size(buffer))
	if(io/=0) stop "Failure to write posterior information."
end do
return
end subroutine posterioroutput
