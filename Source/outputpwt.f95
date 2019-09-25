!	print lower and upper probabilities.
!	id is individual id. prob for model.
!   weightname is name of weight.
!	unitpwt is output unit.
!	eapid indicates if ids are put out.
!	plower is array of lower probabilities.
!   pupper is array of upper probabilities.
subroutine outputpwt(id,weightname,unitpwt,eapid,plower,pupper)
implicit none

character(len=32),intent(in)::id(:)
character(len=32),intent(in)::weightname
integer,intent(in)::unitpwt
logical,intent(in)::eapid
real(kind=8),intent(in)::plower(:),pupper(:)
!	buffers.
character(len=32)::bufferid
character(len=25)::buff,buff1
!	Error flag.  Observation counter
integer::io,obs
write(unitpwt,'(a)',iostat=io) "ID,Lower_probability_"//trim(adjustl(weightname))//",Upper_probability_"//trim(adjustl(weightname))
if(io/=0) stop "Output error for lower and upper probabilities."
do obs=1,size(plower)

	if(eapid)then
		bufferid=id(obs)
	else
		write(bufferid,'(i12)') obs
	end if
    write(buff,'(g25.16e3)') plower(obs)
    write(buff1,'(g25.16e3)') pupper(obs)
	write(unitpwt,'(5a)',iostat=io) trim(adjustl(bufferid)),",",&
        trim(adjustl(buff)),",",trim(adjustl(buff1))
    if(io/=0) stop "Output error for lower and upper probabilities."
end do

end subroutine outputpwt
