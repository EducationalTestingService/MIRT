!	print probabilities.
!	id is individual id. prob for model.
!	unitprob is output unit.
!	eapid indicates if ids are put out.
!	prob is array of probabilities.
subroutine outputprob(id,unitprob,eapid,prob)
implicit none
interface
!	write a label and a floating point number x in csv format to unit.
	subroutine writelabelfl(label,unit,x)
		implicit none
		character(len=*),intent(in)::label
		integer,intent(in)::unit
		real(kind=8),intent(in)::x
	end subroutine writelabelfl
end interface
character(len=*),intent(in),optional::id(:)
integer,intent(in)::unitprob
logical,intent(in)::eapid
real(kind=8),intent(in)::prob(:)
!	buffer.
character(len=32)::bufferid
!	Observation counter
integer::obs
write(unitprob,'(a)') "ID,Probability"
do obs=1,size(prob)
	if(eapid)then
		bufferid=id(obs)
	else
		write(bufferid,'(i12)') obs
	end if
	call writelabelfl(bufferid,unitprob,prob(obs))
end do

end subroutine outputprob
