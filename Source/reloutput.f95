!	Print reliability results for eap's.
!	comment is used to identify output.
!	factorname is variable names.
!	unitrel is unit.
!	covartheta is the covariance matrix of the latent vector.
!	coveap is the covariance matrix of the eaps.
!	meancovtheta is the average conditional covariance of the
!		latent vector given the observations.
!	meaneap is the average eap.
!	releap contains reliabilities.

subroutine reloutput(comment,factorname,unitrel,covartheta,coveap,meancovtheta,meaneap,releap)
implicit none
interface

!	write a label and a floating point number x in csv format to unit.
	subroutine writelabelfl(label,unit,x)
		implicit none
		character(len=*),intent(in)::label
		integer,intent(in)::unit
		real(kind=8),intent(in)::x
	end subroutine writelabelfl
!	print square matrix.
!	method is the type of matrix.
!	rowname provides row names.
!	unitmat is the unit number to use.
!	mat provides the matrix.
	subroutine matoutput(comment,rowname,unitmat,mat)
		implicit none
		character(len=*),intent(in)::comment
		character(len=64),intent(in)::rowname(:)
		integer,intent(in)::unitmat
		real(kind=8),intent(in)::mat(:,:)
	end subroutine matoutput
end interface
character(len=*),intent(in)::comment
character(len=*),intent(in)::factorname(:)
integer,intent(in)::unitrel
real(kind=8),intent(in)::covartheta(:,:),coveap(:,:),meancovtheta(:,:),meaneap(:),releap(:)

!	Count rows.
integer::row
write(unitrel,'(2a)') "Reliability analysis ",comment
write(unitrel,'(a)') "Averages of conditional means"
do row=1,size(releap)
	call writelabelfl(factorname(row),unitrel,meaneap(row))
	
end do

call matoutput("Covariance matrix of conditional means",factorname,unitrel,coveap)
call matoutput("Conditional covariance matrix",factorname,unitrel,meancovtheta) 
call matoutput("Unconditional covariance matrix",factorname,unitrel,covartheta)

write(unitrel,'(a)') "Reliabilities of conditional means"
do row=1,size(releap)
	call writelabelfl(factorname(row),unitrel,releap(row))
end do

return



end subroutine reloutput


