!	Print reliability results for observed scores.
!	comment is used to identify output.
!	factorname is variable names.
!	unitobsrel is unit.
!	covobs is the covariance matrix of the observations vector.
!	covres is the covariance matrix of the residuals.
!   obsave is the average of the scores.
!	rel contains reliabilities.

subroutine obsreloutput(comment,factorname,unitobsrel,&
    covobs,covres,obsave,rel)
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
		character(len=*),intent(in)::rowname(:)
		integer,intent(in)::unitmat
		real(kind=8),intent(in)::mat(:,:)
	end subroutine matoutput
end interface
character(len=*),intent(in)::comment
character(len=32),intent(in)::factorname(:)
integer,intent(in)::unitobsrel
real(kind=8),intent(in)::covobs(:,:),covres(:,:),obsave(:),rel(:)

!	Count rows.
integer::row
write(unitobsrel,'(2a)') "Reliability analysis ",comment
write(unitobsrel,'(a)') "Averages of scores"
do row=1,size(rel)
	call writelabelfl(factorname(row),unitobsrel,obsave(row))
	
end do

call matoutput("Covariance matrix of scores",factorname,unitobsrel,covobs)
call matoutput("Covariance matrix of residuals",factorname,unitobsrel,covres)


write(unitobsrel,'(a)') "Reliabilities of observed scores"
do row=1,size(rel)
	call writelabelfl(factorname(row),unitobsrel,rel(row))
end do

return



end subroutine obsreloutput


