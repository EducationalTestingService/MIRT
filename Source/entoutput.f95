
!	Information output
!   comment is a comment.
!	unitinfo is the unit number.
!	modeldim is the rank of the design matrix.
!	complx indicates complex sampling.
!	akent is corresponding Akaike estimate.
!	ent is entropy per presented item.
!	ghent is coresponding Gilula-Habermann estimate.
!	ghent_complex is corresponding Gilula-Haberman estimate for complex sampling.
!	loglik is log likelihood.
!	sdentropy is estimated asymptotic standard deviation of ent.
!	sdentropy_complex is asymptotic standard error of ent for complex sampling.
subroutine entoutput(comment,modeldim,unitinfo,complx,akent,ent,ghent,ghent_complex,loglik,sdentropy,sdentropy_complex)
implicit none
interface
!	write a label and a floating point number x in csv format to unit.
	subroutine writelabelfl(label,unit,x)
		implicit none
		character(len=*),intent(in)::label
		integer,intent(in)::unit
		real(kind=8),intent(in)::x
	end subroutine writelabelfl
!	write a label and an integer nvin csv format to unit.
	subroutine writelabelint(label,nv,unit)
		implicit none
		character(len=*),intent(in)::label
		integer,intent(in)::nv,unit
	end subroutine writelabelint
end interface
character(len=*),intent(in)::comment
integer,intent(in),optional::modeldim
integer,intent(in)::unitinfo
logical,intent(in)::complx
real(kind=8),intent(in),optional::akent,ent,ghent,ghent_complex,loglik,sdentropy,sdentropy_complex

write(unitinfo,'(a)') comment
if(present(loglik))call writelabelfl("Log_likelihood",unitinfo,loglik)
if(present(modeldim))call writelabelint("Mod_dim",modeldim,unitinfo)
if(present(ent))call writelabelfl("Penalty",unitinfo,ent)
if(present(sdentropy))call writelabelfl("SE_Penalty",unitinfo,sdentropy)
if(present(akent))call writelabelfl("Akaike",unitinfo,akent)
if(present(ghent))call writelabelfl("Gilula-Haberman",unitinfo,ghent)
if(complx.and.present(ghent_complex))then
	call writelabelfl("Complex_Gilula_Haberman",unitinfo,ghent_complex)
	call writelabelfl("Complex_SE",unitinfo,sdentropy_complex)
end if
return



end subroutine entoutput


