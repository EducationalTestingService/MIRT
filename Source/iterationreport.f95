!	Report on iteration progress.
!	comment is used to comment on the iteration progress report.
!	it is the iteration number.
!	unitno is the unit number for the report.
!	printprog indicates that iteration progress is to be printed to unitno.
!	printprogstd indicates that iteration progress is to go to standard output.
!	loglik is the log likelihood.
!	step is the step size. 
subroutine iterationreport(comment,it,unitno,printprog,printprogstd,loglik,step)
use,intrinsic::iso_fortran_env
implicit none
character(len=*),intent(in)::comment
integer,intent(in)::it,unitno
logical,intent(in)::printprog,printprogstd
real(kind=8),intent(in)::loglik,step
!	buff is used to trim it+1.
!	buff1 is used to trim step.
!	buff2 is used to trim loglik.
character(len=12)::buff
character(len=25)::buff1,buff2
write(buff,'(i12)')it+1
write(buff1,'(g25.16e3)')step
write(buff2,'(g25.16e3)')loglik
if(printprog)then
	if(it==0)write(unitno,'(2a)') "Summary of iteration progress ",comment
	if(it==0)write(unitno,'(a)') "Iteration,Step,Log-likelihood"
	write(unitno,'(5A)') trim(adjustl(buff)),",",trim(adjustl(buff1)),",",trim(adjustl(buff2))
    flush(unitno)
end  if
if(printprogstd)then
	if(it==0)write(output_unit,'(2a)') "Summary of iteration progress ",comment
	if(it==0)write(output_unit,'(a)') "Iteration,Step,Log-likelihood"
	write(output_unit,'(5A)') trim(adjustl(buff)),",",trim(adjustl(buff1)),",",trim(adjustl(buff2))
    flush(output_unit)
end if
return
end subroutine iterationreport
