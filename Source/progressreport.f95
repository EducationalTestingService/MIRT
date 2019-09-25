!	Iteration progress reporting.
!	printprogstart indicates that starting iteration progress is to be printed to a file.
!	printprogstartstd indicates that starting iteration progress is to go to standard output.
!	printprog indicates that regular iteration progress is to be printed to a file.
!	printprogstd indicates that regular iteration progress is to go to standard output.
subroutine progressreport(printprogstart,printprogstartstd,printprog,printprogstd)

implicit none
logical,intent(out)::printprogstart,printprogstartstd,printprog,printprogstd
!	Error indicator.
integer::io
namelist/printprogress/printprogstart,printprogstartstd,printprog,printprogstd
printprogstart=.true.
printprogstartstd=.true.
printprog=.true.
printprogstd=.true.
read(*,nml=printprogress,iostat=io)
if(io/=0)stop "Progress reporting not successfully specified."
return
end subroutine progressreport
