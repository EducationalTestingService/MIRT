!	Some basic summary information.
!	filename is the name of the data file.	
!	nitems is the number of items.
!	nobs is the number of observations.
!	npred is the number of predictors.
!	unitno is the unit number used.
!	totalitems is weighted number of items presented.
!	totalobs is weighted number of observations presented.
	
subroutine outputdata(filename,nitems,nobs,npred,unitno,totalitems,totalobs)
implicit none
interface
!	write a label and an integer nv in csv format to unit.
	subroutine writelabelint(label,nv,unit)
		implicit none
		character(len=*),intent(in)::label
		integer,intent(in)::nv,unit
	end subroutine writelabelint
!	write a label and a floating point number x in csv format to unit.
	subroutine writelabelfl(label,unit,x)
		implicit none
		character(len=*),intent(in)::label
		integer,intent(in)::unit
		real(kind=8),intent(in)::x
	end subroutine writelabelfl
end interface
character(len=256)::filename
integer,intent(in)::nitems,nobs,npred,unitno
real(kind=8),intent(in)::totalitems,totalobs
!	io  flag
integer::io
call writelabelint("Num_obs",nobs,unitno)
call writelabelint("Num_items",nitems,unitno)
call writelabelint("Num_pred",npred,unitno)
write(unitno,'(2a)',iostat=io) "File_name,", trim(filename)
call writelabelfl("Tot_items_presented",unitno,totalitems)
call writelabelfl("Tot_obs_presented",unitno,totalobs)
call writelabelfl("Ave_items_presented",unitno,totalitems/totalobs)
if(io/=0)stop "Output failure"
return
end subroutine outputdata
