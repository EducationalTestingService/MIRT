!
!	Obtain predictor names.  Put in predname.
subroutine getprednames(predname)
implicit none
character(len=*),intent(out)::predname(:)

!	Name buffer
!	Replacement for predname
character(len=32),allocatable::pred_name(:)
character(len=4)::buffer
!	io is input flag.
!	pred is predictor number.
integer::io,pred
!	pred_lin_map and pred_quad_map give mappings for predictors and factors.
!	See if custom names needed.

namelist/predictorname/pred_name
allocate(pred_name(size(predname)),stat=io)
if(io/=0)stop "Predictor names not allocated successfully."



!	Default names.
pred_name(1)='Constant'

do pred=2,size(predname)
	write(buffer,'(i4)') pred-1
	pred_name(pred)=trim("Predictor"//adjustl(buffer))
	
end do

!	Get predictors.  Check for successful read.



read(*,nml=predictorname,iostat=io)

if(io/=0)stop "Predictor names not read successfully."
predname=pred_name

return

end subroutine getprednames
