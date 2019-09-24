!	read predictor list.
!	predlabels are labels for predictor vectors.
!	predlist are predictor vectors.
subroutine getpredlist(predlabels,predlist)
implicit none
character(len=*),intent(out)::predlabels(:)
real(kind=8),intent(out)::predlist(:,:)
!	io is flag for allocation or read errors.
!	row counts predictor vectors.
integer::io,row
!	buff is used to generate default labels for predictor vectors.
!	predlabels are predictor labels.
character(len=10)::buff
character(len=32),allocatable::valuelabels(:)
!	predvalues are predictor values.
real(kind=8),allocatable::predvalues(:,:)
namelist/predictors/valuelabels,predvalues
allocate(predvalues(size(predlist,1),size(predlist,2)),valuelabels(size(predlist,2)),stat=io)
if(io/=0)stop 'Storage for predictor list not allocated.'
predvalues=0.0_8
predvalues(1,:)=1.0_8
do row=1,size(predlabels)
	write(buff,'(i10)') row
	valuelabels(row)=trim('Vector'//adjustl(buff))	
end do
read(*,nml=predictors,iostat=io)
if(io/=0) stop 'Predictor list not successfully read.'
predlist=predvalues
predlist(1,:)=1.0_8
predlabels=valuelabels
return
end subroutine getpredlist
