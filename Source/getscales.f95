!	Obtain scale scores for weighted sum.
!	maxscore and minscore give maximum and minimum weighted sums.
!	numberscales is the number of scales.
!	Names are in scalename and scales are in scale.
subroutine getscales(maxscore,minscore,numberscales,scalename,scale)
implicit none
character(len=*),intent(out)::scalename(:)
integer,intent(in)::maxscore,minscore,numberscales
real(kind=8),intent(out)::scale(numberscales,minscore:maxscore)
!	Use a buffer for counts.
character(len=4)::buffer
!	Copy of scalename.
character(len=32),allocatable::scalenames(:)
!	i is an index.
!	io is flag for input error.
integer::i,io
!	scalevalues are values of scale.
real(kind=8),allocatable::scalevalues(:,:)
namelist/scales/scalenames,scalevalues
allocate(scalenames(numberscales),scalevalues(numberscales,minscore:maxscore),stat=io)
if(io/=0)stop 'Cannot allocate space for scales.'
do i=minscore,maxscore
	scalevalues(:,i)=i
end do
do i=1,numberscales
	write(buffer,'(i4)') i
	scalenames(i)='Scale'//trim(adjustl(buffer))
end do
read(*,nml=scales,iostat=io)
if(io/=0) stop 'Unable to read scales.'
scalename=scalenames
scale=scalevalues
return
end subroutine getscales
