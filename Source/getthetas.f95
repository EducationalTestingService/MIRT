!	Find points for item response function computations.  Place in thetacheck.

subroutine getthetas(thetacheck)
implicit none
real(kind=8),intent(out)::thetacheck(:,:)
!	Error indicator for input.
integer::io
real(kind=8),allocatable::thetas(:,:)
namelist/irfpoints/thetas
allocate(thetas(size(thetacheck,1),size(thetacheck,2)),stat=io)
if(io/=0) stop 'Failure to allocate points for reading of item response functions.'
thetas=0.0_8
read(*,nml=irfpoints,iostat=io)
if(io/=0) stop 'Failure to read points for item response functions.'
thetacheck=thetas
return
end subroutine getthetas