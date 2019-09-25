
!	Total observations presented.
!	dat contains response data.  obsweight are observation weights.
!	numcatobs contains observed category counts by item.
real(kind=8) function obspresented(dat,numcatobs,obsweight)
implicit none
integer,intent(in)::dat(:,:),numcatobs(:)
real(kind=8),intent(in)::obsweight(:)
!	obs counts observations.
integer::obs
obspresented=0.0_8
do obs=1,size(obsweight)
	if(any(0<= dat(:,obs).and.dat(:,obs)<numcatobs))obspresented=obspresented+obsweight(obs)
end do

end function obspresented