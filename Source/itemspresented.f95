!	Total items presented.
!	dat contains response data.  obsweight are observation weights.
!	numcatobs contains observed category counts by item.
real(kind=8) function itemspresented(dat,numcatobs,obsweight)
implicit none
integer,intent(in)::dat(:,:),numcatobs(:)
real(kind=8),intent(in)::obsweight(:)
!	obs counts observations.
integer::obs
itemspresented=sum((/(obsweight(obs)*count(0<= dat(:,obs).and.dat(:,obs)<numcatobs),obs=1,size(dat,2))/))

end function itemspresented
