!	Multinomial mean.
!	density is vector of probabilities
!	quadpoint is array of weights.
function mmeanth(density,quadpoint)
implicit none
real(kind=8),intent(in)::density(:),quadpoint(:,:)
real(kind=8)::mmeanth(size(quadpoint,1))
mmeanth=matmul(quadpoint,density)
return


end function mmeanth
