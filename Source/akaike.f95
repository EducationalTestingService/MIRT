!
!	Find the Akaike correction.
!	modeldim is the model dimension.
!	loglik is the log likelilhood.
!	totalitemsis the total number of presented items.
real(kind=8) function akaike(modeldim,loglik,totalitems)
implicit none
integer,intent(in)::modeldim
real(kind=8),intent(in)::loglik,totalitems
akaike=(-loglik+modeldim)/totalitems
return
end function akaike
