!	Standard error for entropy per item under standard sampling.
!	loglik is log likelihood.
!	obsweight is  observation weights.
!	prob is array of observed probabilities.
!	totalitems is weighted sum of items presented.
real(kind=8) function sdent(loglik,obsweight,prob,totalitems)
implicit none
real(kind=8),intent(in)::loglik,obsweight(:),prob(:),totalitems
!	Count observations.
integer::obs
!	Average item score.
real(kind=8)::loglikitem
loglikitem=loglik/sum(obsweight,prob>0.0_8)
sdent=0.0_8
do obs=1,size(prob)
	if(prob(obs)>0.0_8)sdent=sdent+obsweight(obs)*(log(prob(obs))-loglikitem)**2
end do
sdent=sqrt(sdent)/totalitems
return
end function sdent


