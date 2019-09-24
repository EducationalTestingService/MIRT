!
!  Find covariance matrix components for linear and quadratic parts of beta.
!	covbeta is the estimated covariance matrix of beta.
!	thpr is a collection of predictor vectors.
!   covlin is the estimated covariance matrix of linear part.
!	covlinquad is the estimated compact covariance matrix of linear and quadratic parts.
!	covquad is the estimated compact covariance matrix of quadratic parts.


subroutine covarlinquad(covbeta,thpr,covlin,covlinquad,covquad)
implicit none

real(kind=8),intent(in)::covbeta(:,:),thpr(:)
real(kind=8),intent(out)::covlin(:,:),covlinquad(:,:),covquad(:,:)
!   dimcovlat is the number of possibly distinct elements in a dimlatin by dimlatin covariance matrix.
!   dimlatin is the dimension of the underlying covariance matrix.
!   nlin points to linear components.
!   npred is the number of predictors.
!	nquadratic points to quadratic components.
!   position, position1, position2, and position3 point to elements of the parameter vector.
!   position1 points to elements of the parameter vector.
!   pred counts predictors.
!   pred1 counts predictors.




integer::dimcovlat,dimlatin,nlin,npred,nquadratic,position,position1,position2,position3,pred,pred1

!   Find sizes.
dimlatin=size(covlin,1)
dimcovlat=dimlatin*(dimlatin+1)/2
npred=size(thpr)
nquadratic=size(covbeta,1)-npred*dimcovlat+1
nlin=nquadratic-npred*dimlatin
!	Covariance matrix of linear part.
position=nlin

covlin=0.0_8
do pred=1,npred
	position1=position+dimlatin-1
	position2=nlin
	
	do pred1=1,npred
		position3=position2+dimlatin-1
		covlin=covlin+covbeta(position:position1,position2:position3)*thpr(pred)*thpr(pred1)
		position2=position3+1
	end do
	position=position1+1
	
end do

!	Covariance matrix of quadratic part.
position=nquadratic

covquad=0.0_8
do pred=1,npred
	position1=position+dimcovlat-1
	position2=nquadratic
	do pred1=1,npred
		position3=position2+dimcovlat-1
		covquad=covquad+covbeta(position:position1,position2:position3)*thpr(pred)*thpr(pred1)
		position2=position3+1
	end do
	position=position1+1
end do
!	Covariance matrix of linear and quadratic part.
covlinquad=0.0_8
position=nlin
do pred=1,npred
	position1=position+dimlatin-1
	position2=nquadratic
	do pred1=1,npred
		position3=position2+dimcovlat-1
		covlinquad=covlinquad+covbeta(position:position1,position2:position3)*thpr(pred)*thpr(pred1)
		position2=position3+1
	end do

	position=position1+1
end do
return

end subroutine covarlinquad
