!	Find estimated conditional expectations and conditional covariance matrices for transformations of latent vector.
!	dimtrans is the dimension of the transformation.
!	postdensity contains posterior weights.
!	theta contains posterior quadrature points.

!	covtran is the conditional covariance of the transformation of the latent vector.
!	meantran is the conditional mean of the transformation of the latent vector.


subroutine eaptrans(dimtrans,postdensity,theta,covtran,meantran)
implicit none
interface
!	Multinomial covariance matrix.
!	density is vector of probabilities.
!	meanth is vector of means.
!	quadpoint is array of weights.
	function mcovth(density,meanth,quadpoint)
		implicit none
		real(kind=8),intent(in)::density(:),meanth(:),quadpoint(:,:)
		real(kind=8)::mcovth(size(meanth),size(meanth))
	end function mcovth
!	Multinomial mean.
!	density is vector of probabilities
!	quadpoint is array of weights.
	function mmeanth(density,quadpoint)
		implicit none
		real(kind=8),intent(in)::density(:),quadpoint(:,:)
		real(kind=8)::mmeanth(size(quadpoint,1))
	end function mmeanth
!	Define the transformation
!	m is the dimension of the result.
!	theta is the input.
	function transform(m,theta)
		implicit none
		integer,intent(in)::m
		real (kind=8),intent(in)::theta(:)
		real(kind=8)::transform(m)

	end function transform
	
end interface
integer,intent(in)::dimtrans
real(kind=8),intent(in)::postdensity(:,:),theta(:,:,:)
real(kind=8),intent(out)::covtran(:,:,:),meantran(:,:)
!	dimlatout is the dimension of the transformed latent vector.
!	item is an item.
!	ncat is the total number of underlying categories.
!	ncatobs is the total number of observed categories.
!	nitems is the number of items.
!	nquad is the number of quadrature points.
!	nscale is the location in beta of the initial item discrimination.
!	nscale1 is the location in beta of the last item discrimination.
!	obs is observation number.
!	quad counts quadrature points.
!	resp is a single response.
!	row is a row counter.

integer::obs,quad
!	scales is for scale parameters.

real(kind=8)::meanquad(dimtrans,size(postdensity,1))
	
	







!	Obtain results for each observation.
do obs=1,size(postdensity,2)
	
	
	

		
!	Cycle through quadrature points.

	do quad=1,size(postdensity,1)
		
		
		
! Probability computation.


		meanquad(:,quad)=transform(dimtrans,theta(:,quad,obs))
		
	end do

!	mean and covariance.
	meantran(:,obs)=mmeanth(postdensity(:,obs),meanquad)
    covtran(:,:,obs)=mcovth(postdensity(:,obs),meantran(:,obs),meanquad)
    
	
end do
return
end subroutine eaptrans
