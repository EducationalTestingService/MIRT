!	Find estimated conditional expectations and conditional covariance matrices for latent vectors.
!	postdensity is the densities.
!	theta is the latent vectors.
!	covtheta is the array of conditional covariance matrices.
!	meantheta is the array of conditional means. 
subroutine eap(postdensity,theta,covtheta,meantheta)
implicit none
interface

!	Multinomial covariance matrix.
	function mcovth(density,meanth,quadpoint)
		implicit none
		real(kind=8),intent(in)::density(:),meanth(:),quadpoint(:,:)
		real(kind=8)::mcovth(size(meanth),size(meanth))
	end function mcovth
	!	Multinomial mean.
	function mmeanth(density,quadpoint)
		implicit none
		real(kind=8),intent(in)::density(:),quadpoint(:,:)
		real(kind=8)::mmeanth(size(quadpoint,1))
	end function mmeanth
end interface
real(kind=8),intent(in)::postdensity(:,:),theta(:,:,:)
real(kind=8),intent(out)::covtheta(:,:,:),meantheta(:,:)
!	obs is observation number.
!	quad counts quadrature points.
integer::obs,quad
!	Obtain results for each observation.
do obs=1,size(postdensity,2)
	
! Observation obs.
	

!	mean and covariance.
	meantheta(:,obs)=mmeanth(postdensity(:,obs),theta(:,:,obs))
	covtheta(:,:,obs)=mcovth(postdensity(:,obs),meantheta(:,obs),theta(:,:,obs))

end do
return
end subroutine eap
