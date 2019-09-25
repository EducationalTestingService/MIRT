!	Find estimated reliabilities of eaps of latent vectors.
!	Input:

!	covtheta is the array of conditional covariance matrices.
!	meantheta is the array of conditional means.

!	Output:

!	obsweight provides observation weights.
!	coveap is the covariance matrix of the eaps.
!	covartheta is the covariance matrix of the latent vector.
!	meancovtheta is the average conditional covariance of theta given the observations.
!	meaneap is the average eap.
!	releap is the reliability for each element of the latent vector.

subroutine reltheta(covtheta,meantheta,obsweight,coveap,covartheta,meancovtheta,meaneap,releap)
implicit none
interface
!	Multinomial mean.
	function mmeanth(density,quadpoint)
		implicit none
		real(kind=8),intent(in)::density(:),quadpoint(:,:)
		real(kind=8)::mmeanth(size(quadpoint,1))
	end function mmeanth
!	Multinomial covariance matrix.
	function mcovth(density,meanth,quadpoint)
		implicit none
		real(kind=8),intent(in)::density(:),meanth(:),quadpoint(:,:)
		real(kind=8)::mcovth(size(meanth),size(meanth))
	end function mcovth




end interface
real(kind=8),intent(in)::covtheta(:,:,:),meantheta(:,:),obsweight(:)
real(kind=8),intent(out)::coveap(:,:),covartheta(:,:),meancovtheta(:,:),meaneap(:),releap(:)
!	row is row number.
integer::row
!	Average conditional variance, mean eap, variance of eaps, total weight
real(kind=8)::condvar,mean,var,total
total=sum(obsweight)
releap=0.0_8
meaneap=mmeanth(obsweight,meantheta)/total
coveap=mcovth(obsweight,meaneap,meantheta)/total
do row=1,size(meaneap)
	meancovtheta(:,row)=mmeanth(obsweight,covtheta(:,row,:))/total
	covartheta(:,row)=meancovtheta(:,row)+coveap(:,row)
	if(covartheta(row,row)>0.0_8)releap(row)=coveap(row,row)/covartheta(row,row)
end do
	


return
end subroutine reltheta
