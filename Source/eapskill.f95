!	Find estimated conditional expectations and conditional covariance matrices for transformed latent vectors.
!	covtheta is the array of conditional covariance matrices.
!	lintran is the linear transformation for the latent vector.
!	meantheta is the array of conditional means.
!	covnewtheta is the array of transformed conditional covariances.
!	meannewtheta is the  array of transformed conditional means. 
subroutine eapskill(covtheta,lintran,meantheta,covnewtheta,meannewtheta)
implicit none
real(kind=8),intent(in)::covtheta(:,:,:),lintran(:,:),meantheta(:,:)
real(kind=8),intent(out)::covnewtheta(:,:,:),meannewtheta(:,:)
!	obs is observation number.
integer::obs
!	Set up parameter arrays for simplified processing.
do obs=1,size(meantheta,2)
	covnewtheta(:,:,obs)=matmul(lintran,matmul(covtheta(:,:,obs),transpose(lintran)))
	meannewtheta(:,obs)=matmul(lintran,meantheta(:,obs))
	
end do

return
end subroutine eapskill
