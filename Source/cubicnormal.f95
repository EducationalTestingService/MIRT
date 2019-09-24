!	Covariances of products of elements of latent vector and of elements of latent vector with normal
!	distribution with mean mean and  covariance cov.
function cubicnormal(cov,mean)
implicit none
real(kind=8),intent(in)::cov(:,:),mean(:)
real(kind=8)::cubicnormal(size(mean),size(mean),size(mean))
!	Counters
integer::col,row,row2
do row=1,size(mean)
	do col=1,row
		do row2=1,size(mean)
			cubicnormal(row,col,row2)=mean(col)*cov(row,row2)+mean(row)*cov(col,row2)
			
		end do
	end do
end do
return
end function cubicnormal

