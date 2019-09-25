!	Covariances of products of elements of latent vector with normal distribution with covariance matrix cov and mean mean.
function quarticnormal(cov,mean)
implicit none
real(kind=8),intent(in)::cov(:,:),mean(:)
real(kind=8)::quarticnormal(size(mean),size(mean),size(mean),size(mean))
!	Counters
integer::col,col2,row,row2
do row=1,size(mean)
	do col=1,row
		do row2=1,row
			col2=row2
			if(row==row2)col2=col
			quarticnormal(row,col,row2,1:col2)=cov(row,row2)*cov(col,1:col2)+cov(row,1:col2)*cov(row2,col)&
					+mean(row)*mean(row2)*cov(col,1:col2)&
					+mean(col)*mean(row2)*cov(row,1:col2)&
					+mean(row)*mean(1:col2)*cov(col,row2)&
					+mean(col)*mean(1:col2)*cov(row,row2)
			if(row/=row2)quarticnormal(row2,1:col2,row,col)=quarticnormal(row,col,row2,1:col2)
			
		end do
	end do
end do
return
end function quarticnormal

