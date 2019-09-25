!
!	Find estimated reliabilities of eap's of transformed latent vectors.
!	Input:
!	covartheta is the covariance matrix of the latent vector.
!	coveap is the covariance matrix of the eaps.
!	lintran is the transformation for the latent vector.
!	meancovtheta is the average conditional covariance of the latent vector given the observations.
!	meaneap is the average eap.

!	covarthetaskill is the covariance matrix of the transformed latent vector.
!	coveapskill is the covariance matrix of the eaps for the transformed latent vector.
!	meancovthetaskill is the average conditional covariance of the transformed latent vector
!		given the observations.
!	meaneapskill is the average eap for the transformed latent vector. 
!	releapskill is the reliability of each
!		element of the transformed latent vector.
subroutine relthetaskill(covartheta,coveap,lintran,meancovtheta,meaneap,&
	covarthetaskill,coveapskill,meancovthetaskill,meaneapskill,releapskill)
implicit none
real(kind=8),intent(in)::covartheta(:,:),coveap(:,:),&
	lintran(:,:),meancovtheta(:,:),meaneap(:)
	
real(kind=8),intent(out)::covarthetaskill(:,:),coveapskill(:,:),&
	meancovthetaskill(:,:),meaneapskill(:),releapskill(:)
!	row counter.
integer::row
meaneapskill=matmul(lintran,meaneap)
coveapskill=matmul(lintran,matmul(coveap,transpose(lintran)))
meancovthetaskill=matmul(lintran,matmul(meancovtheta,transpose(lintran)))
covarthetaskill=coveapskill+meancovthetaskill
releapskill=0.0_8
do row=1,size(meaneapskill)
	if(covarthetaskill(row,row)>0.0_8)releapskill(row)=coveapskill(row,row)/covarthetaskill(row,row)
end do

end subroutine relthetaskill