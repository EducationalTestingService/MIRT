!
!  Find statistics of mean and covariance matrix of underlying latent vector.
!

!   covinverse is the estimated inverse of the covariance matrix of the underlying latent vector corresponding to beta.

!	covlin is the estimated covariance matrix of the linear component of beta.
!	covlinquad is the estimated covariance matrix of the linear and quadratic components of beta.
!   covmatrix is the estimated covariance matrix of the underlying latent vector corresponding to beta.
!	covquad is the estimated covariance matrix of the quadratic component of beta.
!	meanvec is the estimated mean of the underlying latent vector.

!   covcovmatrix is the estimated covariance matrix of covmatrix.
!   covmeanvec is the estimated covariance matrix of meanvec.
!   Note that this routine is for the case of an underlying normal latent vector.
subroutine covmatstat(covinverse,covlin, covlinquad,covmatrix,covquad,meanvec,covcovmatrix,covmeanvec)
implicit none
interface

!	Asymptotic covariance matrix of matrix inverse of symmetric matrix.
!	Input is matrix mat, its inverse matinv and its covariance matrix covmat.
!	Output is the asymptotic covariance matrix of matinv.
	function covmatinv(mat,matinv,covmat)
		implicit none
		real(kind=8),intent(in)::mat(:,:),matinv(:,:),covmat(:,:,:,:)
		real(kind=8)::covmatinv(size(mat,1),size(mat,1),size(mat,1),size(mat,1))
	end function covmatinv
!   Expand covariance matrix of matrix from compact format.
	function expandcovmat(n,covmat)
		implicit none
		integer,intent(in)::n
		real(kind=8),intent(in)::covmat(:,:)
		real(kind=8)::expandcovmat(n,n,n,n)
	end function expandcovmat
!	convert an n by n symmetric matrix from compact to regular form.
	function expandmat(n,compmat)
		implicit none
		integer,intent(in)::n
		real(kind=8),intent(in)::compmat(:)
		real(kind=8)::expandmat(n,n)
	end function expandmat

!		This program uses the modified Cholesky decomposition
!		in lu to solve the equation lusolve = b.

	subroutine makesym(matrix)
!	Symmetrize n by n matrix based on lower triangle.
		implicit none
		real(kind=8),intent(inout)::matrix(:,:)
	end subroutine makesym


end interface
real(kind=8),intent(in)::covinverse(:,:),covlin(:,:),covlinquad(:,:),covmatrix(:,:),covquad(:,:),meanvec(:)
real(kind=8),intent(out)::covcovmatrix(:,:,:,:),covmeanvec(:,:)
!   dimlatin is the dimension of the underlying covariance matrix.

!   row counts elements of latent vectors.
!   row1 counts elements of latent vectors.

integer::dimlatin,row,row1
!   covmeanvec1,covmeanvec2, and covmeanvec3 are components of covmeanvec.
!   covquade is the expanded covariance matrix of the quadratic component of beta.

real(kind=8)::covmeanvec1(size(covmatrix,1),size(covmatrix,1)),&
	covmeanvec2(size(covmatrix,1),size(covmatrix,1),size(covmatrix,1)),&
	covmeanvec3(size(covmatrix,1),size(covmatrix,1)),&
	covquade(size(covmatrix,1),size(covmatrix,1),size(covmatrix,1),size(covmatrix,1))
	
!   Find dimension of underlying latent vector.
dimlatin=size(covmatrix,1)



!	Covariance matrix of linear part.


covquade=expandcovmat(dimlatin,covquad)
covcovmatrix=covmatinv(covinverse,covmatrix,covquade)
do row=1,dimlatin
	
	covmeanvec1(row,:)=matmul(expandmat(dimlatin,covlinquad(row,:)),meanvec)
	do row1=1,dimlatin
		covmeanvec2(row,:,row1)=matmul(covquade(row,:,row1,:),meanvec)
		covmeanvec3(row,row1)=sum(covmeanvec2(row,:,row1)*meanvec)
	end do
end do

do row=1,dimlatin
	do row1=1,row
		covmeanvec(row,row1)=covlin(row,row1)-covmeanvec1(row,row1)-covmeanvec1(row1,row)+covmeanvec3(row,row1)
	end do
end do
call makesym(covmeanvec)
covmeanvec=matmul(covinverse,matmul(covmeanvec,covinverse))		
end subroutine covmatstat