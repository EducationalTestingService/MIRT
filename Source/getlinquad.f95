!	Extract the linear component lintheta and the quadratic component quadtheta from the parameter vector beta.
subroutine getlinquad(beta,lintheta,quadtheta)
implicit none
interface
!	convert an n by n symmetric matrix from compact to regular form.
	function expandmat(n,compmat)
		implicit none
		integer,intent(in)::n
		real(kind=8),intent(in)::compmat(:)
		real(kind=8)::expandmat(n,n)
	end function expandmat

end interface

real(kind=8),intent(in)::beta(:)
real(kind=8),intent(out)::lintheta(:,:),quadtheta(:,:,:)
!	dimcovlat is the number of elements per predictor of the quadratic component.
!	dimlatin is the dimension of the underlying latent vector.
!	nlin points to linear components.
!	npred is the number of predictors.

!	nquadratic points to quadratic components.
!	position and position1 point to elements of the parameter vector.
!	pred counts predictors.
integer::dimcovlat,dimlatin,nlin,npred,nquadratic,position,position1,pred

dimlatin=size(lintheta,1)
npred=size(lintheta,2)
dimcovlat=dimlatin*(dimlatin+1)/2
nquadratic=size(beta)-dimcovlat*npred+1
nlin=nquadratic-dimlatin*npred
lintheta=reshape(beta(nlin:nquadratic-1),(/dimlatin,npred/))
position=nquadratic
do pred=1,npred
	position1=position+dimcovlat-1
	quadtheta(:,:,pred)=expandmat(dimlatin,beta(position:position1))
	position=position1+1
end do
return
end subroutine getlinquad