!	Convert from adaptive quadrature scale to normalized latent vector.
!	Use nhchol for modified Cholesky decomposition of negative hessian.
!	Use oldalpha for location.
!	Use scalez for scaling.
!	Quadrature point is z.
function scaleadapt(nhchol,oldalpha,scalez,z)
implicit none

real(kind=8),intent(in)::nhchol(:,:),oldalpha(:),scalez(:),z(:)
real(kind=8)::scaleadapt(size(z))
!	dimlat is dimension of z.
!	row is row number
integer::dimlat,row
dimlat=size(z)

scaleadapt(dimlat)=z(dimlat)*scalez(dimlat)
if (dimlat>1) then
	do row=dimlat-1,1,-1
		scaleadapt(row)=z(row)*scalez(row)&
			-dot_product(nhchol(row+1:dimlat,row),scaleadapt(row+1:dimlat))
	end do
end if

scaleadapt=scaleadapt+oldalpha
return




end function scaleadapt
