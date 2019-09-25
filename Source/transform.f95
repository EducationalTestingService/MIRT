!	Scale function
!	dimtrans is the dimension of the transformation.
!	theta is the argument.
function transform(dimtrans,theta)
use transset
implicit none
integer,intent(in)::dimtrans
real (kind=8),intent(in)::theta(:)
real(kind=8)::transform(dimtrans)
!	elementin is the source element.
!	elementout is the output element.
!	Component counter is i.
!	j counts degrees.
!	k and kk are used to find coefficients of polynomials.
integer::elementin,elementout,i,j,k,kk
!	bas is the base.
!	pol is the polynomial value.
!	th is the element of theta used in the component.
!	wt is the component weight.
real(kind=8)::bas,pol,th,wt
transform=0.0_8
kk=1

do i=1,components
	k=kk
	kk=kk+degrees(i)+1
	
	elementout=elementsout(i)
	th=dot_product(inweights(:,i),theta)
	wt=1.0_8
	!	See if th is within bounds.
	if(lbmasks(i)) then
		if(lbs(i)>th) cycle
	end if
	if(ubmasks(i))then
		if(ubs(i)<th)cycle
	end if
	if(th==lbs(i).or.th==ubs(i))wt=0.5_8
!	Polynomial computation.
	pol=0.0_8
	bas=bases(i)
	do j=0,degrees(i)
		pol=pol*(th-bas)+coeff(kk-j-1)
	end do
!	Add component.
	transform(elementout)=transform(elementout)+wt*pol
end do

end function transform

