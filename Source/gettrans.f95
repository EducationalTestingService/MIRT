!	Obtain transformation.
subroutine gettrans(dimlatin,transname)
use transset
implicit none
!	Names of elements of transformation.
character(len=*),intent(out)::transname(:)
character(len=32),allocatable::transnames(:)
!	Dimension of input.
integer,intent(in)::dimlatin
!	Buffer.
character(len=4)::buff
!	degree of polynomial.

!	elementout is the output element.
!	i counts components and output elements.
!	io is an error flag.
!	j counts coefficients.
integer::degree,elementout,i,io,j
!	Existence of lower or upper bounds of polynomial range.
logical::lbmask,ubmask
!	Reference base of polynomial.
!	lb is the lower bound.
!	ub is the upper bound
real(kind=8)::base,lb,ub
!	Polynomial coefficients
real(kind=8),allocatable::coefficients(:),inweight(:)
namelist/transformation/transnames,components
namelist/transcomp/degree,elementout,inweight,lbmask,ubmask,base,lb,ub
namelist/polynomial/coefficients
!	Space needed for names.
allocate(transnames(size(transname)),inweight(dimlatin),stat=io)
if(io/=0) stop 'Allocation of transformation names failed.'
!	Default names.
do i=1,size(transname)
	write(buff,'(i4)') i
	transnames(i)='Transformation'//trim(adjustl(buff))
end do
!	Component count.
components=size(transname)

read(*,nml=transformation,iostat=io)
if(io/=0)stop 'Transformation names and number of components not successfully read.'
transname=transnames
!   Run over components.
if(components<1) components=size(transname)
allocate(degrees(components),inweights(dimlatin,components),elementsout(components),&
    lbmasks(components),ubmasks(components),bases(components),lbs(components),ubs(components),&
    stat=io)
if(io/=0) stop 'Allocation of transformation data failed.'
!	Set defaults and fix input.
do i=1,components
    degree=1
    elementout=i
    inweight=0.0_8
	inweight(min(i,dimlatin))=1.0_8
    lbmask=.FALSE.
    ubmask=.FALSE.
    lb=0.0_8
    ub=0.0_8
    base=0.0_8
    read(*,nml=transcomp,iostat=io)
    if(io/=0) stop 'Component data not successfully read for transformation.'
	degree=max(0,degree)
    degrees(i)=degree
	inweights(:,i)=inweight
	elementout=max(1,min(elementout,size(transnames)))
    elementsout(i)=elementout
    lbmasks(i)=lbmask
    ubmasks(i)=ubmask
    lbs(i)=lb
    ubs(i)=ub
    bases(i)=base
end do
allocate(coeff(sum(degrees)+components),stat=io)
if(io/=0) stop 'Allocation of transformation definition not successful.'
j=1
do i=1,components
    allocate(coefficients(0:degrees(i)),stat=io)
    coefficients(0:degrees(i))=0.0_8
    if(degrees(i)>0)coefficients(1)=1.0_8
    read(*,nml=polynomial,iostat=io)
    if(io/=0) stop 'Coefficients of polynomial not read correctly.'
    coeff(j:j+degrees(i))=coefficients
    deallocate(coefficients)
    j=j+degrees(i)+1
end do
return
end subroutine gettrans	