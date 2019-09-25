!
!  Transformation definition.
module transset
integer::components
integer,allocatable::degrees(:),elementsout(:)
logical,allocatable::lbmasks(:),ubmasks(:)
real(kind=8),allocatable::bases(:),coeff(:),inweights(:,:),lbs(:),ubs(:)
end module transset