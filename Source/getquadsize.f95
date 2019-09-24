
!	Find the size specifications for the quadrature points.
!	dimpoints gives the number  of points per diemnsion if a
!	grid is used.
!	If a grid is used, then the grid has ngrid points.
!	In all cases, the number of quadrature points is nquad.
!	cross-polytope quadrature is used if cross is .true.
!	Each grid dimension has the same number of points if
!	equalpoints is .true.
!	If even is .true. and a grid is used, then
!		evenly spaced quadrature points are used.
!	A complete grid is used if a grid is used and fullgrid is .true.
!	If guasshermite is .true., even is .false., and a grid is used,
!		then Gauss-Hermite quadrature points and weights are
!		produced.
!	A grid is used if grid is .true. and cross and simplex are .true.
!	normal indicates that the latent vector has a normal distribution.
!	Simplex quadrature is used if simplex is .true. and cross is .false.
!	The subroutine returns grid and fullgrid as .false. if  simplex  or
!		cross is .true.
!	The subroutine returns simplex as .false. if cross is .true.
!	In the default case, cross is .false., even is .false., fullgrid is .true.,
!		gausshermite is .true., grid is .true., and simplex is .false.
!	The default value of each element of dimpoints is 5 if Gauss-Hermite
!		quadrature is used and 6 if a grid is used and Gauss-Hermite
!		quadrature is not used.  If Gauss-Hermite quadrature is used, then no
!		element of dimpoints is permitted to exceed 30.
!		Any element of dimpoints is forced to be at least 2 if a grid is used.
!	If a grid is used but nquad is less than 2, then fullgrid is set to .true.
!	If nquad exceeds ngrid and a grid is used, then nquad is set to ngrid.
!	If a full grid is used, then nquad is set to ngrid.
!	If cross is .true., then nquad is set to twice size(dimpoints).
!	If simplex is returned as .true., then nquad is size(dimpoints)+1.
!	If simplex, cross, and grid are .false., then nquad must be given and
!		at least 2.  Otherwise, processing ends with an error message.


subroutine getquadsize(dimpoints,ngrid,nquad,&
	cross,equalpoints,even,fullgrid,gausshermite,grid,normal,simplex,&
	pspread)

implicit none

interface
!	Obtain number of quadrature points for each dimension.
!	The result is placed in dimpoints.
	subroutine getdimpoints(dimpoints)
		implicit none
		integer,intent(out)::dimpoints(:)
	end subroutine getdimpoints
end interface


integer,intent(out)::dimpoints(:),ngrid,nquad
logical,intent(out)::cross,equalpoints,even,fullgrid,gausshermite,grid,normal,simplex
real(kind=8),intent(out)::pspread
!	dimsize is the number of points per dimension.
!	io is the flag for input errors.
integer,parameter::maxquad=30
integer::dimsize,io

namelist/quadsize/dimsize,nquad,cross,equalpoints,even,fullgrid,gausshermite,grid,normal,simplex,&
		pspread
!	Default values.
dimsize=5
nquad=0
cross=.false.
equalpoints=.true.
even=.false.
fullgrid=.true.
gausshermite=.true.
grid=.true.
normal=.true.
simplex=.false.

pspread=0.0_8


read(*,nml=quadsize,iostat=io)
if(io/=0) stop "Quadrature points not successfully specified."
!	Resolve discrepancies.
if(even)gausshermite=.false.
if(cross.or.simplex)grid=.false.
if(cross)simplex=.false.
if(.not.grid)then
	even=.false.
	fullgrid=.false.
	gausshermite=.false.
!	Count points.
	if(cross) nquad=2*size(dimpoints)
	if(simplex) nquad=size(dimpoints)+1
	
else
	if(equalpoints)then
		if(dimsize<1)dimsize=5
		if(gausshermite)dimsize=min(dimsize,maxquad)
		dimpoints=dimsize
	else
		call getdimpoints(dimpoints)
		if(gausshermite)dimpoints=min(dimpoints,maxquad)
	end if
	ngrid=product(dimpoints)
	if(fullgrid.or.nquad<2)then
		nquad=ngrid
	else
		nquad=min(nquad,ngrid)
		if(nquad==ngrid)fullgrid=.true.
	end if
end if



if(nquad<1) stop "Insufficient quadrature points"
return
end subroutine getquadsize
