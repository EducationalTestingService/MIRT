!	Set quadrature points.

!	coord specifies the list of products to use.
!	dimpoints is the number of quadrature points per dimension.
!	cross indicates a cross-polytope.
!	even is for even spacing.
!	gausshermite is the indicator for Gauss-Hermite quadrature.
!	grid indicates a grid of points.
!	simplex indicates a simplex.
!	pspread is for point spread for even.
!	quadpoint provides quadrature points.
!	quadweight provides quadrature weights.
subroutine getquad(coord,dimpoints,cross,equalpoints,&
    even,gausshermite,grid,simplex,cweight,pspread,&
	quadpoint,quadweight)
implicit none
interface
!	Gauss-Hermite weight i of n.
	real(kind=8) function gaussweight(i,n)
		implicit none
		integer,intent(in)::i,n
	end function gaussweight
!	Gauss-Hermite point i of n.
	real(kind=8) function gausspoint(i,n)
		implicit none
		integer,intent(in)::i,n
	end function gausspoint
!	Normal weights for even spacing with centering at 0.
!	pspread is specified as the range of the quadrature points.
!		It is computed automatically  if the input is not positive.
!	Quadrature points are returned in points.
!	Quadrature weights are returned in weights.
	subroutine normalweight(pspread,points,weights)
		implicit none
		real(kind=8),intent(in)::pspread
		real(kind=8),intent(out)::points(:),weights(:)
	end subroutine normalweight
    subroutine readquaddim(points,weights)
        implicit none
        real(kind=8),intent(inout)::points(:),weights(:)
    end subroutine readquaddim
    subroutine readquadrature(vecpoints,weights)
        implicit none
real(kind=8),intent(inout)::vecpoints(:,:),weights(:)
    end subroutine readquadrature

end interface
integer,intent(in)::coord(:,:),dimpoints(:)
logical,intent(in)::cross,equalpoints,even,gausshermite,grid,simplex
real(kind=8),intent(in)::cweight(:),pspread
real(kind=8),intent(out)::quadpoint(:,:),quadweight(:)

!	dimno counts dimensions.
!	evenend ends a specification for grid components.
!	evenstart begins a specification for grid components.
!	io is an indicator for reading error.

!	quad is used for counting quadrature points.


integer::dimno,evenend,evenstart,io,quad
!	evenpoints has location components for products.
!	evenweights has weight components for products.
!	scale, scale 1, and scale 2 are scale factors.
real(kind=8):: evenpoints(sum(dimpoints)),&
	evenweights(sum(dimpoints)),scale,scale1,scale2
!	points for quadrature.
!	vecpoints for quadrature.
!	weights for quadrature.
real(kind=8),allocatable::points(:),vecpoints(:,:),weights(:)

!	Gauss-Hermite
if(grid) then
	if(gausshermite)then
	
		do quad=1,size(quadweight)
			scale=1.0_8
			do dimno=1,size(dimpoints)
				quadpoint(dimno,quad)=gausspoint(coord(dimno,quad),dimpoints(dimno))
				scale=scale*gaussweight(coord(dimno,quad),dimpoints(dimno))
			end do
			quadweight(quad)=scale

		end do
	else
!	Even.
		if(even)then
			evenstart=1
			
			do dimno=1,size(dimpoints)
				evenend=evenstart+dimpoints(dimno)-1
				call normalweight(pspread,&
					evenpoints(evenstart:evenend),evenweights(evenstart:evenend))
				evenstart=evenend+1
			end do
		else
!	User-supplied grid points.
			evenstart=1
			do dimno=1,size(dimpoints)
				allocate(points(dimpoints(dimno)),weights(dimpoints(dimno)),&
					stat=io)
				if(io/=0)stop "Quadrature points and weights not allocated."
                weights=1.0_8
                scale2=size(weights)
                scale1=(scale2+1.0_8)/2.0_8
                do quad=1,size(weights)
                    points(quad)=quad-scale1
                end do
                scale=sqrt(6.0_8/(scale1*(scale2-1.0_8)))
                points=scale*points
                call readquaddim(points,weights)
				if(minval(weights)<=0.0_8) &
					stop "Error in reading quadrature weights."
				evenend=evenstart+dimpoints(dimno)-1
				evenpoints(evenstart:evenend)=points
				evenweights(evenstart:evenend)=weights
				evenstart=evenend+1
				deallocate(points,weights)
			end do	
			
		end if
		do quad=1,size(quadweight)
			scale=1.0_8
			evenstart=0
			do dimno=1,size(quadpoint,1)
				evenend=evenstart+dimpoints(dimno)
				evenstart=evenstart+coord(dimno,quad)
				quadpoint(dimno,quad)=evenpoints(evenstart)
				scale=scale*evenweights(evenstart)
				evenstart=evenend
			end do
			quadweight(quad)=scale

		end do
	end if
	quadweight=quadweight*cweight
else
!	Simplex
	if(simplex)then
		quadpoint=0.0_8
        quadweight=1.0_8
        scale=size(quadweight)
        do quad=1,size(quadweight)
            scale1=quad+1
            do dimno=quad,size(dimpoints)
                scale1=dimno
                quadpoint(dimno,quad)=sqrt(scale/(scale1*(scale1+1.0_8)))
            end do
            scale1=quad-1
            scale2=quad
            if(quad>1)quadpoint(quad-1,quad)=-sqrt(scale*scale1/scale2)
		
		end do
        
!	Cross-polytope
	else
		if(cross)then
			quadweight=1.0_8
			quadpoint=0.0_8
			scale=size(dimpoints)
			scale=sqrt(scale)
			quad=1
			do dimno=1,size(dimpoints)
				quadpoint(dimno,quad)=-scale
				quad=quad+1
				quadpoint(dimno,quad)=scale
				quad=quad+1
			end do

		else
!	Custom
			allocate(vecpoints(size(dimpoints),size(quadweight)),&
				weights(size(quadweight)),stat=io)
			if(io/=0) stop "Quadrature points and weights not allocated."
			weights=1.0_8
			vecpoints=0.0_8
            call readquadrature(vecpoints,weights)
			quadpoint=vecpoints
			quadweight=weights
			if(minval(quadweight)<=0.0_8) stop "Error in reading quadrature weights."			
			
			
			
			
		end if
	end if
end if
!	Renormalize if needed.

quadweight=quadweight/sum(quadweight)
!	Reweight for standard normal case.
do quad=1,size(quadweight)
	quadweight(quad)=quadweight(quad)*&
		exp(dot_product(quadpoint(:,quad),quadpoint(:,quad))/2.0_8)
end do


return
end subroutine getquad
subroutine readquaddim(points,weights)
    implicit none
    real(kind=8),intent(inout)::points(:),weights(:)
    integer::io
    namelist/quaddim/points,weights
    read(*,nml=quaddim,iostat=io)
    if(io/=0)stop "Error in reading quadrature grid."

    return
end subroutine readquaddim
subroutine readquadrature(vecpoints,weights)
    implicit none
    real(kind=8),intent(inout)::vecpoints(:,:),weights(:)
    integer::io
    namelist/quadrature/vecpoints,weights
    read(*,nml=quadrature,iostat=io)
    if(io/=0)stop "Quadrature points and weights not read successfully."

    return
end subroutine readquadrature
