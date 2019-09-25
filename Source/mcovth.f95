!	Multinomial covariance matrix.
!	density is vector of probabilities.
!	meanth is vector of means.
!	quadpoint is array of weights.
function mcovth(density,meanth,quadpoint)
implicit none
interface
!
	subroutine makesym(matrix)
!	Symmetrize n by n matrix based on lower triangle.
		real(kind=8),intent(inout)::matrix(:,:)
	end subroutine makesym
end interface
real(kind=8),intent(in)::density(:),meanth(:),quadpoint(:,:)
real(kind=8)::mcovth(size(meanth),size(meanth))
!	Count columns, quadrature points, and rows.
integer::col,quad,row
mcovth=0.0_8
if(size(density)>1)then 
    do row=1,size(meanth)
        do col=1,row
            do quad=1,size(density)
                mcovth(row,col)=mcovth(row,col)+(quadpoint(row,quad)-meanth(row))*&
                    (quadpoint(col,quad)-meanth(col))*density(quad)
            end do
        end do
    end do
    call makesym(mcovth)
end if

return
end function mcovth
