!	Multinomial third moments.  Covariance of product and single term.
!	density contains weights.
!	meanth contains means.
!	mean2th contains means of products.
!	quadpoint contains points.

function mcubth(density,meanth,mean2th,quadpoint)
implicit none
real(kind=8),intent(in)::density(:),meanth(:),mean2th(:,:),quadpoint(:,:)
real(kind=8)::mcubth(size(meanth),size(meanth),size(meanth))
!	Counts blocks, columns, rows, and quadrature points.
integer::block,col,row,quad
mcubth=0.0_8
do row=1,size(meanth)
	do col=1,row
		do block=1,size(meanth)
			do quad=1,size(density)
				mcubth(row,col,block)=mcubth(row,col,block)+&
					(quadpoint(row,quad)*quadpoint(col,quad)-mean2th(row,col))*&
					(quadpoint(block,quad)-meanth(block))*&
					density(quad)
			end do
		end do
	end do
end do
return

end function mcubth
