!	Multinomial fourth moments.  Covariances of cross-products.
!	density provides weights.
!	points are in quadpoint.
!	Product means are in mean2th.
function mquarth(density,mean2th,quadpoint)
implicit none
real(kind=8),intent(in)::density(:),mean2th(:,:),quadpoint(:,:)
real(kind=8)::mquarth(size(mean2th,1),size(mean2th,1),size(mean2th,1),size(mean2th,1))
!	Counters
integer::col,col1,row,row1,quad,topcol
!	Product
mquarth=0.0_8
do row=1,size(mean2th,1)
	do col=1,row
		do row1=1,row
			topcol=row1
			if(row==row1)topcol=col
			do col1=1,topcol
				do quad=1,size(density)
					mquarth(row,col,row1,col1)=mquarth(row,col,row1,col1)+&
						(quadpoint(row,quad)*quadpoint(col,quad)-mean2th(row,col))*&
						(quadpoint(row1,quad)*quadpoint(col1,quad)-mean2th(row1,col1))*&
						density(quad)
				end do
				
				if(row/=row1.or.col/=col1)mquarth(row1,col1,row,col)=mquarth(row,col,row1,col1)
			end do
		end do
	end do
end do
return
end function mquarth
