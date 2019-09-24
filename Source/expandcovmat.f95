!   Expand covariance matrix of matrix from compact format.
function expandcovmat(n,covmat)
implicit none
integer,intent(in)::n
real(kind=8),intent(in)::covmat(:,:)
real(kind=8)::expandcovmat(n,n,n,n)
!   h, i, j, k, l, m, and n are indices

integer::h,i,j,k,l,m
l=1
do h=1,n
    do i=1,h
        m=1
        do j=1,h
            do k=1,j
                
                if(k<j.and.h==i)then
					expandcovmat(h,i,j,k)=2.0_8*covmat(l,m)
					expandcovmat(h,i,k,j)=expandcovmat(h,i,j,k)
				endif
				if(i<h.and.k==j)then
					expandcovmat(h,i,j,k)=2.0_8*covmat(l,m)
					expandcovmat(i,h,j,k)=expandcovmat(h,i,j,k)
				end if
                if(i<h.and.k<j)then
					expandcovmat(h,i,j,k)=4.0_8*covmat(l,m)
					expandcovmat(h,i,k,j)=expandcovmat(h,i,j,k)
					expandcovmat(i,h,j,k)=expandcovmat(h,i,j,k)
					expandcovmat(i,h,k,j)=expandcovmat(h,i,j,k)
				end if
				if(i==h.and.j==k)expandcovmat(h,i,j,k)=covmat(l,m)
                m=m+1
            end do
            
        end do
        l=l+1
    end do
end do

return
end function expandcovmat