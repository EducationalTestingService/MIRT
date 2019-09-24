!	Asymptotic covariance matrix of matrix inverse of symmetric matrix.
!	Input is matrix mat, its inverse matinv and its covariance matrix covmat.
!	Output is the asymptotic covariance matrix of matinv.
function covmatinv(mat,matinv,covmat)
implicit none
real(kind=8),intent(in)::mat(:,:),matinv(:,:),covmat(:,:,:,:)
real(kind=8)::covmatinv(size(mat,1),size(mat,1),size(mat,1),size(mat,1))
!   Indices i,j,k,l,m,n,o,p and array size s.
integer::i,j,k,l,m,n,o,p,s
s=size(mat,1)
covmatinv=0.0_8
do i=1,s
    do j=1,s
        do k=1,s
            do l=1,s
                do m=1,s
                    do n=1,s
						do o=1,s
							do p=1,s
								covmatinv(i,j,k,l)=covmatinv(i,j,k,l)+matinv(i,m)*matinv(k,n)*matinv(j,o)*matinv(l,p)*covmat(m,o,n,p)
							end do
						end do
                    end do
                end do
                
            end do
        
        end do
		
    end do

end do
return
end function covmatinv
