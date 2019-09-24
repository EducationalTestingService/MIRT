!
!	Find dimension of the non-negative definite matrix mat.
!	tolsing is the singularity tolerance.

integer function dimfind(mat,tolsing)
implicit none
real(kind=8),intent(in)::mat(:,:),tolsing
interface
!	This function computes the modified Cholesky decomposition
!	of the n by n nonnegative definite symmecholc macholx sym.
!	tolsing is the singularity tolerance.	
	function chol(sym,tolsing)
		implicit none
		real(kind=8),intent(in)::sym(:,:),tolsing
	 	real(kind=8)::chol(size(sym,1),size(sym,1))
	end function chol
end interface
!	row counts rows.
integer::row
!	Modified Cholesky decomposition is in matchol.
real(kind=8)::matchol(size(mat,1),size(mat,2))
matchol=chol(mat,tolsing)
dimfind=count((/(matchol(row,row)>0.0_8,row=1,size(mat,1))/))
return
end function dimfind