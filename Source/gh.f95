!
!	Find the Gilula-Haberman correction.  
!	eacov is the estimated asymptotic covariance matrix under the model.
!	inform is the estimated information obtained by the Louis approach.
!	loglik is the log likelihood.
!	totalitems is the total weighted number of presented items.
!
!
real(kind=8) function gh(eacov,inform,loglik,totalitems)
implicit none
real(kind=8),intent(in)::eacov(:,:),inform(:,:),loglik,totalitems
!	col counts columns.
!	obs is an observation counter.
!	row counts rows.
integer::col,obs,row
gh=-loglik
do row=1,size(eacov,1)
	do col=1,row
		if(col<row)then
			gh=gh+2.0_8*eacov(row,col)*inform(col,row)
		else
			gh=gh+eacov(row,row)*inform(row,row)
		end if
	end do
end do
gh=gh/totalitems
return
end function gh
