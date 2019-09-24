!	beta output
!	betaname is the vector of betar names.
!	unitbeta is the unit number.
!	complex indicates complex sampling.
!	beta is the vector of estimates.
!	eacovbeta is the array of estimated covariances of beta estimates under the model.
!	eacovbeta_louis is the array of estimated covariances of beta estimates under the model with the Louis approach.
!	eacovbeta_sandwich is the array of estimated covariances of beta estimates without assuming the model.
!	eacovbeta_complex is the array of estimated covariances of beta estimates under complex sampling.

subroutine betaoutput(betaname,unitbeta,complx,beta,eacovbeta,eacovbeta_complex,eacovbeta_louis,eacovbeta_sandwich)
implicit none
character(len=64)::betaname(:)
integer,intent(in)::unitbeta
logical,intent(in)::complx
real(kind=8),intent(in)::beta(:),eacovbeta(:,:),eacovbeta_complex(:,:),eacovbeta_louis(:,:),eacovbeta_sandwich(:,:)
!	buffers for output.
character(len=25)::buffe,buffse,buffsec,buffsel,buffses
!	Count rows.
integer::row

write(unitbeta,'(a)') "Beta estimates and standard errors"

if(complx)then
	write(unitbeta,'(6a)') "Beta,", "Estimate,", "Std. Err.,", "Louis S.E.,", "Sandwich S.E.,","Complex S.E."
	do row=1,size(beta)
		write(buffe,'(g25.16e3)') beta(row)
		write(buffse,'(g25.16e3)') sqrt(max(0.0_8,eacovbeta(row,row)))
		write(buffsec,'(g25.16e3)') sqrt(max(0.0_8,eacovbeta_complex(row,row)))
		write(buffsel,'(g25.16e3)') sqrt(max(0.0_8,eacovbeta_louis(row,row)))
		write(buffses,'(g25.16e3)') sqrt(max(0.0_8,eacovbeta_sandwich(row,row)))
		write(unitbeta,'(11a)') trim(adjustl(betaname(row))),",",&
			trim(adjustl(buffe)),",",&
			trim(adjustl(buffse)),",",&
			trim(adjustl(buffsel)),",",&
			trim(adjustl(buffses)),",",&
			trim(adjustl(buffsec))
	end do
else
	
	write(unitbeta,'(5a)') "Beta,", "Estimate,", "Std. Err.,", "Louis S.E.,", "Sandwich S.E."
	do row=1,size(beta)
		write(buffe,'(g25.16e3)') beta(row)
		write(buffse,'(g25.16e3)') sqrt(max(0.0_8,eacovbeta(row,row)))
		write(buffsel,'(g25.16e3)') sqrt(max(0.0_8,eacovbeta_louis(row,row)))
		write(buffses,'(g25.16e3)') sqrt(max(0.0_8,eacovbeta_sandwich(row,row)))
		write(unitbeta,'(9a)') trim(adjustl(betaname(row))),",",&
			trim(adjustl(buffe)),",",&
			trim(adjustl(buffse)),",",&
			trim(adjustl(buffsel)),",",&
			trim(adjustl(buffses))
	end do
end if
return

end subroutine betaoutput

