!	gamma output
!	paramname is the vector of parameter names.
!	unitparam is the unit number.
!	complex indicates complex sampling.
!	eacovgam is the array of estimated covariances of gamma estimates under the model.
!	eacovgam_louis is the array of estimated covariances of gamma estimates under the model with the Louis approach.
!	eacovgam_sandwich is the array of estimated covariances of gamma estimates without assuming the model.
!	eacovgam_complex is the array of estimated covariances of gamma estimates under complex sampling.
!	gamma is the vector of estimates.
subroutine gammaoutput(paramname,unitparam,complx,eacovgam,eacovgam_complex,eacovgam_louis,eacovgam_sandwich,gamma)
implicit none
character(len=64)::paramname(:)
integer,intent(in)::unitparam
logical,intent(in)::complx
real(kind=8),intent(in)::eacovgam(:,:),eacovgam_complex(:,:),eacovgam_louis(:,:),eacovgam_sandwich(:,:),gamma(:)
!	buffers for output.
character(len=25)::buffe,buffse,buffsec,buffsel,buffses
!	Count rows.
integer::row
write(unitparam,'(a)') "Parameter estimates and standard errors"

if(complx)then
	write(unitparam,'(6a)') "Parameter,", "Estimate,", "Std. Err.,", "Louis S.E.,", "Sandwich S.E.,","Complex S.E."
	do row=1,size(gamma)
		write(buffe,'(g25.16e3)') gamma(row)
		write(buffse,'(g25.16e3)') sqrt(max(0.0_8,eacovgam(row,row)))
		write(buffsec,'(g25.16e3)') sqrt(max(0.0_8,eacovgam_complex(row,row)))
		write(buffsel,'(g25.16e3)') sqrt(max(0.0_8,eacovgam_louis(row,row)))
		write(buffses,'(g25.16e3)') sqrt(max(0.0_8,eacovgam_sandwich(row,row)))
		write(unitparam,'(11a)') trim(adjustl(paramname(row))),",",&
			trim(adjustl(buffe)),",",&
			trim(adjustl(buffse)),",",&
			trim(adjustl(buffsel)),",",&
			trim(adjustl(buffses)),",",&
			trim(adjustl(buffsec))
	end do
else
	write(unitparam,'(5a)') "Parameter,", "Estimate,", "Std. Err.,", "Louis S.E.,", "Sandwich S.E."
	do row=1,size(gamma)
		write(buffe,'(g25.16e3)') gamma(row)
		write(buffse,'(g25.16e3)') sqrt(max(0.0_8,eacovgam(row,row)))
		write(buffsel,'(g25.16e3)') sqrt(max(0.0_8,eacovgam_louis(row,row)))
		write(buffses,'(g25.16e3)') sqrt(max(0.0_8,eacovgam_sandwich(row,row)))
		write(unitparam,'(9a)') trim(adjustl(paramname(row))),",",&
			trim(adjustl(buffe)),",",&
			trim(adjustl(buffse)),",",&
			trim(adjustl(buffsel)),",",&
			trim(adjustl(buffses))
	end do
end if
return

end subroutine gammaoutput

