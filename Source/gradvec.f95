!	gradient of log conditional probability of an observation.
!	rbetamap is the reduced beta map.
!	mask is the response mask.
!	condprobcat contains conditional item probabilities.
!	meanth is the mean of theta.
!	mean2th is the mean of products of elements of theta.
!	newtheta is the transformed latent vector.
!	probcat contains item probabilities.
!	theta is the latent vector.
!	thpr is the independent vector.
function gradvec(rbetamap,mask,condprobcat,meanth,mean2th,newtheta,probcat,theta,thpr)

implicit none
integer,intent(in)::rbetamap(:,:)
logical,intent(in)::mask(:)
real(kind=8),intent(in)::condprobcat(:),meanth(:),mean2th(:,:),newtheta(:),probcat(:),theta(:),thpr(:)
real(kind=8)::gradvec(size(rbetamap,2))
!	row is a counter.
integer::row


gradvec=0.0_8

do row=1,size(gradvec)

	if(rbetamap(1,row)>0)then
		if(mask(rbetamap(1,row)))then
			if(rbetamap(3,row)==0)then
!	Intercept
				gradvec(row)=condprobcat(rbetamap(2,row))-probcat(rbetamap(2,row))
			else
!	Slope
				gradvec(row)=newtheta(rbetamap(3,row))*(condprobcat(rbetamap(2,row))-probcat(rbetamap(2,row)))
			endif
		endif
	else
		if(rbetamap(4,row)==0)then
!	Linear
			gradvec(row)=thpr(rbetamap(5,row))*(theta(rbetamap(3,row))-meanth(rbetamap(3,row)))
		else
!	Quadratic
			gradvec(row)=thpr(rbetamap(5,row))*&
				(theta(rbetamap(3,row))*theta(rbetamap(4,row))-mean2th(rbetamap(3,row),rbetamap(4,row)))
		end if
	end if
end do
return
end function gradvec
