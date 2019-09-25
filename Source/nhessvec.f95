!	negative hessian of log conditional probability of an observation.
!	rbetamap is the reduced beta map
!	mask indicates the items presented.
!	condprobcat contains conditional item probabilities.
!	covth is the covariance matrix of the latent vector.
!	cubth provides covariances of products of latent vector elements and of latent vector elements.
!	probcat contains item probabilities.
!	quarth provides covariances of products of latent vector elements.
!	theta is the latent vector.
!	thpr is the independent vector.

function nhessvec(rbetamap,mask,condprobcat,covth,cubth,newtheta,probcat,quarth,theta,thpr)

implicit none
integer,intent(in)::rbetamap(:,:)
logical,intent(in)::mask(:)
real(kind=8),intent(in)::condprobcat(:),covth(:,:),cubth(:,:,:),newtheta(:),probcat(:),quarth(:,:,:,:),theta(:),thpr(:)
real(kind=8)::nhessvec(size(rbetamap,2),size(rbetamap,2))





!	col is a column.
!	row is a row, 
integer::row,col


nhessvec=0.0_8
do row=1,size(rbetamap,2)
	do col=1,row

		if(rbetamap(1,row)>0)then
			if(mask(rbetamap(1,row)))then
				if(rbetamap(1,row)==rbetamap(1,col))then
!	Intercepts.
					nhessvec(row, col)=condprobcat(rbetamap(2,row))*condprobcat(rbetamap(2,col))&
						-probcat(rbetamap(2,row))*probcat(rbetamap(2,col))

					if(rbetamap(2,row)==rbetamap(2,col))nhessvec(row,col)=&
						nhessvec(row,col)+probcat(rbetamap(2,row))-condprobcat(rbetamap(2,row))
					if(rbetamap(3,row)>0)then
!	Slope and intercept
						nhessvec(row,col)=newtheta(rbetamap(3,row))*&
							nhessvec(row,col)					
					end if
!	Slope and slope
					if(rbetamap(3,col)>0)then
						nhessvec(row, col)=newtheta(rbetamap(3,col))*nhessvec(row,col)
					end if
				end if
			end if
		else
			if(rbetamap(4,row)==0)then
!	Linear and linear
				if(rbetamap(1,col)==0)then
					nhessvec(row,col)=thpr(rbetamap(5,row))*thpr(rbetamap(5,col))*&
						covth(rbetamap(3,row),rbetamap(3,col))
				end if
			else
				if(rbetamap(1,col)==0)then
					if(rbetamap(4,col)==0)then
!	Quadratic by linear						
						nhessvec(row,col)=thpr(rbetamap(5,row))*thpr(rbetamap(5,col))*&
							cubth(rbetamap(3,row),rbetamap(4,row),rbetamap(3,col))	
					else
!	Quadratic by quadratic
						nhessvec(row,col)=thpr(rbetamap(5,row))*thpr(rbetamap(5,col))*&
							quarth(rbetamap(3,row),rbetamap(4,row),rbetamap(3,col),rbetamap(4,col))
					end if
				end if
			end if
		end if
	end do
end do
return
end function nhessvec
