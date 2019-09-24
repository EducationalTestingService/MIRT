!	Find estimated conditional expectations and conditional covariance matrices for weighted sums.
!	Responses are in dat.
!	numcat provides ranges for underlying observations.
!	numcatobs provides ranges for observations.
!	mask indicates items to use for computations.
!	beta is the parameter vector.
!	lintran is the linear transformation for the latent vector.
!	postdensity contains posterior weights.
!	theta contains posterior quadrature points.
!	wtsum is a matrix used to compute linear combinations of the observations.
!	covsum is the conditional covariance of wtsum.
!	meansum is the conditional mean of wtsum.


subroutine eapwtsum(dat,numcat,numcatobs,mask,beta,lintran,postdensity,theta,wtsum,covsum,meansum)
implicit none
interface




!	Multinomial covariance matrix.
!	density is vector of probabilities.
!	meanth is vector of means.
!	quadpoint is array of weights.
	function mcovth(density,meanth,quadpoint)
		implicit none
		real(kind=8),intent(in)::density(:),meanth(:),quadpoint(:,:)
		real(kind=8)::mcovth(size(meanth),size(meanth))
	end function mcovth
!	Multinomial mean.
!	density is vector of probabilities
!	quadpoint is array of weights.
	function mmeanth(density,quadpoint)
		implicit none
		real(kind=8),intent(in)::density(:),quadpoint(:,:)
		real(kind=8)::mmeanth(size(quadpoint,1))
	end function mmeanth
!	probvec is used to compute underlying item probabilities.
!   numcat indicates the number of underlying categories per item.
!	mask indicates which items were presented.
!	locations provides location parameters and scales provides scale factors.
!	theta is the value of the latent vector.
	function probvec(numcat,mask,locations,scales,theta)
		implicit none
		integer,intent(in)::numcat(:)
		logical,intent(in)::mask(:)
		real(kind=8),intent(in)::locations(:),scales(:,:),theta(:)
		real(kind=8)::probvec(size(locations))
	end function probvec

end interface
integer,intent(in)::dat(:,:),numcat(:),numcatobs(:)
logical,intent(in)::mask(:)
real(kind=8),intent(in)::beta(:),lintran(:,:),postdensity(:,:),theta(:,:,:),&
	wtsum(:,:)
real(kind=8),intent(out)::covsum(:,:,:),meansum(:,:)
!	item is an item.
!	ncat is the total number of underlying categories.
!	obs is observation number.
!	quad counts quadrature points.
!	resp is a single response.
!	row is a row counter.

integer::item,ncat,obs,quad,resp(size(mask)),row




!	locations is for location parameters.

!	meanquad is for means for quadrature points.
!	newtheta is lintran(theta).
!	probcat is used for underlying marginal conditional probabilities.
!	scales is for scale parameters.

real(kind=8)::locations(sum(numcat)),&
	meanquad(size(wtsum,1),size(postdensity,1)),newtheta(size(lintran,1)),&
	probcat(sum(numcat)),scales(size(lintran,1),sum(numcat))
	

!	Set up parameter arrays for simplified processing.
ncat=sum(numcat)
!	Item scale parameters.
locations=beta(1:ncat)
scales=reshape(beta(ncat+1:ncat*(1+size(lintran,1))),(/size(lintran,1),ncat/))




!	Obtain results for each observation.
do obs=1,size(postdensity,2)
	

	

		
!	Cycle through quadrature points.
	meanquad=0.0_8
	do quad=1,size(postdensity,1)
		newtheta=matmul(lintran,theta(:,quad,obs))

		
! Probability computation.


		probcat=probvec(numcat,mask,locations,scales,newtheta)
        
		if(size(wtsum,1)>0)then
			row=1
        
			do item=1,size(mask)
				if(mask(item)) meanquad(:,quad)=meanquad(:,quad)+&
                    mmeanth(probcat(row:row+numcat(item)-1),wtsum(:,row:row+numcat(item)-1))
				row=row+numcat(item)
			end do
			
		end if
	end do

!	mean and covariance.
	meansum(:,obs)=mmeanth(postdensity(:,obs),meanquad)
    covsum(:,:,obs)=mcovth(postdensity(:,obs),meansum(:,obs),meanquad)
    
	
end do
return
end subroutine eapwtsum
