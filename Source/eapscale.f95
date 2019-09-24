!	Find estimated conditional expectations and conditional covariance matrices for transformations of integer-weighted sums.
!	catobsrange maps underlying to observed categories.
!	maxscore is the maximum score.
!	maxw provides maximum item weights.
!	minscore is the minimum score.
!	minw provides minimum item weights.
!	numberscales is the number of scale scores.
!	numcat provides ranges for underlying observations.
!	numcatobs provides ranges for observations.
!	weight is the weight.
!	beta is the parameter vector.
!	lintran is the linear transformation for the latent vector.
!	postdensity contains posterior weights.
!	scale is the matrix of scale transformations.
!	theta contains posterior quadrature points.

!	covscale is the conditional covariance of the scale transformation of the weighted sum.
!	meanscale is the conditional mean of the scale transformation of the weighted sum.


subroutine eapscale(catobsrange,maxscore,maxw,minscore,&
    minw,numberscales,numcat,numcatobs,weight,beta,&
	lintran,postdensity,scale,theta,covscale,meanscale)
implicit none
interface
!	Obtain the frequency distribution of a weighted sum of responses.
!	The weights are integers.
!	maxnw is the array of minimum scores for items.
!	minw is the array of maximum scores for items.
!	numcatobs is the number of observed item categories.

!	weight is the category weights.
!	prob is the vector of observed category probabilities.
!	dist is the probability distribution of the weighted sum.

	subroutine distwtsum(maxw,minw,numcatobs,weight,prob,dist)
		implicit none
		integer,intent(in)::maxw(:),minw(:),numcatobs(:),weight(:)

		real(kind=8),intent(in)::prob(:)
		real(kind=8),intent(out)::dist(sum(minw):sum(maxw))
	end subroutine distwtsum



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
!	probvecobsall is used to compute observed item probabilities given the latent vector.
!	catobsrange is the table of ranges of underlying categories per observed category.
!	numcat provides the number of underlying categories per item, and numcatobs provides the number
!		of observed categories per item.
!	mask indicates which items were presented.
!	probcat is the vector of underlying item probabilities.
	function probvecobsall(catobsrange,numcatobs,mask,probcat)
		implicit none
		integer,intent(in)::catobsrange(:,:),numcatobs(:)
		logical,intent(in)::mask(:)
		real(kind=8),intent(in)::probcat(:)
		real(kind=8)::probvecobsall(sum(numcatobs))
	end function probvecobsall
end interface
integer,intent(in)::catobsrange(:,:),maxscore,maxw(:),minscore,minw(:),numberscales,&
	numcat(:),numcatobs(:),weight(:)

real(kind=8),intent(in)::beta(:),lintran(:,:),postdensity(:,:),scale(numberscales,minscore:maxscore),theta(:,:,:)
real(kind=8),intent(out)::covscale(:,:,:),meanscale(:,:)
!	dimlatout is the dimension of the transformed latent vector.
!	ncat is the total number of underlying categories.
!	ncatobs is the total number of observed categories.
!	nscale is the location in beta of the initial item discrimination.
!	nscale1 is the location in beta of the last item discrimination.
!	obs is observation number.
!	quad counts quadrature points.
!	row is a row counter.

integer::dimlatout,item,ncat,ncatobs,nitems,nquad,nscale,nscale1,obs,quad,row
!   Figure out which items to include.
   logical::datamask(size(numcat))
!	dist is the array of probabilities of raw scores.
!	locations is for location parameters.


!	newtheta is lintran(theta).
!	probcat is used for underlying marginal conditional probabilities.
!	scales is for scale parameters.

real(kind=8)::dist(minscore:maxscore),&
	locations(sum(numcat)),&
	meanquad(numberscales,size(postdensity,1)),&
	newtheta(size(lintran,1)),&
	probcat(sum(numcat)),probcatobs(sum(numcatobs)),scales(size(lintran,1),sum(numcat))
	
	

!	Set up parameter arrays for simplified processing.
dimlatout=size(lintran,1)
ncat=sum(numcat)
ncatobs=sum(numcatobs)
!	Item scale parameters.
locations=beta(1:ncat)
nscale=ncat+1
nscale1=ncat*(1+dimlatout)
scales=reshape(beta(nscale:nscale1),(/dimlatout,ncat/))
datamask=.false.
do row=1,size(numcat)
    if(maxw(row)>minw(row)) datamask(row)=.TRUE.
end do




!	Obtain results for each observation.
do obs=1,size(postdensity,2)
	
	
	

		
!	Cycle through quadrature points.

	do quad=1,size(postdensity,1)
		newtheta=matmul(lintran,theta(:,quad,obs))
		
		
! Probability computation.


		probcat=probvec(numcat,datamask,locations,scales,newtheta)
        probcatobs=probvecobsall(catobsrange,numcatobs,datamask,probcat)

		call distwtsum(maxw,minw,numcatobs,weight,probcatobs,dist)
		meanquad(:,quad)=mmeanth(dist,scale)
		
	end do

!	mean and covariance.
	meanscale(:,obs)=mmeanth(postdensity(:,obs),meanquad)
    covscale(:,:,obs)=mcovth(postdensity(:,obs),meanscale(:,obs),meanquad)
    
	
end do
return
end subroutine eapscale
