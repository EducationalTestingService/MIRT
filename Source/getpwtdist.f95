!	Obtain lower and upper probabilities of a weighted distractor sum.
!	catobsrange maps underlying to observed categories.
!   choices indicates the number of distractors per item.
!	dat provides responses.
!   distdata provides distractor data.
!   distmap is the distractor map.
!	maxscore is the maximum score.
!	maxw provides maximum item weights.
!	minscore is the minimum score.
!	minw provides minimum item weights.
!	numcat provides the underlying number of categories per item.
!	numcatobs provides the observed number of categories per item.
!	weightdist is the weights.
!	beta provides parameter estimates.
!	lintran transforms the latent vector.
!	postdensity is the array of posterior densities.
!	theta is the array of quadrature points.


!	plower is the lower probabilities.
!   pupper is the upper probabilities.

!
subroutine getpwtdist(catobsrange,choices,dat,distdat,distmap,&
    maxscore,maxw,minscore,minw,&
	numcat,numcatobs,weightdist,&
    beta,distprob,&
	lintran,&
	postdensity,theta,&
	plower,pupper)

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

	


!	probvec is used to compute underlying item probabilities.
!   	numcat indicates the number of underlying categories per item.
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

!	probvecobsdist is used to compute observed item distractor probabilities given the latent vector.
!	choices is the number of choices per item.
!	numcatobs provides the number of observed categories per item.
!	mask indicates which items were presented.
!	distprob provides distractor probabilities.
!	probcatobs is the vector of item probabilities.
    function probvecobsdist(choices,distmap,numcatobs,mask,distprob,probcatobs)
        implicit none
        integer,intent(in)::choices(:),distmap(:),numcatobs(:)
        logical,intent(in)::mask(:)
        real(kind=8),intent(in)::distprob(:),probcatobs(:)
        real(kind=8)::probvecobsdist(sum(choices))
    end function probvecobsdist

end interface
integer,intent(in)::catobsrange(:,:),choices(:),dat(:,:),distdat(:,:),distmap(:)
integer,intent(in)::maxscore,maxw(:),minscore,minw(:)
integer,intent(in)::numcat(:),numcatobs(:),weightdist(:)
real(kind=8),intent(in)::beta(:),distprob(:),lintran(:,:),&
	postdensity(:,:),theta(:,:,:)
real(kind=8),intent(out)::plower(:),pupper(:)
	
!	col counts columns.
!	counter finds underlying item codes.
!	dimdesign is the design dimension.
!	dimlatout is the number of skills.
!	item is an item.
!	ncat is number of underlying categories.
!	ncatobs is number of observed categories.
!	nitems is the number of items.

!	nquad is the number of quadrature points.
!	nscale is the position in beta of slopes.
!	nscale1 is the position in beta of the last slope.
!	obs counts observations.
!	quad counts quadrature points.
!	resp is a response vector.

!	score is a weighted sum
integer::counter,dimlatout,item,&
	ncat,ncatobs,nitems,nquad,nscale,nscale1,&
	obs,quad,resp(size(dat,1)),respdist(size(dat,1)),score
!	datamask is a mask for items with variable weights.
!	skip indicates that an item is to be skipped.
logical::datamask(size(dat,1)),skip

!	dist is the conditional distribution.
!	locations gives location parameters.
!	newtheta is the transformed latent vector.

!	probcat is the vector of underlying probabilities.
!	probcatobs is the vector of observed probabilities.
!	respprobvec is the vector of fitted probabilites for an observation.

!	respvec is the corresponding vector of observations.
!	scales gives scale parameters.


real(kind=8)::dist(minscore:maxscore),&
	locations(sum(numcat)),newtheta(size(lintran,1)),&
	probcat(sum(numcat)),probcatobs(sum(numcatobs)),&
    probcatobsdist(sum(choices)),&
	respprobvec(minscore:maxscore),&
	scales(size(lintran,1),sum(numcat))
!	Set up parameter arrays for simplified processing.


dimlatout=size(lintran,1)
ncat=sum(numcat)
ncatobs=sum(numcatobs)
nitems=size(dat,1)

nquad=size(theta,2)
nscale=ncat+1
nscale1=ncat*(dimlatout+1)
locations=beta(1:ncat)




scales=reshape(beta(nscale:nscale1),(/dimlatout,ncat/))


!	Establish mask for items with variable weights.
datamask=.false.

do item=1,nitems
	if(maxw(item)>minw(item))datamask(item)=.true.
end do

plower=0.0_8
pupper=0.0_8
!	Marginal totals.
do obs=1,size(dat,2)

! Observation obs.
	resp=dat(:,obs)
    respdist=distdat(:,obs)
	respprobvec=0.0_8
	skip=.false.
	
	do item=1,nitems
!	Verify if weighted sum observed.
		if(datamask(item))then
			if(resp(item)<0.or.resp(item)>=numcatobs(item))then
				skip=.true.
				exit
			end if
		end if
	end do
	if(skip)cycle
	
	score=0
	counter=1
	do item=1,nitems
		if(datamask(item))then
			score=score+weightdist(counter+respdist(item))
		else
			score=score+minw(item)
		end if
		counter=counter+choices(item)
	end do





!	Cycle through quadrature points.
	do quad=1,nquad
!		Probability computation.
		newtheta=matmul(lintran,theta(:,quad,obs))
		probcat=probvec(numcat,datamask,locations,scales,newtheta)
        probcatobs=probvecobsall(catobsrange,numcatobs,datamask,probcat)
        probcatobsdist=probvecobsdist(choices,distmap,numcatobs,datamask,&
            distprob,probcatobs)




        call distwtsum(maxw,minw,choices,weightdist,&
            probcatobsdist,dist)
		respprobvec=respprobvec+postdensity(quad,obs)*dist
	end do

!   Upper and lower probabilities.
	plower(obs)=sum(respprobvec(minscore:score))
	pupper(obs)=sum(respprobvec(score:maxscore))
end do




return

end subroutine getpwtdist

