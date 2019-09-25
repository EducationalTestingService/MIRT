!	Obtain interactions of item values and weighted sums.
!	catobsrange maps underlying to observed categories.
!	dat provides responses.

!	maxw provides maximum item weights.

!	minw provides minimum item weights.
!	npsu is the number of psus's per stratum.
!	nstratum is the number of strata.
!	numcat provides the underlying number of categories per item.
!	numcatobs provides the observed number of categories per item.
!	psu indicates the psu of an observation.
!	stratum indicates the stratum of an observation.
!	weight provides weighted sum.
!	complx indicates complex sampling.
!	resid is for residuals.
!	stratify is .true. for stratified sampling.
!	usepsu is .true. for primary sampling units within strata.
!	beta provides parameter estimates.

!	eacovgaminv_louis is the inverse of the Louis information plus the
!		constraint component of the negative hessian.
!	gradsdes provides gradients for observations.
!	lintran transforms the latent vector.
!	obsweight provides weights.
!	postdensity is the array of posterior densities.
!	theta is the array of quadrature points.
!	tolres is the residual tolerance.
!	fitpwtitem is the average of fitted product of item indicator and weighted sum.
!	fitwtitem is the total of fitted products of item indicator and weighted sum.
!	obspwtitem is the observed average product of item indicator and weighted sum.
!	obswtitem is the observed total products of item indicator and weighted sum.
!	presentedwtitem counts weighted items presented with defined weighted sums.
!	residawtitem is the adjusted residual of total of products.
!	residpwtitem is the residual average of products.
!	residwtitem is the residual total of products.
!	stdobspwtitem is the asymptotic standard deviation of the observed average of products.
!	stdobswtitem is the asymptotic standard deviation of the observed total of products.
!	stdresidpwtitem is the asymptotic standard deviation of residual average of products.
!	stdresidwtitem is the asymptotic standard deviation of the residual total of products.



subroutine marginwtitem(catobsrange,dat,maxw,minw,npsu,nstratum,numcat,numcatobs,&
	psu,stratum,weight,complx,resid,stratify,usepsu,&
	beta,eacovgaminv_louis,gradsdes,lintran,&
	obsweight,postdensity,theta,tolres,&
	fitpwtitem,fitwtitem,obspwtitem,obswtitem,presentedwtitem,&
	residawtitem,residpwtitem,residwtitem,stdobspwtitem,stdobswtitem,stdresidpwtitem,stdresidwtitem)
	

implicit none
interface


!	Estimate the cross information matrix.
!	adjust indicates adjustment for means.
!	grads is the array of observation gradients.
!	grads1 is another array of observation gradients.
!	obsweight contains observation weights.

	function crossinformation(adjust,grads,grads1,obsweight)
		implicit none
		logical,intent(in)::adjust
		real(kind=8),intent(in)::grads(:,:),grads1(:,:),obsweight(:)
		real(kind=8)::crossinformation(size(grads,1),size(grads1,1))
	end function crossinformation


	
!	The following description is for the function's basic purpose.
!	Its use in this subroutine is a bit different.
!	Estimate the information matrix for a complex sample.
!	This version is designed for random sampling of psu's within strata.
!	Sampling is with replacement.
!	npsu is the number of psus's per stratum.
!	nstratum is the number of strata.
!	psu indicates the psu of an observation.
!	stratum indicates the stratum of an observation.
!	stratify is .true. for stratified sampling.
!	usepsu is .true. for primary sampling units within strata.
!	grads is the array of observation gradients.
!	obsweight contains observation weights.

	function information_complex(npsu,nstratum,psu,stratum,stratify,usepsu,grads,obsweight)
		implicit none
		integer,intent(in)::nstratum
		integer,intent(in),optional::npsu(:),psu(:),stratum(:)
		logical,intent(in)::stratify,usepsu
		real(kind=8),intent(in)::grads(:,:),obsweight(:)
		real(kind=8)::information_complex(size(grads,1),size(grads,1))
	end function information_complex

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



end interface
integer,intent(in)::catobsrange(:,:),dat(:,:),maxw(:),minw(:),&
    npsu(:),nstratum,numcat(:),numcatobs(:),&
	psu(:),stratum(:),weight(:)
logical,intent(in)::complx,resid,stratify,usepsu
real(kind=8),intent(in)::beta(:),eacovgaminv_louis(:,:),&
	gradsdes(:,:),lintran(:,:),obsweight(:),&
	postdensity(:,:),theta(:,:,:),tolres
real(kind=8),intent(out)::fitpwtitem(:),fitwtitem(:),&
	obspwtitem(:),obswtitem(:),presentedwtitem(:),residawtitem(:),residpwtitem(:),residwtitem(:),&
    stdobspwtitem(:),stdobswtitem(:),stdresidpwtitem(:),stdresidwtitem(:)
	

!	al is the allocation indicator.
!	counter finds observed item codes.
!	counter1 also finds observed item codes.
!	dimdesign is the design dimension.
!	dimlatout is the number of skills.
!	item is an item.
!	ncat is number of underlying categories.
!	ncatobs is number of observed categories.
!	nitems is the number of items.
!	nobs is the number of observations.
!	nquad is the number of quadrature points.
!	nscale is the initial position in beta of slopes.
!	nscale1 is the end position in beta of slopes.
!	obs counts observations.
!	quad counts quadrature points.
!	resp is a response vector.
!	score is the weighted sum for the observation.
integer::al,counter,counter1,dimdesign,dimlatout,item,ncat,ncatobs,nitems,nobs,nquad,nscale,nscale1,&
	obs,quad,resp(size(dat,1)),score,score1
!	datamask is a mask for items needed for weighted sum.
!	datamask1 is a mask for items presented.


logical::datamask(size(dat,1)),datamask1(size(dat,1))

!   avediff is an average difference.
!	covobs is the covariance matrix for adjusted observed values.
!	covsum is used to estimate asymptotic variances for fitted marginals.
!   diff is a difference.


!   fr is a fraction.
!	locations gives location parameters.

!	newtheta is the transformed latent vector.
!	obsmat is the array of response indicators and fits.
!		The first row is the observed indicator, the second is the fit,
!		the third is 1 if presented and 0 otherwise.
!	probcat is the vector of underlying conditional probabilities.
!	probcatobs is the vector of observed conditional probabilities.
!   probm is the array of observed marginal probabilities.
!	probs is the array of observed conditional probabilities.
!	respprobvec is the vector of fitted probabilites for an observation.
!	respresid is the residual component for an observation.
!	respvec is the corresponding vector of observations.
!	scales gives scale parameters.


!	slopes is for regression of errors on gradients. 
!	totalweight is sum of all observation weights.
!   x is a real number.

real(kind=8)::avediff,covobs(1,1),covsum(size(gradsdes,1)),diff,&
	fr,locations(sum(numcat)),newtheta(size(lintran,1)),&
	probcat(sum(numcat)),&
	probm(sum(numcatobs)),probs(sum(numcatobs),size(theta,2)),&
    respprobvec(sum(numcatobs)),respvec(sum(numcatobs)),&
	scales(size(lintran,1),sum(numcat)),&
	slopes(size(gradsdes,1)),totalweight,x
real(kind=8),allocatable::obsmat(:,:,:)	

!	Set up parameter arrays for simplified processing.
dimdesign=size(gradsdes,1)

dimlatout=size(lintran,1)

ncat=sum(numcat)
ncatobs=sum(numcatobs)
nitems=size(dat,1)
nobs=size(obsweight)
nquad=size(theta,2)
nscale=ncat+1
nscale1=ncat*(dimlatout+1)
datamask=.false.

fitwtitem=0.0_8
fitpwtitem=0.0_8
locations=beta(1:ncat)
allocate(obsmat(3,nobs,ncatobs),stat=al)
if(al/=0)stop 'Allocation failed for interactions of items and weighted sums.'
obsmat=0.0_8
obswtitem=0.0_8
obspwtitem=0.0_8
presentedwtitem=0.0_8
stdobswtitem=0.0_8
stdobspwtitem=0.0_8
totalweight=sum(obsweight)
if(resid)then
	residwtitem=0.0_8
	residpwtitem=0.0_8
	residawtitem=0.0_8
	stdresidwtitem=0.0_8
	stdresidpwtitem=0.0_8
end if
scales=reshape(beta(nscale:nscale1),(/dimlatout,ncat/))

do item=1,nitems
	if(maxw(item)>minw(item))datamask(item)=.true.
end do
!	Totals of products.
! Observation obs.

do obs=1,nobs
	if(obsweight(obs)<=0.0_8) cycle
!	Establish mask for items presented.
	datamask1=.true.
	resp=dat(:,obs)
	respvec=0.0_8
	respprobvec=0.0_8
!
	
	

	
	

	
	
	
	do item=1,nitems
		
		if(resp(item)<0.or.resp(item)>=numcatobs(item)) datamask1(item)=.false.
			
	end do
!	Check for weighted sum that is not observed.
	if(any(datamask.and.(.not.datamask1)))cycle
!	Weighted sum for observation.
	score=0
	counter=1
	!	Cycle through quadrature points.
	probm=0.0_8
	do quad=1,nquad
!		Probability computation.
		newtheta=matmul(lintran,theta(:,quad,obs))
		probcat=probvec(numcat,datamask1,locations,scales,newtheta)
		probs(:,quad)=probvecobsall(catobsrange,numcatobs,datamask1,probcat)
		probm=probm+postdensity(quad,obs)*probs(:,quad)
	end do
!	Observed and expected item scores and total score.
	
	
	do item=1,nitems
		if(datamask(item))then
			score=score+weight(counter+resp(item))
	!	Expected scores for items.
			
			
		else
			score=score+minw(item)
			
			
		end if
		
		
		counter=counter+numcatobs(item)
	end do


!	Go through all items.
	counter=1
	respvec=0.0_8
	respprobvec=0.0_8
	do item=1,nitems
		
!	Skip item if missing.
		if(datamask1(item))then
            presentedwtitem(item)=presentedwtitem(item)+obsweight(obs)
            obsmat(3,obs,counter:counter+numcatobs(item)-1)=1.0_8
			score1=score-weight(counter+resp(item))
            respvec(counter+resp(item))=score1
            obswtitem(counter+resp(item))=obswtitem(counter+resp(item))+obsweight(obs)*score1
            do counter1=counter,counter+numcatobs(item)-1
				respprobvec(counter1)=score1*probm(counter1)
				fitwtitem(counter1)=fitwtitem(counter1)+respprobvec(counter1)*obsweight(obs)
				
            end do
		end if
		counter=counter+numcatobs(item)
	end do
    
    obsmat(1,obs,:)=respvec
	if(resid)obsmat(2,obs,:)=respvec-respprobvec

end do
if(resid) residwtitem=obswtitem-fitwtitem
!
!	Fractions and variances.
counter=1

do item=1,nitems
    fr=presentedwtitem(item)/totalweight
	if(presentedwtitem(item)>0.0_8) then
		do counter1=counter,counter+numcatobs(item)-1
			fitpwtitem(counter1)= fitwtitem(counter1)/presentedwtitem(item)
			obspwtitem(counter1)= obswtitem(counter1)/presentedwtitem(item)
			if(resid)residpwtitem(counter1)= residwtitem(counter1)/presentedwtitem(item)
				
			
				

	
!	Standard deviation of total of products depends on sampling.
			if(.not.complx)then
                x=obswtitem(counter1)/totalweight
                do obs=1,nobs
                    diff=obsmat(1,obs,counter1)-x
                    stdobswtitem(counter1)=stdobswtitem(counter1)+obsweight(obs)*diff*diff
                    diff=diff-obspwtitem(counter1)*(obsmat(3,obs,counter1)-fr)
                    stdobspwtitem(counter1)=stdobspwtitem(counter1)+obsweight(obs)*diff*diff

                end do
                stdobswtitem(counter1)=sqrt(stdobswtitem(counter1))
                stdobspwtitem(counter1)=sqrt(stdobspwtitem(counter1))/presentedwtitem(item)
				if(resid)then

                    
                    covsum=reshape(crossinformation(.true.,obsmat(2:2,:,counter1),gradsdes,obsweight),&
                        (/dimdesign/))
                       
                    slopes=matmul(eacovgaminv_louis,covsum)
                    avediff=0.0_8
                    do obs=1,nobs
                        obsmat(2,obs,counter1)=obsmat(2,obs,counter1)-dot_product(gradsdes(:,obs),slopes)
                        avediff=avediff+obsmat(2,obs,counter1)*obsweight(obs)
                    end do
                    avediff=avediff/totalweight
                    do obs=1,nobs
                        diff=obsmat(2,obs,counter1)-avediff
                        stdresidwtitem(counter1)=stdresidwtitem(counter1)+diff*diff*obsweight(obs)
                    end do

                    
                end if
				
					
				
			else
                covobs=information_complex(npsu,nstratum,psu,stratum,stratify,usepsu,&
                    obsmat(1:1,:,counter1),obsweight)
                stdobswtitem(counter1)=sqrt(covobs(1,1))

                obsmat(1,:,counter1)=obsmat(1,:,counter1)-obspwtitem(counter1)*obsmat(3,:,counter1)
                covobs=information_complex(npsu,nstratum,psu,stratum,stratify,usepsu,&
                    obsmat(1:1,:,counter1),obsweight)
                stdobspwtitem(counter1)=sqrt(covobs(1,1)/presentedwtitem(item))
                if(resid)then
					covsum=reshape(crossinformation(.true.,obsmat(2:2,:,counter1),gradsdes,obsweight),&
                        (/dimdesign/))
                       
                    slopes=matmul(eacovgaminv_louis,covsum)
				

					
                    do obs=1,nobs
                        obsmat(2,obs,counter1)=obsmat(2,obs,counter1)-dot_product(gradsdes(:,obs),slopes)
                    end do


                    covobs=information_complex(npsu,nstratum,psu,stratum,stratify,usepsu,&
                        obsmat(2:2,:,counter1),obsweight)

                    stdresidwtitem(counter1)=covobs(1,1)


				end if
				
			end if
			if(resid)then
				stdresidwtitem(counter1)=sqrt(max(0.0_8,stdresidwtitem(counter1)))
				stdresidpwtitem(counter1)=stdresidwtitem(counter1)/presentedwtitem(item)
                if(stdresidwtitem(counter1)>tolres*stdobswtitem(counter1).and.stdobswtitem(counter1)>0.0_8)&
                    residawtitem(counter1)=residwtitem(counter1)/stdresidwtitem(counter1)
			end if
			
		
			
		end do
	end if




	counter=counter+numcatobs(item)

end do 

return
deallocate(obsmat)
end subroutine marginwtitem
