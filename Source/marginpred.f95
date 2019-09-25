!	Obtain interactions of category indicators and predictors.
!	catobsrange maps underlying to observed categories.
!	dat provides responses.
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
!	extvar provides external variables.
!	gradsdes provides gradients for observations.
!	indvar provides predictor variables.
!	lintran transforms the latent vector.
!	obsweight provides weights.
!	postdensity is the array of posterior densities.
!	theta is the array of quadrature points.
!	tolres is the residual tolerance.
!	fitppred is the average of fitted product of item indicator and predictors.
!	fitpred is the total of fitted products of item indicator and predictors.
!	obsppred is the observed average product of item indicator and predictors.
!	obspred is the observed total product of item indicator and predictors.
!	presented counts weighted items presented.
!	residapred is the adjusted residual of average of products.
!	residppred is the residual of average of products.
!	residpred is the residual of total of products.
!	stdobsppred is the asymptotic standard deviation of the observed average of products.
!	stdobspred is the asymptotic standard deviation of the observed total of products.
!	stdresidppred is the asymptotic standard deviation of residual average of products.
!	stdresidpred is the asymptotic standard deviation of the residual total of products.



subroutine marginpred(catobsrange,dat,npsu,nstratum,numcat,numcatobs,&
	psu,stratum,complx,resid,stratify,usepsu,&
	beta,eacovgaminv_louis,extvar,indvar,gradsdes,lintran,&
	obsweight,postdensity,theta,tolres,&
	fitppred,fitpred,obsppred,obspred,presented,&
    residapred,residppred,residpred,stdobsppred,stdobspred,stdresidppred,stdresidpred)
	

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
integer,intent(in)::catobsrange(:,:),dat(:,:),&
    npsu(:),nstratum,numcat(:),numcatobs(:),&
	psu(:),stratum(:)
logical,intent(in)::complx,resid,stratify,usepsu
real(kind=8),intent(in)::beta(:),eacovgaminv_louis(:,:),&
	gradsdes(:,:),extvar(:,:),indvar(:,:),lintran(:,:),obsweight(:),&
	postdensity(:,:),theta(:,:,:),tolres
real(kind=8),intent(out)::fitppred(:,:),fitpred(:,:),&
	obsppred(:,:),obspred(:,:),presented(:),residapred(:,:),residppred(:,:),residpred(:,:),&
    stdobsppred(:,:),stdobspred(:,:),stdresidppred(:,:),stdresidpred(:,:)
	

!	al is an allocation indicator.
!	counter finds observedlying item codes.
!	counter1 also finds observed item codes.
!	dimdesign is the design dimension.
!	dimlatout is the number of skills.
!	item is an item.
!	ncat is number of underlying categories.
!	ncatobs is number of observed categories.
!	nitems is the number of items.
!	nlin is the position in beta of linear terms.
!	nobs is the number of observations.
!	nprede is the number of external predictors.
!	npredi is the number of internal predictors.
!	npred is the number of predictors.
!	nquad is the number of quadrature points.
!	nscale is the position in beta of slopes.
!	obs counts observations.
!	pred is a predictor number.
!	quad counts quadrature points.
!	resp is a response vector.

integer::al,counter,counter1,dimdesign,dimlatout,item,ncat,ncatobs,nitems,nlin,nobs,&
	npred,nprede,npredi,nquad,nscale,&
	obs,pred,quad,resp(size(dat,1))
!	datamask is a mask for items presented.

logical::datamask(size(dat,1))

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
!	score is a predictor value.

!	slopes is for regression of errors on gradients. 
!	totalweight is the sum of all observation weights.
!   x is a real number.

real(kind=8)::avediff,covobs(1,1),covsum(size(gradsdes,1)),diff,escore(size(numcatobs),size(theta,2)),&
	fitscore(size(numcatobs)),fitscoret,fr,locations(sum(numcat)),newtheta(size(lintran,1)),&
	probcat(sum(numcat)),&
	probm(sum(numcatobs)),probs(sum(numcatobs)),&
    respprobvec(sum(numcatobs)),respvec(sum(numcatobs)),&
	scales(size(lintran,1),sum(numcat)),score,&
	slopes(size(gradsdes,1)),totalweight,x
	
real(kind=8),allocatable::obsmat(:,:,:)
!	Set up parameter arrays for simplified processing.
dimdesign=size(gradsdes,1)

dimlatout=size(lintran,1)

ncat=sum(numcat)
ncatobs=sum(numcatobs)
nitems=size(dat,1)
nlin=ncat*(dimlatout+1)+1
nobs=size(obsweight)
npredi=size(indvar,1)
nprede=size(extvar,1)
npred=nprede+npredi
nquad=size(theta,2)
nscale=ncat+1

fitppred=0.0_8
fitpred=0.0_8
locations=beta(1:ncat)
allocate(obsmat(3,nobs,ncatobs),stat=al)
if(al/=0)stop 'Allocation failed for interactions of items and predictors.'
obsmat=0.0_8

obsppred=0.0_8
obspred=0.0_8

stdobsppred=0.0_8
stdobspred=0.0_8
totalweight=sum(obsweight)
if(resid)then
	residapred=0.0_8
	residppred=0.0_8
	residpred=0.0_8
	stdresidppred=0.0_8
	stdresidpred=0.0_8
end if
scales=reshape(beta(nscale:nlin-1),(/dimlatout,ncat/))
presented=0.0_8
!	Totals of products.
do pred=2,npred

 do obs=1,nobs
	if(obsweight(obs)<=0.0_8) cycle
! Observation obs.
	resp=dat(:,obs)
    datamask=.false.
    do item=1,nitems

!	Skip item if missing.
		if(resp(item)>=0.and.resp(item)<numcatobs(item))then
            datamask(item)=.true.
            if(pred==2) presented(item)=presented(item)+obsweight(obs)
		end if
	end do
    
	
	if(pred<=npredi)then
		score=indvar(pred,obs)
	else
		score=extvar(pred-npredi,obs)
	end if
	
	!	Cycle through quadrature points.
	probm=0.0_8
	do quad=1,nquad
!		Probability computation.
		newtheta=matmul(lintran,theta(:,quad,obs))
		probcat=probvec(numcat,datamask,locations,scales,newtheta)
		probs=probvecobsall(catobsrange,numcatobs,datamask,probcat)
		probm=probm+postdensity(quad,obs)*probs
	end do

!	Go through all items.
	counter=1
	respvec=0.0_8
	respprobvec=0.0_8
	do item=1,nitems

!	Skip item if missing.
		if(datamask(item))then
			
            
            obsmat(3,obs,counter:counter+numcatobs(item)-1)=1.0_8
            respvec(counter+resp(item))=score
            obspred(counter+resp(item),pred)=obspred(counter+resp(item),pred)+obsweight(obs)*score
            do counter1=counter,counter+numcatobs(item)-1
                respprobvec(counter1)=probm(counter1)*score
                fitpred(counter1,pred)=fitpred(counter1,pred)+respprobvec(counter1)*obsweight(obs)
            end do
		end if
		counter=counter+numcatobs(item)
	end do
    
    obsmat(1,obs,:)=respvec
	if(resid)obsmat(2,obs,:)=respvec-respprobvec

end do

if(resid) residpred(:,pred)=obspred(:,pred)-fitpred(:,pred)
!
!	Fractions and variances.
counter=1

do item=1,nitems

    fr=presented(item)/totalweight
	if(presented(item)>0.0_8) then
		do counter1=counter,counter+numcatobs(item)-1
			fitppred(counter1,pred)= fitpred(counter1,pred)/presented(item)
			obsppred(counter1,pred)= obspred(counter1,pred)/presented(item)
			if(resid)residppred(counter1,pred)= residpred(counter1,pred)/presented(item)
				
			
				

	
!	Standard deviation of total of products depends on sampling.
			if(.not.complx)then
                x=obspred(counter1,pred)/totalweight
                do obs=1,nobs
                    diff=obsmat(1,obs,counter1)-x
                    stdobspred(counter1,pred)=stdobspred(counter1,pred)+obsweight(obs)*diff*diff
                    diff=diff-obsppred(counter1,pred)*(obsmat(3,obs,counter1)-fr)
                    stdobsppred(counter1,pred)=stdobsppred(counter1,pred)+obsweight(obs)*diff*diff

                end do
                stdobspred(counter1,pred)=sqrt(stdobspred(counter1,pred))
                stdobsppred(counter1,pred)=sqrt(stdobspred(counter1,pred))/presented(item)
				if(resid)then

                    
                    covsum=reshape(crossinformation(.true.,obsmat(2:2,:,counter1),gradsdes,obsweight),&
                        (/dimdesign/))
                   
                    slopes=matmul(eacovgaminv_louis,covsum)
					avediff=0.0_8
                do obs=1,nobs
                    if(obsweight(obs)<=0.0_8)cycle
                    obsmat(2,obs,counter1)=obsmat(2,obs,counter1)-dot_product(gradsdes(:,obs),slopes)
                    avediff=avediff+obsmat(2,obs,counter1)*obsweight(obs)
                end do
                avediff=avediff/totalweight
                do obs=1,nobs
                    diff=obsmat(2,obs,counter1)-avediff
                    stdresidpred(counter1,pred)=stdresidpred(counter1,pred)+diff*diff*obsweight(obs)
                end do

                                      
                    
                end if
				
					
				
			else
                covobs=information_complex(npsu,nstratum,psu,stratum,stratify,usepsu,&
                    obsmat(1:1,:,counter1),obsweight)
                stdobspred(counter1,pred)=sqrt(covobs(1,1))

                obsmat(1,:,counter1)=obsmat(1,:,counter1)-obsppred(counter1,pred)*obsmat(3,:,counter1)
                covobs=information_complex(npsu,nstratum,psu,stratum,stratify,usepsu,&
                    obsmat(1:1,:,counter1),obsweight)
                stdobsppred(counter1,pred)=sqrt(covobs(1,1)/presented(item))
                if(resid)then
					covsum=reshape(crossinformation(.true.,obsmat(2:2,:,counter1),gradsdes,obsweight),&
                        (/dimdesign/))
                   
                    slopes=matmul(eacovgaminv_louis,covsum)
					

					
                    do obs=1,nobs
                        obsmat(2,obs,counter1)=obsmat(2,obs,counter1)-dot_product(gradsdes(:,obs),slopes)
                    end do
                    covobs=information_complex(npsu,nstratum,psu,stratum,stratify,usepsu,&
                        obsmat(2:2,:,counter1),obsweight)

					stdresidpred(counter1,pred)=covobs(1,1)
				end if
				
			end if
			if(resid)then
				
				stdresidpred(counter1,pred)=sqrt(max(0.0_8,stdresidpred(counter1,pred)))
				stdresidppred(counter1,pred)=stdresidpred(counter1,pred)/presented(item)
				if(stdresidpred(counter1,pred)>tolres*stdobspred(counter1,pred).and.stdobspred(counter1,pred)>0.0_8)&
                    residapred(counter1,pred)=residpred(counter1,pred)/stdresidpred(counter1,pred)
			end if
			
		
			
		end do
	end if




	counter=counter+numcatobs(item)
	
end do 

end do

return
deallocate(obsmat)
end subroutine marginpred
