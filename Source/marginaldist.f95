!	Obtain observed and fitted marginal distributions of items and examine adjusted residuals.
!	catobsrange maps underlying to observed categories.
!	dat provides responses.
!	npsu is the number of psus's per stratum.
!	nstratum is the number of strata.
!	numcat provides the underlying number of categories per item.
!	numcatobs provides the observed number of categories per item.
!	psu indicates the psu of an observation.
!	stratum indicates the stratum of an observation.
!	complx indicates complex sampling.
!	resid is for residuals.
!	stratify is .true. for stratified sampling.
!	usepsu is .true. for primary sampling units within strata.
!	beta provides parameter estimates.


!	eacovgaminv_louis is the inverse of the Louis information plus the
!		constraint component of the negative hessian.
!	gradsdes provides gradients for observations.
!	lintran transforms the latent vector.!	obsweight provides weights.
!	postdensity is the array of posterior densities.
!	theta is the array of quadrature points.
!	tolres is the rsidual tolerance.
!	fitmarg is the fitted marginal total.
!	fitpmarg is the  fitted marginal proportion.
!	obsmarg is the observed marginal total.
!	obspmarg is the observed marginal proportion.
!	presented counts weighted items presented.
!	relmarg is the reliability of the item category indicator.
!	residmarg is the residual for the marginal total.
!	residamarg is the adjusted residual.
!	residpmarg is the residual for the marginal proportion.
!	seobsmarg is the standard error for the observed indicator.
!	stdobsmarg is the standard deviation of the observed marginal.
!	stdobspmarg is the asymptotic standard deviation of the observed 
!		marginal fraction.
!	stdresidmarg is the asymptotic standard deviation of the residual marginal.
!	stdresidpmarg is the asymptotic standard deviation of the residual 
!		marginal fraction.
!	stobsmarg is the standard deviation of the indicator.


subroutine marginaldist(catobsrange,dat,npsu,nstratum,numcat,numcatobs,&
	psu,stratum,complx,resid,stratify,usepsu,&
	beta,eacovgaminv_louis,gradsdes,lintran,&
	obsweight,postdensity,theta,tolres,&
	fitmarg,fitpmarg,obsmarg,obspmarg,presented,&
	relmarg,residamarg,residmarg,residpmarg,&
	seobsmarg,stdobsmarg,stdobspmarg,&
	stdresidmarg,stdresidpmarg,stobsmarg)

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
integer,intent(in)::catobsrange(:,:),dat(:,:),npsu(:),nstratum,numcat(:),numcatobs(:),&
	psu(:),stratum(:)
logical,intent(in)::complx,resid,stratify,usepsu
real(kind=8),intent(in)::beta(:),eacovgaminv_louis(:,:),&
	gradsdes(:,:),lintran(:,:),obsweight(:),&
	postdensity(:,:),theta(:,:,:),tolres
real(kind=8),intent(out)::fitmarg(:),fitpmarg(:),obsmarg(:),obspmarg(:),presented(:),&
	relmarg(:),residamarg(:),residmarg(:),residpmarg(:),&
	seobsmarg(:),stdobsmarg(:),stdobspmarg(:),&
	stdresidmarg(:),stdresidpmarg(:),stobsmarg(:)
!	al is allocation flag.
!	col counts columns.
!	counter finds underlying item codes.
!	dimdesign is the design dimension.
!	dimlatout is the number of skills.
!	item is an item.
!	ncat is number of underlying categories.
!	ncatobs is number of observed categories.
!	nitems is the number of items.
!	nlin is the position in beta of linear terms.
!	nobs is number of observations.
!	nquad is the number of quadrature points.
!	nscale is the position in beta of slopes.
!	obs counts observations.
!	quad counts quadrature points.
!	resp is a response vector.
!	row counts rows.
integer::al,col,counter,dimdesign,dimlatout,item,ncat,ncatobs,nitems,nlin,nobs,nquad,nscale,&
	obs,quad,resp(size(dat,1)),row
!	datamask is a mask for items presented.
logical::datamask(size(dat,1))
!   avediff is an average difference.
!	covobs is the covariance matrix for adjusted observed values.
!	covsum is used to estimate asymptotic variances for fitted marginals.

!	locations gives location parameters.
!	newtheta is the transformed latent vector.
!	obsmat is the array of response indicators and fits.
!		The first row is the observed indicator, the second is the fit,
!		the third is 1 if presented and 0 otherwise.
!	probcat is the vector of underlying probabilities.
!	probcatobs is the vector of observed probabilities.
!	respprobvec is the vector of fitted probabilites for an observation.
!	respresid is the residual component for an observation.
!	respvec is the corresponding vector of observations.
!	scales gives scale parameters.
!	slopes is for regression of errors on gradients. 
!	totalobs is the sum of the observation weights. 

real(kind=8)::avediff,covobs(1,1),covsum(size(gradsdes,1)),diff,&
	locations(sum(numcat)),newtheta(size(lintran,1)),&
	probcat(sum(numcat)),probcatobs(sum(numcatobs)),&
	respprobvec(sum(numcatobs)),respresid,respvec(sum(numcatobs)),&
	scales(size(lintran,1),sum(numcat)),slopes(size(gradsdes,1)),totalobs
real(kind=8),allocatable::obsmat(:,:,:)
!	Set up parameter arrays for simplified processing.
dimdesign=size(gradsdes,1)

dimlatout=size(lintran,1)
ncat=sum(numcat)
ncatobs=sum(numcatobs)
nitems=size(dat,1)
nlin=ncat*(dimlatout+1)+1
nobs=size(obsweight)
nquad=size(theta,2)
nscale=ncat+1
if(resid.or.complx)then
	allocate(obsmat(3,nobs,ncatobs),stat=al)
	if(al/=0)stop 'Allocation of array failed for marginal distributions.'
	obsmat=0.0_8
endif
fitmarg=0.0_8
fitpmarg=0.0_8
locations=beta(1:ncat)
obsmarg=0.0_8
obspmarg=0.0_8
presented=0.0_8
relmarg=0.0_8
if(resid)then
	residmarg=0.0_8
	residpmarg=0.0_8
	residamarg=0.0_8
end if

scales=reshape(beta(nscale:nlin-1),(/dimlatout,ncat/))
seobsmarg=0.0_8
stobsmarg=0.0_8
stdobsmarg=0.0_8
stdobspmarg=0.0_8
if(resid)then
	stdresidmarg=0.0_8
	stdresidpmarg=0.0_8
end if
totalobs=sum(obsweight)

!	Marginal totals.

do obs=1,nobs
	if(obsweight(obs)<=0.0_8) cycle
! Observation obs.
	resp=dat(:,obs)
	respvec=0.0_8
	respprobvec=0.0_8
	datamask=.false.
	counter=1
	do item=1,nitems
!	Verify if item was presented.
		if(resp(item)>=0.and.resp(item)<numcatobs(item))then
			datamask(item)=.true.
			respvec(counter+resp(item))=1.0_8
			presented(item)=presented(item)+obsweight(obs)
			if(resid.or.complx)obsmat(3,obs,counter:counter+numcatobs(item)-1)=1.0_8
		end if
		counter=counter+numcatobs(item)
	end do

	if(count(datamask)==0)cycle
	if(resid.or.complx)obsmat(1,obs,:)=respvec
	


!	Cycle through quadrature points.
	do quad=1,nquad
!		Probability computation.
		newtheta=matmul(lintran,theta(:,quad,obs))
		probcat=probvec(numcat,datamask,locations,scales,newtheta)
		probcatobs=probvecobsall(catobsrange,numcatobs,datamask,probcat)
		respprobvec=respprobvec+postdensity(quad,obs)*probcatobs
	end do
!	Set up obsmat and observed and fitted marginal totals.
	
	obsmarg=obsmarg+obsweight(obs)*respvec
	fitmarg=fitmarg+obsweight(obs)*respprobvec
!	Get standard error of measurement of indicator
	
	do row=1,ncatobs
		respresid=respvec(row)-respprobvec(row)
		if(respresid/=0.0_8)seobsmarg(row)=seobsmarg(row)+obsweight(obs)*respresid*respresid
		if(resid)obsmat(2,obs,row)=respresid
	end do	
	

end do

!	The residual for marginal totals.
if(resid) residmarg=obsmarg-fitmarg
!
!	Fractions and variances.
counter=1

do item=1,nitems
	if(presented(item)>0.0_8) then
    
		do row=counter,counter+numcatobs(item)-1
			fitpmarg(row)= fitmarg(row)/presented(item)
			obspmarg(row)= obsmarg(row)/presented(item)
			if(resid)then
				residpmarg(row)= residmarg(row)/presented(item)
				
			end if
				
!	Final standard error of measurement.
			seobsmarg(row)=seobsmarg(row)/presented(item)
			stobsmarg(row) = obspmarg(row)&
				*(1.0_8-obspmarg(row))
!	Reliability.
			if(stobsmarg(row)>0.0_8)relmarg(row)=1.0_8-seobsmarg(row)/stobsmarg(row)
			seobsmarg(row)=sqrt(seobsmarg(row))
			stobsmarg(row)=sqrt(stobsmarg(row))
	
!	Standard deviation of count depends on sampling.
			if(.not.complx)then
				stdobsmarg(row)=obsmarg(row)*(1.0_8-obsmarg(row)/totalobs)
				stdobsmarg(row)=sqrt(stdobsmarg(row))
				stdobspmarg(row)=sqrt(obspmarg(row)*(1.0_8-obspmarg(row))/presented(item))
				if(resid)then
					covsum=reshape(crossinformation(.true.,obsmat(2:2,:,row),gradsdes,obsweight),(/dimdesign/))
					slopes=matmul(eacovgaminv_louis,covsum)
                    avediff=0.0_8
                    do obs=1,nobs
                        if(obsweight(obs)<=0.0_8)cycle
                        obsmat(2,obs,row)=obsmat(2,obs,row)-dot_product(gradsdes(:,obs),slopes)
                        avediff=avediff+obsmat(2,obs,row)*obsweight(obs)
                    end do
                    avediff=avediff/totalobs
                    do obs=1,nobs
                        diff=obsmat(2,obs,row)-avediff
                        stdresidmarg(row)=stdresidmarg(row)+diff*diff*obsweight(obs)
					end do
					
					
				end if
					
				
			else
				
				covobs=information_complex(npsu,nstratum,psu,stratum,stratify,usepsu,&
					obsmat(1:1,:,row),obsweight)
				stdobsmarg(row)=sqrt(covobs(1,1))
				
				obsmat(1,:,row)=obsmat(1,:,row)-obspmarg(row)*obsmat(3,:,row)
				covobs=information_complex(npsu,nstratum,psu,stratum,stratify,usepsu,&
					obsmat(1:1,:,row),obsweight)
				stdobspmarg(row)=sqrt(covobs(1,1)/presented(item))
				if(resid)then
					covsum=reshape(crossinformation(.true.,obsmat(2:2,:,row),gradsdes,obsweight),(/dimdesign/))
					slopes=matmul(eacovgaminv_louis,covsum)
					
                    do obs=1,nobs
                        obsmat(2,obs,row)=obsmat(2,obs,row)-dot_product(gradsdes(:,obs),slopes)
                    end do
                    covobs=information_complex(npsu,nstratum,psu,stratum,stratify,usepsu,&
                        obsmat(2:2,:,row),obsweight)
					stdresidmarg(row)=covobs(1,1)
				end if
				
			end if
			if(resid)then
				stdresidmarg(row)=sqrt(stdresidmarg(row))
				stdresidpmarg(row)=stdresidmarg(row)/presented(item)
				if(stdresidmarg(row)>tolres*stdobsmarg(row).and.stdobsmarg(row)>0.0_8)residamarg(row)=residmarg(row)/stdresidmarg(row)
			end if
			
		
			
		end do
	end if




	counter=counter+numcatobs(item)

end do 
if(resid.or.complx)deallocate(obsmat)
return
end subroutine marginaldist
