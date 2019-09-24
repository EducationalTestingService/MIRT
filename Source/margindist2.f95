!	Obtain observed and fitted two-way marginal distributions of items and examine adjusted residuals.
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
!	lintran transforms the latent vector.
!	obsweight provides weights.
!	postdensity is the array of posterior densities.
!	theta is the array of quadrature points.
!	tolres is the rsidual tolerance.
!	fitmarg2 is the fitted marginal total.
!	fitpmarg2 is the  fitted marginal proportion.
!	obsmarg2 is the observed marginal total.
!	obspmarg2 is the observed marginal proportion.
!	presented2 counts weighted items presented in pairs.
!	residmarg2 is the residual for the marginal total.
!	residamarg2 is the adjusted residual.
!	residpmarg2 is the residual for the marginal proportion.
!	stdobsmarg2 is the standard deviation of the observed marginal.
!	stdobspmarg2 is the asymptotic standard deviation of the observed 
!		marginal fraction.
!	stdresidmarg2 is the asymptotic standard deviation of the residual marginal.
!	stdresidpmarg2 is the asymptotic standard deviation of the residual 
!		marginal fraction.



subroutine margindist2(catobsrange,dat,npsu,nstratum,numcat,numcatobs,&
	psu,stratum,complx,resid,stratify,usepsu,&
	beta,eacovgaminv_louis,gradsdes,lintran,&
	obsweight,postdensity,theta,tolres,&
	fitmarg2,fitpmarg2,obsmarg2,obspmarg2,presented2,&
	residamarg2,residmarg2,residpmarg2,&
	stdobsmarg2,stdobspmarg2,&
	stdresidmarg2,stdresidpmarg2)

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
real(kind=8),intent(out)::fitmarg2(:,:),fitpmarg2(:,:),obsmarg2(:,:),obspmarg2(:,:),presented2(:,:),&
	residamarg2(:,:),residmarg2(:,:),residpmarg2(:,:),&
	stdobsmarg2(:,:),stdobspmarg2(:,:),&
	stdresidmarg2(:,:),stdresidpmarg2(:,:)
!	al is an error indicator for allocation.
!	col counts columns.
!	counter, counter0, and counter1 find underlying item codes.
!	dimdesign is the design dimension.
!	dimlatout is the number of skills.
!	item is an item.
!	item1 is another item
!	ncat is number of underlying categories.
!   ncat0 is number of underlying categories for first item1 items.
!	ncatobs is number of observed categories.
!	nitems is the number of items.
!	nobs is number of observations.
!	nquad is the number of quadrature points.
!	nscale is the position in beta of slopes.
!	nscale1 is the end position in beta of slopes.
!	obs counts observations.
!	quad counts quadrature points.
!	resp is a response vector.
!	row counts rows.
!	row1 also counts rows.
integer::al,col,counter,counter0,counter1,dimdesign,dimlatout,item,item1,&
	ncat,ncat0,ncatobs,nitems,nobs,nquad,nscale,nscale1,&
	obs,quad,resp(size(dat,1)),row,row1
!	datamask is a mask for items presented.
logical::datamask(size(dat,1))
!   avediff is an average difference.
!	covobs is the covariance matrix for adjusted observed values.
!	covsum is used to estimate asymptotic variances for fitted marginals.
!   diff is a difference.
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
nobs=size(obsweight)
nquad=size(theta,2)
nscale=ncat+1
nscale1=ncat*(dimlatout+1)
if(resid.or.complx)then
	allocate(obsmat(3,nobs,ncatobs),stat=al)
	if(al/=0)stop 'Allocation failure for two-way marginals.'
	
end if
fitmarg2=0.0_8
fitpmarg2=0.0_8
locations=beta(1:ncat)

obsmarg2=0.0_8
obspmarg2=0.0_8
presented2=0.0_8
if(resid)then
	residmarg2=0.0_8
	residpmarg2=0.0_8
	residamarg2=0.0_8
end if
scales=reshape(beta(nscale:nscale1),(/dimlatout,ncat/))
stdobsmarg2=0.0_8
stdobspmarg2=0.0_8
if(resid)then
	stdresidmarg2=0.0_8
	stdresidpmarg2=0.0_8
end if
totalobs=sum(obsweight)
do obs=1,nobs
	if(obsweight(obs)<=0.0_8)cycle
	datamask=.false.
	resp=dat(:,obs)
	do item1=1,nitems
		if(0<=resp(item1).and.resp(item1)<numcatobs(item1))datamask(item1)=.true.
		if(item1>1.and.datamask(item1))then
			do item=1,item1-1
				if(datamask(item))presented2(item,item1)=presented2(item,item1)+obsweight(obs)
			end do
		end if
	end do
end do
!	Marginal totals.
counter0=numcatobs(1)
ncat0=numcat(1)
counter1=counter0+1
do item1=2,nitems
ncat0=ncat0+numcat(item1)
do row1=counter1,counter1+numcatobs(item1)-1
if(resid.or.complx)obsmat=0.0_8
do obs=1,nobs
	if(obsweight(obs)<=0.0_8) cycle
	
! Observation obs.
	resp=dat(:,obs)
	if(resp(item1)<0.or.resp(item1)>=numcatobs(item1))cycle
    
    respvec(1:counter0)=0.0_8
    
    respprobvec(1:counter0)=0.0_8
    datamask(1:item1-1)=.false.
	datamask(item1)=.true.
	counter=1
	do item=1,item1-1
		
		
!	Verify if item was presented.
		if(resp(item)>=0.and.resp(item)<numcatobs(item))then
			datamask(item)=.true.
			if(resp(item1)==row1-counter1)respvec(counter+resp(item))=1.0_8
			
			if(resid.or.complx)obsmat(3,obs,counter:counter+numcatobs(item)-1)=1.0_8
		end if
		
		counter=counter+numcatobs(item)
	end do

	if(count(datamask(1:item1-1))==0)cycle
	if(resid.or.complx)obsmat(1,obs,1:counter0)=respvec(1:counter0)
	


!	Cycle through quadrature points.
	do quad=1,nquad
!		Probability computation.
		newtheta=matmul(lintran,theta(:,quad,obs))
        probcat(1:ncat0)=&
            probvec(numcat(1:item1),datamask(1:item1),locations(1:ncat0),scales(:,1:ncat0),newtheta)
        probcatobs(1:counter0+numcatobs(item1))=&
            probvecobsall(catobsrange(:,1:item1),numcatobs(1:item1),datamask(1:item1),probcat(1:ncat0))
		
		respprobvec(1:counter0)=respprobvec(1:counter0)+postdensity(quad,obs)*probcatobs(1:counter0)*probcatobs(row1)
	end do
!	Set up obsmat and observed and fitted marginal totals.
    obsmarg2(1:counter0,row1)=obsmarg2(1:counter0,row1)+obsweight(obs)*respvec(1:counter0)
	
    fitmarg2(1:counter0,row1)=fitmarg2(1:counter0,row1)+obsweight(obs)*respprobvec(1:counter0)
!	Get standard deviation of indicator.
	
	
	if(resid)obsmat(2,obs,1:counter0)=respvec(1:counter0)-respprobvec(1:counter0)
	
	
end do
!	The residual for marginal totals.
if(resid) residmarg2(1:counter0,row1)=obsmarg2(1:counter0,row1)-fitmarg2(1:counter0,row1)
!
!	Fractions and variances.
counter=1

do item=1,item1-1
	if(presented2(item,item1)>0.0_8) then
		do row=counter,counter+numcatobs(item)-1
			fitpmarg2(row,row1)= fitmarg2(row,row1)/presented2(item,item1)
			obspmarg2(row,row1)= obsmarg2(row,row1)/presented2(item,item1)
			if(resid)residpmarg2(row,row1)= residmarg2(row,row1)/presented2(item,item1)
			
				

			
	
!	Standard deviation of count depends on sampling.
			if(.not.complx)then
				stdobsmarg2(row,row1)=obsmarg2(row,row1)*(1.0_8-obsmarg2(row,row1)/totalobs)
				stdobsmarg2(row,row1)=sqrt(stdobsmarg2(row,row1))
				stdobspmarg2(row,row1)=sqrt(obspmarg2(row,row1)*(1.0_8-obspmarg2(row,row1))/presented2(item,item1))
				if(resid)then
					covsum=reshape(crossinformation(.true.,obsmat(2:2,:,row),gradsdes,obsweight),(/dimdesign/))
					slopes=matmul(eacovgaminv_louis,covsum)
                    avediff=0.0_8
                    do obs=1,nobs
                        if(obsweight(obs)<=0.0_8) cycle
                        obsmat(2,obs,row)=obsmat(2,obs,row)-dot_product(gradsdes(:,obs),slopes)
                        avediff=avediff+obsmat(2,obs,row)*obsweight(obs)
                    end do
                    avediff=avediff/totalobs
                    do obs=1,nobs
                        diff=obsmat(2,obs,row)-avediff
                        stdresidmarg2(row,row1)=stdresidmarg2(row,row1)+diff*diff*obsweight(obs)
                    end do
					
					
					
				end if
			else
				covobs=information_complex(npsu,nstratum,psu,stratum,stratify,usepsu,&
					obsmat(1:1,:,row),obsweight)
					
				stdobsmarg2(row,row1)=sqrt(covobs(1,1))
				
				obsmat(1,:,row)=obsmat(1,:,row)-obspmarg2(row,row1)*obsmat(3,:,row)
				covobs=information_complex(npsu,nstratum,psu,stratum,stratify,usepsu,&
					obsmat(1:1,:,row),obsweight)
				stdobspmarg2(row,row1)=sqrt(covobs(1,1)/presented2(item,item1))
				if(resid)then
					covsum=reshape(crossinformation(.true.,obsmat(2:2,:,item),gradsdes,obsweight),(/dimdesign/))
					slopes=matmul(eacovgaminv_louis,covsum)
                
					do obs=1,nobs
						if(obsweight(obs)<=0.0_8) cycle
						obsmat(2,obs,item)=obsmat(2,obs,item)-dot_product(gradsdes(:,obs),slopes)
                     
					end do
					covobs=information_complex(npsu,nstratum,psu,stratum,stratify,usepsu,&
						obsmat(2:2,:,item),obsweight)
				
					stdresidmarg2(row,row1)=covobs(1,1)
				end if
			end if
			if(resid)then
				stdresidmarg2(row,row1)=sqrt(stdresidmarg2(row,row1))
				stdresidpmarg2(row,row1)=stdresidmarg2(row,row1)/presented2(item,item1)
				if(stdresidmarg2(row,row1)>tolres*stdobsmarg2(row,row1).and.stdobsmarg2(row,row1)>0.0_8)&
                    residamarg2(row,row1)=residmarg2(row,row1)/stdresidmarg2(row,row1)
			end if
			
		
			
		end do
	end if




	counter=counter+numcatobs(item)

end do 
end do
counter0=counter0+numcatobs(item1)
counter1=counter0+1


end do
if(resid.or.complx)deallocate(obsmat)
return
end subroutine margindist2
