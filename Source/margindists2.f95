!	Obtain observed and fitted products of item scores and examine adjusted residuals.
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
!	scores provides item scores.
!	theta is the array of quadrature points.
!	tolres is the rsidual tolerance.
!	fitmargs2 is the fitted marginal total.
!	fitpmargs2 is the  fitted marginal proportion.
!	obsmargs2 is the observed marginal total.
!	obspmargs2 is the observed marginal proportion.
!	presenteds2 counts weighted items presented in pairs.
!	residmargs2 is the residual for the marginal total.
!	residamargs2 is the adjusted residual.
!	residpmargs2 is the residual for the marginal proportion.
!	stdobsmargs2 is the standard deviation of the observed marginal.
!	stdobspmargs2 is the asymptotic standard deviation of the observed 
!		marginal fraction.
!	stdresidmargs2 is the asymptotic standard deviation of the residual marginal.
!	stdresidpmargs2 is the asymptotic standard deviation of the residual 
!		marginal fraction.



subroutine margindists2(catobsrange,dat,npsu,nstratum,numcat,numcatobs,&
	psu,stratum,complx,resid,stratify,usepsu,&
	beta,eacovgaminv_louis,gradsdes,lintran,&
	obsweight,postdensity,scores,theta,tolres,&
	fitmargs2,fitpmargs2,obsmargs2,obspmargs2,presenteds2,&
	residamargs2,residmargs2,residpmargs2,&
	stdobsmargs2,stdobspmargs2,&
	stdresidmargs2,stdresidpmargs2)

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
	postdensity(:,:),scores(:),theta(:,:,:),tolres
real(kind=8),intent(out)::fitmargs2(:,:),fitpmargs2(:,:),obsmargs2(:,:),obspmargs2(:,:),presenteds2(:,:),&
	residamargs2(:,:),residmargs2(:,:),residpmargs2(:,:),&
	stdobsmargs2(:,:),stdobspmargs2(:,:),&
	stdresidmargs2(:,:),stdresidpmargs2(:,:)
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
!	obs counts observations.
!	quad counts quadrature points.
!	resp is a response vector.
!	row, row1, and row2 count rows.

integer::al,col,counter,counter0,counter1,dimdesign,dimlatout,item,item1,&
	ncat,ncat0,ncatobs,nitems,nobs,nquad,nscale,nscale1,&
	obs,quad,resp(size(dat,1)),row,row1,row2
!	datamask is a mask for items presented.
logical::datamask(size(dat,1))
!   avediff is an average difference.
!	covobs is the covariance matrix for adjusted observed values.
!	covsum is used to estimate asymptotic variances for fitted marginals.
!   diff is a difference.
!	es is an expected score.
!	esc is a conditional expected score.
!	locations gives location parameters.
!	newtheta is the transformed latent vector.
!	obsmat is the array of response indicators and fits.
!		The first row is the observed indicator, the second is the fit,
!		the third is 1 if presented and 0 otherwise.
!	probcat is the vector of underlying probabilities.
!	probcatobs is the vector of observed probabilities.
!	respresid is the residual product component for an observation.
!	respvec is the corresponding vector of products observations.
!	scales gives scale parameters.
!	slopes is for regression of errors on gradients. 
!	totalobs is the sum of the observation weights. 

real(kind=8)::avediff,covobs(1,1),covsum(size(gradsdes,1)),diff,es(size(numcat)),&
	esc(size(numcat)),locations(sum(numcat)),newtheta(size(lintran,1)),&
	probcat(sum(numcat)),probcatobs(sum(numcatobs)),&
	respresid,respvec(size(numcatobs)),&
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
allocate(obsmat(3,nobs,nitems),stat=al)
if(al/=0) stop 'Allocation failed for marginal two-way scores.'

fitmargs2=0.0_8
fitpmargs2=0.0_8
locations=beta(1:ncat)

obsmargs2=0.0_8
obspmargs2=0.0_8
presenteds2=0.0_8
if(resid)then
	residmargs2=0.0_8
	residpmargs2=0.0_8
	residamargs2=0.0_8
end if

scales=reshape(beta(nscale:nscale1),(/dimlatout,ncat/))
stdobsmargs2=0.0_8
stdobspmargs2=0.0_8
if(resid)then
	stdresidmargs2=0.0_8
	stdresidpmargs2=0.0_8
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
				if(datamask(item))presenteds2(item,item1)=presenteds2(item,item1)+obsweight(obs)
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
obsmat=0.0_8
do obs=1,nobs
	if(obsweight(obs)<=0.0_8) cycle
	
! Observation obs.
	resp=dat(:,obs)
	if(resp(item1)<0.or.resp(item1)>=numcatobs(item1))cycle
    
    respvec(1:item1-1)=0.0_8
    
    es(1:item1-1)=0.0_8
    datamask(1:item1-1)=.false.
	datamask(item1)=.true.
	counter=1
	respvec(item1)=scores(counter1+resp(item1))
	do item=1,item1-1
		
		
!	Verify if item was presented.
		if(resp(item)>=0.and.resp(item)<numcatobs(item))then
			datamask(item)=.true.
			respvec(item)=scores(counter+resp(item))*respvec(item1)
			
			obsmat(3,obs,item)=1.0_8
		end if
		
		counter=counter+numcatobs(item)
	end do

	if(count(datamask)==0)cycle
	obsmat(1,obs,1:item1-1)=respvec(1:item1)
	


!	Cycle through quadrature points.
	do quad=1,nquad
!		Probability computation.
		newtheta=matmul(lintran,theta(:,quad,obs))
        probcat(1:ncat0)=&
            probvec(numcat(1:item1),datamask(1:item1),locations(1:ncat0),scales(:,1:ncat0),newtheta)
        probcatobs(1:counter0+numcatobs(item1))=&
            probvecobsall(catobsrange(:,1:item1),numcatobs(1:item1),datamask(1:item1),probcat(1:ncat0))
		row1=1
		row2=0
		do row=1,item1
            row2=row2+numcatobs(row)
			esc(row)=sum(probcatobs(row1:row2)*scores(row1:row2))
			row1=row2+1

		end do

		
		es(1:item1-1)=es(1:item1-1)+postdensity(quad,obs)*esc(1:item1-1)*esc(item1)
	end do
!	Set up obsmat and observed and fitted marginal totals.
    obsmargs2(1:item1-1,item1)=obsmargs2(1:item1-1,item1)+obsweight(obs)*respvec(1:item1-1)
	
    fitmargs2(1:item1-1,item1)=fitmargs2(1:item1-1,item1)+obsweight(obs)*es(1:item1-1)
!	Get standard deviation of indicator.
	
	
	if(resid)obsmat(2,obs,1:item1-1)=respvec(1:item1-1)-es(1:item1-1)
	
	
end do

!	The residual for marginal totals.
if(resid) residmargs2(1:item1-1,item1)=obsmargs2(1:item1-1,item1)-fitmargs2(1:item1-1,item1)
!
!	Fractions and variances.


do item=1,item1-1
	if(presenteds2(item,item1)>0.0_8) then
		fitpmargs2(item,item1)= fitmargs2(item,item1)/presenteds2(item,item1)
		obspmargs2(item,item1)= obsmargs2(item,item1)/presenteds2(item,item1)
		if(resid)residpmargs2(item,item1)= residmargs2(item,item1)/presenteds2(item,item1)
			
				

			
	
!	Standard deviation of count depends on sampling.
		if(.not.complx)then
			stdobsmargs2(item,item1)=sum((obsmat(1,:,item)-obspmargs2(item,item1))**2*obsweight)
			obsmat(1,:,item)=obsmat(1,:,item)-obspmargs2(item,item1)*obsmat(3,:,item)	
			stdobspmargs2(item,item1)=sum(obsmat(1,:,item)**2*obsweight)
			stdobsmargs2(item,item1)=sqrt(stdobsmargs2(item,item1))
			stdobspmargs2(item,item1)=sqrt(stdobsmargs2(item,item1))
			if(resid)then
				covsum=reshape(crossinformation(.true.,obsmat(2:2,:,item),gradsdes,obsweight),(/dimdesign/))
				slopes=matmul(eacovgaminv_louis,covsum)
                avediff=0.0_8
                do obs=1,nobs
                    if(obsweight(obs)<=0.0_8) cycle
                     obsmat(2,obs,item)=obsmat(2,obs,item)-dot_product(gradsdes(:,obs),slopes)
                     avediff=avediff+obsmat(2,obs,item)*obsweight(obs)
               end do
               avediff=avediff/totalobs
               do obs=1,nobs
				if(obsweight(obs)<=0.0_8) cycle
                   diff=obsmat(2,obs,item)-avediff
                   stdresidmargs2(item,item1)=stdresidmargs2(item,item1)+diff*diff*obsweight(obs)
               end do
					
			end if
			
		else
			covobs=information_complex(npsu,nstratum,psu,stratum,stratify,usepsu,&
				obsmat(1:1,:,item),obsweight)
					
			stdobsmargs2(item,item1)=sqrt(covobs(1,1))
				
			obsmat(1,:,item)=obsmat(1,:,item)-obspmargs2(item,item1)*obsmat(3,:,item)
			covobs=information_complex(npsu,nstratum,psu,stratum,stratify,usepsu,&
				obsmat(1:1,:,item),obsweight)
			stdobspmargs2(item,item1)=sqrt(covobs(1,1)/presenteds2(item,item1))
			if(resid)then
				covsum=reshape(crossinformation(.true.,obsmat(2:2,:,item),gradsdes,obsweight),(/dimdesign/))
				slopes=matmul(eacovgaminv_louis,covsum)
                
                do obs=1,nobs
                    if(obsweight(obs)<=0.0_8) cycle
                     obsmat(2,obs,item)=obsmat(2,obs,item)-dot_product(gradsdes(:,obs),slopes)
                     
               end do
				covobs=information_complex(npsu,nstratum,psu,stratum,stratify,usepsu,&
                    obsmat(2:2,:,item),obsweight)
				stdresidmargs2(item,item1)=covobs(1,1)
			end if
		end if
		if(resid)then
			stdresidmargs2(item,item1)=sqrt(stdresidmargs2(item,item1))
			stdresidpmargs2(item,item1)=stdresidmargs2(item,item1)/presenteds2(item,item1)
			if(stdresidmargs2(item,item1)>tolres*stdobsmargs2(item,item1).and.stdobsmargs2(item,item1)>0.0_8)&
				residamargs2(item,item1)=residmargs2(item,item1)/stdresidmargs2(item,item1)
		end if
	end if




end do 
counter0=counter0+numcatobs(item1)
counter1=counter0+1


end do
if(resid.or.complx)deallocate(obsmat)
return
end subroutine margindists2
