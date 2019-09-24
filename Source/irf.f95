!	Obtain observed and fitted item response functions and examine adjusted residuals.
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
!	indvar is array of predictors.
!	lintran transforms the latent vector.
!	obsweight provides weights.
!	postdensity is the array of posterior densities.
!   prob is the arrayof marginal probabilities.
!	theta is the array of quadrature points.
!   thetacheck is the array of points to check.
!	tolres is the rsidual tolerance.
!	fitirf is the unconditional estimate of the tem response function.
!	obsirf is the conditional estimate of the item response function.
!	presentedirf counts weighted items presented.
!	residirf is the residual for the item response function.
!	residairfarg is the adjusted residual.
!	stdresidirf is the asymptotic standard deviation of the residual item response function.



subroutine irf(catobsrange,dat,npsu,nstratum,numcat,numcatobs,&
	psu,stratum,complx,resid,stratify,usepsu,&
	beta,eacovgaminv_louis,gradsdes,indvar,lintran,&
	obsweight,postdensity,prob,theta,thetacheck,tolres,&
	fitirf,obsirf,presentedirf,&
	residairf,residirf,&
	stdresidirf)

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


!	convert an n by n symmetric matrix from compact to regular form.
	function expandmat(n,compmat)
		implicit none
		integer,intent(in)::n
		real(kind=8),intent(in)::compmat(:)
		real(kind=8)::expandmat(n,n)
	end function expandmat
	
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


!	probvecobs is used to compute observed item probabilities given the latent vector.
!	catobsrange is the table of ranges of underlying categories per observed category.
!	numcatobs provides the number of observed categories per item.
!	resp is the response vector.
!	mask indicates which items were presented.
!	probcat is the vector of underlying item probabilities.
	function probvecobs(catobsrange,numcatobs,resp,mask,probcat)
		implicit none
		integer,intent(in)::catobsrange(:,:),numcatobs(:),resp(:)
		logical,intent(in)::mask(:)
		real(kind=8),intent(in)::probcat(:)
		real(kind=8)::probvecobs(size(mask))
	end function probvecobs



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
	gradsdes(:,:),indvar(:,:),lintran(:,:),obsweight(:),&
	postdensity(:,:),prob(:),theta(:,:,:),thetacheck(:,:),tolres
real(kind=8),intent(out)::fitirf(:,:),obsirf(:,:),presentedirf(:),&
	residairf(:,:),residirf(:,:),&
	stdresidirf(:,:)
!	al is an allocation flag.
!	col counts columns.
!	counter finds underlying item codes.
!	dimdesign is the design dimension.
!	dimlatin is the number of factors.
!	dimlatout is the number of skills.
!	item is an item.
!	ncat is number of underlying categories.
!	ncatobs is number of observed categories.
!	nitems is the number of items.
!	nlin is the position in beta of linear terms.
!	nobs is number of observations.
!	npred is number of predictors.
!	nquad is the number of quadrature points.
!	nscale is the position in beta of slopes.
!   nthetacheck is the number of points to check.
!	obs counts observations.
!   point counts points to check.
!	pred counts predictors.
!	quad counts quadrature points.
!	resp is a response vector.
!	row counts rows.
integer::al,col,counter,dimdesign,dimlatin,dimlatout,item,ncat,ncatobs,nitems,nlin,&
	nobs,npred,nquad,nscale,&
	nthetacheck,obs,point,pred,quad,resp(size(dat,1)),row
!	amask is a mask for items presented.
logical::mask(size(dat,1))
!   avediff is an average difference.
!	covobs is the covariance matrix for adjusted observed values.
!	covsum is used to estimate asymptotic variances for fitted marginals.
!   diff is a difference.
!	locations gives location parameters.
!	newtheta is the transformed latent vector.
!	obsmat is the array of residuals for observations.
!	pratio is the ratio of conditional and unconditional probabilities.
!	probcat is the vector of underlying probabilities.
!	probcatobs is the vector of observed probabilities.
!	respprobvec is the vector of fitted probabilites for an observation.
!	respresid is the residual component for an observation.
!	respvec is the corresponding vector of observations.
!	scales gives scale parameters.
!	slopes is for regression of errors on gradients.
!   sumirf is for total latent weights.
!	totalobs is the sum of the observation weights. 

real(kind=8)::avediff,covobs(1,1),covsum(size(gradsdes,1)),diff,&
	locations(sum(numcat)),newtheta(size(lintran,1)),&
	pratio,probcat(sum(numcat)),&
	probcatobs(size(dat,1)),&
	respprobvec(sum(numcatobs)),respresid,respvec(sum(numcatobs)),&
	scales(size(lintran,1),sum(numcat)),slopes(size(gradsdes,1)),&
	sumirf(size(dat,1)),totalobs
real(kind=8),allocatable::obsmat(:,:,:)
!	Set up parameter arrays for simplified processing.
dimdesign=size(gradsdes,1)
dimlatin=size(theta,1)
dimlatout=size(lintran,1)
ncat=sum(numcat)
ncatobs=sum(numcatobs)
nitems=size(dat,1)
nlin=ncat*(dimlatout+1)+1
nobs=size(obsweight)
npred=size(indvar,1)
nquad=size(theta,2)
nscale=ncat+1
nthetacheck=size(thetacheck,2)
fitirf=0.0_8
locations=beta(1:ncat)
obsirf=0.0_8


presentedirf=0.0_8

if(resid)then
	allocate(obsmat(1,nobs,ncatobs),stat=al)
	if(al/=0) stop 'Allocation of array failed for residuals of item response functions.'
	obsmat=0.0_8
	residirf=0.0_8
	
	residairf=0.0_8
end if
scales=reshape(beta(nscale:nlin-1),(/dimlatout,ncat/))


if(resid)stdresidirf=0.0_8
sumirf=0.0_8
totalobs=sum(obsweight)
!	Item response functions.
do point=1,nthetacheck
newtheta=matmul(lintran,thetacheck(:,point))
do obs=1,size(obsweight)
	if(obsweight(obs)<=0.0_8) cycle
! Observation obs.
	resp=dat(:,obs)
	respvec=0.0_8
	respprobvec=0.0_8
	mask=.false.
	counter=1
	do item=1,nitems
!	Verify if item was presented.
		if(resp(item)>=0.and.resp(item)<numcatobs(item))then
			mask(item)=.true.
			respvec(counter+resp(item))=1.0_8
			if(point==1)presentedirf(item)=presentedirf(item)+obsweight(obs)
			
		end if
		counter=counter+numcatobs(item)
	end do

	if(count(mask)==0)cycle
	



			

!	Probability computation.

	
	probcat=probvec(numcat,mask,locations,scales,newtheta)

	probcatobs=probvecobs(catobsrange,numcatobs,resp,mask,probcat)

	pratio=product(probcatobs,mask)/prob(obs)
	respprobvec=probvecobsall(catobsrange,numcatobs,mask,probcat)
	
!	Set up obsmat and observed and fitted marginal totals.
	
	obsirf(:,point)=obsirf(:,point)+obsweight(obs)*pratio*respvec
	fitirf(:,point)=fitirf(:,point)+obsweight(obs)*pratio*respprobvec
!	Get residual variance and residual if needed.
	if(resid) obsmat(1,obs,:)=pratio*(respvec-respprobvec)
	
	

end do
!	Adjust numbers.
counter=1
do item=1,nitems
	if(presentedirf(item)>0.0_8)then
		sumirf(item)=sum(obsirf(counter:counter+numcatobs(item)-1,point))
		obsirf(counter:counter+numcatobs(item)-1,point)=&
			obsirf(counter:counter+numcatobs(item)-1,point)/sumirf(item)
		fitirf(counter:counter+numcatobs(item)-1,point)=&
			fitirf(counter:counter+numcatobs(item)-1,point)/sumirf(item)
		if(resid)obsmat(1,:,counter:counter+numcatobs(item)-1)=&
			obsmat(1,:,counter:counter+numcatobs(item)-1)/sumirf(item)
	end if
	counter=counter+numcatobs(item)
end do
if(resid)then
	residirf(:,point)=obsirf(:,point)-fitirf(:,point)
	do row=1,ncatobs
!	Standard deviations depends on sampling.
		if(.not.complx)then
			covsum=reshape(crossinformation(.true.,obsmat(:,:,row),gradsdes,obsweight),(/dimdesign/))
			slopes=matmul(eacovgaminv_louis,covsum)
			avediff=0.0_8
			do obs=1,size(obsweight)
				if(obsweight(obs)<=0.0_8)cycle
				obsmat(1,obs,row)=obsmat(1,obs,row)-dot_product(gradsdes(:,obs),slopes)
				avediff=avediff+obsmat(1,obs,row)*obsweight(obs)
			end do
			avediff=avediff/totalobs
			do obs=1,size(obsweight)
				diff=obsmat(1,obs,row)-avediff
				stdresidirf(row,point)=stdresidirf(row,point)+diff*diff*obsweight(obs)
			end do
					
					
		else
			
			covsum=reshape(crossinformation(.true.,obsmat(:,:,row),gradsdes,obsweight),(/dimdesign/))
			slopes=matmul(eacovgaminv_louis,covsum)
            do obs=1,size(obsweight)
				obsmat(1,obs,row)=obsmat(1,obs,row)-dot_product(gradsdes(:,obs),slopes)
            end do
            covobs=information_complex(npsu,nstratum,psu,stratum,stratify,usepsu,&
                obsmat(1:1,:,row),obsweight)
			stdresidirf(row,point)=covobs(1,1)
		end if
		stdresidirf(row,point)=sqrt(stdresidirf(row,point))
				
		if(stdresidirf(row,point)>0.0_8)residairf(row,point)&
			=residirf(row,point)/stdresidirf(row,point)
	end do
end if			
		
			
end do
if(resid)deallocate(obsmat)


return
end subroutine irf