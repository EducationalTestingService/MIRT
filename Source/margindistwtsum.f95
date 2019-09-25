!	Obtain marginal distribution of a weighted sum.
!	catobsrange maps underlying to observed categories.
!	dat provides responses.
!	maxscore is the maximum score.
!	maxw provides maximum item weights.
!	minscore is the minimum score.
!	minw provides minimum item weights.
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
!	weight is the weights.
!	fitcummargwtsum is the fitted cumulative marginal
!		frequency distribution.
!	fitmargwtsum is the fitted marginal frequency distribution.
!	fitpcummargwtsum is the fitted cumulative marginal
!		probability distribution.
!	fitpmargwtsum is the fitted marginal probability distribution.
!	obscummargwtsum is the observed cumulative marginal
!		frequency distribution.
!	obsmargwtsum is the observed marginal frequency distribution.
!	obspcummargwtsum is the observed cumulative marginal
!		probability distribution.
!	obspmargwtsum is the observed marginal
!		probability distribution.
!	presentedwtsum counts weighted items presented
!		for the weighted sum.
!	residacummargwtsum is the adjusted residual for the cumulative
!		marginal distribution.
!	residamargwtsum is the adjusted residual for the marginal
!		distribution.
!	residcummargwtsum is the residual for the cumulative marginal
!		frequency distribution.
!	residmargwtsum is the residual for the marginal frequency
!		distribution.
!	residpcummargwtsum is the residual for the cumulative marginal
!		 probability distribution.
!	residpmargwtsum is the residual for the marginal probability
!		distribution.
!	stdobscummargwtsum is the standard deviation of the observed
!		cumulative marginal frequency distribution.
!	stdobsmargwtsum is the standard deviation of the observed
!		marginal frequency distribution.
!	stdobspcummargwtsum is the asymptotic standard deviation of
!		the observed cumulative marginal probability distribution.
!	stdobspmargwtsum is the asymptotic standard deviation of the
!		observed marginal probability distribution.
!	stdresidcummargwtsum is the asymptotic standard deviation
!		of the residual cumulative marginal frequency
!		distribution.
!	stdresidmargwtsum is the asymptotic standard deviation of the
!		residual marginal frequency distribution.
!	stdresidpcummargwtsum is the asymptotic standard deviation of
!		the residual cumulative marginal probability distribution.
!	stdresidpmargwtsum is the asymptotic standard deviation of the
!		residual marginal probability distribution.

!
subroutine margindistwtsum(catobsrange,dat,maxscore,maxw,minscore,minw,&
	npsu,nstratum,&
	numcat,numcatobs,&
	psu,stratum,complx,resid,stratify,usepsu,&
	beta,eacovgaminv_louis,gradsdes,&
	lintran,&
	obsweight,postdensity,theta,tolres,weight,&
	fitcummargwtsum,fitmargwtsum,fitpcummargwtsum,fitpmargwtsum,&
	obscummargwtsum,obsmargwtsum,obspcummargwtsum,obspmargwtsum,&
	presentedwtsum,&
	residacummargwtsum,residamargwtsum,&
	residcummargwtsum,residmargwtsum,residpcummargwtsum,&
	residpmargwtsum,stdobscummargwtsum,stdobsmargwtsum,&
	stdobspcummargwtsum,stdobspmargwtsum,&
	stdresidcummargwtsum,stdresidmargwtsum,&
	stdresidpcummargwtsum,stdresidpmargwtsum)

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
integer,intent(in)::catobsrange(:,:),dat(:,:),maxscore,maxw(:),minscore,minw(:),&
	npsu(:),nstratum,numcat(:),numcatobs(:),&
	psu(:),stratum(:),weight(:)
logical,intent(in)::complx,resid,stratify,usepsu
real(kind=8),intent(in)::beta(:),eacovgaminv_louis(:,:),&
	gradsdes(:,:),lintran(:,:),obsweight(:),&
	postdensity(:,:),theta(:,:,:),tolres
real(kind=8),intent(out)::fitcummargwtsum(minscore:maxscore),fitmargwtsum(minscore:maxscore),&
	fitpcummargwtsum(minscore:maxscore),fitpmargwtsum(minscore:maxscore),&
	obscummargwtsum(minscore:maxscore),obsmargwtsum(minscore:maxscore),&
	obspcummargwtsum(minscore:maxscore),obspmargwtsum(minscore:maxscore),presentedwtsum,&
	residacummargwtsum(minscore:maxscore),residamargwtsum(minscore:maxscore),&
	residcummargwtsum(minscore:maxscore),residmargwtsum(minscore:maxscore),&
	residpcummargwtsum(minscore:maxscore),residpmargwtsum(minscore:maxscore),&
	stdobscummargwtsum(minscore:maxscore),stdobsmargwtsum(minscore:maxscore),&
	stdobspcummargwtsum(minscore:maxscore),stdobspmargwtsum(minscore:maxscore),&
	stdresidcummargwtsum(minscore:maxscore),stdresidmargwtsum(minscore:maxscore),&
	stdresidpcummargwtsum(minscore:maxscore),stdresidpmargwtsum(minscore:maxscore)	
!	al is allocation indicator.	
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
!	score is a weighted sum
integer::al,col,counter,dimdesign,dimlatout,item,&
	ncat,ncatobs,nitems,nlin,nobs,nquad,nscale,&
	obs,quad,resp(size(dat,1)),row,score
!	datamask is a mask for items with variable weights.
!	skip indicates that an item is to be skipped.
logical::datamask(size(dat,1)),skip
!   avediff is an average difference.
!	covobs is the covariance matrix for adjusted observed values.
!	covsum is used to estimate asymptotic variances for fitted marginals.
!   diff is a difference.
!	dist is the conditional distribution.
!	locations gives location parameters.
!	newtheta is the transformed latent vector.
!	nobs is number of observations.
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

real(kind=8)::avediff,covobs(1,1),covsum(size(gradsdes,1)),diff,dist(minscore:maxscore),&
	locations(sum(numcat)),newtheta(size(lintran,1)),&
	probcat(sum(numcat)),probcatobs(sum(numcatobs)),&
	respprobvec(minscore:maxscore),respresid,respvec(minscore:maxscore),&
	scales(size(lintran,1),sum(numcat)),slopes(size(gradsdes,1)),totalobs
!	Set up parameter arrays for simplified processing.
real(kind=8),allocatable::obsmat(:,:,:)
dimdesign=size(gradsdes,1)
dimlatout=size(lintran,1)
ncat=sum(numcat)
ncatobs=sum(numcatobs)
nitems=size(dat,1)
nlin=ncat*(dimlatout+1)+1
nobs=size(obsweight)
nquad=size(theta,2)
nscale=ncat+1
fitcummargwtsum=0.0_8
fitmargwtsum=0.0_8
fitpcummargwtsum=0.0_8
fitpmargwtsum=0.0_8
locations=beta(1:ncat)
if(resid.or.complx)then
	allocate(obsmat(3,nobs,minscore:maxscore),stat=al)
	if(al/=0) stop 'Allocation of array failed for marginal distributions.'
	obsmat=0.0_8
end if
obscummargwtsum=0.0_8
obsmargwtsum=0.0_8
obspcummargwtsum=0.0_8
obspmargwtsum=0.0_8
presentedwtsum=0.0_8

if(resid)then
	residacummargwtsum=0.0_8
	residamargwtsum=0.0_8
	residcummargwtsum=0.0_8
	residmargwtsum=0.0_8
	residpcummargwtsum=0.0_8
	residpmargwtsum=0.0_8

end if

scales=reshape(beta(nscale:nlin-1),(/dimlatout,ncat/))
stdobscummargwtsum=0.0_8

stdobsmargwtsum=0.0_8
stdobspcummargwtsum=0.0_8
stdobspmargwtsum=0.0_8
if(resid)then
	stdresidcummargwtsum=0.0_8
	stdresidmargwtsum=0.0_8
	stdresidpcummargwtsum=0.0_8
	stdresidpmargwtsum=0.0_8
end if
totalobs=sum(obsweight)

!	Establish mask for items with variable weights.
datamask=.false.

do item=1,nitems
	if(maxw(item)>minw(item))datamask(item)=.true.
	
end do

!	Marginal totals.
do obs=1,nobs
	if(obsweight(obs)<=0.0_8) cycle
! Observation obs.
	resp=dat(:,obs)
	respvec=0.0_8
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
	presentedwtsum=presentedwtsum+obsweight(obs)
	score=0
	counter=1
	do item=1,nitems
		if(datamask(item))then
			score=score+weight(counter+resp(item))
		else
			score=score+minw(item)
		end if
		counter=counter+numcatobs(item)
	end do
	if(resid.or.complx)obsmat(3,obs,:)=1.0_8
	respvec(score)=1.0_8
	if(resid.or.complx)obsmat(1,obs,score)=1.0_8

	


!	Cycle through quadrature points.
	do quad=1,nquad
!		Probability computation.
		newtheta=matmul(lintran,theta(:,quad,obs))
		probcat=probvec(numcat,datamask,locations,scales,newtheta)
		probcatobs=probvecobsall(catobsrange,numcatobs,datamask,probcat)
		
		call distwtsum(maxw,minw,numcatobs,weight,&
			probcatobs,dist)
		
		respprobvec=respprobvec+postdensity(quad,obs)*dist
	end do
!	Set up obsmat and observed and fitted marginal totals.
	
	obsmargwtsum(score)=obsmargwtsum(score)+obsweight(obs)
	fitmargwtsum=fitmargwtsum+obsweight(obs)*respprobvec
	if(resid)obsmat(2,obs,:)=respvec-respprobvec
		
	

end do

!	The residual for marginal totals.
if(resid) residmargwtsum=obsmargwtsum-fitmargwtsum
!
!	Cumulative values.
do row=minscore,maxscore
	if(row==minscore)then
		fitcummargwtsum(row)=fitmargwtsum(row)
		obscummargwtsum(row)=obsmargwtsum(row)
	else
		if(row==maxscore)then
			fitcummargwtsum(row)=presentedwtsum
			obscummargwtsum(row)=presentedwtsum
		else
			fitcummargwtsum(row)=fitmargwtsum(row)+fitcummargwtsum(row-1)
			obscummargwtsum(row)=obsmargwtsum(row)+obscummargwtsum(row-1)
		end if
	end if
end do
if(resid)residcummargwtsum=obscummargwtsum-fitcummargwtsum
!	Fractions and variances.
if(presentedwtsum>0.0_8)then
	fitpcummargwtsum=fitcummargwtsum/presentedwtsum
	fitpmargwtsum=fitmargwtsum/presentedwtsum
	obspcummargwtsum=obscummargwtsum/presentedwtsum
	obspmargwtsum=obsmargwtsum/presentedwtsum
	if(resid)then
		residpcummargwtsum=residcummargwtsum/presentedwtsum
		residpmargwtsum=residmargwtsum/presentedwtsum
	end if
end if
if(presentedwtsum==0.0_8)return
!	Standard errors.
!	Simple random sampling.
if(.not.complx)then
	stdobscummargwtsum=sqrt(obscummargwtsum*(1.0_8-obscummargwtsum/totalobs))
	stdobsmargwtsum=sqrt(obsmargwtsum*(1.0_8-obsmargwtsum/totalobs))
	stdobspcummargwtsum=sqrt(obspcummargwtsum*(1.0_8-obspcummargwtsum)/presentedwtsum)
	stdobspmargwtsum=sqrt(obspmargwtsum*(1.0_8-obspmargwtsum)/presentedwtsum)			
	if(resid)then
		do row=minscore,maxscore
            covsum=reshape(crossinformation(.true.,obsmat(2:2,:,row),gradsdes,obsweight),(/dimdesign/))
            slopes=matmul(eacovgaminv_louis,covsum)
            avediff=0.0_8
			do obs=1,nobs
				if(obsweight(obs)<-0.0_8)cycle

                obsmat(2,obs,row)=obsmat(2,obs,row)-dot_product(gradsdes(:,obs),slopes)
                avediff=avediff+obsmat(2,obs,row)*obsweight(obs)
            end do

            avediff=avediff/totalobs

            do obs=1,nobs
                diff=obsmat(2,obs,row)-avediff
                obsmat(2,obs,row)=diff
                stdresidmargwtsum(row)=stdresidmargwtsum(row)+diff*diff*obsweight(obs)
                if(row==maxscore)cycle
                if(row>minscore)then
                    obsmat(2,obs,row)=obsmat(2,obs,row)+obsmat(2,obs,row-1)
                    diff=obsmat(2,obs,row)
                    stdresidcummargwtsum(row)=stdresidcummargwtsum(row)+diff*diff*obsweight(obs)
                 end if
            end do
            if(row==minscore)stdresidcummargwtsum(row)=stdresidmargwtsum(row)
        end do
	end if
else
!	Complex sampling.
	do row=minscore,maxscore
	
		covobs=information_complex(npsu,nstratum,psu,stratum,stratify,usepsu,&
					obsmat(1:1,:,row),obsweight)
		stdobsmargwtsum(row)=sqrt(covobs(1,1))
		
		
		covobs=information_complex(npsu,nstratum,psu,stratum,stratify,usepsu,&
			obsmat(1:1,:,row)-obspmargwtsum(row)*obsmat(3:3,:,row),obsweight)
		stdobspmargwtsum(row)=sqrt(covobs(1,1)/presentedwtsum)
		
		if(row>minscore.and.row<maxscore)then
				
			obsmat(1,:,row)=obsmat(1,:,row)+obsmat(1,:,row-1)
				
			covobs=information_complex(npsu,nstratum,psu,stratum,stratify,usepsu,&
				obsmat(1:1,:,row),obsweight)
			stdobscummargwtsum(row)=sqrt(covobs(1,1))
		
		
			covobs=information_complex(npsu,nstratum,psu,stratum,stratify,usepsu,&
				obsmat(1:1,:,row)-obspcummargwtsum(row)*obsmat(3:3,:,row),obsweight)
			stdobspcummargwtsum(row)=sqrt(covobs(1,1)/presentedwtsum)
		end if
		if(row==minscore) then!
			stdobscummargwtsum(row)=stdobsmargwtsum(row)
			stdobspcummargwtsum(row)=stdobsmargwtsum(row)
		
		end if
			

		
		if(resid)then
			covsum=reshape(crossinformation(.true.,obsmat(2:2,:,row),gradsdes,obsweight),(/dimdesign/))
            slopes=matmul(eacovgaminv_louis,covsum)
			
            do obs=1,nobs
                obsmat(2,obs,row)=obsmat(2,obs,row)-dot_product(gradsdes(:,obs),slopes)
            end do
            covobs=information_complex(npsu,nstratum,psu,stratum,stratify,usepsu,&
				obsmat(2:2,:,row),obsweight)
            stdresidmargwtsum(row)=covobs(1,1)
            if(row==maxscore)cycle
            if(row>minscore)then
                obsmat(2,:,row)=obsmat(2,:,row)+obsmat(2,:,row-1)


                covobs=information_complex(npsu,nstratum,psu,stratum,stratify,usepsu,&
                    obsmat(2:2,:,row),obsweight)
                stdresidcummargwtsum(row)=covobs(1,1)
            else
                stdresidcummargwtsum(row)=stdresidmargwtsum(row)
            end if

		end if

		
	end do
end if
if(resid)then
	do row=minscore,maxscore
		stdresidcummargwtsum(row)=sqrt(max(0.0_8,stdresidcummargwtsum(row)))
		stdresidpcummargwtsum(row)=stdresidcummargwtsum(row)/presentedwtsum
		if(stdresidcummargwtsum(row)>tolres*stdobscummargwtsum(row).and.stdobscummargwtsum(row)>0.0_8)&
				residacummargwtsum(row)=residcummargwtsum(row)/stdresidcummargwtsum(row)
		stdresidmargwtsum(row)=sqrt(max(0.0_8,stdresidmargwtsum(row)))
		stdresidpmargwtsum(row)=stdresidmargwtsum(row)/presentedwtsum
		if(stdresidmargwtsum(row)>tolres*stdobsmargwtsum(row).and.stdobsmargwtsum(row)>0.0_8)&
			residamargwtsum(row)=residmargwtsum(row)/stdresidmargwtsum(row)
	end do
end if	
return
if(resid.or.complx)deallocate(obsmat)
end subroutine margindistwtsum

