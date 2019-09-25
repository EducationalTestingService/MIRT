!	Find posterior distribution for latent vectors.

!	catobsrange is used to relate observed and underlying categories.
!	Responses are in dat.
!	maxita is used to determine the number of iterations to use to find the maximum posterior density.
!	numcat provides ranges for underlying observations,
!	numcatobs provides ranges of observations.
!	mask indicates items to use for computations.
!	normal determines if a normal or multinomial analysis is used.
!	beta is the parameter vector.
!	changemin is used for minimum change.
!	Predictors are in indvar.
!	lintran is the linear transformation for the latent vector.
!	maxdalpha is used to limit changes in alpha approximations.
!	quadpoint provides quadrature points.
!	quadweight provides quadrature weights.
!	tau controls shrinking of step size.
!	tola is used to determine the desired accuracy of the maximization.
!	tolsing is used to determine the tolerance for the Cholesky decomposition.
!	alpha contains locations of maxima.
!	cholnhess contains modified Cholesky decompositions.
!	postdensity gives weights at quadrature points.
!	theta are quadrature points.
subroutine posterior(catobsrange,dat,maxita,numcat,numcatobs,mask,normal,beta,changemin,indvar,lintran,&
	maxdalpha,quadpoint,quadweight,&
	tau,tola,tolsing,alpha,cholnhess,postdensity,theta)
implicit none
interface

	function chol(sym,tolsing)
!	chol is used to compute a modified Cholesky decomposition of the n by n matrix sym.
!	The singularity tolerance is tolsing.
		implicit none
		real(kind=8),intent(in)::sym(:,:),tolsing
	 	real(kind=8)::chol(size(sym,1),size(sym,1))
	end function chol

!	Get normal density for scaled latent vector.  Place in density.
!	Input involves linear component linth, quadrature points quadpoint, quadratic components quadth, quadrature weights quadweight, and scale
!	factor scaletheta.
	function density(linth,quadpoint,quadth,quadweight,scaletheta)
		implicit none
		real(kind=8),intent(in)::linth(:),quadpoint(:,:),quadth(:,:),quadweight(:),scaletheta
		real(kind=8)::density(size(quadweight))
	end function density		
!	Get multinomial probabilities for latent vector.  Place in density.
!	Input involves linear component linth, quadrature points quadpoint, quadratic components quadth, and quadrature weights quadweight.
	function densitym(linth,quadpoint,quadth,quadweight)
		real(kind=8),intent(in)::linth(:),quadpoint(:,:),quadth(:,:),quadweight(:)
		real(kind=8)::densitym(size(quadweight))
	end function densitym	

!	convert an n by n symmetric matrix from compact to regular form.
	function expandmat(n,compmat)
		implicit none
		integer,intent(in)::n
		real(kind=8),intent(in)::compmat(:)
		real(kind=8)::expandmat(n,n)
	end function expandmat

!	Extract the linear component lintheta and the quadratic component quadtheta from the parameter vector beta.
	subroutine getlinquad(beta,lintheta,quadtheta)
		implicit none
		real(kind=8),intent(in)::beta(:)
		real(kind=8),intent(out)::lintheta(:,:),quadtheta(:,:,:)
	end subroutine getlinquad

!	Extract the linear component linth and the quadratic component quadth from lintheta and quadtheta for a given predictor thpr.
	subroutine getlinquadth(lintheta,quadtheta,thpr,linth,quadth)
		implicit none
		real(kind=8),intent(in)::lintheta(:,:),quadtheta(:,:,:),thpr(:)
		real(kind=8),intent(out)::linth(:),quadth(:,:)
	end subroutine getlinquadth
!	
!	Find location of maximum posterior density of the latent vector.
!	catobsrange is used to relate observed and underlying categories.
!	maxita is used to determine the number of iterations to use to find the maximum posterior density.
!	maxitb is used for subiterations.
!	npred is the number of predictors.
!	numcat provides ranges for underlying observations,
!	numcatobs provides ranges of observations, mask indicates items to use for computations.
!	obsmask is an observation mask that can be used to restrict attention to a portion of the response.
!	beta is the parameter vector.
!	changemin is used for minimum change.
!	linth is the linear component of the log density.
!	lintran is the linear transformation of the latent vector.
!	maxdtheta is the maximum permitted change in the approximation of the maximum in one step.
!	quadth is the quadratic component of the log density.
!	resp is the response.
!	tau is used for minimum step sizes.
!	tol is used to determine the desired accuracy of the maximization.
!	tolsing is used to determine the tolerance for the Cholesky decomposition.
!	nhchol is the modified Cholesky decomposition of the negative hessian.
!	theta is the maximum.
	subroutine maxpost(catobsrange,maxita,maxitb,npred,numcat,numcatobs,resp,obsmask,beta,changemin,&
		linth,lintran,maxdtheta,quadth,tau,tol,tolsing,&
		nhchol,theta)
		implicit none
		integer,intent(in)::catobsrange(:,:),maxita,maxitb,npred,numcat(:),numcatobs(:),resp(:)
		logical,intent(in)::obsmask(:)
		real(kind=8),intent(in)::beta(:),changemin,lintran(:,:),linth(:),maxdtheta,quadth(:,:),tau,tol,tolsing
		real(kind=8),intent(inout)::theta(:)
		real(kind=8),intent(out)::nhchol(:,:)
	end subroutine maxpost

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

!	Convert from adaptive quadrature scale to normalized latent vector.
!	Use nhchol for modified Cholesky decomposition of negative hessian.
!	Use oldalpha for location.
!	Use scalez for scaling.
!	Quadrature point is z.
	function scaleadapt(nhchol,oldalpha,scalez,z)
		implicit none
		real(kind=8),intent(in)::nhchol(:,:),oldalpha(:),scalez(:),z(:)
		real(kind=8)::scaleadapt(size(z))
	end function scaleadapt

!		This program uses the modified Cholesky decomposition
!		in lu to solve the equation lusolve = b.
	function solve(lu,b)
		implicit none
		
		real(kind=8),intent(in)::lu(:,:)
		real(kind=8),intent(in)::b(:)
		real(kind=8)::solve(size(b))
	end function solve
end interface
integer,intent(in)::catobsrange(:,:),dat(:,:),maxita,numcat(:),numcatobs(:)
logical,intent(in)::mask(:),normal
real(kind=8),intent(in)::beta(:),changemin,indvar(:,:),lintran(:,:),&
	maxdalpha,quadpoint(:,:),quadweight(:),&
	tau,tola,tolsing
real(kind=8),intent(inout)::alpha(:,:),cholnhess(:,:,:)
real(kind=8),intent(out)::postdensity(:,:),theta(:,:,:)


!	dimlatin is the dimension of the latent vector.
!	dimlatout is the dimension of the transformed latent vector.

!	ncat is the total number of underlying categories.
!	nitems is the number of items, and item is an item.
!	nobs is number of observations.
!	npred is the number of predictors.
!	nquad is the number of quadrature points.
!	nscale is the pointer for scale parameters.
!	nscale1 points to the end of the scale parameters.
!	obs is observation number.
!	pred is a predictor counter.
!	quad counts quadrature points.
!	resp is a single response.
!	row is a row counter.
integer::dimlatin,dimlatout,item,ncat,nitems,&
	nobs,npred,nquad,nscale,&
	nscale1,obs,pred,quad,&
	resp(size(mask)),row
!	Data mask.
logical::datamask(size(mask))
!	alpha is used to find the maximum posterior density.
!	locations is for location parameters.
!	linth is a linear component of the log density of a latent vector.
!	lintheta is for linear parameters for the log density.
!	meanth is the mean of the latent vector.
!	newtheta is lintran(theta).
!	nhchol is a Cholesky decomposition of -2*quadth.
!	priordensity provides a prior density.
!	probcat is used for underlying marginal conditional probabilities.
!	probcatobs is used for observed marginal conditional probabilities.
!	probobs is the probability of an observation.
!	quadth is a quadratic component of the log density of a latent vector.
!	quadtheta is for quadratic parameters for the log density.
!	scales is for scale parameters.
!	scaletheta is used for scaling.
!	scalez is used for scaling.
!	thpr is the predictor vector for an observation.
!	z is a transformation of a latent vector value.
real(kind=8)::linth(size(quadpoint,1)),&
	lintheta(size(quadpoint,1),size(indvar,1)),locations(sum(numcat)),&
	newtheta(size(lintran,1)),meanth(size(quadpoint,1)),&
	nhchol(size(quadpoint,1),size(quadpoint,1)),&
	priordensity(size(quadweight)),probcat(sum(numcat)),probcatobs(size(mask)),probobs,&
	quadth(size(quadpoint,1),size(quadpoint,1)),&
	quadtheta(size(quadpoint,1),size(quadpoint,1),size(indvar,1)),&
	scales(size(lintran,1),sum(numcat)),scaletheta,&
	scalez(size(quadpoint,1)),&
	thpr(size(indvar,1)),z(size(quadpoint,1))
	

!	Set up parameter arrays for simplified processing.
dimlatin=size(quadpoint,1)

dimlatout=size(lintran,1)
ncat=sum(numcat)
nitems=size(mask)
npred=size(indvar,1)

nobs=size(dat,2)
!	Item scale parameters.
nscale=ncat+1

nscale1=ncat*(1+dimlatout)
nquad=size(quadweight)
!	Parameters for quadratic component of log density.

locations=beta(1:ncat)
scales=reshape(beta(nscale:nscale1),(/dimlatout,ncat/))

call getlinquad(beta,lintheta,quadtheta)



!	Obtain results for each observation.
do obs=1,nobs
	
! Observation obs.
	resp=dat(:,obs)

!	Density specification for observation.
	thpr=indvar(:,obs)

	call getlinquadth(lintheta,quadtheta,thpr,linth,quadth)


!	Mask
	datamask=.false.
	do item=1,nitems
		if(resp(item)>=0.and.resp(item)<numcatobs(item).and.mask(item))datamask(item)=.true.
	end do	
	
!	The normal versus multinomial cases differ here.
	if(normal)then
!	Mean and covariance matrix of latent vector.

		nhchol=chol(-2.0_8*quadth,tolsing)
!	Check that covariance matrix is positive definite.
		if(minval((/(nhchol(row,row),row=1,dimlatin)/))<=0.0_8)stop "Unacceptable covariance matrix"
		meanth=solve(nhchol,linth)
		scaletheta=sqrt(product((/(nhchol(row,row),row=1,dimlatin)/)))
		scaletheta=scaletheta*exp(-0.5_8*dot_product(linth,meanth))		
		

!	See if adaptive approach to be used.
		if(maxita>0) then

			call maxpost(catobsrange,maxita,maxita,npred,numcat,numcatobs,resp,datamask,beta,changemin,&
				linth,lintran,maxdalpha,quadth,tau,tola,tolsing,&
				nhchol,alpha(:,obs))
			cholnhess(:,:,obs)=nhchol
		end if
		nhchol=cholnhess(:,:,obs)
		do row=1,dimlatin
			scalez(row)=1.0_8/sqrt(nhchol(row,row))
		end do	
		scaletheta=scaletheta*product(scalez)
		do quad=1,nquad
				theta(:,quad,obs)=scaleadapt(nhchol,alpha(:,obs),scalez,quadpoint(:,quad))
		end do
		priordensity=density(linth,theta(:,:,obs),quadth,quadweight,scaletheta)
	else
		theta(:,:,obs)=quadpoint
		priordensity=densitym(linth,quadpoint,quadth,quadweight)
	end if		
!	Cycle through quadrature points.
	
	do quad=1,nquad
		
		
				


			
! Linear transformation of latent vector.
		newtheta=matmul(lintran,theta(:,quad,obs))

		
! Probability computation.


		probcat=probvec(numcat,datamask,locations,scales,newtheta)
		probcatobs=probvecobs(catobsrange,numcatobs,resp,datamask,probcat)
		postdensity(quad,obs)=priordensity(quad)
		if(any(datamask))postdensity(quad,obs)=postdensity(quad,obs)*product(probcatobs,datamask)
	end do
	probobs=sum(postdensity(:,obs))
	postdensity(:,obs)=postdensity(:,obs)/probobs	
end do
return
end subroutine posterior
