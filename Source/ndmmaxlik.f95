!	Maximum likelihood for log-linear MIRT models with possible indirect observation.
!	Summary of variables.
!	comment is used to comment on the iteration progress report.
!	catobsrange is the range of underlying categories per for each observed category.
!	constdim1 is the number of linear constraints not treated via least squares.
!	dat is the data matrix.
!	maxit is the maximum number of iterations.
!	maxita is the maximum number of iterations used to find alp.
!	maxitb is the maximum number of iterations used to find the step size.
!	numcat is the number of underlying categories per item.
!	numcatobs is the number of observed categories per item.
!	rbetamap is the reduced beta map.
!	unitno is the unit number for reporting iteration progress.
!	normal is for a multivariate normal distribution of the latent vector.  The alternative is multinomial.
!	nr is for Newton-Raphson.  The alternative uses an approximation to Newton-Raphson based on the Louis
!		approximation.
!	printprog indicates that iteration progress is to be printed to unitno.
!	printprogstd indicates that iteration progress is to go to standard output.
!	proport is .true. if the constraint sum of squares is proportional to the sum of the observation weights.
!	changemin is the criterion for minimum change.
!	constmat is the constraint matrix.
!	constvec is the constraint vector.
!	design is the design matrix.
!	indvar is the matrix of independent variables.
!	kappa is the maximum permitted step size.
!	lintran is the linear transformation from basic quadrature points to values of the latent vector.
!	maxdalpha is the maximum change in thee alpha maximum norm permitted per iteration.
!	obsweight is the vector of observation weights.
!	offset is the offset vector.
!	quadpoint is the quadrature points.
!	quadweight is the quadrature weights.
!	rdesign is the reduced design matrix.
!	tau is a step multiplier.
!	tol is the convergence criterion for ent change.
!	tola is the convergence criterion for alpha change.
!   tolc is the criterion for appropriate maximum ratio of diagonal terms in modified Cholesky decomposition.
!	tolsing is a criterion for singularity.
!	alpha is the matrix of approximate maxima of the individual log likelihood.
!	beta is the parameter vector.
!	beta=design(gamma).
!	cholnhess is the modified Cholesky decomposition of the negative hessian at the corresponding alpha.
!	grad is the gradient vector relative to beta=design(gamma).
!	graddes is the gradient vector relative to gamma.
!	grads is a matrix of individual gradient vectors.
!	loglik is the log likelihood.
!	nhess is the negative hessian matrix relative to beta.
!	nhessdes is the negative hessian matrix relative to gamma.
!	nhessdeschol is the modified Cholesky decomposition of nhessdesplus.
!	nhessdesplus is nhessdes plus the component constquad for constraints.
!	postdensity is the array of posterior densities at adaptive quadrature points.
!	prob contains item probabilities.
!	theta is the array of adaptive quadrature points.
subroutine ndmmaxlik(comment,catobsrange,constdim1,dat,maxit,maxita,maxitb,numcat,numcatobs,rbetamap,unitno,&
	normal,nr,printprog,printprogstd,proport,&
	changemin,constmat,constvec,design,indvar,kappa,lintran,&
	maxdalpha,obsweight,offset,quadpoint,quadweight,rdesign,&
	tau,tol,tola,tolc,tolsing,&
	alpha,beta,cholnhess,gamma,constquad,grad,graddes,grads,&
	loglik,nhess,nhessdes,nhessdeschol,nhessdesplus,postdensity,prob,theta)
implicit none
interface
!	
	function chol(sym,tolsing)
!	chol is used to compute a modified Cholesky decomposition of the n by n matrix sym.
!	The singularity tolerance is tolsing.
		implicit none
		real(kind=8),intent(in)::sym(:,:),tolsing
	 	real(kind=8)::chol(size(sym,1),size(sym,1))
	end function chol
!	condprobvec is used to compute observed item conditional probabilities of underlying item categories given
!	observed item categories.
!	catobsrange provides ranges of underlying item categories that correspond to observed item categories.
!	numcatobs provides the number of observed categories per item.
!	resp is the item response.
!	mask indicates which items were presented.
!	probcat is the vector of underlying category probabilities.
!	probcatobs is the vector of observed category probabilities.
	function condprobvec(catobsrange,numcatobs,resp,mask,probcat,probcatobs)
		implicit none
		integer,intent(in)::catobsrange(:,:),numcatobs(:),resp(:)
		logical,intent(in)::mask(:)
		real(kind=8),intent(in)::probcat(:),probcatobs(:)
		real(kind=8)::condprobvec(size(probcat))
	end function condprobvec

!	Covariances of products of elements of latent vector and of elements of latent vector with normal
!	distribution with mean mean and  covariance cov.
	function cubicnormal(cov,mean)
		implicit none
		real(kind=8),intent(in)::cov(:,:),mean(:)
		real(kind=8)::cubicnormal(size(mean),size(mean),size(mean))
	end function cubicnormal
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

!	gradient of log conditional probability of an observation.
!	rbetamap is the reduced beta map.
!	mask is the response mask.
!	condprobcat contains conditional item probabilities.
!	meanth is the mean of theta.
!	mean2th is the mean of products of elements of theta.
!	newtheta is the transformed latent vector.
!	probcat contains item probabilities.
!	theta is the latent vector.
!	thpr is the independent vector.
	function gradvec(rbetamap,mask,condprobcat,meanth,mean2th,newtheta,probcat,theta,thpr)
		implicit none
		integer,intent(in)::rbetamap(:,:)
		logical,intent(in)::mask(:)
		real(kind=8),intent(in)::condprobcat(:),meanth(:),mean2th(:,:),newtheta(:),probcat(:),theta(:),thpr(:)
		real(kind=8)::gradvec(size(rbetamap,2))
	end function gradvec

!	Estimate the information matrix with the Louis approach.
!	adjust is .true. if the average of grads is to be subtracted.
!	grads is the array of observation gradients.
!	obsweight contains observation weights.
	function information(adjust,grads,obsweight)
		implicit none
		logical,intent(in)::adjust
		real(kind=8),intent(in)::grads(:,:),obsweight(:)
		real(kind=8)::information(size(grads,1),size(grads,1))	
	end function information
!

!	Compute inverse for n by n matrix with modified Cholesky
!	decomposition tri.
	function invert(tri)
		implicit none
		real(kind=8),intent(in)::tri(:,:)
		real(kind=8)::invert(size(tri,1),size(tri,1))
	end function invert
	
!	Report on iteration progress.
!	comment is used to comment on the iteration progress report.
!	it is the iteration number.
!	unitno is the unit number for the report.
!	printprog indicates that iteration progress is to be printed to unitno.
!	printprogstd indicates that iteration progress is to go to standard output.
!	loglik is the log likelihood.
!	step is the step size. 
	subroutine iterationreport(comment,it,unitno,printprog,printprogstd,loglik,step)
		implicit none
		character(len=*),intent(in)::comment
		integer,intent(in)::it,unitno
		logical,intent(in)::printprog,printprogstd
		real(kind=8),intent(in)::loglik,step
	end subroutine iterationreport	

	subroutine makesym(matrix)
!	Symmetrize n by n matrix based on lower triangle.
		implicit none
		real(kind=8),intent(inout)::matrix(:,:)
	end subroutine makesym

!	Find location of maximum posterior density of the latent vector.
!	catobsrange is used to relate observed and underlying categories.
!	maxita is used to determine the number of iterations to use to find the maximum posterior density.
!	maxitb is used for subiterations.
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
!	thpr is the predictor.
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
	
!	Multinomial covariance matrix.
!	density is vector of probabilities.
!	meanth is vector of means.
!	quadpoint is array of weights.
	function mcovth(density,meanth,quadpoint)
		implicit none
		real(kind=8),intent(in)::density(:),meanth(:),quadpoint(:,:)
		real(kind=8)::mcovth(size(meanth),size(meanth))	
	end function mcovth
	
!	Multinomial third moments.  Covariance of product and single term.
!	density contains weights.
!	meanth contains means.
!	mean2th contains means of products.
!	quadpoint contains points.

	function mcubth(density,meanth,mean2th,quadpoint)
		implicit none
		real(kind=8),intent(in)::density(:),meanth(:),mean2th(:,:),quadpoint(:,:)
		real(kind=8)::mcubth(size(meanth),size(meanth),size(meanth))
	end function mcubth
	
!	Multinomial mean.
!	density is vector of probabilities
!	quadpoint is array of weights.
	function mmeanth(density,quadpoint)
		implicit none
		real(kind=8),intent(in)::density(:),quadpoint(:,:)
		real(kind=8)::mmeanth(size(quadpoint,1))	
	end function mmeanth

!	Means of products.  meanth contains means.
!	covth contains covariances.
!	quadpoint contains points.  density contains weights.
	function mprodth(covth,meanth)
		implicit none
		real(kind=8),intent(in)::covth(:,:),meanth(:)
		real(kind=8)::mprodth(size(meanth),size(meanth))
	end function mprodth	

!	Multinomial fourth moments.  Covariances of cross-products.
!	density provides weights.
!	points are in quadpoint.
!	Product means are in mean2th.
	function mquarth(density,mean2th,quadpoint)
		implicit none
		real(kind=8),intent(in)::density(:),mean2th(:,:),quadpoint(:,:)
		real(kind=8)::mquarth(size(mean2th,1),size(mean2th,1),size(mean2th,1),size(mean2th,1))
	end function mquarth

!	negative hessian of log conditional probability of an observation.
!	rbetamap is the reduced beta map
!	mask indicates the items presented.
!	condprobcat contains conditional item probabilities.
!	covth is the covariance matrix of the latent vector.
!	cubth provides covariances of products of latent vector elements and of latent vector elements.
!	probcat contains item probabilities.
!	quarth provides covariances of products of latent vector elements.
!	theta is the latent vector.
!	thpr is the independent vector.

	function nhessvec(rbetamap,mask,condprobcat,covth,cubth,newtheta,probcat,quarth,theta,thpr)
		implicit none
		integer,intent(in)::rbetamap(:,:)
		logical,intent(in)::mask(:)
		real(kind=8),intent(in)::condprobcat(:),covth(:,:),cubth(:,:,:),newtheta(:),probcat(:),quarth(:,:,:,:),theta(:),thpr(:)
		real(kind=8)::nhessvec(size(rbetamap,2),size(rbetamap,2))
	end function nhessvec

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



!	Covariances of products of a normal random vector with covariance matric cov and mean mean.
	function quarticnormal(cov,mean)
		implicit none
		real(kind=8),intent(in)::cov(:,:),mean(:)
		real(kind=8)::quarticnormal(size(mean),size(mean),size(mean),size(mean))
	end function quarticnormal
	
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
!

!		This program uses the modified Cholesky decomposition
!		in lu to solve the equation lusolve = b.
	function solve(lu,b)
		implicit none
		
		real(kind=8),intent(in)::lu(:,:)
		real(kind=8),intent(in)::b(:)
		real(kind=8)::solve(size(b))
	end function solve





end interface
character(len=*),intent(in)::comment
integer,intent(in)::catobsrange(:,:),constdim1,dat(:,:),maxit,maxita,maxitb,numcat(:),numcatobs(:),rbetamap(:,:),unitno
logical,intent(in)::normal,nr,printprog,printprogstd,proport
real(kind=8),intent(in)::changemin,constmat(:,:),constvec(:),design(:,:),indvar(:,:),kappa,lintran(:,:),maxdalpha,&
	obsweight(:),offset(:),quadpoint(:,:),quadweight(:),&
	rdesign(:,:),tau,tol,tola,tolc,tolsing
real(kind=8),intent(inout)::alpha(:,:),beta(:),cholnhess(:,:,:),gamma(:)
real(kind=8),intent(out)::constquad(:,:),grad(:),graddes(:),grads(:,:),loglik,&
	nhess(:,:),nhessdes(:,:),nhessdeschol(:,:),nhessdesplus(:,:),postdensity(:,:),prob(:),theta(:,:,:)
!	al is allocation error flag.
!	col is the column number.
!	constdim is the total number of constraints.
!	dimlatin is the dimension of the latent vector.
!	dimlatout is the dimension of the transformed latent vector.

!	it is the iteration number.
!	ita is the iteration number for determination of the adaptive quadrature.
!	item is the item number.
!	ncat is the total number of underlying categories.
!	ncatobs is the total number of observed categories.
!	nitems is the number of items.
!	nobs is the total number of observations.
!	npred is the number of independent variables.

!	nquad is the number of quadrature points.

!	nreduce is the dimension of the reduced version of beta.
!	nscale is the pointer for scale parameters.
!	nscale1 points to the last scale parameter.
!	obs is the  observation number.
!	pred is the predictor number.

!	quad is the quadrature point index.
!	resp is the response vector for the observation under study.
!	row is the row number.
!	rowpt is the row pointer.


integer::al,col,constdim,dimlatin,dimlatout,it,ita,item,ncat,&
	nobs,nitems,npred,nquad,nreduce,nscale,nscale1,&
	obs,pred,quad,resp(size(dat,1)),row,rowpt

!	error is used to indicate a problem with a covariance matrix.	
!	finish indicates whether iterations end.
!	linflag indicates if linear terms vary.
!	mask is .true. if the  item is  present and .false. otherwise.
!	nrskip is used to skip Newton-Raphson.
!	quadflag indicates if  quadratic terms vary.
logical::error,finish,linflag,mask(size(dat,1)),nrskip,quadflag

!	loc is the location vector for the latent vector.
!	lowtri is the lower triangular matrix in the Cholesky decomposition of the
!	covariance matrix of the latent vector.
!	probcat is a vector of conditional underlying item category probabilities.
!	probcatobs is a vector of conditional observed item category probabilities.
!	condprobcat is a vector of conditional probabilities of underlying item
!	categories given observed category probabilities.
!	nhchol is a matrix used in Cholesky decompositions for transformation of latent vectors.
!	gradth is a gradient with respect to the standardized latent vector value, and nhessth
!	is a negative hessian with respect
!	to the standardized latent vector value.
!	alphastep is a proposed change in alpha.
!	betastep is the change in beta proposed.
!	change is the change in weighted log likelihood.
!	condprob is a probability component.
!	consterr is the constraint error.
!	constweight is the multiplier of the sum of squares of constraint errors.
!	covth is the covariance matrix of the latent vector.
!	cubth is the matrix of third central moments of the latent vector.
!	der is the directional derivative of the weighted log likelihood.
!	gammastep is the proposed change in gamma.
!	gradlog is an array of gradients of logs of posterior densities at the latent vector.
!	gradobs is the gradient with respect to beta for a single observation.
!	lintheta is a matrix used for the log density of theta.  It provides a linear contribution.  It is 
!	used in conjunction wih the vector linth.
!	locations is the vector of location parameters.
!	logcondprob is a logarithm of a probability component.
!	loglikn is the old log likelihood within the iteration.
!	loglikold is the previous log likelihood from the last iteration.
!	maxalphastep is the maximum norm of alphastep.
!	maxbeta is the maximum norm of betastep.
!	mean2th is the mean of products of elements of the latent vector.
!	meanth is the mean of the latent vector.
!   maxchol is the maximum diagonal in a Cholesky decomposition.
!   minchol is the minimum diagonal in a Cholesky decomposition.
!	newtheta is lintran(theta)
!	nhessobs is the negative hessian with respect to beta for a single observation.
!	nhessquad is a negative hessian of the log of a posterior density at the latent vector.
!	oldalpha is the previous alpha.
!	oldbeta is the previous beta.
!	oldconsterr is  the last constraint error.
!	oldlogcondprob is the previous value of logcondprob.
!	priordensity is the prior density of the latent vector at quadrature points.
!	quadtheta is ann array of matrices used for the log density of theta.  It is used  in conjunction
!	with quadth.
!	quarth is the matrix of fourth central moments of the latent vector minus products of covariances.
!	rowdiv is a divisor for a variable.
!	scale is the vector of scale factors for the latent vector.
!	scales is the matrix of scale factors.
!	scaletheta is the square root of the inverse of the determinant of covth.
!	scalez is used to scale latent vector values for adaptive quadrature.
!	scalezp is the product of the elements of scalez.
!	step is the step size.
!	stepmin is the minimum step size.
!	thpr is the vector of theta predictors for a particular observation.
!	v is a matrix used for constraints not based on least squares.
!	vv is a second matrix used for constraints not based on least squares.
!	xitemcount is a counter for the number of items presented for an observation.
!	z is a latent vector value in standard units.
real(kind=8)::alphastep(size(quadpoint,1)),&
	betastep(size(beta)),&
	change,condprob,condprobcat(sum(numcat)),consterr(size(constvec)),constweight,&
	covth(size(quadpoint,1),size(quadpoint,1)),scaletheta,&
	cubth(size(quadpoint,1),size(quadpoint,1),size(quadpoint,1)),&	
	der,&
	gammastep(size(gamma)),gradlog(size(rbetamap,2)),gradobs(size(rdesign,1)),&
	gradth(size(quadpoint,1)),&
	linth(size(quadpoint,1)),lintheta(size(quadpoint,1),size(indvar,1)),&
	locations(sum(numcat)),loglikn,loglikold,logcondprob,&
	maxalphastep,maxbeta,maxchol,&
    mean2th(size(quadpoint,1),size(quadpoint,1)),meanth(size(quadpoint,1)),&
	minchol,newtheta(size(lintran,1)),nhchol(size(quadpoint,1),size(quadpoint,1)),&
	nhessobs(size(rdesign,1),size(rdesign,1)),&
	nhessquad(size(rbetamap,2),size(rbetamap,2)),nhessth(size(quadpoint,1),size(quadpoint,1)),&
	oldalpha(size(quadpoint,1)),oldbeta(size(beta)),oldconsterr(size(constvec)),&
	oldlogcondprob,&
	priordensity(size(quadweight)),probcat(sum(numcat)),&
	probcatobs(size(dat,1)),&
	quadth(size(quadpoint,1),size(quadpoint,1)),&
	quadtheta(size(quadpoint,1),size(quadpoint,1),size(indvar,1)),&
	quarth(size(quadpoint,1),size(quadpoint,1),size(quadpoint,1),size(quadpoint,1)),&
	scales(size(lintran,1),sum(numcat)),step,stepmin,&
	scalez(size(quadpoint,1)),scalezp(size(obsweight)),&
	thpr(size(indvar,1)),&
	xitemcount,&
	
	z(size(quadpoint,1))
real(kind=8),allocatable::v(:,:),vv(:,:)
!	Some sizes.
constdim=size(constvec)
dimlatin=size(quadpoint,1)
dimlatout=size(lintran,1)
ncat=sum(numcat)
nitems=size(dat,1)
nobs=size(obsweight)

nquad=size(quadweight)

nreduce=size(rdesign,1)
npred=size(indvar,1)
if(constdim1>0)then
	allocate(v(size(gamma),constdim1),vv(constdim1,constdim1),stat=al)
	if(al/=0)stop 'Allocation of constraint arrays failed.'
end if
!	Set up pointers for beta elements.
!	Item scale parameters.
nscale=ncat+1

nscale1=ncat*(1+dimlatout)
!	Parameters for quadratic component of log density.

!	Check for linear and quadratic terms
linflag=.false.
quadflag=.false.
if(any(rbetamap(3,:)>0.and.rbetamap(5,:)==0))linflag=.true.
if(any(rbetamap(3,:)>0.and.rbetamap(4,:)>0))quadflag=.true.
if(constdim1<constdim) then
	constquad=2.0_8*matmul(constmat(:,constdim1+1:constdim),transpose(constmat(:,constdim1+1:constdim)))
	constweight=1.0_8
	if(proport)then
		constweight=sum(obsweight)
		constquad=constweight*constquad
		
	end if

	
end if




!	Up to maxit iterations.
finish=.false.

do it=0,maxit



	if(it>=maxit)finish=.true.
!	Initialize gradient and hessian.
	grad=0.0_8
	nhess=0.0_8
	prob=0.0_8

!	Initialize information measure.


	loglik=0.0_8
	if(constdim>0)then
		consterr=constvec-matmul(transpose(constmat),gamma)
		if(constdim1<constdim)loglik=-constweight*sum(consterr(constdim1+1:constdim)*consterr(constdim1+1:constdim))
		
	end if
	
!	Set up parameter arrays for simplified processing.

	locations=beta(1:ncat)

	scales=reshape(beta(nscale:nscale1),(/dimlatout,ncat/))

	call getlinquad(beta,lintheta,quadtheta)
	

	do obs=1,nobs			
		
		if(obsweight(obs)<=0.0_8) cycle
		
! Observation obs.
		resp=dat(:,obs)

		mask=.false.
		do item=1,nitems
			if(resp(item)>=0.and.resp(item)<numcatobs(item))mask(item)=.true.
		end do

		
!	Density specification for observation.
		thpr=indvar(:,obs)

		call getlinquadth(lintheta,quadtheta,thpr,linth,quadth)
!	The normal versus multinomial cases differ here.
		if(normal)then
!	Mean and covariance matrix of latent vector.

			nhchol=chol(-2.0_8*quadth,tolsing)
!	Check that covariance matrix is positive definite.
			if(minval((/(nhchol(row,row),row=1,dimlatin)/))<=0.0_8)stop "Unacceptable covariance matrix"
			scaletheta=sqrt(product((/(nhchol(row,row),row=1,dimlatin)/)))
			

			meanth=solve(nhchol,linth)
			scaletheta=scaletheta*exp(-0.5_8*dot_product(linth,meanth))
			covth=invert(nhchol)
			mean2th=mprodth(covth,meanth)
			if(linflag.and.quadflag)cubth=cubicnormal(covth,meanth)
			if(quadflag)quarth=quarticnormal(covth,meanth)

! 	Location maximum of likelihood contribution for examinee k.
			

!	See if adaptive approach to be used.
			if(maxita>0) then

				call maxpost(catobsrange,maxita,maxita,npred,numcat,numcatobs,resp,mask,beta,changemin,&
					linth,lintran,maxdalpha,quadth,tau,tol,tolsing,&
					nhchol,alpha(:,obs))
				cholnhess(:,:,obs)=nhchol

			end if
		end if
!	Initialize gradient and negative hessian.  Initialize probability.

		
		
		gradobs=0.0_8

		nhessobs=0.0_8


!	Scale factors and priordensities.
		if(normal)then
			nhchol=cholnhess(:,:,obs)
			oldalpha=alpha(:,obs)
			


			do row=1,dimlatin
				scalez(row)=1.0_8/sqrt(nhchol(row,row))
			end do
			scalezp(obs)=product(scalez)
			scaletheta=scaletheta*scalezp(obs)
			do quad=1,nquad
				theta(:,quad,obs)=scaleadapt(nhchol,oldalpha,scalez,quadpoint(:,quad))
			end do
			priordensity=density(linth,theta(:,:,obs),quadth,quadweight,scaletheta)

		else
			theta(:,:,obs)=quadpoint
			priordensity=densitym(linth,quadpoint,quadth,quadweight)
			if(linflag.or.quadflag)then
				meanth=mmeanth(priordensity,quadpoint)
				covth=mcovth(priordensity,meanth,quadpoint)
				if(quadflag)then
					mean2th=mprodth(covth,meanth)
					cubth=mcubth(priordensity,meanth,mean2th,quadpoint)
					quarth=mquarth(priordensity,mean2th,quadpoint)
				end if
			end if


		end if
!	Cycle through quadrature points.

		do quad=1,nquad
			
			

!	Probability computation.

			newtheta=matmul(lintran,theta(:,quad,obs))
			probcat=probvec(numcat,mask,locations,scales,newtheta)

			probcatobs=probvecobs(catobsrange,numcatobs,resp,mask,probcat)

			postdensity(quad,obs)=product(probcatobs,mask)*priordensity(quad)

		
			prob(obs)=prob(obs)+postdensity(quad,obs)

			condprobcat=condprobvec(catobsrange,numcatobs,resp,mask,probcat,probcatobs)

!	Gradient and negative hessian computations.
			gradlog= gradvec(rbetamap,mask,condprobcat,meanth,mean2th,newtheta,probcat,theta(:,quad,obs),thpr)
			do row=1,nreduce
				if(gradlog(row)/=0.0_8)gradobs(row)=gradobs(row)+postdensity(quad,obs)*gradlog(row)
			end do
			if(nr)then
				nhessquad=nhessvec(rbetamap,mask,condprobcat,covth,cubth,newtheta,&
					probcat,quarth,theta(:,quad,obs),thpr)
				do row=1,nreduce
					do col=1,row
						if((gradlog(row)/=0.0.and.gradlog(col)/=0.0_8).or.nhessquad(row,col)/=0.0_8)&
							nhessobs(row,col)=nhessobs(row,col)+postdensity(quad,obs)*(nhessquad(row,col)&
							-gradlog(row)*gradlog(col))
					end do
				end do
			end if
		end do
		do row=1,nreduce
			if(gradobs(row)/=0.0_8)	gradobs(row)=gradobs(row)/prob(obs)
		end do
		postdensity(:,obs)=postdensity(:,obs)/prob(obs)
		if(nr)then 
			do row=1,nreduce
				do col=1,nreduce
					if(nhessobs(row,col)/=0.0_8.or.(gradobs(row)/=0.0_8.and.gradobs(col)/=0.0_8))&
						nhessobs(row,col)=nhessobs(row,col)/prob(obs)+gradobs(row)*gradobs(col)
				end do
			end do
		end if
		do row=1,nreduce
			if(gradobs(row)/=0.0_8)grad(row)=grad(row)+obsweight(obs)*gradobs(row)
		end do
		grads(:,obs)=gradobs
		if(nr)then
			do row=1,nreduce
				do col=1,row
					if(nhessobs(row,col)/=0.0_8)&
						nhess(row,col)=nhess(row,col)+obsweight(obs)*nhessobs(row,col)
				end do
			end do
		end if	
		
		

		loglik=loglik+obsweight(obs)*log(prob(obs))


	end do
!	Check for adequate progress.
	if(nr)call makesym(nhess)
	oldbeta=beta
	if(constdim>0)oldconsterr=consterr
	step=1.0_8
!	Save log likelihood.
	loglikn=loglik
!	Gradients and hessians.
	graddes=matmul(transpose(rdesign),grad)	

	if(constdim1<constdim) graddes=graddes+2.0_8*constweight*matmul(constmat(:,constdim1+1:constdim),consterr)
	nrskip=.false.
	if(nr)then
		nhessdes=matmul(transpose(rdesign),matmul(nhess,rdesign))
		nhessdesplus=nhessdes
		if(constdim1<constdim)nhessdesplus=nhessdesplus+constquad
		nhessdeschol=chol(nhessdesplus,tolsing)
!	If not nonnegative definite, then use Louis approach.
!		if(minval((/(nhessdeschol(row,row),row=1,size(gamma))/))<0.0_8)nrskip=.true.
        do ita=1,maxitb
            minchol=minval((/(nhessdeschol(row,row),row=1,size(gamma))/))
            maxchol=maxval((/(nhessdeschol(row,row),row=1,size(gamma))/))
            if(minchol<tolc*maxchol)then
                do row=1,size(gamma)
                    nhessdesplus(row,row)=nhessdesplus(row,row)+ita*maxchol
                end do
                nhessdeschol=chol(nhessdesplus,tolsing)
            else
                exit
            end if
        end do
    end if
	if(.not.nr)then
		nhess=information(.false.,grads,obsweight)
		nhessdes=matmul(transpose(rdesign),matmul(nhess,rdesign))
		nhessdesplus=nhessdes
		if(constdim1<constdim)nhessdesplus=nhessdesplus+constquad	
		nhessdeschol=chol(nhessdesplus,tolsing)
        do ita=1,maxitb
            minchol=minval((/(nhessdeschol(row,row),row=1,size(gamma))/))
            maxchol=maxval((/(nhessdeschol(row,row),row=1,size(gamma))/))
            if(minchol<tolc*maxchol)then
                do row=1,size(gamma)
                    nhessdesplus(row,row)=nhessdesplus(row,row)+ita*maxchol
                end do
                nhessdeschol=chol(nhessdesplus,tolsing)
            else
                exit
            end if
        end do
    end if

	
	
	if(finish)return
!	Tentative steps.

	gammastep=solve(nhessdeschol,graddes)
	if(constdim1>0)then
		do col=1,constdim1
			v(:,col)=solve(nhessdeschol,constmat(:,col))
		end do
		vv=matmul(transpose(v),matmul(nhessdesplus,v))
		vv=chol(vv,tolsing)
		gammastep=gammastep&
            +matmul(v,solve(vv,consterr(1:constdim1)&
            -matmul(transpose(constmat(:,1:constdim1)),gammastep)))

	end if
	betastep=matmul(design,gammastep)
	maxbeta=maxval(abs(betastep))
!	See if step is too large.
	if (maxbeta>kappa) step=kappa/maxbeta	
!	Predict improvement.	
	der=dot_product(graddes,gammastep)


	do ita=1,maxitb
!	Take step and check progress.
		beta=oldbeta+step*betastep
		locations=beta(1:ncat)
		scales=reshape(beta(nscale:nscale1),(/dimlatout,ncat/))
		
		call getlinquad(beta,lintheta,quadtheta)
		
!	Error flag
		error=.false.	


		loglik=0.0_8
	
		if(constdim>0) consterr=oldconsterr-step*matmul(transpose(constmat),gammastep)
		if(constdim1<constdim) loglik=-constweight*sum(consterr(constdim1+1:constdim)*consterr(constdim1+1:constdim))
			
		
		prob=0.0_8
!	New log likelihood.
		do obs=1,nobs
		
			if(obsweight(obs)<=0.0_8) cycle
					
! Observation obs.
			resp=dat(:,obs)
			
			mask=.false.
			do item=1,nitems
				if(resp(item)>=0.and.resp(item)<numcatobs(item))mask(item)=.true.
			end do
			
!	Density specification for observation.
			thpr=indvar(:,obs)	
			call getlinquadth(lintheta,quadtheta,thpr,linth,quadth)
				
			if(normal)then
!	Mean and covariance matrix of latent vector.
				nhchol=chol(-2.0_8*quadth,tolsing)
!	Check that covariance matrix is positive definite.
				if(minval((/(nhchol(row,row),row=1,dimlatin)/))<=0.0_8)then
					error=.TRUE.
					exit
				end if
				meanth=solve(nhchol,linth)

				
				scaletheta=sqrt(product((/(nhchol(row,row),row=1,dimlatin)/)))
				scaletheta=scaletheta*exp(-0.5_8*dot_product(linth,meanth))

				
			

! 	Location maximum of likelihood contribution for examinee k.
			

!	See if adaptive approach to be used.
				if(maxita>0) then

					call maxpost(catobsrange,maxita,maxita,npred,numcat,numcatobs,resp,mask,beta,changemin,&
						linth,lintran,maxdalpha,quadth,tau,tola,tolsing,&
						nhchol,alpha(:,obs))
					cholnhess(:,:,obs)=nhchol

				end if

				oldalpha=alpha(:,obs)



				do row=1,dimlatin
					scalez(row)=1.0_8/sqrt(nhchol(row,row))
				end do
				scalezp(obs)=product(scalez)
				scaletheta=scaletheta*scalezp(obs)
				do quad=1,nquad
					theta(:,quad,obs)=scaleadapt(nhchol,oldalpha,scalez,quadpoint(:,quad))
				end do
				priordensity=density(linth,theta(:,:,obs),quadth,quadweight,scaletheta)


				
					

!	 Scale factor
				
		

				

			else
				priordensity=densitym(linth,quadpoint,quadth,quadweight)
				
			end if				
!	Cycle through quadrature points.
			do quad=1,nquad
				

				

! Probability computation.
				newtheta=matmul(lintran,theta(:,quad,obs))
				probcat=probvec(numcat,mask,locations,scales,newtheta)
				probcatobs=probvecobs(catobsrange,numcatobs,resp,mask,probcat)
				postdensity(quad,obs)=product(probcatobs,mask)*priordensity(quad)

				prob(obs)=prob(obs)+postdensity(quad,obs)
				
			end do

			postdensity(:,obs)=postdensity(:,obs)/prob(obs)
			loglik=loglik+obsweight(obs)*log(prob(obs))
		end do
!	Error in covariance matrix requires step reduction.

		if(error)then
			loglik=loglikn
			step=step/2.0_8
		else



!	See if change in log likelihood acceptable.

			change=loglik-loglikn


			if(change>step*der*changemin) exit
			stepmin = step*tau
			step = step*der/(2.0_8*(der-change/step))
			if(step<stepmin) step = stepmin
		end if
	end do
	if(printprog.or.printprogstd)call  iterationreport(comment,it,unitno,printprog,printprogstd,loglik,step)
	gamma=gamma+step*gammastep
!	Finish up.
	if(change<(-tol*loglikn))finish=.true.
!	The following check arises due to quadrature errors.
	if(it>0)then
		if(loglik<loglikold)finish=.true.
	end if
	loglikold=loglik

end do
return
end subroutine  ndmmaxlik

