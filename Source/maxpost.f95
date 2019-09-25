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
interface

	function chol(sym,tolsing)
!	chol is used to compute a modified Cholesky decomposition of the n by n matrix sym.
!	The singularity tolerance is tolsing.
		implicit none
		real(kind=8),intent(in)::sym(:,:),tolsing
	 	real(kind=8)::chol(size(sym,1),size(sym,1))
	end function chol

!	Obtain conditional covariance matrix of item scale parameters for an observation.
!	catobsrange contains the range of underlying categories per observed category.
!	numcatobs is a vector with the number of observed  categories per item.
!	resp is the observation vector.
!	The items presented are indicated by mask.
!	condmeans contains the conditional means.
!	The underlying item conditional probabilities are in condprobcat.
!	The scale parameters are in scales.
	function condcovmat(catobsrange,numcatobs,resp,mask,condmeans,condprobcat,scales)
		implicit none
		integer,intent(in)::numcatobs(:),catobsrange(:,:),resp(:)
		logical,intent(in)::mask(:)
		real(kind=8),intent(in)::condmeans(:,:),condprobcat(:),scales(:,:)
		real(kind=8)::condcovmat(size(scales,1),size(scales,1),size(resp))!	Obtain conditional covariance matrix of item scale parameters for an observation resp.

	end function condcovmat
!	Obtain conditional mean item scale parameters.
!	catobsrange contains the range of underlying categories per observed category.
!	numcatobs is a vector with the number of observed categories per item.
!	The observation is resp.
!	The items presented are indicated by mask.
!	The underlying conditional item probabilities are in condprobcat.
!	scales is the matrix of scale parameters. 
	function condmean(catobsrange,numcatobs,resp,mask,condprobcat,scales)
		implicit none
		integer,intent(in)::catobsrange(:,:),numcatobs(:),resp(:)
		logical,intent(in)::mask(:)
		real(kind=8),intent(in)::condprobcat(:),scales(:,:)
		real(kind=8)::condmean(size(scales,1),size(resp))
	end function condmean	

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
	
!	Obtain covariance matrix of item scale parameters for an observation.
!	numcat is a vector with the number of underlying categories per category.
!	The items presented are indicated by mask.
!	means contains the corresponding means.
!	The underlying item probabilities are in probcat.
!	scales contains the scale parameters.   
	function covmat(numcat,mask,means,probcat,scales)
		implicit none
		integer,intent(in)::numcat(:)
		logical,intent(in)::mask(:)
		real(kind=8),intent(in)::means(:,:),probcat(:),scales(:,:)
		real(kind=8)::covmat(size(scales,1),size(scales,1),size(mask))
	end function covmat
	
!	
	
!	Find the gradient of the conditional log posterior.
!	mask indicates which items are present.
!	condmeans is the vector of conditional means of scale factors.
!	linth is the linear component of the log density.
!	lintran is the linear transformation of the latent vector.
!	means is the vector of means of scale factors.
!	quadth is the quadratic component of the log density.
!	theta is the value of the latent vector.
	function gradtheta(mask,condmeans,linth,lintran,means,quadth,theta)
		implicit none
		logical,intent(in)::mask(:)
		real(kind=8),intent(in)::condmeans(:,:),linth(:),lintran(:,:),means(:,:),quadth(:,:),theta(:)
		real(kind=8)::gradtheta(size(theta))
	end function gradtheta
	

!	Compute inverse for n by n matrix with modified Cholesky
!	decomposition tri.
	function invert(tri)
		implicit none
		real(kind=8),intent(in)::tri(:,:)
		real(kind=8)::invert(size(tri,1),size(tri,1))
	end function invert

!	Obtain mean item scale parameters for an observation.
!	numcat is a vector with the number of
!		underlying categories per item.
!	The items presented are indicated by mask.
!	The underlying item probabilities are in probcat.
!	scales contains the scale parameters.
	function mean(numcat,mask,probcat,scales)
		implicit none
		integer,intent(in)::numcat(:)
		logical,intent(in)::mask(:)
		real(kind=8),intent(in)::probcat(:),scales(:,:)
		real(kind=8)::mean(size(scales,1),size(mask))	
	end function mean

!	Find the negative hessian of the conditional log posterior of an observation
!	relative to the latent vector theta.
!	mask indicates which items are present.
!	condcovs is the arrayr of conditional covariance matrics of scale factors.
!	covs is the array of covariance matrices of scale factors.
!	lintran is the linear transformation of the latent vector.
!	quadth is the matrix parameter for the log density. 

	function nhesstheta(mask,condcovs,covs,lintran,quadth)
		implicit none
		logical,intent(in)::mask(:)
		real(kind=8),intent(in)::condcovs(:,:,:),covs(:,:,:),lintran(:,:),quadth(:,:)
		real(kind=8)::nhesstheta(size(quadth,1),size(quadth,1))
	end function nhesstheta

!	Find a substitute negative hessian nhessth of the conditional log posterior of an observation
!	relative to the latent vector theta.
!	mask indicates which items are present.
!	condmeans is the vector of conditional means of scale factors.
!	gradth is the gradient.
!	lintran is the linear transformation of the latent vector.
!	means is the vector of means of scale factors.
!	quadth is the matrix of parameters used in the log density.
	function nhessthsub(mask,condmeans,gradth,lintran,means,quadth)
		implicit none
		logical,intent(in)::mask(:)
		real(kind=8),intent(in)::condmeans(:,:),gradth(:),lintran(:,:),means(:,:),quadth(:,:)
		real(kind=8)::nhessthsub(size(means,1),size(means,1))
	end function nhessthsub
	
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


!		This program uses the modified Cholesky decomposition
!		in lu to solve the equation lusolve = b.
	function solve(lu,b)
		implicit none
		
		real(kind=8),intent(in)::lu(:,:)
		real(kind=8),intent(in)::b(:)
		real(kind=8)::solve(size(b))
	end function solve
end interface
integer,intent(in)::catobsrange(:,:),maxita,maxitb,npred,numcat(:),numcatobs(:),resp(:)
logical,intent(in)::obsmask(:)
real(kind=8),intent(in)::beta(:),changemin,linth(:),lintran(:,:),maxdtheta,quadth(:,:),tau,tol,tolsing
real(kind=8),intent(inout)::theta(:)
real(kind=8),intent(out)::nhchol(:,:)

!	dimlatin is the dimension of the latent vector.
!	dimlatout is the dimension of the transformed latent vector.
!	ita counts main iterations.
!	itb counts subiterations.
!	item is an item counter.
!	ncat is the total number of underlying categories.
!	nitems is the number of items.

!	nscale is the pointer for scale parameters.
!	nscale1 is the pointer for the last scale parameter.
!	row is a row counter.

integer::dimlatin,dimlatout,ita,itb,item,ncat,nitems,&
	nscale,nscale1,row
!	Data mask.
logical::mask(size(obsmask))
!	change is a change in logcondprob.
!	condcovs is used to find conditional covariance  matrices of scale parameters.
!	condmeans is used to find conditional means of scale parameters.
!	condprobcat is used for conditional probabilities.
!	covs  is used to find covariance matrices of scale parameters.
!	covth is the covariance matrix of the latent vector.
!	der is a directional derivative.
!	gradth is the gradient of the logarithm of the posterior.
!	lintheta is for linear parameters for the log density.
!	locations is for location parameters.
!	logcondprob is a log conditional probability plus log density.
!	maxthetastep is the maximum norm of thetastep.
!	means is used to find means of scale parameters.
!	meanth is the mean of the latent vector.
!	newtheta is the transformed value of the latent vector.
!	nhchol is a Cholesky decomposition of -2*quadth.
!	nhessth is the negative heessian of the logarithm of the posterior.
!	oldlogcondprob is the previous value of logcondprob.
!	oldtheta is a previous approximation of the location of the maximum.
!	probcat is used for underlying marginal conditional probabilities.
!	probcatobs is used for observed marginal conditional probabilities.
!	quadtheta is for quadratic parameters for the log density.
!	scales is for scale parameters.
!	step is a step size.
!	stepmin is a minimum step size.
!	thetastep is the proposed change in the approximation of the location of the maximum.

real(kind=8)::change,condcovs(size(lintran,1),size(lintran,1),size(obsmask)),&
	condmeans(size(lintran,1),size(obsmask)),condprobcat(sum(numcat)),&
	covs(size(lintran,1),size(lintran,1),size(obsmask)),&
	der,gradth(size(theta)),&
	locations(sum(numcat)),logcondprob,&
	maxthetastep,means(size(lintran,1),size(obsmask)),&
	newtheta(size(lintran,1)),nhessth(size(theta),size(theta)),&
	oldlogcondprob,oldtheta(size(theta)),&
	probcat(sum(numcat)),probcatobs(size(obsmask)),scales(size(lintran,1),sum(numcat)),&
	step,stepmin,thetastep(size(theta))
	
	

!	Set up parameter arrays for simplified processing.
ncat=sum(numcat)
nitems=size(obsmask)

dimlatin=size(theta)
dimlatout=size(lintran,1)


!	Item scale parameters.
nscale=ncat+1
nscale1=nscale+size(scales)-1
locations=beta(1:ncat)
scales=reshape(beta(nscale:nscale1),(/dimlatout,ncat/))






!	Mask
mask=.false.
do item=1,nitems
	if(resp(item)>=0.and.resp(item)<numcatobs(item))mask(item)=.true.
end do
mask=mask.and.obsmask

! 	Location maximum of likelihood contribution for examinee k.


! 	Probabilities.

!	Individual probabilities
newtheta=matmul(lintran,theta)
probcat=probvec(numcat,mask,locations,scales,newtheta)
probcatobs=probvecobs(catobsrange,numcatobs,resp,mask,probcat)
logcondprob=sum(log(probcatobs),mask)+dot_product(linth+matmul(quadth,theta),theta)


do ita=1,maxita

!	Needed  means and covariance matrices.
	condprobcat=condprobvec(catobsrange,numcatobs,resp,mask,probcat,probcatobs)

	means=mean(numcat,mask,probcat,scales)
	covs=covmat(numcat,mask,means,probcat,scales)

	condmeans=condmean(catobsrange,numcatobs,resp,mask,condprobcat,scales)


	condcovs=condcovmat(catobsrange,numcatobs,resp,mask,condmeans,condprobcat,scales)

!	Gradient and negative hessian.
	gradth=gradtheta(mask,condmeans,linth,lintran,means,quadth,theta)
	nhessth=nhesstheta(mask,condcovs,covs,lintran,quadth)

! 	Step computation.
	nhchol=chol(nhessth,tolsing)


!	If not nonnegative definite, then use Louis approach.
	if(minval((/(nhchol(row,row),row=1,dimlatin)/))<0.0_8)then
		nhessth=nhessthsub(mask,condmeans,gradth,lintran,means,quadth)

		nhchol=chol(nhessth,tolsing)

	end if
	thetastep=solve(nhchol,gradth)

			
! 	Step reduction if necessary.
	maxthetastep=maxval(abs(thetastep))
	if(maxthetastep==0.0_8)exit
	step=1.0_8
	if(maxthetastep>maxdtheta)step=maxdtheta/maxthetastep
	der=dot_product(gradth,thetastep)
	oldtheta=theta
	oldlogcondprob=logcondprob
	do itb=1,maxita
		theta=oldtheta+step*thetastep

! 	Probabilities.

		newtheta=matmul(lintran,theta)
		probcat=probvec(numcat,mask,locations,scales,newtheta)
		probcatobs=probvecobs(catobsrange,numcatobs,resp,mask,probcat)
		logcondprob=sum(log(probcatobs),mask)+dot_product(linth+matmul(quadth,theta),theta)


!	Test  for adequate progress.

		change=logcondprob-oldlogcondprob
		if(change>step*der*changemin) exit
		stepmin = step*tau
		step = step*der/(2.0_8*(der-change/step))
		if(step<stepmin) step = stepmin
	end do
	if(change<tol) exit
end do
!	Save results.
condprobcat=condprobvec(catobsrange,numcatobs,resp,mask,probcat,probcatobs)
means=mean(numcat,mask,probcat,scales)
covs=covmat(numcat,mask,means,probcat,scales)
condmeans=condmean(catobsrange,numcatobs,resp,mask,condprobcat,scales)
condcovs=condcovmat(catobsrange,numcatobs,resp,mask,condmeans,condprobcat,scales)
!	Gradient and negative hessian.
gradth=gradtheta(mask,condmeans,linth,lintran,means,quadth,theta)
nhessth=nhesstheta(mask,condcovs,covs,lintran,quadth)


nhchol=chol(nhessth,tolsing)
!	If not nonnegative definite, then use Louis approach.
if(minval((/(nhchol(row,row),row=1,dimlatin)/))<0.0_8)then
	nhessth=nhessthsub(mask,condmeans,gradth,lintran,means,quadth)
	nhchol=chol(nhessth,tolsing)
end if
return
end subroutine maxpost
