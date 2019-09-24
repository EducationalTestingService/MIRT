!	Basic maximum likelihood for item response models with multivariate normal latent vectors.
!	Responses are polytomous and are assumed to be integers.  Predictors of latent variables
!	are floating point.  Weights may be employed for observations.  Responses are conditionally
!	independent given latent variables.  Predictors and responses and conditionally independent
!	given latent vectors.  Responses given  latent vectors satisfy multinomial logit models that
!	have logarithms of probabilities quadratic functions of latent vectors.
!	Responses that are observed may contain more than one possible underlying response.
!	There is an assumption that considerable variation may exist in input and output procedures
!	for different applications.  As a consequence, input and output are found in external subroutines.   
program irtnorm
use unitdef
implicit none
interface
!
!
!	Find the Akaike correction.
!	modeldim is the model dimension.
!	loglik is the log likelilhood.
!	totalitemsis the total number of presented items.
	real(kind=8) function akaike(modeldim,loglik,totalitems)
		implicit none
		integer,intent(in)::modeldim
		real(kind=8),intent(in)::loglik,totalitems
	end function akaike
!
!	Get beta names.
!	Input:
!
!	factorname is the array of factor names.
!	itemname is the arrray of item names.
!	predname is the array of predictor names.
!	skillname is the array of skill names.
!	numcat is the number of underlying categories per item.
!   beta is the beta vector.
!
!	Output:
!
!	betanames is the array of beta names.

!
	function betanames(factorname,itemname,predname,skillname,numcat,beta)
		implicit none

		character(len=32),intent(in)::itemname(:),skillname(:),factorname(:),predname(:)
		integer,intent(in)::numcat(:)
		real(kind=8)::beta(:)


		character(len=64)::betanames(size(beta))
	end function betanames
!	beta output
!	betaname is the vector of betar names.
!	unitbeta is the unit number.
!	complex indicates complex sampling.
!	beta is the vector of estimates.
!	eacovbeta is the array of estimated covariances of beta estimates under the model.
!	eacovbeta_louis is the array of estimated covariances of beta estimates under the model with the Louis approach.
!	eacovbeta_sandwich is the array of estimated covariances of beta estimates without assuming the model.
!	eacovbeta_complex is the array of estimated covariances of beta estimates under complex sampling.

	subroutine betaoutput(betaname,unitbeta,complx,beta,eacovbeta,&
        eacovbeta_complex,eacovbeta_louis,eacovbeta_sandwich)
		implicit none
		character(len=64)::betaname(:)
		integer,intent(in)::unitbeta
		logical,intent(in)::complx
		real(kind=8),intent(in)::beta(:),eacovbeta(:,:),&
            eacovbeta_complex(:,:),eacovbeta_louis(:,:),eacovbeta_sandwich(:,:)
	end subroutine betaoutput
!	buffers for output.
!	This function computes the modified Cholesky decomposition
!	of the n by n nonnegative definite symmecholc macholx sym.
!	tolsing is the singularity tolerance.	
	function chol(sym,tolsing)
		implicit none
		real(kind=8),intent(in)::sym(:,:),tolsing
	 	real(kind=8)::chol(size(sym,1),size(sym,1))
	end function chol

!   See if unit to be closed.
    subroutine closeunit(unit)

        use unitdef
        implicit none
        integer,intent(in)::unit
    end subroutine closeunit

!	Get coordinates for quadrature points and weights.
!	dimpoints contains number of points per dimension.
!	coord is the coordinate vectors.
!	cweight is used to modify quadrature weights.

	subroutine coordinates(dimpoints,coord,cweight)
		implicit none
		integer,intent(in)::dimpoints(:)
		integer,intent(out)::coord(:,:)
		real(kind=8),intent(out)::cweight(:)
	end subroutine coordinates
!	Design specifications.
!	Input:

!	intdim provides the number of columns in itemmatrix for intercepts associated with each item. 
!	slopedim provides the number of columns in itemmatrix for slopes associated with an item.
!	constguess is for a constant guessing parameter 
!	fixdiag is for fixed diagonals for the constant
!	fixedguess is for fixed guessing parameters.
!	guess is for guessing parameters.
!	predlinmap is the mapping of predictors to factors for the linear term.
!	predquadmap is the mapping of predictors to factors for the quadratic term.
!	rasch is for a Rasch model for a skill
!	raschslope1 is for a Rasch model for a skill with slope 1.
!
!	Output:
!
!	constdim is the number of linear constraints.
!	constdim1 is the number of linear constraints not treated via least squares.
!	dimdesign is the number of independent parameters.
!	proport is true if the constraint sum of squares is proportional 
!		to sum of observation weights.
!	specialtrans is true if a special transition matrix is
!	used.
	subroutine designspecs(intdim,slopedim,constguess,&
		fixdiag,fixedguess,guess,&
		predlinmap,predquadmap,rasch,raschslope1,constdim,constdim1,dimdesign,proport,specialtrans)
		implicit none
		integer,intent(in)::intdim(:),slopedim(:,:)

		integer,intent(out)::constdim,constdim1,dimdesign
		logical,intent(in)::constguess,fixdiag(:),fixedguess(:),guess(:),&
			predlinmap(:,:),predquadmap(:,:,:),rasch(:),raschslope1(:)
		logical,intent(out)::proport,specialtrans	
	end subroutine designspecs



!	Find estimated conditional expectations and conditional covariance matrices for latent vectors.
!	postdensity is the densities.
!	theta is the latent vectors.
!	covtheta is the array of conditional covariance matrices.
!	meantheta is the array of conditional means. 
	subroutine eap(postdensity,theta,covtheta,meantheta)
		implicit none
		real(kind=8),intent(in)::postdensity(:,:),theta(:,:,:)
		real(kind=8),intent(out)::covtheta(:,:,:),meantheta(:,:)
	end subroutine eap
	
!	print eap data for latent vector.
!	comment is used for identification.
!	factorname gives names of elements of the latent vector.
!	id is individual id.
!	unit uniteap is the unit for output.
!	eapid indicates if ids are put out.
!	covtheta are conditional covariance matrices.
!	meantheta are eaps. 

	subroutine eapoutput(comment,factorname,id,uniteap,eapid,covtheta,meantheta)
		implicit none
		character(len=*),intent(in)::comment
		character(len=32),intent(in)::factorname(:)
		character(len=32),intent(in),optional::id(:)
		integer,intent(in)::uniteap
		logical,intent(in)::eapid
		real(kind=8),intent(in)::covtheta(:,:,:),meantheta(:,:)
	end subroutine eapoutput
!	Find estimated conditional expectations and conditional covariance matrices
!   for transformations of integer-weighted sums.
!	catobsrange maps underlying to observed categories.
!	Responses are in dat.
!	maxscore is the maximum score.
!	maxw provides maximum item weights.
!	minscore is the minimum score.
!	minw provides minimum item weights.
!	numberscales is the number of scale scores.
!	numcat provides ranges for underlying observations.
!	numcatobs provides ranges for observations.
!	weight is the weight.

!	beta is the parameter vector.
!	predictors are in indvar.
!	lintran is the linear transformation for the latent vector.
!	postdensity contains posterior weights.
!	scale is the matrix of scale transformations.
!	theta contains posterior quadrature points.

!	covscale is the conditional covariance of the scale transformation of the weighted sum.
!	meanscale is the conditional mean of the scale transformation of the weighted sum.


	subroutine eapscale(catobsrange,maxscore,maxw,minscore,&
		minw,numberscales,numcat,numcatobs,weight,beta,&
		lintran,postdensity,scale,theta,covscale,meanscale)
		implicit none
		integer,intent(in)::catobsrange(:,:),maxscore,maxw(:),minscore,minw(:),numberscales,&
			numcat(:),numcatobs(:),weight(:)

		real(kind=8),intent(in)::beta(:),lintran(:,:),postdensity(:,:),&
    scale(numberscales,minscore:maxscore),theta(:,:,:)
		real(kind=8),intent(out)::covscale(:,:,:),meanscale(:,:)
	end subroutine eapscale

!
!	Find estimated conditional expectations and conditional covariance matrices
!   for transformations of latent vector.
!	m is the dimension of the transformation.
!	postdensity contains posterior weights.
!	theta contains posterior quadrature points.

!	covtran is the conditional covariance of the transformation of the latent vector.
!	meantran is the conditional mean of the transformation of the latent vector.


	subroutine eaptrans(m,postdensity,theta,covtran,meantran)
		implicit none
		integer,intent(in)::m
		real(kind=8),intent(in)::postdensity(:,:),theta(:,:,:)
		real(kind=8),intent(out)::covtran(:,:,:),meantran(:,:)
	end subroutine eaptrans

!	Find estimated conditional expectations and conditional covariance matrices
!   for transformed latent vectors.
!	covtheta is the array of conditional covariance matrices.
!	lintran is the linear transformation for the latent vector.
!	meantheta is the array of conditional means.
!	covnewtheta is the array of transformed conditional covariances.
!	meannewtheta is the  array of transformed conditional means. 
	subroutine eapskill(covtheta,lintran,meantheta,covnewtheta,meannewtheta)
		implicit none
		real(kind=8),intent(in)::covtheta(:,:,:),lintran(:,:),meantheta(:,:)
		real(kind=8),intent(out)::covnewtheta(:,:,:),meannewtheta(:,:)	
	end subroutine eapskill
	
!	Find estimated conditional expectations and conditional covariance matrices
!   for weighted sums.
!	Responses are in dat.
!	numcat provides ranges for underlying observations.
!	numcatobs provides ranges for observations.
!	mask indicates items to use for computations.
!	beta is the parameter vector.
!	predictors are in indvar,
!	lintran is the linear transformation for the latent vector.
!	postdensity contains posterior weights.
!	theta contains posterior quadrature points.
!	wtsum is a matrix used to compute linear combinations of the observations.
!	covsum is the conditional covariance of wtsum.
!	meansum is the conditional mean of wtsum.
	subroutine eapwtsum(dat,numcat,numcatobs,mask,beta,lintran,postdensity,&
        theta,wtsum,covsum,meansum)
		implicit none
		integer,intent(in)::dat(:,:),numcat(:),numcatobs(:)
		logical,intent(in)::mask(:)
		real(kind=8),intent(in)::beta(:),lintran(:,:),postdensity(:,:),theta(:,:,:),&
			wtsum(:,:)
		real(kind=8),intent(out)::covsum(:,:,:),meansum(:,:)
	end subroutine eapwtsum


!	Information output
!   comment is a comment.
!   unitinfo is the unit number.
!	modeldim is the rank of the design matrix.
!	complx indicates complex sampling.
!	akent is corresponding Akaike estimate.
!	ent is entropy per presented item.
!	ghent is coresponding Gilula-Habermann estimate.
!	ghent_complex is corresponding Gilula-Haberman estimate for complex sampling.
!	loglik is log likelihood.
!	sdentropy is estimated asymptotic standard deviation of ent.
!	sdentropy is asymptotic standard error of ent for complex sampling.
	subroutine entoutput(comment,modeldim,unitinfo,complx,akent,&
        ent,ghent,ghent_complex,loglik,sdentropy,sdentropy_complex)
		implicit none
        character(len=*),intent(in)::comment
		integer,intent(in)::unitinfo
		integer,intent(in),optional::modeldim
		logical,intent(in)::complx
		real(kind=8),intent(in),optional::akent,ent,ghent,ghent_complex,&
            loglik,sdentropy,sdentropy_complex
	end subroutine entoutput
	
!	gamma output
!	paramname is the vector of parameter names.
!	unitparam is the unit number.
!	complex indicates complex sampling.
!	eacovgam is the array of estimated covariances of gamma estimates under the model.
!	eacovgam_louis is the array of estimated covariances of gamma estimates
!   under the model with the Louis approach.
!	eacovgam_sandwich is the array of estimated covariances of gamma estimates
!   without assuming the model.
!	eacovgam_complex is the array of estimated covariances of gamma estimates under complex sampling.
!	gamma is the vector of estimates.
	subroutine gammaoutput(paramname,unitparam,complx,eacovgam,eacovgam_complex,eacovgam_louis,eacovgam_sandwich,gamma)
		implicit none
		integer,intent(in)::unitparam
		logical,intent(in)::complx
		character(len=64)::paramname(:)
		real(kind=8),intent(in)::eacovgam(:,:),eacovgam_complex(:,:),eacovgam_louis(:,:),eacovgam_sandwich(:,:),gamma(:)
	end subroutine gammaoutput

!       Initialize alpha and cholnhess
!       alpha is location of the maximum of the conditional probability
!		of a response given the latent vector.
!       cholnhess is the modified Cholesky decomposition of the
!       negative hessian of the log conditional probability with respect
!       to the latent vector at the maximum.
	subroutine getalpha(alpha,cholnhess)
		implicit none
		real(kind=8),intent(out)::alpha(:,:),cholnhess(:,:,:)
	end subroutine getalpha


!	Constraint specification.
!	constmat is the matrix .
!	constvec is the vector.
!	constmat(gamma) is to equal constvec.

	subroutine getconstraint(constmat,constvec)
		implicit none
		real(kind=8),intent(out)::constmat(:,:),constvec(:)
	end subroutine getconstraint		

!	Obtain data.
!	Input:
!
!	fileformat is the file format.
!	nexternal is the number of external variables.
!	nitems is the number of items.
!	nobs is the number of observations to be read.
!	npred is the number of predictors, including the constant predictor 1.
!	numcodes is the number of recodes per item.
!	numdistcodes is the number of distractor recodes per item.
!	recodetab is the table of recodes for each item.
!	recodedisttab is the table of distractor recodes for each item.
!	distract is .true. if distractors are read.
!	recode indicates if recodes are present.
!	recodedist indicates if distractor recodes are present.
!	stratify is .true. if sampling is stratified.
!	If true, useid causes reading of an examinee identification code.
!	usepsu is .true. is primary sampling units are considered.
!	weight is .true. if observation weights are used.
!
!	Output:
!
!	id stores examinee identifications.
!	dat is the response matrix.
!	distdat is the distractor matrix.
!	psu is the vector of psus for observations.
!	stratum is the vector of strata for observations.
!	extvar stores external variables.
!	indvar is the matrix of independent variables.
!	obsweight is the vector of observation weights.
	subroutine getdata(fileformat,nexternal,nitems,nobs,npred,&
		numcodes,numdistcodes,recodetab,recodedisttab,distract,&
        recode,recodedist,stratify,useid,usepsu,weight,&
		id,dat,distdat,psu,stratum,extvar,indvar,obsweight)
		implicit none
		character(len=256),intent(in)::fileformat
		integer,intent(in)::nexternal,nitems,nobs,npred
		integer,intent(in),optional::numcodes(:),numdistcodes(:),&
        recodetab(:,:),recodedisttab(:,:)
		logical,intent(in)::distract,recode,recodedist,stratify,useid,usepsu,weight

		character(len=32),intent(out),optional::id(:)
		integer,intent(out)::dat(:,:),distdat(:,:)
		integer,intent(out),optional::psu(:),stratum(:)
		real(kind=8),intent(out),optional::extvar(:,:)
		real(kind=8),intent(out)::indvar(:,:),obsweight(:)
	end subroutine getdata

!	Basic input specifications.
!	fileformat is the Fortran format for reading the data.
!	filename is the file name.
!	nexternal is the number of external variables.
!	nitems is the number of items.
!	nobs is the number of observations,nitems is the number of items.
!	nobsstart is the initial sample size for starting calculations.
!	npred is the number of predictors.
!	complx is stratify.or.usepsu.or.weight.
!	distract is .true. if distractors are read.
!	recode is used for item recoding.
!	recodedist indicates if distractor recodes are present.
!	stratify is the indicator of whether sampling is stratified.
!	If true, useid causes reading of an examinee identification code.
!	usepsu is the indicator of whether primary sampling units are involved.
!	weight indicates that weighting is used for observations.
	subroutine getdataspec(fileformat,filename,nexternal,nitems,nobs,nobsstart,npred,&
		complx,distract,recode,recodedist,stratify,useid,usepsu,weight)
		implicit none

		character(len=256),intent(out)::fileformat
        character(len=256),intent(out)::filename
		integer,intent(out)::nexternal,nitems,nobs,nobsstart,npred
		logical,intent(out)::complx,distract,recode,recodedist,stratify,useid,usepsu,weight
	end subroutine getdataspec

!	Get design settings.
!	Input:
!
!	factorname is the array of factor names.
!	itemname is the arrray of item names.
!	predname is the array of predictor names.
!	skillname is the array of skill names.
!	choices indicates the number of options for an item.
!	intdim gives the number intercept components per item.
!	numcat is the number of underlying categories per item.
!	slopedim gives the number of slope components per dimension and per item.
!	constguess is the indicator for a  constant guessing parameter.
!	fixdiag is the indicator for a fixed diagonal for the first predictor.

!	fixedguess is the indicator for a fixed guessing parameter.
!	guess is the indicator for guessing.
!	predlinmap is the mapping of predictors to factors for the linear term.
!	predquadmap is the mapping of predictors to factors for the quadratic term.

!	rasch is the indicator for a Rasch model.
!	raschslope1 is the indicator for a Rasch model with a slope of 1.
!	specialtrans is the indicator for a custom transition matrix.
!	itemmatrix is the item matrix.
!	setguess indicates the fixed value of the guessing parameter.
!	setslope specifies initialization of slope parameters.
!
!	Output:
!
!	paramname is the array of parameter names.
!	design is the final design matrix.
!	gamma is the parameter vector.
!	offset is the final offset.
!
	subroutine getdesign(factorname,itemname,predname,skillname,choices,intdim,&
		numcat,slopedim,&
		constguess,fixdiag,fixedguess,guess,predlinmap,predquadmap,rasch,raschslope1,&
		specialtrans,itemmatrix,setguess,setslope,&
		paramname,design,gamma,offset)
		implicit none

		character(len=32),intent(in)::itemname(:),skillname(:),factorname(:),predname(:)
		integer,intent(in)::choices(:),intdim(:),numcat(:),slopedim(:,:)
		logical,intent(in)::constguess,fixdiag(:),fixedguess(:),guess(:),&
			predlinmap(:,:),predquadmap(:,:,:),rasch(:),raschslope1(:),specialtrans
		real(kind=8),intent(in)::itemmatrix(:,:),setguess(:),setslope(:,:)
		character(len=64),intent(out)::paramname(:)
		real(kind=8),intent(out)::design(:,:),gamma(:),offset(:)
	end subroutine getdesign

!	Define the latent vector dimension dimlatin
!   and the transformed latent vector dimension dimlatout.
!	custom is for a special linear transformation.
	subroutine getdim(dimlatin,dimlatout,custom)
		implicit none
		integer,intent(out)::dimlatin,dimlatout
		logical,intent(out)::custom
	end subroutine getdim

!	Obtain distractor recodes.
!	numdistcodes gives number of distractor recodes per item.
!	recodedisttab gives distractor recodes.
	subroutine getdistrecodes(numdistcodes,recodedisttab)
		implicit none
		integer,intent(in)::numdistcodes(:)
		integer,intent(out)::recodedisttab(:,:)
	end subroutine getdistrecodes


!	Obtain factor specifications.
!	normal is .true. if normal distribution and .false. otherwise.
!	factorname contains factor names.
!	fixdiag indicates that diagonal terms of quadrature component
!   of log density are fixed.
!	fixquad indicates that quadrature component of log density is fixed.
!	independentf indicates that elements of latent vector are independent
!   of previous elements.
!	no_lin indicates that the linear component of the log density has an intercept of 0.
!	noquad indicates that the quadratic componet of the log density is 0.	
!	predlinmap is map of predictors to factors for linear term.
!	predquadmap is map of predictors to factors for quadratic term.
	subroutine getfactorspecs(normal,factorname,fixdiag,fixquad,&
        independentf,no_lin,noquad,predlinmap,predquadmap)
		implicit none
		character(len=32),intent(out)::factorname(:)
		logical,intent(in)::normal
        logical,intent(out)::fixdiag(:),fixquad,independentf(:),&
        no_lin(:),noquad,predlinmap(:,:),predquadmap(:,:,:)
	end subroutine getfactorspecs
!   Obtain conditional probabilities for distractors.
!   choices provides the number of choices for each item.
!   data provides item data.
!   distdat provides distractor data.
!   numcatobs provides the number of categories for item scores.
!   obsweight gives observation weights.
!   distfreq provides distractor frequencies.
!   distprob provides conditional distractor probabilities.
!   scorefreq provides item score frequencies.
    subroutine getfreq(choices,dat,distdat,numcatobs,obsweight,&
        distmap,distfreq,distprob,scorefreq)
        implicit none
        integer,intent(in)::choices(:),dat(:,:),distdat(:,:),numcatobs(:)
        real(kind=8),intent(in)::obsweight(:)
        integer,intent(out)::distmap(:)
        real(kind=8),intent(out)::distfreq(:),distprob(:),scorefreq(:)
    end subroutine getfreq


!	Initialize gamma.
	subroutine getgamma(gamma)
		implicit none
		real(kind=8),intent(inout)::gamma(:)
	end subroutine getgamma

!	Obtain basic item data.
!	itemname contains item names.
!	choices indicates the number of choices in a right-scored multiple-choice item.
!	intdim  is the dimension of the intercept matrix.
!	numcat is the number of underlying categories per item.
!	numcatobs is the number of observed categories per item.
!	skillnumbers specifies skills that correspond to items in a between-item case.
!	slopedim is the dimension of the slope matrix.
!	betweenitem indicates a between-item model.
!	catmap is for category maps.
!	constguess indicates a constant guessing parameter.
!	fixedguess indicates a fixed guessing parameter.
!	guess indicates a guessing parameter.
!	nominal indicates that a nominal model is used for the item.
!	specialitemint is used to specify items that need special treatment for intercepts.
!	specialitemslope is used to specify items that need special treatment for slopes.
!	setguess is the fixed setting of a guessing parameter.

	subroutine getitemdata(itemname,choices,intdim,numcat,numcatobs,&
        skillnumbers,slopedim,&
		betweenitem,catmap,constguess,fixedguess,guess,nominal,&
        specialitemint,specialitemslope,setguess)
		implicit none
		character(len=32),intent(out)::itemname(:)
		integer,intent(out)::choices(:),intdim(:),numcat(:),numcatobs(:),&
    skillnumbers(:),slopedim(:,:)
		logical,intent(out)::betweenitem(:),catmap(:),constguess,fixedguess(:),guess(:),&
			nominal(:),specialitemint(:),specialitemslope(:)
		real(kind=8),intent(out)::setguess(:)
	end subroutine getitemdata
	


!	Get item details such as category scores and maps of observed
!	categories to underlying categories.
!	intdim provides intercept dimensions for the item matrix.
!	numcat is the category counts for observed categories.
!	numcatobs is the category counts for underlying categories.
!	slopedim provides slope dimensions for the item matrix.
!	catmap is for category maps.
!	guess specifies guessing.
!	nominal is the indicator for a nominal model.
!	specialitemint specifies special treatment of an item intercept.
!	specialitemslope specifies special treatment of an item slope.
!	catobsrange provides the category map.
!	itemmatrix is the item matrix.
!	setslope helps initialize slope parameters.
	subroutine getitemdetail(intdim,numcat,numcatobs,slopedim,&
		catmap,guess,nominal,specialitemint,&
		specialitemslope,catobsrange,itemmatrix,setslope)
		implicit none
		integer,intent(in)::intdim(:),numcat(:),numcatobs(:),slopedim(:,:)
		integer,intent(out)::catobsrange(:,:)
		logical,intent(in)::catmap(:),guess(:),nominal(:),specialitemint(:),specialitemslope(:)
		real(kind=8),intent(out)::itemmatrix(:,:),setslope(:,:)
	end subroutine getitemdetail
!
!	Find linear transformation of latent vector.
!	custom specifies a transformation other than the default.
	subroutine getlintran(custom,lintran)
		implicit none
		logical,intent(in)::custom
		real(kind=8):: lintran(:,:)
	end subroutine getlintran
!	Find means and covariance matrices corresponding to a predictor list for the normal case.

!	beta is the parameter vector.
!	predlist is the predictor list.
!	tolsing is the tolerance for modified Cholesky decomposition.
!	covs are the covariance matrices.
!	means are the means.

	subroutine getmeancov(beta,predlist,tolsing,covs,means)
		implicit none
		real(kind=8),intent(in)::beta(:),predlist(:,:),tolsing
		real(kind=8),intent(out)::covs(:,:,:),means(:,:)
	end subroutine getmeancov
!	Find means and covariance matrices corresponding to a predictor list for the multinomial case.

!	beta is the parameter vector.
!	predlist is the predictor list.
!	quadpoint is the array of quadrature points.
!	quadweight is the array of quadrature weights.

!	covs are the covariance matrices.
!	means are the means.

	subroutine getmeancovm(beta,predlist,quadpoint,quadweight,covs,means)
		implicit none
		real(kind=8),intent(in)::beta(:),predlist(:,:),quadpoint(:,:),quadweight(:)
		real(kind=8),intent(out)::covs(:,:,:),means(:,:)
	end subroutine getmeancovm

	!
!  See what means of latent vectors are needed.
!	npredlist is the number of predicting vectors to be used.  It is relevant if allmeans is .FALSE.
!	allmeans is .TRUE. only if a mean is desired for each observed predicting vector.

	subroutine getmeans(npredlist,allmeans)
		implicit none 
		integer,intent(out)::npredlist
		logical,intent(out)::allmeans
	end subroutine getmeans

!	Get number numcodes of recodes per item.
	subroutine getnumcodes(numcodes)

		implicit none
		integer,intent(out)::numcodes(:)
	end subroutine getnumcodes
!	Get number numdistcodes of distractor recodes per item.
	subroutine getnumdistcodes(numdistcodes)
		implicit none
		integer,intent(out)::numdistcodes(:)
	end subroutine getnumdistcodes	
!	Get the number of weights for distributions.
!	slopedim gives the dimension per item.
!	The default dimension is the dimension of the transformed latent vector.
	integer function getnumweights(slopedim)
		implicit none
		integer,intent(in)::slopedim(:,:)
	end function getnumweights
!	Get number of weights for distractor distributions.
!	slopedim gives the dimension per item.
!	The default dimension is the dimension of the transformed latent vector.
    integer function getnumweightsdist(slopedim)
        implicit none
        integer,intent(in)::slopedim(:,:)
    end function getnumweightsdist

!	Some parameter settings.
!	maxit is the number of iterations for the main cycle.
!	maxita is the number of iterations for secondary cycles.
!	maxitb is the maximum number of iterations used to find the step size.
!	twostages indicates whether two rounds of iterations are used.
!	nr is for Newton-Raphson.  The alternative uses an approximation
!   to Newton-Raphson based on the Louis approximation.
!	changemin is the criterion for minimum change.
!	kappa is a maximum step size.
!	maxdalpha is the maximum permitted alpha change.
!	tau is a maximum shrinkage of a step size.
!	tol is the tolerance for the information change per item per effective observation.
!	tola is the tolerance for the posterior change for an item.
!   tolc is the criterion for appropriate maximum ratio of diagonal terms in modified Cholesky decomposition.
!	tolres is the tolerance for standard errors of residuals.
!	tolsing is the tolerance for the modified Cholesky decomposition.
	subroutine getparam(maxit,maxita,maxitb,nr,twostages,&
		changemin,kappa,maxdalpha,tau,tol,tola,tolc,tolres,tolsing)
		implicit none
		integer,intent(out)::maxit,maxita,maxitb


		logical,intent(out)::nr,twostages
		real(kind=8),intent(out)::changemin,kappa,maxdalpha,&
			tau,tol,tola,tolc,tolres,tolsing
			
	end subroutine getparam

!	read predictor list.
!	predlabels are labels for predictor vectors.
!	predlist are predictor vectors.
	subroutine getpredlist(predlabels,predlist)
		implicit none
		character(len=32),intent(out)::predlabels(:)
		real(kind=8),intent(out)::predlist(:,:)	
	end subroutine getpredlist

!	Obtain predictor names.  Put in predname.
	subroutine getprednames(predname)
		implicit none
		character(len=32),intent(out)::predname(:)
		
	end subroutine getprednames
!   Get distractor probabilities.
!   dat is response vector.
!   distdat is distractor array.
!   choices gives number of distractors per item.
!   numcatobs is number of observed categories for each item.

!   distprob gives conditional probabilities of distractors.
!   prob gives probabilities of item scores.
    function getprobdist(dat,distdat,choices,numcatobs,distprob,prob)
        implicit none
        integer,intent(in)::dat(:,:),distdat(:,:),choices(:),numcatobs(:)
        real(kind=8),intent(in)::distprob(:),prob(:)
        real(kind=8)::getprobdist(size(prob))
    end function getprobdist

!	
!	Get data on psu counts per stratum for the nstratum strata.
	subroutine getpsucount(psu,stratum,npsu)
		implicit none
		integer,intent(in)::psu(:),stratum(:)
		integer,intent(out)::npsu(:)
	end subroutine getpsucount
	
!	Set quadrature points.
!	Obtain lower and upper probabilities of a weighted sum.
!	catobsrange maps underlying to observed categories.
!	dat provides responses.
!	maxscore is the maximum score.
!	maxw provides maximum item weights.
!	minscore is the minimum score.
!	minw provides minimum item weights.
!	numcat provides the underlying number of categories per item.
!	numcatobs provides the observed number of categories per item.
!	weight is the weights.
!	beta provides parameter estimates.
!	lintran transforms the latent vector.
!	postdensity is the array of posterior densities.
!	theta is the array of quadrature points.


!	plower is the lower probabilities.
!   pupper is the upper probabilities.

!
    subroutine getpwt(catobsrange,dat,maxscore,maxw,minscore,minw,&
        numcat,numcatobs,weight,&
        beta,&
        lintran,&
        postdensity,theta,&
        plower,pupper)

        implicit none
        integer,intent(in)::catobsrange(:,:),dat(:,:),maxscore,maxw(:),minscore,minw(:),&
            numcat(:),numcatobs(:),weight(:)
        real(kind=8),intent(in)::beta(:),lintran(:,:),&
            postdensity(:,:),theta(:,:,:)
        real(kind=8),intent(out)::plower(:),pupper(:)
    end subroutine getpwt
!	Obtain lower and upper probabilities of a weighted distractor sum.
!	catobsrange maps underlying to observed categories.
!   choices indicates the number of distractors per item.
!	dat provides responses.
!   distdata provides distractor data.
!   distmap is the distractor map.
!	maxscore is the maximum score.
!	maxw provides maximum item weights.
!	minscore is the minimum score.
!	minw provides minimum item weights.
!	numcat provides the underlying number of categories per item.
!	numcatobs provides the observed number of categories per item.
!	weightdist is the weights.
!	beta provides parameter estimates.
!	lintran transforms the latent vector.
!	postdensity is the array of posterior densities.
!	theta is the array of quadrature points.


!	plower is the lower probabilities.
!   pupper is the upper probabilities.

!
    subroutine getpwtdist(catobsrange,choices,dat,distdat,distmap,&
        maxscore,maxw,minscore,minw,&
        numcat,numcatobs,weightdist,&
        beta,distprob,&
        lintran,&
        postdensity,theta,&
        plower,pupper)

        implicit none
        integer,intent(in)::catobsrange(:,:),choices(:),dat(:,:),distdat(:,:),distmap(:)
        integer,intent(in)::maxscore,maxw(:),minscore,minw(:)
        integer,intent(in)::numcat(:),numcatobs(:),weightdist(:)
        real(kind=8),intent(in)::beta(:),distprob(:),lintran(:,:),&
            postdensity(:,:),theta(:,:,:)
        real(kind=8),intent(out)::plower(:),pupper(:)
    end subroutine getpwtdist




!	coord specifies the list of products to use.
!	dimpoints is the number of quadrature points per dimension.
!	cross indicates a cross-polytope.
!	even is for even spacing.
!	gausshermite is the indicator for Gauss-Hermite quadrature.
!	grid indicates a grid of points.
!	simplex indicates a simplex.
!	cweight modifies quadrature weights.
!	pspread is for point spread for even.
!	quadpoint provides quadrature points.
!	quadweight provides quadrature weights.
	subroutine getquad(coord,dimpoints,cross,equalpoints,&
		even,gausshermite,grid,simplex,cweight,pspread,&
		quadpoint,quadweight)
		implicit none
		integer,intent(in)::coord(:,:),dimpoints(:)
		logical,intent(in)::cross,equalpoints,even,gausshermite,grid,simplex
		real(kind=8),intent(in)::cweight(:),pspread
		real(kind=8),intent(out)::quadpoint(:,:),quadweight(:)
	end subroutine getquad

!	Find the size specifications for the quadrature points.
!	dimpoints gives the number  of points per diemnsion if a
!	grid is used.
!	If a grid is used, then the grid has ngrid points.
!	In all cases, the number of quadrature points is nquad.
!	cross-polytope quadrature is used if cross is .true.
!	Each grid dimension has the same number of points if
!	equalpoints is .true.
!	If even is .true. and a grid is used, then
!		evenly spaced quadrature points are used.
!	A complete grid is used if a grid is used and fullgrid is .true.
!	If guasshermite is .true., even is .false., and a grid is used,
!		then Gauss-Hermite quadrature points and weights are
!		produced.
!	A grid is used if grid is .true. and cross and simplex are .true.
!	normal indicates that the latent vector has a normal distribution.
!	Simplex quadrature is used if simplex is .true. and cross is .false.
!	The subroutine returns grid and fullgrid as .false. if  simplex  or
!		cross is .true.
!	The subroutine returns simplex as .false. if cross is .true.
!	In the default case, cross is .false., even is .false., fullgrid is .true.,
!		gausshermite is .true., grid is .true., and simplex is .false.
!	The default value of each element of dimpoints is 5 if Gauss-Hermite
!		quadrature is used and 6 if a grid is used and Gauss-Hermite
!		quadrature is not used.  If Gauss-Hermite quadrature is used, then no
!		element of dimpoints is permitted to exceed 30.
!		Any element of dimpoints is forced to be at least 2 if a grid is used.
!	If a grid is used but nquad is less than 2, then fullgrid is set to .true.
!	If nquad exceeds ngrid and a grid is used, then nquad is set to ngrid.
!	If a full grid is used, then nquad is set to ngrid.
!	If cross is .true., then nquad is set to twice size(dimpoints).
!	If simplex is returned as .true., then nquad is size(dimpoints)+1.
!	If simplex, cross, and grid are .false., then nquad must be given and
!		at least 2.  Otherwise, processing ends with an error message.


	subroutine getquadsize(dimpoints,ngrid,nquad,&
		cross,equalpoints,even,fullgrid,gausshermite,grid,normal,simplex,&
		pspread)

		implicit none


		integer,intent(out)::dimpoints(:),nquad,ngrid
		logical,intent(out)::cross,equalpoints,even,fullgrid,gausshermite,grid,normal,simplex
		real(kind=8),intent(out)::pspread
	end subroutine getquadsize
	
!	Obtain recodes.
!	numcodes gives number of recodes per item.
!	recodetab gives recodes.
	subroutine getrecodes(numcodes,recodetab)
		implicit none
		integer,intent(in)::numcodes(:)
		integer,intent(out)::recodetab(:,:)
	end subroutine getrecodes	
!
!	Obtain counts numberscales of scale scores per weighted sum.
	subroutine getscalecount(numberscales)
		implicit none
		integer,intent(out)::numberscales
	end subroutine getscalecount	
!	Obtain scale scores for weighted sum.
!	maxscore and minscore give maximum and minimum weighted sums.
!	numberscales is the number of scales.
!	Names are in scalename and scales are in scale.
	subroutine getscales(maxscore,minscore,numberscales,scalename,scale)
		implicit none
		character(len=32),intent(out)::scalename(:)
		integer,intent(in)::maxscore,minscore,numberscales
		real(kind=8),intent(out)::scale(minscore:maxscore,numberscales)
	end subroutine getscales
!
!	Get item scores.  Use numcatobs to find dimension and defaults.
!	numcatobs is the number of observed categories per item
	function getscores(numcatobs)
		implicit none
		integer,intent(in)::numcatobs(:)
		real(kind=8)::getscores(sum(numcatobs))
	end function getscores

!	Obtain skill specifications
!	factorname is the array of factor names.
!	skillname is the array of skill names.
!	rasch is the array of indicators for constant slopes for skills.
!	raschslope1 is the array of indicators for constant slopes of 1 for skills.		
	subroutine getskillspecs(factorname,skillname,rasch,raschslope1)
		implicit none
		character(len=32),intent(in)::factorname(:)
		character(len=32),intent(out)::skillname(:)
		logical,intent(out)::rasch(:),raschslope1(:)
	end subroutine getskillspecs

!	Find points for item response function computations.  Place in thetacheck.

	subroutine getthetas(thetacheck)
		implicit none
		real(kind=8),intent(out)::thetacheck(:,:)
	end subroutine getthetas
!	Check on number of points for item response functions and whether quadrature points are used.
!	quadpoint provides regular quadrature points.
!	irfcount is the desired number of points for item response functions.
!   The default is the number of quadrature points for
!	quadrature.
!	custompoints is only true if special points are to be read.
	subroutine getthetasize(quadpoint,irfcount,custompoints)
		implicit none
		integer,intent(out)::irfcount
		logical,intent(out)::custompoints
		real(kind=8),intent(in)::quadpoint(:,:)
	end subroutine getthetasize
	
	
!	Get  title
	subroutine gettitle(title)
		implicit none
		character(len=80),intent(out)::title

	end subroutine gettitle
!	Obtain transformation.
!	dimlatin is input dimension.
!	transname is output labels.
	subroutine gettrans(dimlatin,transname)
		implicit none
		character(len=32),intent(out)::transname(:)
		integer,intent(in)::dimlatin
	end subroutine gettrans	

!	Get integer weights for weighted sum.
!	skillname is the array of skill names.
!	numcatobs is the array of numbers of categories in items.
!	slopedim is the array of dimensions for each factor and item.
!	weightnames are names of weights and weights are weights.
	subroutine getweights(skillname,numcatobs,slopedim,weightnames,weights)
		implicit none
		character(len=32),intent(in)::skillname(:)
		integer,intent(in)::numcatobs(:),slopedim(:,:)
		character(len=32),intent(out)::weightnames(:)
		integer,intent(out)::weights(:,:)
	end subroutine getweights
!	Get integer weights for weighted sum.
!	skillname is the array of factor names.
!	choices is the array of numbers of distractors of items.
!	distmap is the distractor map.
!   numcatobs counts categories per observed item.
!	slopedim is the array of dimensions for each factor and item.
!	weightnamesdist are names of weights and weights are weights.
    subroutine getweightsdist(skillname,choices,distmap,numcatobs,&
        slopedim,weightnamesdist,weightsdist)
        implicit none
        character(len=32),intent(in)::skillname(:)
        integer,intent(in)::choices(:),distmap(:),numcatobs(:),slopedim(:,:)
        character(len=32),intent(out)::weightnamesdist(:)
        integer,intent(out)::weightsdist(:,:)
    end subroutine getweightsdist
!	Get weighted sum wtsum with row names wtname.
!	skillname is array of skill names.
!	catobsrange is array of ranges of underlying categories.
!	numcatobs is array of observed categories per item.
!	slopedim gives factors assigned to items.
	subroutine getwtsum(skillname,catobsrange,numcatobs,slopedim,wtname,wtsum)
		implicit none
		character(len=32),intent(in)::skillname(:)
		integer,intent(in)::catobsrange(:,:),numcatobs(:),slopedim(:,:)
		character(len=32),intent(out)::wtname(:)
		real(kind=8),intent(out)::wtsum(:,:)
	end subroutine getwtsum
!	Get weighted sum for response choices.
!	skillname is array of skill names.
!	choices are counts for choices for each item.
!	distmap is the distractor map.
!   numcatobs counts categories per observed item.
!	slopedim gives factors assigned to items.
!	wtnamedist gives name of sum.
!	wtsumdist specifies the sum.
    subroutine getwtsumdist(skillname,choices,distmap,numcatobs,slopedim,wtnamedist,wtsumdist)
        implicit none
        character(len=32),intent(in)::skillname(:)
        integer,intent(in)::choices(:),distmap(:),numcatobs(:),slopedim(:,:)
        character(len=32),intent(out)::wtnamedist(:)
        real(kind=8),intent(out)::wtsumdist(:,:)
    end subroutine getwtsumdist
!   Perform guessing tests.
!	dat provides responses.
!	npsu is the number of psus's per stratum.
!	nstratum is the number of strata.
!	numcat provides the underlying number of categories per item.
!	psu indicates the psu of an observation.
!	stratum indicates the stratum of an observation.
!	complx indicates complex sampling.
!	stratify is .true. for stratified sampling.
!	usepsu is .true. for primary sampling units within strata.

!	eacovgaminv_louis is the inverse of the Louis information plus the
!		constraint component of the negative hessian.
!	gradsdes provides gradients for observations.
!	lintran transforms the latent vector.
!	obsweight provides weights.
!	postdensity is the array of posterior densities.
!	theta is the array of quadrature points.
!   guessres gives the residual total.
!   guessresa is the adjusted residual.
!   stdguessres is the standard error of guessres.

    subroutine guesstest(dat,npsu,nstratum,numcat,psu,stratum,&
        complx,stratify,usepsu,beta,eacovgaminv_louis,&
        gradsdes,lintran,obsweight,postdensity,theta,guessres,guessresa,stdguessres)
        implicit none
        integer,intent(in)::dat(:,:),npsu(:),nstratum,numcat(:),&
            psu(:),stratum(:)
        logical,intent(in)::complx,stratify,usepsu
        real(kind=8),intent(in)::beta(:),eacovgaminv_louis(:,:),&
            gradsdes(:,:),lintran(:,:),obsweight(:),postdensity(:,:),theta(:,:,:)
        real(kind=8),intent(out)::guessres(:),guessresa(:),stdguessres(:)
    end subroutine guesstest
		
!
!
!	Find the Gilula-Haberman correction.  
!	eacov is the estimated asymptotic covariance matrix under the model.
!	inform is the estimated information obtained by the Louis approach.
!	loglik is the log likelihood.
!	totalitems is the total weighted number of presented items.
!
!
	real(kind=8) function gh(eacov,inform,loglik,totalitems)
		implicit none
		real(kind=8),intent(in)::eacov(:,:),inform(:,:),loglik,totalitems
	end function gh
	
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

!	Set intercepts by a regression.    


!	catobsrange maps observed to underlying categories.
!	dat is the array of responses.
!	numcat is the number of underlying categories per item.
!	numcatobs is the number of observed categories per item.
!	design is the design matrix.
!	obsweight is the vector of observation weights.
!	offset is the offset for beta.
!	tolsing is a tolerance for modified Cholesky decomposition.
!	gamma is the vector of starting value.
	subroutine interceptstart(catobsrange,dat,numcat,numcatobs,&
        design,obsweight,offset,tolsing,gamma)
		implicit none
		integer,intent(in)::catobsrange(:,:),dat(:,:),numcat(:),numcatobs(:)
		real(kind=8),intent(in)::design(:,:),obsweight(:),offset(:),tolsing
		real(kind=8),intent(out)::gamma(:)
	end subroutine interceptstart
	
!	Compute inverse for n by n matrix with modified Cholesky
!	decomposition tri.
	function invert(tri)
		implicit none
		real(kind=8),intent(in)::tri(:,:)
		real(kind=8)::invert(size(tri,1),size(tri,1))
	end function invert	
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
!   prob is the array of marginal probabilities.
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
		integer,intent(in)::catobsrange(:,:),dat(:,:),&
            npsu(:),nstratum,numcat(:),numcatobs(:),&
			psu(:),stratum(:)
		logical,intent(in)::complx,resid,stratify,usepsu
		real(kind=8),intent(in)::beta(:),eacovgaminv_louis(:,:),&
			gradsdes(:,:),indvar(:,:),lintran(:,:),obsweight(:),&
			postdensity(:,:),prob(:),theta(:,:,:),thetacheck(:,:),tolres
		real(kind=8),intent(out)::fitirf(:,:),obsirf(:,:),presentedirf(:),&
			residairf(:,:),residirf(:,:),&
			stdresidirf(:,:)
	end subroutine irf
	
!	Total items presented.
!	dat contains response data.  obsweight are observation weights.
!	numcatobs contains observed category counts by item.
	real(kind=8) function itemspresented(dat,numcatobs,obsweight)
		implicit none
		integer,intent(in)::dat(:,:),numcatobs(:)
		real(kind=8),intent(in)::obsweight(:)	
	end function itemspresented
	
!	Map beta.  First row is item, second is category, third is first dimension,
!   fourth is second dimension, and fifth is
!	predictor number.
!	dimlatin is dimension of latent vector.
!	dimlatout is dimension of transformed latent vector.
!	npred is number of predictors.
!	numcat is number of underlying categories per item.
	function mapbeta(dimlatin,dimlatout, npred,numcat)
		implicit none
		integer,intent(in)::dimlatin,dimlatout,npred,numcat(:)
		integer::mapbeta(5,sum(numcat)*(dimlatout+1)+npred*dimlatin*(dimlatin+3)/2)
	end function mapbeta
	
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
!	lintran transforms the latent vector.
!	obsweight provides weights.
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
		integer,intent(in)::catobsrange(:,:),dat(:,:),&
			npsu(:),nstratum,numcat(:),numcatobs(:),&
			psu(:),stratum(:)
		logical,intent(in)::complx,resid,stratify,usepsu
		real(kind=8),intent(in)::beta(:),eacovgaminv_louis(:,:),&
			gradsdes(:,:),lintran(:,:),obsweight(:),&
			postdensity(:,:),theta(:,:,:),tolres
		real(kind=8),intent(out)::fitmarg(:),fitpmarg(:),obsmarg(:),obspmarg(:),presented(:),&
			relmarg(:),residamarg(:),residmarg(:),residpmarg(:),&
			seobsmarg(:),stdobsmarg(:),stdobspmarg(:),&
			stdresidmarg(:),stdresidpmarg(:),stobsmarg(:)
	end subroutine marginaldist



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
		integer,intent(in)::catobsrange(:,:),dat(:,:),npsu(:),nstratum,numcat(:),numcatobs(:),&
			psu(:),stratum(:)
		logical,intent(in)::complx,resid,stratify,usepsu
		real(kind=8),intent(in)::beta(:),eacovgaminv_louis(:,:),&
			gradsdes(:,:),lintran(:,:),obsweight(:),&
			postdensity(:,:),scores(:),theta(:,:,:),tolres
		real(kind=8),intent(out)::fitmargs2(:,:),fitpmargs2(:,:),&
            obsmargs2(:,:),obspmargs2(:,:),presenteds2(:,:),&
			residamargs2(:,:),residmargs2(:,:),residpmargs2(:,:),&
			stdobsmargs2(:,:),stdobspmargs2(:,:),&
			stdresidmargs2(:,:),stdresidpmargs2(:,:)
	end subroutine margindists2
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
!	lintran transforms the latent vector.!	obsweight provides weights.
!	postdensity is the array of posterior densities.
!	theta is the array of quadrature points.
!	tolres is the rsidual tolerance.
!	fitmarg2 is the fitted marginal total.
!	fitpmarg2 is the  fitted marginal proportion.
!	obsmarg2 is the observed marginal total.
!	obspmarg2 is the observed marginal proportion.
!	presented2 counts weighted items presented in pairs.
!	relmarg2 is the reliability of the item category indicator pairs.
!	residmarg2 is the residual for the marginal total.
!	residamarg2 is the adjusted residual.
!	residpmarg2 is the residual for the marginal proportion.
!	seobsmarg2 is the standard error for the observed indicator.
!	stdobsmarg2 is the standard deviation of the observed marginal.
!	stdobspmarg2 is the asymptotic standard deviation of the observed 
!		marginal fraction.
!	stdresidmarg2 is the asymptotic standard deviation of the residual marginal.
!	stdresidpmarg2 is the asymptotic standard deviation of the residual 
!		marginal fraction.
!	stobsmarg2 is the standard deviation of the indicator.
	subroutine margindist2(catobsrange,dat,npsu,nstratum,numcat,numcatobs,&
		psu,stratum,complx,resid,stratify,usepsu,&
		beta,eacovgaminv_louis,gradsdes,lintran,&
		obsweight,postdensity,theta,tolres,&
		fitmarg2,fitpmarg2,obsmarg2,obspmarg2,presented2,&
		residamarg2,residmarg2,residpmarg2,&
		stdobsmarg2,stdobspmarg2,&
		stdresidmarg2,stdresidpmarg2)
		integer,intent(in)::catobsrange(:,:),dat(:,:),&
			npsu(:),nstratum,numcat(:),numcatobs(:),&
			psu(:),stratum(:)
		logical,intent(in)::complx,resid,stratify,usepsu
		real(kind=8),intent(in)::beta(:),eacovgaminv_louis(:,:),&
			gradsdes(:,:),lintran(:,:),obsweight(:),&
			postdensity(:,:),theta(:,:,:),tolres
		real(kind=8),intent(out)::fitmarg2(:,:),&
            fitpmarg2(:,:),obsmarg2(:,:),obspmarg2(:,:),presented2(:,:),&
			residamarg2(:,:),residmarg2(:,:),residpmarg2(:,:),&
			stdobsmarg2(:,:),stdobspmarg2(:,:),&
			stdresidmarg2(:,:),stdresidpmarg2(:,:)
	end subroutine margindist2
    
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
	end subroutine margindistwtsum
	
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
!	presentedpred counts weighted items presented.
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
		fitppred,fitpred,obsppred,obspred,presentedpred,&
		residapred,residppred,residpred,stdobsppred,stdobspred,stdresidppred,stdresidpred)
		implicit none
		integer,intent(in)::catobsrange(:,:),dat(:,:),&
			npsu(:),nstratum,numcat(:),numcatobs(:),&
			psu(:),stratum(:)
		logical,intent(in)::complx,resid,stratify,usepsu
		real(kind=8),intent(in)::beta(:),eacovgaminv_louis(:,:),&
			gradsdes(:,:),extvar(:,:),indvar(:,:),lintran(:,:),obsweight(:),&
			postdensity(:,:),theta(:,:,:),tolres
		real(kind=8),intent(out)::fitppred(:,:),fitpred(:,:),&
			obsppred(:,:),obspred(:,:),presentedpred(:),residapred(:,:),residppred(:,:),residpred(:,:),&
			stdobsppred(:,:),stdobspred(:,:),stdresidppred(:,:),stdresidpred(:,:)
	end subroutine marginpred

!	Obtain interactions of category indicators and transformed latent vectors.
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
!		constraint component of the negative hessian..
!	gradsdes provides gradients for observations.
!	lintran transforms the latent vector.
!	obsweight provides weights.
!	postdensity is the array of posterior densities.
!	theta is the array of quadrature points.
!	tolres is the residual tolerance.
!	fitptheta is the average of fitted product of item indicator and latent vectors.
!	fittheta is the total of fitted products of item indicator and latent vectors.
!	obsptheta is the observed average product of item indicator and latent vectors.
!	obstheta is the observed total product of item indicator and predictors.
!	presentedtheta counts weighted items presented.
!	residatheta is the adjusted residual of average of products.
!	residptheta is the residual of average of products.
!	residtheta is the residual of total of products.
!	stdobsptheta is the asymptotic standard deviation of the observed average of products.
!	stdobstheta is the asymptotic standard deviation of the observed total of products.
!	stdresidptheta is the asymptotic standard deviation of residual average of products.
!	stdresidtheta is the asymptotic standard deviation of the residual total of products.



    subroutine margintheta(catobsrange,dat,npsu,nstratum,numcat,numcatobs,&
        psu,stratum,complx,resid,stratify,usepsu,&
        beta,eacovgaminv_louis,gradsdes,lintran,&
        obsweight,postdensity,theta,tolres,&
        fitptheta,fittheta,obsptheta,obstheta,presentedtheta,&
        residatheta,residptheta,residtheta,stdobsptheta,&
        stdobstheta,stdresidptheta,stdresidtheta)


        implicit none
        integer,intent(in)::catobsrange(:,:),dat(:,:),&
        npsu(:),nstratum,numcat(:),numcatobs(:),&
        psu(:),stratum(:)
        logical,intent(in)::complx,resid,stratify,usepsu
        real(kind=8),intent(in)::beta(:),eacovgaminv_louis(:,:),&
            gradsdes(:,:),lintran(:,:),obsweight(:),&
            postdensity(:,:),theta(:,:,:),tolres
        real(kind=8),intent(out)::fitptheta(:,:),fittheta(:,:),&
            obsptheta(:,:),obstheta(:,:),presentedtheta(:),residatheta(:,:),&
            residptheta(:,:),residtheta(:,:),&
            stdobsptheta(:,:),stdobstheta(:,:),stdresidptheta(:,:),stdresidtheta(:,:)
    end subroutine margintheta

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
!	lintran transforms the latent vector.!	obsweight provides weights.
!	postdensity is the array of posterior densities.
!	theta is the array of quadrature points.
!	tolres is the residual tolerance.
!	fitwtitem is the total of fitted products of item indicator and weighted sum.
!	fitpwtitem is the average of fitted product of item indicator and weighted sum.
!	obswtitem is the observed total products of item indicator and weighted sum.
!	obspwtitem is the observed average product of item indicator and weighted sum.
!	presentedwtitem counts weighted items presented with defined weighted sums.
!	residwtitem is the residual total of products.
!	residawtitem is the adjusted residual of total of products.
!	residpwtitem is the residual average of products.
!	stdobswtitem is the asymptotic standard deviation of the observed total of products.
!	stdobspwtitem is the asymptotic standard deviation of the observed average of products.
!	stdresidwtitem is the asymptotic standard deviation of the residual total of products.
!	stdresidpwtitem is the asymptotic standard deviation of residual average of products.


	subroutine marginwtitem(catobsrange,dat,maxw,minw,npsu,nstratum,numcat,numcatobs,&
		psu,stratum,weight,complx,resid,stratify,usepsu,&
		beta,eacovgaminv_louis,gradsdes,lintran,&
		obsweight,postdensity,theta,tolres,&
		fitpwtitem,fitwtitem,obspwtitem,obswtitem,presentedwtitem,&
		residawtitem,residpwtitem,residwtitem,stdobspwtitem,stdobswtitem,stdresidpwtitem,stdresidwtitem)
		implicit none
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
	end subroutine marginwtitem
	

!	print square matrix.
!	comment is the type of matrix.
!	rowname provides row names.
!	unitmat is the unit number to use.
!	mat provides the matrix.
	subroutine matoutput(comment,rowname,unitmat,mat)
		implicit none
		character(len=*),intent(in)::comment
		character(len=*),intent(in)::rowname(:)
		integer,intent(in)::unitmat
		real(kind=8),intent(in)::mat(:,:)
	end subroutine matoutput
!   Print maximum likelihood estimates of theta or a function of theta.

!    factorname gives names of elements of the latent vector.
!    id is individual id.
!    transname is used if dimtrans is positive.

!    catobsrange is used to relate observed and underlying categories.
!    dat are the responses.
!    dimtrans is the transformation dimension.
!    unitmp is the unit for output.
!    maxita is used to determine the number of iterations to use to find the maximum posterior density.
!    maxitb is used for subiterations.
!    npred is the number of predictors.
!    numcat provides ranges for underlying observations,
!    numcatobs provides ranges of observations, mask indicates items to use for computations.
!    eapid indicates if ids are put out.
!    obsmask is an observation mask that can be used to restrict attention to a portion of the response.
!    beta is the parameter vector.
!    changemin is used for minimum change.
!    lintran is the linear transformation of the latent vector.
!    maxdtheta is the maximum permitted change in the approximation of the maximum in one step.

!    tau is used for minimum step sizes.
!    tol is used to determine the desired accuracy of the maximization.
!    tolsing is used to determine the tolerance for the Cholesky decomposition.

    subroutine maxl(factorname,id,transname,catobsrange,dat,dimtrans,maxita,maxitb,&
        npred,numcat,numcatobs,unitml,eapid,obsmask,beta,changemin,&
        lintran,maxdtheta,tau,tol,tolsing)
        implicit none
        character(len=32),intent(in)::factorname(:),id(:),transname(:)
        integer,intent(in)::catobsrange(:,:),dat(:,:),dimtrans,maxita,maxitb,npred,numcat(:),numcatobs(:),unitml
        logical,intent(in)::eapid,obsmask(:)
        real(kind=8),intent(in)::beta(:),changemin,lintran(:,:),maxdtheta,tau,tol,tolsing
    end subroutine maxl


	
!	print maximum posterior likelihood data for latent vector.
!	factorname gives names of elements of the latent vector.
!	id is individual id.
!	unitmp is the unit for output.
!	eapid indicates if ids are put out.
!	alpha are locations of maxima.
!	cholnhess is the array of modified Cholesky decompositions.
	subroutine mpoutput(factorname,id,unitmp,eapid,alpha,cholnhess)
		implicit none
		character(len=32),intent(in)::factorname(:)
		character(len=32),intent(in),optional::id(:)
		integer,intent(in)::unitmp
		logical,intent(in)::eapid
		real(kind=8),intent(in)::alpha(:,:),cholnhess(:,:,:)
	end subroutine mpoutput

!	Maximum likelihood for log-linear MIRT models with possible indirect observation.
!	Summary of variables.
!	comment is used to comment on the iteration progress report.
!	catobsrange is the range of underlying categories per for each observed category.
!	constdim1 is the number of linear constraints not treated via least squares.
!	dat is the data matrix.
!	maxit is the maximum number of iterations.
!	maxita is the maximum number of iterations used to find alpha.
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
!	tola is the convergence criterion for alp change.
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
!	nhessdeschol is the modified Cholesky decomposition of nhessdes.
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
		character(len=*),intent(in)::comment
		integer,intent(in)::catobsrange(:,:),constdim1,dat(:,:),maxit,maxita,maxitb,numcat(:),numcatobs(:),rbetamap(:,:),unitno
		logical,intent(in)::normal,nr,printprog,printprogstd,proport
		real(kind=8),intent(in)::changemin,constmat(:,:),constvec(:),design(:,:),indvar(:,:),kappa,&
			lintran(:,:),maxdalpha,obsweight(:),offset(:),&
			quadpoint(:,:),quadweight(:),&
			rdesign(:,:),tau,&
			tol,tola,tolc,tolsing
		
		real(kind=8),intent(inout)::alpha(:,:),beta(:),cholnhess(:,:,:),gamma(:)
		real(kind=8),intent(out)::constquad(:,:),grad(:),graddes(:),grads(:,:),&
			loglik,nhess(:,:),nhessdes(:,:),nhessdeschol(:,:),&
			nhessdesplus(:,:),postdensity(:,:),prob(:),theta(:,:,:)
	end subroutine ndmmaxlik
	
!	Total observations presented.
!	dat contains response data.  obsweight are observation weights.
!	numcatobs contains observed category counts by item.
	real(kind=8) function obspresented(dat,numcatobs,obsweight)
		implicit none
		integer,intent(in)::dat(:,:),numcatobs(:)
		real(kind=8),intent(in)::obsweight(:)
	end function obspresented

!	Find observed scale score for weighted sums.
!	Responses are in dat.
!	maxscore is the maximum raw score.
!	maxw provides maximum item weights.
!	minscore is the minimum raw score.
!	minw provides minimum item weights.
!	numberscales is the number of scale scores.
!	numcatobs provides ranges for observations.
!	weight is the weight.
!	scale is the matrix of scale transformations.
    subroutine obsscale(dat,maxscore,maxw,minscore,&
        minw,numberscales,numcatobs,weight,presence,&
        scale,observedscales)
        implicit none
        integer,intent(in)::dat(:,:),maxscore,maxw(:),minscore,minw(:),numberscales,&
            numcatobs(:),weight(:)
        logical,intent(out)::presence(:)

        real(kind=8),intent(in)::scale(numberscales,minscore:maxscore)
        real(kind=8),intent(out)::observedscales(:,:)
    end subroutine obsscale

!	print observed scale score.
!	comment is used for identification.
!	factorname gives names of elements of the vector.
!	id is individual id.
!	unitobs unitobs is the unit for output.
!	obsid indicates if ids are put out.
!   presence indicates what is missing.
!	observed gives observations.
    subroutine obsoutput(comment,factorname,id,unitobs,obsid,presence,observed)
        implicit none
        character(len=*),intent(in)::comment
        character(len=32),intent(in)::factorname(:)
        character(len=32),intent(in),optional::id(:)
        integer,intent(in)::unitobs
        logical,intent(in)::obsid,presence(:)
        real(kind=8),intent(in)::observed(:,:)
    end subroutine obsoutput

!   Print observed scale summaries.
!   scalenames gives names of scale scores.
!   unitobscale is unit for output.
!   printobscaleres is .true. if residual information is to be given.
!   aresid gives adjusted residuas.
!   covobsave gives the covariance matrix of the observed averages.
!   covobssum gives the covariance matrix of the observed sums.
!   covresave gives the covariance matrix of the average residuals.
!   covressum gives the covariance matrix of the residual sums.
!   fitave gives the average fitted values.
!   fitsum gives the sum of fitted values.
!   obsave gives the average observed scores.
!   obssum gives the sums of observed scores.
!   resave gives the average residual.
!   ressum gives the sum of residuals.
!   sdobsave gives the standard deviation of the average observed score.
!   sdobssum gives the standard deviation of the sum of observed scores.
!   sdresave gives the standard deviation of the average residual.
!   sdressum gives the standard deviation of the sum of residuals.
!   totsum is the total weighted sample sizee for observations.

    subroutine obsscaleoutput(scalename,unitobsscale,printobsscaleres,&
        aresid,covobsave,covobssum,covresave,covressum,&
        fitave,fitsum,obsave,obssum,resave,&
        ressum,sdobsave,sdobssum,sdresave,sdressum,totsum)
        implicit none
        character(len=32),intent(in)::scalename(:)
        integer,intent(in)::unitobsscale
        logical,intent(in)::printobsscaleres
        real(kind=8),intent(in)::aresid(:),covobsave(:,:),covobssum(:,:),&
            covresave(:,:),covressum(:,:),&
            fitave(:),fitsum(:),&
            obsave(:),obssum(:),&
            resave(:),ressum(:),&
            sdobsave(:),sdobssum(:),sdresave(:),sdressum(:),totsum
        end subroutine obsscaleoutput

!	Print reliability results for observed scores.
!	comment is used to identify output.
!	factorname is variable names.
!	unitobsrel is unit.
!	covobs is the covariance matrix of the observations vector.
!	covres is the covariance matrix of the residuals.
!   obsave is the average of the scores.
!	rel contains reliabilities.

    subroutine obsreloutput(comment,factorname,unitobsrel,&
        covobs,covres,obsave,rel)
        implicit none
        character(len=*),intent(in)::comment
        character(len=32),intent(in)::factorname(:)
        integer,intent(in)::unitobsrel
        real(kind=8),intent(in)::covobs(:,:),covres(:,:),obsave(:),rel(:)
    end subroutine obsreloutput

!	Some basic summary information.
!	filename is the name of the data file.
!	nitems is the number of items.
!	nobs is the number of observations.
!	npred is the number of predictors.
!	unitno is the unit number used.
!	totalitems is weighted number of items presented.
!	totalobs is weighted number of observations presented.


	subroutine outputdata(filename,nitems,nobs,npred,unitno,totalitems,totalobs)
		implicit none
		character(len=32)::filename
		integer,intent(in)::nitems,nobs,npred,unitno
		real(kind=8),intent(in)::totalitems,totalobs
	end subroutine outputdata
	
!	Get output units.  See unitdef for definitions.

	subroutine outputfiles()
	end subroutine outputfiles
!	Print gradients.
!	id is the observation id.
!	paramname is the array of parameter names.
!	unitgrad is the output unit.
!	useid indicates if individual identifications are present in id.
!	gradsdes is the array of gradients.
		subroutine outputgrad(id,paramname,unitgrad,useid,gradsdes)
			implicit none
			character(len=32),intent(in)::id(:)
			character(len=64),intent(in)::paramname(:)
			integer,intent(in)::unitgrad
			logical,intent(in)::useid
			real(kind=8),intent(in)::gradsdes(:,:)
		end subroutine outputgrad

!	print probabilities.
!	id is individual id. prob for model.
!	unitprob is output unit.
!	eapid indicates if ids are put out.
!	prob is array of probabilities.
	subroutine outputprob(id,unitprob,eapid,prob)
		implicit none
		character(len=32),intent(in),optional::id(:)
		integer,intent(in)::unitprob
		logical,intent(in)::eapid
		real(kind=8),intent(in)::prob(:)
	end subroutine outputprob
!	print lower and upper probabilities.
!	id is individual id. prob for model.
!   weightname is name of weight.
!	unitpwt is output unit.
!	eapid indicates if ids are put out.
!	plower is array of lower probabilities.
!   pupper is array of upper probabilities.
    subroutine outputpwt(id,weightname,unitpwt,eapid,plower,pupper)
        implicit none

        character(len=32),intent(in)::id(:)
        character(len=32),intent(in)::weightname
        integer,intent(in)::unitpwt
        logical,intent(in)::eapid
		real(kind=8),intent(in)::plower(:),pupper(:)
    end subroutine outputpwt
	
!	Output the title on unitno
	subroutine outputtitle(title,unitno)
		implicit none
		character(len=80),intent(in)::title
		integer,intent(in)::unitno
	end subroutine outputtitle


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
		integer,intent(in)::catobsrange(:,:),dat(:,:),maxita,numcat(:),numcatobs(:)
		logical,intent(in)::mask(:),normal
		real(kind=8),intent(in)::beta(:),changemin,indvar(:,:),lintran(:,:),&
			maxdalpha,quadpoint(:,:),quadweight(:),&
			tau,tola,tolsing
		real(kind=8),intent(inout)::alpha(:,:),cholnhess(:,:,:)
		real(kind=8),intent(out)::postdensity(:,:),theta(:,:,:)
	end subroutine posterior

!	print posterior distribution data for latent vector.
!	factorname gives names of elements of the latent vector.
!	id is ids.
!	Output unit is unitpost.
!	eapid indicates if ids are put out.
!	postdensity are weights.
!	theta are adaptive quadrature points.  
	subroutine posterioroutput(factorname,id,unitpost,eapid,postdensity,theta)
		implicit none
		character(len=32),intent(in),optional::id(:)
		character(len=32),intent(in)::factorname(:)
		integer,intent(in)::unitpost
		logical,intent(in)::eapid
		real(kind=8),intent(in)::postdensity(:,:),theta(:,:,:)
	end subroutine posterioroutput
!   Output conditional probabilities for distractors.
!   choices provides the number of choices for each item.
!   numcatobs provides the number of categories for item scores.
!   distfreq provides distractor frequencies.
!   distprob provides conditional distractor probabilities.
!   scorefreq provides item score frequencies.
!	sdfreq gives standard errors for frequencies.
!	sdprob gives standard errors for probabilities.
    subroutine printfreqs(itemname,choices,numcatobs,unitfreq,&
        distmap,distfreq,distprob,scorefreq,sdfreq,sdprob)
        implicit none
		character(len=32),intent(in)::itemname(:)

        integer,intent(in)::choices(:),distmap(:),numcatobs(:),unitfreq
        
       
        real(kind=8),intent(in)::distfreq(:),distprob(:),scorefreq(:),sdfreq(:),sdprob(:)
    end subroutine printfreqs
!
!	Print guessing tests.
!	itemname contains item names.
!	numcat counts underlying categories.
!	unitguess provides the unit for output.
!	guessres is the residual totals.
!	guessresa is the adjusted residuals.
!	stdguessres is the standard errors of residual totals.


    subroutine printguesstest(itemname,numcat,unitguess,&
        guessres,guessresa,stdguessres)
        implicit none
        character(len=32),intent(in)::itemname(:)
        integer,intent(in)::numcat(:),unitguess

        real(kind=8),intent(in)::guessres(:),guessresa(:),stdguessres(:)
    end subroutine printguesstest
!
!	Print estimated item response functions for selected points.
!	factorname
!	itemname contains item names.
!	numcatobs counts observed categories.
!	unitirf provides the unit for output.
!	resid indicates if residuals are displayed.
!	fitirf is the unconditional estimated item response functions.
!	obsirf is the conditional estimated item response functions.
!	presentedirf counts weighted items presented.
!	residairf is the adjusted residuals.
!	residirf is the residuals.
!	thetacheck is the array of selected points

	subroutine printirfs(factorname,itemname,numcatobs,unitirf,resid,&
		fitirf,obsirf,presentedirf,&
		residairf,residirf,stdresidirf,thetacheck)
		implicit none
		character(len=32),intent(in)::factorname(:)
		character(len=32),intent(in)::itemname(:)

		integer,intent(in)::numcatobs(:),unitirf
		logical,intent(in)::resid
		real(kind=8),intent(in)::fitirf(:,:),obsirf(:,:),&
			presentedirf(:),residairf(:,:),residirf(:,:),&
			stdresidirf(:,:),thetacheck(:,:)
	end subroutine printirfs
!
!	Print marginal distribution.
!	itemname contains item names.
!	numcatobs counts observed categories.
!	unitmargin provides the unit for output.
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
!	stdresidmarg is the asymptotic standard error for the residual frequency.
!	stdresidpmarg is the asymptotic standard error for the residual proportion.
!	stobsmarg is the standard error of the indicator.

	subroutine printmarginaldist(itemname,numcatobs,unitmargin,resid,&
		fitmarg,fitpmarg,obsmarg,obspmarg,presented,&
		relmarg,residamarg,residmarg,residpmarg,&
		seobsmarg,stdobsmarg,stdobspmarg,stdresidmarg,stdresidpmarg,stobsmarg)
		implicit none
		character(len=32),intent(in)::itemname(:)
		integer,intent(in)::numcatobs(:),unitmargin
		logical,intent(in)::resid
		real(kind=8),intent(in)::fitmarg(:),fitpmarg(:),obsmarg(:),&
			obspmarg(:),presented(:),relmarg(:),residamarg(:),residmarg(:),&
			residpmarg(:),seobsmarg(:),stdobsmarg(:),stdobspmarg(:),stdresidmarg(:),&
			stdresidpmarg(:),stobsmarg(:)
	end subroutine printmarginaldist

!
!	Print two-way cross products of item scores.
!	itemname contains item names.

!	unitmargins2 provides the unit for output.
!	resid indicates if residuals are displayed.
!	fitmargs2 is the fitted total cross products.
!	fitpmargs2 is the fitted average cross products.
!	obsmargs2 is the observed total cross products.
!	obspmargs2 is the observed average cross products.
!	presenteds2 counts weighted pairs of items presented.
!	residamargs2 is the adjusted residual.
!	residmargs2 is the residual for the total cross products.
!	residpmargs2 is the residual for the average cross products.
!	stdobsmargs2 is the standard deviation of the observed total cross products.
!	stdobspmargs2 is the asymptotic standard deviation of the observed 
!		average cross products.
!	stdresidmargs2 is the asymptotic standard error for the residual total cross products.
!	stdresidpmargs2 is the asymptotic standard error for the residual average cross products.


	subroutine printmarginals2(itemname,unitmargins2,resid,&
		fitmargs2,fitpmargs2,obsmargs2,obspmargs2,presenteds2,&
		residamargs2,residmargs2,residpmargs2,&
		stdobsmargs2,stdobspmargs2,stdresidmargs2,stdresidpmargs2)
		implicit none
		character(len=32),intent(in)::itemname(:)
		integer,intent(in)::unitmargins2
		logical,intent(in)::resid
		real(kind=8),intent(in)::fitmargs2(:,:),fitpmargs2(:,:),obsmargs2(:,:),&
			obspmargs2(:,:),presenteds2(:,:),residamargs2(:,:),residmargs2(:,:),&
			residpmargs2(:,:),stdobsmargs2(:,:),stdobspmargs2(:,:),&
			stdresidmargs2(:,:),stdresidpmargs2(:,:)
	end subroutine printmarginals2
	
!
!	Print two-way marginal distribution.
!	itemname contains item names.
!	numcatobs counts observed categories.
!	unitmargin provides the unit for output.
!	fitmarg2 is the fitted two-way marginal total.
!	fitpmarg2 is the fitted two-way marginal proportion.
!	obsmarg2 is the observed two-way marginal total.
!	obspmarg2 is the observed two-way marginal proportion.
!	presented2 counts weighted pairs of items presented.
!	residamarg2 is the adjusted residual.
!	residmarg2 is the residual for the two-way marginal total.
!	residpmarg2 is the residual for the two-way marginal proportion.
!	stdobsmarg2 is the standard deviation of the observed two-way marginal.
!	stdobspmarg2 is the asymptotic standard deviation of the observed 
!		two-way marginal fraction.
!	stdresidmarg2 is the asymptotic standard error for the residual frequency.
!	stdresidpmarg2 is the asymptotic standard error for the residual proportion.


	subroutine printmarginal2(itemname,numcatobs,unitmargin2,resid,&
		fitmarg2,fitpmarg2,obsmarg2,obspmarg2,presented2,&
		residamarg2,residmarg2,residpmarg2,&
		stdobsmarg2,stdobspmarg2,stdresidmarg2,stdresidpmarg2)
		implicit none
		character(len=32)::itemname(:)
		integer,intent(in)::numcatobs(:),unitmargin2
		logical,intent(in)::resid
		real(kind=8),intent(in)::fitmarg2(:,:),fitpmarg2(:,:),obsmarg2(:,:),&
			obspmarg2(:,:),presented2(:,:),residamarg2(:,:),residmarg2(:,:),&
			residpmarg2(:,:),stdobsmarg2(:,:),stdobspmarg2(:,:),&
			stdresidmarg2(:,:),stdresidpmarg2(:,:)	
	end subroutine printmarginal2
	


!
!	Print marginal distribution of weighted sum.
!	weightname contains sum name.
!	maxscore is the maximum score.
!	minscore is the minimum score.
!	unitmarginwtsum provides the unit for output.
!	resid indicates if residuals are displayed.
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

	subroutine printmargindistwtsum(weightname,maxscore,minscore,unitmarginwtsum,resid,&
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
		character(len=32)::weightname
		integer,intent(in)::maxscore,minscore,unitmarginwtsum
		logical,intent(in)::resid
		real(kind=8),intent(in)::fitcummargwtsum(minscore:maxscore),fitmargwtsum(minscore:maxscore),&
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
	end subroutine printmargindistwtsum

!	Print totals and averages of products of predictors and indicators of item categories.
!	itemname contains item names.
!	predname contains predictor names.
!	numcatobs counts observed categories.
!	unitmargpred provides the unit for output.
!	resid indicates if residuals are displayed.
!	fitppred is the  fitted average of products.
!	fitpred is the fitted total of products.
!	obsppred is the observed average of products.
!	obspred is the observed total of products.
!	presentedpred counts weighted items presented.
!	residapred is the adjusted residual.
!	residppred is the residual for the average of products.
!	residpred is the residual for the total of products.
!	stdobsppred is the standard error of the observed average of products.
!	stdobspred is the standard error of the observed 
!		total of products.
!	stdresidppred is the asymptotic standard error for the residual for average
!       of products.
!	stdresidpred is the asymptotic standard error for the residual for total
!       of products.

	subroutine printmarginpred(itemname,predname,numcatobs,unitpreditem,&
		resid,fitppred,fitpred,obsppred,obspred,presentedpred,&
		residapred,residppred,residpred,&
		stdobsppred,stdobspred,stdresidppred,stdresidpred)
		implicit none
		character(len=32),intent(in)::itemname(:)
		character(len=32),intent(in)::predname(:)
		integer,intent(in)::numcatobs(:),unitpreditem
		logical,intent(in)::resid
		real(kind=8),intent(in)::fitppred(:,:),fitpred(:,:),&
			obsppred(:,:),obspred(:,:),presentedpred(:),&
			residapred(:,:),residppred(:,:),residpred(:,:),&
			stdobsppred(:,:),stdobspred(:,:),stdresidppred(:,:),stdresidpred(:,:)
	end subroutine printmarginpred

!	Print totals and averages of products of transformed latent variables and indicators of item categories.
!	itemname contains item names.
!	skillname contains skill\r names.
!	numcatobs counts observed categories.
!	unitthetaitem provides the unit for output.
!	resid indicates if residuals are displayed.
!	fitptheta is the  fitted average of products.
!	fittheta is the fitted total of products.
!	obsptheta is the observed average of products.
!	obstheta is the observed total of products.
!	presentedtheta counts weighted items presented.
!	residatheta is the adjusted residual.
!	residptheta is the residual for the average of products.
!	residtheta is the residual for the total of products.
!	stdobsptheta is the standard error of the observed average of products.
!	stdobstheta is the standard error of the observed
!		total of products.
!	stdresidptheta is the asymptotic standard error for the residual for average
!       of products.
!	stdresidtheta is the asymptotic standard error for the residual for total
!       of products.


    subroutine printmargintheta(itemname,skillname,numcatobs,unitthetaitem,&
        resid,fitptheta,fittheta,obsptheta,obstheta,presentedtheta,&
        residatheta,residptheta,residtheta,&
        stdobsptheta,stdobstheta,stdresidptheta,stdresidtheta)
        implicit none
        character(len=32),intent(in)::itemname(:)
        character(len=32),intent(in)::skillname(:)
        integer,intent(in)::numcatobs(:),unitthetaitem
        logical,intent(in)::resid
        real(kind=8),intent(in)::fitptheta(:,:),fittheta(:,:),&
            obsptheta(:,:),obstheta(:,:),presentedtheta(:),&
            residatheta(:,:),residptheta(:,:),residtheta(:,:),&
            stdobsptheta(:,:),stdobstheta(:,:),stdresidptheta(:,:),stdresidtheta(:,:)
    end subroutine printmargintheta

!
!	Print weighted sum totals by item categories.
!	itemname contains item names.
!	wtname is the weight name.
!	numcatobs counts observed categories.
!	unitmarginwtitem provides the unit for output.
!	resid indicates if residuals are displayed.
!	fitpwtitem is the  fitted average of products.
!	fitwtitem is the fitted total of products.
!	obspwtitem is the observed average of products.
!	obsmargwtitem is the observed total of products.
!	presentedwtitem counts weighted items presented with items.
!	residawtitem is the adjusted residual.
!	residpwtitem is the residual for the average of products.
!	residwtitem is the residual for the total of products.
!	stdobspwtitem is the standard deviation of the observed average of products.
!	stdobswtitem is the asymptotic standard deviation of the observed 
!		total of products.
!	stdresidpwtitem is the asymptotic standard error for the residual for average
!       of products.
!	stdresidwtitem is the asymptotic standard error for the residual for total
!       of products.
	subroutine printmarginwtitem(itemname,wtname,numcatobs,unitmarginwtitem,&
		resid,fitpwtitem,fitwtitem,obspwtitem,obswtitem,presentedwtitem,&
		residawtitem,residpwtitem,residwtitem,&
		stdobspwtitem,stdobswtitem,stdresidpwtitem,stdresidwtitem)
		implicit none
		character(len=32),intent(in)::itemname(:)
		character(len=32),intent(in)::wtname
		integer,intent(in)::numcatobs(:),unitmarginwtitem
		logical,intent(in)::resid
		real(kind=8),intent(in)::fitpwtitem(:),fitwtitem(:),&
			obspwtitem(:),obswtitem(:),presentedwtitem(:),&
			residawtitem(:),residpwtitem(:),residwtitem(:),&
			stdobspwtitem(:),stdobswtitem(:),stdresidpwtitem(:),stdresidwtitem(:)
	end subroutine printmarginwtitem

!	Iteration progress reporting.
!	printprogstart indicates that starting iteration progress is to be printed to a file.
!	printprogstartstd indicates that starting iteration progress is to go to standard output.
!	printprog indicates that regular iteration progress is to be printed to a file.
!	printprogstd indicates that regular iteration progress is to go to standard output.
	subroutine progressreport(printprogstart,printprogstartstd,printprog,printprogstd)
		implicit none
		logical,intent(out)::printprog,printprogstart,printprogstartstd,printprogstd
	end subroutine  progressreport


!	Print reliability results for eap's.
!	comment is used to identify output.
!	factorname is variable names.
!	unitrel is unit.
!	covartheta is the covariance matrix of the latent vector.
!	coveap is the covariance matrix of the eaps.
!	meancovtheta is the average conditional covariance of the
!		latent vector given the observations.
!	meaneap is the average eap.
!	releap contains reliabilities.
	subroutine reloutput(comment,factorname,unitrel,covartheta,coveap,meancovtheta,meaneap,releap)
		implicit none
		character(len=*),intent(in)::comment
		character(len=32),intent(in)::factorname(:)
		integer,intent(in)::unitrel
		real(kind=8),intent(in)::covartheta(:,:),coveap(:,:),meancovtheta(:,:),meaneap(:),releap(:)
	end subroutine reloutput

!	Find estimated reliabilities of eaps of latent vectors.
!	Input:

!	covtheta is the array of conditional covariance matrices.
!	meantheta is the array of conditional means.

!	Output:

!	obsweight provides observation weights.
!	coveap is the covariance matrix of the eaps.
!	covartheta is the covariance matrix of the latent vector.
!	meancovtheta is the average conditional covariance of theta given the observations.
!	meaneap is the average eap.
!	releap is the reliability for each element of the latent vector.

	subroutine reltheta(covtheta,meantheta,obsweight,coveap,covartheta,meancovtheta,meaneap,releap)
		implicit none
		real(kind=8),intent(in)::covtheta(:,:,:),meantheta(:,:),obsweight(:)
		real(kind=8),intent(out)::coveap(:,:),covartheta(:,:),meancovtheta(:,:),meaneap(:),releap(:)
	end subroutine reltheta

!
!	Find estimated reliabilities of eap's of transformed latent vectors.
!	Input:
!	covartheta is the covariance matrix of the latent vector.
!	coveap is the covariance matrix of the eaps.
!	lintran is the transformation for the latent vector.
!	meancovtheta is the average conditional covariance of the latent vector given the observations.
!	meaneap is the average eap.

!	covarthetaskill is the covariance matrix of the transformed latent vector.
!	coveapskill is the covariance matrix of the eaps for the transformed latent vector.
!	meancovthetaskill is the average conditional covariance of the transformed latent vector
!		given the observations.
!	meaneapskill is the average eap for the transformed latent vector. 
!	releapskill is the reliability of each
!		element of the transformed latent vector.
	subroutine relthetaskill(covartheta,coveap,lintran,meancovtheta,meaneap,&
		covarthetaskill,coveapskill,meancovthetaskill,meaneapskill,releapskill)	
		implicit none
		real(kind=8),intent(in)::covartheta(:,:),coveap(:,:),&
			lintran(:,:),meancovtheta(:,:),meaneap(:)
		real(kind=8),intent(out)::covarthetaskill(:,:),coveapskill(:,:),&
			meancovthetaskill(:,:),meaneapskill(:),releapskill(:)
	end subroutine relthetaskill

!	print individual residual data.
!	comment is used for identification.
!	resname gives names of residual elements.
!	unit unitres is the unit for output.
!	useid indicates if ids are put out.
!	covres is residual covariance matrix.
!   fit gives fitted values.
!   observed gives observed values.
!	res  gives residual.
!   sdres gives estimated residual standard errors.
!   zres gives normalized residual.

    subroutine resoutput(comment,resname,id,unitres,useid,&
        covres,fit,observed,res,sdres,zres)
        implicit none
        character(len=*),intent(in)::comment
        character(len=32),intent(in)::resname(:)
        character(len=32),intent(in),optional::id(:)
        integer,intent(in)::unitres
        logical,intent(in)::useid
        real(kind=8),intent(in)::covres(:,:,:),fit(:,:),observed(:,:),res(:,:),sdres(:,:),zres(:,:)
    end subroutine resoutput

!	Find individual residuals for weighted sums.
!	catobsrange provides ranges of underlying item categories that correspond to
!   observed item categories.
!	Responses are in dat.
!	numcat provides ranges for underlying observations.
!	numcatobs provides ranges for observations.

!	beta is the parameter vector.
!	lintran is the linear transformation for the latent vector.
!	postdensity contains posterior weights.
!	theta contains posterior quadrature points.
!	wtsum is a matrix used to compute linear combinations of the observations.

!   covreswt is the estimated covariance matrix of reswt.
!   fitwt is the estimated mean of the weighted sum.
!   obswt is the observed weighted sum.

!   reswt is the raw individual residual for the weighted sum.
!   sdreswt is the estimated standard error of the elements of reswt.
!   zreswt is the ratio of elements of reswt and sdreswt.

    subroutine reswtsum(catobsrange,dat,numcat,numcatobs,beta,lintran,postdensity,theta,wtsum,&
        covreswt,fitwt,obswt,reswt,sdreswt,zreswt)
        implicit none
        integer,intent(in)::catobsrange(:,:),dat(:,:),numcat(:),numcatobs(:)
   

        real(kind=8),intent(in)::beta(:),lintran(:,:),postdensity(:,:),theta(:,:,:),wtsum(:,:)
        real(kind=8),intent(out)::covreswt(:,:,:),fitwt(:,:),obswt(:,:),reswt(:,:),sdreswt(:,:),zreswt(:,:)
    end subroutine reswtsum

!	Find individual residuals for weighted sums.
!	catobsrange provides ranges of underlying item categories that correspond to
!	choices provides number of choices per item.
!	distmap connects distractors to item scores.
!   observed item categories.
!	Response scores are in dat.
!	Raw responses are in distdat.
!	numcat provides ranges for underlying observations.
!	numcatobs provides ranges for observed scores.

!	beta is the parameter vector.
!	distprob gives conditional distractor probabilities.
!	lintran is the linear transformation for the latent vector.
!	postdensity contains posterior weights.
!	theta contains posterior quadrature points.
!	wtsumdist is a matrix used to compute linear combinations of the observations.

!   covreswtdist is the estimated covariance matrix of reswtdist.
!   fitwtdist is the estimated mean of the weighted sum.
!   obswtdist is the observed weighted sum.

!   reswtdist is the raw individual residual for the weighted sum.
!   sdreswtdist is the estimated standard error of the elements of reswtdist.
!   zreswtdist is the ratio of elements of reswtdist and sdreswtdist.

	subroutine reswtsumdist(catobsrange,choices,dat,distdat,&
        distmap,numcat,numcatobs,beta,distprob,lintran,postdensity,theta,wtsumdist,&
		covreswtdist,fitwtdist,obswtdist,reswtdist,sdreswtdist,zreswtdist)
		implicit none
		integer,intent(in)::catobsrange(:,:),choices(:),dat(:,:),distdat(:,:),distmap(:),numcat(:),numcatobs(:)
	

		real(kind=8),intent(in)::beta(:),distprob(:),lintran(:,:),postdensity(:,:),theta(:,:,:),wtsumdist(:,:)
		real(kind=8),intent(out)::covreswtdist(:,:,:),fitwtdist(:,:),obswtdist(:,:),reswtdist(:,:),sdreswtdist(:,:),zreswtdist(:,:)
	end subroutine reswtsumdist



!       Save alpha and cholnhess on unit specified by unitalpha
!       alpha is location of the maximum of the conditional probability
!		of a response given the latent vector.
!       cholnhess is the modified Cholesky decomposition of the
!      		negative hessian of the log conditional probability with respect
!      		to the latent vector at the maximum.
	subroutine savealpha(unitalpha,alpha,cholnhess)
		implicit none

		integer,intent(in)::unitalpha
		real(kind=8),intent(in)::alpha(:,:),cholnhess(:,:,:)
	end subroutine savealpha


!
!	Standard error for entropy per item under standard sampling.
!	loglik is log likelihood.
!	obsweight is  observation weights.
!	prob is array of observed probabilities.
!	totalitems is weighted sum of items presented.
	real(kind=8) function sdent(loglik,obsweight,prob,totalitems)
		implicit none
		real(kind=8),intent(in)::loglik,obsweight(:),prob(:),totalitems
	end function sdent



!	Standard error for entropy per item under complex sampling.
!	This version is designed for random sampling of psu's within strata.
!	Sampling is with replacement.
!	npsu is the number of psus's per stratum.
!	nstratum is the number of strata.
!	stratum is observation stratum.
!	psu is observation psu.
!	stratify is .true. for stratified sampling.
!	usepsu is .true. for primary sampling units within strata.
!	obsweight is  observation weights.
!	prob is array of observed probabilities.
!	totalitems is weighted sum of items presented.
	real(kind=8) function sdent_complex(npsu,nstratum,psu,stratum,stratify,usepsu,obsweight,prob,totalitems)
		implicit none
		integer,intent(in)::nstratum
		integer,intent(in),optional::npsu(:),psu(:),stratum(:)
		logical,intent(in)::stratify,usepsu
		real(kind=8),intent(in)::obsweight(:),prob(:),totalitems
	end function sdent_complex

!	Find standard errors of distractor frequencies and relative frequencies.
!   choices provides the number of choices for each item.
!   dat provides item data.
!   distdat provides distractor data.
!   distmap maps distractors to item scores.
!	npsu is the number of psus's per stratum.
!	nstratum is the number of strata.
!   numcatobs provides the number of categories for item scores.
!	psu indicates the psu of an observation.
!	stratum indicates the stratum of an observation.
!	complx indicates complex sampling.
!	stratify is .true. for stratified sampling.
!	usepsu is .true. for primary sampling units within strata.

!   distfreq provides distractor frequencies.

!   distprob provides conditional distractor probabilities.
!   obsweight gives observation weights.
!   scorefreq provides item score frequencies.
!	sdfreq gives standard errors for frequencies.
!	sdprob gives standard errors for probabilities.
    subroutine sdfreqs(choices,dat,distdat,distmap,npsu,nstratum,numcatobs,obsweight,&
        psu,stratum,complx,stratify,usepsu,distfreq,distprob,scorefreq,sdfreq,sdprob)
        implicit none
        integer,intent(in)::choices(:),dat(:,:),distdat(:,:),distmap(:),npsu(:),nstratum,numcatobs(:),&
            psu(:),stratum(:)
        logical,intent(in)::complx,stratify,usepsu
        real(kind=8),intent(in)::distfreq(:),distprob(:),obsweight(:),scorefreq(:)
        real(kind=8),intent(out)::sdfreq(:),sdprob(:)
    end subroutine sdfreqs

        
!	set parameters for eap output.
!	dimtrans is the transformation dimension.
!	dimwtsum is the dimension of the desired linear combination of observations.
!	eapmask indicates the items to include.
!	altbeta indicates the beta vector to use.
	subroutine seteapoutput(dimtrans,dimwtsum,dimwtsumdist,eapmask,altbeta)
		implicit none
		integer,intent(out)::dimtrans,dimwtsum,dimwtsumdist
		logical,intent(out)::eapmask(:)
		real(kind=8),intent(inout)::altbeta(:)
	end subroutine seteapoutput

!	Set output options.  See unitdef for definitions.

	subroutine setoutput()
	end subroutine setoutput
	
!	Find weight limits for each item.
!	numcatobs gives number of observed categories per item.
!	weights are category weights.
!	maxw gives the maximum weights per item.
!	minw gives the minimum weights per item.
	subroutine setwtlimits(numcatobs,weights,maxw,minw)
		implicit none
		integer,intent(in)::numcatobs(:),weights(:)
		integer,intent(out)::maxw(:),minw(:)
	end subroutine setwtlimits
!
!	sortstratum is used to sort the strata and/or psus and then to
!	recode them so that strata are numbered from 1 to nstratum,
!	the number of strata and psus within strata are numbered from
!	1 to the number of psus within the stratum.
!	psu is the array of psu codes.
!	stratum is the array of stratum codes.
!	stratify is .true. for stratified sampling.
!	usepsu is .true. if psu's are used.
	subroutine sortstratum(psu,stratum,stratify,usepsu,nstratum)
		implicit none
		integer,intent(inout)::psu(:),stratum(:)
		logical,intent(in)::stratify,usepsu
		integer,intent(out)::nstratum
	end subroutine sortstratum

!	Generic procedure for residuals.
!	npsu is the number of psus's per stratum.
!	nstratum is the number of strata.
!	psu indicates the psu of an observation.
!	stratum indicates the stratum of an observation.
!	complx indicates complex sampling.
!   presence indicates if observation present.
!   irel is .true. if reliability is desired.
!   resid is .true. if residuals are desired.
!   stats is .true. if summary statistics are desired.
!	stratify is .true. for stratified sampling.
!	usepsu is .true. for primary sampling units within strata.
!	eacovgaminv_louis is the inverse of the Louis information plus the
!		constraint component of the negative hessian.
!	fitted is the estimated expectation for an observation.
!	gradsdes provides gradients for observations
!	observed is the observed value for an observation.
!   obsweight is the observation weight.
!	tolres is the residual tolerance.
!   aresid is the adjusted residual.
!   covobs is the estimated covariance matrix of an observation.
!   covobsave is the estimated covariance matrix of obsave.
!   covobssum is the estimated covariance matrix of obssum.
!   covres is the estimated covvariance matrix of a residual vector.
!   covresave is the estimated covariance matrix of resave.
!   covressum is the estimated covariance matrix of ressum.
!   fitsum is the expected sum.
!   fitave is the expected average.
!   obsave is the observed average.
!   obssum is the observed sum.
!   rel is a reliability estimate.
!   resave is the average residual.
!   ressum is the sum of residuals.
!   sdobs is the estimated standard deviation of an observation.
!   sdobsave is the estimated standard deviation of obsave.
!   sdobssum is the estimated standard deviation of obssum.
!   sdres is the estimated standard deviation of a residual vector
!   sdresave is the estimated asymptotic standard deviation of
!       ressum.
!   sdressum is the estimated asymptotic standard deviation of resave.
!   totsum is the sum of observations weight for present observations.



    subroutine stdresid(npsu,nstratum,&
        psu,stratum,complx,irel,presence,resid,stats,stratify,usepsu,&
        eacovgaminv_louis,fitted,gradsdes,&
        observed,obsweight,tolres,&
        aresid,covobs,covobsave,covobssum,covres,covresave,covressum,fitave,fitsum,obsave,obssum,rel,resave,&
        ressum,sdobs,sdobsave,sdobssum,sdres,sdresave,sdressum,totsum)

        implicit none
        integer,intent(in)::npsu(:),nstratum,psu(:),stratum(:)
        logical,intent(in)::complx,irel,presence(:),resid,stats,stratify,usepsu
        real(kind=8),intent(in)::eacovgaminv_louis(:,:),&
            fitted(:,:),gradsdes(:,:),&
            observed(:,:),obsweight(:),tolres
        real(kind=8),intent(out)::aresid(:),covobs(:,:),covobsave(:,:),covobssum(:,:),&
            covres(:,:),covresave(:,:),covressum(:,:),&
            fitave(:),fitsum(:),&
            obsave(:),obssum(:),&
            rel(:),resave(:),ressum(:),&
            sdobs(:),sdobsave(:),sdobssum(:),sdres(:),sdresave(:),sdressum(:),totsum
    end subroutine stdresid

end interface

!	Input file and input format.  Optional title
character(len=256)::filename,fileformat
character(len=80)::title
!	Factor names.  Ids. Item names.  Predictor names.  scale names. skill names.  Transformations names. Weight sum names.
character(len=32),allocatable::factorname(:)
character(len=32),allocatable::id(:),itemname(:)
character(len=32),allocatable::predlabels(:),predname(:)
character(len=32),allocatable::scalename(:),skillname(:)
character(len=32),allocatable::transname(:)
character(len=32),allocatable::weightnames(:),weightnamesdist(:),wtname(:),wtnamedist(:)
!   Gamma names.  Beta names.
character(len=64),allocatable::betaname(:),paramname(:)
!	minobs is the minimum observation count.  It must be at least 2.
integer,parameter::minobs=2
!	al is used to detect allocation errors.
!	constdim is the constraint dimension.
!	constdim1 is the number of linear constraints not treated via least squares.
!	dimdesign is the dimension of the design matrix.
!	dimitemmatrix is the dimension of the item matrix.
!	dimlatin is the dimension of the latent vector.
!	dimlatout is the dimension of the transformed latent vector.
!	dimno is a dimension.
!	dimtrans is dimension of transformation of latent vector.
!	dimwtsum is the dimension of the desired linear combination of observations.
!	dimwtsumdist is the dimension of the desired linear combination of distractor observations.
!	intconst is the number of linear constraints on intercepts.
!	intdesdim is the number of columns in the matrix defining the model for
!	intercept parameters.
!	irfcount is the desired number of points for item response functions.  
!	linconst is the number of linear constraints on linear components of the density.
!	lindesdim is the number of columns in the matrix defining the model for
!	the linear component of the log latent density.
!	maxit is the maximum number of iterations.
!	maxita is the maximum number of iterations to find an alpha.
!	maxitb is the maximum number of iterations used to find the step size.
!	maxite is the maximum number of iterations with the evenly spaced option.
!	maxscore is the minimum score.
!	maxscore is the maximum score.
!	modeldim is the model dimension.
!	modeldistdim is the model dimension with distractors.
!	ncat is the total number of underlying categories.
!	ncatobs is the total number of observed categories.
!	nexternal is the number of external variables.
!	ngrid is the number of points in the grid used in selection of quadrature points.
!	nitems is the number of items. 
!	nobs is the number of observations.
!	nobsstart is the initial sample size for starting calculations.
!	nparam is the dimension of beta.
!	npred is the number of predictors.
!	npredlist is the number of vectors in predlist.
!	npredtot is the sum of npred and nexternal.
!	nquad is the number of quadrature points.
!	nreduce is the dimension of the reduced version of beta.
!	numberscales counts scales per weighted sum.
!	nstratum is the number of strata.
!   numchoices is the total number of choices for item responses.
!	numweights is the number of weighted sums to examine for distributions.
!	numweightsdist is the number of weighted sums to examine for distractor distributions.
!	obs is an observation.
!	quadconst is the number of linear constraints on quadratic components of the density.
!	quaddesdim is the number of columns in the matrix defining the model for
!	the quadratic component of the log latent density.
!	row is an array element.
!	slopeconst is the number of linear constraints on slopes.
!	slopedesdim is the number of columns in the matrix defining the model for
!	slope parameters.

integer::al
integer::constdim,constdim1
integer::dimdesign,dimitemmatrix,dimlatin,dimlatout,dimno,dimtrans,dimwtsum,dimwtsumdist
integer::intconst,intdesdim,irfcount
integer::linconst,lindesdim
integer::maxit,maxita,maxitb,maxite,maxscore,minscore,modeldim,modeldistdim
integer::ncat,ncatobs,nexternal,ngrid,nitems,nobs,nobsstart,nparam,npred,npredlist,npredtot
integer::nquad,nreduce,numberscales,nstratum,numchoices,numweights,numweightsdist
integer::obs
integer::quadconst,quaddesdim
integer::row
integer::slopeconst,slopedesdim
	
!	betamap maps beta elements to items, categories, and dimensions.
!	catobsrange has ranges of underlying categories for observed categories.
!	choices is the number of choices for an item.
!	coord specifies coordinates of quadrature points for regular integration.
!	dat contains data.
!	dimpoints is the number of integration points for each dimension.
!	distdat is the distractor matrix.
!   distmap maps distractors to item scores.
!	intdim provides dimensions of intercept matrix.
!	maxw gives minimum item weights.
!	minw gives maximum item weights.
!	npsu lists the number of psu's in each stratum.

!	numcat has observed categories per item.
!	numcatobs has underlying categories per item.
!	numcodes provides the number of recodes per item.
!	numdistcodes is the number of distractor recodes per item.
!	psu is the psu number of the observation.
!	rbetamap is the reduced map for beta.
!	recodetab is the table of recodes for each item.
!	recodedisttab is the table of distractor recodes for each item.
!	skillnumbers specifies skills that correspond to items in a simple-structure case.
!	slopedim provides dimensions of slope matrix.
!	stratum is the stratum number of the observation.
!	weights is for integer item weights.
!	weightsdist is for integer item distractor weights.
integer,allocatable::betamap(:,:)
integer,allocatable::catobsrange(:,:),choices(:),coord(:,:)
integer,allocatable::dat(:,:),dimpoints(:),distdat(:,:),distmap(:)
integer,allocatable::intdim(:)
integer,allocatable::maxw(:),minw(:)
integer,allocatable::numcat(:),numcatobs(:),numcodes(:),numdistcodes(:),npsu(:)
integer,allocatable::psu(:)
integer,allocatable::rbetamap(:,:),recodetab(:,:),recodedisttab(:,:)
integer,allocatable::skillnumbers(:),slopedim(:,:),stratum(:)
integer,allocatable::weights(:,:),weightsdist(:,:)

!	allmeans if .true. if all means and covariances of underlying latent vectors are to be supplied.	
!	altpost is .true. if an alternative posterior distribution is desired.

!	complx is true if either usepsu or stratify or weight is true.
!	constguess indicates if guessing parameters are all the same.
!	cross is .true. for the cross-polytope version of quadrature.
!	custom is for a special linear transformation.
!	custompoints is .true. only if special points are used for item response functions.
!	distract is .true. if distractors are read.
!	error is the error flag for a bad covariance matrix.
!	even is for even spacing.
!	fixquad indicates that the quadratic component of the logarithm of the latent density is fixed.
!	for all predictor values.
!	gausshermite is .true. if Gauss-Hermite quadrature is used.
!	grid is .true. if quadrature points are on a grid.
!	initparam is used to decide if default gamma and offset are needed.
!	meanid indicates that labels for vectors are used for printing of means and covariances of underlying latent vectors.
!	mp is used for maximum posterior likelihood.

!	noquad indicates that the quadratic terms are 0 for the log density.
!	normal is for a multivariate normal distribution of the latent vector.
!	nr is for Newton-Raphson.  The alternative uses an approximation to Newton-Raphson based on the Louis
!		approximation.

!	proport is true if constraint sum of squares is proportional to sum of observation weights.
!	recode is used for item recoding.
!	recodedist indicates if distractor recodes are present.
!	simplex is .true. if a simplex is used for quadrature points.
!	specialtrans specifies a custom transition matrix.
!	stratify indicates the presence of stratification if true.
!	twostages is used for two rounds of iterations.
!	If true, useid causes reading of an examinee identification code.
!	usepsu indicates the presence of psu's if true.
!	weight indicates whether observation weights are used.
logical::allmeans,altpost
logical::complx,constguess,cross,custom,custompoints
logical::distract
logical::equalpoints,error,even
logical::fixquad,fullgrid
logical::gausshermite,grid
logical::initparam,meanid,mp
logical::noquad,normal,nr
logical::proport,recode,recodedist
logical::simplex,specialtrans,stratify
logical::twostages,useid,usepsu,weight
!	betweenitem indicates a between-item model.
!	catmap is for category maps.
!	eapmaskis used to determine items used in eap computations.
!	fixdiag indicates that the diagonal terms of the quadratic component of the
!	logarithm of the latent density are fixed.
!	fixedguess specifies that guessing parameters have a specified value.
!	guess specifies if a guessing parameter is used.
!	independentf indicates which latent vector coordinates are conditionally
!	independent given the independent variables.
!	no_lin indicates that the base linear terms are 0 for the log density.
!	nominal specifies that a nominal model is used for the item.
!	predlinmap is map for predictors for linear component.
!	predquadmap is map for predictors for quadratic component.
!   presence is for existence of weighted sum.
!	rasch specifies that all slope parameters for a skill are equal.
!	raschslope1 specifies that all slope parameters associated with a skill are 1.


!	specialitemint specifies an item which needs special treatment for intercept parameters.
!	specialitemslope specifies an item which needs special treatment for slope parameters.
logical,allocatable::betweenitem(:),catmap(:),eapmask(:)
logical,allocatable::fixdiag(:),fixedguess(:),guess(:),independentf(:),no_lin(:),nominal(:)
logical,allocatable::predlinmap(:,:),predquadmap(:,:,:),presence(:)
logical,allocatable::rasch(:),raschslope1(:),specialitemint(:),specialitemslope(:)

!	akent is the Akaike adjustment of ent.
!	akentdist is the Akaike adjustment of entdist.
!	changemin is  used as a criterion for minimum change.
!	ent is the information per item per unit observation.
!	entdist is the information per item per unit observation with distractors added.
!	gap is used to specify the selection of observations for
!	starting sampling.
!	ghent is the Gilula-Haberman adjustment of ent.
!	ghentdist is the Gilula-Haberman adjustment of entdist.
!	ghent_complex is the Gilula-Haberman adjustment of ent for complex sampling.
!	ghentdist_complex is the Gilula-Haberman adjustment of entdist for complex sampling.
!	kappa is a maximum step size.
!	loglik is log likelihood.
!	loglikdist is log likelihoo with distractors.
!	maxdalpha is the maximum permitted alpha change.
!	presentedwtsum counts weighted items presented
!		for the weighted sum.
!	pspread indicates the width of the evenly spaced grid.
!	sdentropy is estimated asymptotic standard error of ent.
!	sdentropydist is estimated asymptotic standard error of entdist.
!	sdentropy_complex is estimated asymptotic standard error of ent for complex sampling.
!	sdentropydist_complex is estimated asymptotic standard error of entdist for complex sampling.
!	tau is a maximum shrinkage of a step size.
!	tol is the tolerance for the information change per item per effective observation.
!	tola is the tolerance for the posterior change for an item.
!   tolc is the criterion for appropriate maximum ratio of diagonal terms in modified Cholesky decomposition.
!	tole is the tolerance used with the evenly spaced option.
!	tolres is the tolerance for standard errors of residuals.
!	tolsing is the tolerance for the modified Cholesky decomposition.
!	totalitems is weighted total number of items presented.
!	totalobs is weighted total number of observations presented.
real(kind=8)::akent,akentdist,changemin,ent,entdist
real(kind=8)::gap,ghent,ghentdist,ghent_complex,ghentdist_complex
real(kind=8)::kappa,loglik,loglikdist,maxdalpha
real(kind=8)::presentedwtsum,pspread
real(kind=8)::sdentropy, sdentropydist,sdentropy_complex,sdentropydist_complex
real(kind=8)::tau,tol,tola,tolc,tole,tolres,tolsing
real(kind=8)::totalitems,totalobs,totsum


!	alpha is a matrix of locations of maxima of posterior densities.
!	altalpha is an alternative matrix of locations of maxima of posterior densities.
!	altbeta is an alternative beta vector.
!	altcholnhess is an alternative array of modified Cholesky decompositions for each observation.
!	altpostdensity is the alternate array of posterior densities at adaptive quadrature points.	
!	alttheta is the alternate array of adaptive quadrature points.
!	beta is an alternative version of beta for EAP calculations.
!	catscore contains category scores.
!	cholnhess is an array of modified Cholesky decompositions for each observation.
!	constmat is the constraint matrix.
!	constvec is the constraint vector.
!	covmeanscale is the covariance matrix of the eaps for scale means.
!	covmeantrans is the covariance matrix of the eaps for transformations of latent vectors.
!	covartheta is the covariance matrix of the latent vector.
!	covarthetaskill is the covariance matrix of the transformed latent vector.
!	covartrans is the covariance matrix of the transformations of latent vector.
!	covarwt is the covariance matrix of the weighted sum.
!	coveap is the covariance matrix of the eaps.
!	coveapskill is the covariance matrix of the eaps for skills.
!	coveaptrans is the covariance matrix of the eaps for transformations.
!	coveapwt is the covariance matrix of the eaps for weighted sums.
!   covmeanobsscale is the residual component of scaleobscov.
!	covnewtheta is the array of transformed conditional covariance matrices.
!   covobs is  a covariance matrix of an observation.
!   covobsave is a covariance matrix of an observed average.
!   covobssum is a covariance matrix of an observed sum.
!   covres is a covariance matrix of a residual.
!   covresave is a covariance matrix of an average residual.
!   covressum is a covariance matrix of a residual sum.
!   covreswt is the estimated covariance matrix of reswt.
!   covreswtdist is the estimated covariance matrix of reswtdist.
!	covs are covariance matrices of underlying latent vectors.
!	covscale is the covariance matrix of the eap for the scale scores.
!	covsum is the conditional covariance of wtsum.
!	covtheta is the array of conditional covariance matrices.
!	covtrans is the array of conditional covariance matrices of transformations of latent vectors.

!	cweight is used to modify quadrature weights.
!	design is a design matrix.
!   distfreq gives frequencies of distractors.
!   distprob contains conditional distractor probabilities.
!	eacovbeta is the estimated asymptotic covariance matrix of beta.
!	eacovbeta_complex is the estimated asymptotic covariance matrix of beta based on complex sampling.
!	eacovbeta_louis is the estimated asymptotic covariance matrix of beta based on the
!		Louis approximation.
!	eacovbeta_sandwich is the estimated asymptotic covariance matrix of beta based on the sandwich
!		approach.
!	eacovgam is the estimated asymptotic covariance matrix of gamma.
!	eacovgam_complex is the estimated asymptotic covariance matrix of gamma based on complex sampling.
!	eacovgam_louis is the estimated asymptotic covariance matrix of gamma based on the
!		Louis approximation.
!	eacovgam_sandwich is the estimated asymptotic covariance matrix of gamma based on the sandwich
!		approach.
!	eacovgaminv is the inverse of nhessdesplus.
!	
!	eacovgaminv_louis is the inverse of the information plus the constraint component of the 
!		negative hessian.
!	extvar is used to store external variables.
!   fitave is a fitted average.
!	fitcummargwtsum is the fitted cumulative marginal
!		frequency distribution.
!	fitirf is for unconditional item response functions.
!	fitmarg is for fitted marginal frequencies.
!	fitmargs2 is for fitted total cross products of item scores.
!	fitmargwtsum is the fitted marginal frequency distribution.
!	fitmarg2 is for fitted two-way marginal frequencies.
!	fitpcummargwtsum is the fitted cumulative marginal
!		probability distribution.
!	fitpmarg is for fitted marginal proportions.
!	fitpmargs2 is for fitted average cross products of item scores.
!	fitpmargwtsum is the fitted marginal probability distribution.
!	fitpmarg2 is for fitted two-way marginal proportions.
!	fitppred is the average of fitted products of categority indicators and predictors.
!	fitptheta is the average of fitted product of item indicator and latent vectors.
!	fitpwtitem is the average of fitted products of category indicators and weighted sums.
!   fitsum is a fitted sum.
!	fittheta is the total of fitted products of item indicator and latent vectors.

!   fitwt is the individual weighted sum for item scores.
!   fitwtdist is the individual weighted sum for item distractors.
!	fitwtitem is the total of fitted products of item indicator and weighted sum.
!	gamma is a parameter vector such that beta=design*gamma.
!	grad is a gradient in the beta scale.
!	graddes is a gradient in the gamma scale.
!	grads is a matrix of gradients per observation in the beta scale.
!   guessres gives the residual total.
!   guessresa is the adjusted residual.

!	indvar is a matrix of independent variables.
!	informgam is the estimated information matrix for gamma.
!	informgam_complex is the estimated information matrix for gamma for complex sampling.
!	itemmatrix is the basic design matrix without parameter restrictions.
!	lintran is the linear transformation of the latent vector.
!   meancovobsscale is the EAP component of scaleobscov.
!	meancovscale is the average conditional covariance of the scale scores          
!	meancovtheta is the average conditional covariance of the latent vector given the observations.
!	meancovthetaskill is the average conditional covariance of the transformed latent vector given the observations.
!	meancovwt is the average conditional covariance of the weighted sum given the observations.
!   meanscale is the eap of the scales scores.
!	meaneap is the average eap.
!	meaneapskill is the average eap for a skill.
!	meaneaptrans is the average eap for a transformation.
!	meaneapwt is the average eap for a weighted sum.
!	meannewtheta is the array of transformed conditional means.
!	meansum is the conditional mean of wtsum.
!	means are means of underlying latent vectors.
!	meantheta is the array of conditional means.
!	meantrans is the array of conditional means of transformations.
!	nhess is a hessian in the beta scale.
!	nhessdes is a hessian in the gamma scale.
!	nhessdeschol is the modified Cholesky decomposition of nhessdes.
!	nhessdesplus is nhessdes plus the component for constraints.
!   obsave is an observed average.
!	obscummargwtsum is the observed cumulative marginal
!		frequency distribution.
!   observescale is the observed scale score.
!	obsirf is for conditional item response functions.
!	obsmarg is for observed marginal frequencies.
!	obsmargs2 is for totals of cross products of item scores.
!	obsmargwtsum is the observed marginal frequency distribution.
!	obsmarg2 is for observed two-way marginal frequencies.
!	obspcummargwtsum is the observed cumulative marginal
!		probability distribution.
!	obspmarg is for observed marginal proportions.
!	obspmargs2 is for average cross products of item scores.
!	obspmargwtsum is the observed marginal
!		probability distribution.
!	obspmarg2 is for observed two-way marginal proportions.
!	obsppred is the observed average products of category indicators and predictors.
!	obspred is the observed total of products of category indicators and predictors.
!	obspwtitem is the observed average products of category indicators and weighted sum.
!	obsptheta is the observed average product of item indicator and latent vectors.
!	obstheta is the observed total product of item indicator and predictors.
!   obsscalerel is the reliability of the observed scale scores.
!   obssum is an observed sum.
!	obsweight is a vector of observation weights.
!	obsweightsel is a vector of selected observation weights.
!   obswt is the individual weighted sum of item scores.
!   obswtdist is the individual weighted sum of distractors.
!	obswtitem is the observed total products of item indicator and weighted sum.
!	offset is an offset vector.
!	offsettran is the offset vector corresponding to the transition matrix.
!   plower is array of lower probabilities.
!	postdensity is the array of posterior densities at adaptive quadrature points.
!	predlist is the predictor list.
!	presented counts items presented to examinees.
!	presented2 and presenteds2 count pairs of items presented to examinees.
!	presentedirf is for items presented for item response functions.
!	presentedpred counts items presented to examinees for predictor analysis.
!	presentedtheta counts weighted items presented.
!	presentedwtitem counts weighted items presented with defined weighted sums.
!	prob contains item probabilities.
!   probdist contains distractor probabilities.
!   pupper is array of upper probabilities.
!	quadpoint is the array of quadrature points.
!	quadweight is the vector of weights of quadrature points.
!	rdesign is the reduced design matrix.
!	releap is the reliabilities of eap's for the latent vector.
!	releapskill is the reliabilities of eap's for the transformed latent vector.
!	relmarg is for reliabilities of individual category indicators.
!	reltrans is for reliabilities of eap's for transformations.
!   resave is an average residual.
!	residacummargwtsum is the adjusted residual for the cumulative
!		marginal distribution.
!	residairf is for adjusted residuals for item response functions.
!	residamarg is adjusted residuals for marginal distributions.
!	residamargs2 is for adjusted residuals for cross products of item scores.
!	residamargwtsum is the adjusted residual for the marginal
!		distribution.
!	residamarg2 is adjusted residuals for two-way marginal distributions.
!	residapred is adjusted residuals for total of products of category indicators and predictors.
!	residatheta is the adjusted residual of average of products of category indicators and transformed latent variables.
!	residawtitem is the adjusted residual of total of products of category indicators and predictors.
!	residcummargwtsum is the residual for the cumulative marginal
!		frequency distribution.
!	residirf is for residuals for item response functions.
!	residmarg is residuals for marginal frequencies.
!	residmargs2 is for residuals for totals of cross products of item scores.
!	residmargwtsum is the residual for the marginal frequency
!		distribution.
!	residmarg2 is residuals for two-way marginal frequencies.
!	residpcummargwtsum is the residual for the cumulative marginal
!		 probability distribution.
!	residpmarg is residuals for marginal proportions. 
!	residpmargs2 is for residuals for averages of cross products of item scores.
!	residpmargwtsum is the residual for the marginal probability
!		distribution.
!	residpmarg2 is residuals for two-way marginal proportions.
!	residppred is the residual average of products of category indicators and predictors.
!	residpred is the residual total of products of category indicators and predictors.
!	residpwtitem is the residual average of products of category indicators and weighted sums.
!	residptheta is the residual of average of products of category indicators and transformed latent variables.
!	residtheta is the residual of total of products of category indicators and transformed latent variables.
!	residwtitem is the residual total of products of category indicators and weighted sums.
!   ressum is a residual sum.
!   reswt is the raw individual residual for the weighted sum.
!   reswtdist is the raw individual residual for the weighted distractor sum.
!   scalecov is the covariance matrix of the EAP scale scores.
!   scalemean is the mean of the EAP scale scores.
!   scaleobscov is the covariance matrix of the observed scale scores.
!   scaleobsmean is the mean of the observedn scale scores.
!   scaleobsrel is the reliability of the observed scale scores.
!   scalerel is the reliability of the EAP scale scores.
!   scorefreq gives frequencies of item scores.
!	scores are item scores.
!	sdfreq gives standard errors for distractor frequencies.
!   sdobs is standard deviation of an observed residual.
!   sdobsave is the standard deviation of an average residual.
!   sdobssum is the standard deviation of a residual sum.
!	sdprob gives standard errors for distractor probabilities.
!   sdreswt is the estimated standard error of the elements of reswt.
!   sdreswtdist is the estimated standard error of the elements of reswtdist.
!	seobsmarg is standard errors for indicators of categories.
!	setguess is the fixed setting of a guessing parameter.
!	setslope helps initialize slope parameters.
!   stdguessres is the standard error of guessres.
!	stdobscummargwtsum is the standard deviation of the observed
!		cumulative marginal frequency distribution.
!	stdobsmarg is standard errors for observed marginal frequencies.
!	stdobsmargs2 is for standard errors for observed totals of cross products of item scores.
!	stdobsmargwtsum is the standard deviation of the observed
!		marginal frequency distribution.
!	stdobsmarg2 is standard errors for observed two-way marginal frequencies.
!	stdobspcummargwtsum is the asymptotic standard deviation of
!		the observed cumulative marginal probability distribution.
!	stdobspmarg is standard errors for observed marginal proportions.
!	stdobspmargs2 is for standard errors for observed averages of cross products of item scores.
!	stdobspmargwtsum is the asymptotic standard deviation of the
!		observed marginal probability distribution.
!	stdobspmarg2 is standard errors for observed two-way marginal proportions.
!	stdobsppred is the asymptotic standard deviation of the observed average
!		of products of category indicators and predictors.
!	stdobspred is the asymptotic standard deviation of the observed average of products
!		of category indicators and predictors.
!	stdobsptheta is the asymptotic standard deviation of the observed average of products
!        of category indicators and transformed latent variables.
!	stdobspwtitem is the asymptotic standard deviation of the observed average of products
!		of category indicators and weighted sums.
!	stdobspwtitem is the asymptotic standard deviation of the observed sum of products
!		of category indicators and weighted sums.
!	stdobstheta is the asymptotic standard deviation of the observed total of products
!        of category indicators and transformed latent variables.
!	stdresidcummargwtsum is the asymptotic standard deviation
!		of the residual cumulative marginal frequency
!		distribution.
!	stdresidirf is for standard error of residual for item response function.
!	stdresidmarg is standard error of residual for marginal frequencies.
!	stdresidmargs2 is standard error of residual of total of cross products of item scores.
!	stdresidmargwtsum is the asymptotic standard deviation of the
!		residual marginal frequency distribution.
!	stdresidmarg2 is standard error of residual for two-way marginal frequencies.
!	stdresidpcummargwtsum is the asymptotic standard deviation of
!		the residual cumulative marginal probability distribution.
!	stdresidpmarg is standard error of residual for marginal proportions.
!	stdresidpmargs2 is standard error of residual of average of cross products of item scores.
!	stdresidpmargwtsum is the asymptotic standard deviation of the
!		residual marginal probability distribution.
!	stdresidpmarg2 is standard error of residual for two-way marginal proportions.
!	stdresidppred is the asymptotic standard deviation of the residual for the average of products
!		of category indicators and predictors.
!	stdresidpred is the asymptotic standard deviation of the residual of the sum of products
!		of category indicators and predictors.
!	stdresidptheta is the asymptotic standard deviation of residual average of products
!        of category indicators and transformed latent variables.
!	stdresidpwtitem is the asymptotic standard deviation of the residual for the average of products
!		of category indicators and weighted sums.
!	stdresidtheta is the asymptotic standard deviation of the residual total of products
!        of category indicators and transformed latent variables.
!	stdresidwtitem is the asymptotic standard deviation of the residual of the sum of products
!		of category indicators and weighted sums.

!	theta is the array of adaptive quadrature points.
!	thetacheck contains points at which item response functions are to be computed.
!	transition is the transition matrix.
!	v and vv are used in computations of estimated asymptotic covariance matrices if linear constraints are used
!		not associated with least squares.
!	wtsum is a vector used to compute a linear combination of the observations.
!   zreswt is the ratio of elements of reswt and sdreswt.
!   zreswtdist is the ratio of elements of reswtdist and sdreswtdist.
real(kind=8),allocatable::alpha(:,:),altalpha(:,:),altbeta(:),altpostdensity(:,:),&
    altcholnhess(:,:,:),alttheta(:,:,:),aresid(:)
real(kind=8),allocatable::beta(:)
real(kind=8),allocatable::catscore(:),cholnhess(:,:,:),constmat(:,:),constquad(:,:),constvec(:)
real(kind=8),allocatable::covartheta(:,:),covarthetaskill(:,:),covartrans(:,:),covarwt(:,:)
real(kind=8),allocatable::coveap(:,:),coveapskill(:,:),coveaptrans(:,:),coveapwt(:,:)
real(kind=8),allocatable::covmeanobsscale(:,:)
real(kind=8),allocatable::covmeanscale(:,:),covobs(:,:),covobsave(:,:),covobssum(:,:),&
    covnewtheta(:,:,:),covres(:,:),covresave(:,:),covressum(:,:),covreswt(:,:,:),covreswtdist(:,:,:)
real(kind=8),allocatable::covs(:,:,:),covscale(:,:,:),covsum(:,:,:)
real(kind=8),allocatable::covtheta(:,:,:),covtrans(:,:,:)
real(kind=8),allocatable::cweight(:)
real(kind=8),allocatable::design(:,:),distfreq(:),distprob(:)
real(kind=8),allocatable::eacovbeta(:,:),eacovbeta_complex(:,:),eacovbeta_louis(:,:),eacovbeta_sandwich(:,:)
real(kind=8),allocatable::eacovgam(:,:),eacovgam_complex(:,:),eacovgam_louis(:,:),eacovgam_sandwich(:,:)
real(kind=8),allocatable::eacovgaminv(:,:),eacovgaminv_complex(:,:),eacovgaminv_louis(:,:)
real(kind=8),allocatable::extvar(:,:)
real(kind=8),allocatable::fitave(:),fitsum(:),fitcummargwtsum(:),fitirf(:,:)
real(kind=8),allocatable::fitmarg(:),fitmargwtsum(:),fitmarg2(:,:),fitmargs2(:,:)
real(kind=8),allocatable::fitpcummargwtsum(:),fitpmarg(:),fitpmargwtsum(:),fitpmarg2(:,:),fitpmargs2(:,:)
real(kind=8),allocatable::fitppred(:,:),fitpred(:,:),fitptheta(:,:),fitpwtitem(:)
real(kind=8),allocatable::fittheta(:,:),fitwt(:,:),fitwtdist(:,:),fitwtitem(:)
real(kind=8),allocatable::gamma(:),grad(:),graddes(:),grads(:,:),gradsdes(:,:)
real(kind=8),allocatable::guessres(:),guessresa(:)
real(kind=8),allocatable::indvar(:,:),informgam(:,:),informgam_complex(:,:),itemmatrix(:,:)
real(kind=8),allocatable::lintran(:,:)
real(kind=8),allocatable::meancovobsscale(:,:),meancovscale(:,:)
real(kind=8),allocatable::meancovtheta(:,:),meancovthetaskill(:,:),meancovtrans(:,:),meancovwt(:,:)
real(kind=8),allocatable::meaneap(:),meaneapskill(:),meaneaptrans(:),meaneapwt(:)
real(kind=8),allocatable::meannewtheta(:,:),meanscale(:,:),means(:,:),meansum(:,:),meantheta(:,:),meantrans(:,:)
real(kind=8),allocatable::nhess(:,:),nhessdes(:,:),nhessdeschol(:,:),nhessdesplus(:,:)
real(kind=8),allocatable::obsave(:),obscummargwtsum(:),observedscales(:,:),obsirf(:,:)
real(kind=8),allocatable::obsmarg(:),obsmargs2(:,:),obsmargwtsum(:),obsmarg2(:,:)
real(kind=8),allocatable::obspcummargwtsum(:),obspmarg(:),obspmargs2(:,:),obspmargwtsum(:),obspmarg2(:,:)
real(kind=8),allocatable::obsppred(:,:),obspred(:,:),obsptheta(:,:),&
    obspwtitem(:),obsscalerel(:),obsscaleres(:,:),obssum(:)
real(kind=8),allocatable::obstheta(:,:),obsweight(:),obsweightsel(:),obswt(:,:),obswtdist(:,:),obswtitem(:)
real(kind=8),allocatable::offset(:),offsettran(:)
real(kind=8),allocatable::plower(:),postdensity(:,:),predlist(:,:)
real(kind=8),allocatable::presented(:),presentedirf(:),presentedpred(:)
real(kind=8),allocatable::presenteds2(:,:),presentedtheta(:),presentedwtitem(:),presented2(:,:)
real(kind=8),allocatable::prob(:),probdist(:),pupper(:)
real(kind=8),allocatable::quadweight(:),quadpoint(:,:)
real(kind=8),allocatable::rdesign(:,:),releap(:),releapskill(:),reltrans(:),relwt(:),relmarg(:)
real(kind=8),allocatable::resave(:),ressum(:)
real(kind=8),allocatable::residacummargwtsum(:),residairf(:,:)
real(kind=8),allocatable::residamarg(:),residamargs2(:,:),residamargwtsum(:),residamarg2(:,:)
real(kind=8),allocatable::residapred(:,:),residatheta(:,:),residawtitem(:)
real(kind=8),allocatable::residcummargwtsum(:),residirf(:,:)
real(kind=8),allocatable::residmarg(:),residmargs2(:,:),residmargwtsum(:),residmarg2(:,:)
real(kind=8),allocatable::residpcummargwtsum(:),residpmarg(:)
real(kind=8),allocatable::residpmargs2(:,:),residpmargwtsum(:),residpmarg2(:,:)
real(kind=8),allocatable::residppred(:,:),residpred(:,:),residptheta(:,:),residpwtitem(:)
real(kind=8),allocatable::residtheta(:,:),residwtitem(:),reswt(:,:),reswtdist(:,:)
real(kind=8),allocatable::scale(:,:),scalecov(:,:),scalemean(:),scaleobsmean(:),scaleobscov(:,:)
real(kind=8),allocatable::scalerel(:),scalevar(:,:)
real(kind=8),allocatable::scorefreq(:),scores(:)
real(kind=8),allocatable::sdfreq(:),sdobs(:),sdobsave(:),sdobssum(:),sdprob(:)
real(kind=8),allocatable::sdres(:),sdresave(:),sdressum(:),sdreswt(:,:),sdreswtdist(:,:)
real(kind=8),allocatable::seobsmarg(:),setguess(:),setslope(:,:)
real(kind=8),allocatable::stdguessres(:)
real(kind=8),allocatable::stdobscummargwtsum(:)
real(kind=8),allocatable::stdobsmarg(:),stdobsmargs2(:,:),stdobsmargwtsum(:),stdobsmarg2(:,:)
real(kind=8),allocatable::stdobspcummargwtsum(:),stdobspmarg(:),stdobspmargs2(:,:),stdobspmargwtsum(:),stdobspmarg2(:,:)
real(kind=8),allocatable::stdobsppred(:,:),stdobspred(:,:),stdobsptheta(:,:),stdobspwtitem(:)
real(kind=8),allocatable::stdobstheta(:,:),stdobswtitem(:)
real(kind=8),allocatable::stdresidcummargwtsum(:),stdresidirf(:,:)
real(kind=8),allocatable::stdresidmarg(:),stdresidmargs2(:,:),stdresidmargwtsum(:),stdresidmarg2(:,:)
real(kind=8),allocatable::stdresidpcummargwtsum(:),stdresidpmarg(:),stdresidpmargs2(:,:),stdresidpmargwtsum(:),stdresidpmarg2(:,:)
real(kind=8),allocatable::stdresidppred(:,:),stdresidptheta(:,:),stdresidpred(:,:),stdresidpwtitem(:)
real(kind=8),allocatable::stdresidtheta(:,:),stdresidwtitem(:)
real(kind=8),allocatable::stobsmarg(:)
real(kind=8),allocatable::theta(:,:,:),thetacheck(:,:),transition(:,:)
real(kind=8),allocatable::v(:,:),vv(:,:)
real(kind=8),allocatable::wtsum(:,:),wtsumdist(:,:)
real(kind=8),allocatable::zreswt(:,:),zreswtdist(:,:)


!	Name the run.
call gettitle(title)
!	Get assignment of output.
call outputfiles()

!	Output the run.
call outputtitle(title,unittitle)
!	Find the basic structure of the input data.

call getdataspec(fileformat,filename,nexternal,nitems,nobs,nobsstart,npred,&
		complx,distract,recode,recodedist,stratify,useid,usepsu,weight)
npredtot=npred+nexternal
if(recode)then
	allocate(numcodes(nitems),stat=al)
	if(al/=0)stop "Recoding allocation failed."
	call getnumcodes(numcodes)
	if(sum(numcodes)>0)then
		allocate(recodetab(2,sum(numcodes)),stat=al)
		if(al/=0)stop "Recoding allocation failed."
		call getrecodes(numcodes,recodetab)
	end if
end if
if(recodedist)then
	allocate(numdistcodes(nitems),stat=al)
	if(al/=0)stop "Distractor recoding allocation failed."
	call getnumdistcodes(numdistcodes)
	if(sum(numdistcodes)>0)then
		allocate(recodedisttab(2,sum(numdistcodes)),stat=al)
		if(al/=0)stop "Distractor recoding allocation failed."
		call getdistrecodes(numdistcodes,recodedisttab)
	end if
end if

!	Obtain space for data.
allocate(dat(nitems,nobs),indvar(npred,nobs),extvar(nexternal,nobs),obsweight(nobs),&
	prob(nobs),stat=al)
if(al/=0) stop "Unable to allocate space for data storage"
if(distract)then
	allocate(distdat(nitems,nobs),stat=al)
	if(al/=0) stop "Unable to allocate space for storage of original item codes"
end if
if(stratify)then
	allocate(stratum(nobs),stat=al)
	if(al/=0) stop "Unable to allocate  space for storage of data on stratrum"
end if
if(useid)then
	allocate(id(nobs),stat=al)
	if(al/=0) stop "Unable to allocate  space for storage of examinee identifications"
end if
if(usepsu) then


	allocate(psu(nobs),stat=al)
	if(al/=0) stop "Unable to allocate space for storage of data on primary sampling units."
	
end if

!	Obtain the data

call getdata(fileformat,nexternal,nitems,nobs,npred,&
	numcodes,numdistcodes,recodetab,recodedisttab,distract,recode,recodedist,stratify,useid,usepsu,weight,&
	id,dat,distdat,psu,stratum,extvar,indvar,obsweight)

if(stratify.or.usepsu)call sortstratum(psu,stratum,stratify,usepsu,nstratum)

if(usepsu)then
	allocate(npsu(nstratum),stat=al)
	if(al/=0) stop "Unable to allocate space for storage of data on primary sampling units."
	call getpsucount(psu,stratum,npsu)
end if


!	Set the parameters for maximum-likelihood estimation.
call getparam(maxit,maxita,maxitb,nr,twostages,&
		changemin,kappa,maxdalpha,tau,tol,tola,tolc,tolres,tolsing)
!	Define the latent vector, quadrature weights, and quadrature points.
!	Select sample of starting observations by systematic sampling if starting iterations are used.
if(nobsstart<minobs)twostages=.false.

if(twostages)then
	allocate(obsweightsel(nobs),stat=al)
	if(al/=0) stop "Unable to allocate space for storage of weights for starting iterations."
	obsweightsel=0.0_8
	gap=(nobs+0.5_8)/nobsstart
	do obs=1,nobsstart
		row=obs*gap
		obsweightsel(row)=obsweight(row)
	end do

end if
!	Get number of factors, number of skills, and transition matrix.
call getdim(dimlatin,dimlatout,custom)

!	Allocate arrays dependent only on number of factors and number of skills.
allocate(factorname(dimlatin),skillname(dimlatout),dimpoints(dimlatin),&
	fixdiag(dimlatin),independentf(dimlatin),no_lin(dimlatin),rasch(dimlatout),raschslope1(dimlatout),&
	lintran(dimlatout,dimlatin),stat=al)
if(al/=0) stop "Arrays for factors and skills not allocated."
!	Allocate arrays dependent on number of factors and number of predictors.
allocate(predlinmap(dimlatin,npred),predquadmap(dimlatin,dimlatin,npred),stat=al)
if(al/=0) stop "Arrays for predictor maps to factors not allocated."
!	Get linear transformation of latent vector.
call getlintran(custom,lintran)

!	Quadrature data.

call getquadsize(dimpoints,ngrid,nquad,&
	cross,equalpoints,even,fullgrid,gausshermite,grid,normal,simplex,&
	pspread)



allocate(quadpoint(dimlatin,nquad),quadweight(nquad),coord(dimlatin,nquad),&
	cweight(nquad),stat=al)
if(al/=0) stop "Unable to allocate storage for quadrature."
if(grid) call coordinates(dimpoints,coord,cweight)

call getquad(coord,dimpoints,cross,equalpoints,&
		even,gausshermite,grid,simplex,cweight,pspread,&
		quadpoint,quadweight)

!	factor specifications.
call getfactorspecs(normal,factorname,fixdiag,fixquad,independentf,no_lin,noquad,predlinmap,predquadmap)



call getskillspecs(factorname,skillname,rasch,raschslope1)		
	
	


!	Obtain basic item  data.
allocate(choices(nitems),intdim(nitems),itemname(nitems),numcat(nitems),numcatobs(nitems),&
	skillnumbers(nitems),slopedim(dimlatout,nitems),&
	betweenitem(nitems),catmap(nitems),fixedguess(nitems),guess(nitems),nominal(nitems),&
	specialitemint(nitems),specialitemslope(nitems),&
	setguess(nitems),setslope(dimlatout,nitems),stat=al)
if(al/=0) stop "Unable to allocate storage for item information."

call getitemdata(itemname,choices,intdim,numcat,numcatobs,skillnumbers,slopedim,&
		betweenitem,catmap,constguess,fixedguess,guess,nominal,specialitemint,specialitemslope,setguess)

ncat=sum(numcat)
ncatobs=sum(numcatobs)
!	Beta dimension
nparam=ncat*(dimlatout+1)+npred*dimlatin*(dimlatin+3)/2
!
!	Rank of itemmatrix.
dimitemmatrix=sum(intdim)+sum(slopedim)+npred*dimlatin*(dimlatin+3)/2

!	Set up itemmatrix
allocate(itemmatrix(nparam,dimitemmatrix),betamap(5,nparam),stat=al)
if(al/=0) stop "Unable to allocate storage for item design information."
totalitems=itemspresented(dat,numcatobs,obsweight)
totalobs=obspresented(dat,numcatobs,obsweight)

!	Set betamap
betamap=mapbeta(dimlatin,dimlatout,npred,numcat)

if(totalitems<=0.0_8) stop "No items presented."
call outputdata(filename,nitems,nobs,npred,unitoutdata,totalitems,totalobs)
!	Set up category mapping.
allocate(catobsrange(2,ncatobs),stat=al)

call getitemdetail(intdim,numcat,numcatobs,slopedim,&
		catmap,guess,nominal,specialitemint,&
		specialitemslope,catobsrange,itemmatrix,setslope)
!	Get predictor names.

allocate(predname(npredtot),stat=al)
if(al/=0) stop "Unable to allocate storage for predictor names."

call getprednames(predname)


!	Get design specifications.



call designspecs(intdim,slopedim,constguess,&
		fixdiag,fixedguess,guess,&
		predlinmap,predquadmap,rasch,raschslope1,constdim,constdim1,dimdesign,proport,specialtrans)

allocate(constmat(dimdesign,constdim),constvec(constdim),&
	design(nparam,dimdesign),offset(nparam),paramname(dimdesign),&
	beta(nparam),gamma(dimdesign), stat=al)
if(al/=0) stop "Unable to allocate storage for model design and model parameters."

call getdesign(factorname,itemname,predname,skillname,choices,intdim,&
	numcat,slopedim,&
	constguess,fixdiag,fixedguess,guess,&
	predlinmap,predquadmap,rasch,raschslope1,specialtrans,&
	itemmatrix,setguess,setslope,&
	paramname,design,gamma,offset)

if(constdim>0)call getconstraint(constmat,constvec)

!	Read initial gamma and beta settings.

if(twostages) then
	call interceptstart(catobsrange,dat,numcat,numcatobs,design,obsweightsel,offset,tolsing,gamma)

else
	
	call interceptstart(catobsrange,dat,numcat,numcatobs,design,obsweight,offset,tolsing,gamma)
end if

call getgamma(gamma)

beta=offset+matmul(design,gamma)



!       Set alpha and cholnhess.
if(normal) then
	allocate(alpha(dimlatin,nobs),cholnhess(dimlatin,dimlatin,nobs),stat=al)
	if(al/=0) stop "Unable to allocate storage for adaptive quadrature."

	call getalpha(alpha,cholnhess)
end if
!	Reduce model.
nreduce=0
do row=1,nparam
	if(maxval(abs(design(row,:)))>0.0_8) nreduce=nreduce+1
	
end do
!	Parameter mapping.
allocate(rbetamap(5,nreduce),rdesign(nreduce,dimdesign),stat=al)
if(al/=0)stop "Parameter mappings not allocated successfully."
dimno=1
do row=1,nparam
	if(maxval(abs(design(row,:)))>0_8)then
		
		rbetamap(:,dimno)=betamap(:,row)
		rdesign(dimno,:)=design(row,:)
		dimno=dimno+1
	end if
end do

!	Iterate
allocate(constquad(dimdesign,dimdesign),grad(nreduce),graddes(dimdesign),grads(nreduce,nobs),&
	nhess(nreduce,nreduce),nhessdes(dimdesign,dimdesign),nhessdeschol(dimdesign,dimdesign),&
	nhessdesplus(dimdesign,dimdesign),postdensity(nquad,nobs),theta(dimlatin,nquad,nobs),&
	stat=al)
if(al/=0) stop "Unable to allocate storage for gradients and hessians."
!	Check on progress reporting.
	call progressreport(printprogstart,printprogstartstd,printprog,printprogstd)

!	First round.
if(twostages) then




	
	call ndmmaxlik("for starting iterations",catobsrange,constdim1,dat,&
		maxit,maxita,maxitb,numcat,numcatobs,rbetamap,unititerationstart,&
		normal,nr,printprogstart,printprogstartstd,proport,&
		changemin,constmat,constvec,design,indvar,kappa,lintran,&
		maxdalpha,obsweightsel,offset,quadpoint,quadweight,rdesign,&
		tau,tol,tola,tolc,tolsing,&
		alpha,beta,cholnhess,gamma,constquad,grad,graddes,grads,&
		loglik,nhess,nhessdes,nhessdeschol,nhessdesplus,&
		postdensity,prob,theta)
else
    printprogstart=.false.
end if
!	Second round.


call ndmmaxlik("for regular iterations",catobsrange,constdim1,dat,&
		maxit,maxita,maxitb,numcat,numcatobs,rbetamap,unititeration,&
		normal,nr,printprog,printprogstd,proport,&
		changemin,constmat,constvec,design,indvar,kappa,lintran,&
		maxdalpha,obsweight,offset,quadpoint,quadweight,rdesign,&
		tau,tol,tola,tolc,tolsing,&
		alpha,beta,cholnhess,gamma,constquad,grad,graddes,grads,&
		loglik,nhess,nhessdes,nhessdeschol,nhessdesplus,&
		postdensity,prob,theta)




!	Estimated asymptotic covariance matrices for gamma.
allocate(eacovgam(dimdesign,dimdesign),eacovgam_louis(dimdesign,dimdesign),&
	eacovgam_sandwich(dimdesign,dimdesign),&
	eacovgaminv(dimdesign,dimdesign),eacovgaminv_louis(dimdesign,dimdesign),&
	informgam(dimdesign,dimdesign),stat=al)
if(al/=0) stop "Unable to allocate storage for estimated asymptotic covariances"
!	get  output options



call setoutput()
if(maxita<=0.or.(.not.normal))then
    if(printml)then
        call closeunit(unitml)
        printml=.FALSE.
    end if
    if(printmp)then
        call closeunit(unitmp)
        printmp=.FALSE.
    end if
    if(printalpha)then
        call closeunit(unitalpha)
        printalpha=.FALSE.
    end if
end if

if(printprogstart)call closeunit(unititerationstart)

if(printprog)call closeunit(unititeration)
!	Covariance matrices.
informgam=matmul(transpose(rdesign),matmul(information(.true.,grads,obsweight),rdesign))
eacovgam_louis=chol(informgam,tolsing)
modeldim=dimdesign-constdim
if(constdim1<constdim)eacovgam_louis=chol(informgam+constquad,tolsing)
eacovgaminv_louis=invert(eacovgam_louis)
if(constdim1>0) then
	allocate(v(dimdesign,constdim1),vv(constdim1,constdim1),stat=al)
	if(al/=0) stop 'Allocation for estimated asymptotic covariances failed for constraint case.'
	v=matmul(eacovgaminv_louis,constmat(:,1:constdim1))
	vv=matmul(transpose(v),constmat(:,1:constdim1))
	vv=chol(vv,tolsing)
	vv=invert(vv)
	eacovgaminv_louis=eacovgaminv_louis-matmul(v,matmul(vv,transpose(v)))
end if
if(constdim1<constdim)then
	eacovgam_louis=matmul(eacovgaminv_louis,matmul(informgam,eacovgaminv_louis))
else
	
	eacovgam_louis=eacovgaminv_louis
end if
if(.not.nr) then
	eacovgam=eacovgam_louis
	eacovgam_sandwich=eacovgam
	eacovgaminv=eacovgaminv_louis
else
	eacovgaminv=invert(nhessdeschol)
	if(constdim1>0) then
		v=matmul(eacovgaminv,constmat(:,1:constdim1))
		vv=matmul(transpose(v),constmat(:,1:constdim1))
		vv=chol(vv,tolsing)
		vv=invert(vv)
		eacovgaminv=eacovgaminv-matmul(v,matmul(vv,transpose(v)))
	end if
	if(constdim1<constdim)then
		eacovgam=matmul(eacovgaminv,matmul(nhessdes,eacovgaminv))
	else
		eacovgam=eacovgaminv
	end if
	eacovgam_sandwich=matmul(eacovgaminv,matmul(informgam,eacovgaminv))
end if
if(printbeta.or.printbetacov.or.printbetacov_sandwich.or.printbetacov_complex) then
	allocate(betaname(nparam),stat=al)
    if(al/=0) stop "Unable to allocate storage for beta names."
	betaname=betanames(factorname,itemname,predname,skillname,numcat,beta)
end if
if(printbeta.or.printbetacov)then
    allocate(eacovbeta(nparam,nparam),stat=al)
    if(al/=0) stop "Unable to allocate storage for estimated asymptotic covariance matrix of beta."
	
    eacovbeta=matmul(design,matmul(eacovgam,transpose(design)))
end if
if(printbeta.or.printbetacov_louis)then
    allocate(eacovbeta_louis(nparam,nparam),stat=al)
    if(al/=0) stop "Unable to allocate storage for estimated Louis asymptotic covariance matrix of beta."
    eacovbeta_louis=matmul(design,matmul(eacovgam_louis,transpose(design)))
end if
if(printbeta.or.printbetacov_sandwich)then
    allocate(eacovbeta_sandwich(nparam,nparam),stat=al)
    if(al/=0) stop "Unable to allocate storage for estimated sandwich asymptotic covariance matrix of beta."
    eacovbeta_sandwich=matmul(design,matmul(eacovgam_sandwich,transpose(design)))
end if
    
	


if(complx)then
	allocate(eacovgam_complex(dimdesign,dimdesign),eacovgaminv_complex(dimdesign,dimdesign),&
		informgam_complex(dimdesign,dimdesign),stat=al)
	if(al/=0) stop "Unable to allocate storage for estimated asymptotic covariances for complex sampling."
	
	informgam_complex=matmul(transpose(rdesign),&
		matmul(information_complex(npsu,nstratum,psu,stratum,stratify,usepsu,grads,obsweight),rdesign))
		
	
		
		
		
	eacovgam_complex=matmul(matmul(eacovgaminv,informgam_complex),eacovgaminv)
    if(printbeta.or.printbetacov_complex)then
        allocate(eacovbeta_complex(nparam,nparam),stat=al)
        if(al/=0) stop "Unable to allocate storage for estimated complex asymptotic covariance matrix of beta."
        eacovbeta_complex=matmul(design,matmul(eacovgam_complex,transpose(design)))
    end if	
end if
!	Save probabilities.
if(printprob)then
	call outputprob(id,unitprob,useid,prob)
	call closeunit(unitprob)
end if
!   Distractor probabilities.
if(distract)then
    numchoices=sum(choices)
    allocate(distfreq(numchoices),distmap(numchoices),distprob(numchoices),scorefreq(ncatobs),stat=al)
    if(al/=0) stop "Unable to allocate storage for probabilities of distractor choices."
    call getfreq(choices,dat,distdat,numcatobs,obsweight,distmap,distfreq,distprob,scorefreq)
    if(printprobdist.or.printfreq.or.printentdist)then
        allocate(sdfreq(numchoices),sdprob(numchoices),stat=al)
        if(al/=0) stop 'Unable to allocate storage for standard errors for distractor frequencies and proportions.'
        call sdfreqs(choices,dat,distdat,distmap,npsu,nstratum,numcatobs,obsweight,&
            psu,stratum,complx,stratify,usepsu,distfreq,distprob,scorefreq,sdfreq,sdprob)
    end if
	if(printfreq) call printfreqs(itemname,choices,numcatobs,unitfreq,&
		distmap,distfreq,distprob,scorefreq,sdfreq,sdprob)
    if(printprobdist.or.printentdist)then

        allocate(probdist(nobs),stat=al)
        if(al/=0) stop "Unable to allocate storage for probabilities of distractor vectors."
        probdist=getprobdist(dat,distdat,choices,numcatobs,distprob,prob)
        if(printprobdist)then
			call outputprob(id,unitprobdist,useid,probdist)
			call closeunit(unitprobdist)
		end if
    end if
end if


!	Entropy output
if(printent) then

	ent=-loglik/totalitems
	sdentropy=sdent(loglik,obsweight,prob,totalitems)
	akent=akaike(modeldim,loglik,totalitems)
	ghent=gh(eacovgaminv,informgam,loglik,totalitems)
	if(complx) then
	
		ghent_complex=gh(eacovgaminv,informgam_complex,loglik,totalitems)
		
		sdentropy_complex=sdent_complex(npsu,nstratum,psu,stratum,stratify,usepsu,obsweight,prob,totalitems)
	end if
	
	call entoutput('Model prediction',modeldim,unitinfo,complx,akent,ent,ghent,ghent_complex,loglik,sdentropy,sdentropy_complex)
    call closeunit(unitinfo)
end if

if(printentdist) then
    loglikdist=sum(obsweight*log(probdist),probdist>0.0_8)
	entdist=-loglikdist/totalitems
	sdentropydist=sdent(loglikdist,obsweight,probdist,totalitems)
	modeldistdim=modeldim+count(distfreq>0.0_8)-count(scorefreq>0.0_8)
	akentdist=akaike(modeldistdim,loglikdist,totalitems)
	ghentdist=gh(eacovgaminv,informgam,loglikdist-modeldistdim+modeldim,totalitems)
	if(complx) then
	
		ghentdist_complex=gh(eacovgaminv,informgam_complex,loglikdist,totalitems)&
            +sum(sdprob*sdprob/distprob,distprob>0.0_8)/totalitems
		
		sdentropydist_complex=sdent_complex(npsu,nstratum,psu,stratum,stratify,usepsu,obsweight,probdist,totalitems)
	end if
	
	call entoutput('Distractor prediction',modeldistdim,unitinfodist,complx,&
        akentdist,entdist,ghentdist,ghentdist_complex,loglikdist,sdentropydist,&
        sdentropydist_complex)
    call closeunit(unitinfodist)
end if
!	gamma output
if(allocated(probdist))deallocate(probdist)
if(printparam)then
    call gammaoutput(paramname,unitparam,complx,eacovgam,eacovgam_complex,eacovgam_louis,eacovgam_sandwich,gamma)
    call closeunit(unitparam)
end if
if(printparamcov) then
    call matoutput("Model-based estimated asymptotic covariance matrix",&
        paramname,unitparamcov,eacovgam)
    call closeunit(unitparamcov)
end if
if(printparamcov_louis) then
    call matoutput("Louis estimate of asymptotic covariance matrix",&
        paramname,unitparamcov_louis,eacovgam_louis)
    call closeunit(unitparamcov_louis)
end if
if(printparamcov_sandwich) then
    call matoutput("Sandwich estimate of asymptotic covariance matrix",&
        paramname,unitparamcov_sandwich,eacovgam_sandwich)
    call closeunit(unitparamcov_sandwich)
end if
if(printparamcov_complex.and.complx) then
    call matoutput("Complex sampling estimate of asymptotic covariance matrix",&
        paramname,unitparamcov_complex,eacovgam_complex)
    call closeunit(unitparamcov_complex)
end if
    
!	beta output

if(printbeta)then
    call betaoutput(betaname,unitbeta,complx,beta,eacovbeta,eacovbeta_complex,eacovbeta_louis,eacovbeta_sandwich)
    call closeunit(unitbeta)
end if

if(printbetacov) then
    call matoutput("Model-based estimated asymptotic covariance matrix for beta",&
        betaname,unitbetacov,eacovbeta)
    call closeunit(unitbetacov)
end if
if(printbetacov_louis) then
    call matoutput("Louis estimate of asymptotic covariance matrix for beta",&
        betaname,unitbetacov_louis,eacovbeta_louis)
    call closeunit(unitbetacov_louis)
end if
if(printbetacov_sandwich) then
    call matoutput("Sandwich estimate of asymptotic covariance matrix for beta",&
        betaname,unitbetacov_sandwich,eacovbeta_sandwich)
    call closeunit(unitbetacov_sandwich)
end if
if(printbetacov_complex.and.complx) then
    call matoutput("Complex sampling estimate of asymptotic covariance matrix for beta",&
        betaname,unitbetacov_complex,eacovbeta_complex)
    call closeunit(unitbetacov_complex)
end if
!	Covariance matrices
if(printmeans)then
	call getmeans(npredlist,allmeans)
	if(allmeans.or.npredlist>0) then
		meanid=.true.
		if(allmeans)npredlist=nobs
		allocate(predlist(npred,npredlist),covs(dimlatin,dimlatin,npredlist),means(dimlatin,npredlist),stat=al)
		if(al/=0) stop 'Unable to allocate list of predictors.'
		if(.not.allmeans)then
			allocate(predlabels(npredlist),stat=al)
			if(al/=0) stop 'Unable to allocate list of predictors.'
			call getpredlist(predlabels,predlist)
		else
			if(useid) then
				allocate(predlabels(nobs),stat=al)
				if(al/=0) stop 'Unable to allocate list of predictors.'
				predlabels=id
			end if
			meanid=useid
			predlist=indvar
		end if
		if(normal)then
			
			call getmeancov(beta,predlist,tolsing,covs,means)
		else
			call getmeancovm(beta,predlist,quadpoint,quadweight,covs,means)
		end if
		
		call  eapoutput('of underlying latent vectors given predictors',factorname,predlabels,unitmeans,meanid,covs,means)
		deallocate(predlist,covs,means)
        call closeunit(unitmeans)
	end if
end if
!	eap output

if(printpost.or.printeap.or.printeapscale.or.printeapskill.or.printeaptrans.or.printml.or.printmp.or. &
	printrel.or.printrelskill.or.printreltrans.or.printeapwt.or.printrelwt.or. &
    printreswt.or.printreswtdist.or.printscalerel.or.printpwt.or. &
    printpwtdist.or.printobsscaleres.or.printobsscalerel)then
	allocate(eapmask(nitems),altbeta(size(beta)),stat=al)
	if(al/=0) stop "EAP arrays not allocated successfully."

	altbeta=beta
	call seteapoutput(dimtrans,dimwtsum,dimwtsumdist,eapmask,altbeta)
	altpost=.false.
	if(any(beta-altbeta/=0.0_8).or.(.not.all(eapmask))) altpost=.true.
	if(dimwtsum>0) then
		allocate(wtname(dimwtsum),wtsum(dimwtsum,ncat),&
			meansum(dimwtsum,nobs),covsum(dimwtsum,dimwtsum,nobs),&
			stat=al)
		if(al/=0) stop "Allocation failed for eap's for weighted sums."
		call getwtsum(skillname,catobsrange,numcatobs,slopedim,wtname,wtsum)
	end if
	if(altpost)then
		allocate(alttheta(dimlatin,nquad,nobs),altpostdensity(nquad,nobs),stat=al)
		if(al/=0) stop "EAP arrays not allocated successfully."
		if(normal)then
			allocate(altalpha(dimlatin,nobs),altcholnhess(dimlatin,dimlatin,nobs),stat=al)
			if(al/=0) stop "Arrays for maximum posterior likelihood not allocated successfully."
			altalpha=alpha
			altcholnhess=cholnhess
		end if
        if(printpost.or.printeap.or.printeapscale.or.printeapskill.or.printeaptrans.or. &
            printmp.or. &
            printrel.or.printrelskill.or.printreltrans.or.printeapwt.or.printrelwt.or. &
            printreswt.or.printreswtdist.or.printscalerel.or.printpwt.or. &
            printpwtdist.or.printobsscaleres.or.printobsscalerel)then

		    call posterior(catobsrange,dat,maxita,numcat,numcatobs,eapmask,normal,altbeta,changemin,indvar,lintran,&
			    maxdalpha,quadpoint,quadweight,&
			    tau,tola,tolsing,altalpha,altcholnhess,altpostdensity,alttheta)
        end if
	end if

end if

!   Posterior distributions.

if(printpost)then
	if(.not.altpost) then	
		call posterioroutput(factorname,id,unitpost,useid,postdensity,theta)
	else
		call posterioroutput(factorname,id,unitpost,useid,altpostdensity,alttheta)
	end if
    call closeunit(unitpost)
end if
!   EAP results.
if(printeap.or.printeapskill.or.printrel.or.printrelskill)then
	allocate(meantheta(dimlatin,nobs),covtheta(dimlatin,dimlatin,nobs),stat=al)
	if(al/=0) stop "EAP arrays not allocated successfully."	
	if(.not.altpost)then
		call eap(postdensity,theta,covtheta,meantheta)
	else
		call eap(altpostdensity,alttheta,covtheta,meantheta)
	end if
	if(printeap)then
        call eapoutput("for latent vector",factorname,id,uniteap,useid,covtheta,meantheta)
        call closeunit(uniteap)
    end if
	
	if(printrel.or.printrelskill)then
		allocate(releap(dimlatin),meaneap(dimlatin),coveap(dimlatin,dimlatin),&
			meancovtheta(dimlatin,dimlatin),covartheta(dimlatin,dimlatin),&
			stat=al)
		if(al/=0) stop "Reliability array not allocated successfully"

		call reltheta(covtheta,meantheta,obsweight,coveap,covartheta,meancovtheta,meaneap,releap)

		if(printrel)then
            call reloutput("for latent vector",factorname,unitrel,covartheta,&
                coveap,meancovtheta,meaneap,releap)
            call closeunit(unitrel)
        end if
		
		
	end if
	if(printeapskill.or.printrelskill)then

		allocate(meannewtheta(dimlatout,nobs),covnewtheta(dimlatout,dimlatout,nobs),stat=al)
		if(al/=0) stop "EAP arrays not allocated successfully."

		call eapskill(covtheta,lintran,meantheta,covnewtheta,meannewtheta)

		if(printeapskill)then
            call eapoutput("for transformed latent vector",skillname,id,uniteapskill,useid,covnewtheta,meannewtheta)
            call closeunit(uniteapskill)
        end if
		if(printrelskill)then
			allocate(releapskill(dimlatout),meaneapskill(dimlatout),coveapskill(dimlatout,dimlatout),&
				meancovthetaskill(dimlatout,dimlatout),covarthetaskill(dimlatout,dimlatout),stat=al)
			if(al/=0) stop "Reliability array not allocated successfully"

			call relthetaskill(covartheta,coveap,lintran,meancovtheta,meaneap,&
				covarthetaskill,coveapskill,meancovthetaskill,meaneapskill,releapskill)

			call reloutput("for transformed latent vector",skillname,unitrelskill,covarthetaskill,&
				coveapskill,meancovthetaskill,meaneapskill,releapskill)
            call closeunit(unitrelskill)
			deallocate(releapskill,meaneapskill,coveapskill,&
				meancovthetaskill,covarthetaskill)
		end if
		deallocate(meannewtheta,covnewtheta)
	end if
	deallocate(covtheta,meantheta)
	if(printrel)deallocate(covartheta,coveap,meancovtheta,meaneap,releap)
end if
!	EAP results for transformed latent vectors.
if(dimtrans>0.and.(printeaptrans.or.printreltrans.or.printml)) then
	
	
	allocate(transname(dimtrans),stat=al)
    call gettrans(dimlatin,transname)
	if(al/=0) stop "EAP arrays not allocated successfully."
    if(printeaptrans.or.printreltrans)then
        allocate(covtrans(dimtrans,dimtrans,nobs),meantrans(dimtrans,nobs),stat=al)

	    if(.not.altpost)then
		    call eaptrans(dimtrans,postdensity,theta,covtrans,meantrans)
	    else
		    call eaptrans(dimtrans,altpostdensity,alttheta,covtrans,meantrans)
	    end if
	    if(printeaptrans)then
	       call eapoutput("for transformation of latent vector",transname,id,uniteaptrans,useid,covtrans,meantrans)
	       call closeunit(uniteaptrans)
	    end if
	    if(printreltrans)then
		    allocate(covartrans(dimtrans,dimtrans),coveaptrans(dimtrans,dimtrans),&
			    meancovtrans(dimtrans,dimtrans),meaneaptrans(dimtrans),reltrans(dimtrans),stat=al)
		    if(al/=0) stop "Reliability array not allocated successfully"
		    call reltheta(covtrans,meantrans,obsweight,coveaptrans,covartrans,meancovtrans,meaneaptrans,reltrans)
		    call reloutput("for transformation of latent vector",transname,unitreltrans,covartrans,coveaptrans,&
			    meancovtrans,meaneaptrans,reltrans)
	        call closeunit(unitreltrans)
		    deallocate(coveaptrans,covartrans,meancovtrans,meaneaptrans,reltrans)
        end if
        deallocate(covtrans,meantrans)
	end if

	
end if


!   EAP results for weighted sums.
if((printeapwt.or.printrelwt).and.dimwtsum>0)then

	
	if(.not.altpost)then
		call eapwtsum(dat,numcat,numcatobs,eapmask,beta,lintran,postdensity,theta,wtsum,covsum,meansum)
	else
		call eapwtsum(dat,numcat,numcatobs,eapmask,beta,lintran,altpostdensity,alttheta,wtsum,covsum,meansum)
	end if
	if(printeapwt) then
		call eapoutput("for weighted sum",wtname,id,uniteapwt,useid,covsum,meansum)
		call closeunit(uniteapwt)
	end if
	if(printrelwt) then
		allocate(relwt(dimwtsum),meaneapwt(dimwtsum),coveapwt(dimwtsum,dimwtsum),&
			meancovwt(dimwtsum,dimwtsum),covarwt(dimwtsum,dimwtsum),stat=al)
		if(al/=0) stop "Reliability array not allocated successfully"
		call reltheta(covsum,meansum,obsweight,coveapwt,covarwt,meancovwt,meaneapwt,relwt)
!	releap contains reliabilities.
		call reloutput("for weighted sum",wtname,unitrelwt,covarwt,coveapwt,meancovwt,meaneapwt,relwt)
        call closeunit(unitrelwt)
		deallocate(covarwt,coveapwt,meancovwt,meaneapwt,relwt)
	end if
	deallocate(wtname,wtsum,meansum,covsum)
end if
!   Individual residuals for weighted sums.
if(printreswt.and.dimwtsum>0)then
    allocate(covreswt(dimwtsum,dimwtsum,nobs),fitwt(dimwtsum,nobs),obswt(dimwtsum,nobs),&
        reswt(dimwtsum,nobs),sdreswt(dimwtsum,nobs),zreswt(dimwtsum,nobs),&
        stat=al)
    if(al/=0) stop "Unable to allocate space for individual residuals for weighted sums."
	if(.not.altpost)then
		call  reswtsum(catobsrange,dat,numcat,numcatobs,beta,lintran,postdensity,theta,wtsum,&
			covreswt,fitwt,obswt,reswt,sdreswt,zreswt)
	else
		call reswtsum(catobsrange,dat,numcat,numcatobs,beta,lintran,altpostdensity,alttheta,wtsum,&
			covreswt,fitwt,obswt,reswt,sdreswt,zreswt)
	end if
    
    call resoutput("for individual weighted sums",wtname,id,unitreswt,useid,&
        covreswt,fitwt,obswt,reswt,sdreswt,zreswt)
    call closeunit(unitreswt)
    deallocate(covreswt,fitwt,obswt,reswt,sdreswt,zreswt)
end if
!	Individual distractor residuals.
if(distract.and.printreswtdist.and.dimwtsumdist>0)then
    allocate(wtnamedist(dimwtsumdist),wtsumdist(dimwtsumdist,numchoices),&
        covreswtdist(dimwtsumdist,dimwtsumdist,nobs),&
        fitwtdist(dimwtsumdist,nobs),obswtdist(dimwtsumdist,nobs),reswtdist(dimwtsumdist,nobs),&
		sdreswtdist(dimwtsumdist,nobs),zreswtdist(dimwtsumdist,nobs),&
        stat=al)
    if(al/=0) stop "Allocation failed for individual residuals for distractor sums."

    call getwtsumdist(skillname,choices,distmap,numcatobs,slopedim,wtnamedist,wtsumdist)
	if(.not.altpost)then
		call reswtsumdist(catobsrange,choices,dat,distdat,distmap,numcat,numcatobs,&
			beta,distprob,lintran,postdensity,theta,wtsumdist,&
			covreswtdist,fitwtdist,obswtdist,reswtdist,sdreswtdist,zreswtdist)
	else
		call reswtsumdist(catobsrange,choices,dat,distdat,distmap,numcat,numcatobs,&
			beta,distprob,lintran,altpostdensity,alttheta,wtsumdist,&
			covreswtdist,fitwtdist,obswtdist,reswtdist,sdreswtdist,zreswtdist)
	end if
		
	call resoutput("for individual weighted distractor sums",wtnamedist,id,&
        unitreswtdist,useid,covreswtdist,fitwtdist,obswtdist,reswtdist,sdreswtdist,zreswtdist)
    call closeunit(unitreswtdist)
    deallocate(covreswtdist,fitwtdist,obswtdist,reswtdist,sdreswtdist,zreswtdist)
end if
if(printmp)then
	if(.not.altpost) then

		call mpoutput(factorname,id,unitmp,useid,alpha,cholnhess)
	else

		call mpoutput(factorname,id,unitmp,useid,altalpha,altcholnhess)
	end if
    call closeunit(unitmp)
end if

if(printgrad.or.printguess.or.printirfres.or.printmarginres.or.printmargin2res.or.printmargins2res.or.&
	printpreditemres.or.printmarginwtsumres.or.&
    printthetaitemres.or.printwtitemres.or.printobsscaleres)then
	allocate(gradsdes(dimdesign,nobs),stat=al)
	if(al/=0)stop "Allocation of gradients for residual analysis was not successful."
	do obs=1,nobs
		gradsdes(:,obs)=matmul(transpose(rdesign),grads(:,obs))
	end do
	if(printgrad)then
		call outputgrad(id,paramname,unitgrad,useid,gradsdes)
		call closeunit(unitgrad)
	end if
else
	allocate(gradsdes(0,0))
end if	

if(printalpha)then
	call savealpha(unitalpha,alpha,cholnhess)
	call closeunit(unitalpha)
end if
!   Guessing tests.
if(minval(numcat)==2.and.printguess)then
    allocate(guessres(nitems),guessresa(nitems),stdguessres(nitems),stat=al)
    if(al/=0) stop 'Allocation failed for guessing tests'
    
    call guesstest(dat,npsu,nstratum,numcat,psu,stratum,&
        complx,stratify,usepsu,beta,eacovgaminv_louis,&
        gradsdes,lintran,obsweight,postdensity,theta,guessres,guessresa,stdguessres)

    call printguesstest(itemname,numcat,unitguess,&
        guessres,guessresa,stdguessres)
    call closeunit(unitguess)
    deallocate(guessres,guessresa,stdguessres)
end if
!	Marginals for items.
if(printirf)then

	call getthetasize(quadpoint,irfcount,custompoints)

	if(irfcount>0) then
	allocate(thetacheck(size(quadpoint,1),irfcount),stat=al)
	
	if(al/=0) stop 'Allocation failed for points for item response functions'
	if(custompoints)then
		call getthetas(thetacheck)
	else
		thetacheck=quadpoint
	end if

	allocate(fitirf(ncatobs,size(thetacheck,2)),obsirf(ncatobs,size(thetacheck,2)),presentedirf(nitems),stat=al)
	if(al/=0) stop 'Allocation failed for item response functions'
	if(printirfres)then
		allocate(residairf(ncatobs,size(thetacheck,2)),residirf(ncatobs,size(thetacheck,2)),&
			stdresidirf(ncatobs,size(thetacheck,2)),stat=al)
			if(al/=0) stop 'Allocation failed for residuals for item response functions'
	end if

	call irf(catobsrange,dat,npsu,nstratum,numcat,numcatobs,&
		psu,stratum,complx,printirfres,stratify,usepsu,&
		beta,eacovgaminv_louis,gradsdes,indvar,lintran,&
		obsweight,postdensity,prob,theta,thetacheck,tolres,&
		fitirf,obsirf,presentedirf,&
		residairf,residirf,&
		stdresidirf)

	call printirfs(factorname,itemname,numcatobs,unitirf,printirfres,&
		fitirf,obsirf,presentedirf,&
		residairf,residirf,stdresidirf,thetacheck)

    call closeunit(unitirf)

	deallocate(fitirf,obsirf,presentedirf)
	if(printirfres)then
        deallocate(residairf,residirf,stdresidirf)
        call closeunit(unitirf)
    end if
	end if
end if
if(printmargin)then

	allocate(fitmarg(ncatobs),fitpmarg(ncatobs),obsmarg(ncatobs),&
		obspmarg(ncatobs),presented(nitems),relmarg(ncatobs),&
		seobsmarg(ncatobs),stdobsmarg(ncatobs),stdobspmarg(ncatobs),&
		stobsmarg(ncatobs),stat=al)
	if(al/=0)stop "Allocation failed for marginal distributions."

	if(printmarginres)then
		allocate(residamarg(ncatobs),residmarg(ncatobs),&
			residpmarg(ncatobs),stdresidmarg(ncatobs),stdresidpmarg(ncatobs),stat=al)
		if(al/=0)stop "Allocation failed for marginal distributions."
	end if	
	
	call marginaldist(catobsrange,dat,npsu,nstratum,numcat,numcatobs,&
		psu,stratum,complx,printmarginres,stratify,usepsu,&
		beta,eacovgaminv_louis,gradsdes,lintran,&
		obsweight,postdensity,theta,tolres,&
		fitmarg,fitpmarg,obsmarg,obspmarg,presented,&
		relmarg,residamarg,residmarg,residpmarg,&
		seobsmarg,stdobsmarg,stdobspmarg,&
		stdresidmarg,stdresidpmarg,stobsmarg)

	call printmarginaldist(itemname,numcatobs,unitmargin,printmarginres,&
		fitmarg,fitpmarg,obsmarg,obspmarg,presented,&
		relmarg,residamarg,residmarg,residpmarg,&
		seobsmarg,stdobsmarg,stdobspmarg,stdresidmarg,stdresidpmarg,stobsmarg)

    call closeunit(unitmargin)


	deallocate(fitmarg,fitpmarg,obsmarg,&
		obspmarg,presented,relmarg,&
		seobsmarg,stdobsmarg,stdobspmarg,&
		stobsmarg)
	if(printmarginres) then
	
        deallocate(residamarg,residmarg,&
            residpmarg,stdresidmarg,stdresidpmarg)
        call closeunit(unitmargin)
	
    end if
end if

!	Marginal distributions for weighted sums.
if(printmarginwtsum.or.printwtitem.or.printeapscale.or.printscalerel&
    .or.printpwt.or.printobservedscale.or.printobsscale.or.printobsscalerel)then
	numweights=getnumweights(slopedim)
	
	if(numweights>0)then

		allocate(weightnames(numweights),maxw(nitems),minw(nitems),weights(ncatobs,numweights),stat=al)
		if(al/=0)stop "Weights for distributions not allocated successfully."


		call getweights(skillname,numcatobs,slopedim,weightnames,weights)
		if(printpwt) then

            allocate(plower(nobs),pupper(nobs),stat=al)
            if(al/=0)stop "Upper and lower probabilities not allocated successfully."
        end if
		if(printeapscale.or.printscalerel.or.printobservedscale.or.printobsscale.or.printobsscalerel)then
			
			call getscalecount(numberscales)
		end if
		
		do row=1,numweights
			call setwtlimits(numcatobs,weights(:,row),maxw,minw)
			maxscore=sum(maxw)
			minscore=sum(minw)

			if(printmarginwtsum)then
				allocate(fitcummargwtsum(minscore:maxscore),fitmargwtsum(minscore:maxscore),&
					fitpcummargwtsum(minscore:maxscore),fitpmargwtsum(minscore:maxscore),&
					obscummargwtsum(minscore:maxscore),obsmargwtsum(minscore:maxscore),&
					obspcummargwtsum(minscore:maxscore),obspmargwtsum(minscore:maxscore),&
					stdobscummargwtsum(minscore:maxscore),&
					stdobsmargwtsum(minscore:maxscore),&
					stdobspcummargwtsum(minscore:maxscore),stdobspmargwtsum(minscore:maxscore),&
					stat=al)
				if(al/=0)stop "Allocation failure for distribution of weighted sums."
			end if
			if(printmarginwtsumres)then
				allocate(residacummargwtsum(minscore:maxscore),&
					residamargwtsum(minscore:maxscore),&
					residcummargwtsum(minscore:maxscore),residmargwtsum(minscore:maxscore),&
					residpcummargwtsum(minscore:maxscore),&
					residpmargwtsum(minscore:maxscore),&
					stdresidcummargwtsum(minscore:maxscore),stdresidmargwtsum(minscore:maxscore),&
					stdresidpcummargwtsum(minscore:maxscore),stdresidpmargwtsum(minscore:maxscore),&
					stat=al)
				if(al/=0)stop "Allocation failure for distribution of weighted sums."
			end if
			
			if(printmarginwtsum)then
				call margindistwtsum(catobsrange,dat,maxscore,maxw,minscore,minw,&
					npsu,nstratum,&
					numcat,numcatobs,&
					psu,stratum,complx,printmarginwtsumres,stratify,usepsu,&
					beta,eacovgaminv_louis,gradsdes,&
					lintran,&
					obsweight,postdensity,theta,tolres,weights(:,row),&
					fitcummargwtsum,fitmargwtsum,fitpcummargwtsum,fitpmargwtsum,&
					obscummargwtsum,obsmargwtsum,obspcummargwtsum,obspmargwtsum,&
					presentedwtsum,&
					residacummargwtsum,residamargwtsum,&
					residcummargwtsum,residmargwtsum,residpcummargwtsum,&
					residpmargwtsum,stdobscummargwtsum,stdobsmargwtsum,&
					stdobspcummargwtsum,stdobspmargwtsum,&
					stdresidcummargwtsum,stdresidmargwtsum,&
					stdresidpcummargwtsum,stdresidpmargwtsum)
			
				call printmargindistwtsum(weightnames(row),maxscore,minscore,&
					unitmarginwtsum,printmarginwtsumres,&
					fitcummargwtsum,fitmargwtsum,fitpcummargwtsum,fitpmargwtsum,&
					obscummargwtsum,obsmargwtsum,obspcummargwtsum,obspmargwtsum,&
					presentedwtsum,&
					residacummargwtsum,residamargwtsum,&
					residcummargwtsum,residmargwtsum,residpcummargwtsum,&
					residpmargwtsum,stdobscummargwtsum,stdobsmargwtsum,&
					stdobspcummargwtsum,stdobspmargwtsum,&
					stdresidcummargwtsum,stdresidmargwtsum,&
					stdresidpcummargwtsum,stdresidpmargwtsum)
				deallocate(fitcummargwtsum,fitmargwtsum,&
					fitpcummargwtsum,fitpmargwtsum,&
					obscummargwtsum,obsmargwtsum,&
					obspcummargwtsum,obspmargwtsum,&
					stdobscummargwtsum,&
					stdobsmargwtsum,&
					stdobspcummargwtsum,stdobspmargwtsum)
				if(printmarginwtsumres)deallocate(residacummargwtsum,&
					residamargwtsum,&
					residcummargwtsum,residmargwtsum,&
					residpcummargwtsum,&
					residpmargwtsum,&
					stdresidcummargwtsum,stdresidmargwtsum,&
					stdresidpcummargwtsum,stdresidpmargwtsum)
			end if
            if(printpwt)then
               if(.not.altpost)then
                    call getpwt(catobsrange,dat,maxscore,maxw,minscore,minw,&
                        numcat,numcatobs,weights(:,row),&
                        beta,&
                        lintran,&
                        postdensity,theta,&
                        plower,pupper)

               else

                      call getpwt(catobsrange,dat,maxscore,maxw,minscore,minw,&
                         numcat,numcatobs,weights(:,row),&
                         altbeta,lintran,&
                         altpostdensity,alttheta,&
                         plower,pupper)

                end if
                call outputpwt(id,weightnames(row),unitpwt,useid,plower,pupper)

            end if
			if(printwtitem.and.row==1)then
				allocate(fitpwtitem(ncatobs),fitwtitem(ncatobs),obspwtitem(ncatobs),obswtitem(ncatobs),&
					presentedwtitem(nitems),stdobspwtitem(ncatobs),stdobswtitem(ncatobs),stat=al)
				if(al/=0)stop "Allocation failure for totals of products of item indicators and weighted sums."
			end if
			if(printwtitemres.and.row==1)then
				allocate(residawtitem(ncatobs),residpwtitem(ncatobs),residwtitem(ncatobs),&
					stdresidpwtitem(ncatobs),stdresidwtitem(ncatobs),stat=al)
				if(al/=0)stop "Allocation failure for residuals for totals of products of item indicators and weighted sums."
			end if
			if(printwtitem)then
				call marginwtitem(catobsrange,dat,maxw,minw,npsu,nstratum,numcat,numcatobs,&
					psu,stratum,weights(:,row),complx,printwtitemres,stratify,usepsu,&
					beta,eacovgaminv_louis,gradsdes,lintran,&
					obsweight,postdensity,theta,tolres,&
					fitpwtitem,fitwtitem,obspwtitem,obswtitem,presentedwtitem,&
					residawtitem,residpwtitem,residwtitem,stdobspwtitem,stdobswtitem,stdresidpwtitem,stdresidwtitem)
				call printmarginwtitem(itemname,weightnames(row),numcatobs,unitwtitem,&
					printwtitemres,fitpwtitem,fitwtitem,obspwtitem,obswtitem,presentedwtitem,&
					residawtitem,residpwtitem,residwtitem,&
					stdobspwtitem,stdobswtitem,stdresidpwtitem,stdresidwtitem)
				
			end if

			if((printeapscale.or.printscalerel.or.printobservedscale.or.printobsscale &
                .or.printobsscalerel).and.numberscales>0)then

				allocate(scalename(numberscales),scale(numberscales,minscore:maxscore),stat=al)
                if(al/=0)stop 'Unable to allocate space for scale scores.'
                call getscales(maxscore,minscore,numberscales,scalename,scale)
                if(printeapscale.or.printscalerel.or.printobsscalerel.or.printobsscaleres)then
                    allocate(meanscale(numberscales,nobs),covscale(numberscales,numberscales,nobs),stat=al)
                    if(al/=0)stop 'Unable to allocate space for scale scores.'
                end if
				if(printobservedscale.or.printobsscale.or.printobsscalerel)then
                    allocate(observedscales(numberscales,nobs),presence(nobs),stat=al)
                    if(al/=0)stop 'Unable to allocate space for scale scores.'
                end if
				if(printscalerel) then
					allocate(covmeanscale(numberscales,numberscales),&
						meancovscale(numberscales,numberscales),&
						scalemean(numberscales),scalecov(numberscales,numberscales),&
                        scalerel(numberscales), stat=al)
					if(al/=0)stop 'Unable to allocate space for scale scores.'
					
				end if

                if(printeapscale.or.printscalerel.or.printobsscalerel.or.printobsscaleres)then
                    if(.not.altpost)then
                        call eapscale(catobsrange,maxscore,maxw,minscore,&
                            minw,numberscales,numcat,numcatobs,weights(:,row),beta,&
                            lintran,postdensity,scale,theta,covscale,meanscale)
                    else
                        call eapscale(catobsrange,maxscore,maxw,minscore,&
                            minw,numberscales,numcat,numcatobs,weights(:,row),altbeta,&
                            lintran,altpostdensity,scale,alttheta,covscale,meanscale)
                    end if
                end if


                if(printobservedscale.or.printobsscale.or.printobsscalerel)&
                    call obsscale(dat,maxscore,maxw,minscore,&
                        minw,numberscales,numcatobs,weights(:,row),presence,&
                        scale,observedscales)

                if(printobservedscale)call obsoutput('Observed scale scores',scalename,&
                    id,unitobservedscale,useid,presence,observedscales)
!
                if(printobsscale.or.printobsscalerel)then

                    allocate(obsave(numberscales),obssum(numberscales),stat=al)
                    if(al/=0)stop 'Unable to allocate average and total scale scores'
                end if

                if(printobsscale)then
                    allocate(covobsave(numberscales,numberscales),covobssum(numberscales,numberscales),&
                        sdobsave(numberscales),sdobssum(numberscales),stat=al)
                    if(al/=0) stop 'Unable to allocate space for statistics on scale scores'
                end if
                if(printobsscalerel.or.printobsscaleres)then

                    allocate(fitave(numberscales),fitsum(numberscales),stat=al)
                    if(al/=0) stop 'Unable to allocate space for fitted averages.'
                end if
                if(printobsscalerel)then

                    allocate(covobs(numberscales,numberscales),&
                        covres(numberscales,numberscales),&
                        obsscalerel(numberscales),sdobs(numberscales),&
                        sdres(numberscales),stat=al)

                    if(al/=0) stop 'Unable to allocate space for reliability of observed scale scores'
                end if
                if(printobsscaleres)then
                    allocate(aresid(numberscales),covresave(numberscales,numberscales),&
                        covressum(numberscales,numberscales),&
                        resave(numberscales),ressum(numberscales),sdresave(numberscales),&
                        sdressum(numberscales),stat=al)
                    if(al/=0) stop 'Unable to allocate space for residuals for scale scores'
                end if
                if(printobsscale.or.printobsscalerel)then

                    call stdresid(npsu,nstratum,&
                        psu,stratum,complx,printobsscalerel,presence,&
                        printobsscaleres,printobsscale,stratify,usepsu,&
                        eacovgaminv_louis,meanscale,gradsdes,&
                        observedscales,obsweight,tolres,&
                        aresid,covobs,covobsave,covobssum,covres,covresave,covressum,&
                        fitave,fitsum,obsave,obssum,obsscalerel,resave,&
                        ressum,sdobs,sdobsave,sdobssum,sdres,sdresave,sdressum,totsum)

                end if
                if(printobsscalerel)then
                    call obsreloutput(' for observed scale scores',scalename,unitobsscalerel,&
                        covobs,covres,obsave,obsscalerel)
                    deallocate(covobs,covres,obsscalerel,sdobs)
                end if
                if(printobsscale)&
                    call obsscaleoutput(scalename,unitobsscale,printobsscaleres,&
                        aresid,covobsave,covobssum,covresave,covressum,&
                        fitave,fitsum,obsave,obssum,resave,&
                        ressum,sdobsave,sdobssum,sdresave,sdressum,totsum)
                if(printobsscaleres)deallocate(aresid,covresave,fitave,fitsum,&
                    covressum,obsave,obssum,resave,ressum,sdresave,&
                    sdressum)
                if(printobsscale)deallocate(covobsave,covobssum,&
                    sdobsave,sdobssum)



				if(printeapscale) call eapoutput("for scaled score",scalename,id,uniteapscale,&
                    useid,covscale,meanscale)
				if(printscalerel)then
					call reltheta(covscale,meanscale,obsweight,covmeanscale,&
                        scalecov,meancovscale,scalemean,scalerel)
                    call reloutput("for scaled score ",scalename,unitscalerel,scalecov,covmeanscale,&
                        meancovscale,scalemean,scalerel)
					deallocate(covmeanscale,meancovscale,scalecov,scalemean,scalerel)
				end if
				deallocate(scalename,scale)
                if(printeapscale.or.printscalerel) deallocate(covscale,meanscale)

			end if
		end do
		if(printwtitem)then
            deallocate(fitpwtitem,fitwtitem,obspwtitem,obswtitem,&
                presentedwtitem,stdobspwtitem,stdobswtitem)
            call closeunit(unitwtitem)
        end if
		if(printwtitemres)then
            deallocate(residawtitem,residpwtitem,residwtitem,&
					stdresidpwtitem,stdresidwtitem)
            call closeunit(unitwtitem)
        end if
        if(printmarginwtsum)call closeunit(unitmarginwtsum)
        if(printmarginwtsumres)call closeunit(unitmarginwtsum)

        if(printpwt)then
            deallocate(plower,pupper)
            call closeunit(unitpwt)
        end if
        if(printscalerel)call closeunit(unitscalerel)
        if(printeapscale)call closeunit(uniteapscale)
        if(printobsscale)call closeunit(unitobsscale)
        if(printobservedscale)call closeunit(unitobservedscale)
        if(printobsscalerel)call closeunit(unitobsscalerel)
        deallocate(maxw,minw,weightnames,weights)
	end if
	
end if

if(distract.and.printpwtdist)then
	allocate(plower(nobs),pupper(nobs),stat=al)
    if(al/=0)stop "Upper and lower probabilities not allocated successfully."
    
    numweightsdist=getnumweightsdist(slopedim)
    if(numweightsdist>0)then

        allocate(weightnamesdist(numweightsdist),maxw(nitems),minw(nitems),weightsdist(numchoices,numweightsdist),stat=al)
        if(al/=0)stop "Weights for distractor distributions not allocated successfully."

        call getweightsdist(skillname,choices,distmap,numcatobs,slopedim,weightnamesdist,weightsdist)

        
        do row=1,numweightsdist
            call setwtlimits(choices,weightsdist(:,row),maxw,minw)
            maxscore=sum(maxw)
            minscore=sum(minw)


            if(.not.altpost)then


                call getpwtdist(catobsrange,choices,dat,distdat,distmap,&
                    maxscore,maxw,minscore,minw,&
                    numcat,numcatobs,weightsdist(:,row),&
                    beta,distprob,&
                    lintran,&
                    postdensity,theta,&
                    plower,pupper)


             else

                call getpwtdist(catobsrange,choices,dat,distdat,distmap,&
                     maxscore,maxw,minscore,minw,&
                     numcat,numcatobs,weightsdist(:,row),&
                     altbeta,distprob,&
                     lintran,&
                     altpostdensity,alttheta,&
                     plower,pupper)
            end if


            call outputpwt(id,weightnamesdist(row),unitpwtdist,useid,plower,pupper)
        end do
    end if
    call closeunit(unitpwtdist)
    deallocate(maxw,minw,plower,pupper,weightnamesdist,weightsdist)
end if



!	Marginals for item pairs.
if(printmargin2)then
	allocate(fitmarg2(ncatobs,ncatobs),fitpmarg2(ncatobs,ncatobs),obsmarg2(ncatobs,ncatobs),&
		obspmarg2(ncatobs,ncatobs),presented2(nitems,nitems),&
		stdobsmarg2(ncatobs,ncatobs),stdobspmarg2(ncatobs,ncatobs),&
		stat=al)
	if(al/=0)stop "Allocation failed for two-way marginal distributions."
	if(printmargin2res)then
		allocate(residamarg2(ncatobs,ncatobs),residmarg2(ncatobs,ncatobs),&
			residpmarg2(ncatobs,ncatobs),stdresidmarg2(ncatobs,ncatobs),&
			stdresidpmarg2(ncatobs,ncatobs),stat=al)
		if(al/=0)stop "Allocation failed for two-way marginal distributions."
	end if	
	
	call margindist2(catobsrange,dat,npsu,nstratum,numcat,numcatobs,&
		psu,stratum,complx,printmargin2res,stratify,usepsu,&
		beta,eacovgaminv_louis,gradsdes,lintran,&
		obsweight,postdensity,theta,tolres,&
		fitmarg2,fitpmarg2,obsmarg2,obspmarg2,presented2,&
		residamarg2,residmarg2,residpmarg2,&
		stdobsmarg2,stdobspmarg2,&
		stdresidmarg2,stdresidpmarg2)
	
	call printmarginal2(itemname,numcatobs,unitmargin2,printmargin2res,&
		fitmarg2,fitpmarg2,obsmarg2,obspmarg2,presented2,&
		residamarg2,residmarg2,residpmarg2,&
		stdobsmarg2,stdobspmarg2,stdresidmarg2,stdresidpmarg2)
	call closeunit(unitmargin2)
	deallocate(fitmarg2,fitpmarg2,obsmarg2,&
		obspmarg2,presented2,&
		stdobsmarg2,stdobspmarg2)
	if(printmargin2res) then
        deallocate(residamarg2,residmarg2,&
            residpmarg2,stdresidmarg2,stdresidpmarg2)
        call closeunit(unitmargin2)
    end if
end if
if(printmargins2)then
	allocate(fitmargs2(nitems,nitems),fitpmargs2(nitems,nitems),obsmargs2(nitems,nitems),&
		obspmargs2(nitems,nitems),presenteds2(nitems,nitems),scores(ncatobs),&
		stdobsmargs2(nitems,nitems),stdobspmargs2(nitems,nitems),&
		stat=al)
	if(al/=0)stop "Allocation failed for cross products of item scores."
	if(printmargins2res)then
		allocate(residamargs2(nitems,nitems),residmargs2(nitems,nitems),&
			residpmargs2(nitems,nitems),stdresidmargs2(nitems,nitems),&
			stdresidpmargs2(nitems,nitems),stat=al)
		if(al/=0)stop "Allocation failed for cross products of item scores."
	end if

	scores=getscores(numcatobs)

	call margindists2(catobsrange,dat,npsu,nstratum,numcat,numcatobs,&
		psu,stratum,complx,printmargins2res,stratify,usepsu,&
		beta,eacovgaminv_louis,gradsdes,lintran,&
		obsweight,postdensity,scores,theta,tolres,&
		fitmargs2,fitpmargs2,obsmargs2,obspmargs2,presenteds2,&
		residamargs2,residmargs2,residpmargs2,&
		stdobsmargs2,stdobspmargs2,&
		stdresidmargs2,stdresidpmargs2)
	
	call printmarginals2(itemname,unitmargins2,printmargins2res,&
		fitmargs2,fitpmargs2,obsmargs2,obspmargs2,presenteds2,&
		residamargs2,residmargs2,residpmargs2,&
		stdobsmargs2,stdobspmargs2,stdresidmargs2,stdresidpmargs2)
	
	deallocate(fitmargs2,fitpmargs2,obsmargs2,&
		obspmargs2,presenteds2,&
		stdobsmargs2,stdobspmargs2)
    call closeunit(unitmargins2)
	if(printmargin2res) then
        deallocate(residamargs2,residmargs2,&
            residpmargs2,stdresidmargs2,stdresidpmargs2)
        call closeunit(unitmargins2)
    end if

end if
if(printpreditem.and.npredtot>1)then
	allocate(fitppred(ncatobs,npredtot),fitpred(ncatobs,npredtot),&
		obsppred(ncatobs,npredtot),obspred(ncatobs,npredtot),presentedpred(nitems),&
		stdobsppred(ncatobs,npredtot),stdobspred(ncatobs,npredtot),stat=al)
	if(al/=0)stop "Allocation failed for sums and averages of products of category indicators and predictors."
	if(printpreditemres)then
		allocate(residapred(ncatobs,npredtot),&
			residppred(ncatobs,npredtot),residpred(ncatobs,npredtot),&
			stdresidppred(ncatobs,npredtot),&
			stdresidpred(ncatobs,npredtot),stat=al)
		if(al/=0)stop "Allocation failed for residuals for sums and averages of products of category indicators and predictors."
	end if

	call marginpred(catobsrange,dat,npsu,nstratum,numcat,numcatobs,&
		psu,stratum,complx,printpreditemres,stratify,usepsu,&
		beta,eacovgaminv_louis,extvar,indvar,gradsdes,lintran,&
		obsweight,postdensity,theta,tolres,&
		fitppred,fitpred,obsppred,obspred,presentedpred,&
		residapred,residppred,residpred,stdobsppred,stdobspred,stdresidppred,stdresidpred)
	
	call printmarginpred(itemname,predname,numcatobs,unitpreditem,&
		printpreditemres,fitppred,fitpred,obsppred,obspred,presentedpred,&
		residapred,residppred,residpred,&
		stdobsppred,stdobspred,stdresidppred,stdresidpred)

	deallocate(fitppred,fitpred,obsppred,obspred,presentedpred,stdobsppred,stdobspred)
	call closeunit(unitpreditem)
	
	if(printpreditemres) then
        deallocate(residapred,residppred,&
            residpred,stdresidppred,stdresidpred)
        call closeunit(unitpreditem)
    end if
	
end if

if(printthetaitem)then
    allocate(fitptheta(ncatobs,dimlatout),fittheta(ncatobs,dimlatout),&
        obsptheta(ncatobs,dimlatout),obstheta(ncatobs,dimlatout),presentedtheta(nitems),&
        stdobsptheta(ncatobs,dimlatout),stdobstheta(ncatobs,dimlatout),stat=al)
    if(al/=0)stop "Allocation failed for sums and averages of products of category indicators and latent variable."
    if(printthetaitemres)then
        allocate(residatheta(ncatobs,dimlatout),&
            residptheta(ncatobs,dimlatout),residtheta(ncatobs,dimlatout),&
            stdresidptheta(ncatobs,dimlatout),&
            stdresidtheta(ncatobs,dimlatout),stat=al)
        if(al/=0)stop "Allocation failed for sums and averages of products of category indicators and latent variables."
    end if

    call margintheta(catobsrange,dat,npsu,nstratum,numcat,numcatobs,&
        psu,stratum,complx,printthetaitemres,stratify,usepsu,&
        beta,eacovgaminv_louis,gradsdes,lintran,&
        obsweight,postdensity,theta,tolres,&
        fitptheta,fittheta,obsptheta,obstheta,presentedtheta,&
        residatheta,residptheta,residtheta,stdobsptheta,stdobstheta,stdresidptheta,stdresidtheta)

    call printmargintheta(itemname,skillname,numcatobs,unitthetaitem,&
        printthetaitemres,fitptheta,fittheta,obsptheta,obstheta,presentedtheta,&
        residatheta,residptheta,residtheta,&
        stdobsptheta,stdobstheta,stdresidptheta,stdresidtheta)
    call closeunit(unitthetaitem)
    deallocate(fitptheta,fittheta,obsptheta,&
        obstheta,presentedtheta,&
        stdobsptheta,stdobstheta)

        if(printthetaitemres) then
            deallocate(residatheta,residptheta,&
                residtheta,stdresidptheta,stdresidtheta)
        call closeunit(unitthetaitem)
    end if
        
end if
if(printml)then
    if(altpost)then
        call maxl(factorname,id,transname,catobsrange,dat,dimtrans,maxita,maxitb,&
            npred,numcat,numcatobs,unitml,useid,eapmask,altbeta,changemin,&
            lintran,maxdalpha,tau,tol,tolsing)
    else
        
        call maxl(factorname,id,transname,catobsrange,dat,dimtrans,maxita,maxitb,&
            npred,numcat,numcatobs,unitml,useid,eapmask,beta,changemin,&
            lintran,maxdalpha,tau,tol,tolsing)
    end if
    call closeunit(unitml)
end if
end program irtnorm
