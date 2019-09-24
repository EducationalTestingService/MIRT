!	Design specifications.
!	Input:

!	intdim provides the number of columns in itemmatrix for intercepts associated with each item. 

!	slopedim provides the number of columns in itemmatrix for slopes associated with an item.
!	constguess is for a constant guessing parameter 
!	fixdiag is for fixed diagonals for the constant
!	fixedguess is for fixed guessing parameters.
!	fixquad is for constant quadratic parameters
!	guess is for guessing parameters.
!	predlinmap is map of predictors to factors for linear term.
!	predquadmap is map of predictors to factors for quadratic term.	
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
	predlinmap,predquadmap,rasch,raschslope1,constdim,constdim1,dimdesign,proport,&
	specialtrans)
implicit none
integer,intent(in)::intdim(:),slopedim(:,:)

integer,intent(out)::constdim,constdim1,dimdesign
logical,intent(in)::constguess,fixdiag(:),fixedguess(:),guess(:),&
		predlinmap(:,:),predquadmap(:,:,:),rasch(:),raschslope1(:)
logical,intent(out)::proport,specialtrans
	
!	Design specifications.
!	dimlatcov counts covariance parameters.
!	dimlatin is the dimension of the latent vector.
!	dimno is the count for independent true.
!	intdesdim is the dimension of the intercepts.
!	io is the flag for input error.
!	lindesdim is the dimension of the linear components of the density
!	quaddesdim is the dimension of the quadratic components of the density
!	slopedesdim is the dimension of the slopes.
integer::dimcount,dimlatcov,dimlatin,dimno,dimnoold,&
	intdesdim,io,lindesdim,quaddesdim,slopedesdim


namelist/designparameters/constdim,constdim1,dimdesign,proport,specialtrans
dimlatin=size(fixdiag)
!	Defaults.
constdim=0
constdim1=0
if(constguess)then
	intdesdim=sum(intdim)-count(guess)+1
else
	intdesdim=sum(intdim)-count(fixedguess)
end if
slopedesdim=sum(sum(slopedim,2),mask=(.not.rasch))+count(rasch.and.(.not.raschslope1))
lindesdim=count(predlinmap)
quaddesdim=count(predquadmap)-count(fixdiag)

dimdesign=intdesdim+slopedesdim+lindesdim+quaddesdim

proport=.true.
specialtrans=.false.

read(*,nml=designparameters,iostat=io)
if(io/=0)stop "Design parameters not read successfully."
if(constdim<0) constdim=0
if(constdim1<0) constdim1=0
if(constdim1>constdim) constdim1=constdim
if(dimdesign<1) stop "At least one independent parameter is needed."







return
end subroutine designspecs


