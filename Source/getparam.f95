!	Some parameter settings.
!	maxit is the number of iterations for the main cycle.
!	maxita is the number of iterations for finding alpha.
!	maxitb is the maximum number of iterations used to find the step size.
!	nr is for Newton-Raphson.  The alternative uses an approximation to Newton-Raphson based on the Louis
!		approximation.
!	twostages indicates whether two rounds of iterations are used.
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
	
!	io is flag for read error.
integer::io
namelist/paramspec/maxit,maxita,maxitb,nr,twostages,&
	changemin,kappa,maxdalpha,tau,tolres,tolsing,tol,tola
	
!	Default values.
maxit=50
maxita=10
maxitb=10
nr=.true.
twostages=.true.
changemin=0.0625_8
kappa=2.0_8
maxdalpha=3.0_8
tau=0.1_8
tol=1.0d-5
tola=1.0d-4
tolc=1.0d-5
tolres=0.01_8
tolsing=1.0d-6

!	Read values.
read(*,nml=paramspec,iostat=io)
if(io/=0)stop "Parameter specifications not satisfactory"
if(maxit<0)maxit=0
if(maxita<0)maxita=0
if(maxitb<0)maxitb=0
if(changemin<0.0_8.or.changemin>=0.5_8)changemin=0.0625_8
if(kappa<=0.0_8) kappa=2.0_8

if(maxdalpha<=0.0_8)maxdalpha=3.0_8
if(tau<0.0_8.or.tau>=1.0_8)tau=0.1_8
if(tol<=0.0_8) tol=1.0d-5
if(tola<=0.0_8)tola=1.0d-4
if(tolc<=0.0_8) tolc=1.0d-5
if(tolres<=0.0_8) tolres=0.01_8
if(tolsing<=0.0_8) tolsing=1.0d-9
return
end subroutine getparam
