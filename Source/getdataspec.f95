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

!	idcheck is a character variable used in the input testing.
character(len=32),allocatable::idcheck(:)

!	At least minobs observations are needed.
integer,parameter::minobs=2
!	counter counts records.
!	io indicates a reading issue.
!	nfloat is the number of floating-point numbers to be read.
!	nid is the number of character variables to read.
!	nitems1 is nitems if distractors are not analyzed and 2 times nitems of distractors are analyzed.
!	nsamp is the number of variables related to sampling.
integer::counter,io,nfloat,nid,nitems1,nsamp
!	csamp is the specification for a complex sample.
!	resp is a response.
integer,allocatable::csamp(:),resp(:)
!	x is for floating point portions of a record.
real(kind=8),allocatable::x(:)
namelist/dataspec/fileformat,filename,nexternal,nitems,nobs,nobsstart,npred,&
	complx,distract,recode,recodedist,stratify,useid,usepsu,weight


!	Default values
fileformat="*"
filename="data.txt"
nexternal=0
nitems=0
nobs=0
nobsstart=1000
npred=1

complx=.false.
distract=.false.
recode=.false.
recodedist=.false.
stratify=.false.
useid=.false.
usepsu=.false.
weight=.false.

nid=0
nsamp=0


read(*,nml=dataspec,iostat=io)
if(io/=0)stop "Data specification information not valid"
if(nexternal<0) nexternal=0
if(nitems<1) stop "Number of items not positive"
if(npred<1) npred=1
if(.not.distract)recodedist=.false.
complx=stratify.or.usepsu.or.(weight.and.complx)
open(8,file=filename,iostat=io)
if(io/=0) stop "Data file cannot be opened."
if(stratify)nsamp=nsamp+1
if(usepsu)nsamp=nsamp+1
if(weight)then
	nfloat=npred+nexternal
else
	nfloat=npred+nexternal-1
end if
if(useid)nid=1
!	Count up number of observations.
counter=0
if(nobs<=0) nobs=huge(nobs)
nitems1=nitems
if(distract)nitems1=nitems1+nitems1
allocate(idcheck(nid),resp(nitems1),x(nfloat),csamp(nsamp))
do while(counter<nobs)
	
	if(fileformat/="*")then
		read(8,fmt=fileformat,iostat=io) idcheck,resp,x,csamp
	else
		read(8,*,iostat=io) idcheck,resp,x,csamp
	end if
	if(io>0) stop "Read error encountered"
	if(io<0) exit
	counter=counter+1
end do

nobs=counter
if(nobs<minobs) stop "Insufficient data for analysis."
!	Initial sample cannot be larger than final sample.
if(nobsstart<=0) nobsstart=1000
nobsstart=min(nobs,nobsstart)
rewind(8)
return
end subroutine getdataspec
