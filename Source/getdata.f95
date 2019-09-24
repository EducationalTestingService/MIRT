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
!	recode indicates if recodes are present.
!	recodedisttab is the table of distractor recodes for each item.
!	distract is .true. if distractors are read.
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
	numcodes,numdistcodes,recodetab,recodedisttab,distract,recode,recodedist,stratify,useid,usepsu,weight,&
	id,dat,distdat,psu,stratum,extvar,indvar,obsweight)
implicit none
character(len=*),intent(in)::fileformat
integer,intent(in)::nexternal,nitems,nobs,npred
integer,intent(in),optional::numcodes(:),numdistcodes(:),recodetab(:,:),recodedisttab(:,:)
logical,intent(in)::distract,recode,recodedist,stratify,useid,usepsu,weight

character(len=32),intent(out),optional::id(:)
integer,intent(out)::dat(:,:),distdat(:,:)
integer,intent(out),optional::psu(:),stratum(:)
real(kind=8),intent(out)::extvar(:,:),indvar(:,:),obsweight(:)
!	ident is for the id.
character(len=32),allocatable::ident(:)
!	code and code1 count codes.
!	io is flag for correct input.
!	item is an item number.
!	ndist gives the size of the distractor respdist.
!	nident gives the size of the id array ident.
!	nps gives the size of the psu array ps.
!	nstr gives the size of the stratum array str.
!	nwgt gives the size of the weight array wgt.
!	obs is counter for observations.
!	resp is the response array.

integer::code,code1,io,item,nident,ndist,nps,nstr,nwgt,obs,resp(nitems)
!	ps is for the psu.
!	respdist is the distractor response array.
!	str is for the stratum variable.
integer,allocatable::ps(:),respdist(:),str(:)
!	ext is the external array.
!	pr is the predictor array.
real(kind=8)::pr(npred-1),ext(size(extvar,1))
!	wgt is for the weight variable.
real(kind=8),allocatable::wgt(:)

!	Allocate needed arrays.
ndist=0
nident=0
nps=0
nstr=0
nwgt=0

if(useid)nident=1
if(usepsu)nps=1
if(stratify)nstr=1
if(weight)nwgt=1
if(distract)ndist=nitems

allocate(ident(nident),ps(nps),respdist(ndist),str(nstr),wgt(nwgt),stat=io)
if(io/=0)stop "Allocation of arrays for input failed."

obsweight=1.0_8
do obs=1,nobs
	if(fileformat/="*")then
		read(8,fmt=fileformat,iostat=io) ident,resp,respdist,pr,wgt,ext,str,ps
	else
		read(8,*,iostat=io) ident,resp,respdist,pr,wgt,ext,str,ps
	end if
	if(io>0) stop "Read error encountered"
	if(weight)then
		if(wgt(1)<0.0_8)  stop "Read error encountered"
	end if
	
!	Recode if necessary.
	if(recode)then
		code1=0
		do item=1,size(resp)
			if(numcodes(item)>0)then
				do code=1,numcodes(item)
					if(resp(item)==recodetab(1,code1+code))then
						resp(item)=recodetab(2,code1+code)
						exit
					end if
				end do
				code1=code1+numcodes(item)
			end if
		end do
	end if
	if(recodedist)then
		code1=0
		do item=1,size(resp)
			if(numdistcodes(item)>0)then
				do code=1,numdistcodes(item)
					if(respdist(item)==recodedisttab(1,code1+code))then
						respdist(item)=recodedisttab(2,code1+code)
						exit
					end if
				end do
				code1=code1+numdistcodes(item)
			end if
		end do
	end if
	dat(:,obs)=resp
	if(distract)distdat(:,obs)=respdist
	indvar(1,obs)=1.0_8
	indvar(2:npred,obs)=pr
	extvar(:,obs)=ext
	if(weight)obsweight(obs)=wgt(1)
	if(stratify)stratum(obs)=str(1)
	if(usepsu) psu(obs)=ps(1)
	if(useid) id(obs)=ident(1)
	
end do
close(8)
return
end subroutine getdata
