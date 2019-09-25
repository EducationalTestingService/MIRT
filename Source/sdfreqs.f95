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
interface
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
end interface


integer,intent(in)::choices(:),dat(:,:),distdat(:,:),distmap(:),npsu(:),nstratum,numcatobs(:),&
	psu(:),stratum(:)
logical,intent(in)::complx,stratify,usepsu
real(kind=8),intent(in)::distfreq(:),distprob(:),obsweight(:),scorefreq(:)
real(kind=8),intent(out)::sdfreq(:),sdprob(:)
!   a1 is an allocation indicator.
!   col is a distractor for an item.
!   counter marks distractors.

!	item is an item.


!	nitems is the number of items.

!	nobs is number of observations.
!   numchoices is total number of distractors.

!	obs counts observations.
!   position is a position of an item score
!   positiondist is a position of a distractor score.
!	resp is a response.
!   resdidst is a distractor response.
!   row is an item marker

integer::al,col,counter,item,nitems,nobs,numchoices,&
	obs,position,positiondist,resp,respdist,row


!	covobs is the covariance matrix for adjusted observed values.



!	obsmat is the array of response indicators and fits.
!		The first row is the observed indicator.
!		The second is 1 if consistent with the item score and 0 otherwise.


real(kind=8)::covobs(1,1)
real(kind=8),allocatable::obsmat(:,:)
!

!	Set up parameter arrays for simplified processing.



nitems=size(dat,1)
nobs=size(obsweight)
numchoices=sum(choices)
sdfreq=0.0_8
sdprob=0.0_8
!   Complex sampling requires an examimation of individual cases.
if(complx)then
	allocate(obsmat(2,nobs),stat=al)
	if(al/=0) stop 'Allocation of array failed for standard errors of distractor frequencies and proportions.'
    counter=1
    positiondist=1
    row=1
    do item=1,nitems
        do col=1,choices(item)
            if(distmap(counter)>=0) then
                obsmat=0.0_8
                do obs=1,nobs

                    if(obsweight(obs)<=0.0_8) cycle
                    resp=dat(item,obs)
                    respdist=distdat(item,obs)
                    if(resp<0.or.resp>=numcatobs(item))cycle

!	Verify if item score was observed.
                    position=distmap(positiondist+respdist)

                    if(resp==position-row)then
                        obsmat(2,obs)=1.0_8
                        if(respdist==col-1)obsmat(1,obs)=1.0_8

                    end if
                end do
                row=distmap(counter)
                covobs=information_complex(npsu,nstratum,psu,stratum,stratify,usepsu,&
					obsmat(1:1,:),obsweight)
                sdfreq(counter)=sqrt(covobs(1,1))

                obsmat(1,:)=obsmat(1,:)-distprob(counter)*obsmat(2,:)
                covobs=information_complex(npsu,nstratum,psu,stratum,stratify,usepsu,&
                    obsmat(1:1,:),obsweight)
                sdprob(counter)=sqrt(covobs(1,1)/scorefreq(row))
            end if
            counter=counter+1
        end do
        row=row+numcatobs(item)
        positiondist=positiondist+choices(item)
    end do
else


    do counter=1,numchoices
        row=distmap(counter)
        if(row<0) cycle
		sdfreq(counter)=distfreq(counter)*(1.0_8-distprob(counter))
        sdprob(counter)=sqrt(sdfreq(counter)/scorefreq(row))
        sdfreq(counter)=sqrt(sdfreq(counter))
    end do
end if

if(complx)deallocate(obsmat)
return
end subroutine sdfreqs