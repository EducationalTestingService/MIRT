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


integer,intent(in)::nstratum
integer,intent(in),optional::npsu(:),psu(:),stratum(:)
logical,intent(in)::stratify,usepsu
real(kind=8),intent(in)::obsweight(:),prob(:),totalitems


!	al indicates allocation error.
!   nobs is the number of observations.
!	obs is an observation counter.


integer::al,nobs,obs
!   Variance goes to covobs.
!	Log probabilities.
real(kind=8)::covobs(1,1)
real(kind=8),allocatable::obsmat(:,:)
nobs=size(prob)
allocate(obsmat(1,nobs),stat=al)
if(al/=0) stop "Unable to allocate arrays for asymptotic variances computations for complex sampling."
obsmat=0.0_8
do obs=1,nobs
    if(prob(obs)>0.0_8)obsmat(1,obs)=log(prob(obs))
end do
covobs=information_complex(npsu,nstratum,psu,stratum,stratify,usepsu,&
    obsmat(1:1,:),obsweight)
sdent_complex=sqrt(covobs(1,1))/totalitems
deallocate(obsmat)

return
end function sdent_complex


