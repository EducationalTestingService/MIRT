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
interface


!	Estimate the cross information matrix.
!	adjust indicates adjustment for means.
!	grads is the array of observation gradients.
!	grads1 is another array of observation gradients.
!	obsweight contains observation weights.

    function crossinformation(adjust,grads,grads1,obsweight)
        implicit none
        logical,intent(in)::adjust
        real(kind=8),intent(in)::grads(:,:),grads1(:,:),obsweight(:)
        real(kind=8)::crossinformation(size(grads,1),size(grads1,1))
    end function crossinformation


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

    function information_complex(npsu,nstratum,psu,stratum,stratify,&
        usepsu,grads,obsweight)
        implicit none
        integer,intent(in)::nstratum
        integer,intent(in),optional::npsu(:),psu(:),stratum(:)
        logical,intent(in)::stratify,usepsu
        real(kind=8),intent(in)::grads(:,:),obsweight(:)
        real(kind=8)::information_complex(size(grads,1),size(grads,1))
    end function information_complex

!	probvec is used to compute underlying item probabilities.
!   numcat indicates the number of underlying categories per item.
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




end interface
integer,intent(in)::dat(:,:),npsu(:),nstratum,numcat(:),&
    psu(:),stratum(:)
logical,intent(in)::complx,stratify,usepsu
real(kind=8),intent(in)::beta(:),eacovgaminv_louis(:,:),&
    gradsdes(:,:),lintran(:,:),obsweight(:),postdensity(:,:),theta(:,:,:)
real(kind=8),intent(out)::guessres(:),guessresa(:),stdguessres(:)
!	al is allocation indicator.
!	counter counts 
integer::al,counter,dimdesign,dimlatout,item,ncat,&
    nitems,nlin,nobs,nscale,obs,quad,resp(size(guessres))
logical::datamask(size(guessres))
real(kind=8)::avediff,diff,covobs(1,1),covsum(size(gradsdes,1)),newtheta(size(lintran,1)),&
    locations(sum(numcat)),&
    probcat(sum(numcat)),respvec(size(guessres)),scales(size(lintran,1),sum(numcat)),&
    slopes(size(gradsdes,1)),totalobs
real(kind=8),allocatable::obsmat(:,:,:)
dimdesign=size(gradsdes,1)
dimlatout=size(lintran,1)
guessres=0.0_8
guessresa=0.0_8
ncat=sum(numcat)
locations=beta(1:ncat)
nitems=size(numcat)
nlin=ncat*(dimlatout+1)+1
nobs=size(obsweight)
allocate(obsmat(1,nobs,nitems),stat=al)
if(al/=0)stop 'Allocation failed for guessing test'
nscale=ncat+1
scales=reshape(beta(nscale:nlin-1),(/dimlatout,ncat/))
stdguessres=0.0_8
totalobs=sum(obsweight)
do obs=1,nobs
    if(obsweight(obs)<=0.0_8)cycle
    resp=dat(:,obs)
    respvec=0.0_8
    datamask=.false.
    do item=1,nitems

        
        if(resp(item)>=0.and.resp(item)<2.and.numcat(item)==2)datamask(item)=.TRUE.





        
    end do
    if(count(datamask)==0)cycle
!	Cycle through quadrature points.
    do quad=1,size(theta,2)
!		Probability computation.
        newtheta=matmul(lintran,theta(:,quad,obs))
        probcat=probvec(numcat,datamask,locations,scales,newtheta)
        counter=2
        do item=1,size(numcat)
            if(datamask(item))respvec(item)=&
                respvec(item)+postdensity(quad,obs)*(resp(item)-probcat(counter))/probcat(counter)
            counter=counter+numcat(item)
        end do
    end do
    obsmat(1,obs,:)=respvec
!   Residual
   guessres=guessres+obsweight(obs)*respvec

end do
!   Standard error
do item=1,nitems
    if(numcat(item)==2)then
        if(.not.complx)then
            covsum=reshape(crossinformation(.true.,obsmat(:,:,item),gradsdes,obsweight),&
                (/dimdesign/))
            slopes=matmul(eacovgaminv_louis,covsum)
            avediff=0.0_8
            do obs=1,nobs
                if(obsweight(obs)<=0.0_8)cycle
                obsmat(1,obs,item)=obsmat(1,obs,item)-dot_product(gradsdes(:,obs),slopes)
                avediff=avediff+obsmat(1,obs,item)*obsweight(obs)
            end do
            avediff=avediff/totalobs
            do obs=1,nobs
                diff=obsmat(1,obs,item)-avediff
                stdguessres(item)=stdguessres(item)+diff*diff*obsweight(obs)
            end do
        else
        
            covsum=reshape(crossinformation(.true.,obsmat(:,:,item),gradsdes,obsweight),(/dimdesign/))
            slopes=matmul(eacovgaminv_louis,covsum)
            do obs=1,nobs
                obsmat(1,obs,item)=obsmat(1,obs,item)-dot_product(gradsdes(:,obs),slopes)
            end do
            covobs=information_complex(npsu,nstratum,psu,stratum,stratify,usepsu,&
                obsmat(:,:,item),obsweight)
            stdguessres(item)=covobs(1,1)
        end if
        stdguessres(item)=sqrt(stdguessres(item))
        if(stdguessres(item)>0.0_8)guessresa(item)=guessres(item)/stdguessres(item)
    end if
    
end do
deallocate(obsmat)
return
end subroutine guesstest