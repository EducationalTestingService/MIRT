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

subroutine reswtsumdist(catobsrange,choices,dat,distdat,distmap,numcat,numcatobs,beta,&
	distprob,lintran,postdensity,theta,wtsumdist,&
    covreswtdist,fitwtdist,obswtdist,reswtdist,sdreswtdist,zreswtdist)
implicit none
interface




!	Multinomial covariance matrix.
!	density is vector of probabilities.
!	meanth is vector of means.
!	quadpoint is array of weights.
    function mcovth(density,meanth,quadpoint)
        implicit none
        real(kind=8),intent(in)::density(:),meanth(:),quadpoint(:,:)
        real(kind=8)::mcovth(size(meanth),size(meanth))
    end function mcovth
!	Multinomial mean.
!	density is vector of probabilities
!	quadpoint is array of weights.
    function mmeanth(density,quadpoint)
        implicit none
        real(kind=8),intent(in)::density(:),quadpoint(:,:)
        real(kind=8)::mmeanth(size(quadpoint,1))
    end function mmeanth
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
!	probvecobsall is used to compute observed item probabilities given the latent vector.
!	catobsrange is the table of ranges of underlying categories per observed category.
!	numcat provides the number of underlying categories per item, and numcatobs provides the number
!		of observed categories per item.
!	mask indicates which items were presented.
!	probcat is the vector of underlying item probabilities.
	function probvecobsall(catobsrange,numcatobs,mask,probcat)
		implicit none
		integer,intent(in)::catobsrange(:,:),numcatobs(:)
		logical,intent(in)::mask(:)
		real(kind=8),intent(in)::probcat(:)
		real(kind=8)::probvecobsall(sum(numcatobs))
	end function probvecobsall

!	probvecobsdist is used to compute observed item distractor probabilities given the latent vector.
!	choices is the number of choices per item.
!	numcatobs provides the number of observed categories per item.
!	mask indicates which items were presented.
!	distprob provides distractor probabilities.
!	probcatobs is the vector of item probabilities.
    function probvecobsdist(choices,distmap,numcatobs,mask,distprob,probcatobs)
        implicit none

        integer,intent(in)::choices(:),distmap(:),numcatobs(:)
        logical,intent(in)::mask(:)
        real(kind=8),intent(in)::distprob(:),probcatobs(:)
        real(kind=8)::probvecobsdist(sum(choices))
    end function probvecobsdist
end interface
integer,intent(in)::catobsrange(:,:),choices(:),dat(:,:),distdat(:,:),distmap(:),numcat(:),numcatobs(:)


real(kind=8),intent(in)::beta(:),distprob(:),lintran(:,:),postdensity(:,:),theta(:,:,:),wtsumdist(:,:)
real(kind=8),intent(out)::covreswtdist(:,:,:),fitwtdist(:,:),obswtdist(:,:),reswtdist(:,:),sdreswtdist(:,:),zreswtdist(:,:)
!	item is an item.
!	ncat is the total number of underlying categories.
!	obs is observation number.
!   quad counts quadrature points.
!	resp is a single response.
!	respdist is a single distractor response.
!	row is a row counter.

integer::item,ncat,obs,quad,resp(size(numcatobs)),respdist(size(choices)),row
!   datamask is for identification of item responses for observations.
logical::datamask(size(choices))

!   covar is the covariance matrix of the differences of observed and duplicate sums.
!   covquad is the covariance matrix of sums given the latent vector.
!   covquaditem is the covariance matrix of weighted sums for items given the latent vector.

!	locations is for location parameters.

!	meanquad is for means for quadrature points.

!	newtheta is lintran(theta).
!	probcat is used for underlying marginal conditional probabilities.
!   probcatobs is used for marginal conditional probabilities of observations.
!	probcatobsdist is used for marginal conditional probabilities of original observations.
!	scales is for scale parameters.

real(kind=8)::covar(size(wtsumdist,1),size(wtsumdist,1))
real(kind=8)::covquad(size(wtsumdist,1),size(wtsumdist,1),size(postdensity,1)),&
    covquaditem(size(wtsumdist,1),size(wtsumdist,1),size(numcatobs),size(postdensity,1))
real(kind=8)::locations(sum(numcat))
real(kind=8)::meanquad(size(wtsumdist,1),size(postdensity,1)),meanquaditem(size(wtsumdist,1),size(numcatobs),size(postdensity,1))
real(kind=8)::newtheta(size(lintran,1))

real(kind=8)::probcat(sum(numcat)),probcatobs(sum(numcatobs)),probcatobsdist(sum(choices))
real(kind=8)::scales(size(lintran,1),sum(numcat))





!	Set up parameter arrays for simplified processing.
ncat=sum(numcat)
!	Item scale parameters.
locations=beta(1:ncat)
scales=reshape(beta(ncat+1:ncat*(1+size(lintran,1))),(/size(lintran,1),ncat/))

covar=0.0_8
covquad=0.0_8
covquaditem=0.0_8

covreswtdist=0.0_8
fitwtdist=0.0_8
obswtdist=0.0_8

reswtdist=0.0_8
sdreswtdist=0.0_8
zreswtdist=0.0_8




!	Obtain results for each observation.
obsloop:do obs=1,size(postdensity,2)

    resp=dat(:,obs)
	respdist=distdat(:,obs)

    datamask=.false.
    do item=1,size(numcat)
        if(resp(item)>=0.and.resp(item)<numcatobs(item))datamask(item)=.true.
    end do

!	Cycle through quadrature points.
    covar=0.0_8
    covquad=0.0_8
    covquaditem=0.0_8

    meanquad=0.0_8
    meanquaditem=0.0_8



    do quad=1,size(postdensity,1)
        newtheta=matmul(lintran,theta(:,quad,obs))


! Probability computation.


        probcat=probvec(numcat,datamask,locations,scales,newtheta)
        probcatobs=probvecobsall(catobsrange,numcatobs,datamask,probcat)

		

		probcatobsdist=probvecobsdist(choices,distmap,numcatobs,datamask,distprob,probcatobs)
		
		
        if(size(wtsumdist,1)>0)then
            row=1

            do item=1,size(choices)
                if(.not.datamask(item).and.any(wtsumdist(:,row:row+choices(item)-1)/=0.0_8))&
                    cycle obsloop
                if(datamask(item)) then

                    meanquaditem(:,item,quad)=mmeanth(probcatobsdist(row:row+choices(item)-1),wtsumdist(:,row:row+choices(item)-1))
                    meanquad(:,quad)=meanquad(:,quad)+meanquaditem(:,item,quad)
                    covquaditem(:,:,item,quad)=&
                        mcovth(probcatobsdist(row:row+choices(item)-1),meanquaditem(:,item,quad),&
                        wtsumdist(:,row:row+choices(item)-1))

                    covquad(:,:,quad)=covquad(:,:,quad)+covquaditem(:,:,item,quad)
                    if(quad==1)obswtdist(:,obs)=obswtdist(:,obs)+wtsumdist(:,row+respdist(item))


				end if
                row=row+choices(item)
            end do
            covar=covar+postdensity(quad,obs)*covquad(:,:,quad)

        end if
    end do


!	mean and covariance.

    fitwtdist(:,obs)=mmeanth(postdensity(:,obs),meanquad)


    reswtdist(:,obs)=obswtdist(:,obs)-fitwtdist(:,obs)


    covreswtdist(:,:,obs)=mcovth(postdensity(:,obs),fitwtdist(:,obs),meanquad)+covar
    do row=1,size(wtsumdist,1)
        sdreswtdist(row,obs)=sqrt(covreswtdist(row,row,obs))
        if(sdreswtdist(row,obs)>0.0_8)zreswtdist(row,obs)=reswtdist(row,obs)/sdreswtdist(row,obs)
    end do
end do obsloop
return
end subroutine reswtsumdist

