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
interface
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
end interface
integer,intent(in)::catobsrange(:,:),dat(:,:),numcat(:),numcatobs(:)


real(kind=8),intent(in)::beta(:),lintran(:,:),postdensity(:,:),theta(:,:,:),wtsum(:,:)
real(kind=8),intent(out)::covreswt(:,:,:),fitwt(:,:),obswt(:,:),reswt(:,:),sdreswt(:,:),zreswt(:,:)
!	item is an item.
!	ncat is the total number of underlying categories.
!	obs is observation number.
!   quad counts quadrature points.
!	resp is a single response.
!	row is a row counter.
!   row1 is a second row counter.

integer::item,ncat,obs,quad,resp(size(numcatobs)),row
!   datamask is for identification of item responses for observations.
logical::datamask(size(numcat))
!   condprob is for conditional probabilities of underlying responses given observed ones.
!   covar is the covariance matrix of the differences of observed and duplicate sums.
!   covquad is the covariance matrix of sums given the latent vector.
!   covquaditem is the covariance matrix of weighted sums for items given the latent vector.
!	locations is for location parameters.

!	meanquad is for means of sums for quadrature points.
!	meanquaditem is for means of sums for items and for quadrature points.
!	obsquad is for conditional means of sums for quadrature points.
!	obsquaditem is for conditional means of sums for items and for quadrature points.
!	newtheta is lintran(theta).
!	probcat is used for underlying marginal conditional probabilities.
!   probcatobs is used for marginal conditional probabilities of observations.
!	scales is for scale parameters.

real(kind=8)::condprob(sum(numcat)),covar(size(wtsum,1),size(wtsum,1))
real(kind=8)::covquad(size(wtsum,1),size(wtsum,1),size(postdensity,1)),&
    covquaditem(size(wtsum,1),size(wtsum,1),size(numcat),size(postdensity,1))
real(kind=8)::locations(sum(numcat))
real(kind=8)::meanquad(size(wtsum,1),size(postdensity,1)),meanquaditem(size(wtsum,1),size(numcat),size(postdensity,1))
real(kind=8)::newtheta(size(lintran,1))
real(kind=8)::obsquad(size(wtsum,1),size(postdensity,1)),obsquaditem(size(wtsum,1),size(numcat),size(postdensity,1))
real(kind=8)::probcat(sum(numcat)),probcatobs(size(numcatobs))
real(kind=8)::scales(size(lintran,1),sum(numcat))





!	Set up parameter arrays for simplified processing.
ncat=sum(numcat)
!	Item scale parameters.
locations=beta(1:ncat)
scales=reshape(beta(ncat+1:ncat*(1+size(lintran,1))),(/size(lintran,1),ncat/))
covar=0.0_8
covquad=0.0_8
covquaditem=0.0_8
covreswt=0.0_8
fitwt=0.0_8
obswt=0.0_8
reswt=0.0_8
sdreswt=0.0_8
zreswt=0.0_8



!	Obtain results for each observation.
obsloop:do obs=1,size(postdensity,2)

    resp=dat(:,obs)

    datamask=.false.
    do item=1,size(numcat)
        if(resp(item)>=0.and.resp(item)<numcatobs(item))datamask(item)=.true.
    end do

!	Cycle through quadrature points.
    covquad=0.0_8
    covquaditem=0.0_8
    meanquad=0.0_8
    meanquaditem=0.0_8
    obsquad=0.0_8
    obsquaditem=0.0_8

	covar=0.0_8
    do quad=1,size(postdensity,1)
        newtheta=matmul(lintran,theta(:,quad,obs))


! Probability computation.


        probcat=probvec(numcat,datamask,locations,scales,newtheta)
        probcatobs=probvecobs(catobsrange,numcatobs,resp,datamask,probcat)
        condprob=condprobvec(catobsrange,numcatobs,resp,datamask,probcat,probcatobs)
        if(size(wtsum,1)>0)then
            row=1
            

            do item=1,size(numcat)
                if(.not.datamask(item).and.any(wtsum(:,row:row+numcat(item)-1)/=0.0_8))&
                    cycle obsloop
                if(datamask(item)) then
                    meanquaditem(:,item,quad)=mmeanth(probcat(row:row+numcat(item)-1),wtsum(:,row:row+numcat(item)-1))
                    meanquad(:,quad)=meanquad(:,quad)+meanquaditem(:,item,quad)
                    obsquaditem(:,item,quad)=mmeanth(condprob(row:row+numcat(item)-1),wtsum(:,row:row+numcat(item)-1))
                    obsquad(:,quad)=obsquad(:,quad)+obsquaditem(:,item,quad)
                    covquaditem(:,:,item,quad)=&
                        mcovth(probcat(row:row+numcat(item)-1),meanquaditem(:,item,quad),wtsum(:,row:row+numcat(item)-1))&
                        +mcovth(condprob(row:row+numcat(item)-1),obsquaditem(:,item,quad),wtsum(:,row:row+numcat(item)-1))
                    covquad(:,:,quad)=covquad(:,:,quad)+covquaditem(:,:,item,quad)


                end if

                row=row+numcat(item)
            end do
            covar=covar+postdensity(quad,obs)*covquad(:,:,quad)

        end if
    end do

!	mean and covariance.
    fitwt(:,obs)=mmeanth(postdensity(:,obs),meanquad)
    obswt(:,obs)=mmeanth(postdensity(:,obs),obsquad)

    reswt(:,obs)=obswt(:,obs)-fitwt(:,obs)
    covreswt(:,:,obs)=mcovth(postdensity(:,obs),fitwt(:,obs),meanquad)+covar
    do row=1,size(wtsum,1)
        sdreswt(row,obs)=sqrt(covreswt(row,row,obs))
        if(sdreswt(row,obs)>0.0_8)zreswt(row,obs)=reswt(row,obs)/sdreswt(row,obs)
    end do
end do obsloop
return
end subroutine reswtsum
