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

!   fitave is the expected average.
!   fitsum is the expected sum.
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

	function information_complex(npsu,nstratum,psu,stratum,stratify,usepsu,grads,obsweight)
		implicit none
		integer,intent(in)::nstratum
		integer,intent(in),optional::npsu(:),psu(:),stratum(:)
		logical,intent(in)::stratify,usepsu
		real(kind=8),intent(in)::grads(:,:),obsweight(:)
		real(kind=8)::information_complex(size(grads,1),size(grads,1))
	end function information_complex

    subroutine makesym(matrix)
!	Symmetrize n by n matrix based on lower triangle.
        implicit none
        real(kind=8),intent(inout)::matrix(:,:)
    end subroutine makesym


end interface
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
!	al is allocation flag.
!   col is a column
!	dimdesign is the design dimension.
!   dimo is the dimension of the observation vectors.
!	obs counts observations.
!   row is a row.
integer::col,dimo,dimdesign,nobs,obs,row

!   avediff is an average difference.
!   covsum is the cross-product matrix for gradsdes and diff.
!	slopes is for regression of errors on gradients.
!   totsum2 is the square of totsum.

real(kind=8)::avediff(size(observed,1)),diff(size(observed,1),size(obsweight)),&
    covsum(size(observed,1),size(gradsdes,1)),&
    slopes(size(observed,1),size(gradsdes,1)),totsum2

!	Set up parameter arrays for simplified processing.
dimdesign=size(gradsdes,1)
dimo=size(observed,1)


nobs=size(obsweight)
totsum=sum(obsweight,mask=presence)

if(totsum<=0.0_8) then
    obsave=0.0_8
    obssum=0.0_8
    if(stats)then
        covobsave=0.0_8
        covobssum=0.0_8
        sdobsave=0.0_8
        sdobssum=0.0_8
    end if
    if(irel.or.resid)then
        fitave=0.0_8
        fitsum=0.0_8
    end if
    if(resid)then
        aresid=0.0_8

        covresave=0.0_8
        covressum=0.0_8
        resave=0.0_8
        ressum=0.0_8
        sdresave=0.0_8
        sdressum=0.0_8


    end if
    if(irel)then
        covobs=0.0_8
        covres=0.0_8
        sdres=0.0_8
        rel=0.0_8
    end if
    return
end if

!   Basic statistics
totsum2=totsum*totsum
do row=1,dimo
    obssum(row)=sum(obsweight*observed(row,:),mask=presence)

    obsave(row)=obssum(row)/totsum

end do





if(irel.or.resid)then
    do row=1,dimo
        fitsum(row)=sum(obsweight*fitted(row,:),mask=presence)
        fitave(row)=fitsum(row)/totsum
    end do

    avediff=obsave-fitave
    do obs=1,nobs
        if(presence(obs))then
            diff(:,obs)=observed(:,obs)-fitted(:,obs)-avediff
        else
            diff(:,obs)=0.0_8
        end if
    end do

end if
!   Reliability
if(irel)then
    do row=1,dimo
        do col=1,row
            covobs(row,col)=sum(obsweight*(observed(row,:)-obsave(row))&
                *(observed(col,:)-obsave(col)),mask=presence)/totsum
            covres(row,col)=sum(obsweight*diff(row,:)*diff(col,:),mask=presence)/totsum
        end do
        sdobs(row)=sqrt(covobs(row,row))
        sdres(row)=sqrt(covres(row,row))
        if(sdobs(row)>0.0_8)then
            rel(row)=sum(obsweight*(fitted(row,:)-fitave(row)) &
                *(observed(row,:)-obsave(row)),mask=presence)/(totsum*covobs(row,row))
            rel(row)=max(0.0_8,min(1.0_8,rel(row)))
        else
            rel(row)=0.0_8
        end if
    end do
    if(dimo>1)then
        call makesym(covobs)
        call makesym(covres)
    end if


end if
!   Covariance matrices and standard errors for sums and averages.
if(stats)then
    if(.not.complx)then
        if(irel)then
            covobsave=covobs/totsum
            covobssum=totsum*covobs
            sdobsave=sdobs/sqrt(totsum)
            sdobssum=totsum*sdobsave
        else
            do row=1,dimo
                do col=1,row
                    covobssum(row,col)=sum(obsweight*(observed(row,:)-obsave(row))*(observed(col,:)-obsave(col)),mask=presence)
                    covobsave(row,col)=covobssum(row,col)/totsum2
                end do
                sdobssum(row)=sqrt(covobssum(row,row))
                sdobsave(row)=sdobssum(row)/totsum

            end do
            if(dimo>1)then
                call makesym(covobsave)
                call makesym(covobssum)
            end if
        endif
    else
        covobssum=information_complex(npsu,nstratum,psu,stratum,stratify,usepsu,&
            diff,obsweight)
        covobsave=covobssum/totsum2
        do row=1,dimo
            sdobssum(row)=sqrt(covobssum(row,row))
            sdobsave(row)=sdobssum(row)/totsum
        end do
    end if
end if
!   Covariance matrices and standard errors for residuals.
if(resid)then
!   Regression residuals.
    ressum=obssum-fitsum
    resave=ressum/totsum
    covsum=crossinformation(.true.,diff,gradsdes,obsweight)
    slopes=matmul(covsum,eacovgaminv_louis)
    diff=diff-matmul(slopes,gradsdes)
    do row=1,dimo
        avediff(row)=sum(obsweight*diff(row,:))/totsum
        diff(row,:)=diff(row,:)-avediff(row)
    end do
    if(.not.complx)then
        do row=1,dimo
            do col=1,row
                covressum(row,col)=sum(obsweight*diff(row,:)*diff(col,:),mask=presence)
                covresave(row,col)=covressum(row,col)/totsum2
            end do
            sdressum(row)=sqrt(covressum(row,row))
            sdresave(row)=sdressum(row)/totsum
        end do
        if(dimo>1)then
            call makesym(covressum)
            call makesym(covresave)
        end if

    else
        covressum=information_complex(npsu,nstratum,psu,stratum,stratify,usepsu,&
            diff,obsweight)
        covresave=covressum/totsum2
        do row=1,dimo
            sdressum(row)=sqrt(covressum(row,row))
            sdresave(row)=sdressum(row)/totsum
        end do
   end if



    do row=1,dimo
        if(sdressum(row)>tolres*sdobssum(row).and.sdobssum(row)>0.0_8)then
            aresid(row)=ressum(row)/sdressum(row)
        else
            aresid(row)=0.0_8
        end if

    end do
end if


return
end subroutine stdresid
