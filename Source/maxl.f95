!   Print maximum likelihood estimates of theta or a function of theta.

!    factorname gives names of elements of the latent vector.
!    id is individual id.
!    transname is used if dimtrans is positive.

!    catobsrange is used to relate observed and underlying categories.
!    dat are the responses.
!    dimtrans is the transformation dimension.
!    unitmp is the unit for output.
!    maxita is used to determine the number of iterations to use to find the maximum posterior density.
!    maxitb is used for subiterations.
!    npred is the number of predictors.
!    numcat provides ranges for underlying observations,
!    numcatobs provides ranges of observations, mask indicates items to use for computations.
!    eapid indicates if ids are put out.
!    obsmask is an observation mask that can be used to restrict attention to a portion of the response.
!    beta is the parameter vector.
!    changemin is used for minimum change.
!    lintran is the linear transformation of the latent vector.
!    maxdtheta is the maximum permitted change in the approximation of the maximum in one step.

!    tau is used for minimum step sizes.
!    tol is used to determine the desired accuracy of the maximization.
!    tolsing is used to determine the tolerance for the Cholesky decomposition.
subroutine maxl(factorname,id,transname,catobsrange,dat,dimtrans,maxita,maxitb,&
    npred,numcat,numcatobs,unitml,eapid,obsmask,beta,changemin,&
    lintran,maxdtheta,tau,tol,tolsing)
    implicit none


interface
!    Compute inverse for n by n matrix with modified Cholesky
!    decomposition tri.
    function invert(tri)
        implicit none
        real(kind=8),intent(in)::tri(:,:)
        real(kind=8)::invert(size(tri,1),size(tri,1))
    end function invert
!    Find location of maximum posterior density of the latent vector.
!    catobsrange is used to relate observed and underlying categories.
!    maxita is used to determine the number of iterations to use to find the maximum posterior density.
!    maxitb is used for subiterations.
!    npred is the number of predictors.
!    numcat provides ranges for underlying observations,
!    numcatobs provides ranges of observations, mask indicates items to use for computations.
!    obsmask is an observation mask that can be used to restrict attention to a portion of the response.
!    beta is the parameter vector.
!    changemin is used for minimum change.
!    linth is the linear component of the log density.
!    lintran is the linear transformation of the latent vector.
!    maxdtheta is the maximum permitted change in the approximation of the maximum in one step.
!    quadth is the quadratic component of the log density.
!    resp is the response.
!    tau is used for minimum step sizes.
!    tol is used to determine the desired accuracy of the maximization.
!    tolsing is used to determine the tolerance for the Cholesky decomposition.
!    nhchol is the modified Cholesky decomposition of the negative hessian.
!    theta is the maximum.

    subroutine maxpost(catobsrange,maxita,maxitb,npred,numcat,numcatobs,resp,obsmask,beta,changemin,&
        linth,lintran,maxdtheta,quadth,tau,tol,tolsing,&
        nhchol,theta)
        implicit none
        integer,intent(in)::catobsrange(:,:),maxita,maxitb,npred,numcat(:),numcatobs(:),resp(:)
        logical,intent(in)::obsmask(:)
        real(kind=8),intent(in)::beta(:),changemin,linth(:),lintran(:,:),maxdtheta,quadth(:,:),tau,tol,tolsing
        real(kind=8),intent(inout)::theta(:)
        real(kind=8),intent(out)::nhchol(:,:)
    end subroutine maxpost
!    Scale function
!    dimtrans is the dimension of the transformation.
!    theta is the argument.
    function transform(dimtrans,theta)
        use transset
        implicit none
        integer,intent(in)::dimtrans
        real (kind=8),intent(in)::theta(:)
        real(kind=8)::transform(dimtrans)
    end function transform

end interface


character(len=32),intent(in)::factorname(:),transname(:)
character(len=32),intent(in)::id(:)
integer,intent(in)::catobsrange(:,:),dat(:,:),dimtrans,maxita,maxitb,npred,numcat(:),numcatobs(:),unitml
logical,intent(in)::eapid,obsmask(:)
real(kind=8),intent(in)::beta(:),changemin,lintran(:,:),maxdtheta,tau,tol,tolsing



!	Buffers for writing.
character(len=4)::buffer
character(len=32)::bufferid
character(len=25),allocatable::buffm(:),buffc(:,:)
!	format
character(len=7)::writefmt
!	col is number of entries.
!	io is error flag.
!	obs is observation number.
!   resp is an observation.
!	row is row number.	
integer::col,io,obs,resp(size(dat,1)),row
!   linth is linear component and quadth is quadratic component for maxpost.
!   theta is MLE.
!   tran is transformed MLE
real(kind=8)::inv(size(lintran,2),size(lintran,2)),linth(size(lintran,2)),&
    nhchol(size(lintran,2),size(lintran,2)),quadth(size(lintran,2),size(lintran,2)),&
    theta(size(lintran,2)),tran(dimtrans)
!	Check output dimension.

if(dimtrans>0) then
    col=size(transname)
else
    col=size(factorname)
    allocate(buffc(col,col),stat=io)
    if(io/=0) stop "Allocation failure for arrays for maximum likelihood printing."
end if


allocate(buffm(col),stat=io)
if(io/=0) stop "Allocation failure for arrays for maximum likelihood printing."

if(dimtrans>0)then
    write(unitml,'(a)') "Maximum likelihood estimate"
    write(buffer,'(i4)') 1+2*col
else
    write(unitml,'(a)') "Maximum likelihood estimate and inverse information matrix"
    write(buffer,'(i4)') 1+2*col*(col+1)
end if
writefmt='('//trim(adjustl(buffer))//'a)'
if(dimtrans>0)then
    write(unit=unitml,fmt=writefmt,iostat=io) "ID",(",",trim(adjustl(transname(row))),row=1,col)
    
else
    write(unit=unitml,fmt=writefmt,iostat=io) "ID",&
        (",","Estimate_"//trim(adjustl(factorname(row))),row=1,size(factorname)),&
        ((",","Inverse_"//trim(adjustl(factorname(row)))//"_"//trim(adjustl(factorname(col))),&
        col=1,size(factorname)),row=1,size(factorname))
endif
if(io/=0) stop "Writing of maximum likelihood estimates failed."
!   No prior.
linth=0.0_8
quadth=0.0_8
do obs=1,size(dat,2)
    resp=dat(:,obs)
    if(eapid)then
        bufferid=id(obs)
    else
        write(bufferid,'(i12)') obs
    end if
    theta=0.0_8
    call maxpost(catobsrange,maxita,maxitb,npred,numcat,numcatobs,resp,obsmask,beta,changemin,&
        linth,lintran,maxdtheta,quadth,tau,tol,tolsing,&
        nhchol,theta)
    if(dimtrans>0)then
        tran=transform(dimtrans,theta)
        do row=1,col
                write(buffm(row),'(g25.16e3)')tran(row)
        end do
        write(unit=unitml,fmt=writefmt,iostat=io) trim(adjustl(bufferid)),(",",trim(adjustl(buffm(row))),row=1,col)
    else
        inv=invert(nhchol)
        do row=1,size(factorname)
            write(buffm(row),'(g25.16e3)')theta(row)
            do col=1,size(factorname)
                write(buffc(row,col),'(g25.16e3)')inv(row,col)
            end do
        end do
        do row=1,size(factorname)
            write(buffm(row),'(g25.16e3)')theta(row)
        end do
        write(unit=unitml,fmt=writefmt,iostat=io) &
            trim(adjustl(bufferid)),(",",trim(adjustl(buffm(row))),row=1,size(factorname)),&
            ((",",trim(adjustl(buffc(row,col))),col=1,size(factorname)),row=1,size(factorname))


    endif

	
	if(io/=0) stop "Writing of maximum likelihood estimates failed."
end do
return
end subroutine maxl
