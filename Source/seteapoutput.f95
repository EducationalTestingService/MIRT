!	set parameters for eap output.
!	dimtrans is the dimension of the transformation of the latent vector.
!	dimwtsum is the dimension of the desired linear combination of observations.
!	dimwtsum is the dimension of the desired linear combination of distractor observations.
!	eapmask indicates the items to include.
!	altbeta indicates the beta vector to use.
subroutine seteapoutput(dimtrans,dimwtsum,dimwtsumdist,eapmask,altbeta)
implicit none

integer,intent(out)::dimtrans,dimwtsum,dimwtsumdist
logical,intent(out)::eapmask(:)
real(kind=8),intent(inout)::altbeta(:)
!   The input file.
character(len=256)::betafile
character(len=64)::paramname
!   i is a counter.
!	io error flag
!   skip is amount to skip in betafile.
integer::i,io,skip
!	eapmask clone
logical::usebetafile
logical,allocatable::eap_mask(:)
!	altbeta clone
real(kind=8),allocatable::alt_beta(:)
namelist/eapoutput/betafile,dimtrans,dimwtsum,dimwtsumdist,skip,eap_mask,usebetafile,alt_beta
dimtrans=0
dimwtsum=0
dimwtsumdist=0
allocate(eap_mask(size(eapmask)),alt_beta(size(altbeta)),stat=io)
if(io/=0) stop "Allocation of EAP input specification failed."
betafile='betas.csv'
eap_mask=.true.
usebetafile=.false.
skip=2
alt_beta=altbeta

read(*,nml=eapoutput,iostat=io)
if(io/=0)stop "EAP specifications in error."
if(usebetafile) then
    open(8,file=betafile,iostat=io)

    if(io/=0) stop "Gamma parameters not successfully read."
    skip=max(0,skip)
    do i=1,skip

        read(8,*)
        if(io/=0) stop "Beta parameters not successfully read."

    end do

    do i=1,size(altbeta)

        read(8,fmt=*,iostat=io) paramname,alt_beta(i)

        if(io/=0) stop "Beta parameters not successfully read."
    end do



    close(8)
end if

if(dimtrans<0) dimtrans=0
if(dimwtsum<0) dimwtsum=0
if(dimwtsumdist<0) dimwtsumdist=0

eapmask=eap_mask
altbeta=alt_beta
return
end subroutine seteapoutput
