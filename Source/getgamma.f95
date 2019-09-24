!	Initialize gamma.  Read any settings.
subroutine getgamma(gamma)
implicit none
interface
    subroutine readgammas(usegammafile,gammafile,skip,gammas)
        character(len=256),intent(inout)::gammafile
        integer,intent(inout)::skip
        logical,intent(inout)::usegammafile
        real(kind=8),intent(inout)::gammas(:)

    end subroutine readgammas
end interface
!	The input file and input format.
character(len=256)::gammafile
character(len=64)::paramname
!	io is the flag for read error.
integer::i,io,skip
!	Flag for using external file.
logical::usegammafile
real(kind=8),intent(inout)::gamma(:)
real(kind=8),allocatable::gammas(:)

allocate(gammas(size(gamma)),stat=io)
if(io/=0) stop "Gamma parameters not allocated for reading."
gammafile='gammas.csv'

usegammafile=.FALSE.
gammas=gamma
skip=2
call readgammas(usegammafile,gammafile,skip,gammas)

if(usegammafile) then
	
	open(8,file=gammafile,iostat=io)

    if(io/=0) stop "Gamma parameters not successfully read."
    skip=max(0,skip)
    do i=1,skip

        read(8,*,iostat=io)
        if(io/=0) stop "Gamma parameters not successfully read."

    end do
	
	do i=1,size(gamma)

        read(8,fmt=*,iostat=io) paramname,gammas(i)

        if(io/=0) stop "Gamma parameters not successfully read."
    end do
	

	
	close(8)
end if
gamma=gammas
return
end subroutine getgamma
subroutine readgammas(usegammafile,gammafile,skip,gammas)
    character(len=256),intent(inout)::gammafile
    integer,intent(inout)::skip
    logical,intent(inout)::usegammafile
    real(kind=8),intent(inout)::gammas(:)
    integer::io
    namelist/readgamma/usegammafile,gammafile,skip,gammas

    read(*,nml=readgamma,iostat=io)
    if(io/=0) stop "Gamma parameters not successfully read."
    return
end subroutine readgammas
