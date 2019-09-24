!       Initialize alpha and cholnhess
!       alpha is location of the maximum of the conditional probability
!	of a response given the latent vector.
!       cholnhess is the modified Cholesky decomposition of the
!       negative hessian of the log conditional probability with respect
!       to the latent vector at the maximum.
subroutine getalpha(alpha,cholnhess)
implicit none
real(kind=8),intent(out)::alpha(:,:),cholnhess(:,:,:)
!	fileformat is the format for reading.
!	filename is the filename for reading.
character(len=256)::fileformat,filename
!       io checks for input errors.
!       obs counts observations.
!       row is used to count rows of cholnhess.
integer::io,obs,row
!       readalpha checks if alpha and cholnhess are to be read.
logical::readalpha
namelist/inputinformation/fileformat,filename,readalpha
fileformat="*"
filename="alpha.csv"
alpha=0_8
cholnhess=0.0_8
readalpha=.false.
do row=1,size(alpha,1)
        cholnhess(row,row,:)=1.0_8
end do
read(*,nml=inputinformation,iostat=io)
if(io/=0.or.(filename=="".and.readalpha)) stop &
	 "Input settings for initial values for adaptive quadrature not read successfully."
IF(readalpha)then
        open(unit=8,file=filename,iostat=io)
        IF(io/=0) stop "Input file for initial values for adaptive quadrature not opened succesfully."

       	do obs=1,size(alpha,2)
       		if(fileformat=="*")then
                	read(unit=8,fmt=*,iostat=io) alpha(:,obs),cholnhess(:,:,obs)
                else
                	read(unit=8,fmt=fileformat,iostat=io) alpha(:,obs),cholnhess(:,:,obs)
                end if
                if(io/=0) stop "Input data for initial values for adaptive quadrature not read successfully."
        end do
	close(8)
end if

return

end subroutine getalpha
