!       Save alpha and cholnhess on unit specified by unitalpha
!       alpha is location of the maximum of the conditional probability
!	of a response given the latent vector.
!       cholnhess is the modified Cholesky decomposition of the
!       negative hessian of the log conditional probability with respect
!       to the latent vector at the maximum.
subroutine savealpha(unitalpha,alpha,cholnhess)
implicit none

integer,intent(in)::unitalpha
REAL(KIND=8),INTENT(in)::alpha(:,:),cholnhess(:,:,:)
!	buffers.
character(len=4)::buff
character(len=25)::buffer(size(alpha,1)*(size(alpha,1)+1))
!	writefmt is used for producing the output.
character(len=7)::writefmt
!	col is used to count columns.
!	counter is a counter.
!       io checks for output errors.
!       obs counts observations.
!       row is used to count rows.
integer::col,counter,io,obs,row

write(buff,'(i4)') 2*size(alpha,1)*(size(alpha,1)+1)-1
writefmt='('//trim(adjustl(buff))//'a)'
do obs=1,size(alpha,2)
	counter=size(alpha,1)+1
	do row=1,size(alpha,1)
		write(buffer(row),'(g25.15e3)') alpha(row,obs)
	end do
	do col=1,size(alpha,1)
		do row=1,size(alpha,1)
			write(buffer(counter),'(g25.16e3)') cholnhess(row,col,obs)
			counter=counter+1
		end do
	end do
	write(unitalpha,writefmt) (trim(adjustl(buffer(row))),',',row=1,size(buffer)-1),buffer(size(buffer))
       	
end do
return

end subroutine savealpha
