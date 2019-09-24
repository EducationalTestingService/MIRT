!	Estimate the information matrix with the Louis approach.
!	adjust is .true. if the average of grads is to be subtracted.
!	grads is the array of observation gradients.
!	obsweight contains observation weights.
function information(adjust,grads,obsweight)
implicit none
!
interface
	subroutine makesym(matrix)
!	Symmetrize n by n matrix based on lower triangle.
		real(kind=8),intent(inout)::matrix(:,:)
	end subroutine makesym
end interface

logical,intent(in)::adjust
real(kind=8),intent(in)::grads(:,:),obsweight(:)
real(kind=8)::information(size(grads,1),size(grads,1))
!	avegrad is the average  gradient and diff is the deviation between the gradient and
!	the average gradient.
real(kind=8)::avegrad(size(grads,1)),diff(size(grads,1))
!	col is a column number.
!	obs is an observation counter.
!	row is a row  number.
integer::col,obs,row
if(adjust)then
	avegrad=0.0_8
	do obs=1,size(obsweight)
		if(obsweight(obs)<=0.0_8)cycle
		do row=1,size(grads,1)
			if(grads(row,obs)/=0.0_8)avegrad(row)=avegrad(row)+obsweight(obs)*grads(row,obs)
		end do
	end do
	avegrad=avegrad/sum(obsweight)
end if
information=0.0_8
do obs=1,size(obsweight)
	if(obsweight(obs)<=0.0_8) cycle
	do row=1,size(grads,1)
		diff(row)=grads(row,obs)
		if(adjust)diff(row)=diff(row)-avegrad(row)
		if(diff(row)/=0.0_8)then
			do col=1,row
				if(diff(col)/=0.0_8)information(row,col)=information(row,col)+obsweight(obs)*diff(row)*diff(col)
			end do
		end if
	end do
end do
call makesym(information)
return
end
