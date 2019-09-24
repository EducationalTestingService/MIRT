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
!	al indicates allocation error.
!	col counts columns.
!	col1 counts columns.
!	counts count strata.
!	locpsu is psu locator. 
!	obs is an observation counter.
!	point indicates the location of a psu.
!	row counts rows.
!	row1 counts rows.
integer::al,col,obs,point,row,row1
!	avegrad is used for weighted summations for grads.
!	avegrad1 is used for weighted summations for grads1.
!	diff is for differences for grads.
!	diff1 is for differences for grads1.

!	totalweight is the sum of weights.


real(kind=8)::avegrad(size(grads,1)),avegrad1(size(grads1,1)),&
	diff(size(grads,1)),diff1(size(grads1,1)),totalweight
	
crossinformation=0.0_8
!	Something must be present.
if(adjust)then
	totalweight=sum(obsweight)
	avegrad=0.0_8
	avegrad1=0.0_8
	do obs=1,size(obsweight)
		if(obsweight(obs)<=0.0_8)cycle
		do row=1,size(grads,1)
			if(grads(row,obs)/=0.0_8)avegrad(row)=avegrad(row)+obsweight(obs)*grads(row,obs)
		end do
		do row=1,size(grads1,1)
			if(grads1(row,obs)/=0.0_8)avegrad1(row)=avegrad1(row)+obsweight(obs)*grads1(row,obs)
		end do
	end do
	avegrad=avegrad/totalweight
	avegrad1=avegrad1/totalweight
end if
do obs=1,size(obsweight)
	if(obsweight(obs)<=0.0_8) cycle
	do row=1,size(grads,1)
		diff(row)=grads(row,obs)
		if(adjust)diff(row)=diff(row)-avegrad(row)
	end do
	do row=1,size(grads1,1)	
		diff1(row)=grads1(row,obs)
		if(adjust)diff1(row)=diff1(row)-avegrad1(row)
	end do
	do row=1,size(grads,1)
		if(diff(row)/=0.0_8)then
			do col=1,size(grads1,1)
				if(diff1(col)/=0.0_8)crossinformation(row,col)&
					=crossinformation(row,col)+obsweight(obs)*diff(row)*diff1(col)
			end do
		end if
	end do
end do

return
end function crossinformation