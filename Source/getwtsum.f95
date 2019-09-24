!	Get weighted sum.
!	skillname is array of skill names.
!	catobsrange is array of ranges of underlying categories.
!	numcatobs is array of observed categories per item.
!	slopedim gives factors assigned to items.
!	wtname gives name of sum.
!	wtsum specifies the sum.
subroutine getwtsum(skillname,catobsrange,numcatobs,slopedim,wtname,wtsum)
implicit none
interface
    subroutine readwtsum(weight_name,weight_sum)
        implicit none
        character(len=32),intent(inout)::weight_name
        real(kind=8),intent(inout)::weight_sum(:)
    end subroutine readwtsum
end interface
character(len=32),intent(in)::skillname(:)
integer,intent(in)::catobsrange(:,:),numcatobs(:),slopedim(:,:)
character(len=32),intent(out)::wtname(:)
real(kind=8),intent(out)::wtsum(:,:)
!	Buffer
character(len=4)::buff
!	Clone of wtname
character(len=32)::weight_name
!	cat counts categories.
!	count counts dimensions.
!	counter counts positions in catobsrange.
!	io is the error flag.
!	item counts items.
!	row counts rows.
integer::cat,count,counter,io,item,row
!	Clone of wtsum
real(kind=8),allocatable::weight_sum(:)
allocate(weight_sum(size(wtsum,2)),stat=io)
if(io/=0) stop "Allocation of arrays for weighted sums failed."
wtsum=0.0_8
if(size(skillname)>=size(wtname))then
	wtname=skillname(1:size(wtname))
	if(size(skillname)>size(wtname))then
		if(size(wtname)>1)then
			wtname(size(wtname))='Remainder'
		else
			wtname(1)='Total'
		end if
	end if
	counter=0
	do item=1,size(numcatobs)
		do count=1,size(skillname)
			do cat=1,numcatobs(item)
				if(slopedim(count,item)>0)wtsum(min(count,size(wtname)),catobsrange(1,counter+cat):catobsrange(2,counter+cat))=cat-1
			end do
			
		end do
		counter=counter+numcatobs(item)
	end do
else
	wtname(1:size(skillname))=skillname
	counter=0
	do item=1,size(numcatobs)
		do count=1,size(skillname)
			do cat=1,numcatobs(item)
				if(slopedim(count,item)>0)wtsum(count,catobsrange(1,counter+cat):catobsrange(2,counter+cat))=cat-1
			end do
			
		end do
		counter=counter+numcatobs(item)
	end do
	do count=size(skillname)+1,size(wtname)
		if(count==size(skillname)+1.and.count>2)then
			wtname(count)='Total'
			counter=0
			do item=1,size(numcatobs)
				do cat=1,numcatobs(item)
					wtsum(count,catobsrange(1,counter+cat):catobsrange(2,counter+cat))=cat-1
				end do
				counter=counter+numcatobs(item)
			end do
		else
			write(buff,'(i4)') count
			wtname(count)='Sum_'//trim(adjustl(buff))
		end if
	end do
end if



do row=1,size(wtname)
	weight_name=wtname(row)
	weight_sum=wtsum(row,:)
	call readwtsum(weight_name,weight_sum)

	wtname(row)=weight_name
	wtsum(row,:)=weight_sum

end do
return





end subroutine getwtsum
subroutine readwtsum(weight_name,weight_sum)
    implicit none
    character(len=32),intent(inout)::weight_name
    real(kind=8),intent(inout)::weight_sum(:)
    integer::io
    namelist/weightedsum/weight_name,weight_sum
    read(*,nml=weightedsum,iostat=io)
    if(io/=0) stop "Weighted sum not successfully read."

    return
end subroutine readwtsum
