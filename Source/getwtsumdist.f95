!	Get weighted sum for response choices.
!	skillname is array of skill names.
!	choices are counts for choices for each item.
!	distmap is the distractor map.
!   numcatobs counts categories per observed item.
!	slopedim gives factors assigned to items.
!	wtnamedist gives name of sum.
!	wtsumdist specifies the sum.
subroutine getwtsumdist(skillname,choices,distmap,numcatobs,slopedim,wtnamedist,wtsumdist)
implicit none
interface
    subroutine readwtsumdist(weight_name_dist,weight_sum_dist)
        implicit none
        character(len=32)::weight_name_dist
        real(kind=8)::weight_sum_dist(:)
    end subroutine readwtsumdist
end interface
character(len=32),intent(in)::skillname(:)
integer,intent(in)::choices(:),distmap(:),numcatobs(:),slopedim(:,:)
character(len=32),intent(out)::wtnamedist(:)
real(kind=8),intent(out)::wtsumdist(:,:)
!	Buffer
character(len=4)::buff
!	Clone of wtnamedist
character(len=32)::weight_name_dist
!	cat counts categories.
!	count counts dimensions.
!	counter counts positions in catobsrange.
!	io is the error flag.
!	item counts items.
!	row counts rows.
integer::cat,count,counter,io,item,row
!	Clone of wtsumdist
real(kind=8),allocatable::weight_sum_dist(:)
allocate(weight_sum_dist(size(wtsumdist,2)),stat=io)
if(io/=0) stop "Allocation of arrays for weighted sums of distractors failed."
wtsumdist=0.0_8
if(size(skillname)>=size(wtnamedist))then
	wtnamedist=skillname(1:size(wtnamedist))
	if(size(skillname)>size(wtnamedist))then
		if(size(wtnamedist)>1)then
			wtnamedist(size(wtnamedist))='Remainder'
		else
			wtnamedist(1)='Total'
		end if
	end if
	counter=0
    row=1
	do item=1,size(choices)
		do count=1,size(skillname)
			do cat=1,choices(item)
				if(slopedim(count,item)>0.and.distmap(counter+cat)>=0)&
					wtsumdist(min(count,size(wtnamedist)),counter+cat)=distmap(counter+cat)-row
			end do
			
		end do
		counter=counter+choices(item)
        row=row+numcatobs(item)
	end do
else
	wtnamedist(1:size(skillname))=skillname
	counter=0
    row=1
	do item=1,size(choices)
		do count=1,size(skillname)
			do cat=1,choices(item)
				if(slopedim(count,item)>0.and.distmap(counter+cat)>=0)&
					wtsumdist(count,counter+cat)=distmap(counter+cat)-row
			end do
			
		end do
		counter=counter+choices(item)
        row=row+numcatobs(item)
	end do
	do count=size(skillname)+1,size(wtnamedist)
		if(count==size(skillname)+1.and.count>2)then
			wtnamedist(count)='Total'
			counter=0
            row=1
			do item=1,size(choices)
				do cat=1,choices(item)
					if(distmap(counter+cat)>=0)wtsumdist(count,counter+cat)=distmap(counter+cat)-row
				end do
				counter=counter+choices(item)
                row=row+numcatobs(item)
			end do
		else
			write(buff,'(i4)') count
			wtnamedist(count)='Sum_'//trim(adjustl(buff))
		end if
	end do
end if



do row=1,size(wtnamedist)
	weight_name_dist=wtnamedist(row)
	weight_sum_dist=wtsumdist(row,:)
    call readwtsumdist(weight_name_dist,weight_sum_dist)


	wtnamedist(row)=weight_name_dist
	wtsumdist(row,:)=weight_sum_dist

end do
return





end subroutine getwtsumdist
subroutine readwtsumdist(weight_name_dist,weight_sum_dist)
    implicit none
    character(len=32)::weight_name_dist
    real(kind=8)::weight_sum_dist(:)
    integer::io
    namelist/weightedsumdist/weight_name_dist,weight_sum_dist

    read(*,nml=weightedsumdist,iostat=io)
    if(io/=0) stop "Weighted distractor sum not successfully read."

    return
end subroutine readwtsumdist
