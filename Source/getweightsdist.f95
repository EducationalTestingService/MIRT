!	Get integer weights for weighted sum.
!	skillname is the array of factor names.
!	choices is the array of numbers of distractors of items.
!	distmap is the distractor map.
!   numcatobs counts categories per observed item.
!	slopedim is the array of dimensions for each factor and item.
!	weightnamesdist are names of weights and weights are weights.
subroutine getweightsdist(skillname,choices,distmap,numcatobs,slopedim,weightnamesdist,weightsdist)
implicit none
interface
    subroutine readweightsdist(weightnamedist,weightdist)
        implicit none
        character(len=32),intent(inout)::weightnamedist
        integer,intent(inout)::weightdist(:)
    end subroutine readweightsdist
end interface
character(len=32),intent(in)::skillname(:)
integer,intent(in)::choices(:),distmap(:),numcatobs(:),slopedim(:,:)
character(len=32),intent(out)::weightnamesdist(:)
integer,intent(out)::weightsdist(:,:)
!	buff is a buffer.
!	weightnamedist is used for a nameless group.
character(len=4)::buff
character(len=32)::weightnamedist
!	cat counts categories.
!	count weights.
!	counter finds position in weight.
!	io is error flag.
!	item counts items.
!   row counts rows.
!	
integer::cat,count,counter,io,item,row
!	weight to be read. 
integer,allocatable::weightdist(:)


allocate(weightdist(size(weightsdist,1)),stat=io)
if(io/=0)stop 'Allocation for distractor weights failed.'

weightsdist=0
if(size(skillname)>=size(weightnamesdist))then
	weightnamesdist=skillname(1:size(weightnamesdist))
	if(size(skillname)>size(weightnamesdist))then
		if(size(weightnamesdist)>1)then
			weightnamesdist(size(weightnamesdist))='Remainder'
		else
			weightnamesdist(1)='Total'
		end if
	end if
	counter=0
    row=1
	do item=1,size(choices)
		do count=1,size(skillname)
			do cat=1,choices(item)
				if(slopedim(count,item)>0.and.distmap(counter+cat)>=0)&
                    weightsdist(counter+cat,min(count,size(weightnamesdist)))=distmap(counter+cat)-row
			end do
			
		end do
		counter=counter+choices(item)
        row=row+numcatobs(item)
	end do
else
	weightnamesdist(1:size(skillname))=skillname
	counter=0
    row=1
	do item=1,size(choices)
		do count=1,size(skillname)
			do cat=1,choices(item)
				if(slopedim(count,item)>0.and.distmap(counter+cat)>=0)&
                    weightsdist(counter+cat,count)=distmap(counter+cat)-row
			end do
			
		end do
		counter=counter+choices(item)
        row=row+numcatobs(item)
	end do
	do count=size(skillname)+1,size(weightnamesdist)
		if(count==size(skillname)+1.and.count>2)then
			weightnamesdist(count)='Total'
			counter=0
            row=1
			do item=1,size(choices)
				do cat=1,choices(item)
					if(distmap(counter+cat)>=0)weightsdist(counter+cat,count)=distmap(counter+cat)-row
				end do
				counter=counter+choices(item)
                row=row+numcatobs(item)
			end do
		else
			write(buff,'(i4)') count
			weightnamesdist(count)='Sum_'//trim(adjustl(buff))
		end if
	end do
end if

do count=1,size(weightnamesdist)
	weightnamedist=weightnamesdist(count)
	weightdist=weightsdist(:,count)
	call readweightsdist(weightnamedist,weightdist)
	weightnamesdist(count)=weightnamedist
	weightsdist(:,count)=weightdist
end do


return
end subroutine getweightsdist
subroutine readweightsdist(weightnamedist,weightdist)
    implicit none
    character(len=32),intent(inout)::weightnamedist
    integer,intent(inout)::weightdist(:)
    integer::io
    namelist/readweightdist/weightnamedist,weightdist
    read(*,nml=readweightdist,iostat=io)
    if(io/=0)stop "Weights not read successfully."

    return
end subroutine readweightsdist
