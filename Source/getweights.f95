!	Get integer weights for weighted sum.
!	skillname is the array of factor names.
!	numcatobs is the array of numbers of categories in items.
!	slopedim is the array of dimensions for each factor and item.
!	weightnames are names of weights and weights are weights.
subroutine getweights(skillname,numcatobs,slopedim,weightnames,weights)
implicit none
interface
    subroutine readweights(weightname,weight)
        implicit none
        character(len=32),intent(inout)::weightname
        integer,intent(inout)::weight(:)
    end subroutine readweights
end interface
character(len=32),intent(in)::skillname(:)
integer,intent(in)::numcatobs(:),slopedim(:,:)
character(len=32),intent(out)::weightnames(:)
integer,intent(out)::weights(:,:)
!	buff is a buffer.
!	weightname is used for a nameless group.
character(len=4)::buff
character(len=32)::weightname
!	cat counts categories.
!	count weights.
!	counter finds position in weight.
!	io is error flag.
!	item counts items.
!	
integer::cat,count,counter,io,item
!	weight to be read. 
integer,allocatable::weight(:)

allocate(weight(size(weights,1)),stat=io)
if(io/=0)stop 'Allocation for category weights failed.'

weights=0
if(size(skillname)>=size(weightnames))then
	weightnames=skillname(1:size(weightnames))
	if(size(skillname)>size(weightnames))then
		if(size(weightnames)>1)then
			weightnames(size(weightnames))='Remainder'
		else
			weightnames(1)='Total'
		end if
	end if
	counter=0
	do item=1,size(numcatobs)
		do count=1,size(skillname)
			do cat=1,numcatobs(item)
				if(slopedim(count,item)>0)weights(counter+cat,min(count,size(weightnames)))=cat-1
			end do
			
		end do
		counter=counter+numcatobs(item)		
	end do
else
	weightnames(1:size(skillname))=skillname
	counter=0
	do item=1,size(numcatobs)
		do count=1,size(skillname)
			do cat=1,numcatobs(item)
				if(slopedim(count,item)>0)weights(counter+cat,count)=cat-1
			end do
			
		end do
		counter=counter+numcatobs(item)
	end do
	do count=size(skillname)+1,size(weightnames)
		if(count==size(skillname)+1.and.count>2)then
			weightnames(count)='Total'
			counter=0
			do item=1,size(numcatobs)
				do cat=1,numcatobs(item)
					weights(counter+cat,count)=cat-1
				end do
				counter=counter+numcatobs(item)
			end do
		else
			write(buff,'(i4)') count
			weightnames(count)='Sum_'//trim(adjustl(buff))
		end if
	end do
end if

do count=1,size(weightnames)
	weightname=weightnames(count)
	weight=weights(:,count)
	call readweights(weightname,weight)
	weightnames(count)=weightname
	weights(:,count)=weight
end do


return
end subroutine getweights
subroutine readweights(weightname,weight)
    implicit none
    character(len=32),intent(inout)::weightname
    integer,intent(inout)::weight(:)
    integer::io
    namelist/readweight/weightname,weight
    read(*,nml=readweight,iostat=io)
    if(io/=0)stop "Weights not read successfully."

    return
end subroutine readweights
