!	Obtain recodes.
!	numcodes gives number of recodes per item.
!	recodetab gives recodes.
subroutine getrecodes(numcodes,recodetab)
implicit none
interface
    subroutine readrecodetab(recode_tab)
        implicit none
        integer,intent(inout)::recode_tab(:,:)
    end subroutine readrecodetab
end interface
integer,intent(in)::numcodes(:)
integer,intent(out)::recodetab(:,:)
!	counter finds recodetab position.
!	io is the error indicator.
!	item counts items
integer::counter,io,item
!	recode_tab is the allocation table for an item.
integer,allocatable::recode_tab(:,:)

counter=0
do item=1,size(numcodes)
	if(numcodes(item)>0)then
		allocate(recode_tab(2,numcodes(item)),stat=io)
		if(io/=0)stop "Allocation error for recodes."
        call readrecodetab(recode_tab)

		recodetab(:,counter+1:counter+numcodes(item))=&
			recode_tab
		deallocate(recode_tab)
	end if
	counter=counter+numcodes(item)
end do

return
end subroutine getrecodes
subroutine readrecodetab(recode_tab)
    implicit none
    integer,intent(inout)::recode_tab(:,:)
    integer::io
    namelist/recodetable/recode_tab
    read(*,nml=recodetable,iostat=io)
    if(io/=0)stop "Recodes not read successfully."
    return
end subroutine readrecodetab
