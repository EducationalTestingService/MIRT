!	Obtain distractor recodes.
!	numdistcodes gives number of distractor recodes per item.
!	recodedisttab gives distractor recodes.
subroutine getdistrecodes(numdistcodes,recodedisttab)
implicit none
interface
    subroutine readdistrecodes(recodedist_tab)
        implicit none
integer,intent(inout)::recodedist_tab(:,:)
    end subroutine readdistrecodes
end interface
integer,intent(in)::numdistcodes(:)
integer,intent(out)::recodedisttab(:,:)
!	counter finds recodetab position.
!	io is the error indicator.
!	item counts items
integer::counter,io,item
!	recodedist_tab is the allocation table for an item.
integer,allocatable::recodedist_tab(:,:)

counter=0
do item=1,size(numdistcodes)
	if(numdistcodes(item)>0)then
		allocate(recodedist_tab(2,numdistcodes(item)),stat=io)
		if(io/=0)stop "Allocation error for distractor recodes."
        call readdistrecodes(recodedist_tab)

		recodedisttab(:,counter+1:counter+numdistcodes(item))=&
			recodedist_tab
		deallocate(recodedist_tab)
	end if
	counter=counter+numdistcodes(item)
end do

return
end subroutine getdistrecodes
subroutine readdistrecodes(recodedist_tab)
    implicit none
integer,intent(inout)::recodedist_tab(:,:)
    integer::io
    namelist/recodedisttable/recodedist_tab
    read(*,nml=recodedisttable,iostat=io)
    if(io/=0)stop "Recodes not read successfully."
    return
end subroutine readdistrecodes
