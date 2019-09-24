!
!	Get item scores.  Use numcatobs to find dimension and defaults.
!	numcatobs is the number of observed categories per item
function getscores(numcatobs)
implicit none
integer,intent(in)::numcatobs(:)
real(kind=8)::getscores(sum(numcatobs))
!	Error flag is io and category counter is cat.  Item count is item.  Position is placement of
!	entry for item score for a specific item and category.
interface
    subroutine readscores(scores)
        implicit none
        real(kind=8),intent(inout)::scores(:)
    end subroutine readscores
end interface
integer::cat,io,item,position
!	copy of getscores.
real(kind=8),allocatable::scores(:)

!	Default setting.
position=1
allocate(scores(sum(numcatobs)),stat=io)
if(io/=0) stop 'Allocation failed for item scores.'
do item=1,size(numcatobs)
	do cat=1,numcatobs(item)
		scores(position)=cat-1
		position=position+1
	end do
end do
call readscores(scores)

getscores=scores
return
end function getscores
subroutine readscores(scores)
implicit none
integer::io
real(kind=8),intent(inout)::scores(:)
namelist/itemscores/scores
read(*,nml=itemscores,iostat=io)
if(io/=0) stop 'Item scores not read successfully.'
return
end subroutine readscores
