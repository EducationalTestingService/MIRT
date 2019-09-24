!	Get number numcodes of recodes per item.
subroutine getnumcodes(numcodes)
implicit none
interface
    subroutine readnumbercodes(numbercodes)
        implicit none
        integer,intent(inout)::numbercodes(:)

    end subroutine readnumbercodes
end interface
integer,intent(out)::numcodes(:)
!	The error flag
integer::io
!	The clone.
integer,allocatable::numbercodes(:)

allocate(numbercodes(size(numcodes)),stat=io)
if(io/=0)stop "Allocation failure for recodes."
numbercodes=0
call readnumbercodes(numbercodes)
numbercodes=max(numbercodes,0)
numcodes=numbercodes
return
end subroutine getnumcodes
subroutine readnumbercodes(numbercodes)
    implicit none
    integer,intent(inout)::numbercodes(:)
    integer::io
    namelist/numberrecodes/numbercodes
    read(*,nml=numberrecodes,iostat=io)
    if(io/=0) stop "Failure to read number of recodes"
    return
end subroutine readnumbercodes
