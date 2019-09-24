!	Get number numdistcodes of distractor recodes per item.
subroutine getnumdistcodes(numdistcodes)
implicit none
interface
    subroutine readnumberdistcodes(numberdistcodes)
        implicit none
        integer,intent(inout)::numberdistcodes(:)
    end subroutine readnumberdistcodes
end interface
integer,intent(out)::numdistcodes(:)
!	The error flag
integer::io
!	The clone.
integer,allocatable::numberdistcodes(:)

allocate(numberdistcodes(size(numdistcodes)),stat=io)
if(io/=0)stop "Allocation failure for distractor recodes."
numberdistcodes=0
call readnumberdistcodes(numberdistcodes)

numberdistcodes=max(numberdistcodes,0)
numdistcodes=numberdistcodes
return
end subroutine getnumdistcodes
subroutine readnumberdistcodes(numberdistcodes)
    implicit none
    integer,intent(inout)::numberdistcodes(:)
    integer::io
    namelist/numberdistrecodes/numberdistcodes
    read(*,nml=numberdistrecodes,iostat=io)
    if(io/=0) stop "Failure to read number of recodes"
    return
end subroutine readnumberdistcodes
