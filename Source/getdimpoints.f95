!	Obtain number of quadrature points for each dimension.
!	The result is placed in dimpoints.
!	The default for quadrature points per dimension is 5.

subroutine getdimpoints(dimpoints)
implicit none
interface
    subroutine readdimpoints(Q)
        implicit none
        integer,intent(inout)::Q(:)
    end subroutine readdimpoints
end interface
integer,intent(out)::dimpoints(:)
!	io is an input flag.
!	row is a counter.
integer::io,row
!	Q is used for dimpoints.
integer,allocatable::Q(:)

allocate(Q(size(dimpoints)),stat=io)
if(io/=0)stop "Allocation for quadrature points failed."
Q=5
call readdimpoints(Q)
dimpoints=Q
do row=1,size(dimpoints)
	if(dimpoints(row)<2)dimpoints(row)=5
end do



return
end subroutine getdimpoints

subroutine readdimpoints(Q)
    implicit none
    integer,intent(inout)::Q(:)
    integer::io
    namelist/nquadperdim/Q
    read(*,nml=nquadperdim,iostat=io)
    if(io/=0) stop "Quadrature points per dimension not read properly."
    return
end subroutine readdimpoints
