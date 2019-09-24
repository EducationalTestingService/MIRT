!	Constraint specification.
!	constmat is the matrix .
!	constvec is the vector.
!	constmat(gamma) is to equal constvec.
subroutine getconstraint(constmat,constvec)
implicit none
interface
    subroutine readconstraint(const_mat,const_vec)
        implicit none
real(kind=8),intent(inout)::const_mat(:,:),const_vec(:)
    end subroutine readconstraint
end interface
real(kind=8),intent(out)::constmat(:,:),constvec(:)
!	Error flag and counter.
integer::io,i
!	Clones.
real(kind=8),allocatable::const_mat(:,:),const_vec(:)

allocate(const_mat(size(constmat,1),size(constmat,2)),const_vec(size(constvec)),stat=io)
if(io/=0)stop "Constraint allocation failed."
const_mat=0.0_8
do i=1,min(size(constmat,1),size(constmat,2))
	const_mat(i,i)=1.0_8
end do
const_vec=0.0_8
call readconstraint(const_mat,const_vec)
constmat=const_mat
constvec=const_vec
return

end subroutine getconstraint		
subroutine readconstraint(const_mat,const_vec)
    implicit none
    real(kind=8),intent(inout)::const_mat(:,:),const_vec(:)
    integer::io
    namelist/constraints/const_mat,const_vec
    read(*,nml=constraints,iostat=io)
    if(io/=0) stop "Constraints not read successfully."
    return
end subroutine readconstraint
