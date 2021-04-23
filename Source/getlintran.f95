!	Obtain the linear transformation lintran from the underlying
!	latent vector to the latent vector for skills.
!	custom indicates whether the default transformation 
!	is used.
subroutine getlintran(custom,lintran)
implicit none
interface
    subroutine readlintran(lin_tran)
        implicit none
        real(kind=8)::lin_tran(:,:)

    end subroutine readlintran
end interface
logical,intent(in)::custom
real(kind=8),intent(out):: lintran(:,:)

!	diffdim is a dimension difference.
!	io is indicator for an input error.
!	row counts rows.
integer::diffdim,io,row
!	lin_tran is for input of lintran
real(kind=8),allocatable::lin_tran(:,:)


diffdim=size(lintran,2)-size(lintran,1)
lintran=0.0_8
if(diffdim>=0)then
	lintran(:,1:diffdim)=1.0_8
	do row=1,size(lintran,1)
		lintran(row,row+diffdim)=1.0_8
	end do
else
	lintran(1:-diffdim,:)=1.0_8
	do row=1,size(lintran,2)
		lintran(row-diffdim,row)=1.0_8
	end do
end if

if(custom)then
	allocate(lin_tran(size(lintran,1),size(lintran,2)),stat=io)
	if(io/=0) stop "Allocation of linear transformation failed."
	lin_tran=lintran
    call readlintran(lin_tran)
    lintran=lin_tran
end if
return
end subroutine getlintran
subroutine readlintran(lin_tran)
    implicit none
    real(kind=8)::lin_tran(:,:)
    integer::io

    namelist/linearspec/lin_tran
    read(*,nml=linearspec,iostat=io)
    if(io/=0) stop "Custom linear transformationom linear transformation not read successfully."
    return
end subroutine readlintran
