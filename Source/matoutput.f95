!	print square matrix.
!	comment is the type of matrix.
!	rowname provides row names.
!	unitmat is the unit number to use.
!	mat provides the matrix.

subroutine matoutput(comment,rowname,unitmat,mat)
implicit none
character(len=*),intent(in)::comment
character(len=*),intent(in)::rowname(:)
integer,intent(in)::unitmat
real(kind=8),intent(in)::mat(:,:)
!	buffers
character(len=4)::buffer
character(len=25)::buff(size(rowname))
!	write format
character(len=9)::writefmt
!	col counts columns.
!	io indicates write error.
!	row counts rows.
integer::col,io,row
write(buffer,'(i4)')2*size(rowname)
writefmt='(a,'//trim(adjustl(buffer))//'a)'
write(unitmat,'(a)') comment
write(unit=unitmat,fmt=writefmt,iostat=io) "Row name",(",",trim(adjustl(rowname(col))),col=1,size(rowname))
if(io/=0) stop "Failure to write matrix."
do row=1,size(rowname)
	do col=1,size(rowname)
		write(buff(col),'(g25.16e3)') mat(row,col)
	end do
	write(unit=unitmat,fmt=writefmt,iostat=io)&
				trim(adjustl(rowname(row))),&
				(",",trim(adjustl(buff(col))),&
				col=1,size(rowname))
	if(io/=0) stop "Failure to write matrix."
end do		
return
end subroutine matoutput

