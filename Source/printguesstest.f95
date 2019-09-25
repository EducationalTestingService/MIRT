!
!	Print guessing tests.
!	itemname contains item names.
!	numcat counts underlying categories.
!	unitguess provides the unit for output.
!	guessres is the residual totals.
!	guessresa is the adjusted residuals.
!	stdguessres is the standard errors of residual totals.


subroutine printguesstest(itemname,numcat,unitguess,&
    guessres,guessresa,stdguessres)
implicit none
character(len=32),intent(in)::itemname(:)
integer,intent(in)::numcat(:),unitguess

real(kind=8),intent(in)::guessres(:),guessresa(:),stdguessres(:)
!	buff are buffers.

character(len=25)::buff(3)


!	io is the flag for output error.
!	item counts items.
!   row counts buffer elements.

integer::io,item,row

write(unit=unitguess,fmt='(a)',iostat=io) 'Guessing tests'
if(io/=0) stop "Printing of marginal distribution failed."
write(unit=unitguess,fmt='(33a)',iostat=io) "Item_name",",","Residual",",","Std_err",",",&
	"Adjusted_residual"


do item=1,size(itemname)
	if(numcat(item)>2)cycle
	write(buff(1),'(g25.16e3)') guessres(item)
	write(buff(2),'(g25.16e3)') stdguessres(item)
	write(buff(3),'(g25.16e3)') guessresa(item)
    write(unit=unitguess,fmt='(7a)',iostat=io)trim(adjustl(itemname(item))),&
		(",",trim(adjustl(buff(row))),row=1,3)
	if(io/=0) stop "Printing of guessing tests failed."
	
	
end do
return
end subroutine printguesstest
