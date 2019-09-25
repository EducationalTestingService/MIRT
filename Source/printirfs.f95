!
!	Print estimated item response functions for selected points.
!	factorname
!	itemname contains item names.
!	numcatobs counts observed categories.
!	unitirf provides the unit for output.
!	resid indicates if residuals are displayed.
!	fitirf is the unconditional estimated item response functions.
!	obsirf is the conditional estimated item response functions.
!	presentedirf counts weighted items presented.
!	residairf is the adjusted residuals.
!	residirf is the residuals.
!	thetacheck is the array of selected points

subroutine printirfs(factorname,itemname,numcatobs,unitirf,resid,&
	fitirf,obsirf,presentedirf,&
	residairf,residirf,stdresidirf,thetacheck)
implicit none
character(len=32),intent(in)::factorname(:)
character(len=32),intent(in)::itemname(:)

integer,intent(in)::numcatobs(:),unitirf
logical,intent(in)::resid
real(kind=8),intent(in)::fitirf(:,:),obsirf(:,:),&
		presentedirf(:),residairf(:,:),residirf(:,:),&
		stdresidirf(:,:),thetacheck(:,:)
!	buff1 and buff are buffers.  fmt1 is for formats.
character(len=12)::buff1
character(len=25),allocatable::buff(:)
character(len=25)::fmt1
!	cat counts categories.
!	counter counts positions.
!	io is the flag for output error.
!	item counts items.
!	point counts quadrature points.
!	position counts positions.
!	row and row1 count rows.
integer::cat,counter,io,item,point,position,row,row1

counter=1
row=3+size(factorname)
position=row-1
if(resid)row=row+3
allocate(buff(row),stat=io)
if(io/=0)stop 'Allocation failed for printing item response function.'

write(fmt1,'(a,i2,a)') '(',row+row+3,'a)'
write(unit=unitirf,fmt=fmt1,iostat=io) 'Item response functions'
if(io/=0) stop "Printing of item response functions failed."
if(resid)then
	write(unit=unitirf,fmt=fmt1,iostat=io) "Item_name",",",&
		"Item_freq",",","Cat_no",",",(trim(adjustl(factorname(row1))),',',row1=1,size(factorname)),"Cond_irf",",",&
		"Uncond_irf",",",&
		"Res_irf",",","Std_err_res_irf",",",&
		"Adjusted_residual"
else
	write(unit=unitirf,fmt=fmt1,iostat=io) "Item_name",",",&
		"Item_freq",",","Cat_no",",",(trim(adjustl(factorname(row1))),',',row1=1,size(factorname)),"Cond_irf",",",&
		"Uncond_irf"
end if

do item=1,size(itemname)
	do cat=1,numcatobs(item)
		write(buff(1),'(g25.16e3)') presentedirf(item)
		write(buff1,'(i12)') cat-1
		
		do point=1,size(thetacheck,2)
			do row1=1,size(factorname)
				write(buff(row1+1),'(g25.16e3)') thetacheck(row1,point)
			end do
			write(buff(position),'(g25.16e3)') obsirf(counter,point)
			write(buff(position+1),'(g25.16e3)') fitirf(counter,point)
			if(resid)then
				write(buff(position+2),'(g25.16e3)') residirf(counter,point)
				write(buff(position+3),'(g25.16e3)') stdresidirf(counter,point)
				write(buff(position+4),'(g25.16e3)') residairf(counter,point)
			end if
			write(unit=unitirf,fmt=fmt1,iostat=io)trim(adjustl(itemname(item))),&
				",",trim(adjustl(buff(1))),",",trim(adjustl(buff1)),&
				(",",trim(adjustl(buff(row1))),row1=2,row)
		
			if(io/=0) stop "Printing of marginal distribution failed."
			
		end do
		counter=counter+1
	end do
end do
return
end subroutine printirfs
