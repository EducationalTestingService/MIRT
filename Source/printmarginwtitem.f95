!
!	Print weighted sum totals by item categories.
!	itemname contains item names.
!	wtname is the weight name.
!	numcatobs counts observed categories.
!	unitmarginwtitem provides the unit for output.
!	resid indicates if residuals are displayed.
!	fitpwtitem is the  fitted average of products.
!	fitwtitem is the fitted total of products.
!	obspwtitem is the observed average of products.
!	obsmargwtitem is the observed total of products.
!	presentedwtitem counts weighted items presented with items.
!	residawtitem is the adjusted residual.
!	residpwtitem is the residual for the average of products.
!	residwtitem is the residual for the total of products.
!	stdobspwtitem is the standard deviation of the observed average of products.
!	stdobswtitem is the asymptotic standard deviation of the observed 
!		total of products.
!	stdresidpwtitem is the asymptotic standard error for the residual for average
!       of products.
!	stdresidwtitem is the asymptotic standard error for the residual for total
!       of products.


subroutine printmarginwtitem(itemname,wtname,numcatobs,unitmarginwtitem,&
    resid,fitpwtitem,fitwtitem,obspwtitem,obswtitem,presentedwtitem,&
    residawtitem,residpwtitem,residwtitem,&
    stdobspwtitem,stdobswtitem,stdresidpwtitem,stdresidwtitem)
implicit none
character(len=*),intent(in)::itemname(:)
character(len=*),intent(in)::wtname
integer,intent(in)::numcatobs(:),unitmarginwtitem
logical,intent(in)::resid
real(kind=8),intent(in)::fitpwtitem(:),fitwtitem(:),&
    obspwtitem(:),obswtitem(:),presentedwtitem(:),&
    residawtitem(:),residpwtitem(:),residwtitem(:),&
    stdobspwtitem(:),stdobswtitem(:),stdresidpwtitem(:),stdresidwtitem(:)
!	buff1 and buff are buffers.
character(len=12)::buff1
character(len=25)::buff(15)

!	cat counts categories.
!	counter counts positions.
!	io is the flag for output error.
!	item counts items.
!	row counts rows.
integer::cat,counter,io,item,row
!	var is a variance.
real(kind=8)::var
counter=1
write(unit=unitmarginwtitem,fmt='(2a)',iostat=io) 'Totals and averages of products for ',wtname

if(io/=0) stop "Printing of marginal distribution failed."
if(resid)then
	write(unit=unitmarginwtitem,fmt='(27a)',iostat=io) "Item_name",",",&
		"Cat_no",",","Item_freq",",","Obs_tot_prod",",","Std_err_obs_tot_prod",",",&
		"Fit_tot_prod",",","Obs_ave_prod",",","Std_err_obs_ave_prod",",",&
		"Fit_ave_prod",",",&
		"Res_tot_prod",",","Std_err_res_tot_prod",",",&
		"Res_ave_prod",",","Std_err_res_ave_prod",",",&
		"Adjusted_residual"
else
	write(unit=unitmarginwtitem,fmt='(17a)',iostat=io) "Item_name",",",&
        "Cat_no",",","Item_freq",",","Obs_tot_prod",",",&
        "Std_err_obs_tot_prod",",",&
        "Fit_tot_prod",",","Obs_ave_prod",",","Std_err_obs_ave_prod",",",&
        "Fit_ave_prod"
end if

do item=1,size(itemname)
	do cat=1,numcatobs(item)
		write(buff1,'(i12)') cat-1
		write(buff(1),'(g25.16e3)') presentedwtitem(item)
		write(buff(2),'(g25.16e3)') obswtitem(counter)
		write(buff(3),'(g25.16e3)') stdobswtitem(counter)
		write(buff(4),'(g25.16e3)') fitwtitem(counter)
		write(buff(5),'(g25.16e3)') obspwtitem(counter)
		write(buff(6),'(g25.16e3)') stdobspwtitem(counter)
		write(buff(7),'(g25.16e3)') fitpwtitem(counter)
		if(resid)then
			write(buff(8),'(g25.16e3)') residwtitem(counter)
			write(buff(9),'(g25.16e3)') stdresidwtitem(counter)
			write(buff(10),'(g25.16e3)') residpwtitem(counter)
			write(buff(11),'(g25.16e3)') stdresidpwtitem(counter)
			write(buff(12),'(g25.16e3)') residawtitem(counter)
			write(unit=unitmarginwtitem,fmt='(27a)',iostat=io)trim(adjustl(itemname(item))),",",trim(adjustl(buff1)),&
				(",",trim(adjustl(buff(row))),row=1,12)
		else
			write(unit=unitmarginwtitem,fmt='(17a)',iostat=io)trim(adjustl(itemname(item))),",",trim(adjustl(buff1)),&
				(",",trim(adjustl(buff(row))),row=1,7)
		end if
		if(io/=0) stop "Printing of product summaries failed."
		counter=counter+1
	end do
end do
return
end subroutine printmarginwtitem
