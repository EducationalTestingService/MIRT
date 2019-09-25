!
!	Print two-way cross products of item scores.
!	itemname contains item names.

!	unitmargins2 provides the unit for output.
!	resid indicates if residuals are displayed.
!	fitmargs2 is the fitted total cross products.
!	fitpmargs2 is the fitted average cross products.
!	obsmargs2 is the observed total cross products.
!	obspmargs2 is the observed average cross products.
!	presenteds2 counts weighted pairs of items presented.
!	residamargs2 is the adjusted residual.
!	residmargs2 is the residual for the total cross products.
!	residpmargs2 is the residual for the average cross products.
!	stdobsmargs2 is the standard deviation of the observed total cross products.
!	stdobspmargs2 is the asymptotic standard deviation of the observed 
!		average cross products.
!	stdresidmargs2 is the asymptotic standard error for the residual total cross products.
!	stdresidpmargs2 is the asymptotic standard error for the residual average cross products.


subroutine printmargins2(itemname,unitmargins2,resid,&
	fitmargs2,fitpmargs2,obsmargs2,obspmargs2,presenteds2,&
	residamargs2,residmargs2,residpmargs2,&
	stdobsmargs2,stdobspmargs2,stdresidmargs2,stdresidpmargs2)
implicit none
character(len=*),intent(in)::itemname(:)
integer,intent(in)::unitmargins2
logical,intent(in)::resid
real(kind=8),intent(in)::fitmargs2(:,:),fitpmargs2(:,:),obsmargs2(:,:),&
		obspmargs2(:,:),presenteds2(:,:),residamargs2(:,:),residmargs2(:,:),&
		residpmargs2(:,:),stdobsmargs2(:,:),stdobspmargs2(:,:),&
		stdresidmargs2(:,:),stdresidpmargs2(:,:)
!	buff is a  buffer.

character(len=25)::buff(12)


!	io is the flag for output error.
!	item counts items, as does item1.
!	row counts rows.
integer::io,item,item1,row
!	var is a variance.
real(kind=8)::var

write(unit=unitmargins2,fmt='(a)',iostat=io) 'Cross products of item scores'
if(io/=0) stop "Printing of cross products of item scores failed."
if(resid)then
	write(unit=unitmargins2,fmt='(27a)',iostat=io) "Item_name",",","Item_name",",",&
		"Item_pair_freq",",","Obs_tot_cross_prod",",","Std_err_obs_tot_cross_prod",",",&
		"Fit_tot_cross_prod",",","Obs_ave_cross_prod",",","Std_err_ave_cross_prod",",",&
		"Fit_ave_cross_prod",",",&
		"Res_tot_cross_prod",",","Std_err_res_tot_cross_prod",",",&
		"Res_ave_cross_prod",",","Std_err_res_ave_cross_prod",",",&
		"Adjusted_residual"
else
	write(unit=unitmargins2,fmt='(17a)',iostat=io) "Item_name",",",&
		"Item_name",",",&
		"Item_pair_freq",",","Obs_tot_cross_prod",",","Std_err_obs_tot_cross_prod",",",&
		"Fit_tot_cross_prod",",","Obs_ave_cross_prod",",","Std_err_ave_cross_prod",",",&
		"Fit_ave_cross_prod"
		
end if
do item1=2,size(itemname)
	
	do item=1,item1-1
	
		
		
		write(buff(1),'(g25.16e3)') presenteds2(item,item1)
		write(buff(2),'(g25.16e3)') obsmargs2(item,item1)
		write(buff(3),'(g25.16e3)') stdobsmargs2(item,item1)
		write(buff(4),'(g25.16e3)') fitmargs2(item,item1)
		write(buff(5),'(g25.16e3)') obspmargs2(item,item1)
		write(buff(6),'(g25.16e3)') stdobspmargs2(item,item1)
		write(buff(7),'(g25.16e3)') fitpmargs2(item,item1) 
		
		if(resid)then
			write(buff(8),'(g25.16e3)') residmargs2(item,item1)
			write(buff(9),'(g25.16e3)') stdresidmargs2(item,item1)
			write(buff(10),'(g25.16e3)') residpmargs2(item,item1)
			write(buff(11),'(g25.16e3)') stdresidpmargs2(item,item1)
			write(buff(12),'(g25.16e3)') residamargs2(item,item1)
			write(unit=unitmargins2,fmt='(27a)',iostat=io)trim(adjustl(itemname(item1))),",",&
				trim(adjustl(itemname(item))),&
				(",",trim(adjustl(buff(row))),row=1,12)
		else
			write(unit=unitmargins2,fmt='(17a)',iostat=io)trim(adjustl(itemname(item1))),",",&
				trim(adjustl(itemname(item))),&
				(",",trim(adjustl(buff(row))),row=1,7)				
		end if
		if(io/=0) stop "Printing of marginal distribution failed."
		
	end do
	
end do


return
end subroutine printmargins2


