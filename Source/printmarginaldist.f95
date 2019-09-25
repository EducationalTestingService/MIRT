!
!	Print marginal distribution.
!	itemname contains item names.
!	numcatobs counts observed categories.
!	unitmargin provides the unit for output.
!	resid indicates if residuals are displayed.
!	fitmarg is the fitted marginal total.
!	fitpmarg is the  fitted marginal proportion.
!	obsmarg is the observed marginal total.
!	obspmarg is the observed marginal proportion.
!	presented counts weighted items presented.
!	relmarg is the reliability of the item category indicator.
!	residamarg is the adjusted residual.
!	residmarg is the residual for the marginal total.
!	residpmarg is the residual for the marginal proportion.
!	seobsmarg is the standard error for the observed indicator.
!	stdobsmarg is the standard deviation of the observed marginal.
!	stdobspmarg is the asymptotic standard deviation of the observed 
!		marginal fraction.
!	stdresidmarg is the asymptotic standard error for the residual frequency.
!	stdresidmarg is the asymptotic standard error for the residual proportion.
!	stobsmarg is the standard error of the indicator.

subroutine printmarginaldist(itemname,numcatobs,unitmargin,resid,&
	fitmarg,fitpmarg,obsmarg,obspmarg,presented,&
	relmarg,residamarg,residmarg,residpmarg,&
	seobsmarg,stdobsmarg,stdobspmarg,stdresidmarg,stdresidpmarg,stobsmarg)
implicit none
character(len=32),intent(in)::itemname(:)
integer,intent(in)::numcatobs(:),unitmargin
logical,intent(in)::resid
real(kind=8),intent(in)::fitmarg(:),fitpmarg(:),obsmarg(:),&
		obspmarg(:),presented(:),relmarg(:),residamarg(:),residmarg(:),&
		residpmarg(:),seobsmarg(:),stdobsmarg(:),stdobspmarg(:),stdresidmarg(:),&
		stdresidpmarg(:),stobsmarg(:)
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
write(unit=unitmargin,fmt='(a)',iostat=io) 'Marginal distribution of items'
if(io/=0) stop "Printing of marginal distribution failed."
if(resid)then
	write(unit=unitmargin,fmt='(33a)',iostat=io) "Item_name",",",&
		"Cat_no",",","Item_freq",",","Obs_cat_freq",",","Std_err_cat_freq",",",&
		"Fit_cat_freq",",","Obs_cat_prop",",","Std_err_cat_prop",",",&
		"Fit_cat_prop",",","Std_dev_cat_ind",",",&
		"Std_err_cat_ind",",","Rel_cat_ind",",",&
		"Res_cat_freq",",","Std_err_res_cat_freq",",",&
		"Res_cat_prop",",","Std_err_res_cat_prop",",",&
		"Adjusted_residual"
else
	write(unit=unitmargin,fmt='(23a)',iostat=io) "Item_name",",",&
		"Cat_no",",","Item_freq",",","Obs_cat_freq",",","Std_err_cat_freq",",",&
		"Fit_cat_freq",",","Obs_cat_prop",",","Std_err_cat_prop",",",&
		"Fit_cat_prop",",","Std_dev_cat_ind",",",&
		"Std_err_cat_ind",",","Rel_cat_ind"
end if

do item=1,size(itemname)
	do cat=1,numcatobs(item)
		write(buff1,'(i12)') cat-1
		write(buff(1),'(g25.16e3)') presented(item)
		write(buff(2),'(g25.16e3)') obsmarg(counter)
		write(buff(3),'(g25.16e3)') stdobsmarg(counter)
		write(buff(4),'(g25.16e3)') fitmarg(counter)
		write(buff(5),'(g25.16e3)') obspmarg(counter)
		write(buff(6),'(g25.16e3)') stdobspmarg(counter)
		write(buff(7),'(g25.16e3)') fitpmarg(counter) 
		write(buff(8),'(g25.16e3)') stobsmarg(counter)
		write(buff(9),'(g25.16e3)') seobsmarg(counter)
		write(buff(10),'(g25.16e3)') relmarg(counter)
		if(resid)then
			write(buff(11),'(g25.16e3)') residmarg(counter)
			write(buff(12),'(g25.16e3)') stdresidmarg(counter)
			write(buff(13),'(g25.16e3)') residpmarg(counter)
			write(buff(14),'(g25.16e3)') stdresidpmarg(counter)
			write(buff(15),'(g25.16e3)') residamarg(counter)
			write(unit=unitmargin,fmt='(33a)',iostat=io)trim(adjustl(itemname(item))),",",trim(adjustl(buff1)),&
				(",",trim(adjustl(buff(row))),row=1,15)
		else
			write(unit=unitmargin,fmt='(23a)',iostat=io)trim(adjustl(itemname(item))),",",trim(adjustl(buff1)),&
				(",",trim(adjustl(buff(row))),row=1,10)				
		end if
		if(io/=0) stop "Printing of marginal distribution failed."
		counter=counter+1
	end do
end do
return
end subroutine printmarginaldist
