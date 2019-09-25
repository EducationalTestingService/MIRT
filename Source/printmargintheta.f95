!	Print totals and averages of products of transformed latent variables and indicators of item categories.
!	itemname contains item names.
!	skillname contains skill\r names.
!	numcatobs counts observed categories.
!	unitthetaitem provides the unit for output.
!	resid indicates if residuals are displayed.
!	fitptheta is the  fitted average of products.
!	fittheta is the fitted total of products.
!	obsptheta is the observed average of products.
!	obstheta is the observed total of products.
!	presentedtheta counts weighted items presented.
!	residatheta is the adjusted residual.
!	residptheta is the residual for the average of products.
!	residtheta is the residual for the total of products.
!	stdobsptheta is the standard error of the observed average of products.
!	stdobstheta is the standard error of the observed 
!		total of products.
!	stdresidptheta is the asymptotic standard error for the residual for average
!       of products.
!	stdresidtheta is the asymptotic standard error for the residual for total
!       of products.


subroutine printmargintheta(itemname,skillname,numcatobs,unitthetaitem,&
    resid,fitptheta,fittheta,obsptheta,obstheta,presentedtheta,&
    residatheta,residptheta,residtheta,&
    stdobsptheta,stdobstheta,stdresidptheta,stdresidtheta)
implicit none
character(len=*),intent(in)::itemname(:)
character(len=*),intent(in)::skillname(:)
integer,intent(in)::numcatobs(:),unitthetaitem
logical,intent(in)::resid
real(kind=8),intent(in)::fitptheta(:,:),fittheta(:,:),&
    obsptheta(:,:),obstheta(:,:),presentedtheta(:),&
    residatheta(:,:),residptheta(:,:),residtheta(:,:),&
    stdobsptheta(:,:),stdobstheta(:,:),stdresidptheta(:,:),stdresidtheta(:,:)
!	buff1 and buff are buffers.
character(len=12)::buff1
character(len=25)::buff(15)

!	cat counts categories.
!	counter counts positions.
!	io is the flag for output error.
!	item counts items.
!	skill counts skills.
!	row counts rows.
integer::cat,counter,io,item,skill,row
!	var is a variance.
real(kind=8)::var
counter=1
do skill=1,size(skillname)
write(unit=unitthetaitem,fmt='(2a)',iostat=io) 'Totals and averages of products for ',skillname(skill)

if(io/=0) stop "Printing of totals and averages of products failed."
if(resid)then
	write(unit=unitthetaitem,fmt='(27a)',iostat=io) "Item_name",",",&
		"Cat_no",",","Item_freq",",","Obs_tot_prod",",","Std_err_obs_tot_prod",",",&
		"Fit_tot_prod",",","Obs_ave_prod",",","Std_err_obs_ave_prod",",",&
		"Fit_ave_prod",",",&
		"Res_tot_prod",",","Std_err_res_tot_prod",",",&
		"Res_ave_prod",",","Std_err_res_ave_prod",",",&
		"Adjusted_residual"
else
	write(unit=unitthetaitem,fmt='(17a)',iostat=io) "Item_name",",",&
        "Cat_no",",","Item_freq",",","Obs_tot_prod",",",&
        "Std_err_obs_tot_prod",",",&
        "Fit_tot_prod",",","Obs_ave_prod",",","Std_err_obs_ave_prod",",",&
        "Fit_ave_prod"
end if

do item=1,size(itemname)
	do cat=1,numcatobs(item)
		write(buff1,'(i12)') cat-1
		write(buff(1),'(g25.16e3)') presentedtheta(item)
		write(buff(2),'(g25.16e3)') obstheta(counter,skill)
		write(buff(3),'(g25.16e3)') stdobstheta(counter,skill)
		write(buff(4),'(g25.16e3)') fittheta(counter,skill)
		write(buff(5),'(g25.16e3)') obsptheta(counter,skill)
		write(buff(6),'(g25.16e3)') stdobsptheta(counter,skill)
		write(buff(7),'(g25.16e3)') fitptheta(counter,skill)
		if(resid)then
			write(buff(8),'(g25.16e3)') residtheta(counter,skill)
			write(buff(9),'(g25.16e3)') stdresidtheta(counter,skill)
			write(buff(10),'(g25.16e3)') residptheta(counter,skill)
			write(buff(11),'(g25.16e3)') stdresidptheta(counter,skill)
			write(buff(12),'(g25.16e3)') residatheta(counter,skill)
			write(unit=unitthetaitem,fmt='(27a)',iostat=io)trim(adjustl(itemname(item))),",",trim(adjustl(buff1)),&
				(",",trim(adjustl(buff(row))),row=1,12)
		else
			write(unit=unitthetaitem,fmt='(17a)',iostat=io)trim(adjustl(itemname(item))),",",trim(adjustl(buff1)),&
				(",",trim(adjustl(buff(row))),row=1,7)
		end if
		if(io/=0) stop "Printing of product summaries failed."
		counter=counter+1
	end do
end do
end do
return
end subroutine printmargintheta
