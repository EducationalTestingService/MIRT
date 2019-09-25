!
!	Print two-way marginal distribution.
!	itemname contains item names.
!	numcatobs counts observed categories.
!	unitmargin2 provides the unit for output.
!	resid indicates if residuals are displayed.
!	fitmarg2 is the fitted two-way marginal total.
!	fitpmarg2 is the fitted two-way marginal proportion.
!	obsmarg2 is the observed two-way marginal total.
!	obspmarg2 is the observed two-way marginal proportion.
!	presented2 counts weighted pairs of items presented.
!	residamarg2 is the adjusted residual.
!	residmarg2 is the residual for the two-way marginal total.
!	residpmarg2 is the residual for the two-way marginal proportion.
!	stdobsmarg2 is the standard deviation of the observed two-way marginal.
!	stdobspmarg2 is the asymptotic standard deviation of the observed 
!		two-way marginal fraction.
!	stdresidmarg2 is the asymptotic standard error for the residual frequency.
!	stdresidpmarg2 is the asymptotic standard error for the residual proportion.


subroutine printmarginal2(itemname,numcatobs,unitmargin2,resid,&
	fitmarg2,fitpmarg2,obsmarg2,obspmarg2,presented2,&
	residamarg2,residmarg2,residpmarg2,&
	stdobsmarg2,stdobspmarg2,stdresidmarg2,stdresidpmarg2)
implicit none
character(len=32),intent(in)::itemname(:)
integer,intent(in)::numcatobs(:),unitmargin2
logical,intent(in)::resid
real(kind=8),intent(in)::fitmarg2(:,:),fitpmarg2(:,:),obsmarg2(:,:),&
		obspmarg2(:,:),presented2(:,:),residamarg2(:,:),residmarg2(:,:),&
		residpmarg2(:,:),stdobsmarg2(:,:),stdobspmarg2(:,:),&
		stdresidmarg2(:,:),stdresidpmarg2(:,:)
!	buff0, buff1, and buff are buffers.
character(len=12)::buff0,buff1
character(len=25)::buff(12)

!	cat counts categories, as does cat1.
!	counter counts positions, as does counter1.
!	io is the flag for output error.
!	item counts items, as does item1.
!	row counts rows.
integer::cat,cat1,counter,counter1,io,item,item1,row
!	var is a variance.
real(kind=8)::var
counter1=numcatobs(1)+1
write(unit=unitmargin2,fmt='(a)',iostat=io) 'Two-way marginal distribution of items'
if(io/=0) stop "Printing of two-way marginal distribution failed."
if(resid)then
	write(unit=unitmargin2,fmt='(31a)',iostat=io) "Item_name",",",&
		"Cat_no",",","Item_name",",",&
		"Cat_no",",","Item_pair_freq",",","Obs_cat_pair_freq",",","Std_err_obs_cat_pair_freq",",",&
		"Fit_cat_pair_freq",",","Obs_cat_pair_prop",",","Std_err_cat_pair_prop",",",&
		"Fit_cat_pair_prop",",",&
		"Res_cat_pair_freq",",","Std_err_res_cat_pair_freq",",",&
		"Res_cat_pair_prop",",","Std_err_res_cat_pair_prop",",",&
		"Adjusted_residual"
else
	write(unit=unitmargin2,fmt='(21a)',iostat=io) "Item_name",",",&
		"Cat_no",",","Item_name",",",&
		"Cat_no",",","Item_pair_freq",",","Obs_cat_pair_freq",",","Std_err_obs_cat_pair_freq",",",&
		"Fit_cat_pair_freq",",","Obs_cat_pair_prop",",","Std_err_cat_pair_prop",",",&
		"Fit_cat_pair_prop"
		
end if
do item1=2,size(itemname)
do cat1=1,numcatobs(item1)
write(buff0,'(i12)') cat1-1
counter=1
do item=1,item1-1
	
	do cat=1,numcatobs(item)
		write(buff1,'(i12)') cat-1
		write(buff(1),'(g25.16e3)') presented2(item,item1)
		write(buff(2),'(g25.16e3)') obsmarg2(counter,counter1)
		write(buff(3),'(g25.16e3)') stdobsmarg2(counter,counter1)
		write(buff(4),'(g25.16e3)') fitmarg2(counter,counter1)
		write(buff(5),'(g25.16e3)') obspmarg2(counter,counter1)
		write(buff(6),'(g25.16e3)') stdobspmarg2(counter,counter1)
		write(buff(7),'(g25.16e3)') fitpmarg2(counter,counter1) 
		
		if(resid)then
			write(buff(8),'(g25.16e3)') residmarg2(counter,counter1)
			write(buff(9),'(g25.16e3)') stdresidmarg2(counter,counter1)
			write(buff(10),'(g25.16e3)') residpmarg2(counter,counter1)
			write(buff(11),'(g25.16e3)') stdresidpmarg2(counter,counter1)
			write(buff(12),'(g25.16e3)') residamarg2(counter,counter1)
			write(unit=unitmargin2,fmt='(31a)',iostat=io)trim(adjustl(itemname(item1))),",",trim(adjustl(buff0)),",",&
				trim(adjustl(itemname(item))),",",trim(adjustl(buff1)),&
				(",",trim(adjustl(buff(row))),row=1,12)
		else
			write(unit=unitmargin2,fmt='(21a)',iostat=io)trim(adjustl(itemname(item1))),",",trim(adjustl(buff0)),",",&
				trim(adjustl(itemname(item))),",",trim(adjustl(buff1)),&
				(",",trim(adjustl(buff(row))),row=1,7)				
		end if
		if(io/=0) stop "Printing of marginal distribution failed."
		counter=counter+1
	end do
	
end do
counter1=counter1+1
end do
end do
return
end subroutine printmarginal2

