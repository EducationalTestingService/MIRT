!	Print totals and averages of products of predictors and indicators of item categories.
!	itemname contains item names.
!	predname contains predictor names.
!	numcatobs counts observed categories.
!	unitpreditem provides the unit for output.
!	resid indicates if residuals are displayed.
!	fitppred is the  fitted average of products.
!	fitpred is the fitted total of products.
!	obsppred is the observed average of products.
!	obspred is the observed total of products.
!	presented counts weighted items presented.
!	residapred is the adjusted residual.
!	residppred is the residual for the average of products.
!	residpred is the residual for the total of products.
!	stdobsppred is the standard error of the observed average of products.
!	stdobspred is the standard error of the observed 
!		total of products.
!	stdresidppred is the asymptotic standard error for the residual for average
!       of products.
!	stdresidpred is the asymptotic standard error for the residual for total
!       of products.


subroutine printmarginpred(itemname,predname,numcatobs,unitpreditem,&
    resid,fitppred,fitpred,obsppred,obspred,presented,&
    residapred,residppred,residpred,&
    stdobsppred,stdobspred,stdresidppred,stdresidpred)
implicit none
character(len=32),intent(in)::itemname(:)
character(len=32),intent(in)::predname(:)
integer,intent(in)::numcatobs(:),unitpreditem
logical,intent(in)::resid
real(kind=8),intent(in)::fitppred(:,:),fitpred(:,:),&
    obsppred(:,:),obspred(:,:),presented(:),&
    residapred(:,:),residppred(:,:),residpred(:,:),&
    stdobsppred(:,:),stdobspred(:,:),stdresidppred(:,:),stdresidpred(:,:)
!	buff1 and buff are buffers.
character(len=12)::buff1
character(len=25)::buff(15)

!	cat counts categories.
!	counter counts positions.
!	io is the flag for output error.
!	item counts items.
!	pred counts predictors.
!	row counts rows.
integer::cat,counter,io,item,pred,row
!	var is a variance.
real(kind=8)::var

do pred=2,size(predname)

counter=1
write(unit=unitpreditem,fmt='(2a)',iostat=io) 'Totals and averages of products for ',trim(adjustl(predname(pred)))

if(io/=0) stop "Printing of totals and averages of products failed."
if(resid)then
	write(unit=unitpreditem,fmt='(27a)',iostat=io) "Item_name",",",&
		"Cat_no",",","Item_freq",",","Obs_tot_prod",",","Std_err_obs_tot_prod",",",&
		"Fit_tot_prod",",","Obs_ave_prod",",","Std_err_obs_ave_prod",",",&
		"Fit_ave_prod",",",&
		"Res_tot_prod",",","Std_err_res_tot_prod",",",&
		"Res_ave_prod",",","Std_err_res_ave_prod",",",&
		"Adjusted_residual"
else
	write(unit=unitpreditem,fmt='(17a)',iostat=io) "Item_name",",",&
        "Cat_no",",","Item_freq",",","Obs_tot_prod",",",&
        "Std_err_obs_tot_prod",",",&
        "Fit_tot_prod",",","Obs_ave_prod",",","Std_err_obs_ave_prod",",",&
        "Fit_ave_prod"
end if

do item=1,size(itemname)
	do cat=1,numcatobs(item)
		write(buff1,'(i12)') cat-1
		write(buff(1),'(g25.16e3)') presented(item)
		write(buff(2),'(g25.16e3)') obspred(counter,pred)
		write(buff(3),'(g25.16e3)') stdobspred(counter,pred)
		write(buff(4),'(g25.16e3)') fitpred(counter,pred)
		write(buff(5),'(g25.16e3)') obsppred(counter,pred)
		write(buff(6),'(g25.16e3)') stdobsppred(counter,pred)
		write(buff(7),'(g25.16e3)') fitppred(counter,pred)
		if(resid)then
			write(buff(8),'(g25.16e3)') residpred(counter,pred)
			write(buff(9),'(g25.16e3)') stdresidpred(counter,pred)
			write(buff(10),'(g25.16e3)') residppred(counter,pred)
			write(buff(11),'(g25.16e3)') stdresidppred(counter,pred)
			write(buff(12),'(g25.16e3)') residapred(counter,pred)
			write(unit=unitpreditem,fmt='(27a)',iostat=io)trim(adjustl(itemname(item))),",",trim(adjustl(buff1)),&
				(",",trim(adjustl(buff(row))),row=1,12)
		else
			write(unit=unitpreditem,fmt='(17a)',iostat=io)trim(adjustl(itemname(item))),",",trim(adjustl(buff1)),&
				(",",trim(adjustl(buff(row))),row=1,7)
		end if
		if(io/=0) stop "Printing of product summaries failed."
		counter=counter+1
	end do
end do
end do
return
end subroutine printmarginpred
