!   Print observed scale summaries.
!   scalenames gives names of scale scores.
!   unitobscale is unit for output.
!   printobscaleres is .true. if residual information is to be given.
!   aresid gives adjusted residuas.
!   covobsave gives the covariance matrix of the observed averages.
!   covobssum gives the covariance matrix of the observed sums.
!   covresave gives the covariance matrix of the average residuals.
!   covressum gives the covariance matrix of the residual sums.
!   fitave gives the average fitted values.
!   fitsum gives the sum of fitted values.
!   obsave gives the average observed scores.
!   obssum gives the sums of observed scores.
!   resave gives the average residual.
!   ressum gives the sum of residuals.
!   sdobsave gives the standard deviation of the average observed score.
!   sdobssum gives the standard deviation of the sum of observed scores.
!   sdresave gives the standard deviation of the average residual.
!   sdressum gives the standard deviation of the sum of residuals.
!   totsum is the total weighted sample sizee for observations.

subroutine obsscaleoutput(scalename,unitobsscale,printobsscaleres,&
    aresid,covobsave,covobssum,covresave,covressum,&
    fitave,fitsum,obsave,obssum,resave,&
    ressum,sdobsave,sdobssum,sdresave,sdressum,totsum)
implicit none
character(len=32),intent(in)::scalename(:)
integer,intent(in)::unitobsscale
logical,intent(in)::printobsscaleres
real(kind=8),intent(in)::aresid(:),covobsave(:,:),covobssum(:,:),&
    covresave(:,:),covressum(:,:),&
    fitave(:),fitsum(:),&
    obsave(:),obssum(:),&
    resave(:),ressum(:),&
    sdobsave(:),sdobssum(:),sdresave(:),sdressum(:),totsum
!	buff gives buffers.

character(len=25)::buff(12)

!   col counts columns.
!	io is the flag for output error.

!	row counts rows.
integer::col,io,row
!	var is a variance.

write(unit=unitobsscale,fmt='(a)',iostat=io) 'Summary for scale scores'
if(io/=0) stop "Printing of scale score summary failed."
if(printobsscaleres)then
    write(unit=unitobsscale,fmt='(25a)',iostat=io) "Scale name",",",&
    "Weighted_obs_count",",","Obs_sum",",","Std_err_obs_sum",",",&
    "Obs_ave",",","Std_err_obs_ave",",","Fit_sum",",",&
    "Fit_ave",",",&
    "Res_sum",",","Std_err_res_sum",",",&
    "Res_ave",",","Std_err_res_ave",",",&
    "Adj_res"
else
    write(unit=unitobsscale,fmt='(11a)',iostat=io) "Scale name",",",&
        "Weighted_obs_count",",","Obs_sum",",","Std_err_obs_sum",",",&
        "Obs_ave",",","Std_err_obs_ave"
end if

do row=1,size(scalename)
    write(buff(1),'(g25.16e3)') totsum
    write(buff(2),'(g25.16e3)') obssum(row)
    write(buff(3),'(g25.16e3)') sdobssum(row)
    write(buff(4),'(g25.16e3)') obsave(row)
    write(buff(5),'(g25.16e3)') sdobsave(row)
    if(printobsscaleres)then
        write(buff(6),'(g25.16e3)') fitsum(row)
        write(buff(7),'(g25.16e3)') fitave(row)
        write(buff(8),'(g25.16e3)') ressum(row)
        write(buff(9),'(g25.16e3)') sdressum(row)
        write(buff(10),'(g25.16e3)') resave(row)
        write(buff(11),'(g25.16e3)') sdresave(row)
        write(buff(12),'(g25.16e3)') aresid(row)
    end if
    if(printobsscaleres)then
        write(unit=unitobsscale,fmt='(25a)',iostat=io)trim(adjustl(scalename(row))),&
            (",",trim(adjustl(buff(col))),col=1,12)
    else
        write(unit=unitobsscale,fmt='(11a)',iostat=io)trim(adjustl(scalename(row))),&
            (",",trim(adjustl(buff(col))),col=1,5)

    end if
    if(io/=0) stop "Printing of summary of scale scores failed."

end do

return

end subroutine obsscaleoutput

!