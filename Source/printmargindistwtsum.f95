
!
!	Print marginal distribution of weighted sum.
!	weightname contains sum name.
!	maxscore is the maximum score.
!	minscore is the minimum score.
!	unitmarginwtsum provides the unit for output.
!	resid indicates if residuals are displayed.
!	fitcummargwtsum is the fitted cumulative marginal
!		frequency distribution.
!	fitmargwtsum is the fitted marginal frequency distribution.
!	fitpcummargwtsum is the fitted cumulative marginal
!		probability distribution.
!	fitpmargwtsum is the fitted marginal probability distribution.
!	obscummargwtsum is the observed cumulative marginal
!		frequency distribution.
!	obsmargwtsum is the observed marginal frequency distribution.
!	obspcummargwtsum is the observed cumulative marginal
!		probability distribution.
!	obspmargwtsum is the observed marginal
!		probability distribution.
!	presentedwtsum counts weighted items presented
!		for the weighted sum.
!	residacummargwtsum is the adjusted residual for the cumulative
!		marginal distribution.
!	residamargwtsum is the adjusted residual for the marginal
!		distribution.
!	residcummargwtsum is the residual for the cumulative marginal
!		frequency distribution.
!	residmargwtsum is the residual for the marginal frequency
!		distribution.
!	residpcummargwtsum is the residual for the cumulative marginal
!		 probability distribution.
!	residpmargwtsum is the residual for the marginal probability
!		distribution.
!	stdobscummargwtsum is the standard deviation of the observed
!		cumulative marginal frequency distribution.
!	stdobsmargwtsum is the standard deviation of the observed
!		marginal frequency distribution.
!	stdobspcummargwtsum is the asymptotic standard deviation of
!		the observed cumulative marginal probability distribution.
!	stdobspmargwtsum is the asymptotic standard deviation of the
!		observed marginal probability distribution.
!	stdresidcummargwtsum is the asymptotic standard deviation
!		of the residual cumulative marginal frequency
!		distribution.
!	stdresidmargwtsum is the asymptotic standard deviation of the
!		residual marginal frequency distribution.
!	stdresidpcummargwtsum is the asymptotic standard deviation of
!		the residual cumulative marginal probability distribution.
!	stdresidpmargwtsum is the asymptotic standard deviation of the
!		residual marginal probability distribution.

subroutine printmargindistwtsum(weightname,maxscore,minscore,unitmarginwtsum,resid,&
	fitcummargwtsum,fitmargwtsum,fitpcummargwtsum,fitpmargwtsum,&
	obscummargwtsum,obsmargwtsum,obspcummargwtsum,obspmargwtsum,&
	presentedwtsum,&
	residacummargwtsum,residamargwtsum,&
	residcummargwtsum,residmargwtsum,residpcummargwtsum,&
	residpmargwtsum,stdobscummargwtsum,stdobsmargwtsum,&
	stdobspcummargwtsum,stdobspmargwtsum,&
	stdresidcummargwtsum,stdresidmargwtsum,&
	stdresidpcummargwtsum,stdresidpmargwtsum)

implicit none
character(len=32)::weightname
integer,intent(in)::maxscore,minscore,unitmarginwtsum
logical,intent(in)::resid
real(kind=8),intent(in)::fitcummargwtsum(minscore:maxscore),fitmargwtsum(minscore:maxscore),&
	fitpcummargwtsum(minscore:maxscore),fitpmargwtsum(minscore:maxscore),&
	obscummargwtsum(minscore:maxscore),obsmargwtsum(minscore:maxscore),&
	obspcummargwtsum(minscore:maxscore),obspmargwtsum(minscore:maxscore),presentedwtsum,&
	residacummargwtsum(minscore:maxscore),residamargwtsum(minscore:maxscore),&
	residcummargwtsum(minscore:maxscore),residmargwtsum(minscore:maxscore),&
	residpcummargwtsum(minscore:maxscore),residpmargwtsum(minscore:maxscore),&
	stdobscummargwtsum(minscore:maxscore),stdobsmargwtsum(minscore:maxscore),&
	stdobspcummargwtsum(minscore:maxscore),stdobspmargwtsum(minscore:maxscore),&
	stdresidcummargwtsum(minscore:maxscore),stdresidmargwtsum(minscore:maxscore),&
	stdresidpcummargwtsum(minscore:maxscore),stdresidpmargwtsum(minscore:maxscore)
!	buff0 and buff are buffers.
character(len=12)::buff0
character(len=25)::buff(23)


!	io is the flag for output error.
!	row counts entries.
!	score is a score.
integer::io,row,score
write(unit=unitmarginwtsum,fmt='(2a)',iostat=io) 'Marginal distribution of weighted sum ',trim(adjustl(weightname))
if(io/=0) stop "Printing failed for marginal distribution of weighted sum."
if(resid)then
	write(unit=unitmarginwtsum,fmt='(47a)',iostat=io) "Score",",",&
		"Tot_scores",",","Obs_score_freq",",","Std_err_obs_score_freq",",",&
		"Cum_obs_score_freq",",","Std_err_cum_obs_score_freq",",",&
		"Fit_score_freq",",","Cum_fit_score_freq",",",&
		"Obs_score_prop",",","Std_err_obs_score_prop",",",&
		"Cum_obs_score_prop",",","Std_err_cum_obs_score_prop",",",&
		"Fit_score_prop",",",&
		"Cum_fit_score_prop",",",&
		"Res_score_freq",",","Std_err_res_score_freq",",",&
		"Res_cum_score_freq",",","Std_err_res_cum_score_freq",",",&
		"Res_score_prop",",","Std_err_res_score_prop",",",&
		"Res_cum_score_prop",",","Std_err_res_cum_score_prop",",",&
		"Adj_score_res",",","Adj_cum_score_res"
else
	write(unit=unitmarginwtsum,fmt='(27a)',iostat=io) "Score",",",&
		"Tot_scores",",","Obs_score_freq",",","Std_err_obs_score_freq",",",&
		"Cum_obs_score_freq",",","Std_err_cum_obs_score_freq",",",&
		"Fit_score_freq",",","Cum_fit_score_freq",",",&
		"Obs_score_prop",",","Std_err_obs_score_prop",",",&
		"Cum_obs_score_prop",",","Std_err_cum_obs_score_prop",",",&
		"Fit_score_prop",",",&
		"Cum_fit_score_prop"
		
end if
do score=lbound(obsmargwtsum,1),ubound(obsmargwtsum,1)
	write(buff0,'(i12)') score
	write(buff(1),'(g25.16e3)') presentedwtsum
	write(buff(2),'(g25.16e3)') obsmargwtsum(score)
	write(buff(3),'(g25.16e3)') stdobsmargwtsum(score)
	write(buff(4),'(g25.16e3)') obscummargwtsum(score)
	write(buff(5),'(g25.16e3)') stdobscummargwtsum(score)
	write(buff(6),'(g25.16e3)') fitmargwtsum(score)
	write(buff(7),'(g25.16e3)') fitcummargwtsum(score)
	write(buff(8),'(g25.16e3)') obspmargwtsum(score)
	write(buff(9),'(g25.16e3)') stdobspmargwtsum(score)
	write(buff(10),'(g25.16e3)') obspcummargwtsum(score)
	write(buff(11),'(g25.16e3)') stdobspcummargwtsum(score)
	write(buff(12),'(g25.16e3)') fitpmargwtsum(score)
	write(buff(13),'(g25.16e3)') fitpcummargwtsum(score) 
	if(resid)then
		write(buff(14),'(g25.16e3)') residmargwtsum(score)
		write(buff(15),'(g25.16e3)') stdresidmargwtsum(score)
		write(buff(16),'(g25.16e3)') residcummargwtsum(score)
		write(buff(17),'(g25.16e3)') stdresidcummargwtsum(score)
		write(buff(18),'(g25.16e3)') residpmargwtsum(score)
		write(buff(19),'(g25.16e3)') stdresidpmargwtsum(score)
		write(buff(20),'(g25.16e3)') residpcummargwtsum(score)
		write(buff(21),'(g25.16e3)') stdresidpcummargwtsum(score)
		
		write(buff(22),'(g25.16e3)') residamargwtsum(score)
		write(buff(23),'(g25.16e3)') residacummargwtsum(score)
		write(unit=unitmarginwtsum,fmt='(47a)',iostat=io)trim(adjustl(buff0)),&
				(",",trim(adjustl(buff(row))),row=1,23)
	else
		write(unit=unitmarginwtsum,fmt='(27a)',iostat=io)trim(adjustl(buff0)),&
			(",",trim(adjustl(buff(row))),row=1,13)				
	end if
	if(io/=0) stop "Printing failed for marginal distribution of weighted sum."
end do
return
end subroutine printmargindistwtsum

