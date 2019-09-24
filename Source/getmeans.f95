!
!  See what means of latent vectors are needed.
!	npredlist is the number of predicting vectors to be used.  It is relevant if allmeans is .FALSE.
!	allmeans is .TRUE. only if a mean is desired for each observed predicting vector.

subroutine getmeans(npredlist,allmeans)
implicit none 
integer,intent(out)::npredlist
logical,intent(out)::allmeans
namelist/meanprintspecs/npredlist,allmeans
integer::io
allmeans=.false.
npredlist=0
read(*,nml=meanprintspecs,iostat=io)
if(io/=0) stop 'Specifications for printing means not successfully read.'
if (npredlist<0) npredlist=0
return
end subroutine getmeans