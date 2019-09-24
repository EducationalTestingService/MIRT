subroutine gettitle(title)
implicit none
integer::io
character(len=80),intent(out):: title
namelist/runtitle/title
title="GPCM/2PL Model with Normal Latent Distribution"
read(*,nml=runtitle,iostat=io)
if(io/=0)stop " runtitle namelist not read successfully"
end subroutine gettitle
