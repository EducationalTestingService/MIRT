!	print individual residual data.
!	comment is used for identification.
!	resname gives names of residual elements.
!	unit unitres is the unit for output.
!	useid indicates if ids are put out.
!	covres is residual covariance matrix.
!   fit gives fitted values.
!   observed gives observed values.
!	res  gives residual.
!   sdres gives estimated residual standard errors.
!   zres gives normalized residual.

subroutine resoutput(comment,resname,id,unitres,useid,covres,fit,observed,res,sdres,zres)
implicit none
character(len=*),intent(in)::comment
character(len=*),intent(in)::resname(:)
character(len=*),intent(in),optional::id(:)
integer,intent(in)::unitres
logical,intent(in)::useid
real(kind=8),intent(in)::covres(:,:,:),fit(:,:),observed(:,:),res(:,:),sdres(:,:),zres(:,:)

!	Buffers for writing.
character(len=4)::buffer
character(len=32)::bufferid
character(len=25),allocatable::bufff(:),buffo(:),buffr(:),buffs(:),buffz(:),buffc(:,:)
!	format
character(len=7)::writefmt
!	col is column number.
!	io is error flag.
!   nw is number of weights.
!	obs is observation number.
!	row is row number.

integer::col,io,nw,obs,row
nw=size(resname)


allocate(bufff(nw),buffo(nw),buffr(nw),buffs(nw),buffz(nw),buffc(nw,nw),stat=io)
if(io/=0) stop "Allocation failure for arrays for output of individual residuals."
write(unitres,'(2a)') "Residuals ",comment
write(buffer,'(i4)') 1+2*nw*(nw+5)
writefmt='('//trim(adjustl(buffer))//'a)'
write(unit=unitres,fmt=writefmt,iostat=io) "ID",&
    (",","Observed_"//trim(adjustl(resname(row))),row=1,nw),&
    (",","Fitted_"//trim(adjustl(resname(row))),row=1,nw),&
    (",","Residual_"//trim(adjustl(resname(row))),row=1,nw),&
    (",","SE_Res_"//trim(adjustl(resname(row))),row=1,nw),&
    (",","Adj_res__"//trim(adjustl(resname(row))),row=1,nw),&
    ((",","Cov_Res_"//trim(adjustl(resname(row)))//"_"//trim(adjustl(resname(col))),&
    col=1,nw),row=1,nw)
if(io/=0) stop "Writing of individual residuals failed."

do obs=1,size(res,2)
    do row=1,nw
        write(bufff(row),'(g25.16e3)')fit(row,obs)
        write(buffo(row),'(g25.16e3)')observed(row,obs)
        write(buffr(row),'(g25.16e3)')res(row,obs)
        write(buffs(row),'(g25.16e3)')sdres(row,obs)
        write(buffz(row),'(g25.16e3)')zres(row,obs)
        do col=1,nw
            write(buffc(row,col),'(g25.16e3)')covres(row,col,obs)
        end do
    end do
    if(useid)then
        bufferid=id(obs)
    else
        write(bufferid,'(i12)') obs
    end if
    write(unit=unitres,fmt=writefmt,iostat=io) trim(adjustl(bufferid)),&
        (",",trim(adjustl(buffo(row))),row=1,nw),&
        (",",trim(adjustl(bufff(row))),row=1,nw),&
        (",",trim(adjustl(buffr(row))),row=1,nw),&
        (",",trim(adjustl(buffs(row))),row=1,nw),&
        (",",trim(adjustl(buffz(row))),row=1,nw),&
        ((",",trim(adjustl(buffc(row,col))),col=1,nw),&
        row=1,nw)

        if(io/=0) stop "Writing of individual residuals failed."
end do
return



end subroutine resoutput
