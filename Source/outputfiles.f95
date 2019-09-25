!	Output units.  See unitdef for descriptions.

subroutine outputfiles()
use unitdef
implicit none



character(len=256)::filename
!	count is number of units.
!	io is the error indicator for reading.
!	unitno is a unit counter.
integer::count,io,unitno

namelist/numberunits/count
namelist/files/filename
namelist/units/unitalpha,&
    unitbeta,unitbetacov,unitbetacov_complex,unitbetacov_louis,unitbetacov_sandwich,&
	uniteap,uniteapscale,uniteapskill,uniteaptrans,uniteapwt,&
	unitfreq,&
    unitgrad,unitguess,&
    unitinfo,unitinfodist,unitirf,unititeration,unititerationstart,&
    unitmargin,unitmargins2,unitmarginwtsum,unitmargin2,&
    unitmeans,unitml,unitmp,unitobservedscale,unitobsscale,unitobsscalerel,&
    unitoutdata,&
	unitparam,unitparamcov,unitparamcov_complex,unitparamcov_louis,unitparamcov_sandwich,&
    unitpost,unitpreditem,unitprob,unitprobdist,unitpwt,unitpwtdist,&
	unitrel,unitrelskill,unitrelwt,unitreswt,unitreswtdist,&
    unitscalerel,&
    unitthetaitem,unittitle,&
    unitwtitem
count=1


filename="output.csv"

read(*,nml=numberunits,iostat=io)
if(io/=0.or.count<1.or.count>90) stop "Number of output units not successfully read."

do unitno=1,count
	read(*,nml=files,iostat=io)
	open(unit=9+unitno,file=filename,iostat=io)
	if(io/=0) stop "Output file not opened successfully."
end do
read(*,nml=units,iostat=io)
if(io/=0) stop "Output units not successfully read."
if(min(unitalpha,&
unitbeta,unitbetacov,unitbetacov_complex,unitbetacov_louis,unitbetacov_sandwich,&
uniteap,uniteapscale,uniteapskill,uniteaptrans,uniteapwt,&
unitfreq,&
unitgrad,unitguess,&
unitinfo,unitinfodist,unitirf,unititeration,unititerationstart,&
unitmargin,unitmargins2,unitmarginwtsum,unitmargin2,&
unitmeans,unitml,unitmp,unitobservedscale,unitobsscale,unitobsscalerel,&
unitoutdata,&
unitparam,unitparamcov,unitparamcov_complex,unitparamcov_louis,unitparamcov_sandwich,&
unitpost,unitpreditem,unitprob,unitprobdist,unitpwt,unitpwtdist,&
unitrel,unitrelskill,unitrelwt,unitreswt,unitreswtdist,&
unitscalerel,&
unitthetaitem,unittitle,&
unitwtitem)<10.or.&
max(unitalpha,&
unitbeta,unitbetacov,unitbetacov_complex,unitbetacov_louis,unitbetacov_sandwich,&
uniteap,uniteapscale,uniteapskill,uniteaptrans,uniteapwt,&
unitfreq,&
unitgrad,unitguess,&
unitinfo,unitinfodist,unitirf,unititeration,unititerationstart,&
unitmargin,unitmargins2,unitmarginwtsum,unitmargin2,&
unitmeans,unitmp,unitobservedscale,unitobsscale,unitobsscalerel,&
unitoutdata,&
unitparam,unitparamcov,unitparamcov_complex,unitparamcov_louis,unitparamcov_sandwich,&
unitpost,unitpreditem,unitprob,unitprobdist,unitpwt,unitpwtdist,&
unitrel,unitrelskill,unitrelwt,unitreswt,unitreswtdist,&
unitscalerel,&
unitthetaitem,unittitle,&
unitwtitem)>9+count)  stop "Output units not successfully read."
return
end subroutine outputfiles

