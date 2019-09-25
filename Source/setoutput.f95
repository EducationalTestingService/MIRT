!	Set output

subroutine setoutput()
use unitdef
implicit none

!   counter is i.
!	io error.
!   unit is a unit.
integer::i,io,unit
namelist/output/&
    printalpha,&
    printbeta,printbetacov,printbetacov_louis,printbetacov_sandwich,printbetacov_complex,&
    printeap,printeapscale,printeaptrans,printeapskill,printeapwt,printent,printentdist,&
	printfreq,&
    printgrad,printguess,&
	printirf,printirfres,&
    printmargin,printmarginres,&
    printmargins2,printmargins2res,&
	printmarginwtsum,printmarginwtsumres,&
	printmargin2,printmargin2res,printmeans,printml,printmp,&
    printobservedscale,printobsscale,printobsscalerel,printobsscaleres,&
    printparam,printparamcov,printparamcov_complex,printparamcov_louis,printparamcov_sandwich,&
    printpost,printpreditem,printpreditemres,printprob,printprobdist,&
    printpwt,printpwtdist,&
    printrel,printrelskill,printreltrans,printrelwt,printreswt,printreswtdist,&
    printscalerel,&
    printthetaitem,printthetaitemres,&
    printwtitem,printwtitemres
!	Default settings.



read(*,nml=output,iostat=io)
if(io/=0) stop "Output options not read successfully."
if(printirfres)printirf=.true.
if(printmarginres)printmargin=.true.
if(printmargins2res) printmargins2=.true.
if(printmarginwtsumres)printmarginwtsum=.true.
if(printmargin2res)printmargin2=.true.
if(printobsscaleres)printobsscale=.true.
if(printpreditemres) printpreditem=.true.
if(printthetaitemres)printthetaitem=.true.
if(printwtitemres)printwtitem=.true.


if(printalpha)unitcount(unitalpha)=unitcount(unitalpha)+1
if(printbeta)unitcount(unitbeta)=unitcount(unitbeta)+1
if(printbetacov)unitcount(unitbetacov)=unitcount(unitbetacov)+1
if(printbetacov_louis)unitcount(unitbetacov_louis)=unitcount(unitbetacov_louis)+1
if(printbetacov_sandwich)unitcount(unitbetacov_sandwich)=unitcount(unitbetacov_sandwich)+1
if(printbetacov_complex)unitcount(unitbetacov_complex)=unitcount(unitbetacov_complex)+1
if(printeap)unitcount(uniteap)=unitcount(uniteap)+1
if(printeapscale)unitcount(uniteapscale)=unitcount(uniteapscale)+1
if(printeapskill)unitcount(uniteapskill)=unitcount(uniteapskill)+1
if(printeaptrans)unitcount(uniteaptrans)=unitcount(uniteaptrans)+1
if(printeapwt)unitcount(uniteapwt)=unitcount(uniteapwt)+1
if(printent)unitcount(unitinfo)=unitcount(unitinfo)+1
if(printentdist)unitcount(unitinfodist)=unitcount(unitinfodist)+1
if(printfreq)unitcount(unitfreq)=unitcount(unitfreq)+1
if(printgrad)unitcount(unitgrad)=unitcount(unitgrad)+1
if(printguess)unitcount(unitguess)=unitcount(unitguess)+1
if(printirf)unitcount(unitirf)=unitcount(unitirf)+1
if(printmargin)unitcount(unitmargin)=unitcount(unitmargin)+1
if(printmargins2)unitcount(unitmargin2)=unitcount(unitmargin2)+1
if(printmarginwtsum)unitcount(unitmarginwtsum)=unitcount(unitmarginwtsum)+1
if(printmargin2)unitcount(unitmargin2)=unitcount(unitmargin2)+1
if(printmeans)unitcount(unitmeans)=unitcount(unitmeans)+1
if(printmp)unitcount(unitml)=unitcount(unitml)+1
if(printmp)unitcount(unitmp)=unitcount(unitmp)+1
if(printobservedscale)unitcount(unitobservedscale)=unitcount(unitobservedscale)+1
if(printobsscale)unitcount(unitobsscale)=unitcount(unitobsscale)+1
if(printobsscalerel)unitcount(unitobsscalerel)=unitcount(unitobsscalerel)+1
if(printparam)unitcount(unitparam)=unitcount(unitparam)+1
if(printparamcov)unitcount(unitparamcov)=unitcount(unitparamcov)+1
if(printparamcov_complex)unitcount(unitparamcov_complex)=unitcount(unitparamcov_complex)+1
if(printparamcov_louis)unitcount(unitparamcov_louis)=unitcount(unitparamcov_louis)+1
if(printparamcov_sandwich)unitcount(unitparamcov_sandwich)=unitcount(unitparamcov_sandwich)+1
if(printpost)unitcount(unitpost)=unitcount(unitpost)+1
if(printpreditem)unitcount(unitpreditem)=unitcount(unitpreditem)+1
if(printprob)unitcount(unitprob)=unitcount(unitprob)+1
if(printprobdist)unitcount(unitprobdist)=unitcount(unitprobdist)+1
if(printpwt)unitcount(unitpwt)=unitcount(unitpwt)+1
if(printpwtdist)unitcount(unitpwtdist)=unitcount(unitpwtdist)+1
if(printrel)unitcount(unitrel)=unitcount(unitrel)+1
if(printrelskill)unitcount(unitrelskill)=unitcount(unitrelskill)+1
if(printreltrans)unitcount(unitreltrans)=unitcount(unitreltrans)+1
if(printrelwt)unitcount(unitrelwt)=unitcount(unitrelwt)+1
if(printreswt)unitcount(unitreswt)=unitcount(unitreswt)+1
if(printreswtdist)unitcount(unitreswtdist)=unitcount(unitreswtdist)+1
if(printscalerel)unitcount(unitscalerel)=unitcount(unitscalerel)+1
if(printthetaitem)unitcount(unitthetaitem)=unitcount(unitthetaitem)+1
if(printwtitem)unitcount(unitwtitem)=unitcount(unitwtitem)+1


if(printprog)unitcount(unititeration)=unitcount(unititeration)+1
if(printprogstart)unitcount(unititerationstart)=unitcount(unititerationstart)+1
unitcount(unitoutdata)=unitcount(unitoutdata)+1


unitcount(unittitle)=unitcount(unittitle)+1


return
end subroutine setoutput
