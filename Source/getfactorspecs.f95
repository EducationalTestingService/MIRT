!	Obtain factor specifications.
!	normal is .true. if normal distribution and .false. otherwise.
!	factorname contains factor names.
!	fixdiag indicates that diagonal terms of quadrature component of log density are fixed.
!	fixquad indicates that quadrature component of log density is fixed.
!	independentf indicates that elements of latent vector are independent of previous elements.
!	no_lin indicates that the linear component of the log density has an intercept of 0.
!	noquad indicates that the quadratic componet of the log density is 0.
!	predlinmap is map of predictors to factors for linear term.
!	predquadmap is map of predictors to factors for quadratic term.	
subroutine getfactorspecs(normal,factorname,fixdiag,fixquad,independentf,no_lin,noquad,predlinmap,predquadmap)
implicit none
interface
    subroutine readfactorspecs(factor_name,fix_diag,independence,nolin,pred_lin_map,pred_quad_map)
        implicit none
        character(len=32),intent(inout)::factor_name
        logical,intent(inout)::fix_diag,independence,nolin
        logical,intent(inout)::pred_lin_map(:),pred_quad_map(:,:)

    end subroutine readfactorspecs
end interface
character(len=32),intent(out)::factorname(:)
logical,intent(in)::normal
logical,intent(out)::fixdiag(:),fixquad,independentf(:),no_lin(:),noquad,predlinmap(:,:),predquadmap(:,:,:)
!	buffer is a character buffer
character(len=4)::buffer
!	factor_name is for a factor name.
character(len=32)::factor_name
!	factor and factor1 count factors.
!	io is an input flag.

integer::factor,factor1,io
!	If factor_specs is .true., then individual factor information is read.
!	fix_diag specifies fixed values for diagonals for constant part
!		of quadratic component of log density.
!	fixquad indicates that the quadratic component of the logarithm of the latent
!		density is fixed for all independent vector values.
!	independence specifies independence of a component.
!	nolin specifies lack of an intercept for the linear component of
!		the log density.
!	noquad specifies lack of a quadratic component for the log density.
logical::factor_specs,fix_diag,independence,nolin
logical,allocatable::pred_lin_map(:),pred_quad_map(:,:)

namelist/allfactorspecs/factor_specs,fix_diag,fixquad,independence,nolin,noquad

allocate(pred_lin_map(size(predlinmap,2)),pred_quad_map(size(predquadmap,1),size(predquadmap,3)),stat=io)
if(io/=0)stop 'Allocatable of predictor maps failed.'

!	Default names.
do factor=1,size(factorname)
	write(buffer,'(i4)') factor
	factorname(factor)=trim("Factor"//adjustl(buffer))
end do
!	Defaults
factor_specs=.false.
fix_diag=.true.
fixquad=.true.
independence=.false.
nolin=.true.
noquad=.false.
predlinmap=.true.
predquadmap=.true.

if(size(factorname)>1)then
	do factor=2,size(factorname)
		do factor1=1,factor-1
			predquadmap(factor1,factor,:)=.false.
		end do
	end do
end if

!	General factor specifications.
read(*,nml=allfactorspecs,iostat=io)
if(io/=0)stop "General factor specifications cannot be read."
fixdiag=fix_diag
if(normal)noquad=.false.
no_lin=nolin

predlinmap(:,1)=.not.(nolin)

if(noquad)then
	predquadmap=.false.
	fixdiag=.false.
	fixquad=.false.
	independence=.false.
else
	if(fixquad.and.size(predquadmap,3)>1) predquadmap(:,:,2:size(predquadmap,2))=.false.
	if(independence.and.size(factorname)>1)then
		do factor=2,size(factorname)
			do factor1=1,factor-1
				predquadmap(factor,factor1,:)=.false.
			end do
		end do
	end if
end if
independentf=independence

!	Factor-specific setting.
if(factor_specs)then
	do factor=1,size(factorname)
		factor_name=factorname(factor)
		fix_diag=fixdiag(factor)
		nolin=no_lin(factor)
		pred_lin_map=predlinmap(factor,:)
		pred_quad_map=predquadmap(factor,:,:)
		independence=independentf(factor)
		
		
        call readfactorspecs(factor_name,fix_diag,independence,nolin,pred_lin_map,pred_quad_map)

		if(io/=0)stop "Factor specifications cannot be read."
		
		
		factorname(factor)=factor_name
		if(.not.noquad)fixdiag(factor)=fix_diag
		if(no_lin(factor).neqv.nolin) pred_lin_map(1)=.not.nolin
		no_lin(factor)=nolin
		if(.not.independentf(factor).and.independence)independentf(factor)=.true.
		
		predlinmap(factor,:)=pred_lin_map
		if(.not.noquad)predquadmap(factor,1:factor,:)=pred_quad_map(1:factor,:)
		
		
	end do
end if
!	Consistency settings.

if(size(factorname)>1)then
	do factor=2,size(factorname)
		do factor1=1,factor-1
			if(independentf(factor).or.independentf(factor1)) predquadmap(factor,factor1,:)=.false.
			
			
		end do
	end do
end if

if(fixquad.and.size(predquadmap,3)>1) then
	predquadmap(:,:,2:size(predquadmap,3))=.false.
	if(normal)then
		do factor=1,size(factorname)
			predquadmap(factor,factor,1)=.true.
		end do
	end if
end if			


return
end subroutine getfactorspecs		
subroutine readfactorspecs(factor_name,fix_diag,independence,nolin,pred_lin_map,pred_quad_map)
    implicit none
    character(len=32),intent(inout)::factor_name
    logical,intent(inout)::fix_diag,independence,nolin
    logical,intent(inout)::pred_lin_map(:),pred_quad_map(:,:)
    integer::io
    namelist/factorspecs/factor_name,fix_diag,independence,nolin,pred_lin_map,pred_quad_map
    read(*,nml=factorspecs,iostat=io)
    if(io/=0)stop "Factor specifications cannot be read."
    return
end subroutine readfactorspecs
