!	Get design settings.
!	Input:
!
!	factorname is the array of factor names.
!	itemname is the arrray of item names.
!	predname is the array of predictor names.
!	skillname is the array of skill names.
!	choices indicates the number of options for an item.
!	intdim gives the number intercept components per item.
!	numcat is the number of underlying categories per item.
!	slopedim gives the number of slope components per dimension and per item.
!	constguess is the indicator for a  constant guessing parameter.
!	fixdiag is the indicator for a fixed diagonal for the first predictor.
!	fixedguess is the indicator for a fixed guessing parameter.
!	guess is the indicator for guessing.
!	rasch is the indicator for a Rasch model.
!	raschslope1 is the indicator for a Rasch model with a slope of 1.
!	specioltrans specifies a custom transition matrix.!	itemmatrix is the item matrix.
!	setguess indicates the fixed value of the guessing parameter.
!	setslope specifies initialization of slope parameters.

!
!	Output:
!
!	paramname is the array of parameter names.
!	design is the final design matrix.
!	gamma is the parameter vector.
!	offset is the final offset.
!
subroutine getdesign(factorname,itemname,predname,skillname,choices,intdim,&
	numcat,slopedim,&
	constguess,fixdiag,fixedguess,guess,predlinmap,predquadmap,rasch,raschslope1,&
	specialtrans,itemmatrix,setguess,setslope,&
	paramname,design,gamma,offset)
implicit none
interface
    subroutine readdesignspecs(param_name,offsettran,transition)
        character(len=64),intent(inout)::param_name(:)
        real(kind=8),intent(inout)::offsettran(:),transition(:,:)
    end subroutine readdesignspecs
end interface

character(len=32),intent(in)::itemname(:),skillname(:),factorname(:),predname(:)
integer,intent(in)::choices(:),intdim(:),numcat(:),slopedim(:,:)
logical,intent(in)::constguess,fixdiag(:),fixedguess(:),guess(:),&
	predlinmap(:,:),predquadmap(:,:,:),rasch(:),raschslope1(:),specialtrans
real(kind=8),intent(in)::itemmatrix(:,:),setguess(:),setslope(:,:)
character(len=64),intent(out)::paramname(:)
real(kind=8),intent(out)::design(:,:),gamma(:),offset(:)
!	buff is used to create parameter names.
character(len=4)::buff
!	param-name is used for alternative parameter names.
character(len=64),allocatable::param_name(:)
!	cat counts categories.
!	col counts columns.
!	colguess records the location of the guessing parameter for a constant guessing parameter.

!	colrasch records the location of the common slope parameters.
!	dimno is a dimension counters.
!	dimcovlat is the dimension of the covariance matrix of the latent vector.
!	dimlatin is the dimension of the latent vector.
!	dimlatout is the dimension of the transformed latent vector.
!	dim2 is a dimension counter/
!	io is the error flag for input.
!	item counts items.
!	npred is the number of predictors.
!	position records positions.
!	pred counts predictors.
!	row counts rows.
integer::cat,col,colguess,colrasch(size(factorname)),dimno,dimcovlat,dimlatin,dimlatout,dim2,&
	io,item,npred,position,pred,row
!	previous checks if a previous guessing parameter has been encountered, and prevrasch checks if
!	a previous slope parameter has been found.
logical::previous,prevrasch(size(rasch))
!	diff is a difference.

real(kind=8)::diff	
real(kind=8),allocatable::offsettran(:),transition(:,:)

dimlatin=size(factorname)
dimlatout=size(slopedim,1)
npred=size(predlinmap,2)


allocate(offsettran(size(itemmatrix,2)),&
	transition(size(itemmatrix,2),size(design,2)),stat=io)
if(io/=0)stop "Space for transition matrix and offset for this matrix not successfully allocated"
!	Default settings.
gamma=0.0_8
offsettran=0.0_8
transition=0.0_8

do row=1,size(design,2)
	write(buff,'(i4)') row
	paramname(row)="Parameter"//trim(adjustl(buff))
end do
!	Take care of cases with unrestricted intercepts and some straightforward handling of guessing.
row=1
col=1
position=1
if(.not.specialtrans)then
   if(.not.constguess)then
	do item=1,size(intdim)
		if(.not.guess(item))then
			do cat=1,intdim(item)
				transition(row,col)=1.0_8
				if(intdim(item)==1)then
					paramname(col)=trim(itemname(item))//"_intercept"
				else
					
					write(buff,'(i4)') cat
					paramname(col)=trim(itemname(item))//"_intercept"//trim(adjustl(buff))
				end if
				col=col+1
				
				row=row+1
			end do
			
			
		else
			transition(row,col)=1.0_8
			paramname(col)=trim(itemname(item))//"_intercept"
			row=row+1
			col=col+1
			diff=itemmatrix(position+3,row-1)-itemmatrix(position+2,row-1)
			if(.not.fixedguess(item)) then
				transition(row,col)=1.0_8
				paramname(col)=trim(itemname(item))//"_guesslogit"
				if(diff/=0.0_8)gamma(col)=setguess(item)/diff
				col=col+1
			else
				
				if(diff/=0.0_8)offsettran(row)=setguess(item)/diff
			end if
			row=row+1
			
		end if
		position=position+numcat(item)
	end do
   else
	previous=.false.
	do item=1,size(intdim)
		if(.not.guess(item))then
			do cat=1,intdim(item)
				transition(row,col)=1.0_8
				if(intdim(item)==1)then
					paramname(col)=trim(itemname(item))//"_intercept"
				else
					write(buff,'(i4)') cat
					paramname(col)=trim(itemname(item))//"_intercept"//trim(adjustl(buff))
				end if
				col=col+1
				row=row+1
			end do
		else
			transition(row,col)=1.0_8
			paramname(col)=trim(itemname(item))//"_intercept"
			row=row+1
			col=col+1
			if(previous)then
				transition(row,colguess)=1.0_8
			else
				previous=.true.
				colguess=col
				paramname(col)="Guesslogit"
				transition(row,col)=1.0_8
				diff=itemmatrix(position+3,row-1)-itemmatrix(position+2,row-1)
				if(diff/=0.0_8)gamma(col)=setguess(item)/diff
				col=col+1
			end if
			row=row+1
		end if
		position=position+numcat(item)
	end do
   end if		

!	Handle combinations of unrestricted slopes and simple Rasch cases.


   prevrasch=.false.

   do item=1,size(intdim)
	do dimno=1,dimlatout
		do cat=1,slopedim(dimno,item)
			if(rasch(dimno))then
				if(raschslope1(dimno))then
					offsettran(row)=1.0_8
				else
					if(.not.prevrasch(dimno))then
						colrasch(dimno)=col
						prevrasch(dimno)=.true.
						paramname(col)=trim(skillname(dimno))//"_slope"
						gamma(col)=setslope(dimno,item)
						
						col=col+1
					end if
					transition(row,colrasch(dimno))=1.0_8
				end if
			else
				transition(row,col)=1.0_8
				if(slopedim(dimno,item)==1) then
					paramname(col)=&
						trim(skillname(dimno))//"_"//trim(itemname(item))//"_slope"
					gamma(col)=setslope(dimno,item)
				else
					write(buff,'(i4)') cat
					paramname(col)=&
						trim(skillname(dimno))//"_"//trim(itemname(item))//"_slope"//&
						trim(adjustl(buff))
					gamma(col)=setslope(dimno,item)/slopedim(dimno,item)
				end if
				
				col=col+1
			end if
			row=row+1
		end do
		
	end do
   end do
   !	Linear part.
	do pred=1,npred
								
!	Default linear part.
		do dimno=1,dimlatin				
			if(predlinmap(dimno,pred))then
				transition(row,col)=1.0_8
				if(npred==1)then
					paramname(col)="Lin._"//trim(factorname(dimno))
					
					
				else
					
					transition(row,col)=1.0_8
					paramname(col)="Lin._"//trim(predname(pred))//"_"//trim(factorname(dimno))
					
				end if
				col=col+1
			end if
			row=row+1
			
			
		end do
	end do



!	Quadratic part

do pred=1,npred
	do dimno=1,dimlatin
		do dim2=1,dimno
			if(predquadmap(dimno,dim2,pred)) then
				if(pred>1.or.dimno/=dim2.or.(.not.fixdiag(dimno))) then
					transition(row,col)=1.0_8
					if(npred==1) then
						paramname(col)="Quad._"//trim(factorname(dimno))&
							//"_"//trim(factorname(dim2))
					else
						paramname(col)="Quad._"//trim(predname(pred))//"_"//trim(factorname(dimno))&
							//"_"//trim(factorname(dim2))
					end if
					if(pred==1.and.dimno==dim2)gamma(col)=-0.5_8
					col=col+1
				else
					offsettran(row)=-0.5_8	
					
				end if
				
			end if
			row=row+1
		end do
	end do
end do

else
	allocate(param_name(size(design,2)),stat=io)
    if(io/=0)stop "Space for alternative parameter names not successfully allocated"

	do col=1,min(size(transition,1),size(transition,2))
		transition(col,col)=1.0_8
	end do
	
	param_name=paramname
    call readdesignspecs(param_name,offsettran,transition)
    paramname=param_name
end if					

design=matmul(itemmatrix,transition)
offset=matmul(itemmatrix,offsettran)


return
end subroutine getdesign
subroutine readdesignspecs(param_name,offsettran,transition)
    character(len=64),intent(inout)::param_name(:)
    real(kind=8),intent(inout)::offsettran(:),transition(:,:)
    integer::io
    namelist/designspecs/param_name,offsettran,transition
    read(*,nml=designspecs,iostat=io)

    if(io/=0)stop "Parameter names, design matrix, and/or design vector not successfully read."
    return
end subroutine readdesignspecs
