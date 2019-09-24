!	Get beta names.
!	Input:
!
!	factorname is the array of factor names.
!	itemname is the arrray of item names.
!	predname is the array of predictor names.
!	skillname is the array of skill names.
!	numcat is the number of underlying categories per item.
!   beta is the beta vector.
!
!	Output:
!
!	betanames is the array of beta names.

!
function betanames(factorname,itemname,predname,skillname,numcat,beta)
implicit none

character(len=32),intent(in)::itemname(:),skillname(:),factorname(:),predname(:)
integer,intent(in)::numcat(:)
real(kind=8)::beta(:)


character(len=64)::betanames(size(beta))

!	buff is used to create parameter names.
character(len=4)::buff
!	param-name is used for alternative parameter names.

!	cat counts categories.
!	col counts columns.

!	dimno is a dimension counters.
!	dimcovlat is the dimension of the covariance matrix of the latent vector.
!	dimlatin is the dimension of the latent vector.
!	dimlatout is the dimension of the transformed latent vector.
!	dim2 is a dimension counter/
!	item counts items.
!	nitems is the number of items.
!	npred is the number of predictors.
!	position records positions.
!	pred counts predictors.
!	row counts rows.
integer::cat,col,dimno,dimcovlat,dimlatin,dimlatout,dim2,&
	item,nitems,npred,position,pred,row


dimlatin=size(factorname)
dimlatout=size(skillname)


dimcovlat=dimlatin*(dimlatin+1)/2
nitems=size(numcat)
npred=size(predname)
!	Intercepts
row=1
position=1
do item=1,nitems
	do cat=0,numcat(item)-1
		write(buff,'(i4)') cat
		betanames(position)=trim(itemname(item))//"_intercept"//trim(adjustl(buff))
		position=position+1
	end do
end do
!	Slopes
do item=1,nitems
	do cat=0,numcat(item)-1
		do dimno=1,size(skillname)
			write(buff,'(i4)') cat
			betanames(position)=trim(adjustl(skillname(dimno)))//"_"//trim(itemname(item))//"_slope"//trim(adjustl(buff))
			position=position+1
		end do
	end do
end do
					
!	Linear part	
					
do pred=1,npred
	do dimno=1,dimlatin
		betanames(position)="Lin._"//trim(predname(pred))//"_"//trim(factorname(dimno))
		position=position+1
	end do
end do



!	Quadratic part
do pred=1,npred
	do dimno=1,dimlatin
		do dim2=1,dimno			
			betanames(position)="Quad._"//trim(predname(pred))//"_"//trim(factorname(dimno))&
							//"_"//trim(factorname(dimno))
			position=position+1
		end do
	end do
end do
return
end function betanames

