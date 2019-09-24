!	Map beta.  First row is item, second is category, third is first dimension, fourth is second dimension, and fifth is
!	predictor number.
!	dimlatin is dimension of latent vector.
!	dimlatout is dimension of transformed latent vector.
!	npred is number of predictors.
!	numcat is number of underlying categories per item.
function mapbeta(dimlatin,dimlatout, npred,numcat)
implicit none
integer,intent(in)::dimlatin,dimlatout,npred,numcat(:)
integer::mapbeta(5,sum(numcat)*(dimlatout+1)+npred*dimlatin*(dimlatin+3)/2)
!	position, position1, cat, item, dimno, dim1, and pred are used as counters.
integer::cat,dimno,dim1,item,position,position1, pred

position=1
mapbeta=0
!	Intercepts.
do item=1,size(numcat)
	do cat=1,numcat(item)
		mapbeta(1,position)=item
		mapbeta(2,position)=position
		position=position+1
	end do
end do
!	Slopes.
position1=1
do item=1,size(numcat)
	do cat=1,numcat(item)
		do dimno=1,dimlatout
			mapbeta(1,position)=item
			mapbeta(2,position)=position1	
			mapbeta(3,position)=dimno
			position=position+1
		end do
		position1=position1+1
	end do
end do
!	Linear components
do pred=1,npred
	do dimno=1,dimlatin
		mapbeta(3,position)=dimno
		mapbeta(5,position)=pred
		position=position+1
	end do
end do
!	Quadratic components
do pred=1,npred
	do dimno=1,dimlatin
		do dim1=1,dimno
			mapbeta(3,position)=dimno
			mapbeta(4,position)=dim1
			mapbeta(5,position)=pred
			position=position+1
		end do
	end do
end do
return
end function mapbeta
