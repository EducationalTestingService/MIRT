!	Estimate the information matrix for a complex sample.
!	This version is designed for random sampling of psu's within strata.
!	Sampling is with replacement.
!	npsu is the number of psus's per stratum.
!	nstratum is the number of strata.
!	psu indicates the psu of an observation.
!	stratum indicates the stratum of an observation.
!	stratify is .true. for stratified sampling.
!	usepsu is .true. for primary sampling units within strata.
!	grads is the array of observation gradients.
!	obsweight contains observation weights.

function information_complex(npsu,nstratum,psu,stratum,stratify,usepsu,grads,obsweight)
implicit none
!
interface
	subroutine makesym(matrix)
!	Symmetrize n by n matrix based on lower triangle.
		real(kind=8),intent(inout)::matrix(:,:)
	end subroutine makesym
end interface
integer,intent(in)::nstratum
integer,intent(in),optional::npsu(:),psu(:),stratum(:)
logical,intent(in)::stratify,usepsu
real(kind=8),intent(in)::grads(:,:),obsweight(:)
real(kind=8)::information_complex(size(grads,1),size(grads,1))
!	al indicates allocation error.
!	col counts columns.
!	col1 counts columns.
!	counts count strata.
!	locpsu is psu locator. 
!	obs is an observation counter.
!	point indicates the location of a psu.
!	row counts rows.
!	row1 counts rows.
integer::al,col,col1,counts(nstratum),locpsu(nstratum),obs,point,row,row1

!	diff is for differences.
!	gradstrat is used for weighted summations within strata.
!	inform is used for information within a psu.

!	ratio is a ratio.

real(kind=8)::diff(size(grads,1)),gradstrat(size(grads,1)),inform(size(grads,1),size(grads,1)),ratio
!	gradsum is for weighted summations within psu units.
!	informst is used for information within a stratum.
real(kind=8),allocatable::gradsum(:,:),informst(:,:,:)


information_complex=0.0_8
!	Something must be present.
if(count(obsweight>0.0_8)==0)return
counts=0	
gradstrat=0.0_8	

if(usepsu)then

	allocate(gradsum(size(grads,1),sum(npsu)),stat=al)
	if(al/=0) stop "Unable to allocate arrays for asymptotic variances computations for complex sampling."
	gradsum=0.0_8
	
	point=1
	do row=1,nstratum
	
		locpsu(row)=point
		point=point+npsu(row)
	end do
	
	
	
	if(stratify) then
		do obs=1,size(obsweight)
			if(obsweight(obs)<=0.0_8) cycle
			point=locpsu(stratum(obs))+psu(obs)-1
			
			do row=1,size(grads,1)
				if(grads(row,obs)/=0.0_8)gradsum(row,point)=gradsum(row,point)+obsweight(obs)*grads(row,obs)
			end do
			
		end do
		
		point=1
		do row=1,nstratum
			gradstrat=0.0_8
				do row1=1,size(grads,1)
					do col=point,point+npsu(row)-1
						if(gradsum(row1,col)/=0.0_8)gradstrat(row1)=gradstrat(row1)+gradsum(row1,col)
					end do
				end do
			gradstrat=gradstrat/npsu(row)
			inform=0.0_8
			ratio=1.0_8
			if(npsu(row)>1) ratio=npsu(row)/(npsu(row)-1.0_8)	
					
			do col=1,npsu(row)
				diff=gradsum(:,point)-gradstrat
				do row1=1,size(grads,1)
					if(diff(row1)==0.0_8)cycle
					do col1=1,row1
						if(diff(col1)/=0.0_8)inform(row1,col1)=inform(row1,col1)&
							+diff(row1)*diff(col1)
					end do
				end do
				point=point+1
			end do
			do row1=1,size(grads,1)
				do col1=1,row1
					if(inform(row1,col1)/=0.0_8)information_complex(row1,col1)=information_complex(row1,col1)&
						+ratio*inform(row1,col1)
				end do
			end do
			
			
			
		end do
		
	else
	
		do obs=1,size(obsweight)
			if(obsweight(obs)==0.0_8)cycle
			do row=1,size(grads,1)
				
				if(grads(row,obs)/=0.0_8)gradsum(row,psu(obs))=&
					gradsum(row,psu(obs))+obsweight(obs)*grads(row,obs)
			end do
		end do
		gradstrat=0.0_8
		
		do col=1,npsu(1)
			do row=1,size(grads,1)
				if(gradsum(row,col)/=0.0_8) gradstrat(row)=gradstrat(row)+gradsum(row,col)
			end do
		end do
		
		gradstrat=gradstrat/npsu(1)	
		ratio=1
		if(npsu(1)>1)ratio=npsu(1)/(npsu(1)-1.0_8)
		do col=1,npsu(1)
			do row1=1,size(grads,1)
				diff(row1)=gradsum(row1,col)-gradstrat(row1)
				if(diff(row1)==0.0_8)cycle
				do col1=1,row1
					if(diff(col1)/=0.0_8)information_complex(row1,col1)=information_complex(row1,col1)+diff(row1)*diff(col1)
				end do
			end do
		end do
		
		do row=1,size(grads,1)
			information_complex(row,1:row)=ratio*information_complex(row,1:row)
		end do
		
	end if
else
	if(stratify) then
		allocate(gradsum(size(grads,1),nstratum),informst(size(grads,1),size(grads,1),nstratum),stat=al)
		if(al/=0) stop "Unable to allocate arrays for asymptotic variances computations for complex sampling."
		gradsum=0.0_8
		informst=0.0_8
		do obs=1,size(obsweight)
			if(obsweight(obs)<=0.0_8)cycle
			do row=1,size(grads,1)
				if(grads(row,obs)/=0.0_8)gradsum(row,stratum(obs))=&
					gradsum(row,stratum(obs))+obsweight(obs)*grads(row,obs)
			end do
			counts(stratum(obs))=counts(stratum(obs))+1
		end do
		do row=1,nstratum
			if(counts(row)>0)gradsum(:,row)=gradsum(:,row)/counts(row)
		end do
		do obs=1,size(obsweight)
			if(obsweight(obs)<=0.0_8)cycle
			
			do row=1,size(grads,1)
				diff(row)=obsweight(obs)*grads(row,obs)-gradsum(row,stratum(obs))
				if(diff(row)==0.0_8)cycle
				do col=1,row
					if(diff(col)/=0.0_8)informst(row,col,stratum(obs))&
						=informst(row,col,stratum(obs))+diff(row)*diff(col)
				end do
			end do
		end do
		do row1=1,nstratum
			if(counts(row1)>1)then
				ratio=counts(row1)/(counts(row1)-1.0_8)
				do row=1,size(grads,1)
					do col=1,row
						if(informst(row,col,row1)/=0.0_8)information_complex(row,col)&
							=information_complex(row,col)+ratio*informst(row,col,row1)
					end do
				end do
			end if
		end do
	else
		gradstrat=0.0_8
		do obs=1,size(obsweight)
			if(obsweight(obs)==0.0_8)cycle
			do row=1,size(grads,1)
				if(grads(row,obs)/=0.0_8)gradstrat(row)=gradstrat(row)+obsweight(obs)*grads(row,obs)
			end do
		end do
		gradstrat=gradstrat/count(obsweight>0.0_8)
		
		do obs=1,size(obsweight)
			if(obsweight(obs)<=0.0_8)cycle
			
			do row=1,size(grads,1)
				diff(row)=obsweight(obs)*grads(row,obs)-gradstrat(row)
				if(diff(row)==0.0_8)cycle
				do col=1,row
					if(diff(col)/=0.0_8)information_complex(row,col)=information_complex(row,col)+diff(row)*diff(col)
				end do
			end do
		end do
		
		ratio=1.0_8
		if(count(obsweight>0.0_8)>1)ratio=count(obsweight>0.0_8)/(count(obsweight>0.0_8)-1.0_8)
		do row=1,size(grads,1)
			information_complex(row,1:row)=ratio*information_complex(row,1:row)
		end do
	end if
end if

call makesym(information_complex)
return
end function information_complex
