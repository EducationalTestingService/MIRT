!	Estimate the cross information matrix for a complex sample.
!	This version is designed for random sampling of psu's within strata.
!	Sampling is with replacement.
!	npsu is the number of psus's per stratum.
!	nstratum is the number of strata.
!	psu indicates the psu of an observation.
!	stratum indicates the stratum of an observation.
!	stratify is .true. for stratified sampling.
!	usepsu is .true. for primary sampling units within strata.
!	grads is the array of observation gradients.
!	grads1 is another array of observation gradients.
!	obsweight contains observation weights.

function crossinformation_complex(npsu,nstratum,psu,stratum,stratify,usepsu,grads,grads1,obsweight)
implicit none
!
integer,intent(in)::nstratum
integer,intent(in),optional::npsu(:),psu(:),stratum(:)
logical,intent(in)::stratify,usepsu
real(kind=8),intent(in)::grads(:,:),grads1(:,:),obsweight(:)
real(kind=8)::crossinformation_complex(size(grads,1),size(grads1,1))
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

!	diff is for differences for grads.
!	diff1 is for differences for grads1.
!	gradstrat is used for weighted summations within strata for grads.
!	grad1strat is used for weighted summations with strata for grads1
!	crossinform is used for cross information within a psu.

!	ratio is a ratio.

real(kind=8)::diff(size(grads,1)),diff1(size(grads1,1)),&
	gradstrat(size(grads,1)),grad1strat(size(grads1,1)),&
	crossinform(size(grads,1),size(grads1,1)),ratio
!	gradsum is for weighted summations within psu units for grads.
!	grad1sum is for weighted summations within psu units for grads1.
!	crossinformst is used for cross information within a stratum.
real(kind=8),allocatable::gradsum(:,:),grad1sum(:,:),crossinformst(:,:,:)


crossinformation_complex=0.0_8
!	Something must be present.
if(count(obsweight>0.0_8)==0)return
counts=0	
gradstrat=0.0_8
grad1strat=0.0_8	

if(usepsu)then

	allocate(gradsum(size(grads,1),sum(npsu)),grad1sum(size(grads1,1),sum(npsu)),stat=al)
	if(al/=0) stop "Unable to allocate arrays for asymptotic variances computations for complex sampling."
	gradsum=0.0_8
	grad1sum=0.0_8
	point=1
	do row=1,nstratum
		locpsu(row)=point
		point=point+npsu(row)
	end do
	if(stratify) then
		do obs=1,size(stratum)
			if(obsweight(obs)<=0.0_8) cycle
			point=locpsu(stratum(obs))+psu(obs)-1
			gradsum(:,point)=0.0_8
			grad1sum(:,point)=0.0_8
			do row=1,size(grads,1)
				if(grads(row,obs)/=0.0_8)gradsum(row,point)=gradsum(row,point)+obsweight(obs)*grads(row,obs)
			end do
			do row=1,size(grads1,1)
				if(grads1(row,obs)/=0.0_8)grad1sum(row,point)=grad1sum(row,point)+obsweight(obs)*grads1(row,obs)
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
			grad1strat=0.0_8
			do row1=1,size(grads1,1)
				do col=point,point+npsu(row)-1
					if(grad1sum(row1,col)/=0.0_8)grad1strat(row1)=grad1strat(row1)+grad1sum(row1,col)
				end do
			end do
			grad1strat=grad1strat/npsu(row)
			
			crossinform=0.0_8
			ratio=1.0_8
			if(npsu(row)>1)ratio=npsu(row)/(npsu(row)-1.0_8)			
			do col=1,npsu(row)
				diff=gradsum(:,point)-gradstrat
				diff1=grad1sum(:,point)-grad1strat
				do row1=1,size(grads,1)
					crossinform(row1,:)=crossinform(row1,:)&
						+diff(row1)*diff1
				end do
				point=point+1
			end do
			crossinformation_complex=crossinformation_complex+ratio*crossinform
			
			
			
			
		end do
	else
		do obs=1,size(psu)
			if(obsweight(obs)>0.0_8)then
				do row=1,size(grads,1)
					if(grads(row,obs)/=0.0_8)gradsum(row,psu(obs))=&
						gradsum(row,psu(obs))+obsweight(obs)*grads(row,obs)
				end do
				do row=1,size(grads1,1)
					if(grads1(row,obs)/=0.0_8)grad1sum(row,psu(obs))=&
						grad1sum(row,psu(obs))+obsweight(obs)*grads1(row,obs)
				end do
				
			end if
		end do
		gradstrat=0.0_8
		grad1strat=0.0_8
		do col=1,npsu(1)
			do row=1,size(grads,1)
				if(gradsum(row,col)/=0.0_8)gradstrat(row)=&
					gradstrat(row)+gradsum(row,col)
			end do
			do row=1,size(grads,1)
				if(grad1sum(row,col)/=0.0_8)grad1strat(row)=&
					grad1strat(row)+grad1sum(row,col)
			end do
		end do
		gradstrat=gradstrat/npsu(1)
		grad1strat=grad1strat/npsu(1)
		ratio=1.0_8
		if(npsu(1)>1)ratio=npsu(1)/(npsu(1)-1.0_8)
		
		
		do col=1,npsu(1)
			do row1=1,size(grads,1)
				diff(row1)=gradsum(row1,col)-gradstrat(row1)
			end do
			do row1=1,size(grads1,1)
				diff1(row1)=grad1sum(row1,col)-grad1strat(row1)
			end do
			do row=1,size(grads,1)
				if(diff(row)==0.0_8) cycle
				do row1=1,size(grads1,1)
					
					if(diff1(row1)/=0.0_8)crossinformation_complex(row,row1)&
						=crossinformation_complex(row,row1)+diff(row)*diff1(row1)
				end do
			end do
		end do
		crossinformation_complex=ratio*crossinformation_complex
		
	end if
else
	if(stratify) then
		allocate(gradsum(size(grads,1),nstratum),grad1sum(size(grads1,1),nstratum),&
			crossinformst(size(grads,1),size(grads1,1),nstratum),stat=al)
		if(al/=0) stop "Unable to allocate arrays for asymptotic variances computations for complex sampling."
		gradsum=0.0_8
		grad1sum=0.0_8
		crossinformst=0.0_8
		do obs=1,size(stratum)
			if(obsweight(obs)<=0.0_8)cycle
			gradsum(:,stratum(obs))=&
				gradsum(:,stratum(obs))+obsweight(obs)*grads(:,obs)
			grad1sum(:,stratum(obs))=&
				grad1sum(:,stratum(obs))+obsweight(obs)*grads1(:,obs)
			counts(stratum(obs))=counts(stratum(obs))+1
		end do
		do row=1,nstratum
			if(counts(row)>0)then
				gradsum(:,row)=gradsum(:,row)/counts(row)
				grad1sum(:,row)=grad1sum(:,row)/counts(row)
			end if
		end do
		do obs=1,size(stratum)
			if(obsweight(obs)<=0.0_8)cycle
			diff=obsweight(obs)*grads(:,obs)-gradsum(:,stratum(obs))
			diff1=obsweight(obs)*grads1(:,obs)-grad1sum(:,stratum(obs))
			do row=1,size(grads,1)
				if(diff(row)==0.0_8)cycle
				do col=1,size(grads1,1)
					if(diff1(col)/=0.0_8)crossinformst(row,col,stratum(obs))&
						=crossinformst(row,col,stratum(obs))+diff(row)*diff1(col)
				end do
			end do
		end do
		do row1=1,nstratum
			if(counts(row1)>1)then
				ratio=counts(row1)/(counts(row1)-1.0_8)
				crossinformation_complex=crossinformation_complex&
						+ratio*crossinformst(:,:,row1)

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
		grad1strat=0.0_8
		do obs=1,size(obsweight)
			if(obsweight(obs)==0.0_8)cycle
			do row=1,size(grads,1)
				if(grads1(row,obs)/=0.0_8)grad1strat(row)=grad1strat(row)+obsweight(obs)*grads1(row,obs)
			end do
		end do
		grad1strat=grad1strat/count(obsweight>0.0_8)
		
		
		
		do obs=1,size(obsweight)
			if(obsweight(obs)<=0.0_8)cycle
			diff=obsweight(obs)*grads(:,obs)-gradstrat
			diff1=obsweight(obs)*grads1(:,obs)-grad1strat
			do row=1,size(grads,1)
				if(diff(row)==0.0_8)cycle
				do col=1,size(grads1,1)
					if(diff1(col)/=0.0_8)crossinformation_complex(row,col)&
						=crossinformation_complex(row,col)+diff(row)*diff1(col)
				end do
			end do
		end do
		ratio=1.0_8
		if(count(obsweight>0.0_8)>1)ratio=count(obsweight>0.0_8)/(count(obsweight>0.0_8)-1.0_8)
		crossinformation_complex=ratio*crossinformation_complex
	end if
end if

return
end function crossinformation_complex

!
