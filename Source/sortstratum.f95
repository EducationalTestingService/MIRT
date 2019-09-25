!
!	sortstratum is used to sort the strata and/or psus and then to
!	recode them so that strata are numbered from 1 to nstratum,
!	the number of strata and psus within strata are numbered from
!	1 to the number of psus within the stratum.
!	psu is the array of psu codes.
!	stratum is the array of stratum codes.
!	stratify is .true. for stratified sampling.
!	usepsu is .true. if psu's are used.
subroutine sortstratum(psu,stratum,stratify,usepsu,nstratum)
implicit none
interface
!	Integer sort in ascending order.
	recursive subroutine sort(a)
		implicit none
		integer, intent(in out):: a(:)
		integer :: pivot
	end subroutine sort
!	Sort array a of integer pairs in ascending order.
	recursive subroutine sort8(a)
		implicit none
		integer, intent(in out) :: a(:,:)
	end subroutine sort8
!	Remove redundant integers.
	subroutine reduce(na,a)
		implicit none
		integer, intent(in out)::a(:)
		integer,intent(out)::na
	end subroutine
!	reduce a sorted array a of integer pairs to distinct entries.
	subroutine reduce8(na,a)
		implicit none
		integer, intent(in out)::a(:,:)
		integer,intent(out)::na
	end subroutine reduce8
!	Look up integers.	
	integer function lookup(x,a)
		implicit none
		integer, intent(in) :: x,a(:)
	end function lookup
	!	Look up integer pairs.
	integer function lookup8(x,a)
		implicit none
		integer, intent(in) :: x(2),a(:,:)
	end function lookup8
	

end interface

integer,intent(inout)::psu(:),stratum(:)
logical,intent(in)::stratify,usepsu
integer,intent(out)::nstratum
!	al is allocation error flag.
!	na is number of distinct entries.
!	obs is observation number.

!	q is table location.
!	qq is entry to look up.



integer::al,na,obs,q,qq(2)
!	Copies for sorting.
integer,allocatable::psuc(:),stratumc(:)
integer, allocatable::psustr(:,:),b(:,:)

!	Stratified and psus.

if(stratify)then
	if(usepsu)then
		allocate(psustr(2,size(psu)),stat=al)
		if(al/=0) stop "Cannot allocate arrays to sort psu and stratum data."
		
		
		do obs=1,size(stratum)
			 
			psustr(1,obs)=stratum(obs)
			psustr(2,obs)=psu(obs)
			
		end do
		
		call sort8(psustr)
		
		call reduce8(na,psustr)
		
		allocate(b(2,na),stat=al)
		if(al/=0) stop "Cannot allocate arrays to sort psu and stratum data."
		do q=1,na
			if(q==1)then
				b(:,q)=1
			else
				if(psustr(1,q)>psustr(1,q-1))then
					b(1,q)=b(1,q-1)+1
					b(2,q)=1
				else
					b(1,q)=b(1,q-1)
					b(2,q)=b(2,q-1)+1
				end if
			end if
		end do
		nstratum=b(1,na)
		
		do obs=1,size(stratum)
			qq(1)=stratum(obs)
			qq(2)=psu(obs)
			q=lookup8(qq,psustr(:,1:na))
			stratum(obs)=b(1,q)
			psu(obs)=b(2,q)
			
			
			
		end do
		
	else
!	Just stratified.
		allocate(stratumc(size(stratum)),stat=al)
		if(al/=0) stop "Cannot allocate arrays to sort stratum data."
		
		stratumc=stratum
		call sort(stratumc)
		call reduce(na,stratumc)
		do obs=1,size(stratum)
			stratum(obs)=lookup(stratum(obs),stratumc(1:na))
		end do
		nstratum=na
		
	end if
	
else
	nstratum=1
	allocate(psuc(size(psu)),stat=al)
	if(al/=0)stop "Failure to allocate space for sorting primary sampling units."
	psuc=psu
	
	call sort(psuc)
	
	call reduce(na,psuc)
	
	do obs=1,size(psu)
			psu(obs)=lookup(psu(obs),psuc(1:na))
			
	end do
	
	
	
	
	
end if
return
end subroutine sortstratum