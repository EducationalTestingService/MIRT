
!	Get item details such as category scores and maps of observed
!	categories to underlying categories.
!	intdim provides intercept dimensions for the item matrix.
!	numcat is the category counts for observed categories.
!	numcatobs is the category counts for underlying categories.
!	slopedim provides slope dimensions for the item matrix.
!	catmap is for category maps.
!	guess specifies guessing.
!	nominal is the indicator for a nominal model.
!	specialitemint specifies special treatment of an item intercept.
!	specialitemslope specifies special treatment of an item slope.
!	catobsrange provides the category map.
!	itemmatrix is the item matrix.
!	setslope helps initialize slope parameters.
subroutine getitemdetail(intdim,numcat,numcatobs,slopedim,&
	catmap,guess,nominal,specialitemint,&
	specialitemslope,catobsrange,itemmatrix,setslope)
implicit none
interface
    subroutine readcatarray(cat_array)
        implicit none
        integer,intent(inout)::cat_array(:)
    end subroutine readcatarray
    subroutine readintarray(int_array)
        implicit none
        real(kind=8),intent(inout)::int_array(:,:)
    end subroutine readintarray
    subroutine readslopearray(slope_array)
        implicit none
        real(kind=8),intent(inout)::slope_array(:,:)
    end subroutine readslopearray


end interface
integer,intent(in)::intdim(:),numcat(:),numcatobs(:),slopedim(:,:)
integer,intent(out)::catobsrange(:,:)
logical,intent(in)::catmap(:),guess(:),nominal(:),specialitemint(:),specialitemslope(:)
real(kind=8),intent(out)::itemmatrix(:,:),setslope(:,:)
!	cat and cat1 count categories.
!	counter, counter1, and counter2 are general counters.
!	countercat and countercatobs count positions for catobsrange.
!	counterdens and counterdenso are for entries for
!	densities.
!	counterint and counterinto count positions for intercepts.
!	counterslope and counterslopeo count positions for slopes.
!	dimcount is number of dimensions for an item.
!	dimno is a dimension.
!	dimlatout is the dimension of the latent vector.
!	io is the indicator for read error.
!	item is an item.
!	ncat is the total number of underlying categories.
integer::cat,cat1,counter,counter1,counter2,countercat,countercatobs,&
	counterdens,counterdenso,counterint,counterinto,&
	counterslope,counterslopeo,&
	dimcount,dimlatout,dimno,item,io,ncat,nitems
!	cat_array is used to specify category boundaries.
integer,allocatable::cat_array(:)
!	int_array specifies an intercept matrix for an item.
!	slope_array specifies a slope array for a skill
!	and item.
real(kind=8),allocatable::int_array(:,:),slope_array(:,:)


countercat=1
countercatobs=1
counterint=1
counterinto=1
ncat=sum(numcat)
counterslope=ncat+1
counterslopeo=sum(intdim)+1
dimlatout=size(slopedim,1)
counterdens=counterslope+ncat*dimlatout
counterdenso=counterslopeo+sum(slopedim)


nitems=size(numcat)
setslope=0.0_8
itemmatrix=0.0_8

do item=1,nitems
	
!	Default category mapping.

	do cat=1,numcatobs(item)-1
		catobsrange(:,countercatobs)=countercat
		countercat=countercat+1
		countercatobs=countercatobs+1
	end do
	catobsrange(1,countercatobs)=countercat
	countercat=countercat+numcat(item)-numcatobs(item)
	catobsrange(2,countercatobs)=countercat
	countercat=countercat+1
	countercatobs=countercatobs+1
!	Special category map.
	if(catmap(item)) then
		allocate(cat_array(numcatobs(item)-1),stat=io)
		cat_array=catobsrange(1,countercatobs-numcatobs(item)+1:countercatobs-1)&
			-catobsrange(1,countercatobs-numcatobs(item))
		
		if(io/=0) stop "Category mapping array not allocated."
        call readcatarray(cat_array)

!		Force mapping to be feasible.
		cat_array(1)=min(max(1,cat_array(1)),numcat(item)-numcatobs(item)+1)
		if(numcatobs(item)>2)then
			do cat=2,numcatobs(item)-1
				cat_array(cat)=min(max(cat_array(cat-1)+1,cat_array(cat)),numcat(item)-numcatobs(item)+cat)
			end do
		end if
		catobsrange(1,countercatobs-numcatobs(item)+1:countercatobs-1)&
			=catobsrange(1,countercatobs-numcatobs(item))+cat_array
		catobsrange(2,countercatobs-numcatobs(item):countercatobs-2)&
			=catobsrange(1,countercatobs-numcatobs(item)+1:countercatobs-1)-1
		
		
		deallocate(cat_array)
	end if


!	Intercepts


!	Standard generalized partial credit or nominal approach.
	if(guess(item))then
		itemmatrix(counterint+1,counterinto)=1.0_8
		itemmatrix(counterint+3,counterinto)=1.0_8
		itemmatrix(counterint+2:counterint+3,counterinto+1)=1.0_8
		counterint=counterint+4
		counterinto=counterinto+2
	else
		
		do cat=1,intdim(item)
			itemmatrix(counterint+cat:counterint+numcat(item)-1,counterinto)=1.0_8
			counterinto=counterinto+1
		end do
		counterint=counterint+numcat(item)
	end if
!	Custom matrix.
	if(specialitemint(item))then
		allocate(int_array(numcat(item),intdim(item)),stat=io)
		if(io/=0)stop "Item intercept array not allocated."
		int_array=itemmatrix(counterint-numcat(item):counterint-1,counterinto-intdim(item):counterinto-1)
        call readintarray(int_array)

		itemmatrix(counterint-numcat(item):counterint-1,counterinto-intdim(item):counterinto-1)&
			= int_array
		deallocate(int_array)
	end if
!	Slopes
	counter1=counterslope
	counter2=counterslopeo
	if(guess(item))then
		counterslope=counterslope+dimlatout
		dimcount=count(slopedim(:,item)>0)
		do dimno=1,dimlatout
			if(slopedim(dimno,item)==1)then
				itemmatrix(counterslope,counterslopeo)=1.0_8
				itemmatrix(counterslope+2*dimlatout,counterslopeo)=1.0_8
				counterslopeo=counterslopeo+1
				setslope(dimno,item)=1.0_8/dimcount
			end if
			counterslope=counterslope+1
		end do
		counterslope=counterslope+2*dimlatout
	else
		do dimno=1,dimlatout
			dimcount=count(slopedim(:,item)>0)*(numcat(item)-1)
			do counter=1,slopedim(dimno,item)
!	Generalized partial credit default.
				if(slopedim(dimno,item)==1)then
					do cat=2,numcat(item)
						itemmatrix(counterslope+dimlatout*(cat-1)+dimno-1,counterslopeo)=cat-1
					end do
					setslope(dimno,item)=1.0_8/dimcount
					
				else
!	Nominal default
					do cat=counter+1,numcat(item)
						itemmatrix(counterslope+dimlatout*(cat-1)+dimno-1,counterslopeo)=1.0_8
					end do
				end if
				counterslopeo=counterslopeo+1
			end do
			if(slopedim(dimno,item)>0)setslope(dimno,item)=1.0_8/dimcount
		end do
		counterslope=counterslope+dimlatout*numcat(item)
	end if	
!	Custom array
	if(specialitemslope(item))then
		do dimno=1,size(slopedim,1)
			if(slopedim(dimno,item)>0) then
				allocate(slope_array(numcat(item),slopedim(dimno,item)),stat=io)
				if(io/=0) stop "Item slope array not allocated successfully"
				do counter=1,slopedim(dimno,item)
					do cat=1,numcat(item)
						slope_array(cat,counter)=&
						itemmatrix(counter1+dimno-1+(cat-1)*dimlatout,counter2+counter-1)
					end do
					
				end do
                call readslopearray(slope_array)
				do counter=1,slopedim(dimno,item)
					do cat=1,numcat(item)
						itemmatrix(counter1+dimno-1+(cat-1)*dimlatout,&
							counter2+counter-1)&
							=slope_array(cat,counter)
					end do
					
				end do
				deallocate(slope_array)	
			end if
		end do
	end if
end do


!	Fill out itemmatrix with entries for density of latent vector.

do counter=counterdens,size(itemmatrix,1)
	itemmatrix(counter,counterdenso)=1.0_8
	
	counterdenso=counterdenso+1
end do	
return
end subroutine getitemdetail
subroutine readcatarray(cat_array)
    implicit none
    integer,intent(inout)::cat_array(:)
    integer::io
    namelist/catspecs/cat_array


    read(*,nml=catspecs,iostat=io)

    if(io/=0) stop "Categery mapping array not read successfully."


    return
end subroutine readcatarray
subroutine readintarray(int_array)
    implicit none
    real(kind=8),intent(inout)::int_array(:,:)
    integer::io
    namelist/intspecs/int_array
    read(*,nml=intspecs,iostat=io)
    if(io/=0)stop "Item intercept array not read successfully."
    return
end subroutine readintarray
subroutine readslopearray(slope_array)
    implicit none
    real(kind=8),intent(inout)::slope_array(:,:)
    integer::io
    namelist/slopespecs/slope_array
    read(*,nml=slopespecs,iostat=io)
    if(io/=0) stop "Item slope array not read successfully."

    return
end subroutine readslopearray

