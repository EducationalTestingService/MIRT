


!	Obtain basic item data.
!	itemname contains item names.
!	choices indicates the number of choices in a right-scored multiple-choice item.
!	intdim  is the dimension of the intercept matrix.
!	numcat is the number of underlying categories per item.
!	numcatobs is the number of observed categories per item.
!	skillnumbers specifies skills that correspond to items in a between-item case.
!	slopedim is the dimension of the slope matrix.
!	betweenitem indicates a between-item model.
!	catmap is for category maps.
!	constguess indicates a constant guessing parameter.
!	fixedguess indicates a fixed guessing parameter.
!	guess indicates a guessing parameter.
!	nominal indicates that a nominal model is used for the item.
!	specialitemint is used to specify items that need special treatment for intercepts.
!	specialitemslope is used to specify items that need special treatment for slopes.
!	setguess is the fixed setting of a guessing parameter.

subroutine getitemdata(itemname,choices,intdim,numcat,numcatobs,skillnumbers,slopedim,&
	betweenitem,catmap,constguess,fixedguess,guess,nominal,&
	specialitemint,specialitemslope,setguess)
implicit none
interface
    subroutine readallitemspecs(int_dim,num_choices,num_cat,num_cat_obs,&
        slope_dim,between_item,cat_map,constguess,fixguess,guessing,&
        item_specs,nommod,special_item_int,special_item_slope,set_guess)
        implicit none
        integer,intent(inout)::int_dim,num_cat,&
            num_cat_obs,num_choices
        integer,intent(inout)::slope_dim(:)


        logical,intent(inout)::between_item,cat_map,constguess,fixguess,guessing,item_specs,&
            nommod,special_item_int,special_item_slope

        real(kind=8),intent(inout)::set_guess

    end subroutine readallitemspecs
    subroutine readitemspecs(item_name,int_dim,num_choices,num_cat,num_cat_obs,skill_num,&
        slope_dim,between_item,cat_map,fixguess,guessing,special_item_int,&
        special_item_slope,nommod,set_guess)

        implicit none
        character(len=32)::item_name
        integer,intent(inout)::int_dim,num_cat,&
            num_cat_obs,num_choices,skill_num
        integer,intent(inout)::slope_dim(:)


        logical,intent(inout)::between_item,cat_map,fixguess,guessing,&
            nommod,special_item_int,special_item_slope

        real(kind=8),intent(inout)::set_guess

    end subroutine readitemspecs
end interface
character(len=32),intent(out)::itemname(:)
integer,intent(out)::choices(:),intdim(:),numcat(:),numcatobs(:),skillnumbers(:),slopedim(:,:)
logical,intent(out)::betweenitem(:),catmap(:),constguess,&
	fixedguess(:),guess(:),nominal(:),specialitemint(:),specialitemslope(:)
real(kind=8),intent(out)::setguess(:)
!	buffer is used to define the default item parameter.
character(len=4)::buffer
!	item_name is an item name
character(len=32)::item_name
!	int_dim is the dimension for intercepts.
!	io is the input indicator.
!	item counts items.
!	num_cat is the number of underlying categories per item.
!	num_cat_obs is the number of observed categories per item.
!	num_choices is the number of choices per item.
!	skill_num is the skill number.

integer::int_dim,io,item,num_cat,&
	num_cat_obs,num_choices,&
	skill_num
!	between_item indicates a between item model.
!	cat_map is for category maps.
!	fixguess is for fixed guessing parameters.
!	guessing indicates if the item has a guessing parameter.
!	item_specs indicates that item-specific data are provided.
!	nommod is for a nominal model.
!	special_item_int indicates if a special model is used for an item intercept.
!	special_item_slope indicates if a special model is used for an item slope.
!	slope_dim is the dimensions for slopes.
integer,allocatable::slope_dim(:)
logical::between_item,cat_map,fixguess,guessing,item_specs,nommod,special_item_int,special_item_slope
!	set_guess is the standard guessing parameter.
real(kind=8)::set_guess
!	skillnumbers is used for simple structure specification.
!	Default settings for itemspec.
allocate(slope_dim(size(slopedim,1)),stat=io)
if(io/=0)  stop "Allocation of slope dimension specification not successful."
int_dim=-1
num_choices=0
num_cat=2
num_cat_obs=2
slope_dim=1
between_item=.true.
cat_map=.false.
constguess=.true.
fixguess=.false.
guessing=.false.
item_specs=.false.
nommod=.false.
special_item_int=.false.
special_item_slope=.false.
set_guess=-1.0_8
!	Read itemspec group.
call readallitemspecs(int_dim,num_choices,num_cat,num_cat_obs,&
    slope_dim,between_item,cat_map,constguess,fixguess,guessing,&
    item_specs,nommod,special_item_int,special_item_slope,set_guess)

!	Make revisions.
num_cat_obs=max(2,num_cat_obs)
num_cat=max(num_cat_obs,num_cat)
if(num_cat==num_cat_obs)cat_map=.false.

If(num_cat_obs>2)guessing=.false.
if(int_dim<0)int_dim=num_cat-1
if(int_dim>=num_cat)int_dim=num_cat-1
if(int_dim==0)special_item_int=.false.
slope_dim=max(0,slope_dim)
if(nommod)slope_dim=num_cat-1
slope_dim=min(num_cat-1,slope_dim)
if(size(slopedim,1)>1.and.between_item)then
	item_specs=.true.
end if

if(guessing.and.num_cat_obs<3)then
	int_dim=2
	num_cat=4
	cat_map=.false.
	if(set_guess==-1.0_8.and.num_choices>1) set_guess=-log(num_choices-1.0_8)	
end if
if(.not.guessing)fixguess=.false.
if(sum(slope_dim)==0)special_item_slope=.false.
if(special_item_slope.or.special_item_int.or.cat_map) item_specs=.true.
choices=num_choices
intdim=int_dim

numcat=num_cat
numcatobs=num_cat_obs
skillnumbers=1


betweenitem=between_item
catmap=cat_map

fixedguess=fixguess
guess=guessing
nominal=nommod
specialitemint=special_item_int
specialitemslope=special_item_slope
setguess=set_guess

do item=1,size(numcat)
	slopedim(:,item)=slope_dim
	write(buffer,'(i4)') item
	itemname(item)=trim("Item"//adjustl(buffer))
end do

if(item_specs)then
	do item=1,size(numcat)
		item_name=itemname(item)
		between_item=betweenitem(item)
		cat_map=catmap(item)
		num_choices=choices(item)
		if(intdim(item)==numcat(item)-1)then
			int_dim=-1
		else
			int_dim=intdim(item)
		end if
		num_cat=numcat(item)
		num_cat_obs=numcatobs(item)
		skill_num=skillnumbers(item)
		slope_dim=slopedim(:,item)
		fixguess=fixedguess(item)
		guessing=guess(item)
		nommod=nominal(item)
		special_item_int=specialitemint(item)
		special_item_slope=specialitemslope(item)
		set_guess=setguess(item)
        call readitemspecs(item_name,int_dim,num_choices,num_cat,num_cat_obs,skill_num,&
            slope_dim,between_item,cat_map,fixguess,guessing,special_item_int,&
            special_item_slope,nommod,set_guess)


		num_cat_obs=max(1,num_cat_obs)
		num_cat=max(num_cat_obs,num_cat)
		if(num_cat_obs/=2)guessing=.false.
		if(num_cat==num_cat_obs) cat_map=.false.
		if(int_dim<0)int_dim=num_cat-1
		if(int_dim>=num_cat)int_dim=num_cat-1
		if(int_dim==0)special_item_int=.false.
		if(skill_num<1)skill_num=1
		if(skill_num>size(slopedim,1))skill_num=size(slopedim,1)
		slope_dim=max(0,slope_dim)
        if(num_cat==1) nommod=.true.
		if(nommod)slope_dim=num_cat-1
		slope_dim=min(num_cat-1,slope_dim)
		if(nommod)guessing=.false.
		if(guessing.and.num_cat_obs==2)then
			int_dim=2
			num_cat=4
			cat_map=.false.
			slope_dim=min(1,slope_dim)
			if(set_guess==-1.0_8.and.num_choices>0) set_guess=-log(num_choices-1.0_8)
			if(choices(item)>1.and.setguess(item)==-log(choices(item)-1.0_8).and.num_choices>1)	&
				set_guess=-log(num_choices-1.0_8)
		else
			fixguess=.false.
			guessing=.false.
		end if
			
		if(sum(slope_dim)==0)special_item_slope=.false.
		itemname(item)=item_name
		choices(item)=num_choices
		intdim(item)=int_dim
		numcat(item)=num_cat
		numcatobs(item)=num_cat_obs
		skillnumbers(item)=skill_num
		if(between_item)then
			slopedim(:,item)=0
			slopedim(skill_num,item)=slope_dim(skill_num)
		else
			slopedim(:,item)=slope_dim
		end if
		betweenitem(item)=between_item
		catmap(item)=cat_map
		fixedguess(item)=fixguess
		guess(item)=guessing
		nominal(item)=nommod
		specialitemint(item)=special_item_int
		specialitemslope(item)=special_item_slope
		setguess(item)=set_guess
		
	end do
end if
if(any(fixedguess).or.(.not.any(guess)))constguess=.false.
return
end subroutine getitemdata
subroutine readallitemspecs(int_dim,num_choices,num_cat,num_cat_obs,&
    slope_dim,between_item,cat_map,constguess,fixguess,guessing,&
    item_specs,nommod,special_item_int,special_item_slope,set_guess)
    implicit none
    integer,intent(inout)::int_dim,num_cat,&
        num_cat_obs,num_choices
    integer,intent(inout)::slope_dim(:)


    logical,intent(inout)::between_item,cat_map,constguess,fixguess,guessing,item_specs,&
        nommod,special_item_int,special_item_slope

    real(kind=8),intent(inout)::set_guess
    namelist/allitemspecs/int_dim,num_choices,num_cat,num_cat_obs,&
        slope_dim,between_item,cat_map,constguess,fixguess,guessing,&
        item_specs,nommod,special_item_int,special_item_slope,set_guess

    integer::io
    read(*,nml=allitemspecs,iostat=io)
    if(io/=0) stop "Item specifications not successfully read."
    return
end subroutine readallitemspecs
subroutine readitemspecs(item_name,int_dim,num_choices,num_cat,num_cat_obs,skill_num,&
    slope_dim,between_item,cat_map,fixguess,guessing,special_item_int,&
    special_item_slope,nommod,set_guess)
    implicit none
    character(len=32)::item_name
    integer,intent(inout)::int_dim,num_cat,&
        num_cat_obs,num_choices,skill_num
    integer,intent(inout)::slope_dim(:)


    logical,intent(inout)::between_item,cat_map,fixguess,guessing,&
        nommod,special_item_int,special_item_slope

    real(kind=8),intent(inout)::set_guess

    integer::io
    namelist/itemspecs/item_name,int_dim,num_choices,num_cat,num_cat_obs,skill_num,&
        slope_dim,between_item,cat_map,fixguess,guessing,special_item_int,&
        special_item_slope,nommod,set_guess

    read(*,nml=itemspecs,iostat=io)
    if(io/=0) stop "Item specifications not successfully read."
    return
end subroutine readitemspecs
