!	Set intercepts by a regression.    


!	catobsrange maps observed to underlying categories.
!	dat is the array of responses.
!	numcat is the number of underlying categories per item.
!	numcatobs is the number of observed categories per item.
!	design is the design matrix.
!	obsweight is the vector of observation weights.
!	offset is the offset for beta.
!	tolsing is a tolerance for modified Cholesky decomposition.
!	gamma is the vector of starting value.
subroutine interceptstart(catobsrange,dat,numcat,numcatobs,design,obsweight,offset,tolsing,gamma)
implicit none
interface
	function chol(sym,tolsing)
!	chol is used to compute a modified Cholesky decomposition of the n by n matrix sym.
!	The singularity tolerance is tolsing.
		implicit none
		real(kind=8),intent(in)::sym(:,:),tolsing
	 	real(kind=8)::chol(size(sym,1),size(sym,1))
	end function chol
!		This program uses the modified Cholesky decomposition
!		in lu to solve the equation lusolve = b.
	function solve(lu,b)
		implicit none

		real(kind=8),intent(in)::lu(:,:)
		real(kind=8),intent(in)::b(:)
		real(kind=8)::solve(size(b))
		end function solve
end interface
integer,intent(in)::catobsrange(:,:),dat(:,:),numcat(:),numcatobs(:)
real(kind=8),intent(in)::design(:,:),obsweight(:),offset(:),tolsing
real(kind=8),intent(out)::gamma(:)
!	al is the allocation error flag.
!	cat looks at observed categories.
!	colsub is the number of columns in subdesign.
!	obs counts observations.
!	item counts items.
!	position and position1 count observed category groups.
!	sumnumcat is sum(numcat).
integer::al,cat,colsub,item,obs,position,position1,sumnumcat
!	gammas to use.
logical::usegamma(size(gamma))
!	ave is a average.
!	counts are marginal frequencies.
!	logs are adjusted log frequencies.
!	total is a total.
real(kind=8)::ave,counts(sum(numcat)),logs(sum(numcat)),total
!	solution is used for a solution associated with the subdesign.
!	subchol is used for a modified Cholesky decomposition associated with subdesign'subdesign.
!	subdesign is the portion of design used in computations of the gamma starting values for intercepts.
real(kind=8),allocatable::solution(:),subchol(:,:),subdesign(:,:)
counts=0.0_8
sumnumcat=sum(numcat)


!	Count marginal totals.
do obs=1,size(dat,2)
	position=1
	do item=1,size(numcat)
		if(dat(item,obs)>=0.and.dat(item,obs)<numcatobs(item))then
			cat=position+dat(item,obs)
			counts(catobsrange(1,cat):catobsrange(2,cat))=&
				counts(catobsrange(1,cat):catobsrange(2,cat))&
				+obsweight(obs)/(catobsrange(2,cat)-catobsrange(1,cat)+1.0_8)
		end if
		position=position+numcatobs(item)
	end do
end do

!	Obtain adjusted log frequencies.
position=1
do  item=1,size(dat,1)
	
	do cat=position,position+numcat(item)-1
		
		if(counts(cat)>0.0_8) then
			logs(cat)=log(counts(cat))-offset(cat)
		else
			logs(cat)=-log(2.0_8)-offset(cat)
		end if
		
	end do
	total=sum(counts(position:position+numcat(item)-1))
	if(total>0.0_8)then
		ave=dot_product(logs(position:position+numcat(item)-1), &
			counts(position:position+numcat(item)-1))/total
	else
		ave=0.0_8
	end if
	logs(position:position+numcat(item)-1)=logs(position:position+numcat(item)-1)-ave
	position=position+numcat(item)
end do
!	Find a matrix to use if possible.
usegamma=.false.

do position=1,size(gamma)
	if(maxval(abs(design(1:sumnumcat,position)))>0.0_8&
		.and.maxval(abs(design(sumnumcat+1:size(offset),position)))==0.0_8) &
		usegamma(position) = .true.
end do
!	Check for failure.
colsub=count(usegamma)

if(colsub==0)return
!	Set up submatrix.
allocate(subdesign(sumnumcat,colsub),subchol(colsub,colsub),solution(colsub),stat=al)
if(al/=0) stop "Allocation failed for calculation of starting values of intercepts."
position1=1
do position=1,size(gamma)
	if(usegamma(position))then
		subdesign(:,position1)=design(1:sumnumcat,position)
		position1=position1+1
	end if
end do

do position=1,colsub
	position1=1
	do  item=1,size(dat,1)
		total=sum(counts(position1:position1+numcat(item)-1))
		if(total>0.0_8)then
			ave=dot_product(subdesign(position1:position1+numcat(item)-1,position), &
				counts(position1:position1+numcat(item)-1))/total
		else
			ave=0.0_8
		end if
		subdesign(position1:position1+numcat(item)-1,position)=&
			subdesign(position1:position1+numcat(item)-1,position)-ave
		position1=position1+numcat(item)
	end do	
	do position1=1,position
		subchol(position,position1)=sum(counts*subdesign(:,position)*subdesign(:,position1))
	end do
end do

subchol=chol(subchol,tolsing)
do position=1,colsub
	solution(position)=sum(subdesign(:,position)*counts*logs)
end do
solution=solve(subchol,solution)

position1=1
do position=1,size(gamma)
	if(usegamma(position))then
		gamma(position)=solution(position1)
		position1=position1+1
	end if
end do
return
end subroutine interceptstart	
