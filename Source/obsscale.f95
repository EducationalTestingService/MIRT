!	Find observed scale score for weighted sums.
!	Responses are in dat.
!	maxscore is the maximum raw score.
!	maxw provides maximum item weights.
!	minscore is the minimum raw score.
!	minw provides minimum item weights.
!	numberscales is the number of scale scores.
!	numcatobs provides ranges for observations.
!	weight is the weight.
!	scale is the matrix of scale transformations.





subroutine obsscale(dat,maxscore,maxw,minscore,&
    minw,numberscales,numcatobs,weight,presence,&
	scale,observedscales)
implicit none


integer,intent(in)::dat(:,:),maxscore,maxw(:),minscore,minw(:),numberscales,&
	numcatobs(:),weight(:)
logical,intent(out)::presence(:)

real(kind=8),intent(in)::scale(numberscales,minscore:maxscore)
real(kind=8),intent(out)::observedscales(:,:)

!   loc is a pointer to the right place in weight
!	obs is observation number.

!	row is a row counter.
!   sumraw is a cumulated sum.

integer::loc,obs,row,sumraw



	




presence=.true.





!	Obtain results for each observation.
do obs=1,size(dat,2)

    observedscales(:,obs)=0.0_8
    loc=1
    sumraw=0
    do row=1,size(numcatobs)
        if(maxw(row)==minw(row)) then
            sumraw=sumraw+weight(loc)
        else
            if(dat(row,obs)<0.or.dat(row,obs)>=numcatobs(row)) then
                presence(obs)=.false.

                exit
            else
                sumraw=sumraw+weight(loc+dat(row,obs))
            end if
        end if
        loc=loc+numcatobs(row)
    end do
    if(presence(obs)) then
        observedscales(:,obs)=scale(:,sumraw)
    else
        observedscales(:,obs)=0.0_8
        
    end if
	
end do
return
end subroutine obsscale
