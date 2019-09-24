!   Get distractor probabilities.
!   dat is response vector.
!   distdat is distractor array.
!   choices gives number of distractors per item.
!   numcatobs is number of observed categories for each item.
!   distprob gives conditional probabilities of distractors.
!   prob gives probabilities of item scores.
function getprobdist(dat,distdat,choices,numcatobs,distprob,prob)
implicit none
integer,intent(in)::dat(:,:),distdat(:,:),choices(:),numcatobs(:)
real(kind=8),intent(in)::distprob(:),prob(:)
real(kind=8)::getprobdist(size(prob))

integer::item,obs,positiondist,resp(size(choices)),respdist(size(choices))
getprobdist=prob
do obs=1,size(prob)
    resp=dat(:,obs)
    respdist=distdat(:,obs)

    positiondist=1
    do item=1,size(choices)
        if(0<=resp(item).and.resp(item)<numcatobs(item))then
            if(distprob(positiondist+respdist(item))<1.0_8)&
                getprobdist(obs)=getprobdist(obs)&
                *distprob(positiondist+respdist(item))
        end if

        positiondist=positiondist+choices(item)
    end do
  
end do
return
end function getprobdist