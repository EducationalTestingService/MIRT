!   Obtain conditional probabilities for distractors.
!   choices provides the number of choices for each item.
!   data provides item data.
!   distdat provides distractor data.
!   numcatobs provides the number of categories for item scores.
!   obsweight gives observation weights.
!   distfreq provides distractor frequencies.
!   distmap maps distractors to item scores.
!   distprob provides conditional distractor probabilities.
!   scorefreq provides item score frequencies.
subroutine getfreq(choices,dat,distdat,numcatobs,obsweight,&
    distmap,distfreq,distprob,scorefreq)
implicit none
integer,intent(in)::choices(:),dat(:,:),distdat(:,:),numcatobs(:)
real(kind=8),intent(in)::obsweight(:)
integer,intent(out)::distmap(:)
real(kind=8),intent(out)::distfreq(:),distprob(:),scorefreq(:)
!   item counts items.

!   obs counts observations.
!   position locates scored response data.
!   positiondist locates distractor data.
!   resp is a score.
!   respdist is a distractor value.
integer::item,obs,position,positiondist,resp(size(choices)),respdist(size(choices))
distmap=-1
distfreq=0.0_8
scorefreq=0.0_8

positiondist=1
do item=1,size(numcatobs)
    distprob(positiondist:positiondist+choices(item)-1)=1.0_8/choices(item)
    positiondist=positiondist+choices(item)
end do

!   Check observations.
do obs=1,size(obsweight)
    if(obsweight(obs)<=0.0_8)cycle
    resp=dat(:,obs)
    respdist=distdat(:,obs)
    position=1
    positiondist=1
!   Iterate by items.

    do item=1,size(numcatobs)

!   Check for missing responses.
        if(resp(item)>=0.and.resp(item)<numcatobs(item))then
            scorefreq(position+resp(item))=scorefreq(position+resp(item))&
                +obsweight(obs)

            if(respdist(item)>=0.and.respdist(item)<choices(item))then
                distfreq(positiondist+respdist(item))=&
                    distfreq(positiondist+respdist(item))+obsweight(obs)
            else
                stop "Defective distractor input."
            end if
            if(distmap(positiondist+respdist(item))<0)then
                distmap(positiondist+respdist(item))=resp(item)+position
            else
                if(distmap(positiondist+respdist(item))/=resp(item)+position) &
                    stop "Inconsistent combination of item score and distractor."
                
            end if
        end if
        position=position+numcatobs(item)
        positiondist=positiondist+choices(item)
    end do
end do


do positiondist=1,sum(choices)
    if(distmap(positiondist)>=0)then
        if(scorefreq(distmap(positiondist))>0.0_8)distprob(positiondist)=&
            distfreq(positiondist)/scorefreq(distmap(positiondist))
        
    end if
    
end do


return
end subroutine getfreq

