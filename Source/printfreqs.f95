!   Output conditional probabilities for distractors.
!   choices provides the number of choices for each item.

!   numcatobs provides the number of categories for item scores.
!   distfreq provides distractor frequencies.
!   distprob provides conditional distractor probabilities.
!   scorefreq provides item score frequencies.
!	sdfreq gives standard errors for frequencies.
!	sdprob gives standard errors for probabilities.
subroutine printfreqs(itemname,choices,numcatobs,unitfreq,&
    distmap,distfreq,distprob,scorefreq,sdfreq,sdprob)
implicit none
character(len=32),intent(in)::itemname(:)

integer,intent(in)::choices(:),distmap(:),numcatobs(:),unitfreq
        
       
real(kind=8),intent(in)::distfreq(:),distprob(:),scorefreq(:),sdfreq(:),sdprob(:)
character(len=12)::buff1,buff2
character(len=25)::buff(5)

!	cat counts categories.
!	dist counts distractors.
!	io is the flag for output error.
!	item counts items.
!	position is for item scores.
!	positiondist is for distractors.
!	row counts rows.
integer::cat,dist,io,item,position,positiondist,row


write(unit=unitfreq,fmt='(a)',iostat=io) 'Marginal distribution of distractor choices'
if(io/=0) stop "Printing of marginal distribution failed for distractor choices."
write(unit=unitfreq,fmt='(15a)',iostat=io) "Item_name",",","Dist_no",",",&
	"Cat_no",",","Cat_freq",",","Dist_freq",",","SE_dist_freq",",",&
	"Dist_prop",",","SE_dist_prop"

position=1
positiondist=1
do item=1,size(itemname)
	do dist=1,choices(item)
		cat=distmap(positiondist)
		if(cat>=0) then 
			write(buff1,'(i12)') dist-1
		
			
			write(buff2,'(i12)') cat-position
			write(buff(1),'(g25.16e3)') scorefreq(cat)
			write(buff(2),'(g25.16e3)') distfreq(positiondist)
			write(buff(3),'(g25.16e3)') sdfreq(positiondist)
            write(buff(4),'(g25.16e3)') distprob(positiondist)
            write(buff(5),'(g25.16e3)') sdprob(positiondist)
			write(unit=unitfreq,fmt='(15a)',iostat=io)trim(adjustl(itemname(item))),",",trim(adjustl(buff1)),",",trim(adjustl(buff2)),&
				(",",trim(adjustl(buff(row))),row=1,5)
		
			if(io/=0) stop "Printing of marginal distribution failed for distractor choices."
		end if
		positiondist=positiondist+1
	end do
	position=position+numcatobs(item)
end do
return
end subroutine printfreqs
