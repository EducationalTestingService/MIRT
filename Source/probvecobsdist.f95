!	probvecobsdist is used to compute observed distractor probabilities given the latent vector.
!	choices is the number of choices per item.
!	distmap maps distractors to item scores.
!	numcatobs provides the number of observed categories per item.

!	mask indicates which items were presented.
!	distprob provides distractor probabilities.
!	probcatobs is the vector of item probabilities.
function probvecobsdist(choices,distmap,numcatobs,mask,distprob,probcatobs)
implicit none

integer,intent(in)::choices(:),distmap(:),numcatobs(:)
logical,intent(in)::mask(:)
real(kind=8),intent(in)::distprob(:),probcatobs(:)
real(kind=8)::probvecobsdist(sum(choices))

!
!	item is the item number.
	
!	nlocobs is the location of the observed item category.
!	position locates the item response.
integer::item,nlocobs,position
probvecobsdist=0.0_8
nlocobs=1

do item=1,size(mask)
	if(mask(item))then
		do position=nlocobs,nlocobs+choices(item)-1
			if(distmap(position)>=0) probvecobsdist(position)=probcatobs(distmap(position))*distprob(position)
			
			
			
		end do
		
		
	end if
	nlocobs=nlocobs+choices(item)
end do

return
end function probvecobsdist
