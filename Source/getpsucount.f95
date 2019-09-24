!	Get data on psu counts per stratum for the nstratum strata.
!	psu is the array of psu codes.
!	stratum is the array of stratum codes.
!	npsu is the array of psu counts.
subroutine getpsucount(psu,stratum,npsu)
implicit none
integer,intent(in)::psu(:),stratum(:)
integer,intent(out)::npsu(:)
!	s counts strata.
integer::s

if(size(npsu)==1)then
	npsu(1)=maxval(psu)
else
	do s=1,size(npsu)
		npsu(s)=maxval(psu,stratum==s)
	end do
end if

return
end subroutine getpsucount
