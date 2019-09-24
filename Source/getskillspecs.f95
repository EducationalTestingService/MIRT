!	Obtain skill specifications
!	factorname is the array of factor names.
!	skillname is the array of skill names.
!	rasch is the array of indicators for constant slopes for skills.
!	raschslope1 is the array of indicators for constant slopes of 1 for skills.		
subroutine getskillspecs(factorname,skillname,rasch,raschslope1)
implicit none
character(len=32),intent(in)::factorname(:)
character(len=32),intent(out)::skillname(:)
logical,intent(out)::rasch(:),raschslope1(:)
!	buffer is a character buffer
character(len=4)::buffer
!	skill_name is a skill name
character(len=32)::skill_name
!	io is an input flag.
!	skill counts skills.
integer::io,skill
!	rasch_model is for a Rasch-type model.
!	rasch_slope_1 is for a Rasch-type model with a slope parameter of 1.
!	if repeat_names, skill names are obtained from factor names.
!	skill_specs indicates that data are required on individual skills.
logical::rasch_model,rasch_slope_1,repeat_names,skill_specs
namelist/allskillspecs/rasch_model,rasch_slope_1,repeat_names,skill_specs
namelist/skillspecs/skill_name,rasch_model,rasch_slope_1
!	Default names.
do skill=1,size(skillname)
	write(buffer,'(i4)') skill
	skillname(skill)=trim("Skill"//adjustl(buffer))
end do
!	Defaults
rasch_model=.false.
rasch_slope_1=.false.
repeat_names=.false.
if(size(factorname)==size(skillname))repeat_names=.true.
skill_specs=.false.
read(*,nml=allskillspecs,iostat=io)
if(io/=0)stop "Skill specifications cannot be read."
rasch=rasch_model
raschslope1=rasch_slope_1
if(repeat_names)then
	skill=min(size(factorname),size(skillname))
	skillname(1:skill)=factorname(1:skill)
end if	
!	Individual skills.
if(skill_specs)then
	do skill=1,size(skillname)
		rasch_model=rasch(skill)
		rasch_slope_1=raschslope1(skill)
		skill_name=skillname(skill)
		read(*,nml=skillspecs,iostat=io)
		if(io/=0)stop "Skill specifications cannot be read."
		rasch(skill)=rasch_model
		raschslope1(skill)=rasch_slope_1
		skillname(skill)=skill_name
	end do
end if
!	raschslope1 implies rasch.
rasch=rasch.or.raschslope1
return
end subroutine getskillspecs
