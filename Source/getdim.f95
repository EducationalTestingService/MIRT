!	Define the latent vector dimension dimlatin and the transformed latent vector dimension dimlatout.,
!	custom is for a special linear transformation.
subroutine getdim(dimlatin,dimlatout,custom)
implicit none
integer::io
integer,intent(out)::dimlatin,dimlatout
logical,intent(out)::custom
namelist/dimension/dimlatin,dimlatout,custom
dimlatin=1
dimlatout=1

custom=.false.
read(*,nml=dimension,iostat=io)
if(io/=0) stop "Dimension or integration data not read successfully."
dimlatin=max(1,dimlatin)
dimlatout=max(1,dimlatout)
return
end subroutine getdim

