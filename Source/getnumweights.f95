!	Get number of weights for distributions.
!	slopedim gives the dimension per item.
!	The default dimension is the dimension of the transformed latent vector.
integer function getnumweights(slopedim)
implicit none
integer,intent(in)::slopedim(:,:)
!	io is error flag for input and numweights is the number of weights.
integer::io,numweights
namelist/numberweights/numweights
numweights=size(slopedim,1)
read(*,nml=numberweights,iostat=io)
if(io/=0)stop "Number of weights not read successfully."
getnumweights=max(0,numweights)
return
end function getnumweights
