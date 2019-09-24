!	Get number of weights for distractor distributions.
!	slopedim gives the dimension per item.
!	The default dimension is the dimension of the transformed latent vector.
integer function getnumweightsdist(slopedim)
implicit none
integer,intent(in)::slopedim(:,:)
!	io is error flag for input and numweights is the number of weights.
integer::io,numweightsdist
namelist/numberweightsdist/numweightsdist
numweightsdist=size(slopedim,1)
read(*,nml=numberweightsdist,iostat=io)
if(io/=0)stop "Number of distractor weights not read successfully."
getnumweightsdist=max(0,numweightsdist)
return
end function getnumweightsdist
