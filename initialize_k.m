function k0=initialize_k(kA,kB)
%Intial values of the elements of hydraulic conductivity tensor 
%for single layer (type) soil with horizontal stratigraphic planes

global P 

k0=zeros(P,1); %pre allocation of intial values of the element of hydraulic conductivity tensor 

k0(:)=kB; k0(1:floor(P/2))=kA; %horizontal stratigraphic planes of soil(assumed)

end