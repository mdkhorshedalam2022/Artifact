function k0=initialize_k_DoubleLayers(kA,kB,orientation)
%Intial values of the elements of hydraulic conductivity tensor 
%for double layers (type) soil 
global P M N 

k0=zeros(P,1);%pre allocation of intial values of the elements of hydraulic conductivity tensor

if strcmp(orientation,'H') % Two soil layer stacked up horiontally 
    k0(:)=kB; %top layer
    k0(1:floor(P/2))=kA; %bottom layer
else % Two soil layer stacked side by side vertically 
    k0(:)=kB; %right layer
        for ii=1:M
          for jj=1:N
              if jj<N/2 
                  k0((ii-1)*N+jj)=kA; %left layer
              end
          end
        end
end