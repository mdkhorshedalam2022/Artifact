function Nodes(type)
%Define node types based on their distribution
global M N nodetype

for ii=N+2:(M-2)*N+N-1
    nodetype(ii)=0;%'Core';
end

for ii=1
    nodetype(ii)=1; %'bottom left node'
end
for ii=N
    nodetype(ii)=2; %'bottom right node'
end

for ii=(M-1)*N+1
    nodetype(ii)=3;%'top left node'
end
    
for ii=M*N
    nodetype(ii)=4; %'top right node'
end



for ii=2:N-1
    nodetype(ii)=5;%'bottom';
end

for ii= (M-1)*N+2:M*N-1
    nodetype(ii)=6;%'top'; 
end



for ii=N+1:N:(M-2)*N+1
    nodetype(ii)=7;%'left side'; 
end

for ii=2*N:N:(M-2)*N+N
    nodetype(ii)=8;%'right side';
end

if strcmp(type,'A')
    %========Case A Horziontal Gradient: A one-to-three Combiantion======== 
    %left boundary midpoint inlet
    for ii=(floor(M/2)-1)*N+1:N:(floor(M/2))*N+1
       nodetype(ii)=9; %inlet;
    end
    %right boundary midpoint outlet
    for ii=(floor(M/2)-1)*N+N:N:(floor(M/2))*N+N
       nodetype(ii)=10;%oulet;
    end
    %right boundary endpoints outlet
    for ii=[2*N, 3*N,(M-2)*N+N, (M-3)*N+N]
       nodetype(ii)=10;%oulet;
    end
    
elseif strcmp(type,'B') 
    %========Case B Vertical Gradient: A  one-to-two Combiantion=========== 
    % midpoint bottom  inlet
    for ii=(floor(N/2)):(floor(N/2)+1)%[31,41]%
       nodetype(ii)=9; %inlet;
    end
    
    %top boundary midpoint outlet
    for ii=(M-1)*N+(floor(N/2)):(M-1)*N+(floor(N/2)+1)
       nodetype(ii)=10;%oulet;
    end
    %top boundary enddpoints outlet
    for ii=[(M-1)*N+2,(M-1)*N+3,M*N-2,M*N-1]
       nodetype(ii)=10;%oulet;
    end
elseif strcmp(type,'C') 
    %==================Case C: A one-to-one Combiantion====================
    %bottom left corner inlet
    for ii=[1,N+1]
       nodetype(ii)=9; %'inlet';
    end
    %top right corner outlet
    for ii=[M*N-1, M*N]
       nodetype(ii)=10;%'oulet';
    end
    
    
    %top left corner inlet
    for ii=[(M-2)*N+1,(M-1)*N+1]
        nodetype(ii)=9; %inlet;
    end
    %bottom right corner outlet
    for ii=[N-1, N]
        nodetype(ii)=10;%oulet;
    end
else  
    %==========Case D  Concrete Dam: A one-to-two Combiantion==============  
    %upstram lake inlet
    for ii=(((M-1)*N+1)):(((M-1)*N+3)+1)
       nodetype(ii)=9; %inlet;
    end
    %downstream lake outlet
    for ii=(M-1)*N+N:-1:(M-1)*N+N-3
       nodetype(ii)=10;%oulet;
    end
end
