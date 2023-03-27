function [dx_final,dz_final,M,N,P,tfinal]=calc_parameters(m,n, tn)
%Calculate the following parameters:
%horizontal step size: dx
%vertical step size: dz
%number of horizonal nodes: N
%number of vertical nodes: M
%total number of time steps: tn

global L dt T

%forcing no of horizonal and vertical nodes to be even to get symmetry about layers
if mod(n,2)==0
    N=n;
else
    N=n+1;
end

if mod(m,2)==0
    M=m;
else
    M=m+1;
end

P=M*N;%total number of nodes
dx_final=L/(N-1);%recalculate dx to consistent with final N
dz_final=T/(M-1);%recalculate dy to conistent with final M
tfinal=(tn*dt);%total number of time steps
end