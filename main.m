%-------------------------------------------------------------------------
%---------------------Comprehensive Exam Artifact------------------------- 
%-------------------------------------------------------------------------

%++++++++++++++++++++++++++++++++++
% Md Khorshed Alam                +
% mdkhorshedalam@u.boisestate.edu +
% Computing PhD program           +
% Boise State University          +
% Boise, ID                       +
% Date: March.26.2022             +
%++++++++++++++++++++++++++++++++++

%The purpose of this program is solve the following 2D transient seepage
%problem(where X & Z are non principle axes)using finite difference Method
% to simulate the impact non-diagonal hydraulic conductivity tensor.
%==========================================================================
% diF v= -theta'(t),i.e.,
% diF(-k*grad(h))= -mv h'(t) for 0<=x<=L and 0<=z<=T
% IC: h(0)=0 and h(T)=0
% Neumann BC: left side boundary: vx=0; right boundary: vx=0; 
%             top boundary: vz=0; bottom boundary: vz=0.
% Dirichlet BC: inlet(s): h=H1, outlet(s): h=H2.
%==========================================================================

% The output of this program is 
% 1- The polt for the distributed nodes(Fig 1)
% 2- Calculate total head and flow velocity   
%    to generate velocity plot with equipotential lines(Fig 2) 
% 3- Calculate pressure head to generate it's contour plot(Fig 3)

%=============================Program Begin================================
clear; %clear workspace 
close all; % close figures
%=======================Defining Soil Properties===========================
global L M N P T a1 a2 nodetype dt htplot hpplot %declaring global avriables
H1=7; % total head at inlets (m);
H2=3; % total head at outlets (m);
L=2.1; %horizontal length (m)
T=2; %vertical length i.e. thickness of soil (m)
dt=10;%time step size 
tn=100;%total number of time steps 
a1=1; a2=3; %Crank-Nicholson coefficients

fprintf('This program simulate seepage flow through soil\n\n')

prompt = "Please enter D for default soil type otherwise enter ND: ";
soiltype=input(prompt,'s');
valid_inputs = {'D';'ND'};

while ~any(strcmp(soiltype,valid_inputs)) % checking invalid inputs
    f = msgbox(["Invalid Value";"Enter Correct Value"],"Error","error");
    disp("Please Re-Enter")
    prompt = "Please enter D for default soil type otherwise enter ND: ";
    soiltype=input(prompt,'s');
end

if strcmp(soiltype,'D')
    fprintf ('Default soil have the following properties: \n\n')
    disp('Major principal hydraulic conductivity: 5.5e-3 (m/s)')
    disp('Minor principal hydraulic conductivity: 1.5e-3 (m/s)')
    fprintf('Slope of soil stratigraphic plane: 0%s\n\n', char(176))
    k1A=5.5e-3; k1B=k1A; alphaA=0;
    k3A=1.5e-3; k3B=k3A; alphaB=0;
    soil_type=1; % for single type of soil

else
    fprintf('\n                      Define soil(s) properties                           \n\n')
    disp('      Typical values of the hydraulic conductivity of soils               ')
    disp('__________________________________________________________________________')
    disp('Soil Type                                  Hydrauclic Conductivity(m/s)   ')
    disp('__________________________________________________________________________')
    disp('Gravels                                               E+00-E-04           ')
    disp('Sands                                                 E-02-E-03           ')
    disp('Silts                                                 E-01-E-05           ')
    disp('Clays                                                 E-03-E-09           ')
    disp('__________________________________________________________________________')
    fprintf('\nThis program can simulate seepage flow through upto two types of soil\n\n')
    soil_type=input('Enter, 1 for single type of soil, and 2 for two types of soil: ');
    if soil_type==1
        %k1 is major principal orientation of soil stratigraphy
        k1= "Enter the values of the major principal hydraulic conductivity, k1(m/s)= ";
        k1A=input(k1);
        k1B=k1A; %for single soil type 
        %k3 is minor principle orientation of soil stratigraphy
        k3= "Enter the values of the minor principal hydraulic conductivity, k3(m/s)= ";
        k3A=input(k3);
        k3B=k3A; %for single soil type 
        %degree of stratigraphy below horizon(<0) and above horzion (>0)
        deg_strata= "Enter the tilt of soil stratigraphy in degree(below horizon, tilt<0 and above horizon, tilt>0)= ";
        alphaA=input(deg_strata); %degree of stratigraphy
        alphaB=alphaA; %for single soil type 
    else
        disp('Please provide the following details of the two types of soil:')
        disp('How two soil layers stacked up: Horizontally or Vertically?')
        layer_orientation='Enter, H if two soil layer stacked up horiontally, or V if wo soil layer stacked side by side vertically : ';
        orientation=input(layer_orientation,'s');
        valid_inputs = {'H';'V'};
            while ~any(strcmp(orientation,valid_inputs)) % checking invalid inputs
                f = msgbox(["Invalid input";"Give Correct Input"],"Error","error"); %throughing error message for invalid input
                layer_orientation='Enter, H if two soil layer stacked up horiontally, or V if wo soil layer stacked side by side vertically :';
                orientation=input(layer_orientation,'s');
            end
            if strcmp(orientation,'H') %two soil layer stacked up horiontally
    
                %bottom soil type A 
                k1_btm= "Enter the values of the major principal hydraulic conductivity for bottom layer of soil, k1(m/s)= ";
                k1A=input(k1_btm);
                k3_btm= "Enter the values of the minor principal hydraulic conductivity for bottom layer of soil, k3(m/s)= ";
                k3A=input(k3_btm);
                deg_strata_btm= "Enter the tilt of bottom soil layer stratigraphy in degree(below horizon, tilt<0 and above horizon, tilt>0)= ";
                alphaA=input(deg_strata_btm);
    
                %top soil type B
                k1_top= "Enter the values of the major principal hydraulic conductivity for top layer of soil, k1(m/s)= ";
                k1B=input(k1_top);
                k3_top= "Enter the values of the minor principal hydraulic conductivity for top layer of soil, k3(m/s)= ";
                k3B=input(k3_top);
                deg_strata_top= "Enter the tilt of top soil layer stratigraphy in degree(below horizon, tilt<0 and above horizon, tilt>0)= ";
                alphaB=input(deg_strata_top);
    
            else %two soil layer stacked up vertically side by side 
    
                %left pile soil type A
                k1_left= "Enter the values of the major principal hydraulic conductivity for left pile (layer) of soil, k1(m/s)= ";
                k1A=input(k1_left);
                k3_left= "Enter the values of the minor principal hydraulic conductivity for left pile (layer) of soil, k3(m/s)= ";
                k3A=input(k3_left);
                deg_strata_left= "Enter the tilt of left pile soil layer stratigraphy in degree(below horizon, tilt<0 and above horizon, tilt>0)= ";
                alphaA=input(deg_strata_left);
    
                %right pile Soil type B 
                k1_right= "Enter the values of the major principal hydraulic conductivity for right pile (layer) of soil, k1(m/s)= ";
                k1B=input(k1_right);
                k3_right= "Enter the values of the minor principal hydraulic conductivity for right pile (layer) of soil, k3(m/s)= ";
                k3B=input(k3_right);
                deg_strata_top= "Enter the tilt of right pile soil layer stratigraphy in degree(below horizon, tilt<0 and above horizon, tilt>0)= ";
                alphaB=input(deg_strata_top);
            end
    end

end
%===================End of Defining Soil Properties========================

%===========Calculate components of the non-diagonalized tensors===========
[kxxA,kzzA,kxzA,kzxA]=twoD_nondiag_k(k1A,k3A,alphaA); %elemnets of non-diagonalized tensor for soil A
[kxxB,kzzB,kxzB,kzxB]=twoD_nondiag_k(k1B,k3B,alphaB); %elemnets of non-diagonalized tensor for soil B
%=====End of calculation of components of hydraulic conductivity tensor====

%negelcting non-diagonal terms
%kxzA=0;kxzB=0; kzxA=0;kzxB=0;

%=====================Calculate some essential parameters==================
%dx; horizontal step size %dz;vertical step size
%N and M; no of horizonal and vertical nodes respectively (forced to be
%even for symmetry); P=M*N; % total no of nodes
%tn; total number of time step 

zNode="Please enter no of vertical nodes:  ";
m=input(zNode);
xNode= "Please enter no of horizontal nodes:  ";
n=input(xNode);
[dx,dz,M,N,P,tfinal]=calc_parameters(m,n,tn);

while P>=1e+09 % checking maximum no of array size MATLAB can handle
    f = msgbox("Requested array exceeds maximum array size preference","Error","error");
    disp("Please Re-Enter less number of nodes")
    zNode="Enter no of vertical nodes:  ";
    m=input(zNode);
    xNode= "Enter no of horizontal nodes:  ";
    n=input(xNode);
    [dx,dz,M,N,P,tfinal]=calc_parameters(m,n,tn);
end
%======================End of calculation of parameters====================

%=========================Pre allocation of variables =====================
a=zeros(P,P); %coefficient matrix
b=zeros(P,1); %constant matrix
mv=zeros(P,1);%retention of water
x=zeros(1,P); z=zeros(1,P); %coordinate of nodes
he=zeros(P,1); %elevation head
h=zeros(P,1); %total head 
h_old=zeros(P,1);%old values of total head
hplot=zeros(M,N);%total head in 2D
htplot=zeros(M,N,tn); %total head with time step
hpplot=zeros(M,N,tn);% pressure head with time step
vx=zeros(M,N); vz=zeros(M,N);%velociti components in 2D

%elements of hyraulic conductivity tensor:
kxx=zeros(P,1); kzz=zeros(P,1); kxz=zeros(P,1); kzx=zeros(P,1);
%old values of the elements of hyraulic conductivity tensor:
kxxold=zeros(P,1); kzzold=zeros(P,1); kxzold=zeros(P,1); kzxold=zeros(P,1); 
%elements of hyraulic conductivity tensor in 2D:
Kxx=zeros(M,N); Kzz=zeros(M,N); Kxz=zeros(M,N); Kzx=zeros(M,N);
%=========================End of allocation of variables ==================

%===============================Initialization ============================
%Calculate intial values of the elements of hydraulic conductivity tensor 
%based on type(s)of soil and orientation of soil layer(s)

if soil_type==1 %for single type of soil
    kxx0=initialize_k(kxxA,kxxB);
    kzz0=initialize_k(kzzA,kzzB);
    kxz0=initialize_k(kxzA,kxzB);
    kzx0=initialize_k(kzxA,kzxB);
else %for two types of soil
    kxx0=initialize_k_DoubleLayers(kxxA,kxxB,orientation);
    kzz0=initialize_k_DoubleLayers(kzzA,kzzB,orientation);
    kxz0=initialize_k_DoubleLayers(kxzA,kxzB,orientation);
    kzx0=initialize_k_DoubleLayers(kzxA,kzxB,orientation);
end

%Intailization of elevation head (he), pressure head(hp) and total head(h)
for ii=1:P
    [x(ii),z(ii)]=nodes_coordinate(dx,dz,ii); % calculate coordinate of each node
    he(ii)=z(ii); %elevation head (m)at each node initially
    h(ii)=0; %total head (m) at each node initially
end

h_old(:,:)=h; % assiging initial value of h as old value of h
%============================End of Intializations=========================

%==========Taking input for types of inlet and outlet combinations=========
disp("What type of inlet and oulet combinations you have on the boundary walls?")
disp(" type A: Horizontal Gradient(Midpoint left inlet, Midpoint right outlet & Endpoints right outlets) ")
disp(" type B: Vertical Gradient(Midpoint bottom  inlet, Midpoint top outlet & Endpoints top outlet)  ")
disp(" type C: One-to-one Combiantion (Bottom left inlet, Top right outlet & Top left inlet, Bottom right outlet)")
disp(" type D: Concrete Dam (Upstram lake inlet & Downstream lake outlet)")
prompt = "Please enter A for type A and so on: ";
type=input(prompt,'s');
valid_inputs = {'A';'B';'C'; 'D'};
while ~any(strcmp(type,valid_inputs)) % checking invalid inputs
    f = msgbox(["Invalid Value";"Enter Correct Value"],"Error","error");
    prompt = "Please enter A for type A and so on: ";
    type=input(prompt,'s'); 
end
%==============================End of taking input=========================

   
%==============================Plot type of nodes==========================
    Nodes(type); %function which specified the type of nodes and boundaries
    f1=figure;
    hold on
    plotNodes(x,z);%plot nodes distribution
%=============================End of node plot=============================
    
%===============================Solver Begin===============================
    
    
%==============================Time Steps Begin============================
for tt=1:tn
       %calculate components of hydraulic conductivity tensor 
       %based on h from previous time step
        [kxx,kzz,kxz,kzx]=calculate_k(kxx0,kzz0,kxz0,kzx0,h_old,h,he);

        %=======================For Loop Begin to Update k ================
        switchk='diverge'; %keep iterate if the elements of k not converged
        
        for counter=1:1000 %limit the loop to run maximum 1000 times
          
           %=============================Loop for node 1 to P==============
           for ii=1:P 
               
              switch nodetype(ii) %assign to appropriate case based on type of nodes 
    
                %================Calculate core equations==================
                  case 0% Core
                      if h(ii)-he(ii)>=0 %for saturated soil
                          mv(ii)=0.00001;
                      else
                          mv(ii)=0.001;%for unsaturated soil
                      end
                      [a(ii,ii-1),a(ii,ii),a(ii,ii+1),a(ii,ii-N),a(ii,ii+N),a(ii,ii+N+1),a(ii,ii-N+1),a(ii,ii+N-1),a(ii,ii-N-1),b(ii)]=coreEqns(kxx(ii),kxx(ii+1),kzz(ii),kzz(ii+N),kxz(ii),kxz(ii+1),kzx(ii),kzx(ii+N),dx,dz,dt,mv(ii,1),h(ii));                
                %=================End of core equations calculation========
    
    
                %================Calculate Boundary conditions(BC's)=======
    
                %**************Calculate BC's for four corner nodes********
                 case 1 %bottom left node(Nueman BC, vz=0)
                      a(ii,:)=0;
                      a(ii,ii+1)=(kzx(ii)/dx);
                      a(ii,ii+N)=(kzz(ii)/dz);
                      a(ii,ii)=-(a(ii,ii+1)+a(ii,ii+N));
                      b(ii)=0; 
                  case 2 %bottom right node(Nueman BC, vz=0)
                      a(ii,:)=0;
                      a(ii,ii-1)=-(kzx(ii)/dx);
                      a(ii,ii+N)=(kzz(ii)/dz);
                      a(ii,ii)=-(a(ii,ii-1)+a(ii,ii+N));
                      b(ii)=0;
                  case 3 %top left node(Nueman BC, vz=0)
                      a(ii,:)=0;
                      a(ii,ii+1)=(kzx(ii)/dx);
                      a(ii,ii-N)=-(kzz(ii)/dz);
                      a(ii,ii)=-(a(ii,ii+1)+a(ii,ii-N));
                      b(ii)=0;
                  case 4 %top right node(Nueman BC, vz=0)
                      a(ii,:)=0;
                      a(ii,ii-1)=-(kzx(ii)/dx);
                      a(ii,ii-N)=-(kzz(ii)/dz);
                      a(ii,ii)=-(a(ii,ii-1)+a(ii,ii-N));
                      b(ii)=0;
                %**********************************************************               
    
                %****************Calculate BC's for four sides*************
                 case 5 %bottom(Nueman BC, vy=0)
                      a(ii,:)=0;
                      a(ii,ii+1)=(kzx(ii)/(2*dx));
                      a(ii,ii-1)=-(kzx(ii)/(2*dx));
                      a(ii,ii+N)=(kzz(ii)/dz);
                      a(ii,ii)=-(kzz(ii)/dz);
                      b(ii)=0;
                  case 6 %top(Nueman BC, vy=0)
                      a(ii,:)=0;
                      a(ii,ii+1)=(kzx(ii)/(2*dx));
                      a(ii,ii-1)=-(kzx(ii)/(2*dx));
                      a(ii,ii)=(kzz(ii)/dz);
                      a(ii,ii-N)=-(kzz(ii)/dz);
                      b(ii)=0;
                  case 7 %left side(Nueman BC, vx=0)
                      a(ii,:)=0;
                      a(ii,ii+1)=(kxx(ii)/dx);
                      a(ii,ii)=-(kxx(ii)/dx);
                      a(ii,ii+N)=(kxz(ii)/(2*dz));
                      a(ii,ii-N)=-(kxz(ii)/(2*dz));
                      b(ii)=0;
                  case 8 %right side(Nueman BC, vx=0)
                      a(ii,:)=0;
                      a(ii,ii)=(kxx(ii)/dx);
                      a(ii,ii-1)=-(kxx(ii)/dx);
                      a(ii,ii+N)=(kxz(ii)/(2*dz));
                      a(ii,ii-N)=-(kxz(ii)/(2*dz));
                      b(ii)=0;
                  %********************************************************                
    
                  %***********Calculate BC's for inlet and outlet**********
                  
                  case 9 %inlet(Dirichlet BC)
                      a(ii,:)=0;
                      a(ii,ii)=1;
                      b(ii)=H1; 
                  case 10 %outlet(Dirichlet BC)
                      a(ii,:)=0;
                      a(ii,ii)=1;
                      b(ii)=H2;  
                  %********************************************************
              end %end of swtich case
              %===================End of Calculation of BC's===============
              
           end %end of for loop
           %========================End of for Loop for nodes==============
         
            h=(a^-1)*b; %solve for total head h

            %storing previous values of the elements of hydrauclic conductivity tensor(k) as old k 
            %before updating k
            [kxxold, kzzold, kxzold, kzxold]=reassign_k(kxx,kzz,kxz,kzx); 
            
            %==============Upadate hydrauclic conductivity tensor (k)======
            %update elements of k with Crank-Nicholson scheme for
            %calculated values of total head h
             for ii=1:P
                 [mv(ii),kxx(ii),kzz(ii),kxz(ii),kzx(ii)]=update_k(h(ii),he(ii),h_old(ii),kxx0(ii),kzz0(ii),kxz0(ii),kzx0(ii));
             end
            %=========================End of Upadate of k==================
             
            %forcing the count loop to keep iterate until values of the
            %elements of k converge to relative error <0.0001
             if max(abs(kxx-kxxold)./kxx0)<= 0.0001 && max(abs(kzz-kzzold)./kxx0)<= 0.0001 && max(abs(kxz-kxzold)./kxx0)<= 0.0001 && max(abs(kzx-kzxold)./kxx0)<= 0.0001
                swithck='converge';
                break % terminate loop if elements of k converge
             end

        end 
        %========================End of For Loop for k=====================

        h_old=h; % reassign value of h as old value of h 
        
%===============================End of Solver==============================
    
    %==========================Calculation of velocity=====================
    %mapping heads and elements of hydarulic conductivity tersor on a grid 
    [hplot,htplot,hpplot,Kxx,Kzz,Kxz,Kzx]=mapping(h,he,kxx,kzz,kxz,kzx,tt);
    %calculate velocity components
    [vx,vz]= calculate_velocity(hplot,Kxx,Kzz,Kxz,Kzx,dx,dz);
    %=========================End of velocity calculation==================
    
    %**********************************Plotting****************************
     figure(2)
     plotQuiver(htplot(:,:,tt),tt,vx,vz); %flow velocity plot with equipotential lines
 end
 %===============================End of Time Loop==========================
    figure(3)
    for tt=1:tn
     plotContour(hpplot(:,:,tt),tt);%plot contour lines of pressure head;
    end
    %******************************End of plotting*************************
      
  %================================Program End=============================