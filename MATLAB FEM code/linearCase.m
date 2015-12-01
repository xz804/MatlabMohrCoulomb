%  This script file is a simple FE code for a Rectangual element with 
%  4,8 or 9 nodes or for a triangular element with 3 or 6 nodes 
%  under uniform loading.the input parameters are the dimensions,
%  material properties and selection of element types.
%  The out put will be the deformation and stress of the region.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clc
tic;   
feval('setpath')

%   GLOBAL VARIABLES
global node element elemType E nu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp([num2str(toc),'    INPUT PARAMETERS '])

E  = 10e3 ;      % modulus of elasticity, KPa
nu = 0.3 ;       % poissons ratio  
L = 5;           % width of domain
D = 5;           % depth of domain
numx=30;          % number of elements in x-direction
numy=30;          % number of elements in y-direction
gamma=0;         % unit weight of soil
sigmatoy=-10;     % vertical load
sigmatox =0;     % horizontal load
load_edge1=2;    % starting x-coordinate point of load, SHOULD BE AT A NODE
load_edge2=3;    % ending x-coordinate point of load, SHOULD BE AT A NODE
elemType ='Q4';
normal_order = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% elastic compliance matrix 
 %C= E/(1+nu)/(1-2*nu)*[ 1-nu nu 0;nu 1-nu 0;0 0 0.5-nu];
 C= E/(1+nu)/(1-2*nu)*[ 1-nu nu 0 nu;nu 1-nu 0 nu;0 0 0.5-nu 0;nu nu 0 1-nu];
% Four corner points
pt1 = [ 0 -D/2];
pt2 = [ L -D/2];
pt3 = [ L  D/2];
pt4 = [ 0  D/2];

disp([num2str(toc),'    MESH GENERATION....'])

 switch elemType
     case {'Q9','Q8','Q4','T3'}
      [node,element] = mesh_region(pt1, pt2, pt3, pt4,numx,numy,elemType);
     otherwise 
      [node,element] =mesh_t6_elem(L,D,numx,numy);
 end
 [topEdge,topEdge1,dispNodes,dispNodes1,leftNodes1]=supportcond(elemType,numx,numy);  

 
disp([num2str(toc),'    STIFFNESS MATRIX COMPUTATION....'])

K=stiffness_matrix(node,element,elemType,normal_order,C);
numnode = size(node,1);
numelem = size(element,1);
total_unknown = numnode*2;selfwt=selfwt_matrix(elemType,normal_order,gamma,node,element);
  switch elemType
     case {'Q9','Q8','T6'}          
            [f,sctry]=force_matrix689(node,topEdge,topEdge1,sigmatoy,sigmatox,load_edge1,load_edge2);
     case{'Q4','T3'}
            [f,sctry]=force_matrix(node,topEdge,sigmatoy,sigmatox,load_edge1,load_edge2);
  end


disp([num2str(toc),'    CALCULATING DISPLACEMENTS'])

[U,u_x,u_y] =displacements(dispNodes,dispNodes1,numnode,K,f,selfwt);


disp([num2str(toc),'    STRESS COMPUTATION....'])
[stress,strain] =stress_calculation(node,element,elemType,U,normal_order,C);
[r] = internalrxn( node,element,elemType,normal_order,stress);
 s=permute(stress,[2 1 3]);  sigma=permute(stress,[3 2 1]);  
% Plot the FEM mesh 
plot_m(elemType,dispNodes,dispNodes1)
title('Undeformed FE mesh')
%Plot numerical deformed configuration
dispnorm=L/max(sqrt(u_x.^2+u_y.^2));
fac =dispnorm*0.05;                    %magnification factor
plot_def(fac,u_x,u_y,elemType,dispNodes,dispNodes1);
%plot stress and deformation intensity with a colormap
plot_defo(fac,u_x,u_y,elemType)
plot_sig(fac,u_x,u_y,elemType,sigma)

% end of the script file


