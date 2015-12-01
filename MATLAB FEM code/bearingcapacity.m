%  This script file is a FE code for a Rectangual element with 
%  4,8 or 9 nodes or for a triangular element with 3 or 6 nodes 
%  to do a bearing capacity problem for a Mohr Coulomb soil
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clc
tic; 
feval('setpath')
global node element elemType 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    INPUT PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E=10e3;
nu=0.3;
D=6;
L=5;
phi=5;
tsi=5;
c=1;
numx=30;
numy=40;
gamma=0;
df=-0.1;
load_edge1=0;
load_edge2=1;
elemType = 'Q4' ;
normal_order =2; 
nsteps=80; 
maxit=20;
tol=1e-2;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C= E/(1+nu)/(1-2*nu)*[ 1-nu nu 0 nu;nu 1-nu 0 nu;0 0 0.5-nu 0;nu nu 0 1-nu];
dt=4*(1+nu)*(1-2*nu)/(E*(1-2*nu+(sind(phi))^2));
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
 
ko=1-sind(phi); % k0 according to Jacky's equation
sigmatox=0; % set horizontal loading zero 
numnode = size(node,1);
numelem = size(element,1);
nonelm=size(element,2);
total_unknown = 2*numnode;  
udofs  = [(dispNodes.*2)-1;(dispNodes1.*2)-1]; %prescribed disp.in x-dir
vdofs  = dispNodes.*2;                         %prescribed disp. in y-dir
dofs=union(udofs(:),vdofs(:));                 %overall prescribed disp.
unknowndof=setdiff((1:total_unknown)',dofs);
 
disp([num2str(toc),'    STIFFNESS MATRIX COMPUTATION....'])
K=stiffness_matrix(node,element,elemType,normal_order,C);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   disp([num2str(toc),'    INITIAL STRESS WITH K0 PROCEDURE....'])
selfwt=selfwt_matrix(elemType,normal_order,gamma,node,element);sigmatoy1=0;                         
switch elemType
   case {'Q9','Q8','T6'}          
       [f,sctry]=force_matrix689(node,topEdge,topEdge1,sigmatoy1,sigmatox,load_edge1,load_edge2);
   case{'Q4','T3'}
       [f,sctry]=force_matrix(node,topEdge,sigmatoy1,sigmatox,load_edge1,load_edge2);
end
[U,u_x,u_y] =displacements(dispNodes,dispNodes1,numnode,K,f,selfwt);
[stress,strain] =stress_calculation(node,element,elemType,U,normal_order,C);
 stress(1,:,:)=ko*stress(2,:,:); stress(4,:,:)=ko*stress(2,:,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 switch elemType
     case {'Q9','Q8','Q4'}    
           ule=numelem-numx+1;stp=1;
     case {'T3','T6'} 
           ule=(2*numx*numy)-2*numx+2;stp=3;
 end
strainP=zeros(4,nonelm,numelem); % set parameters to zero
stress_tr=zeros(4,nonelm,numelem);
ui=zeros(total_unknown,1);load=0;
force=zeros(total_unknown,1);
f_old=zeros(total_unknown,1);
r=zeros(total_unknown,1);
b=zeros(total_unknown,1);
du=zeros(total_unknown,1);
dgds=zeros(4,nonelm,numelem);
stressule=stress(:,stp,ule); 
[p,q,theta] = invariants(stressule);
                     % prepare space for plotting data
pvq=zeros(2,nsteps+1);     pvq(1,1)=p(1);    pvq(2,1)=q(1);
evsyy=zeros(2,nsteps+1);   evsyy(1,1)=0;     evsyy(2,1)=stress(2,stp,ule);  
epsvq=zeros(2,nsteps+1);   epsvq(1,1)=0;     epsvq(2,1)=q(1);
fvu=zeros(2,nsteps+1);     fvu(1,1)=0;     fvu(2,1)=0;                              
sigmatoy=0;
                    % start load stepping
for steps=1:nsteps
    stepno=steps;
    err=1; nit=0;  
    sigmatoy=df;
    switch elemType
     case {'Q9','Q8','T6'}          
            [f,sctry]=force_matrix689(node,topEdge,topEdge1,sigmatoy,sigmatox,load_edge1,load_edge2);
     case{'Q4','T3'}
            [f]=force_matrix(node,topEdge,sigmatoy,sigmatox,load_edge1,load_edge2);
    end
    
    r=zeros(total_unknown,1);   % reset unknown variables to zero         
    DEPS_PLA=zeros(4,nonelm,numelem);
    Du=zeros(total_unknown,1);
    du_old=zeros(total_unknown,1);
    Deps=zeros(4,nonelm,numelem);
    dsig_pla=zeros(4,nonelm,numelem);
    Dsig=zeros(4,nonelm,numelem);
    deps_pla=zeros(4,nonelm,numelem);
   
                        % start iteration loop
while (err>tol) && (nit<maxit)
        nit=nit+1;
        du(unknowndof)=K(unknowndof,unknowndof)\f(unknowndof);         
        Du=Du+du;
for iel = 1 : numelem       % start looping on all elements
    sctr = element(iel,:); % element connectivity
    nn   = length(sctr);   % number of nodes per element
    eldof =elementdof(elemType,sctr);  
   [W,Q] = gauss_rule(iel,elemType,normal_order);     
       
    for kk = 1:nn          % start looping on all gauss points
        pt = Q(kk,:);                             
        [N,dNdxi] = shape_func(elemType,pt);  
        J0 = node(sctr,:)'*dNdxi;                    
        Bfem4 =Bmatrix4(pt,elemType,iel);
        Deps(:,kk,iel)=Bfem4*du(eldof);       
        Deps(end,:)=0;
        Dsig(:,kk,iel)= C*(Deps(:,kk,iel)-DEPS_PLA(:,kk,iel)); 
        stress_tr(:,kk,iel)=stress(:,kk,iel)+Dsig(:,kk,iel);              
       
       [p,q,theta] =invariants2(kk,iel,stress_tr);
       F=p*sind(phi)+q*((cosd(theta)/sqrt(3))-(sind(theta)*sind(phi)/3))-c*cosd(phi);
      
         if F<0  
             stress(:,kk,iel)=stress_tr(:,kk,iel);
             err=0;
         else      
             [m1,m2,m3] = formm(kk,iel,stress_tr);
             [dg1,dg2,dg3] =formdg(tsi,q,theta);
             dgds(:,kk,iel)=(dg1*m1+dg2*m2+dg3*m3)*stress_tr(:,kk,iel);
             deps_pla(:,kk,iel)=dt*F*dgds(:,kk,iel);         
             DEPS_PLA(:,kk,iel)=DEPS_PLA(:,kk,iel)+deps_pla(:,kk,iel);  
                            if nit==1
                                err=1;
                            else
                               err=max(abs(du_old(2:2:end)-du(2:2:end)));
                            end   
             
         end       
          r(eldof) =r(eldof)+Bfem4'*(C*deps_pla(:,kk,iel))*W(kk)*det(J0);          
    end                 % end of looping on GPs   
   
          f(sctry)=f(sctry)+r(sctry);
end                     % end of looping on elements 
                   
          du_old=du;f_old=f;
end                   %end of iteration
 
     stress(:,kk,iel)=stress_tr(:,kk,iel);       
     s=permute(stress,[2 1 3]);  sigma=permute(stress,[3 2 1]);   
     ui=ui+Du;
     u_x=ui(1:2:2*numnode-1) ;
     u_y=ui(2:2:2*numnode) ;
     strain=strain+Deps;    
     strainP=strainP+DEPS_PLA;    
     xxx=strain(2,stp,ule);
     yyy=stress(2,stp,ule);               
     eps_vol=strainP(1,stp,ule)+strainP(2,stp,ule)+strainP(4,stp,ule);
     stressule=stress(:,stp,ule); 
     [p,q,theta] = invariants(stressule);
     pvq(1,steps+1)=p(end,:);
     pvq(2,steps+1)=q(end,:);
     epsvq(1,steps+1)=eps_vol(end,:);
     epsvq(2,steps+1)=q(end,:);    
     evsyy(1,steps+1)=xxx(end,:);
     evsyy(2,steps+1)=yyy(end,:);
     load=load+df; 
     fvu(1,steps+1)=u_y(end,:);
     fvu(2,steps+1)=load(end,:);
end                     %end of load step

figure
hold on
plot(abs(pvq(1,:)),pvq(2,:),'--b*','linewidth',2);

  xlabel({'P'},'FontSize',16);
  ylabel({'q'},'FontSize',16);

figure
hold on
plot(abs(epsvq(1,:)),epsvq(2,:),'--r*','linewidth',2);

  xlabel({'epsvol'},'FontSize',16);
  ylabel({'q'},'FontSize',16);

  
  figure
hold on
plot(abs(evsyy(1,:)),abs(evsyy(2,:)),'--r*','linewidth',2);

  xlabel({'epsyy'},'FontSize',16);
  ylabel({'sigyy'},'FontSize',16);
  
  figure
hold on
plot(abs(fvu(1,:)),abs(fvu(2,:)),'--r*','linewidth',2);

  xlabel({'u_y'},'FontSize',16);
  ylabel({'load'},'FontSize',16);
  
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







