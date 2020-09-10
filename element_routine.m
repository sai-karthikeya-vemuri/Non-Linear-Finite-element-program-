%%Name: Sai Karthikeya Vemuri.
%Immatrikulation Number : 65124


##Element Routine function
##Functionality : Manages shape functions and computes B, Strains,Jacobian using Gauss QUadrature scheme with 1 point
## Input - Radial distance of inner and outer nodes,parameters,nodal displacements, time step, time instance previous overstress
##Output - Elemental internal force , stress, strain and overstress

function [Fint,Ket,sigma,strain,overstress] = element_routine (inner_node, outer_node,params,nodal_displacements,delDisp,time_step,time_instance,prev_overstress)
  J=(outer_node-inner_node)/2;%Calculating the Jacobian
  r=[inner_node,outer_node];
  xhi=0;%quadrature point
  
  N1=0.5*(1-xhi);
  N2=0.5*(1+xhi); %Shape functions
  N=[N1,N2];
  
  B=[-0.5/J,0.5/J;N1/(N1*r(1)+N2*r(2)),N2/(N1*r(1)+N2*r(2))]; %Calculating the Strain Matrix
  strain = B*nodal_displacements;
  delStrain=B*delDisp;
  NR=N1*inner_node+N2*outer_node;
  [sigma,Mat_tan,overstress]=material_routine(params,strain,delStrain,prev_overstress,time_step,time_instance);%Calling the material routine to get stress, Material Tangent Stiffness Matrix
  
  Fint=2*(B'*sigma*J*NR); %Calcuulating the Internal force vector
  
  Ket=2*(B'*Mat_tan*B*J*NR); % Calculating Stiffness Matrix
  
  
endfunction
