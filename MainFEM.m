 %Non-Linear Finite Element Methods- Assignment 2020 Examiner: Dr.Geralf Hutter
%This program solves the problem of creep of thick walled pipe with internal pressure with assumed plane strain conditions
%Name: Sai Karthikeya Vemuri.
%Immatrikulation Number : 65124

%##### Main Function: Execute this file on Octave.
%Functionality:
%This Script defines the parameters of the problem such as the number of elements , the inner and outer radius of the cylinder, pressure, loading and final times etc.
%It manages and implements essential boundary conditions, implements Newton-Raphson method to get displacements for every time step.
%In the Newton-Raphson scheme , Global stiffness Matrix, Internal and External forces are assembled by calling element routine for each element, the system is solved for displacements until error minimum criteria is met  
%Calculates and Plots stresses, widening history of nodes, displacements of nodes ,overstresses.

######%
clc;
n_el=10;
time_step=0.1;
params=[200000, 0.30, 100000 ,1, 30, 60, 140, 2, 10];%Given parameters
pressure=params(7);
t=[time_step:time_step:params(9)];
nodes=mesh_refinement(n_el,params(5),params(6));%Calling Mesh Refinement Function to get element ends
F_Int=zeros(n_el+1,1);
K_Global=zeros(n_el+1); %Initialization of Stiffness and Internal forces Global Matrices
%initialization of stresses and displacements
radial_stress=zeros(n_el,1);
phi_stress=zeros(n_el,1);
radial_overstress=zeros(n_el,1);
phi_overstress=zeros(n_el,1);
delDisp=zeros(n_el+1,1);
displacement=zeros(n_el+1,1);
load_scaling=(1/params(8));%Load scaling parameter
widening_history=zeros(length(t),1);
counter=1;
%Iterating over time from start to end time
for k=t
  
  iter=0;
if k<params(8)
  %Load Scaling for time less than loading time
  F_Ext=zeros(n_el+1,1);
  F_Ext(1,1)=pressure*params(5)*k*load_scaling;
  
else
  F_Ext=zeros(n_el+1,1);
  F_Ext(1,1)=pressure*params(5);
 end
 %Newton-Raphson Method implementation
while iter <100
  F_Int=0*F_Int;
  K_Global=0*K_Global;
  radial_stress=0*radial_stress;
  phi_stress=0*phi_stress;
for i=1:n_el
  inner=nodes(i);
  outer=nodes(i+1);
  disp1=displacement(i);
  disp2=displacement(i+1);
  nodal_displacements=[disp1;disp2];
  val1=delDisp(i);
  val2=delDisp(i+1);
  nodal_delDisp=[val1;val2];
  temp1=radial_overstress(i);
  temp2=phi_overstress(i);
  prev_overstress=[temp1;temp2];
  %Calling Element routine for the given element
  [Fe_Int,Ket,stress,strain,overstress_element]=element_routine(inner,outer,params,nodal_displacements,nodal_delDisp,time_step,k,prev_overstress);
  %Assembling K_Global and F_Int
  F_Int(i,1)=F_Int(i,1)+Fe_Int(1);
  F_Int(i+1,1)=F_Int(i+1,1)+Fe_Int(2); 
  K_Global(i,i)=K_Global(i,i)+Ket(1,1);
  K_Global(i,i+1)=K_Global(i,i+1)+Ket(1,2);
  K_Global(i+1,i)=K_Global(i+1,i)+Ket(2,1);
  K_Global(i+1,i+1)=K_Global(i+1,i+1)+Ket(2,2);
  radial_stress(i,1)=stress(1,1);
  phi_stress(i,1)=stress(2,1);
  radial_overstress(i,1)=overstress_element(1,1);
  phi_overstress(i,1)=overstress_element(2,1);
endfor
  %Solving the system of equations
  delDisp = linsolve(K_Global,F_Ext-F_Int); 
  displacement = displacement +delDisp;

 %Exit condition for Newton-Raphson
  if ((max(abs(F_Int-F_Ext)))<(.005*(max(abs(F_Int)))) &&(max(abs(delDisp)))<(.005*max(abs(displacement))))
    
    fprintf("converged for time:");
    
    
    disp(k);
    fprintf("converged after iteration:");
    disp(iter);
    break
    
  else
    %fprintf("did'nt converge at:");

    %disp(iter);
    
    
    iter=iter+1;
  end
endwhile
widening_history(counter,1)=displacement(end,1);
counter=counter+1;


endfor

u=analytical_solution(params,nodes);
displacement
radial_stress
phi_stress
radial_overstress
phi_overstress
%plotting the necessary data
%plot(nodes,displacement,"-*","LineWidth",2);
plot(nodes(1:end-1),radial_stress,"m-s","LineWidth",2.5);
%plot(nodes(1:end-1),phi_stress);
%plot(t,widening_history,"-s","LineWidth",2.5);
legend("Stress(r-r) ")
title("Radial Stress after final time");
xlabel("Radius(mm)");
ylabel("Radial Stress(MPa)");





