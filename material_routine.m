%Name: Sai Karthikeya Vemuri.
%Immatrikulation Number : 65124
%Material Routine Function:
%Input :Given parameters,Strain, Overstress of previous time
%Output:Stress, Material Tangent Stiffness, Overstress
%Functionality : Computes Stress and Material Tangent Stiffness Matrix for the given NOn-Linear Visco Elastic case. 

function [stress,Mat_tan,overstress] = material_routine(params,strain,delStrain,prev_overstress,time_step,time_instance)
  E=params(1);
  mu=params(2);
  dt=time_step;
  t=1;
  Q=params(3);
  G=E/((1+mu)*(1-2*mu));
  C=G*[1-mu,mu;mu,1-mu];
  overstress=[0;0];
  dev=[0;0];
  %Hydrostatic and deviatoric parts of strain tensor
  volumetric_strain = (delStrain(1,1)+delStrain(2,1))/3;
  dev(1,1)=delStrain(1,1) -  volumetric_strain;
  dev(2,1)= delStrain(2,1) - volumetric_strain;
  %calculating overstress from previous overstress.
  overstress= (1/(1+(dt/(2*t)))) * ((prev_overstress*(1-(dt/(2*t)))) + Q*dev);
  Mat_tan=[0,0;0,0];
  %calculating the material tangent stiffness using euler modified method.
  Mat_tan(1,1)=C(1,1)+((2/3)*Q/(1+(dt/(2*t))));
  Mat_tan(1,2)=C(1,2)-((1/3)*Q/(1+(dt/(2*t))));
  Mat_tan(2,1)=C(2,1)-((1/3)*Q/(1+(dt/(2*t))));
  Mat_tan(2,2)=C(2,2)+((2/3)*Q/(1+(dt/(2*t))));
  stress= C*strain + overstress;
endfunction
