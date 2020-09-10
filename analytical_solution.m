

##Name: Sai Karthikeya 
##Immatrikulation Number : 65124
## Analytical solution function 
##Input : parameters and nodes
##Output : analytical solution of displacements for linear-elastic case

function u = analytical_solution (params, nodes)
  u=zeros(length(nodes),1);
  mu=params(2);
  p=params(7);
  E=params(1);
  a=params(5);
  b=params(6);
  %Calculate analytical displacements
  u=(1+mu)*(p/E)*(a^2/(a^2-b^2))*(((1-2*mu).*nodes)+(b^2./nodes));
  

endfunction
