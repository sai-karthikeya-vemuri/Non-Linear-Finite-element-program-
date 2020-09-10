%% Generate list of position of nodes according to a geometric series
%    for assignement in "Nonlinear Finite Element Methods" 
%    in summer term 2020
%    lecturer in charge: Dr. Geralf Hütter
%
%Input Parameters
%b=150; %outer radius
%a=100; %inner radius
%nelem=10; %number of elements
function[rnodes]=mesh_refinement(nelem,a,b)
meshrefinementfactor=2; %ratio of element sizes at outer and inner radius
%
%ratio between element sizes of subsequent elements for a geometric series
q=meshrefinementfactor^(1./(nelem-1));
%size of first interval
dr=(b-a)*(1-q)/(1-meshrefinementfactor*q);
rnode=a;
rnodes=[a];
%loop over all elements
for i=1:nelem
        rnode=rnode+dr;
        rnodes=[rnodes;rnode];
        dr=dr*q;
end
end
%visualize location of nodes
%plot(rnodes,zeros(1,nelem+1),'x')
%xlabel('r')
%legend('nodes')
