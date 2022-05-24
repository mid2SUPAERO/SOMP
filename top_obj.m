function [c, dt]=top_obj(x)
global nelx nely vol volfrac ang angle  penal rmin

% define design variables
vol = x(1:length(x)/2); 
ang = x((length(x)/2+1):end); 

vol = reshape(vol,nely,nelx); 
ang = reshape(ang,nely,nelx);


% FE-ANALYSIS
[U]=FE(nelx,nely,vol,ang,penal);
% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
c=0.;
for ely = 1:nely
    for elx = 1:nelx
n1 = (nely+1)*(elx-1)+ely;
n2=(nely+1)*elx +ely;
Ue = U([2*n1-1;2*n1; 2*n2-1;2*n2; 2*n2+1;2*n2+2; 2*n1+1;2*n1+2],1); 
[KE,dKE] =lkOd(ang(ely,elx)); 
c = c + vol(ely,elx)^penal*Ue'*KE*Ue;
dc(ely,elx) = -penal*vol(ely,elx)^(penal-1)*Ue'*KE*Ue;
dca(ely,elx) = -penal*vol(ely,elx)*Ue'*dKE*Ue;
    end
end
% FILTERING OF SENSITIVITIES

 [dc] =check(nelx,nely,rmin,vol,dc);
 [dca] =check(nelx,nely,rmin,ang,dca);
  dt=[dc(:); dca(:)];

