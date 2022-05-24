%%%%%%%%%% FE-ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [U]=FE(nelx,nely,vol,ang,penal)
K = sparse(2*(nelx+1)*(nely+1), 2*(nelx+1)*(nely+1));
F = sparse(2*(nely+1)*(nelx+1),1); U = zeros(2*(nely+1)*(nelx+1),1);
for elx = 1:nelx
    for ely = 1:nely
n1 = (nely+1)*(elx-1)+ely;
n2=(nely+1)*elx +ely;
edof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2]; 
[KE,dKE]=lkOd(ang(ely,elx));
K(edof,edof) = K(edof,edof) + vol(ely,elx)^penal*KE;
K = (K+K')/2;
    end
end
% % DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
%F(2,1) = -1; %DOF 2 NODE 1
%fixeddofs =union([1:2:2*(nely+1)],[2*(nelx+1)*(nely+1)]);
% %DEFINE CANT
%   F(2*(nelx+1)*(nely+1),1) = -1;  %DOF 2*(nelx+1)*(nely+1)the last NODE 1
%   fixeddofs = [1:2*(nely+1)];
% %DEFINE CANT SYM
F(2*((nelx)*(nely+1)+1+nely/2),1) = 1; 
fixeddofs = [1:2*(nely+1)]; %first colum all DOFs blocked
  
  
  
 alldofs = [1:2*(nely+1)*(nelx+1)];
 freedofs = setdiff(alldofs,fixeddofs);


% SOLVING
U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);
U(fixeddofs,:)= 0;


% upper left corner = 1
% bottom left corner = nely+1 
% upper right corner = (nelx)*(nely+1)+1
% bottom right corner = (nelx+1)*(nely+1)
% Finally, considering that nely is even, the middle node of the right edge is
% midpoint = (nelx)*(nely+1)+1+nely/2