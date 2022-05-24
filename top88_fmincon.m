%%%% Compliance topology optimization using FMINCON/L-BFGS
%%%% Based on:
%%%%
%%%% AN 88 LINE TOPOLOGY OPTIMIZATION CODE Nov, 2010 %%%%
%%%%
function top88_fmincon(nelx,nely,volfrac,penal,rmin)
%% MATERIAL PROPERTIES
E0 = 1;
Emin = 1e-9;
nu = 0.3;
%% PREPARE FINITE ELEMENT ANALYSIS
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
F = sparse(2,1,-1,2*(nely+1)*(nelx+1),1);
fixeddofs = union([1:2:2*(nely+1)],[2*(nelx+1)*(nely+1)]);
alldofs = [1:2*(nely+1)*(nelx+1)];
freedofs = setdiff(alldofs,fixeddofs);
%% PREPARE FILTER
iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for i1 = 1:nelx
  for j1 = 1:nely
    e1 = (i1-1)*nely+j1;
    for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
      for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
        e2 = (i2-1)*nely+j2;
        k = k+1;
        iH(k) = e1;
        jH(k) = e2;
        sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
      end
    end
  end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);
%% STORE FILTER AS A SINGLE SPARSE MATRIX
HsH = spdiags(1.0./Hs,0,nelx*nely,nelx*nely)*H;
%% PACK DATA IN ONE STRUCTURE FOR USE IN SUBROUTINES
data.KE = KE; data.iK = iK; data.jK = jK; 
data.nelx = nelx; data.nely = nely; data.F = F;
data.freedofs = freedofs; data.edofMat = edofMat;
data.penal = penal; 
data.Emin = Emin; data.E0 = E0; data.volfrac = volfrac;
data.HsH  = HsH;
%% VOLUME CONSTRAINT
A = (data.HsH'*repmat(1/(nelx*nely),nelx*nely,1))';
%% LOWER AND UPPER DESIGN BOUNDS
xmin = zeros(nelx*nely,1);
xmax = ones (nelx*nely,1);
%% INITIAL DESIGN
x0 = volfrac(ones(nelx*nely,1));
%% OPTIMIZATION OPTIONS (help fmincon)
options = optimset('Algorithm', 'interior-point', ...
                   'GradObj','on', ...
                   'DerivativeCheck', 'off', ...
                   'PlotFcn', {@(x,ov,s)PlotX(x,ov,s,data)}, ...                   
                   'Display','iter', ...
                   'Hessian',{'lbfgs',25});
%% CALL THE OPTIMIZATION SOFTWARE
[x,fval,flag,output] = fmincon( ...
    ... % objective function/gradient:
    @(x) fdf(x,data), ...
    ... % initial guess:
    x0,...
    ... % linear inequality constraints: A*x <= volfrac
    A, volfrac, ...
    ... % linear equality constraints: none
    [],[], ...
    ... % lower/upper bounds
    xmin, xmax, ...
    ... % non-linear constraints: none
    [], ...
    ... % finally, optimization options
    options);
%------------------------------------------------------------------
%% Objective function and its gradient
function [c,dc] = fdf(x,data)
if(nargout>1)
  [c,dc]=FE(x,data);
else
  [c]=FE(x,data);
end
%------------------------------------------------------------------
%% Plot function
function [stop]=PlotX(x,~,~,data)
%% Filtering
xPhys = data.HsH * x; 
%% PLOT DENSITIES
colormap(gray); imagesc(1-reshape(xPhys,data.nely,data.nelx)); 
caxis([0 1]); axis ij; axis equal; axis off; drawnow;
stop=false;
%------------------------------------------------------------------
%% perform FE-analysis
function [c,dc]=FE(x,d)
% store some variables between the calls
persistent U x_old xPhys L s ce
if length(x_old) ~= length(x)
  % pre-allocate memory
  x_old = repmat(NaN,size(x));
  U     = zeros(2*(d.nely+1)*(d.nelx+1),1); 
  A     = zeros(2*(d.nely+1)*(d.nelx+1),1); 
end
%
if any(x_old~=x)
  % need to re-assemble and re-factorize
  x_old = x;
  %% Filtering
  xPhys = d.HsH * x; 
  %% FE-ANALYSIS
  sK = reshape(d.KE(:)*(d.Emin+xPhys(:)'.^d.penal*(d.E0-d.Emin)),...
               64*d.nelx*d.nely,1);
  K = sparse(d.iK,d.jK,sK); K = (K+K')/2;
  % Cholesky factorization
  [L,~,s]=chol(K(d.freedofs,d.freedofs),'lower','vector');
  % Forward/backward substitution
  U(d.freedofs(s))=L'\(L\d.F(d.freedofs(s)));
  %
  ce =  sum((U(d.edofMat)*d.KE).*U(d.edofMat),2);
end
%
% compute outputs
%
%FIXIT: compute objective function (compliance)
c  = 
if nargout > 1
  %FIXIT: compute sensitivities of compliance
  dc = 
  %% MODIFICATION OF SENSITIVITIES
  dc = d.HsH' * dc;
end
