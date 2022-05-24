%%%% A 99 LINE TOPOLOGY OPTIMIZATION CODE BY OLE SIGMUND, JANUARY 2000 %%%
%%%% A 142 LINE CODE MODIFIED FOR USE WITH ORTHOTROPIC MATERIALS %%% 
%function TopO(nelx,nely,volfrac,angle,penal,rmin);
%The results IS HIGHLY SENSITIVE TO X0
% INITIALIZE
clear all; close all;
global nelx nely vol volfrac ang angle  penal rmin

%to use as function remove % and clear all; close all; add end at the end
nelx=30; %number of quad elements in X
nely=20; % number of quad elements in Y
volfrac=0.35; %constraint in volume
angle=0.01; % to avoid 0 fiber aligned with 0 rad
penal=3; %Penalization
rmin=1.5; % filter
filt=3; % convolution filter on orientation
% better results when fibers are driven with limited orientation design bounds
UP=0.33; % max bounds on fiber orientation [-pi/3, pi/3]


% option for fmincon
option = optimoptions('fmincon','Algorithm','interior-point',...
    'GradObj','on',...
'TolX',1E-10,...
'TolFun',1E-10,...
'PlotFcns',@optimplotfval,...
	'GradConstr','on',...
	'MaxFunctionEvaluations',2000);
%
c0 = volfrac;
rho0 = volfrac*ones(nely,nelx); 
%theta0 = angle*ones(nely,nelx); % could be default


% works well for Cantilever Sym Beam
thetau=-1*pi/4*ones(round(nely/2),nelx);
thetas=1*pi/4*ones(round(nely/2),nelx);
 
 
theta0=[thetau;thetas];
nele=nelx*nely;
% define phi (rotation about y axis) 
% offset = 0;
% a = -offset/180*pi; b = offset/180*pi;
% phi = a + (b-a)*rand(nele,1); 
% phi = reshape(phi,nely,nelx); %


%The results IS HIGHLY SENSITIVE TO X0

x0 = [rho0(:);theta0(:)];
lb = [1E-6*ones(length(rho0(:)),1);-pi*UP*ones(length(theta0(:)),1)];
ub = [ones(length(rho0(:)),1);pi*UP*ones(length(theta0(:)),1)];
%
% equality constraint
Aeq = [ones(1,length(rho0(:))) zeros(1,length(theta0(:)))]; 
beq = nelx*nely*c0;
% fmincon % add @top_nonlcon (only for passive element) 
x = fmincon('top_obj',x0,[],[],Aeq,beq,lb,ub,[],option); 


%
rho = x(1:length(x)/2); theta = x(length(x)/2+1:end); 
rho = reshape(rho,nely,nelx); 
theta = reshape(theta,nely,nelx);
% ielem1 = ielem((find(rho>.5)),:); rho1 = rho(find(rho>.5),:); theta1= theta(find(rho>.5),:);
%
[XC,YC]=meshgrid(1:nelx,1:nely);
rho_plot = 1-rho(:,:); % for imagesc
%
[XC,YC]=meshgrid(1:nelx,1:nely);

% % PLOT DENSITIES and ORIENTATION as HEATMAP
figure(2);subplot(211);
heatmap(rho);
subplot(212);
heatmap(rad2deg(theta));
colormap(jet(512));

% % Threshold on density
thres_rho=double(rho>0.3);

% % PLOT DENSITIES and ORIENTATION as IMAGE
figure(3)
imagesc(thres_rho);
colormap(jet(512));
hold on;
quiver(XC, YC, cos((-theta)) , sin((-theta)),0.6,'y','ShowArrowHead', false);axis equal; axis tight; axis off;
 
 
% POSTPROCESSING 
K = (1/(filt.^2))*ones(filt);
theta = conv2(theta,K,'same');

% % PLOT DENSITIES and ORIENTATION as HEATMAP
figure(4);subplot(211);
heatmap(rho);
subplot(212);
heatmap(rad2deg(theta));
colormap(jet(512));

% % PLOT DENSITIES and ORIENTATION as IMAGE with POSTPROCESSING
thres_rho=double(rho>0.3);

figure(5)
imagesc(thres_rho);
colormap(jet(512));
hold on;
quiver(XC, YC, cos((-theta)) , sin((-theta)),0.6,'y','ShowArrowHead', false);axis equal; axis tight; axis off;



% % dim must be odd
% dim = 1;
% sigma = dim/7;
% filter = gaussian2d(dim,sigma);
% dim = size(filter,1);
% thetaOpt = padarray(theta,[floor(dim/2) floor(dim/2)],'symmetric','both');
% thetaOpt = conv2(thetaOpt,filter,'valid');
% %thetaOpt = reshape(thetaOpt,1,[]);
% 
% figure(2)
%  colormap(gray); imagesc(rho_plot); axis equal; axis tight; axis off;
% hold on;
%  quiver(XC, YC, cos(thetaOpt) , sin(thetaOpt),0.6,'r','ShowArrowHead', false);axis equal; axis tight; axis off;
% %





