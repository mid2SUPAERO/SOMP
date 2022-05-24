clear all; close all;
Ex=44.8e+03; % longitudinal Elastic modulus [MPa]
Ey=4.2e+03; % transversal Elastic modulus [MPa]
%Glt=1.9e+03; % Shear Modulus [MPa]
nuxy=0.49; % Poisson ratio
nuyx=nuxy*Ey/Ex;

% Ex=1;
% Ey=5;
% nuxy = 0.3;
% nuyx = 0.3;
h=1;
x1=-1;
y1=-1;
x2=1;
y2=-1;
x3=1;
y3=1;
x4=-1;
y4=1;
ps=2;
%T = angle*3.141592654;
syms s t T;
a = (y1*(s-1)+y2*(-1-s)+y3*(1+s)+y4*(1-s))/4;
b = (y1*(t-1)+y2*(1-t)+y3*(1+t)+y4*(-1-t))/4;
c = (x1*(t-1)+x2*(1-t)+x3*(1+t)+x4*(-1-t))/4;
d = (x1*(s-1)+x2*(-1-s)+x3*(1+s)+x4*(1-s))/4;
B1 = [a*(t-1)/4-b*(s-1)/4 0 ; 0 c*(s-1)/4-d*(t-1)/4 ;
c*(s-1)/4-d*(t-1)/4 a*(t-1)/4-b*(s-1)/4];
B2 = [a*(1-t)/4-b*(-1-s)/4 0 ; 0 c*(-1-s)/4-d*(1-t)/4;
c*(-1-s)/4-d*(1-t)/4 a*(1-t)/4-b*(-1-s)/4];
B3 = [a*(t+1)/4-b*(s+1)/4 0 ; 0 c*(s+1)/4-d*(t+1)/4 ;
c*(s+1)/4-d*(t+1)/4 a*(t+1)/4-b*(s+1)/4];
B4 = [a*(-1-t)/4-b*(1-s)/4 0 ; 0 c*(1-s)/4-d*(-1-t)/4 ;
c*(1-s)/4-d*(-1-t)/4 a*(-1-t)/4-b*(1-s)/4]; 
Bfirst = [B1 B2 B3 B4];
Jfirst = [0 1-t t-s s-1 ; t-1 0 s+1 -s-t ;
s-t -s-1 0 t+1 ; 1-s s+t -t-1 0];
J = [x1 x2 x3 x4]*Jfirst*[y1 ; y2 ; y3 ; y4]/8; 
B = Bfirst/J;
R = [ cos(T) -sin(T) 0; sin(T) cos(T) 0; 0 0 1]; 
if ps==1
D = (Ex/(1-nuxy*nuxy))*[1, nuxy, 0 ; nuxy, 1, 0 ; 0, 0, (1-nuxy)/2]; 
elseif ps == 2
D = [ Ex/(1-(nuxy*nuyx)) (nuyx*Ex)/(1-(nuxy*nuyx)) 0; (nuyx*Ex)/(1-(nuxy*nuyx)) Ey/(1-(nuxy*nuyx)) 0;0 0 Ex/(2*(1+nuxy))];
end
BD = J*transpose(B)*(R)*D*transpose(R)*B;
r = int(int(BD, t, -1, 1), s, -1, 1);
z=h*r;

%T is syms
dw=diff(z,T);

% w = double(z);
% KE=w;
% dKE=double(dw);