%%%%%%%%%% ELEMENT STIFFNESS MATRIX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [KE,dKE]=lkOd(angle);
% BilinearQuadElementStiffness This function returns the element
% stiffness matrix for a bilinear quadrilateral element with modulus
% of elasticity E, Poisson’s ratio nu, thickness h, coordinates of
% node 1 (x1,y1), coordinates node 3 (x3,y3), and coordinates of
% node 4 (x4,y4). Use ps = 1 for cases with isotropic material, and ps = 2 for
% cases with orthotropic material. The size of the element stiffness matrix is 8 x 8. 

%angle is given 0
%Ke is constructed symbolically for 
% Ex=1;
% Ey=5;
% nuxy = 0.3;
% nuyx = 0.3;

%angle belongto [0 1]
%map radians -pi/2 too pi/2
%T = angle*pi-pi/2;
T=angle;



z=[535/273 - (20*cos(T)*sin(T))/91 - (400*cos(T)^2)/273, 5/28 - (15*sin(T)^2)/91 - (50*sin(2*T))/91, (400*cos(T)^2)/273 + (20*cos(T)*sin(T))/91 - 965/546,            (15*cos(T)^2)/91 - (50*sin(2*T))/91 - 5/28, (200*cos(T)^2)/273 + (10*cos(T)*sin(T))/91 - 535/546,            (50*sin(2*T))/91 + (15*sin(T)^2)/91 - 5/28, 215/273 - (10*cos(T)*sin(T))/91 - (200*cos(T)^2)/273,            (50*sin(2*T))/91 - (15*cos(T)^2)/91 + 5/28;   
    (15*cos(T)^2)/91 - (100*cos(T)*sin(T))/91 + 5/364,    (400*cos(T)^2)/273 + (20*cos(T)*sin(T))/91 + 45/91,     (100*cos(T)*sin(T))/91 - (15*cos(T)^2)/91 + 5/28,     (200*cos(T)^2)/273 + (10*cos(T)*sin(T))/91 + 5/91,    (100*cos(T)*sin(T))/91 - (15*cos(T)^2)/91 - 5/364, - (200*cos(T)^2)/273 - (10*cos(T)*sin(T))/91 - 45/182,     (15*cos(T)^2)/91 - (100*cos(T)*sin(T))/91 - 5/28, - (400*cos(T)^2)/273 - (20*cos(T)*sin(T))/91 - 55/182;
    (400*cos(T)^2)/273 + (20*cos(T)*sin(T))/91 - 965/546,            (50*sin(2*T))/91 - (15*cos(T)^2)/91 + 5/28, 535/273 - (20*cos(T)*sin(T))/91 - (400*cos(T)^2)/273,            (50*sin(2*T))/91 + (15*sin(T)^2)/91 - 5/28, 215/273 - (10*cos(T)*sin(T))/91 - (200*cos(T)^2)/273,            (15*cos(T)^2)/91 - (50*sin(2*T))/91 - 5/28, (200*cos(T)^2)/273 + (10*cos(T)*sin(T))/91 - 535/546,            5/28 - (15*sin(T)^2)/91 - (50*sin(2*T))/91;
    (15*cos(T)^2)/91 - (100*cos(T)*sin(T))/91 - 5/28,     (200*cos(T)^2)/273 + (10*cos(T)*sin(T))/91 + 5/91,    (100*cos(T)*sin(T))/91 - (15*cos(T)^2)/91 - 5/364,    (400*cos(T)^2)/273 + (20*cos(T)*sin(T))/91 + 45/91,     (100*cos(T)*sin(T))/91 - (15*cos(T)^2)/91 + 5/28, - (400*cos(T)^2)/273 - (20*cos(T)*sin(T))/91 - 55/182,    (15*cos(T)^2)/91 - (100*cos(T)*sin(T))/91 + 5/364, - (200*cos(T)^2)/273 - (10*cos(T)*sin(T))/91 - 45/182;
    (200*cos(T)^2)/273 + (10*cos(T)*sin(T))/91 - 535/546,            (50*sin(2*T))/91 + (15*sin(T)^2)/91 - 5/28, 215/273 - (10*cos(T)*sin(T))/91 - (200*cos(T)^2)/273,            (50*sin(2*T))/91 - (15*cos(T)^2)/91 + 5/28, 535/273 - (20*cos(T)*sin(T))/91 - (400*cos(T)^2)/273,            5/28 - (15*sin(T)^2)/91 - (50*sin(2*T))/91, (400*cos(T)^2)/273 + (20*cos(T)*sin(T))/91 - 965/546,            (15*cos(T)^2)/91 - (50*sin(2*T))/91 - 5/28;
    (100*cos(T)*sin(T))/91 - (15*cos(T)^2)/91 - 5/364, - (200*cos(T)^2)/273 - (10*cos(T)*sin(T))/91 - 45/182,     (15*cos(T)^2)/91 - (100*cos(T)*sin(T))/91 - 5/28, - (400*cos(T)^2)/273 - (20*cos(T)*sin(T))/91 - 55/182,    (15*cos(T)^2)/91 - (100*cos(T)*sin(T))/91 + 5/364,    (400*cos(T)^2)/273 + (20*cos(T)*sin(T))/91 + 45/91,     (100*cos(T)*sin(T))/91 - (15*cos(T)^2)/91 + 5/28,     (200*cos(T)^2)/273 + (10*cos(T)*sin(T))/91 + 5/91;
    215/273 - (10*cos(T)*sin(T))/91 - (200*cos(T)^2)/273,            (15*cos(T)^2)/91 - (50*sin(2*T))/91 - 5/28, (200*cos(T)^2)/273 + (10*cos(T)*sin(T))/91 - 535/546,            5/28 - (15*sin(T)^2)/91 - (50*sin(2*T))/91, (400*cos(T)^2)/273 + (20*cos(T)*sin(T))/91 - 965/546,            (50*sin(2*T))/91 - (15*cos(T)^2)/91 + 5/28, 535/273 - (20*cos(T)*sin(T))/91 - (400*cos(T)^2)/273,            (50*sin(2*T))/91 + (15*sin(T)^2)/91 - 5/28;
    (100*cos(T)*sin(T))/91 - (15*cos(T)^2)/91 + 5/28, - (400*cos(T)^2)/273 - (20*cos(T)*sin(T))/91 - 55/182,    (15*cos(T)^2)/91 - (100*cos(T)*sin(T))/91 + 5/364, - (200*cos(T)^2)/273 - (10*cos(T)*sin(T))/91 - 45/182,     (15*cos(T)^2)/91 - (100*cos(T)*sin(T))/91 - 5/28,     (200*cos(T)^2)/273 + (10*cos(T)*sin(T))/91 + 5/91,    (100*cos(T)*sin(T))/91 - (15*cos(T)^2)/91 - 5/364,    (400*cos(T)^2)/273 + (20*cos(T)*sin(T))/91 + 45/91];

dw =[(20*sin(T)^2)/91 - (20*cos(T)^2)/91 + (800*cos(T)*sin(T))/273,                   - (100*cos(2*T))/91 - (30*cos(T)*sin(T))/91, (20*cos(T)^2)/91 - (20*sin(T)^2)/91 - (800*cos(T)*sin(T))/273,                   - (100*cos(2*T))/91 - (30*cos(T)*sin(T))/91, (10*cos(T)^2)/91 - (10*sin(T)^2)/91 - (400*cos(T)*sin(T))/273,                     (100*cos(2*T))/91 + (30*cos(T)*sin(T))/91, (10*sin(T)^2)/91 - (10*cos(T)^2)/91 + (400*cos(T)*sin(T))/273,                     (100*cos(2*T))/91 + (30*cos(T)*sin(T))/91;(100*sin(T)^2)/91 - (100*cos(T)^2)/91 - (30*cos(T)*sin(T))/91, (20*cos(T)^2)/91 - (20*sin(T)^2)/91 - (800*cos(T)*sin(T))/273, (100*cos(T)^2)/91 - (100*sin(T)^2)/91 + (30*cos(T)*sin(T))/91, (10*cos(T)^2)/91 - (10*sin(T)^2)/91 - (400*cos(T)*sin(T))/273, (100*cos(T)^2)/91 - (100*sin(T)^2)/91 + (30*cos(T)*sin(T))/91, (10*sin(T)^2)/91 - (10*cos(T)^2)/91 + (400*cos(T)*sin(T))/273, (100*sin(T)^2)/91 - (100*cos(T)^2)/91 - (30*cos(T)*sin(T))/91, (20*sin(T)^2)/91 - (20*cos(T)^2)/91 + (800*cos(T)*sin(T))/273;
    (20*cos(T)^2)/91 - (20*sin(T)^2)/91 - (800*cos(T)*sin(T))/273,                     (100*cos(2*T))/91 + (30*cos(T)*sin(T))/91, (20*sin(T)^2)/91 - (20*cos(T)^2)/91 + (800*cos(T)*sin(T))/273,                     (100*cos(2*T))/91 + (30*cos(T)*sin(T))/91, (10*sin(T)^2)/91 - (10*cos(T)^2)/91 + (400*cos(T)*sin(T))/273,                   - (100*cos(2*T))/91 - (30*cos(T)*sin(T))/91, (10*cos(T)^2)/91 - (10*sin(T)^2)/91 - (400*cos(T)*sin(T))/273,                   - (100*cos(2*T))/91 - (30*cos(T)*sin(T))/91;
    (100*sin(T)^2)/91 - (100*cos(T)^2)/91 - (30*cos(T)*sin(T))/91, (10*cos(T)^2)/91 - (10*sin(T)^2)/91 - (400*cos(T)*sin(T))/273, (100*cos(T)^2)/91 - (100*sin(T)^2)/91 + (30*cos(T)*sin(T))/91, (20*cos(T)^2)/91 - (20*sin(T)^2)/91 - (800*cos(T)*sin(T))/273, (100*cos(T)^2)/91 - (100*sin(T)^2)/91 + (30*cos(T)*sin(T))/91, (20*sin(T)^2)/91 - (20*cos(T)^2)/91 + (800*cos(T)*sin(T))/273, (100*sin(T)^2)/91 - (100*cos(T)^2)/91 - (30*cos(T)*sin(T))/91, (10*sin(T)^2)/91 - (10*cos(T)^2)/91 + (400*cos(T)*sin(T))/273;
    (10*cos(T)^2)/91 - (10*sin(T)^2)/91 - (400*cos(T)*sin(T))/273,                     (100*cos(2*T))/91 + (30*cos(T)*sin(T))/91, (10*sin(T)^2)/91 - (10*cos(T)^2)/91 + (400*cos(T)*sin(T))/273,                     (100*cos(2*T))/91 + (30*cos(T)*sin(T))/91, (20*sin(T)^2)/91 - (20*cos(T)^2)/91 + (800*cos(T)*sin(T))/273,                   - (100*cos(2*T))/91 - (30*cos(T)*sin(T))/91, (20*cos(T)^2)/91 - (20*sin(T)^2)/91 - (800*cos(T)*sin(T))/273,                   - (100*cos(2*T))/91 - (30*cos(T)*sin(T))/91;
    (100*cos(T)^2)/91 - (100*sin(T)^2)/91 + (30*cos(T)*sin(T))/91, (10*sin(T)^2)/91 - (10*cos(T)^2)/91 + (400*cos(T)*sin(T))/273, (100*sin(T)^2)/91 - (100*cos(T)^2)/91 - (30*cos(T)*sin(T))/91, (20*sin(T)^2)/91 - (20*cos(T)^2)/91 + (800*cos(T)*sin(T))/273, (100*sin(T)^2)/91 - (100*cos(T)^2)/91 - (30*cos(T)*sin(T))/91, (20*cos(T)^2)/91 - (20*sin(T)^2)/91 - (800*cos(T)*sin(T))/273, (100*cos(T)^2)/91 - (100*sin(T)^2)/91 + (30*cos(T)*sin(T))/91, (10*cos(T)^2)/91 - (10*sin(T)^2)/91 - (400*cos(T)*sin(T))/273;
    (10*sin(T)^2)/91 - (10*cos(T)^2)/91 + (400*cos(T)*sin(T))/273,                   - (100*cos(2*T))/91 - (30*cos(T)*sin(T))/91, (10*cos(T)^2)/91 - (10*sin(T)^2)/91 - (400*cos(T)*sin(T))/273,                   - (100*cos(2*T))/91 - (30*cos(T)*sin(T))/91, (20*cos(T)^2)/91 - (20*sin(T)^2)/91 - (800*cos(T)*sin(T))/273,                     (100*cos(2*T))/91 + (30*cos(T)*sin(T))/91, (20*sin(T)^2)/91 - (20*cos(T)^2)/91 + (800*cos(T)*sin(T))/273,                     (100*cos(2*T))/91 + (30*cos(T)*sin(T))/91;
    (100*cos(T)^2)/91 - (100*sin(T)^2)/91 + (30*cos(T)*sin(T))/91, (20*sin(T)^2)/91 - (20*cos(T)^2)/91 + (800*cos(T)*sin(T))/273, (100*sin(T)^2)/91 - (100*cos(T)^2)/91 - (30*cos(T)*sin(T))/91, (10*sin(T)^2)/91 - (10*cos(T)^2)/91 + (400*cos(T)*sin(T))/273, (100*sin(T)^2)/91 - (100*cos(T)^2)/91 - (30*cos(T)*sin(T))/91, (10*cos(T)^2)/91 - (10*sin(T)^2)/91 - (400*cos(T)*sin(T))/273, (100*cos(T)^2)/91 - (100*sin(T)^2)/91 + (30*cos(T)*sin(T))/91, (20*cos(T)^2)/91 - (20*sin(T)^2)/91 - (800*cos(T)*sin(T))/273];
w = double(z);
KE=w;
dKE=double(dw);


end