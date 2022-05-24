# SOMP
Solid Orthotropic Material with Penalisation

![SOMP](Images/OutPost.jpg)
The code is based on top99: less efficient, restricted to 2D.
But More readable for beginners ;)

## Tutorial 

main.m : main programm setup the optimization problem solved with fmincode

x0 is the initial design vector x0 = [rho0(:);theta0(:)];
global nelx nely vol volfrac ang angle  penal rmin % global variable

function [dcn]=check(nelx,nely,rmin,x,dc) : top99 MESH-INDEPENDENCY FILTER

[c, dt]=top_obj(x) : output compliance c and dc/drho, dc/dtheta

function [cneq, ceq, gradc, gradceq] = myConstrFcn(x) : output nonlinear constraints and derivative

function [KE,dKE]=lkOd(angle); CLT for 1-layer composite membrane fully integrated Ke (8x8 matrix), and derivative with respect to angle, called in FE.m
with fixed material

Ex=1;
Ey=5;
nuxy = 0.3;
nuyx = 0.3;


function [KE,dKE]=lkOd_laminate(angle); CLT for 1-layer composite membrane fully integrated Ke (8x8 matrix), and derivative with respect to angle, called in FE.m
with fixed material

Ex=44.8e+03; % longitudinal Elastic modulus [MPa]
Ey=4.2e+03; % transversal Elastic modulus [MPa]
%Glt=1.9e+03; % Shear Modulus [MPa]
nuxy=0.49; % Poisson ratio
nuyx=nuxy*Ey/Ex;

integK_laminate.m is the symbolic integration of Ke for a fixed material

function [U]=FE(nelx,nely,vol,ang,penal); output displacement as a function of 

myOutputFcn.m needed for output of objective function






## Bibliography
Topology and printing orientation optimization of orthotropic material for additive manufacturing
https://yorkspace.library.yorku.ca/xmlui/handle/10315/38783


An Anisotropic Topology Optimization Method For Carbon Fiber-Reinforced Fused Filament Fabrication
https://baylor-ir.tdl.org/handle/2104/9821
 

Three dimensional topology optimization with orthotropic material orientation design for additive manufacturing structures.
https://baylor-ir.tdl.org/handle/2104/10163 


Jiang's journal paper
https://www.mdpi.com/2079-6439/7/2/14/htm

