% This script uses the analytical results from Kim and Karilla pg125 to model the
% Jeffery orbit of a torque-free axisymmetric slender body in linear shear flow.

% Rotating body is defined by unit vector d lying along the axis of symmetry of the body.
% r is the aspect ratio of the object.
% C is the orbital constant.
% gamma (representing dot(gamma) ) is the coefficient on the strength of the shear flow.

%% Setup

clear all; close all;

tMin = 0;
tMax = 1e-3;
dt   = 1e-5;
t    = [tMin:dt:tMax];
Nt   = length(t);

r     = 1000;
gamma = -1e4;
d     = zeros(3,Nt);
C     = 1;

d(:,1) = [0;0.5;0];
d(:,1) = d(:,1)/norm(d(:,1));   % initial config.

%% Simulate over time
figure
for n = 1:Nt

    ph   = solvePhi(r,gamma,t(n));
    thet = solveTheta(C,r,ph);

    d1 = sin(thet)*cos(ph);
    d2 = sin(thet)*sin(ph);
    d3 = cos(thet);

    d(:,n)  = [d1;d2;d3];

    rodFig = drawRod(d(:,n));
    pause(0.01)

end

%% Completion

fprintf('Script complete. \n')

%% Functions

function [ theta ] = solveTheta( C,r,phi )
% Analytical solution for the theta component of the slender body orientation vector.

theta = C*r*(r^2*(cos(phi)^2) + sin(phi)^2 )^(-1/2);
theta = atan(theta);

end  % function

function [ phi ] = solvePhi( r,gamma,t )
% Analytical solution for the phi component of the slender body orientation vector.

phi = r*tan( gamma*t/(r+(1/r)) );
phi = atan(phi);

end  % function

function [ rodFig ] = drawRod( d)
% Draws rod given axis of symmetry vector d at time t(n)

clf
box on
hold on

vecx   = linspace(0,d(1),2);
vecy   = linspace(0,d(2),2);
vecz   = linspace(0,d(3),2);

vecx = [-vecx,vecx];
vecy = [-vecy,vecy];
vecz = [-vecz,vecz];

rodFig = plot3(vecx,vecy,vecz,'k','Linewidth',4);

axis([-1.1 1.1 -1.1 1.1 -1 1])
xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
zlabel('$z$','interpreter','latex')
hold off

view(2)

end  % function
