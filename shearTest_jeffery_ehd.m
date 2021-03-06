% This script uses the analytical results from Kim and Karilla pg125 to model the
% Jeffery orbit of a torque-free axisymmetric slender body in linear shear flow.

% Rotating body is defined by unit vector d lying along the axis of symmetry of the body.
% r is the aspect ratio of the object.
% C is the orbital constant.
% gamma (representing dot(gamma) ) is the coefficient on the strength of the shear flow.

%% Setup

clear all; close all;
viewInt = 1000;

tMin = 0;
tMax = 1e-2;
dt   = 1e-9;
t    = [tMin:dt:tMax];
Nt   = length(t);

r     = 1000;
gamma = -1e4;
d     = zeros(3,Nt);
C     = 100000;

% d(:,1) = [0;0.5;0];
% d(:,1) = d(:,1)/norm(d(:,1));   % initial config.

%% Simulate over time
rodFig = figure;
for n = 1:Nt

    ph   = solvePhi(r,gamma,t(n));
    thet = solveTheta(C,r,ph);

    d1 = sin(thet)*cos(ph);
    d2 = sin(thet)*sin(ph);
    d3 = cos(thet);

    d(:,n)  = [d1;d2;d3];

    if mod(n,viewInt)== 0
        clf
        rodFig = drawRod(d(:,n));
        pause(0.01)
    end

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

phi = -r*tan( (gamma*t)/(r+(1/r)) );
phi = atan(phi);

end  % function

function [ rodFig,vecEnd ] = drawRod( d)
% Draws rod given axis of symmetry vector d at time t(n)

vecx   = linspace(-d(1),d(1),2);
vecy   = linspace(-d(2),d(2),2);
vecz   = linspace(-d(3),d(3),2);

vecEnd = [vecx(end);vecy(end);vecz(end)];

subplot(1,2,1)
box on
hold on
rodFig = plot3(vecx,vecy,vecz,'k','Linewidth',4);
axis([-1.1 1.1 -1.1 1.1 -1 1])
axis square
xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
zlabel('$z$','interpreter','latex')
hold off
view(2) % view along z axis.

subplot(1,2,2)
box on
hold on
rodFig = plot3(vecx,vecy,vecz,'k','Linewidth',4);
axis([-1.1 1.1 -1.1 1.1 -1 1])
axis square
xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
zlabel('$z$','interpreter','latex')
hold off
view(0,0) % view along y axis.

end  % function
