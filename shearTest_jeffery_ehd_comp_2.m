% In this scipt, a filament positioned symmetrically about the origin, is placed
% in a shear flow. The filament is allowed to rotate a quarter turn, at which point
% the sim is ended.

% this script compares the previously created shearTest_jeffery.m to
% results from elastohydrodynamics (analytical results from Kim and
% Karilla)

% this script does not draw figures. Use xx_animated.m version for that.

clear all; %close all
figs.viewInt = 100;

%% Setup: bead model.
beads.a    = 0.001;
beads.L    = 1;
beads.calS = 5e4;

beads.Nb      = 11;
beads.epsilon = beads.a;
beads.b0      = 1/(beads.Nb-1);

time.tMin = 0;
time.tMax = 1e-2;
time.dt   = 1e-9;
time.t    = [time.tMin:time.dt:time.tMax];
time.Nt   = length(time.t);

%% Setup: elastohydrodynamics.
ehd.r     = 1000;
ehd.gamma = -1e4;
ehd.d     = zeros(3,time.Nt);
ehd.C     = 1e10;

%% Set initial position.
beads.pos        = zeros(3,beads.Nb);
beads.pos(1,:)   = linspace(-beads.L/2,beads.L/2,beads.Nb)/beads.L;

%% Main
count1 = 10;
count2 = 1;
tic;
for n = 1:time.Nt

    beads.F  = zeros(3,beads.Nb);
    beads.U  = zeros(3,beads.Nb);

    %% bending forces.
    beads.Fb = get_bending_forces(beads.pos);
    beads.F = beads.F + (beads.Nb-1).*beads.Fb;

    %% spring forces.
    beads.Fs = get_spring_forces(beads.pos, beads.b0);
    beads.F  = beads.F + beads.calS.*beads.Fs;

    %% hydrodynamic interactions.
    for p = 1:beads.Nb
        clear stokeslets
        stokeslets = get_stokeslets(beads.pos, beads.pos(:,p),beads.epsilon);
        G = reshape(beads.F',[1,3*beads.Nb])';
        beads.U(:,p) = stokeslets*G ;
    end
    beads.Us = [zeros(1,beads.Nb); -ehd.gamma.*beads.pos(1,:); zeros(1,beads.Nb)];    % planar shear flow
    beads.U  = beads.U + beads.Us;

    %% EHD analytical solution.
    ehd.ph   = solvePhi(ehd.r,ehd.gamma,time.t(n));
    ehd.thet = solveTheta(ehd.C,ehd.r,ehd.ph);

    ehd.d1 = sin(ehd.thet)*cos(ehd.ph);
    ehd.d2 = sin(ehd.thet)*sin(ehd.ph);
    ehd.d3 = cos(ehd.thet);

    ehd.d(:,n)  = [ehd.d1;ehd.d2;ehd.d3];

    %% check arclength.
    beads.s = 0;
    for p = 1:beads.Nb-1
        beads.s = beads.s + norm(beads.pos(:,p+1) - beads.pos(:,p));
        if beads.s > 1.2
            fprintf('Filament inextensibility broken... \n')
            figs.breakPlot = drawBeads(beads.pos(:,n));
        end
    end

    %% animated plot.
    if mod(n,figs.viewInt)== 0
        clf
        hold on
        figs.f1 = figure('visible','off');
        
        [figs.f1,ehd.vecEnd] = drawRod(ehd.d(:,n));
        
        hold off      
        error.theta(count2) = delTheta(ehd.vecEnd,beads.pos(:,end));
        error.t(count2)  = time.t(n);
        count2 = count2+1;
   
    end

    %% update positions.
    beads.pos = beads.pos + beads.U*time.dt;
    
    %% check filament rotation and end code if quarter-turn completed.
%     if abs(xb(1,Nb)) < 0.01 && abs(xb(1,1)) < 0.01      % ie first and last beads close to y axis.
%         time.tFin = time.t(n);
%         fprintf('Filament aligned along x axis; script stopping at t=%g...', time.tFin)
%         return
%     end

    %% script progress counter.
    if mod(n,(time.Nt-1)/10)==0
        fprintf('%g perc. of time steps completed... \n',count1)
        save('workspace_jeffery_ehd_comp_part.mat')
        count1 = count1+10;
    end

end

time.runtime = toc;
save(sprintf('workspace_jeffery_ehd_comp.mat'))
fprintf('Workspace saved. \n')
fprintf('Final time: t=%g. \n',time.t(n))
fprintf('Script completed in %g minutes. \n', time.runtime/60)

%% Function definitions.

function thet = delTheta(x1,x2)
% Calculates angular difference between two lines given by position vector
% from the origin.

nx1   = norm(x1);
nx2   = norm(x2);
dx1x2 = dot(x1,x2);

thet = acos(dx1x2/nx1/nx2);

end % function

function fig = drawBeads(x)
% Modified version of bead figure drawing code excluding CoM plot.

Nb   = size(x,2);
fig1 = scatter3(x(1,:),x(2,:),x(3,:),30,'ro','filled');
axis([-1.1 1.1 -1.1 1.1 0 1])
axis square
view(2)
% axis equal
xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
zlabel('$z$','interpreter','latex')
grid on
hold on
lines = zeros(3,100,Nb-1);
for kk = 1:Nb-1
    lines(1,:,kk) = linspace(x(1,kk),x(1,kk+1));
    lines(2,:,kk) = linspace(x(2,kk),x(2,kk+1));
    lines(3,:,kk) = linspace(x(3,kk),x(3,kk+1));
    fig = plot3(squeeze(lines(1,:,kk)),squeeze(lines(2,:,kk)),squeeze(lines(3,:,kk)),'r','Linewidth',2);
end
end

function [ rodFig,vecEnd ] = drawRod(d)
% Draws rod given axis of symmetry vector d at time t(n).
% Reduced from shearTest_jeffery_ehd.m version of code - only draws xy
% planar view (ie no verification extra view drawn).

vecx   = linspace(-d(1)/2,d(1)/2,2);
vecy   = linspace(-d(2)/2,d(2)/2,2);
vecz   = linspace(-d(3)/2,d(3)/2,2);

vecEnd = [vecx(2);vecy(2);vecz(2)];
box on
hold on
rodFig = plot3(vecx,vecy,vecz,'b','Linewidth',6);
axis([-1.1 1.1 -1.1 1.1 -1 1])
axis square
xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
zlabel('$z$','interpreter','latex')
% hold off
view(2) % view along z axis.

end  % function

function F = get_bending_forces(x)
%% GET_BENDING_FORCES
%   Calculates restorative bending forces at Np 3-dimensional beads to keep
%   a chain aligned.
%   INPUTS: x (3xNp matrix of bead positions)
%   OUTPUTS: F (bending forces)

%% Main.

Np = size(x,2);
F = zeros(3,Np);

% bead 1.
bj   = x(:,1) - x(:,2);
bjp1 = x(:,2) - x(:,3);
F(:,1) = bjp1/norm(bj)/norm(bjp1) + dot(bj,bjp1)*(-bj)/norm(bj)^3/norm(bjp1);

% bead 2.
bjn1 = x(:,1) - x(:,2);
bj   = x(:,2) - x(:,3);
bjp1 = x(:,3) - x(:,4);
Fb = (bjn1 - bj)/norm(bjn1)/norm(bj) + dot(bjn1,bj)*(bjn1/norm(bjn1)^3/norm(bj) - bj/norm(bjn1)/norm(bj)^3);
Fc = bjp1/norm(bj)/norm(bjp1) + dot(bj,bjp1)*(-bj)/norm(bj)^3/norm(bjp1);
F(:,2) = Fb+Fc;

% beads 3 to Np-2.
for kk = 3:(Np-2)
    bjn2 = x(:,kk-2) - x(:,kk-1);
    bjn1 = x(:,kk-1) - x(:,kk);
    bj   = x(:,kk)   - x(:,kk+1);
    bjp1 = x(:,kk+1) - x(:,kk+2);
    Fa = -bjn2/norm(bjn2)/norm(bjn1)   + dot(bjn2,bjn1)*bjn1/norm(bjn2)/norm(bjn1)^3;
    Fb = (bjn1-bj)/norm(bjn1)/norm(bj) + dot(bjn1,bj)*(bjn1/norm(bjn1)^3/norm(bj) - bj/norm(bjn1)/norm(bj)^3);
    Fc = bjp1/norm(bj)/norm(bjp1)      + dot(bj,bjp1)*(-bj)/norm(bj)^3/norm(bjp1);
    F(:,kk) = Fa+Fb+Fc;
end

% bead Np-1.
bjn2 = x(:,Np-3) - x(:,Np-2);
bjn1 = x(:,Np-2) - x(:,Np-1);
bj   = x(:,Np-1) - x(:,Np);
Fa = -bjn2/norm(bjn2)/norm(bjn1)   + dot(bjn2,bjn1)*bjn1/norm(bjn2)/norm(bjn1)^3;
Fb = (bjn1-bj)/norm(bjn1)/norm(bj) + dot(bjn1,bj)*(bjn1/norm(bjn1)^3/norm(bj) - bj/norm(bjn1)/norm(bj)^3);
F(:,Np-1) = Fa+Fb;

% bead Np.
bjn2 = x(:,Np-2) - x(:,Np-1);
bjn1 = x(:,Np-1) - x(:,Np);
F(:,Np) = -bjn2/norm(bjn2)/norm(bjn1) + dot(bjn2,bjn1)*bjn1/norm(bjn2)/norm(bjn1)^3;

end

function F = get_spring_forces(x, b0)
%% GET_SPRING_FORCES
%   Calculates the spring forces on a connected chain of beads in 3-dimensional space.
%   INPUTS: X (3xNp matrix of bead locations)
%           b0 (equilibrium distance between beads)
%   OUTPUTS: F (3xNp matrix of forces)

%% Main.

Np = size(x,2);
F  = zeros(3,Np);

for p = 1:Np

    if p > 1
        bpn1  = x(:,p-1) - x(:,p);
        Nbpn1 = norm(bpn1);
        F1    = -2*(Nbpn1 - b0)*(-bpn1/Nbpn1);
    else
        F1 = 0;
    end
    if p < Np
        bp   = x(:,p) - x(:,p+1);
        Nbp  = norm(bp);
        F2   = -2*(Nbp - b0)*(bp/Nbp);
    else
        F2 = 0;
    end

    F(:,p) = F1 + F2;
end

end

function S = get_stokeslets(x, x0, eps)
%GET_REG_STOKESLETS Calculates the regularised stokeslets for a given list of nodes and regularisation parameter.
%   ______________________________________________________________________
%   INPUTS
%   q_nodes: positions of quadrature nodes along the boundary of object.
%   f_nodes: positions of force nodes on the boundary.
%   eps: the regularisation parameter epsilon.
%   mu: the dynamic viscosity of the fluid in question.
%
%   q_nodes and f_nodes must be inputted as a 3xN/3xQ
%   ______________________________________________________________________
%   OUTPUTS
%   S: the block matrix of stokeslets
%   ______________________________________________________________________

eps2 = eps^2;
Q = size(x0,2);

for p = 1:Q

    rx = x(1,:) - x0(1,p);
    ry = x(2,:) - x0(2,p);
    rz = x(3,:) - x0(3,p);

    r2 = rx.^2 + ry.^2 + rz.^2;
    r_eps = (r2 + eps2).^1.5;
    r_num = (r2 + 2*eps2);

    Sxx(p,:) = (r_num + rx.*rx)./r_eps;
    Sxy(p,:) = (rx.*ry)./r_eps;
    Sxz(p,:) = (rx.*rz)./r_eps;

    Syx(p,:) = (ry.*rx)./r_eps;
    Syy(p,:) = (r_num + ry.*ry)./r_eps;
    Syz(p,:) = (ry.*rz)./r_eps;

    Szx(p,:) = (rz.*rx)./r_eps;
    Szy(p,:) = (rz.*ry)./r_eps;
    Szz(p,:) = (r_num + rz.*rz)./r_eps;

end

A1 = horzcat(Sxx, Sxy, Sxz);
A2 = horzcat(Syx, Syy, Syz);
A3 = horzcat(Szx, Szy, Szz);

S = vertcat(A1, A2, A3);

end

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

