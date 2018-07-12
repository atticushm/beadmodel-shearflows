% In this scipt, a filament positioned symmetrically about the origin, is placed
% in a shear flow. The filament is allowed to rotate a quarter turn, at which point
% the sim is ended.

% this script compares the previously created shearTest_jeffery.m to
% results from elastohydrodynamics (analytical results from Kim and
% Karilla)

clear all; %close all
viewInt = 1000;
figure

%% Setup - bead model.
a    = 0.001;                                   % radius of sphere approximating each bead; non-dimensional.
mu   = 1e-3;                                    % dynamic viscosity of surrounding fluid.
L    = 1;                                       % length scale; dimensionless.
calS = 5e4;                                     % non-dim constant S.

Nb = 11;                                        % number of beads.
epsilon = a;                                    % regularisation parameter in reg. stokeslet method.

tMin = 0;
tMax = 1e-2;
dt   = 1e-9;
t    = [tMin:dt:tMax];
Nt   = length(t);                            

b0 = 1/(Nb-1);                     

%% Set initial position.
xb = zeros(3,Nb);                               % stores the xyz positions of nP particles at nt time-steps.
xc = zeros(3);                                  % stores the xyz positions of the centre of mass at nt time-steps.

xb(2,:,1) = linspace(-L/2,L/2,Nb)/L;            % position filament along y axis.
xc(:,1)   = mean(xb(:,:),2);                    % initial centre of mass.

%% Main
count = 10;
tic;
for n = 1:Nt

    F  = zeros(3,Nb);
    U  = zeros(3,Nb);

    %% bending forces.
    Fb = get_bending_forces(xb);
    F = F + (Nb-1).*Fb;

    %% spring forces.
    Fs = get_spring_forces(xb, b0);
    F  = F + calS.*Fs;

    %% hydrodynamic interactions.
    for p = 1:Nb
        clear stokeslets
        stokeslets = get_stokeslets(xb, xb(:,p),epsilon);
        G = reshape(F',[1,3*Nb])';
        U(:,p) = stokeslets*G ;
    end
    Us = [-1e4.*xb(2,:); zeros(1,Nb); zeros(1,Nb)];    % xy planar shear flow
    U  = U + Us;

    %% check arclength.
    s = 0;
    for p = 1:Nb-1
        s = s + norm(xb(:,p+1) - xb(:,p));
        if s > 1.2
            fprintf('Filament inextensibility broken... \n')
            breakPlot = draw_bm_buckling(x,xc);
        end
    end

    %% animated plot.
    if mod(n,viewInt)== 0
        clf
        hold on
        figBM = draw_bm_buckling(xb,xc);
        xlabel('x')
        ylabel('y')
        pause(0.01)
        hold off
    end

    %% update positions and centre of mass.
    xb = xb + U*dt;
    xc = mean(xb,2);

    %% check filament rotation and end code if quarter-turn completed.
    if abs(xb(2,Nb)) < 0.01 && abs(xb(2,1)) < 0.01      % ie first and last beads close to x axis.
        tFin = t(n);
        fprintf('Filament aligned along x axis; script stopping at t=%g...',tFin)
        return
    end

    %% script progress counter.
    if mod(n,(Nt-1)/10)==0
        fprintf('%g perc. of time steps completed... \n',count)
        save('workspace_multiNb_part.mat')
        count = count+10;
    end

end
runtime = toc;
save(sprintf('workspace_Nb=%g.mat',Nb))
fprintf('Workspace saved. \n')
fprintf('Final time: t=%g. \n',t(n))
fprintf('Script completed in %g minutes. \n', runtime/60)

%% Function definitions.

function fig = draw_bm_buckling(x,xc)

Nb = size(x,2);

fig1 = scatter3(x(1,:),x(2,:),x(3,:),30,'ro','filled');
axis([-1.1 1.1 -1.1 1.1 0 1])
view(2)
% axis equal
xlabel('x')
ylabel('y')
zlabel('z')
grid on
hold on
lines = zeros(3,100,Nb-1);
for kk = 1:Nb-1
    lines(1,:,kk) = linspace(x(1,kk),x(1,kk+1));
    lines(2,:,kk) = linspace(x(2,kk),x(2,kk+1));
    lines(3,:,kk) = linspace(x(3,kk),x(3,kk+1));
    fig = plot3(squeeze(lines(1,:,kk)),squeeze(lines(2,:,kk)),squeeze(lines(3,:,kk)),'r');
end
scatter3(xc(1),xc(2),xc(3),50, 'rx','Linewidth',2)

end

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
