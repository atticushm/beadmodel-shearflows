% this code plots figures outputted from shearTest_jeffery.m simulations runs.

clc; close all;
load('workspace.mat')

%% Initial config.
h = subplot(1,2,1);

intConfig        = zeros(3,Nb);                                                          
intConfig(2,:,1) = linspace(-.5,.5,Nb);            

box on
hold on
intFig = draw_bm(intConfig);
axis square
xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
title('Initial config','interpreter','latex')
str = sprintf('N=%g',Nb);
annotation('textbox','string',str,'FitBoxToText','on')
hold off

%% Final config.
h = subplot(1,2,2);

endConfig = xb;

box on
hold on
endFig = draw_bm(endConfig);
axis square
xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
title('Final config.','interpreter','latex')
hold off


%% Function definitions.

function fig = draw_bm(x)
% Draws a snapshot of the bead model given position data xb.
% This version does not draw the centre of mass (not interesting in jeffery
% orbit experiment).

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
% scatter3(xc(1),xc(2),xc(3),50, 'rx','Linewidth',2)

end