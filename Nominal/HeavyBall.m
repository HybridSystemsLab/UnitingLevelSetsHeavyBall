%--------------------------------------------------------------------------
% Matlab M-file Project: HyEQ Toolbox @  Hybrid Systems Laboratory (HSL), 
% https://hybrid.soe.ucsc.edu/software
% http://hybridsimulator.wordpress.com/
% Filename: HeavyBall.m
%--------------------------------------------------------------------------
% Project: Testing out parameters lambda and gamma for fast, oscillatory
% convergence globally and slow, smooth convergence locally.
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%   See also HYEQSOLVER, PLOTARC, PLOTARC3, PLOTFLOWS, PLOTHARC,
%   PLOTHARCCOLOR, PLOTHARCCOLOR3D, PLOTHYBRIDARC, PLOTJUMPS.
%   Copyright @ Hybrid Systems Laboratory (HSL),
%   Revision: 0.0.0.3 Date: 05/20/2015 3:42:00

clear all

set(0,'defaultTextInterpreter','latex');

% global variables
global gamma lambda gamma_0 gamma_1 lambda_0 lambda_1 c_0 c_10 c_1 V0 V1 delta z1Star 
%%%%%%%%%%%%%%%%%%%%%
% setting the globals
%%%%%%%%%%%%%%%%%%%%%

% Heavy-ball constants: 
lambda = 10.5; % Gravity. 
            % For gamma fixed, "large values of  lambda are seen to give rise to slowly converging 
            % solutions resembling the steepest descent’s while smaller values give 
            % rise to fast solutions with oscillations getting wilder as lambda decreases."
gamma = 1/2; % Viscous friction to mass ratio.

lambda_0 = 10.5;
lambda_1 = 1/5;

gamma_0 = 1/2;
gamma_1 = 1/2;

z1Star = 0;

timeToDeltaSlow = 0;
timeToDeltaOscillate = 0;
timeToDeltaUniting = 0;
timeToDeltaIdxSlow = 1;
timeToDeltaIdxOscillate = 1;
timeToDeltaIdxUniting = 1;

z1deltaSlow = 0;
z1deltaOscillate = 0;
z1deltaUniting = 0;

z2deltaSlow = 0;
z2deltaOscillate = 0;
z2deltaUniting = 0;

delta = 0.01;

c_0 = 12.5; % \mathcal{U}_0
c_10 = 6.3; % \mathcal{T}_{1,0}

% initial conditions
z1_0 = -10;
z2_0 = 0;
q_0 = 0;

% Assign initial conditions to vector
x0 = [z1_0;z2_0;q_0];

% simulation horizon
TSPAN=[0 340];
JSPAN = [0 20];

% rule for jumps
% rule = 1 -> priority for jumps
% rule = 2 -> priority for flows
rule = 1;

options = odeset('RelTol',1e-6,'MaxStep',.1);

% Simulate the slow individual heavy ball controller
[tSlow,jSlow,xSlow] = HyEQsolver(@f,@g,@C,@D,...
    x0,TSPAN,JSPAN,rule,options);

% Update lambda, so we can run the fast, oscillating individual heavy ball controller 
lambda = 1/5;

% Then simulate the fast, oscillating individual heavy ball controller
[tOscillate,jOscillate,xOscillate] = HyEQsolver(@f,@g,@C,@D,...
    x0,TSPAN,JSPAN,rule,options);

% Finally, simulate the hybrid closed-loop heavy ball system
[tUniting,jUniting,xUniting] = HyEQsolver(@fU,@gU,@CU,@DU,...
    x0,TSPAN,JSPAN,rule,options);


% Finding time of convergence for the slow controller:
    for i=2:length(xSlow(:,1))
        if (distance(xSlow(i,1)) <= delta) && (distance(xSlow(i-1,1)) > delta)
            timeToDeltaIdxSlow = i;
            z1deltaSlow = xSlow(i,1);
        end
    end
    z2deltaSlow = xSlow(timeToDeltaIdxSlow,2);
    timeToDeltaSlow = tSlow(timeToDeltaIdxSlow,1);
    
% Finding time of convergence for the fast, oscillatory controller
    for i=2:length(xOscillate(:,1))
        if (distance(xOscillate(i,1)) <= delta) && (distance(xOscillate(i-1,1)) > delta)
            timeToDeltaIdxOscillate = i;
            z1deltaOscillate = xOscillate(i,1);
        end
    end
    z2deltaOscillate = xOscillate(timeToDeltaIdxOscillate,2);
    timeToDeltaOscillate = tOscillate(timeToDeltaIdxOscillate,1);
    
% Finding time of convergence for the hybrid closed-loop system
    for i=2:length(xUniting(:,1))
        if (distance(xUniting(i,1)) <= delta) && (distance(xUniting(i-1,1)) > delta)
            timeToDeltaIdxUniting = i;
            z1deltaUniting = xUniting(i,1);
        end
    end
    z2deltaUniting = xUniting(timeToDeltaIdxUniting,2);
    timeToDeltaUniting = tUniting(timeToDeltaIdxUniting,1);
    
    
minarc = min([length(xUniting),length(xSlow),length(xOscillate)]);
ta = [tUniting(1:minarc),tSlow(1:minarc),tOscillate(1:minarc)];
ja = [jUniting(1:minarc),jSlow(1:minarc),jOscillate(1:minarc)];
xa = [xUniting(1:minarc,1),xSlow(1:minarc,1),xOscillate(1:minarc,1)];
xb = [xUniting(1:minarc,2),xSlow(1:minarc,2),xOscillate(1:minarc,2)];

% Create the level curves for plotting.
x1 = -11:0.005:6;
x2 = -4:0.005:6;

[Z1,Z2] = meshgrid(x1,x2);
V0_0 = gamma_0*(0.25*(Z1-z1Star).^2 - CalculateLStar()) + (1/2)*Z2.^2;
V1_1 = gamma_1*(0.25*(Z1-z1Star).^2 - CalculateLStar()) + (1/2)*Z2.^2;

% subplot of planars
figure(1)
clf
pos1 = [0.1 0.6 0.35 0.35];
subplot('Position', pos1), plotHarcColor(xSlow(:,1),jSlow,xSlow(:,2),tSlow);
xlabel('$\mathrm{z_1}$','FontSize',16)
ylabel('$\mathrm{z_2}$','FontSize',16)
contour(Z1,Z2,V0_0,[c_0 c_0],'-r','ShowText','on')
contour(Z1,Z2,V1_1,[c_10 c_10],'-g','ShowText','on')
plot(z1deltaSlow,z2deltaSlow,'r.','MarkerSize', 14)
strDeltaSlow = [num2str(timeToDeltaSlow),' s'];
text(z1deltaSlow,z2deltaSlow,strDeltaSlow,'HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',12);
axis([-11 6 -0.5 6]);
grid on
pos2 = [0.55 0.6 0.35 0.35];
subplot('Position', pos2), plotHarcColor(xOscillate(:,1),jOscillate,xOscillate(:,2),tOscillate);
xlabel('$\mathrm{z_1}$','FontSize',16)
ylabel('$\mathrm{z_2}$','FontSize',16)
contour(Z1,Z2,V0_0,[c_0 c_0],'-r','ShowText','on')
contour(Z1,Z2,V1_1,[c_10 c_10],'-g','ShowText','on')
plot(z1deltaOscillate,z2deltaOscillate,'r.','MarkerSize', 14)
strDeltaOscillate = [num2str(timeToDeltaOscillate),' s'];
text(z1deltaOscillate,z2deltaOscillate,strDeltaOscillate,'HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',12);
axis([-11 6 -4 6]);
grid on
pos3 = [0.33 0.15 0.35 0.35];
subplot('Position', pos3), plotHarcColor(xUniting(:,1),jUniting,xUniting(:,2),tUniting);
xlabel('$\mathrm{z_1}$','FontSize',16)
ylabel('$\mathrm{z_2}$','FontSize',16)
contour(Z1,Z2,V0_0,[c_0 c_0],'-r','ShowText','on')
contour(Z1,Z2,V1_1,[c_10 c_10],'-g','ShowText','on')
plot(z1deltaUniting,z2deltaUniting,'r.','MarkerSize', 14)
strDeltaUniting = [num2str(timeToDeltaUniting),' s'];
text(z1deltaUniting,z2deltaUniting,strDeltaUniting,'HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',12);
axis([-11 6 -0.5 6]);
grid on
% h = gcf;
% pos = h.PaperPosition;
% h.PaperPositionMode = 'auto';
% set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(h,'Plots\MultiPlane','-dpdf','-r0')
saveas(gcf,'Plots\MultiPlane','epsc')

% subplots of x_1 and x_2, on same plot
figure(2) 
clf
modificatorF{1} = '';
modificatorF{2} = 'LineWidth';
modificatorF{3} = 1.5;
modificatorJ{1} = '*--';
modificatorJ{2} = 'LineWidth';
modificatorJ{3} = 1.5;

subplot(2,1,1), plotHarc(ta,ja,xa,[],modificatorF,modificatorJ);
hold on
plot(timeToDeltaSlow,z1deltaSlow,'k.','MarkerSize', 14)
plot(timeToDeltaOscillate,z1deltaOscillate,'k.','MarkerSize', 14)
plot(timeToDeltaUniting,z1deltaUniting,'k.','MarkerSize', 14)
strDeltaSlow = [num2str(timeToDeltaSlow),' s'];
strDeltaOscillate = [num2str(timeToDeltaOscillate),' s'];
strDeltaUniting = [num2str(timeToDeltaUniting),' s'];
text(timeToDeltaSlow,z1deltaSlow,strDeltaSlow,'HorizontalAlignment','right','VerticalAlignment','top','FontSize',12);
text(timeToDeltaOscillate,z1deltaOscillate,strDeltaOscillate,'HorizontalAlignment','right','VerticalAlignment','top','FontSize',12);
text(timeToDeltaUniting,z1deltaUniting,strDeltaUniting,'HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',12);
grid on
hold off
axis([0 330 -10 7]);
ylabel('$\mathrm{z_1}$','FontSize',16)
subplot(2,1,2), plotHarc(ta,ja,xb,[],modificatorF,modificatorJ);
hold on
plot(timeToDeltaSlow,z2deltaSlow,'k.','MarkerSize', 14)
plot(timeToDeltaOscillate,z2deltaOscillate,'k.','MarkerSize', 14)
plot(timeToDeltaUniting,z2deltaUniting,'k.','MarkerSize', 14)
grid on
hold off
axis([0 330 -3 5]);
ylabel('$\mathrm{z_2}$','FontSize',16)
% h = gcf;
% pos = h.PaperPosition;
% h.PaperPositionMode = 'auto';
% set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(h,'Plots\Multiplots','-dpdf','-r0')
saveas(gcf,'Plots\Multiplots','epsc')