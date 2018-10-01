%--------------------------------------------------------------------------
% Matlab M-file Project: HyEQ Toolbox @  Hybrid Systems Laboratory (HSL), 
% https://hybrid.soe.ucsc.edu/software
% http://hybridsimulator.wordpress.com/
% Filename: HeavyBall.m
%--------------------------------------------------------------------------
% Project: Unifying local and global control, using the Heavy-Ball Method
% for convergence to a global minimum. Different parameters for lambda and
% gamma are used globally and locally. 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%   See also HYEQSOLVER, PLOTARC, PLOTARC3, PLOTFLOWS, PLOTHARC,
%   PLOTHARCCOLOR, PLOTHARCCOLOR3D, PLOTHYBRIDARC, PLOTJUMPS.
%   Copyright @ Hybrid Systems Laboratory (HSL),
%   Revision: 0.0.0.3 Date: 05/20/2015 3:42:00

clear all

set(0,'defaultTextInterpreter','latex');

% global variables
global gamma_0 lambda_0 gamma_1 lambda_1 c_0 c_10 V0 V1 delta z1Star sigma randomsInterp randomsIndex
%%%%%%%%%%%%%%%%%%%%%%%
% setting the globals %
%%%%%%%%%%%%%%%%%%%%%%%

z1Star = 0; %3;

% Heavy-ball constants: Need to play around with these!
lambda_0 = 10.5; % Gravity. 
            % For gamma fixed, "large values of  lambda are seen to give rise to slowly converging 
            % solutions resembling the steepest descent’s while smaller values give 
            % rise to fast solutions with oscillations getting wilder as lambda decreases."
gamma_0 = 1/2; % Viscous friction to mass ratio.

lambda_1 = 1/5;

gamma_1 = 1/2;

% Temporarily using a fixed constant for the level sets. Will find better
% solution soon.
c_0 = 12.5;%31; % \mathcal{U}_0
c_10 = 6.3;%30; % \mathcal{T}_{1,0}

delta = 0.01;
deltaNoise = 0.01;%0.01;

V0 = 0;
V1 = 0;

%%%%%%%%%%%%%%%%%%%%%%
% setting the locals %
%%%%%%%%%%%%%%%%%%%%%%

timeToDelta = 0;
timeToDeltaIdx = 1;
z1delta = 0;
z2delta = 0;
jumpsx1 = [];
jumpsx2 = [];
jumpst = [];
jumpIndex = 1;

timeToDeltab = 0;
timeToDeltaIdxb = 1;
z1deltab = 0;
z2deltab = 0;
jumpsx1b = [];
jumpsx2b = [];
jumpstb = [];
jumpIndexb = 1;

timeToDeltac = 0;
timeToDeltaIdxc = 1;
z1deltac = 0;
z2deltac = 0;
jumpsx1c = [];
jumpsx2c = [];
jumpstc = [];
jumpIndexc = 1;


timeToDeltaNom = 0;
timeToDeltaIdxNom = 1;
z1deltaNom = 0;
z2deltaNom = 0;

% initial conditions
z1_0 = -10;%-9.79;
z2_0 = 0;
q_0 = 0;
% Assign initial conditions to vector
x0 = [z1_0;z2_0;q_0];

% simulation horizon
TSPAN=[0 200];
JSPAN = [0 10];

% rule for jumps
% rule = 1 -> priority for jumps
% rule = 2 -> priority for flows
rule = 1;

options = odeset('RelTol',1e-6,'MaxStep',.1);

% simulate nominal system
[tNom,jNom,xNom] = HyEQsolver(@fNom,@gNom,@CNom,@DNom,...
    x0,TSPAN,JSPAN,rule,options);

nomLen = 100*length(xNom);
sample = 0.1;

% simulate the perturbed system: sigma = 0.1
randoms = randn(1,nomLen);
randInd = 1:nomLen;
randomsInterpIndex = 1:sample:nomLen;
randomsInterp = interp1(randInd,randoms,randomsInterpIndex);
randomsIndex = 1;
sigma = 0.1;
% simulate the perturbed system
[ta,ja,xa] = HyEQsolver(@f,@g,@C,@D,...
    x0,TSPAN,JSPAN,rule,options);


% simulate the perturbed system: sigma = 0.5
% rng('shuffle');
% randoms = randn(1,nomLen);
% randInd = 1:nomLen;
% randomsInterpIndex = 1:sample:nomLen;
% randomsInterp = interp1(randInd,randoms,randomsInterpIndex);
randomsIndex = 1;
sigma = 0.5;
% simulate the perturbed system
[tb,jb,xb] = HyEQsolver(@f,@g,@C,@D,...
    x0,TSPAN,JSPAN,rule,options);

% simulate the perturbed system: sigma = 1
% rng('shuffle');
% randoms = randn(1,nomLen);
% randInd = 1:nomLen;
% randomsInterpIndex = 1:sample:nomLen;
% randomsInterp = interp1(randInd,randoms,randomsInterpIndex);
randomsIndex = 1;
sigma = 1;
% simulate the perturbed system
[tc,jc,xc] = HyEQsolver(@f,@g,@C,@D,...
    x0,TSPAN,JSPAN,rule,options);



% Finding time of convergence for the perturbed system: sigma = 0.1
for i=2:length(xa(:,1))
    if (abs(z1Star - xa(i,1)) <= deltaNoise) && (abs(z1Star - xa(i-1,1)) > deltaNoise)
        timeToDeltaIdx = i;
        z1delta = xa(i,1);
    end
end

z2delta = xa(timeToDeltaIdx,2);
timeToDelta = ta(timeToDeltaIdx,1);

% Find the jumps for sigma = 0.1:
for i=2:length(ja)
    if(ja(i,1) ~= ja(i-1,1))
        jumpsx1(jumpIndex) = xa(i,1);
        jumpsx2(jumpIndex) = xa(i,2);
        jumpst(jumpIndex) = ta(i,1);
        jumpIndex = jumpIndex + 1;
    end
end


% Finding time of convergence for the perturbed system: sigma = 0.5
for i=2:length(xb(:,1))
    if (abs(z1Star - xb(i,1)) <= deltaNoise) && (abs(z1Star - xb(i-1,1)) > deltaNoise)
        timeToDeltaIdxb = i;
        z1deltab = xb(i,1);
    end
end

z2deltab = xb(timeToDeltaIdxb,2);
timeToDeltab = tb(timeToDeltaIdxb,1);

% Find the jumps for sigma = 0.5:
for i=2:length(jb)
    if(jb(i,1) ~= jb(i-1,1))
        jumpsx1b(jumpIndexb) = xb(i,1);
        jumpsx2b(jumpIndexb) = xb(i,2);
        jumpstb(jumpIndexb) = tb(i,1);
        jumpIndexb = jumpIndexb + 1;
    end
end

% Finding time of convergence for the perturbed system: sigma = 1
for i=2:length(xc(:,1))
    if (abs(z1Star - xc(i,1)) <= deltaNoise) && (abs(z1Star - xc(i-1,1)) > deltaNoise)
        timeToDeltaIdxc = i;
        z1deltac = xc(i,1);
    end
end

z2deltac = xc(timeToDeltaIdxc,2);
timeToDeltac = tc(timeToDeltaIdxc,1);

% Find the jumps for sigma = 0.5:
for i=2:length(jc)
    if(jc(i,1) ~= jc(i-1,1))
        jumpsx1c(jumpIndexc) = xc(i,1);
        jumpsx2c(jumpIndexc) = xc(i,2);
        jumpstc(jumpIndexc) = tc(i,1);
        jumpIndexc = jumpIndexc + 1;
    end
end

% Finding time of convergence for the nominal system
for i=2:length(xNom(:,1))
    if (abs(z1Star - xNom(i,1)) <= delta) && (abs(z1Star - xNom(i-1,1)) > delta)
        timeToDeltaIdxNom = i;
        z1deltaNom = xNom(i,1);
    end
end

z2deltaNom = xNom(timeToDeltaIdxNom,2);
timeToDeltaNom = tNom(timeToDeltaIdxNom,1);

x1 = -11:0.005:6;
x2 = -1:0.005:8;
[Z1,Z2] = meshgrid(x1,x2);
V0_0 = gamma_0*(0.25*(Z1-z1Star).^2 - CalculateLStar()) + (1/2)*Z2.^2;
V1_1 = gamma_1*(0.25*(Z1-z1Star).^2 - CalculateLStar()) + (1/2)*Z2.^2;

% plot phase plane
figure(1) % position
clf
pos1 = [0.1 0.6 0.35 0.35];
subplot('Position', pos1), plotHarc(xa(:,1),ja,xa(:,2));
hold on
xlabel('$\mathrm{z_1}$','FontSize',16)
ylabel('$\mathrm{z_2}$','FontSize',16)
contour(Z1,Z2,V0_0,[c_0 c_0],'-r','ShowText','on')
contour(Z1,Z2,V1_1,[c_10 c_10],'-g','ShowText','on')
plot(z1delta,z2delta,'r.','MarkerSize', 14)
strDelta = [num2str(timeToDelta), 's'];
text(z1delta,z2delta,strDelta,'HorizontalAlignment','left','VerticalAlignment','bottom');
axis([-11 6 -0.5 6]);
% axis([-6 -5.95 3.44 3.48]);
grid on

pos2 = [0.55 0.6 0.35 0.35];
subplot('Position', pos2), plotHarc(xb(:,1),jb,xb(:,2));
hold on
xlabel('$\mathrm{z_1}$','FontSize',16)
ylabel('$\mathrm{z_2}$','FontSize',16)
contour(Z1,Z2,V0_0,[c_0 c_0],'-r','ShowText','on')
contour(Z1,Z2,V1_1,[c_10 c_10],'-g','ShowText','on')
plot(z1deltab,z2deltab,'r.','MarkerSize', 14)
strDeltab = [num2str(timeToDeltab), 's'];
text(z1deltab,z2deltab,strDeltab,'HorizontalAlignment','left','VerticalAlignment','bottom');
axis([-11 6 -0.5 6]);
grid on

pos3 = [0.1 0.15 0.35 0.35];
subplot('Position', pos3), plotHarc(xc(:,1),jc,xc(:,2));
hold on
xlabel('$\mathrm{z_1}$','FontSize',16)
ylabel('$\mathrm{z_2}$','FontSize',16)
contour(Z1,Z2,V0_0,[c_0 c_0],'-r','ShowText','on')
contour(Z1,Z2,V1_1,[c_10 c_10],'-g','ShowText','on')
plot(z1deltac,z2deltac,'r.','MarkerSize', 14)
strDeltac = [num2str(timeToDeltac), 's'];
text(z1deltac,z2deltac,strDeltac,'HorizontalAlignment','left','VerticalAlignment','bottom');
axis([-11 6 -0.5 6]);
grid on

pos4 = [0.55 0.15 0.35 0.35];
subplot('Position', pos4), plotHarc(xNom(:,1),jNom,xNom(:,2));
hold on
xlabel('$\mathrm{z_1}$','FontSize',16)
ylabel('$\mathrm{z_2}$','FontSize',16)
%title('Nominal System')
contour(Z1,Z2,V0_0,[c_0 c_0],'-r','ShowText','on')
contour(Z1,Z2,V1_1,[c_10 c_10],'-g','ShowText','on')
plot(z1deltaNom,z2deltaNom,'r.','MarkerSize', 14)
strDeltaNom = [num2str(timeToDeltaNom),'s'];
text(z1deltaNom,z2deltaNom,strDeltaNom,'HorizontalAlignment','left','VerticalAlignment','bottom');
axis([-11 6 -0.5 6]);
grid on
h = gcf;
pos = h.PaperPosition;
h.PaperPositionMode = 'auto';
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'Plots\PhasePlaneMulti','-dpdf','-r0')

% Attempt to plot results on same figure:
figure(2)
clf
plotHarc(xNom(:,1),jNom,xNom(:,2));
hold on
plot(xa(:,1),xa(:,2),'-k');
plot(jumpsx1(1),jumpsx2(1),'*k');
plot(jumpsx1(2),jumpsx2(2),'*k');

plot(xb(:,1),xb(:,2),'-m');
plot(jumpsx1b(1),jumpsx2b(1),'*m');
plot(jumpsx1b(2),jumpsx2b(2),'*m');

plot(xc(:,1),xc(:,2),'-c');
plot(jumpsx1c(1),jumpsx2c(1),'*c');
plot(jumpsx1c(2),jumpsx2c(2),'*c');
xlabel('$\mathrm{z_1}$','FontSize',16)
ylabel('$\mathrm{z_2}$','FontSize',16)
%title('Nominal System')
contour(Z1,Z2,V0_0,[c_0 c_0],'-r','ShowText','on')
contour(Z1,Z2,V1_1,[c_10 c_10],'-g','ShowText','on')
plot(z1deltaNom,z2deltaNom,'r.','MarkerSize', 14)
strDeltaNom = [num2str(timeToDeltaNom),'s'];
text(z1deltaNom,z2deltaNom,strDeltaNom,'HorizontalAlignment','left','VerticalAlignment','bottom');
plot(z1delta,z2delta,'r.','MarkerSize', 14)
strDelta = [num2str(timeToDelta), 's'];
text(z1delta,z2delta,strDelta,'HorizontalAlignment','left','VerticalAlignment','top');
plot(z1deltab,z2deltab,'r.','MarkerSize', 14)
strDeltab = [num2str(timeToDeltab), 's'];
text(z1deltab,z2deltab,strDeltab,'HorizontalAlignment','right','VerticalAlignment','bottom');
plot(z1deltac,z2deltac,'r.','MarkerSize', 14)
strDeltac = [num2str(timeToDeltac), 's'];
text(z1deltac,z2deltac,strDeltac,'HorizontalAlignment','right','VerticalAlignment','top');
axis([-11 6 -0.5 6]);
grid on
h = gcf;
pos = h.PaperPosition;
h.PaperPositionMode = 'auto';
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
saveas(gcf,'Plots\PhasePlaneNomAll','epsc')