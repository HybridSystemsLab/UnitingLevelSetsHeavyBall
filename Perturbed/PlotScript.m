% Rerun on existing .mat file, to make adjustments to plot.
figure(1)
clf
modificatorF{1} = '';
modificatorF{2} = 'LineWidth';
modificatorF{3} = 1.5;
modificatorJ{1} = '*--';
modificatorJ{2} = 'LineWidth';
modificatorJ{3} = 1.5;
pos1 = [0.1 0.15 0.45 0.75];
subplot('Position', pos1),plotHarc(xNom(:,1),jNom,xNom(:,2),[],modificatorF,modificatorJ);
hold on
plot(xa(:,1),xa(:,2),'Color',[0 .5 .0],'LineWidth',1.5);
plot(jumpsx1(1),jumpsx2(1),'*','Color',[0 .5 .0]);
plot(jumpsx1(2),jumpsx2(2),'*','Color',[0 .5 .0]);

plot(xb(:,1),xb(:,2),'-r','LineWidth',1.5);
plot(jumpsx1b(1),jumpsx2b(1),'*r');
plot(jumpsx1b(2),jumpsx2b(2),'*r');

plot(xc(:,1),xc(:,2),'Color',[.44 0 .44],'LineWidth',1.5);
plot(jumpsx1c(1),jumpsx2c(1),'*','Color',[.44 0 .44]);
plot(jumpsx1c(2),jumpsx2c(2),'*','Color',[.44 0 .44]);
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
axis([-11 3 -1 5.5]);
grid on
pos2 = [0.65 0.15 0.25 0.75];
subplot('Position', pos2),plotHarc(xNom(:,1),jNom,xNom(:,2),[],modificatorF,modificatorJ);
hold on
plot(xa(:,1),xa(:,2),'Color',[0 .5 .0],'LineWidth',1.5);
plot(jumpsx1(1),jumpsx2(1),'*','Color',[0 .5 .0]);
plot(jumpsx1(2),jumpsx2(2),'*','Color',[0 .5 .0]);

plot(xb(:,1),xb(:,2),'-r','LineWidth',1.5);
plot(jumpsx1b(1),jumpsx2b(1),'*r');
plot(jumpsx1b(2),jumpsx2b(2),'*r');

plot(xc(:,1),xc(:,2),'Color',[.44 0 .44],'LineWidth',1.5);
plot(jumpsx1c(1),jumpsx2c(1),'*','Color',[.44 0 .44]);
plot(jumpsx1c(2),jumpsx2c(2),'*','Color',[.44 0 .44]);
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
axis([-0.45 0.25 -0.5 4]);
grid on
saveas(gcf,'Plots\PhasePlaneNomAllCloseup','epsc')