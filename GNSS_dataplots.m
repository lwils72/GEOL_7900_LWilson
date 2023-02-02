% Clean Memory
clear;
% Clear screen
a=gcf;
clf(a);

% Download data directly from the source

!curl http://geodesy.unr.edu/gps_timeseries/tenv3/plates/NA/P395.NA.tenv3 > P395.tenv3.txt;
!curl http://geodesy.unr.edu/gps_timeseries/tenv3/plates/NA/P396.NA.tenv3 > P396.tenv3.txt;
!curl http://geodesy.unr.edu/gps_timeseries/tenv3/plates/NA/P404.NA.tenv3 > P404.tenv3.txt;

% File content 
% Either using textscan or fscanf 

fidP395 ='Users/leora/Desktop/LSU23/MATLAB/P395.tenv3.txt';
fidP396 ='Users/leora/Desktop/LSU23/MATLAB/P396.tenv3.txt';
fidP404 ='Users/leora/Desktop/LSU23/MATLAB/P404.tenv3.txt';

A=readmatrix(fidP395);
B=readmatrix(fidP396);
C=readmatrix(fidP404);

%format='%s %s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f','headerlines',1;
%C=textscan(fid395,format);
%fclose(fid395);

%organize what columns are what 
P395t=A(:,3); % time - date
P395x=A(:,9); % eastings
P395y=A(:,11); % northings
P395z=A(:,13); % vertical

station395lat=A(:,21);
stationP395lat=station395lat(1);
station395long=A(:,22);
stationP395long=station395long(1);

P396t=B(:,3); % time - date
P396x=B(:,9); % eastings
P396y=B(:,11); % northings
P396z=B(:,13); % vertical

station396lat=B(:,21);
stationP396lat=station396lat(1);
station396long=B(:,22);
stationP396long=station396long(1);

P404t=C(:,3); % time - date
P404x=C(:,9); % eastings
P404y=C(:,11); % northings
P404z=C(:,13); % vertical

station404lat=C(:,21);
stationP404lat=station404lat(1);
station404long=C(:,22);
stationP404long=station404long(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%plot the figures as seen from GNSS

% P395

figure(1);
hold on;
set(gca)

subplot(311);
title('P395 Station Eastings, Northings, Vertical');
hold on;
pP395E = polyfit(P395t,P395x,1); %fits plots for eastings
vP395E = pP395E(1);
plot(P395t,P395x,'y');
plot(P395t,polyval(pP395E,P395t),'linewidth',3);
grid,xlabel ('Time'),ylabel ('eastings');

subplot(312);
hold on;
pP395N = polyfit(P395t,P395y,1); %fits plots for northings
vP395N = pP395N(1);
plot(P395t,P395y,'r')
plot(P395t,polyval(pP395N,P395t),'linewidth',3);
grid,xlabel ('Time'),ylabel ('northings');

subplot(313);
hold on;
pP395Z = polyfit(P395t,P395z,1); %fits plots for northings
vP395Z = pP395Z(1);
plot(P395t,P395z,'r');
plot(P395t,polyval(pP395Z,P395t),'linewidth',3);
grid,xlabel ('Time'),ylabel ('vertical');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% P396

figure(2);
hold on;
set(gca)

subplot(311);
title('P396 Station Eastings, Northings, Vertical');
pP396E=polyfit(P396t,P396x,1); %fits plots for eastings
vP396E=pP396E(1);
plot(P396t,P396x,'r');
plot(P396t,polyval(pP396E,P396t),'linewidth',3);
grid,xlabel ('Time'),ylabel ('eastings');

subplot(312);
pP396N=polyfit(P396t,P396y,1); %fits plots for northings
vP396N=pP396N(1);
plot(P396t,P396y,'r');
plot(P396t,polyval(pP396N,P396t),'linewidth',3);
grid,xlabel ('Time'),ylabel ('northings');

subplot(313);
pP396Z=polyfit(P396t,P396z,1); %fits plots for vertical
vP396Z=pP396Z(1);
plot(P396t,P396z,'r');
plot(P396t,polyval(pP396Z,P396t),'linewidth',3);
grid,xlabel ('Time'),ylabel ('vertical');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% P404

figure(3);
hold on;
set(gca)

subplot(311);
title('P395 Station Eastings, Northings, Vertical');
pP404E = polyfit(P404t,P404x,1); %fits plots for eastings
vP404E = pP404E(1);
plot(P404t,P404x,'b'); 
plot(P404t,polyval(pP404E,P404t),'linewidth',3);
grid,xlabel ('Time'),ylabel ('eastings');

subplot(312);
pP404N = polyfit(P404t,P404y,1); %fits plots for northings
vP404N = pP404N(1);
plot(P404t,P404y,'b');
plot(P404t,polyval(pP404N,P404t),'linewidth',3);
grid,xlabel ('Time'),ylabel ('northings');

subplot(313);
pP404Z = polyfit(P404t,P404z,1); %fits plots for vertical
vP404Z = pP404Z(1);
plot(P404t,P404z,'r');
plot(P404t,polyval(pP404Z,P404t),'linewidth',3);
grid,xlabel ('Time'),ylabel ('vertical');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plotting coast line with political boundaries 

!curl https://raw.githubusercontent.com/lsugeo/GEOL_7900_Demo/main/coastfile.xy > coastfile.xy
!curl https://raw.githubusercontent.com/lsugeo/GEOL_7900_Demo/main/politicalboundaryfile.xy > politicalboundaryfile.xy

c=load('coastfile.xy');
b=load('politicalboundaryfile.xy');
 
figure(4),clf
plot(c(:,1),c(:,2));
hold on;
plot(b(:,1),b(:,2),'color',[1 1 1]*0.75);

plot(stationP395long,stationP395lat,'^k');
plot(stationP396long,stationP396lat,'^k');
plot(stationP404long,stationP404lat,'^k');

quiver(stationP395long,stationP395lat,vP395E,vP395N,1e2);
quiver(stationP396long,stationP396lat,vP396E,vP396N,1e2);
quiver(stationP404long,stationP404lat,vP404E,vP404N,1e2);

scatter(stationP395long,stationP395lat,50,vP395Z*1000,'filled');
scatter(stationP396long,stationP396lat,50,vP396Z*1000,'filled');
scatter(stationP404long,stationP404lat,50,vP404Z*1000,'filled');

colorbar
clim([-1,1]);
colormap(jet);
xlim([-127,-115]);
ylim([41,50]);


