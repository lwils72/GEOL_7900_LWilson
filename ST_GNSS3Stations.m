clear;
% Clear screen
a=gcf;
clf(a);

% Download data directly from the source

!curl http://geodesy.unr.edu/gps_timeseries/tenv3/plates/NA/P395.NA.tenv3 > P395.tenv3.txt;
!curl http://geodesy.unr.edu/gps_timeseries/tenv3/plates/NA/P396.NA.tenv3 > P396.tenv3.txt;
!curl http://geodesy.unr.edu/gps_timeseries/tenv3/plates/NA/P404.NA.tenv3 > P404.tenv3.txt;

% File content 

fidP395 ='Users/leora/Desktop/LSU23/MATLAB/P395.tenv3.txt';
fidP396 ='Users/leora/Desktop/LSU23/MATLAB/P396.tenv3.txt';
fidP404 ='Users/leora/Desktop/LSU23/MATLAB/P404.tenv3.txt';

A=readmatrix(fidP395);
B=readmatrix(fidP396);
C=readmatrix(fidP404);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%organize what columns are what 
P395t=A(:,3); % time - date
P395x=A(:,9); % eastings
P395y=A(:,11); % northings
P395z=A(:,13); % vertical

station395lat=A(:,21);
stationP395lat=station395lat(1);
station395long=A(:,22);
stationP395long=station395long(1);

[P395a,P395b,P395c]=ll2utm(station395lat,station395long);

P396t=B(:,3); % time - date
P396x=B(:,9); % eastings
P396y=B(:,11); % northings
P396z=B(:,13); % vertical

station396lat=B(:,21);
stationP396lat=station396lat(1);
station396long=B(:,22);
stationP396long=station396long(1);

[P396a,P396b,P396c]=ll2utm(station396lat,station396long);

P404t=C(:,3); % time - date
P404x=C(:,9); % eastings
P404y=C(:,11); % northings
P404z=C(:,13); % vertical

station404lat=C(:,21);
stationP404lat=station404lat(1);
station404long=C(:,22);
stationP404long=station404long(1);

[P404a,P404b,P404c]=ll2utm(station404lat,station404long);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pP395E = polyfit(P395t,P395a,1); %fits plots for eastings
vP395E = pP395E(1);
pP396E = polyfit(P396t,P396a,1); %fits plots for eastings
vP396E = pP396E(1);
pP404E = polyfit(P404t,P404a,1); %fits plots for eastings
vP404E = pP404E(1);

pP395N = polyfit(P395t,P395b,1); %fits plots for northings
vP395N = pP395N(1);
pP396N = polyfit(P396t,P396b,1); %fits plots for northings
vP396N = pP396N(1);
pP404N = polyfit(P404t,P404b,1); %fits plots for northings
vP404N = pP404N(1);

pP395Z = polyfit(P395t,P395z,1); %fits plots for vertical
vP395Z = pP395Z(1);
pP396Z = polyfit(P396t,P396z,1); %fits plots for vertical
vP396Z = pP396Z(1);
pP404Z = polyfit(P404t,P404z,1); %fits plots for vertical
vP404Z = pP404Z(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% estimated vx and vy values at each station

d = [0.0069;0.0106;0.0071;0.0107;0.0052;0.0094];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate the centroid or center of mass of the triangle

meanx1 = sum(P395a)/6148;
meany1 = sum(P395b)/6148;
meanx2 = sum(P396a)/5640;
meany2 = sum(P396b)/5640;
meanx3 = sum(P404a)/3606;
meany3 = sum(P404b)/3606;

% calculate the distance between each station and the centroid

dx1 = P395a-meanx1 ;
dxP395 = dx1(1);
dy1 = P395b-meany1 ;
dyP395 = dy1(1);

dx2 = P396a-meanx2 ;
dxP396 = dx2(1);
dy2 = P396b-meany2 ;
dyP396 = dy2(1);

dx3 = P404a-meanx3 ;
dxP404 = dx3(1);
dy3 = P404b-meany3 ;
dyP404 = dy3(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% g matrix 

g = [1 , 0 ,dxP395 , dyP395 , 0 , -dyP395;...
    0 , 1 , 0 , dxP395 , dyP395 , dxP395;...
    1 , 0 ,dxP396 , dyP396 , 0 , -dyP396;...
    0 , 1 , 0 , dxP396 , dyP396 , dxP396;...
    1 , 0 ,dxP404 , dyP404 , 0 , -dyP404;...
    0 , 1 , 0 , dxP404 , dyP404 , dxP404];

m = g \ d; 

E = [m(3),m(4);m(4),m(5)];

[paxes,pvalues] = eig(E)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plotting the stations and the strain tensor 
 
figure(1),clf
hold on;

plot(P395a,P395b,'^k'); % plotting stations
plot(P396a,P396b,'^k');
plot(P404a,P404b,'^k');

quiver(P395a(1),P395b(1),vP395E,vP395N,1e2); % plot station velocities 
quiver(P396a(1),P396b(1),vP396E,vP396N,1e2);
quiver(P404a(1),P404b(1),vP404E,vP404N,1e2);

scatter(P395a(1),P395b(1),50,vP395Z*1000,'filled','g'); % plot stations based on vertical motion
scatter(P396a(1),P396b(1),50,vP396Z*1000,'filled','g');
scatter(P404a(1),P404b(1),50,vP404Z*1000,'filled','g');
colorbar
colormap(jet);
xlim([4.3e5,4.7e5]);
ylim([4.98e6,5.03e6]);

%plot()
%plot(E(1,2),E(2,1),'-','linewidth',2);







