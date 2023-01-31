% Clean Memory
clear;
% Clear screen
a=gcf;
clf(a);

% Download data directly from the source

!curl http://geodesy.unr.edu/gps_timeseries/tenv3/IGS14/P395.tenv3 > P395.tenv3.txt;
!curl http://geodesy.unr.edu/gps_timeseries/tenv3/IGS14/P396.tenv3 > P396.tenv3.txt;
!curl http://geodesy.unr.edu/gps_timeseries/tenv3/IGS14/P404.tenv3 > P404.tenv3.txt;

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

%P395t=C{3};
%P395x=C{9};
%P395y=C{11};
%P395z=C{13};

%organize what columns are what 
P395t=A(:,3); % time - date
P395x=A(:,9); % eastings
P395y=A(:,11); % northings
P395z=A(:,13); % vertical

P396t=B(:,3); % time - date
P396x=B(:,9); % eastings
P396y=B(:,11); % northings
P396z=B(:,13); % vertical

P404t=B(:,3); % time - date
P404x=B(:,9); % eastings
P404y=B(:,11); % northings
P404z=B(:,13); % vertical

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%plot the figures as seen from GNSS

% P395

figure(1);
hold on;
set(gca)

subplot(311);
plot(P395t,P395x,'b'); 
title('P395 Station Eastings, Northings, Vertical');
grid,xlabel ('Time'),ylabel ('eastings');

subplot(312);
plot(P395t,P395y,'b');
grid,xlabel ('Time'),ylabel ('northings');

subplot(313);
plot(P395t,P395z,'r');
grid,xlabel ('Time'),ylabel ('vertical');

polyfit(P395t,P395x,1) %fits plots for eastings
polyfit(P395t,P395y,1) %fits plots for northings
polyfit(P395t,P395z,1) %fits plots for vertical

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% P396

figure(2);
hold on;
set(gca)

subplot(311);
plot(P396t,P396x,'b'); 
title('P395 Station Eastings, Northings, Vertical');
grid,xlabel ('Time'),ylabel ('eastings');

subplot(312);
plot(P396t,P396y,'b');
grid,xlabel ('Time'),ylabel ('northings');

subplot(313);
plot(P396t,P396z,'r');
grid,xlabel ('Time'),ylabel ('vertical');

polyfit(P396t,P396x,1) %fits plots for eastings
polyfit(P396t,P396y,1) %fits plots for northings
polyfit(P396t,P396z,1) %fits plots for vertical

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% P404

figure(3);
hold on;
set(gca)

subplot(311);
%pP404E = polyfit(P404t,P404x,1); %fits plots for eastings
%fP404E = polyval(pP404E,P404x);
plot(P404t,P404x,'b',P404t,fP404E,'-'); 
title('P395 Station Eastings, Northings, Vertical');
grid,xlabel ('Time'),ylabel ('eastings');

subplot(312);
%pP404N = polyfit(P404t,P404y,1); %fits plots for northings
%fP404N = polyval(pP404N,P404y);
plot(P404t,P404y,'b',P404t,fP404N,'-');
grid,xlabel ('Time'),ylabel ('northings');

subplot(313);
%pP404V = polyfit(P404t,P404z,1); %fits plots for vertical
%fP404V = polyval(pP404V,P404z);
plot(P404t,P404z,'r',P404t,fP404V,'-');
grid,xlabel ('Time'),ylabel ('vertical');

