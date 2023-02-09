% Clean Memory
clear;
% Clear screen
a=gcf;
clf(a);

% Download data directly from the source

!curl http://geodesy.unr.edu/gps_timeseries/tenv3/plates/NA/LAIB.NA.tenv3 > LAIB.tenv3.txt;

% File content 

fidLAIB =fopen('Users/leora/Desktop/LSU23/MATLAB/LAIB.tenv3.txt');
A=textscan(fidLAIB,'%s %s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f','headerlines',1);
fclose(fidLAIB);

%organize what columns are what 

LAIBtd=A{3}; % time - date
LAIBx=A{9}; % eastings
LAIBy=A{11}; % northings
LAIBz=A{13}; % vertical

LAIBt=datenum(A{2},'yyyymmmdd'); % time in "Matlab time", sequential days

stationLAIBLat=A{21}(1);
stationLAIBLong=A{22}(1);

%calculate best fit lines for stations

pLAIBE = polyfit(LAIBt,LAIBx,1); %fits plots for eastings
vLAIBE = pLAIBE(1);
pLAIBN = polyfit(LAIBt,LAIBy,1); %fits plots for northings
vLAIBN = pLAIBN(1);
pLAIBZ = polyfit(LAIBt,LAIBz,1); %fits plots for northings
vLAIBZ = pLAIBZ(1);

%plot the subplots of original data for the eastings, northings, vertical 

figure(1);
hold on;
set(gca)

subplot(311);
title('LAIB Station Eastings, Northings, Vertical');
hold on;
plot(LAIBt,LAIBx,'g');
plot(LAIBt,polyval(pLAIBE,LAIBt),'linewidth',1);
grid,xlabel ('Time'),ylabel ('eastings');
datetick;

subplot(312);
hold on;
plot(LAIBt,LAIBy,'r')
plot(LAIBt,polyval(pLAIBN,LAIBt),'linewidth',1);
grid,xlabel ('Time'),ylabel ('northings');
datetick;

subplot(313);
hold on;
plot(LAIBt,LAIBz,'b');
plot(LAIBt,polyval(pLAIBZ,LAIBt),'linewidth',1);
grid,xlabel ('Time'),ylabel ('vertical');
datetick;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clean data 

dt=1;

[t,x]=filltimegap(LAIBt,LAIBx,dt);
[t,y]=filltimegap(LAIBt,LAIBy,dt);
[t,z]=filltimegap(LAIBt,LAIBz,dt);

LAIBtNaN=t;
LAIBxNaN=x;
LAIByNaN=y;
LAIBzNaN=z; 

% I don't think that there are NaNs inserted into the code

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fix the jump

tjump=datenum([0017 08 04]); % jump happens on date 2017 08 02
iafter=find(LAIBt>tjump); % all values after 1429 and 2.017593400000000e+03
jumpLAIBx= 6.73;
jumpLAIBy= 6.65;
jumpLAIBz= 1.5;
LAIBx(iafter)=LAIBx(iafter)-jumpLAIBx;
LAIBy(iafter)=LAIBy(iafter)+jumpLAIBy;
LAIBz(iafter)=LAIBz(iafter)+jumpLAIBz;

%tsjump=datenum([0017 10 24]); % another jump in eastings at oct 24 2017
%i2after=find(LAIBt>tsjump); % second jump
%jump2LAIBx= 0.7;
%LAIBx(i2after)=LAIBx(i2after)-jump2LAIBx;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot to show the new adjusted jump

figure(3);
hold on;
set(gca)

subplot(311);
title('LAIB Station adjusted for jump');
hold on;
plot(LAIBt,LAIBx,'g');
%plot(LAIBt,polyval(pLAIBE,LAIBt),'linewidth',1);
grid,xlabel ('Time'),ylabel ('eastings');
datetick;

subplot(312);
hold on;
plot(LAIBt,LAIBy,'r')
%plot(LAIBt,polyval(pLAIBN,LAIBt),'linewidth',1);
grid,xlabel ('Time'),ylabel ('northings');
datetick;

subplot(313);
hold on;
plot(LAIBt,LAIBz,'b');
%plot(LAIBt,polyval(pLAIBZ,LAIBt),'linewidth',1);
grid,xlabel ('Time'),ylabel ('vertical');
datetick;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fix outliers in the eastings, northing and vertical components if needed

% eastings

LAIBxcutoff=-1.0; % outliers seem to be below this level
ibad=find(LAIBx<LAIBxcutoff); % these are the points we'll change
LAIBx(ibad)=NaN; % replace with NaN

% northings 

LAIBycutoff=0.0; % outliers seem to be below this level
ibad=find(LAIBy<LAIBycutoff); % these are the points we'll change
LAIBy(ibad)=NaN; % replace with NaN

% vertical 

LAIBzcutoff=-1.0; % outliers seem to be below this level
ibad=find(LAIBz<LAIBzcutoff); % these are the points we'll change
LAIBz(ibad)=NaN; % replace with NaN

% plot again after outliers are taken out 

figure(4);
hold on;
set(gca)

subplot(311);
title('LAIB Station NaN for outliers');
hold on;
plot(LAIBt,LAIBx,'g');
%plot(LAIBt,polyval(pLAIBE,LAIBt),'linewidth',1);
grid,xlabel ('Time'),ylabel ('eastings');
datetick;

subplot(312);
hold on;
plot(LAIBt,LAIBy,'r')
%plot(LAIBt,polyval(pLAIBN,LAIBt),'linewidth',1);
grid,xlabel ('Time'),ylabel ('northings');
datetick;

subplot(313);
hold on;
plot(LAIBt,LAIBz,'b');
%plot(LAIBt,polyval(pLAIBZ,LAIBt),'linewidth',1);
grid,xlabel ('Time'),ylabel ('vertical');
datetick;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Fill in values that are NaN with data that smoothes section

[numel(find(isnan(x))),numel(find(isnan(y))),numel(find(isnan(z)))]

LAIBxs=movmedian(LAIBx,50,'omitnan');
LAIBys=movmedian(LAIBy,100,'omitnan');
LAIBzs=movmedian(LAIBz,50,'omitnan');

% plot again with smoothed sections 

figure(5);
hold on;
set(gca)

subplot(311);
title('LAIB Station smoothed');
hold on;
plot(LAIBt,LAIBxs,'g');
%plot(LAIBt,polyval(pLAIBE,LAIBt),'linewidth',1);
grid,xlabel ('Time'),ylabel ('eastings');
datetick;

subplot(312);
hold on;
plot(LAIBt,LAIBys,'r')
%plot(LAIBt,polyval(pLAIBN,LAIBt),'linewidth',1);
grid,xlabel ('Time'),ylabel ('northings');
datetick;

subplot(313);
hold on;
plot(LAIBt,LAIBzs,'b');
%plot(LAIBt,polyval(pLAIBZ,LAIBt),'linewidth',1);
grid,xlabel ('Time'),ylabel ('vertical');
datetick;

% I followed the order from the demo but am 
% thinking you could change it up (ie outliers first?)
% to get a cleaner data set but it may not be 'accurate'



% I tried to load the function in on its own 
% and could not get it to work /:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fill the gaps in the data with NaNs, and make a continuous time vector
%  - note how the time and data vectors get longer
%  - we'll save a copy of the old version, for our reference

function[t2,d2]=filltimegap(t,d,dt)
%
%  [t2,d2]=filltimegap(t,d,dt)
%     t = time column vector, with gaps
%     d = data column vector, with gaps
%     dt = spacing of time vector
%     t2 = time vector, without gaps
%     d2 = data vector, gaps filled with NaNs
%

%
% are there any gaps?
%    
  difft=diff(t)/dt;
  i=find(difft>1.5); % these are the starts of the gaps

  if numel(i)>0
  
    ngap=round(difft(i));  % this is how many points are missing
    i2=[1;i+1;numel(t)+1]; % these are the starts of the data sections
%
% divide the time and data vectors by the gaps  
%
    bits=cell(2,numel(ngap)*2+1);
    for j=1:numel(ngap)+1
      bits{1,(j*2)-1}=t(i2(j):i2(j+1)-1);  
      bits{2,(j*2)-1}=d(i2(j):i2(j+1)-1);  
    end
%
% fill the gaps with nans
%
    for j=1:numel(ngap)
      bits{1,j*2}=(t(i(j))+dt:dt:t(i(j))+dt*(ngap(j)-1))';
      bits{2,j*2}=ones(ngap(j)-1,1)*NaN;
    end
%
% reconstitute the time series
%
    t2=cell2mat(bits(1,:)');
    d2=cell2mat(bits(2,:)');
  
  else
    t2=t;
    d2=d;
  end

  i=find(difft<0.5);
  if numel(i)>0
    't has repeats';
  end
end



