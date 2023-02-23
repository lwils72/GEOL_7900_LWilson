clear
tic

%
% Mini project on the hydrologic loading in the Honk Kong wetlands
% Link for my presentation : https://docs.google.com/presentation/d/1EVVXRmTVCodxjPbl8RmJxQO4P3DWpIbunHpLs1GGgRQ/edit?usp=sharing
%

%  !curl http://geodesy.unr.edu/gps_timeseries/tenv3/plates/EU/HKWS.EU.tenv3 > HKWS.EU.tenv3.txt;
%  !curl http://geodesy.unr.edu/gps_timeseries/tenv3/plates/EU/HKLT.EU.tenv3 > HKLT.EU.tenv3.txt;
%  !curl http://geodesy.unr.edu/gps_timeseries/tenv3/plates/EU/HKNP.EU.tenv3 > HKNP.EU.tenv3.txt;
%  !curl http://geodesy.unr.edu/gps_timeseries/tenv3/plates/EU/DSMG.EU.tenv3 > DSMG.EU.tenv3.txt;
%  !curl http://geodesy.unr.edu/gps_timeseries/tenv3/plates/EU/HKMW.EU.tenv3 > HKMW.EU.tenv3.txt;
%  !curl http://geodesy.unr.edu/gps_timeseries/tenv3/plates/EU/HKKS.EU.tenv3 > HKKS.EU.tenv3.txt;

  %
  % Load and parse the data
  % The stations from Indochina: CPNM, NIMT, CMUM, HKLM, KMNM, JUNA, PTGG
  % The stations from Hong Kong: HKWS, HKLT, HKOK, DSMG, HKMW
  %
  % Going with Hong Kong GNSS stations for a more localized region with a
  % large amount of stations and to compare to Jiang et al. 2021
  %
    mystations={'HKWS','HKLT','HKNP','DSMG','HKMW','HKKS'};
    nGNSS=numel(mystations); 

    for k=1:nGNSS
      station=mystations{k};
      %system(['!curl http://geodesy.unr.edu/gps_timeseries/tenv3/plates/EU/',station,'.EU.tenv3' > 'MATLAB/',station,'.EU.tenv3']);
      filename=['MATLAB/',station,'.EU.tenv3.txt'];
      fid=fopen(filename);
      C=textscan(fid,'%s %s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f','headerlines',1);
      fclose(fid);

      X{k}=C{9}-C{9}(1);
      Y{k}=C{11}-C{11}(1);
      Z{k}=C{13}-C{13}(1);
     % T{k}=datenum(C{2},'yyMMMdd');
      T{k} = datetime(C{2},'InputFormat','yyMMMdd');

      stationlat(k)=C{21}(1);
      stationlon(k)=mod(C{22}(1)+180,360)-180; % make sure the lon is [-180,180]
      stationelv(k)=C{23}(1);

    end

  %
  % plot of raw data and station locations
  %
  
    shift=-0.05;
    figure(1),clf
    for k=1:nGNSS
      subplot(321),plot(T{k},X{k}-shift*(k-1),'.'),datetick,hold on
      subplot(323),plot(T{k},Y{k}-shift*(k-1),'.'),datetick,hold on
      subplot(325),plot(T{k},Z{k}-shift*(k-1),'.'),datetick,hold on
    end

    c=load('/Users/leorawilson/Desktop/LSU23/MATLAB/DATS/coastfile.xy');
    b=load('/Users/leorawilson/Desktop/LSU23/MATLAB/DATS/politicalboundaryfile.xy');
    subplot(3,2,[2,4,6]),
      plot(c(:,1),c(:,2)),hold on,plot(b(:,1),b(:,2),'color',[1 1 1]*0.75)
      plot(stationlon,stationlat,'k^') % plot stations as triangles 
      text(stationlon,stationlat-0.02,mystations) % since my map is zoomed in, titles are 0.02 away from station
    R=[113,115,21.6,23];
    axis(R)
    subplot(321),legend(mystations,'location','northeast')


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Fill time gaps with NaNs, trim time series, and smooth
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

    t0=datenum([2016,01,01]);
    t1=datenum([2022,12,31]); % time series shared by all stations

    dt=1;
    t=(t0:dt:t1);
    t=t.';

    for k=1:nGNSS;

            [t2,z2]=filltimegap(t,Z{k},dt); 
            [t2,x2]=filltimegap(t,X{k},dt); 
            [t2,y2]=filltimegap(t,Y{k},dt);

      iduring=find(t2>=t0 & t2<=t1); % time between 2014 and 2022

      easting(:,k)=x2(iduring)-x2(iduring(1)); % make the first value in the time series be zero
      northing(:,k)=y2(iduring)-y2(iduring(1));
      vertical(:,k)=z2(iduring)-z2(iduring(1));
       
    end

    dt_smooth=20; % smoothing window 
    sm_easting=movmedian(easting,dt_smooth,2,'omitnan');
    sm_northing=movmedian(northing,dt_smooth,2,'omitnan');
    sm_vertical=movmedian(vertical,dt_smooth,2,'omitnan');

    % verify we got rid of all the NaNs by smoothing...
    [numel(find(isnan(easting))),numel(find(isnan(northing))),numel(find(isnan(vertical)));...
     numel(find(isnan(sm_easting))),numel(find(isnan(sm_northing))),numel(find(isnan(sm_vertical)))]

  %
  % Quick plot of the cleaned and trimmed data of the vertical portions
  % shows a general pattern with all the vertical variations
  % at different elevations of the stations, the pattern changes
  %

    figure(12),clf
        plot(t,vertical(:,1)-shift*(1-1),'.'),datetick,hold on
        plot(t,vertical(:,2)-shift*(2-1),'.'),datetick
        plot(t,vertical(:,3)-shift*(3-1),'.'),datetick
        plot(t,vertical(:,4)-shift*(4-1),'.'),datetick
        plot(t,vertical(:,5)-shift*(5-1),'.'),datetick
        plot(t,vertical(:,6)-shift*(6-1),'.'),datetick
        legend(mystations,'location','northeast'),set(legend,'fontsize',15),


 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Velocities over the set time period
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k=1:nGNSS;
  
      vE(:,k) = polyfit(t,easting(:,k),1); %fits plots for northings
      VE(:,k) = vE(1);
      vN(:,k) = polyfit(t,northing(:,k),1); %fits plots for northings
      VN(:,k) = vN(1);
      vZ(:,k) = polyfit(t,vertical(:,k),1); %fits plots for vertical
      VZ(:,k) = vZ(1);

end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot new time series with best fit, velocities and color scale for
% magnitude
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    figure(2),clf
    for k=1:nGNSS
      subplot(321),plot(t,easting(k)-shift*(k-1),'.'),datetick,hold on,plot(t,polyval(easting(k),t),'linewidth',1);
      subplot(323),plot(t,northing(k)-shift*(k-1),'.'),datetick,hold on,plot(t,polyval(northing(k),t),'linewidth',1);
      subplot(325),plot(t,vertical(k)-shift*(k-1),'.'),datetick,hold on,plot(t,polyval(vertical(k),t),'linewidth',1);
    end

    c=load('/Users/leorawilson/Desktop/LSU23/MATLAB/DATS/coastfile.xy');
    b=load('/Users/leorawilson/Desktop/LSU23/MATLAB/DATS/politicalboundaryfile.xy');
    subplot(3,2,[2,4,6]),hold on
      plot(c(:,1),c(:,2)),plot(b(:,1),b(:,2),'color',[1 1 1]*0.75)
      plot(stationlon,stationlat,'k^') % plot stations as triangles 
      text(stationlon,stationlat-0.02,mystations) % since my map is zoomed in, titles are 0.02 away from station
      quiver(stationlon,stationlat,VE,VN);
      scatter(stationlon,stationlat,50,VZ*1000,'filled');
    R=[113,115,21.6,23];
    axis(R)
    subplot(321),legend(mystations,'location','northeast')
    title('Location of stations with velocity vectors')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PCA and ICA analyses from Luttrell
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %
  % PCA analysis: vertical only
  %
    [PCA_comp2,PCA_weight2]=pca(sm_vertical);
    Ncomponents=size(PCA_comp2,2);

    f3=figure(3);clf,f4.Name='PCA Vertical only';
    for k=1:Ncomponents
      subplot(Ncomponents,2,k*2-1),
      plot(t,PCA_weight2(:,k)),datetick,grid
      title(['PCA component # ',num2str(k)])
    end
    subplot(122),
      imagesc(PCA_weight2'),colorbar
      xticks(1:nGNSS),xticklabels(mystations),xlabel('station')
      yticks(1:Ncomponents),ylabel('principal component number')
      title('PCA weights at each station')

      % After running this multiple times there was only 1 component that
      % was calculated 

  %
  % ICA analysis: vertical only
  %
    Ncomponents=6;

    data=sm_vertical;
    Mdl=rica(data,Ncomponents);
    ICA_comp2=Mdl.TransformWeights;
    ICA_weight2=transform(Mdl,data);
    toc

    f4=figure(4);clf,f7.Name='ICA Vertical only';
    for k=1:Ncomponents
      subplot(Ncomponents,2,k*2-1),
      plot(t,ICA_weight2(:,k)),datetick,grid
      title(['ICA component # ',num2str(k)])
    end
    subplot(122),
      imagesc(ICA_weight2'),colorbar
      xticks(1:nGNSS),xticklabels(mystations),xlabel('station')
      yticks(1:Ncomponents),ylabel('principal component number')
      title('ICA weights at each station')


%
%   Focus just the small annual variations in hydrologic loading in and around
%   Hong Kong, ie plotting only the components at each location together
%   Work on smoothing the data or collected the average for a year?
%

%
% Smoothing for just a year interval ... I tried with various loops but
% they were giving me problems..
%


    t02=datenum([2020,02,01]);
    t12=datenum([2021,02,01]); % time series to look at a single hydroperiod

    dt2=1;
    tt=(t02:dt2:t12);
    tt=tt.';

      iduring2=find(tt>=t02 & tt<=t12); 

ICA1=ICA_weight2(:,1);
ICA2=ICA_weight2(:,2);
ICA3=ICA_weight2(:,3);
ICA4=ICA_weight2(:,4);
ICA5=ICA_weight2(:,5);
ICA6=ICA_weight2(:,6);

      ICAs1=ICA1(iduring2); 
      ICAs2=ICA2(iduring2); 
      ICAs3=ICA3(iduring2); 
      ICAs4=ICA4(iduring2); 
      ICAs5=ICA5(iduring2); 
      ICAs6=ICA6(iduring2); 

    dt_smooth2=50; 
    ICAs1=movmedian(ICAs1,dt_smooth2,2);
    ICAs2=movmedian(ICAs2,dt_smooth2,2);
    ICAs3=movmedian(ICAs3,dt_smooth2,2);
    ICAs4=movmedian(ICAs4,dt_smooth2,2);
    ICAs5=movmedian(ICAs5,dt_smooth2,2);
    ICAs6=movmedian(ICAs6,dt_smooth2,2);


%
% Fit a plot for the ICA_weights, type of 
%
  
      vI1 = polyfit(tt,ICAs1,1); 
      VI1 = vI1(1);
      vI2 = polyfit(tt,ICAs2,1); 
      VI2 = vI2(1);
      vI3 = polyfit(tt,ICAs3,1); 
      VI3 = vI3(1);
      vI4 = polyfit(tt,ICAs4,1);
      VI4 = vI4(1);
      vI5 = polyfit(tt,ICAs5,1); 
      VI5 = vI5(1);
      vI6 = polyfit(tt,ICAs6,1);  
      VI6 = vI6(1);

%
% Plot only the best fit lines of the hydroperiods from the ICA Analysis
% for 2020-2021
%

    figure(5),clf
        subplot(611), hold on, plot(tt,ICAs1), datetick,title('Station HKWS');
        subplot(612), hold on, plot(tt,ICAs2), datetick,title('Station HKLT');
        subplot(613), hold on, plot(tt,ICAs3), datetick,title('Station HKNP');
        subplot(614), hold on, plot(tt,ICAs4), datetick,title('Station DSMG');
        subplot(615), hold on, plot(tt,ICAs5), datetick,title('Station HKMW');
        subplot(616), hold on, plot(tt,ICAs6), datetick,title('Station HKKS');


VI=[VI1,VI2,VI3,VI4,VI5,VI6]


    figure(6),clf
        hold on;
        plot(c(:,1),c(:,2)),plot(b(:,1),b(:,2),'color',[1 1 1]*0.75)
        plot(stationlon,stationlat,'c^') % plot stations  
        text(stationlon,stationlat-0.03,mystations) % since my map is zoomed in, titles are 0.02 away from station
        scatter(stationlon,stationlat,50,VI*1000,'filled');
        quiver(stationlon,stationlat,VE,VN,'b');
    R=[113,115,21.6,23];
    axis(R)
    imagesc(VI'),colorbar
    title('Overall trend of hydrologic period over a year')

%
% Thinking about working out the average amplitude of the hydroperiods from
% ICA analysis 
%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Function 'filltimegap'
% Even with matlab code for this function in my directory
% it still will not work through it without this lower part
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [t2,d2]=filltimegap(t,d,dt)
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







