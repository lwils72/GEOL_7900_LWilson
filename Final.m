clear
tic

%
% Presentation : https://docs.google.com/presentation/d/15VSt46SLCdv89lMab9b21E-GOPVyjO7VBI_dRrRAs7o/edit?usp=sharing
%

%
% Six coseismic files for Panhandle FL sinkholes based on 6 scenes:
% Ascending passes on:
%  - 2022-03-22
%  - 2022-04-15
%  - 2022-04-27
%  - 2022-05-09
%

%
% Load the individual interferograms
%
 
  %
  % There are enough of these that it's simpler to handle as a loop.  So put them all in cell arrays
  %
    filename{1}='Final/S1-GUNW-A-R-019-tops-20220427_20220322-234621-00087W_00030N-PP-7231-v2_0_5.nc';
    filename{2}='Final/S1-GUNW-A-R-019-tops-20220427_20220403-234621-00087W_00030N-PP-d8dd-v2_0_5.nc';
    filename{3}='Final/S1-GUNW-A-R-019-tops-20220427_20220415-234621-00087W_00030N-PP-b347-v2_0_5.nc';
    filename{4}='Final/S1-GUNW-A-R-019-tops-20220509_20220403-234622-00087W_00030N-PP-60bd-v2_0_5.nc';
    filename{5}='Final/S1-GUNW-A-R-019-tops-20220509_20220415-234622-00087W_00030N-PP-aa15-v2_0_5.nc';
    filename{6}='Final/S1-GUNW-A-R-019-tops-20220509_20220427-234622-00087W_00030N-PP-575f-v2_0_5.nc';

    T0{1}='2022-04-27';
    T0{2}='2022-04-27';
    T0{3}='2022-04-27';
    T0{4}='2022-05-09';
    T0{5}='2022-05-09';
    T0{6}='2022-05-04';

    T1{1}='2022-03-22';
    T1{2}='2022-04-03';
    T1{3}='2022-04-15';
    T1{4}='2022-04-03';
    T1{5}='2022-04-15';
    T1{6}='2022-04-27';

  %
  % Load them all into a different cell arrays
  %

    for k=1:6
      full_x{k}=ncread(filename{k},'/science/grids/data/longitude');
      full_y{k}=ncread(filename{k},'/science/grids/data/latitude');
      full_u{k}=ncread(filename{k},'/science/grids/data/unwrappedPhase')'; % unwrapped phase (radians)
      full_c{k}=ncread(filename{k},'/science/grids/data/coherence')';
      full_m{k}=ncread(filename{k},'/science/grids/data/connectedComponents')';
      full_a{k}=ncread(filename{k},'/science/grids/data/amplitude')';
    end
    L=ncread(filename{k},'/science/radarMetaData/wavelength'); % wavelength (m)

%
% Plot the components
%
  figure(1),clf
  for k=1:6
    subplot(6,4,(k-1)*4+1),imagesc(full_x{k},full_y{k},full_u{k}),axis xy,colorbar,title([T0{k},' - ',T1{k},' unwrappedPhase'])
    subplot(6,4,(k-1)*4+2),imagesc(full_x{k},full_y{k},full_c{k}),axis xy,colorbar,title([T0{k},' - ',T1{k},' coherence'])
    subplot(6,4,(k-1)*4+3),imagesc(full_x{k},full_y{k},full_m{k}),axis xy,colorbar,title([T0{k},' - ',T1{k},' connectedComponents'])
    subplot(6,4,(k-1)*4+4),imagesc(full_x{k},full_y{k},full_a{k}),axis xy,colorbar,title([T0{k},' - ',T1{k},' amplitude']),caxis([0,1e4])
  end

%
% Trim them to all be the same size and area
%

  %
  % Determine smaller area that have the same limits (for math)
  %

    figure(21),clf % determine the right coordinates
    imagesc(full_x{1},full_y{1},full_c{1});
    axis xy;
    colorbar;
    title([T0{1},' - ',T1{1},' coherence']);

    figure(23),clf % determine the right coordinates
    imagesc(full_x{3},full_y{3},full_c{3});
    axis xy;
    colorbar;
    title([T0{3},' - ',T1{3},' coherence']);

    x0=-86.6;
    x1=-85.0;
    y0=30.68;
    y1=31.7;

    dx=1/60/20; % how many degrees apart are the points (20 arcseconds)
    dy=dx;

  %
  % make the vectors that will represent the new trimmed edge coordinates
  %
    x=x0:dx:x1;
    y=y0:dy:y1;

  %
  % Loop over each interferogram, trim the data, and store as layers
  % of a 3D array (ny x nx x 4)
  %

    for k=1:6
      ix=find(full_x{k}>=x0-dx/2 & full_x{k}<=x1+dx/2);
      iy=find(full_y{k}>=y0-dy/2 & full_y{k}<=y1+dy/2);

      U(:,:,k)=flipud(full_u{k}(iy,ix));
      C(:,:,k)=flipud(full_c{k}(iy,ix));
      M(:,:,k)=flipud(full_m{k}(iy,ix));
      A(:,:,k)=flipud(full_a{k}(iy,ix));
    end

%
% How different are each of the scenes?
%

  % Convert to LOS displacement
    LOS=U*L/4/pi;

figure(3),clf
subplot(3,2,1),imagesc(x,y,LOS(:,:,1)),axis xy,colorbar,title('LOS displacement'),title([T0{1},' - ',T1{1};
subplot(3,2,2),imagesc(x,y,LOS(:,:,2)),axis xy,colorbar,title('LOS displacement'),title([T0{2},' - ',T1{2};
subplot(3,2,3),imagesc(x,y,LOS(:,:,3)),axis xy,colorbar,title('LOS displacement'),title([T0{3},' - ',T1{3};
subplot(3,2,4),imagesc(x,y,LOS(:,:,4)),axis xy,colorbar,title('LOS displacement'),title([T0{4},' - ',T1{4};
subplot(3,2,5),imagesc(x,y,LOS(:,:,5)),axis xy,colorbar,title('LOS displacement'),title([T0{5},' - ',T1{5};
subplot(3,2,6),imagesc(x,y,LOS(:,:,6)),axis xy,colorbar,title('LOS displacement'),title([T0{6},' - ',T1{6};
%Colorscale(1,:)=[1 1 1]; % white [1 1 1] or black [0 0 0]
colormap(jet);
caxis([-0.2,0.2]); 

  %
  % Plot each LOS on the diagonals, and the difference between each pair on the off-diagnonals
  % (just to see how they differ, not for analysis)
  %
    figure(4),clf
    for k=1:3
      ax((k-1)*3+k)=subplot(3,6,(k-1)*3+k);
      imagesc(x,y,LOS(:,:,k)),axis xy,colorbar, %csym(1)
      title([T0{k},' - ',T1{k},' LOS (m)'])
      set(gca,'dataaspectratio',[1/cosd(35.7),1,1])
      for j=k+1:3
        ax((k-1)*3+j)=subplot(3,6,(k-1)*3+j);
        imagesc(x,y,LOS(:,:,k)-LOS(:,:,j)),axis xy,colorbar, %csym(0.15)
        title('difference')
        set(gca,'dataaspectratio',[1/cosd(35.7),1,1])
      end
    end
    colormap(cpolar);
    linkaxes(ax,'xy');

%
% Wanting to do source separation (atmosphere vs motion)
% 

  %
  % reshape each 2D image into a 1D vector of observations
  %

    nx=numel(x);
    ny=numel(y);

    LOS1d=reshape(LOS,[nx*ny,6,1]);


for k=1:6;
  i_incoherent=find(isnan(LOS1d(k))); % still ocean section defined as '0' values
  LOS(i_incoherent)= 0;

  i_incoherent=find(LOS1d(k)==NaN); % still ocean section defined as '0' values
  LOS(i_incoherent)=0;


%   i_incoherent=find(LOS1d{k}<=0.02); % still ocean section defined as '0' values
%   LOS(i_incoherent)= 0.2;

end

  %
  % Check that this rearranged everything the way we intended: can we get one back? yes!
  %
    figure(99),clf,imagesc(reshape(LOS1d(:,1),ny,nx)),axis xy,colorbar

  %
  % use Principal Component Analysis to identify the different
  % "signal" and "noise" contributions
  %
    [PCAcomp1D,PCAweight,~,~,explained,~]=pca(LOS1d','centered',false);

  %
  % Plot the PCA components ("sources")
  %

    figure(5),clf;
    for k=1:6
      PCAcomp1D(:,:,k)=reshape(PCAcomp1D(:,k),ny,nx);
      subplot(2,2,k),imagesc(x,y,PCAcomp(:,:,k)),axis xy,colorbar,
      csym(quantile(abs(PCAcomp(:,:,k)),0.97,'all')) % automatically choose useful color limits
      title(['this component explains ',num2str(explained(k),'%0.1f'),'% of the total data variance'])
    end
    colormap(cpolar);

  %
  % Plot the PCA weights (how much is each component present in each interferogram?)
  %
    PCAweight % this prints the numbers to the screen

    figure(6),clf
    imagesc(PCAweight),colorbar,csym
    colormap(cpolar)
    xticks(1:numel(filename)),xlabel('PCA component number')
    yticks(1:numel(filename)),ylabel('interferogram number')

    stop 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% three coseismic files for Panhandle FL sinkholes based on 3 scenes:
% Ascending passes on:
%  - 2022-03-22
%  - 2022-04-15
%  - 2022-04-27
%  - 2022-05-09
%

%
% Load the individual interferograms
%
 
  %
  % There are enough of these that it's simpler to handle as a loop.  So put them all in cell arrays
  %
    filename{1}='Final/S1-GUNW-A-R-121-tops-20220504_20220329-233759-00085W_00030N-PP-418f-v2_0_5.nc';
    filename{2}='Final/S1-GUNW-A-R-121-tops-20220504_20220410-233759-00085W_00030N-PP-2375-v2_0_5.nc';
    filename{3}='Final/S1-GUNW-A-R-121-tops-20220504_20220422-233759-00085W_00030N-PP-9aa6-v2_0_5.nc';
    
    T0{1}='2022-05-04';
    T0{2}='2022-05-04';
    T0{3}='2022-05-04';
   
    T1{1}='2022-03-29';
    T1{2}='2022-04-10';
    T1{3}='2022-04-22';
  
  %
  % Load them all into a different cell arrays
  %

    for k=1:3
      full_x{k}=ncread(filename{k},'/science/grids/data/longitude');
      full_y{k}=ncread(filename{k},'/science/grids/data/latitude');
      full_u{k}=ncread(filename{k},'/science/grids/data/unwrappedPhase')'; % unwrapped phase (radians)
      full_c{k}=ncread(filename{k},'/science/grids/data/coherence')';
      full_m{k}=ncread(filename{k},'/science/grids/data/connectedComponents')';
      full_a{k}=ncread(filename{k},'/science/grids/data/amplitude')';
    end
    L=ncread(filename{k},'/science/radarMetaData/wavelength'); % wavelength (m)

%
% Plot the components
%
  figure(1),clf
  for k=1:3
    subplot(3,4,(k-1)*4+1),imagesc(full_x{k},full_y{k},full_u{k}),axis xy,colorbar,title([T0{k},' - ',T1{k},' unwrappedPhase'])
    subplot(3,4,(k-1)*4+2),imagesc(full_x{k},full_y{k},full_c{k}),axis xy,colorbar,title([T0{k},' - ',T1{k},' coherence'])
    subplot(3,4,(k-1)*4+3),imagesc(full_x{k},full_y{k},full_m{k}),axis xy,colorbar,title([T0{k},' - ',T1{k},' connectedComponents'])
    subplot(3,4,(k-1)*4+4),imagesc(full_x{k},full_y{k},full_a{k}),axis xy,colorbar,title([T0{k},' - ',T1{k},' amplitude']),caxis([0,1e4])
  end

   figure(21),clf % determine the right coordinates
    imagesc(full_x{1},full_y{1},full_c{1});
    axis xy;
    colorbar;
    title([T0{1},' - ',T1{1},' coherence']);

    figure(23),clf % determine the right coordinates
    imagesc(full_x{3},full_y{3},full_c{3});
    axis xy;
    colorbar;
    title([T0{3},' - ',T1{3},' coherence']);


    x0=-83.65;
    x1=-82.5;
    y0=30.8;
    y1=30.2;

    dx=1/60/20; % how many degrees apart are the points (20 arcseconds)
    dy=dx;

  %
  % make the vectors that will represent the new trimmed edge coordinates
  %
    x=x0:dx:x1;
    y=y0:dy:y1;

  %
  % Loop over each interferogram, trim the data, and store as layers
  % of a 3D array (ny x nx x 4)
  %

    for k=1:3
      ix=find(full_x{k}>=x0 & full_x{k}<=x1);
      iy=find(full_y{k}>=y0 & full_y{k}<=y1);

      U(:,:,k)=flipud(full_u{k}(iy,ix));
      C(:,:,k)=flipud(full_c{k}(iy,ix));
      M(:,:,k)=flipud(full_m{k}(iy,ix));
      A(:,:,k)=flipud(full_a{k}(iy,ix));
    end

%
% How different are each of the scenes?
%

  % Convert to LOS displacement
    LOS=U*L/4/pi;

figure(3),clf
subplot(3,2,1),imagesc(x,y,LOS(:,:,1)),axis xy,colorbar,title('LOS displacement'),title([T0{1},' - ',T1{1};
subplot(3,2,2),imagesc(x,y,LOS(:,:,2)),axis xy,colorbar,title('LOS displacement'),title([T0{2},' - ',T1{2};
subplot(3,2,3),imagesc(x,y,LOS(:,:,3)),axis xy,colorbar,title('LOS displacement'),title([T0{3},' - ',T1{3};
subplot(3,2,4),imagesc(x,y,LOS(:,:,4)),axis xy,colorbar,title('LOS displacement'),title([T0{4},' - ',T1{4};
subplot(3,2,5),imagesc(x,y,LOS(:,:,5)),axis xy,colorbar,title('LOS displacement'),title([T0{5},' - ',T1{5};
subplot(3,2,6),imagesc(x,y,LOS(:,:,6)),axis xy,colorbar,title('LOS displacement'),title([T0{6},' - ',T1{6};
%Colorscale(1,:)=[1 1 1]; % white [1 1 1] or black [0 0 0]
colormap(jet);
caxis([-0.2,0.2]); 


  %
  % Plot each LOS on the diagonals, and the difference between each pair on the off-diagnonals
  % (just to see how they differ, not for analysis)
  %
    figure(4),clf
    for k=1:3
      ax((k-1)*3+k)=subplot(3,6,(k-1)*3+k);
      imagesc(x,y,LOS(:,:,k)),axis xy,colorbar, %csym(1)
      title([T0{k},' - ',T1{k},' LOS (m)'])
      set(gca,'dataaspectratio',[1/cosd(35.7),1,1])
      for j=k+1:3
        ax((k-1)*3+j)=subplot(3,6,(k-1)*3+j);
        imagesc(x,y,LOS(:,:,k)-LOS(:,:,j)),axis xy,colorbar, %csym(0.15)
        title('difference')
        set(gca,'dataaspectratio',[1/cosd(35.7),1,1])
      end
    end
    colormap(cpolar);
    linkaxes(ax,'xy');


%
% Implementing LiDAR data to look at certain sink hole areas in 2022..
% only found sections in my area for 2018 and data is not consistent and
% not usable
%

    %
    % LiDAR raster differencing and iterative closest point differencing
    %

filename='output_USGS30m';
[array,metadata]=geotiffread(filename);
x=metadata.XWorldLimits;
y=metadata.YWorldLimits;
z=flipud(array);
z(find(z==-9999))=NaN;

%
% Showing shadows and differencing the before and after rasters
%

dz_x=diff(z);
dz_y=diff(z,1,2);
dz_2=-sqrt(dz_x(:,2:end).^2+dz_y(2:end,:).^2);

figure(1),clf
imagesc(x/1e3,y/1e3,dz_x)
axis xy
axis equal
colorbar
caxis([-1,1])
colormap(gray)

