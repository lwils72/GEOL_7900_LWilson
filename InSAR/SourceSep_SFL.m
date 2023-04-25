clear
tic

%
% four coseismic files for Southern FL sinkholes based on 4 scenes:
% Ascending passes on:
%  - 2017-11-21
%  - 2018-10-23
%  - 2018-11-04
%  - 2018-11-16
%  - 2018-11-28
%

%
% Load the individual interferograms
%
 
  %
  % There are enough of these that it's simpler to handle as a loop.  So put them all in cell arrays
  %
    filename{1}='InSAR/S1-GUNW-A-R-048-tops-20181128_20181116-232808-26565N_24669N-PP-f47a-v2_0_4.nc';
    filename{2}='InSAR/S1-GUNW-A-R-048-tops-20181128_20181104-232808-26565N_24669N-PP-b5e5-v2_0_4.nc';
    filename{3}='InSAR/S1-GUNW-A-R-048-tops-20181128_20181023-232808-26565N_24669N-PP-d826-v2_0_4.nc';
    filename{4}='InSAR/S1-GUNW-A-R-048-tops-20181128_20171121-232808-26565N_24669N-PP-a0f0-v2_0_4.nc';

    T0{1}='2018-11-28';
    T0{2}='2018-11-28';
    T0{3}='2018-11-28';
    T0{4}='2018-11-28';

    T1{1}='2018-11-16';
    T1{2}='2018-11-04';
    T1{3}='2018-10-23';
    T1{4}='2017-11-21';

  %
  % Load them all into a different cell arrays
  %

    for k=1:4
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
  for k=1:4
    subplot(4,4,(k-1)*4+1),imagesc(full_x{k},full_y{k},full_u{k}),axis xy,colorbar,title([T0{k},' - ',T1{k},' unwrappedPhase'])
    subplot(4,4,(k-1)*4+2),imagesc(full_x{k},full_y{k},full_c{k}),axis xy,colorbar,title([T0{k},' - ',T1{k},' coherence'])
    subplot(4,4,(k-1)*4+3),imagesc(full_x{k},full_y{k},full_m{k}),axis xy,colorbar,title([T0{k},' - ',T1{k},' connectedComponents'])
    subplot(4,4,(k-1)*4+4),imagesc(full_x{k},full_y{k},full_a{k}),axis xy,colorbar,title([T0{k},' - ',T1{k},' amplitude']),caxis([0,1e4])
  end


%
% Trim them to all be the same size and area
%

  %
  % pick a subset and make smaller versions that have the same limits (so we can do math lol)
  %
    x0=-82.5;
    x1=-79.5;
    y0=24.5;
    y1=26.5;

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

    for k=1:4
      ix=find(full_x{k}>=x0 & full_x{k} <=x1);
      iy=find(full_y{k}>=y0 & full_y{k} <=y1);

      U(:,:,k)=flipud(full_u{k}(iy,ix));
      C(:,:,k)=flipud(full_c{k}(iy,ix));
      M(:,:,k)=flipud(full_m{k}(iy,ix));
      A(:,:,k)=flipud(full_a{k}(iy,ix));
    end

  %
  % Plot the trimmed components, just to make sure everything still looks ok
  %
    figure(2),clf,
    for k=1:4
      subplot(4,4,(k-1)*4+1),imagesc(x,y,U(:,:,k)),axis xy,colorbar,title([T0{k},' - ',T1{k},' unwrappedPhase'])
      subplot(4,4,(k-1)*4+2),imagesc(x,y,C(:,:,k)),axis xy,colorbar,title([T0{k},' - ',T1{k},' coherence'])
      subplot(4,4,(k-1)*4+3),imagesc(x,y,M(:,:,k)),axis xy,colorbar,title([T0{k},' - ',T1{k},' connectedComponents'])
      subplot(4,4,(k-1)*4+4),imagesc(x,y,A(:,:,k)),axis xy,colorbar,title([T0{k},' - ',T1{k},' amplitude']),caxis([0,1e4])
    end

%
% How different are each of the scenes?
%

  % Convert to LOS displacement
    LOS=U*L/4/pi;

  %
  % Plot each LOS on the diagonals, and the difference between each pair on the off-diagnonals
  % (just to see how they differ, not for analysis)
  %
    figure(3),clf
    for k=1:4
      ax((k-1)*4+k)=subplot(4,4,(k-1)*4+k);
      imagesc(x,y,LOS(:,:,k)),axis xy,colorbar, %csym(1)
      title([T0{k},' - ',T1{k},' LOS (m)'])
      set(gca,'dataaspectratio',[1/cosd(35.7),1,1])
      for j=k+1:4
        ax((k-1)*4+j)=subplot(4,4,(k-1)*4+j);
        imagesc(x,y,LOS(:,:,k)-LOS(:,:,j)),axis xy,colorbar, %csym(0.15)
        title('difference')
        set(gca,'dataaspectratio',[1/cosd(35.7),1,1])
      end
    end
    colormap(cpolar)
    linkaxes(ax,'xy')


%    
% Source separation: pull out atmosphere vs motion?
%
  %
  % reshape each 2D image into a 1D vector of observations
  % Since we have 4 images in a 3D array (ny x nx x 4),
  % the result will be a new 2D array (ny*nx x 4)
  %
    nx=numel(x);
    ny=numel(y);
    n=nx*ny;


    LOS1d=reshape(LOS,size(n,4,1));

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
    figure(4),clf
    for k=1:4
      PCAcomp(:,:,k)=reshape(PCAcomp1D(:,k),ny,nx);
      subplot(2,2,k),imagesc(x,y,PCAcomp(:,:,k)),axis xy,colorbar,
      csym(quantile(abs(PCAcomp(:,:,k)),0.97,'all')) % automatically choose useful color limits
      title(['this component explains ',num2str(explained(k),'%0.1f'),'% of the total data variance'])
    end
    colormap(cpolar)






