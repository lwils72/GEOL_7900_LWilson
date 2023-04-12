clear

%
% Processing the files from S1 InSAR in the Florida Everglades & Miami 
%

filename='S1-GUNW-A-R-048-tops-20180929_20180917-232808-26565N_24670N-PP-a173-v2_0_4.nc';

% ncdisp(filename); data looks good, all there (downloaded correctly - YAY!)

% list the components included in the file metadata

A=ncread(filename,'/science/grids/data/amplitude');
y=ncread(filename,'/science/grids/data/latitude');
x=ncread(filename,'/science/grids/data/longitude');
phase=ncread(filename,'/science/grids/data/unwrappedPhase');
coh=ncread(filename,'/science/grids/data/coherence');
concomp=ncread(filename,'/science/grids/data/connectedComponents');
wavelength=ncread(filename,'/science/radarMetaData/wavelength');

%
% Initial plot from the raw data including:
% amplitude, phase, coherence and connected components
%

figure(1),clf
  subplot(2,2,1);
  A = imrotate(A,270);
    imagesc(x,y,A);
    axis xy,colorbar;
    title('amplitude');
    caxis([0,1e4]);
    set(gca, 'XDir','reverse')
  subplot(2,2,2);
  phase = imrotate(phase,270);
    imagesc(x,y,phase);
    axis xy,colorbar;
    title('phase');
    set(gca, 'XDir','reverse')
  subplot(2,2,3);
  coh = imrotate(coh,270);
    imagesc(x,y,coh);
    axis xy,colorbar;
    title('coherence');
    set(gca, 'XDir','reverse')
  subplot(2,2,4);
  concomp = imrotate(concomp,270);
    imagesc(x,y,concomp);
    axis xy, colorbar;
    title('connected components');
    set(gca, 'XDir','reverse');

%
% convert phase to line of sight displacement
%

LOSdisp=phase*wavelength/4/pi;

i_incoherent=find(concomp==0); % still ocean section defined as '0' values
LOSdisp(i_incoherent)= NaN;

i_incoherent=find(coh<0.4); % again regions within the ocean
LOSdisp(i_incoherent)=NaN;

figure(2),clf
imagesc(x,y,LOSdisp),axis xy,colorbar,title('LOS displacement'),
Colorscale=jet;
Colorscale(1,:)=[1 1 1]; % white [1 1 1] or black [0 0 0]
% colormap(jet)
colormap(Colorscale)
caxis([-0.2,0.2])
set(gca, 'XDir','reverse')

% not much movement in southern florida, however there looks to be some
% variaitons throughout the everglades (western portion)
