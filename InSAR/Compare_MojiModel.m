clear

%
% Processing the files from S1 InSAR in the SAmerican Andes
%

filename='S1-GUNW-A-R-048-tops-20181128_20181104-232833-28225N_26167N-PP-1c98-v2_0_4.nc';

%ncdisp(filename); % data looks good, all there (downloaded correctly - YAY!)

% list the components included in the file metadata

A=ncread(filename,'/science/grids/data/amplitude')';
y=ncread(filename,'/science/grids/data/latitude');
x=ncread(filename,'/science/grids/data/longitude');
phase=ncread(filename,'/science/grids/data/unwrappedPhase')';
coh=ncread(filename,'/science/grids/data/coherence')';
concomp=ncread(filename,'/science/grids/data/connectedComponents')';
wavelength=ncread(filename,'/science/radarMetaData/wavelength');

%
% convert phase to line of sight displacement
%

LOSdisp=phase*wavelength/4/pi;

i_incoherent=find(concomp==0); % still ocean section defined as '0' values
LOSdisp(i_incoherent)= NaN;

i_incoherent=find(coh<0.4); % again regions within the ocean
LOSdisp(i_incoherent)=NaN;

figure(1),clf
imagesc(x,y,LOSdisp),axis xy,colorbar,title('LOS displacement'),
Colorscale=jet;
Colorscale(1,:)=[1 1 1]; % white [1 1 1] or black [0 0 0]
% colormap(jet)
colormap(Colorscale)
caxis([-0.2,0.2])


%
% Make MOGI model of displacement due to 
% an inflating sphere embeded within an
% elastic halfspace, and calculate the modeled
% line-of-sight displacement for comparison 
% with unwrapped interferograms
%

%
% Note, sphere is radially symmetric, so surface displacement
% depends only on radial distance from center of sphere
%

%
% 1D demo
%

  R=0:100:50e3; % radial points for sampling the predicted displacement field (m)
  F=-1e2; % depth to center of displacement sphere in (m)
  V=1e3; % sphere volumetric change, in m^3 (assumes a point source)
  nu=0.25; % Poisson's ratio
  %A=4e6; % radius of sphere (m)
  %P=1e4; % hydrostatic pressure change in sphere (Pa)
  %E=40e9; % Youngs modulus (Pa)

  [ur,uz,dt,er,et] = mogi(R,F,V,nu); % assumes a point source with volume change
  %[ur,uz,dt,er,et] = mogi(R,F,A,P,E,nu); % gives radius and pressure change

%
% Figure to visualize: 
% radial displacement 
% tagential displacement 
% ground tilt
% radial strain
% tangential strain
%

%   figure(2),clf
%   subplot(511),plot(R/1e3,ur),xlabel('radial distance from source (km)'),ylabel('m'),title('radial displacement')
%   subplot(512),plot(R/1e3,uz),xlabel('radial distance from source (km)'),ylabel('m'),title('tangential displacement')
%   subplot(513),plot(R/1e3,dt),xlabel('radial distance from source (km)'),ylabel('rad'),title('ground tilt')
%   subplot(514),plot(R/1e3,er),xlabel('radial distance from source (km)'),ylabel('m/m'),title('radial strain')
%   subplot(515),plot(R/1e3,et),xlabel('radial distance from source (km)'),ylabel('m/m'),title('tangential strain')

%
% 2D demo:
% Define the grid of points where we'll sample the displacement field
% (units of km from center of dislocation at 0,0)
%
  xmodel=[-50:.1:50];
  ymodel=[-50:.1:50];
  [xmesh,ymesh]=meshgrid(xmodel,ymodel);

  rmesh=sqrt(xmesh.^2+ymesh.^2);
  azimesh=90-atan2d(ymesh,xmesh);

% check that the grids have the right values, and that they're not backward
%   figure(3),clf,
%   subplot(221),imagesc(xmodel,ymodel,xmesh),axis xy,colorbar,title('x value of grid for calculating displacements')
%   subplot(222),imagesc(xmodel,ymodel,ymesh),axis xy,colorbar,title('y value of grid for calculating displacements')
%   subplot(223),imagesc(xmodel,ymodel,rmesh),axis xy,colorbar,title('r value of grid for calculating displacements')
%   subplot(224),imagesc(xmodel,ymodel,azimesh),axis xy,colorbar,title('azimuth value of grid for calculating displacements')

%
% compute displacement at radius values from across the grid
%

  R=rmesh(:)*1e2; % radial points for sampling the predicted displacement field (m)
  F=-1e6; % depth to center of sphere (m)
  nu=0.50; % Poisson's ratio
  V=2e4; % volume of material injected at this point (m^3). Note: 1e6 = 100 meter cube
  A=4e20; % radius of sphere (m)
  P=1e20; % hydrostatic pressure change in sphere (Pa)
  E=10e9; % Young's modulus (Pa)

  [uR,uZ,dT,eR,eT] = mogi(R,F,A,P,E,nu); % gives radius and pressure change
  %[uR,uT] = mogi(R,F,V,nu); % we just care about displacement

%
% convert radial and tangential dislocations into E and N
%
  uRmesh=reshape(uR,size(xmesh));
  eTmesh=reshape(eT,size(xmesh));

  uE=uRmesh.*sind(azimesh);
  uN=uRmesh.*cosd(azimesh);
  uZ=eTmesh; % does tangential mean vertically? not azimuthally?

  figure(4),clf,
  subplot(221),imagesc(xmodel,ymodel,uRmesh),axis xy,colorbar,colormap(Colorscale),title('radial displacement (m)')
  subplot(222),imagesc(xmodel,ymodel,eTmesh),axis xy,colorbar,colormap(Colorscale),title('tangential (vertical) displacement (m)')
  subplot(223),imagesc(xmodel,ymodel,uE),axis xy,colorbar,colormap(Colorscale),title('East displacement (m)')
  subplot(224),imagesc(xmodel,ymodel,uN),axis xy,colorbar,colormap(Colorscale),title('North displacement (m)')

%
% Resolve displacements from East/North/Vertical directions into 
% line-of-sight direction.
%  - Note: actual radar geometry varies a bit across the scene,
%    but using a single number is close enough for our purposes.
%

  % define radar geometry
  IncAngle=46;   % degrees from vertical - ranges from 32 - 46, so this is a good average number
  HeadAngle=-30; % degrees from north - ascending flight direction, in deg E of N
  % HeadAngle=-170; % degrees from north - descending flight direction, in deg E of N

  % unit vector components for the line-of-sight direction
  px=sind(IncAngle)*cosd(HeadAngle);
  py=-sind(IncAngle)*sind(HeadAngle);
  pz=-cosd(IncAngle);

  % displacement in LOS direction is dot product of 3D displacement with LOS unit vector
  uLOS=uE*px+uN*py+uZ*pz;

  % Just for fun (and learning/demonstration), convert the predicted displacement back to wrapped phase
  wavelength=0.05546576;
  uphase=mod(uLOS/wavelength*4*pi,2*pi);

  figure(5),clf
  subplot(121),imagesc(xmodel,ymodel,uLOS),axis xy,colorbar,title('modeled line of sight displacement (m)'); % csym - was not working
  subplot(122),imagesc(xmodel,ymodel,uphase),axis xy,colorbar,title('modeled wrapped phase')
  colormap(jet)




