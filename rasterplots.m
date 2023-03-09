clear

filename='ridgecrest_after.tif'; % data from Github (after tif) 

[array,metadata]=geotiffread(filename); % read in data and metadata
x=metadata.XWorldLimits ; % access data in metadata 
y=metadata.YWorldLimits ;
z=flipud(array); % flip up and down for map style of data
z(find(z==-9999))=NaN; % 'clean' data

figure (1),clf,hold on % view initial data after inserted NaNs
imagesc(x/1e3,y/1e3,z);
axis xy ;
colorbar;

dz_x=diff(z);
dz_y=diff(z,1,2); % difference of z in the y direction
dz_2=-sqrt(dz_x(:,2:end).^2+dz_y(2:end,:).^2); % differences between the two (dz_x and dz_y)
% ^
% no positive differences?

figure (2),clf, hold on
title('DZ_X') % titled for no confusion 
imagesc(x/1e3,y/1e3,dz_x);
axis xy; axis equal; grid; colorbar;
caxis([-0.75,0.75]); % more defined axis
colormap(gray); % using gray to give good shadows 

figure (3),clf,hold on
title('DZ_Y');
imagesc(x/1e3,y/1e3,dz_y);
axis xy; axis equal; grid; colorbar;
caxis([-0.75,0.75]);
colormap(gray);

figure(4),clf,hold on
title('DZ_2');
imagesc(x/1e3,y/1e3,dz_2);
axis xy; axis equal; grid; colorbar;
caxis([-0.5,0.5]);
colormap(cpolar);

dz=dz_2(find(dz_2<=0.4)); % find all differences <=0.4 
% thinking the large differences could have to do with tectonic movement?

geotiffwrite("ridgecrest_diff",x,y,dz_2); %I know its geotiffwrite but format?




