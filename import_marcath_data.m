%import STDS data 'topography.tif' and lava data 'lava_thickness_flat.tif'
tiff1 = Tiff('topography.tif'); tiff2 = Tiff('lavafield.tif');

%read as matrices, crop (1 arcsec resolution)
topo_post= double(read(tiff1)); lava_thickness_unsmoothed = 11*double(read(tiff2));

[sizey,sizex]=size(topo_post);

%smooth edges of lava over radius 25m
L=25;
lava_thickness_smoothed = imgaussfilt(lava_thickness_unsmoothed,[L/30.87, L/(30.87*0.7826)]);

%subtract to give pre-eruptive topo
topo_pre = topo_post - lava_thickness_smoothed;

%smooth topographies with gaussian radius L metres - change this parameter to
%vary topographic roughness
L=25;
topo_post = imgaussfilt(topo_post,[L/30.87, L/(30.87*0.7826)]);
topo_pre = imgaussfilt(topo_pre,[L/30.87, L/(30.87*0.7826)]);

%x and y axes in metres - 
%30.87m and 30.87*0.7826m are the sizes of one arcsec of latitude and longitude 
spacingx=30.87*0.7826; spacingy=30.87;
xdom=linspace(-1,1,sizex)*sizex*spacingx/2;
ydom=linspace(-1,1,sizey)*sizey*spacingy/2;

%and x and y axes in degrees 
xdms=-115.99167-xdom/(spacingx*3600);
ydms=38.4875+ydom/(spacingy*3600);

%plots
figure; contour(xdom,ydom,topo_pre,50); title('Pre-eruptive topography, Marcath')
%figure; mesh(xdom,ydom,topo_post); title('Post-eruptive topography, Marcath')

%matrices of numerical slopes for transient solver
[dTdx_pre,dTdy_pre] = gradient(topo_pre,spacingx,spacingy);
[dTdx_post,dTdy_post] = gradient(topo_post,spacingx,spacingy);



