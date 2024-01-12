clear ;
close all ;

% Number of voxels in x,y and z directions
xnum = 32 ; 
ynum = 32 ; 
znum = 32 ; 

% Number of header lines
hlines = 20 ;

% Read from file containing distinct grain IDs and corresponding orientations
fname1 = 'orientations.txt' ;
grdat = dlmread(fname1,'',1,0) ;

% Read from file containing all grain IDs
fname1 = 'grainID.txt' ;
micdat = dlmread(fname1,'',[hlines 0 hlines+ynum*znum-1 xnum-1 ]) ;

% Data set containing orientations for all voxels - Rodrigues vector
% convention
sztot = size(micdat,1)*size(micdat,2) ;
micdat = reshape(micdat,[sztot,1]) ;
dats = grdat(micdat(:,1),2:4) ;

%% Specify crystal symmetry. In this case it is Copper
cs = crystalSymmetry('cubic');

%% Specify sample symmetry
ss = specimenSymmetry('triclinic');


%% Compute vector magnitude
mags = sqrt(dats(:,1).*dats(:,1) +  dats(:,2).*dats(:,2) + dats(:,3).*dats(:,3) ) ;

%% Compute angle
angs = 2*atan(mags) ;

%% Compute axis(normalized)
vecs = [ dats(:,1)./mags(:,1) , dats(:,2)./mags(:,1) , dats(:,3)./mags(:,1) ] ;

%% Define object array of type orientation
v(:,1) = vector3d(vecs(:,1),vecs(:,2),vecs(:,3))   ; 
ori(:,1) = orientation('axis',v(:,1),'angle',angs(:,1),cs,ss) ;

%% Plotting convention
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','outOfPlane');

%% Calculate ODF

psi = deLaValeePoussinKernel('halfwidth',5*degree) ;
odf = calcODF(ori,'kernel',psi) ;
%odf = calcODF(ori,'halfwidth',8*degree) ;

%% Calculate pole figure
pf1 = calcPoleFigure(odf,Miller({0,0,1},ori.CS),'resolution',5*degree,'complete') ;
pf2 = calcPoleFigure(odf,Miller({0,1,1},ori.CS),'resolution',5*degree,'complete') ;
pf3 = calcPoleFigure(odf,Miller({1,1,1},ori.CS),'resolution',5*degree,'complete') ;
%%pf4 = calcPoleFigure(odf,Miller({0,0,1},{0,1,1},{1,1,1},ori.CS),'resolution',5*degree,'complete') ;
%% Plot pole figures
figure ;
plot(pf1,'smooth','colorrange','equal') ;
colorbar ;
figure ;
plot(pf2,'smooth','colorrange','equal') ;
colorbar  ;
figure ;
plot(pf3,'smooth','colorrange','equal') ;
colorbar ;
%%figure ;
%%plot(pf4,'smooth','colorrange','equal') ;
%%colorbar ;
