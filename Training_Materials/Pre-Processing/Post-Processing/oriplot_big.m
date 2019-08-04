%% Script to plot pole figures using MTEX toolbox, from PRISMS-CPFE output
clear ; 
close all ;

%% Specify crystal symmetry. In this case it is Magnesium
cs = crystalSymmetry('6/mmm', [3.21 3.21 5.213], 'X||a*', 'Y||b', 'Z||c*');

%% Specify sample symmetry
ss = specimenSymmetry('triclinic');

%% Plotting convention
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','outOfPlane');

%% Filename - Output from PRISMS CPFE
fname = 'orientationsOutput' ;

%% Read in data
dats = dlmread(fname,' ',0,3) ;

%% Separate Rodrigues vector components
dats = dats(:,1:3) ;

%% Compute vector magnitude
mags = sqrt(dats(:,1).*dats(:,1) +  dats(:,2).*dats(:,2) + dats(:,3).*dats(:,3) ) ;

%% Compute angle
angs = 2*atan(mags) ;

%% Compute axis(normalized)
vecs = [ dats(:,1)./mags(:,1) , dats(:,2)./mags(:,1) , dats(:,3)./mags(:,1) ] ;

%% Define object array of type orientation
v(:,1) = vector3d(vecs(:,1),vecs(:,2),vecs(:,3))   ; 
ori(:,1) = orientation('axis',v(:,1),'angle',angs(:,1),cs,ss) ;


%% Calculate ODF

%psi = deLaValeePoussinKernel('halfwidth',2*degree) ;
%odf = calcODF(ori,'kernel',psi) ;
odf = calcODF(ori,'halfwidth',8*degree) ;

%% Calculate pole figure
pf1 = calcPoleFigure(odf,Miller({0,0,0,1},ori.CS),'resolution',2*degree,'complete') ;
pf2 = calcPoleFigure(odf,Miller({-1,0,1,0},ori.CS),'resolution',2*degree,'complete') ;
pf3 = calcPoleFigure(odf,Miller({-2,1,1,0},ori.CS),'resolution',2*degree,'complete') ;

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