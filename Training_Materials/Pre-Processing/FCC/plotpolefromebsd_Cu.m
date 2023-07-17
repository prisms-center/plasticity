clear ;
close all ;

% Read from file containing distinct grain IDs and corresponding orientations
fname1 = 'copperdata1.txt' ;
grdat = dlmread(fname1,'',1,0) ;
dats = grdat(:,1:3);

%% Specify crystal symmetry. In this case it is Magnesium
cs = crystalSymmetry('cubic');


%% Specify sample symmetry
ss = specimenSymmetry('triclinic');


%% Define object array of type orientation
ori(:,1) = orientation('Euler',dats(:,1:3),cs,ss) ;

%% Plotting convention
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','outOfPlane');

%% Calculate ODF

psi = deLaValeePoussinKernel('halfwidth',5*degree) ;
odf = calcODF(ori,'kernel',psi) ;
%odf = calcODF(ori,'halfwidth',2*degree) ;

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
