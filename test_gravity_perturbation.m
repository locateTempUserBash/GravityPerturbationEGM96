 clear all; close all; clc;

 addpath 'data';
 addpath 'test_functions';
 
% This script tests for the potential contribution of the geoid
% in terms of spherical harmonics up to degree and order 10 using
% EGM96
 
% EGM-96 constants used here
        re         = 6378.137e3;         % m
        flat       = 1.0/298.257223563;
        omegaearth = 7.292115E-5;     % rad/s
        mu         = 398600.4418;      % km3/s2
        mum        = 3.986004418e14;   % m3/s2

% Spherical harmonic order and degree 10
lMax = 10;
mMax = 10;
lIndMax = lMax-1;

% Initialisation
r0 = re;
Clm = zeros(lIndMax,lIndMax);
Slm = zeros(lIndMax,lIndMax);
P = zeros(lIndMax,lIndMax);
for l = 0:lMax
lInd = l+1;
    P(lInd,1:lInd) = abs(legendre(l,sind(0)));
end

% Reading EGM96 Coefficients
[Clm Slm] = EGM96(lMax, mMax);

% A vector R is created for testing the gravity potential
xVec = linspace(-1.5*re,1.5*re,100);
yVec = linspace(-1.5*re,1.5*re,100);
z = 0;
[XVEC,YVEC] = meshgrid(xVec,yVec);
result = zeros(size(XVEC,1),size(XVEC,2));
iMax = size(XVEC,1);
for i=1:size(XVEC,1)
    for j=1:size(XVEC,2)
        
        % R vector is expressed in geographic lat/long
        R = sqrt(XVEC(i,j)^2+YVEC(i,j)^2+z^2);
        lat = asind(z/R);
        lon = (180/pi)*atan2(YVEC(i,j)/(R*cosd(lat)),XVEC(i,j)/(R*cosd(lat)));
        
        % Gravity potential due to mass distribution
        U_SH = U_spherical_harmonics(lat,lon,R,mu,r0, Clm, Slm);
        
        % The difference btw. SH representation and point mass assumption
        result(i,j) = abs(U_SH-(mu/R^2)/U_SH);
    end
end

theta = linspace(0,2*pi,100);
circleX = r0*cos(theta);
circleY = r0*sin(theta);

figure
contourf(XVEC,YVEC,result,linspace(0,1,10))
caxis([0 1])
title('Spherical Harmonic Contribution')
xlabel('X (m)');
ylabel('Y (m)');
set(gca(),'fontsize',12)
h = colorbar;
ylabel(h,'Relative Difference (%)')
hold on
plot(circleX,circleY,'--w','LineWidth',1.5)

