function [Clm, Slm] = EGM96 (lmax, mmax)

addpath 'data'
% Spherical Harmonic Gravity Field Coefficients for EGM96

FID = fopen('egm96.txt');
coefficient_line = fgetl(FID);
C=[];
S=[];
l=0;
m=0;
while 1
    F = sscanf(coefficient_line,'%f');
    
    if (size(F)==0), break, end
    
    if (l==lmax && m==mmax), break, end
    
    l = F(1);
    m = F(2);
    lInd = l + 1;
    mInd = m + 1;
    Clm(lInd,mInd) = F(3);
    Slm(lInd,mInd) = F(4);
    sigmaC(lInd,mInd) = F(5);
    sigmaS(lInd,mInd) = F(6);
    
    coefficient_line = fgetl(FID);

end
fclose(FID);
