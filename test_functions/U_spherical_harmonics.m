function U = U_spherical_harmonics( lat, lon, r, mu, r0, Clm, Slm )

% Calculates the value of the spherical harmonic representation of the
% gravity potential for a point given the latitude, longitude and radius.

% Check the degree of the spherical harmonic expression
lIndMax = size(Clm,1);
lMax = lIndMax-1;
% Prepare matrix of potential contribution and legendre polynomials
potContrib = zeros(lIndMax,lIndMax);
P = zeros(lIndMax,lIndMax);
for l = 0:lMax
lInd = l+1;
    P(lInd,1:lInd) = abs(legendre(l,sind(lat)));
end
% Populate matrix of the contribution of each degree and order
for l = 0:lMax
    for m = 0:l
        lInd = l+1;
        mInd = m+1;
        potContrib(lInd,mInd) = (r0/r)^l*P(lInd,mInd)...
          *(Clm(lInd,mInd)*cosd(m*lon)+Slm(lInd,mInd)*sind(m*lon));
    end
end
% Calculate gravity potential
U = (mu/r)*sum(sum(potContrib,2),1);
 
end
