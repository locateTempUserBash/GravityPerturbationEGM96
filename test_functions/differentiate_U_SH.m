
function GradU = differentiate_U_SH (lmax, mmax)

 mu         = 398600.4418;      % km3/s2
[Clm,Slm]=EGM96(lmax, mmax);

syms x y z r m r0 lon lat latdeg londeg
P = sym('P',[lmax mmax]);
potContrib = sym('potContrib',[lmax mmax]);

a = sin(lat);
b = cos(m*lon);
c = sin(m*lon);

lIndMax = size(Clm,1);
lMax = lIndMax-1;
% Prepare matrix of potential contribution and legendre polynomials
%potContrib = zeros(lIndMax,lIndMax);
%P = zeros(lIndMax,lIndMax);
for l = 0:lMax
lInd = l+1;
    polyActual = abs(legendreP(l,sin(lat)));
    P(lInd,1:lInd) = polyActual;
end
% Populate matrix of the contribution of each degree and order
for l = 0:lMax
    for m = 0:l
        lInd = l+1;
        mInd = m+1;
        potContrib(lInd,mInd) = (r0/r)^l*P(lInd,mInd)...
          *(Clm(lInd,mInd)*cos(m*lon)+Slm(lInd,mInd)*sin(m*lon));
    end
end
% Calculate gravity potential
U = (mu/r)*sum(sum(potContrib,2),1);

subs(U,[lat lon],[latdeg*(pi/180) londeg*(pi/180)])

U = subs(U,'r','sqrt(x^2 + y^2 + z^2)');

GradU = [diff(U,x); diff(U,y); diff(U,z)];
disp('Grad U is:')
pretty(GradU)

disp('A Simplified version is:')
GradU   = simplify(GradU);
pretty(subs(GradU,'sqrt(x^2 + y^2 + z^2)','r'))
end


 
%end
