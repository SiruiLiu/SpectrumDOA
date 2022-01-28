function [CRLB] = calCRLB(SNR,T,N,r,Lambda)
    CRLB = ((2.*Lambda^2)./(4.*pi^2.*T.*SNR.*r^2.*N))*180/pi;
end

