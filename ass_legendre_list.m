function P = ass_legendre_list(l,x)
    P = zeros(l+1,1);

    if ~isnumeric(x)
        P = subs(P);
    end

    if l == 0
        P(1) = 1;
    elseif l == 1
        P(1) = x;
        P(2) = -(1 - x^2)^0.5;
    elseif l == 2
        P(1) = 0.5*(3*x^2 - 1);
        P(2) = -3*x*(1 - x)^0.5;
        P(3) = 3*(1 - x^2);
    elseif l == 3
        P(1) = 0.5*x*(5*x^2 - 3);
        P(2) = 1.5*(1 - 5*x^2)*(1 - x^2)^0.5;
        P(3) = 15*x*(1 - x^2);
        P(4) = -15*(1 - x^2)^1.5;
    elseif l == 4
        P(1) = 1/8*(35*x^4 - 30*x^2 + 3);
        P(2) = -5/2*(7*x^3 -3*x)*(1 - x^2)^0.5;
        P(3) = 15/2*(7*x^2 - 1)*(1 - x^2);
        P(4) = -105*x*(1 - x^2)^1.5;
        P(5) = 105*(1 - x^2)^2;
    elseif l == 5
        P(1) = 1/8*(15*x - 70*x^3 + 63*x^5);
        P(2) = -15/8*sqrt(1 - x^2)*(1 - 14*x^2 + 21*x^4);
        P(3) = -105/2*(-1 + x^2)*(-x + 3*x^3);
        P(4) = -105/2*(1 - x^2)^(3/2)*(-1 + 9*x^2);
        P(5) = 945*x*(-1 + x^2)^2;
        P(6) = -945*(1 - x^2)^(5/2);
    elseif l == 6
        P(1) = 1/16*(-5 + 105*x^2 - 315*x^4 + 231*x^6);
        P(2) = -21/8*sqrt(1 - x^2)*(5*x - 30*x^3 + 33*x^5);
        P(3) = -105/8*(-1 + x^2)*(1 - 18*x^2 + 33*x^4);
        P(4) = -315/2*(1 - x^2)^(3/2)*(-3*x + 11*x^3);
        P(5) = 945/2*(-1 + x^2)^2*(-1 + 11*x^2);
        P(6) = -10395*x*(1 - x^2)^(5/2);
        P(7) = -10395*(-1 + x^2)^3;
    elseif l == 7
        P(1) = 1/16*(-35*x + 315*x^3 - 693*x^5 + 429*x^7);
        P(2) = -7/16*sqrt(1 - x^2)*(-5 + 135*x^2 - 495*x^4 + 429*x^6);
        P(3) = -63/8*(-1 + x^2)*(15*x - 110*x^3 + 143*x^5);
        P(4) = -315/8*(1 - x^2)^(3/2)*(3 - 66*x^2 + 143*x^4);
        P(5) = 3465/2*(-1 + x^2)^2*(-3*x + 13*x^3);
        P(6) = -10395/2*(1 - x^2)^(5/2)*(-1 + 13*x^2);
        P(7) = -135135*x*(-1 + x^2)^3;
        P(8) = -135135*(1 - x^2)^(7/2);
    elseif l == 8
        P(1) = 1/128*(35 - 1260*x^2 + 6930*x^4 - 12012*x^6 + 6435*x^8);
        P(2) = -9/16*sqrt(1 - x^2)*(-35*x + 385*x^3 - 1001*x^5 + 715*x^7);
        P(3) = -315/16*(-1 + x^2)*(-1 + 33*x^2 - 143*x^4 + 143*x^6);
        P(4) = -3465/8*(1 - x^2)^(3/2)*(3*x - 26*x^3 + 39*x^5);
        P(5) = 10395/8*(-1 + x^2)^2*(1 - 26*x^2 + 65*x^4);
        P(6) = -135135/2*(1 - x^2)^(5/2)*(-x + 5*x^3);
        P(7) = -135135/2*(-1 + x^2)^3*(-1 + 15*x^2);
        P(8) = -2027025*x*(1 - x^2)^(7/2);
        P(9) = 2027025*(-1 + x^2)^4;
    elseif l == 9
        P(1) = 1/128*(315*x - 4620*x^3 + 18018*x^5 - 25740*x^7 + 12155*x^9);
        P(2) = -45/128*sqrt(1 - x^2)*(7 - 308*x^2 + 2002*x^4 - 4004*x^6 + 2431*x^8);
        P(3) = -495/16*(-1 + x^2)*(-7*x + 91*x^3 - 273*x^5 + 221*x^7);
        P(4) = -3465/16*(1 - x^2)^(3/2)*(-1 + 39*x^2 - 195*x^4 + 221*x^6);
        P(5) = 135135/8*(-1 + x^2)^2*(x - 10*x^3 + 17*x^5);
        P(6) = -135135/8*(1 - x^2)^(5/2)*(1 - 30*x^2 + 85*x^4);
        P(7) = -675675/2*(-1 + x^2)^3*(-3*x + 17*x^3);
        P(8) = -2027025/2*(1 - x^2)^(7/2)*(-1 + 17*x^2);
        P(9) = 34459425*x*(-1 + x^2)^4;
        P(10) = -34459425*(1 - x^2)^(9/2);
    elseif l == 10
        P(1) = 1/256*(-63 + 3465*x^2 - 30030*x^4 + 90090*x^6 - 109395*x^8 + 46189*x^10);
        P(2) = -55/128*sqrt(1 - x^2)*(63*x - 1092*x^3 + 4914*x^5 - 7956*x^7 + 4199*x^9);
        P(3) = -495/128*(-1 + x^2)*(7 - 364*x^2 + 2730*x^4 - 6188*x^6 + 4199*x^8);
        P(4) = -6435/16*(1 - x^2)^(3/2)*(-7*x + 105*x^3 - 357*x^5 + 323*x^7);
        P(5) = 45045/16*(-1 + x^2)^2*(-1 + 45*x^2 - 255*x^4 + 323*x^6);
        P(6) = -135135/8*(1 - x^2)^(5/2)*(15*x - 170*x^3 + 323*x^5);
        P(7) = -675675/8*(-1 + x^2)^3*(3 - 102*x^2 + 323*x^4);
        P(8) = -11486475/2*(1 - x^2)^(7/2)*(-3*x + 19*x^3);
        P(9) = 34459425/2*(-1 + x^2)^4*(-1 + 19*x^2);
        P(10) = -654729075*x*(1 - x^2)^(9/2);
        P(11) = -654729075*(-1 + x^2)^5;
    elseif l == 11
        P(1) = 1/256*(-693*x + 15015*x^3 - 90090*x^5 + 218790*x^7 - 230945*x^9 + 88179*x^11);
        P(2) = -33/256*sqrt(1 - x^2)*(-21 + 1365*x^2 - 13650*x^4 + 46410*x^6 - 62985*x^8 + 29393*x^10);
        P(3) = -2145/128*(-1 + x^2)*(21*x - 420*x^3 + 2142*x^5 - 3876*x^7 + 2261*x^9);
        P(4) = -45045/128*(1 - x^2)^(3/2)*(1 - 60*x^2 + 510*x^4 - 1292*x^6 + 969*x^8);
        P(5) = 135135/16*(-1 + x^2)^2*(-5*x + 85*x^3 - 323*x^5 + 323*x^7);
        P(6) = -135135/16*(1 - x^2)^(5/2)*(-5 + 255*x^2 - 1615*x^4 + 2261*x^6);
        P(7) = -2297295/8*(-1 + x^2)^3*(15*x - 190*x^3 + 399*x^5);
        P(8) = -34459425/8*(1 - x^2)^(7/2)*(1 - 38*x^2 + 133*x^4);
        P(9) = 654729075/2*(-1 + x^2)^4*(-x + 7*x^3);
        P(10) = -654729075/2*(1 - x^2)^(9/2)*(-1 + 21*x^2);
        P(11) = -13749310575*x*(-1 + x^2)^5;
        P(12) = -13749310575*(1 - x^2)^(11/2);
    elseif l == 12
        P(1) = (231 - 18018*x^2 + 225225*x^4 - 1021020*x^6 + 2078505*x^8 - 1939938*x^10 + 676039*x^12)/1024;
        P(2) = -39/256*sqrt(1 - x^2)*(-231*x + 5775*x^3 - 39270*x^5 + 106590*x^7 - 124355*x^9 + 52003*x^11);
        P(3) = -3003/256*(-1 + x^2)*(-3 + 225*x^2 - 2550*x^4 + 9690*x^6 - 14535*x^8 + 7429*x^10);
        P(4) = -15015/128*(1 - x^2)^(3/2)*(45*x - 1020*x^3 + 5814*x^5 - 11628*x^7 + 7429*x^9);
        P(5) = 135135/128*(-1 + x^2)^2*(5 - 340*x^2 + 3230*x^4 - 9044*x^6 + 7429*x^8);
        P(6) = -2297295/16*(1 - x^2)^(5/2)*(-5*x + 95*x^3 - 399*x^5 + 437*x^7);
        P(7) = -2297295/16*(-1 + x^2)^3*(-5 + 285*x^2 - 1995*x^4 + 3059*x^6);
        P(8) = -130945815/8*(1 - x^2)^(7/2)*(5*x - 70*x^3 + 161*x^5);
        P(9) = 654729075/8*(-1 + x^2)^4*(1 - 42*x^2 + 161*x^4);
        P(10) = -4583103525/2*(1 - x^2)^(9/2)*(-3*x + 23*x^3);
        P(11) = -13749310575/2*(-1 + x^2)^5*(-1 + 23*x^2);
        P(12) = -316234143225*x*(1 - x^2)^(11/2);
        P(13) = 316234143225*(-1 + x^2)^6;
    end
end