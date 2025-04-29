function radf = sph_harmonics(sph_w,theta,phi,nf,num_coeff)
    % save Slms in vector
    Slm = zeros(num_coeff,1);

    if ~isnumeric(theta)
        Slm = subs(Slm);
    end
    
    % start loop
    for l = 0:nf
        P = ass_legendre_list(l,cos(theta));
        S = zeros(2*l + 1,1);
        if ~isnumeric(theta)
            S = subs(S);
        end
        for m = -l:l
            mtemp = abs(m);
            Plm = P(mtemp+1);
            Nlm = sqrt((2*l + 1)/(4*pi)*prod(1:l - mtemp)/prod(1:l + mtemp));
            S(m+l+1) = Nlm*Plm*sqrt(2)*cos(mtemp*phi);
        end
        Slm((l+1)^2-2*(l+1)+2:(l+1)^2) = S;
    end
    radf = sph_w'*Slm;
end