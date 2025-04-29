function R = sphDFS_half(theta,phi,shape,nphi,nth)    
    % get coefficients
    as1 = shape(1:nth+1)';
    as2 = shape(nth+2:end)';
    % first sum
    R1 = sum(as1.*cos((0:nth)*theta));
    
    % second sum
    R2 = 0;
    for m = 1:nphi
        asm = as2((m-1)*nth+1:m*nth);
        if mod(m,2) == 0
            s = 1;
        else
            s = 0;
        end
        for l = 1:nth
            R2 = R2 + asm(l)*sin(theta)^s*sin(l*theta)*cos(m*phi);
        end
    end
	
    % add sums
    R = R1 + R2;
end