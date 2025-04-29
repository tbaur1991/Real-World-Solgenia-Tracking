function ps = as_polygon_super(X,nPoly)
    % Returns the polygon chain that corresponds to the given x.
    a = c1(X(7),0,'lower'); 
    b = c1(X(8),0,'lower'); 
    e = c1(X(10),1,'lower');
    
    as = (0:360/(nPoly-1):360) * pi / 180;
    ps = [a;b].*sign([cos(as); sin(as)]).*abs([cos(as); sin(as)]).^(2/e);
end