function y = motion_model(x,T)
    % Memory Allocation
    y = zeros(size(x));
    
    if ~x(6) == 0
        y(1) = x(1) + 2*x(4)/x(6)*sin(x(6)*T/2)*cos(x(5) + x(6)*T/2); 
        y(2) = x(2) + 2*x(4)/x(6)*sin(x(6)*T/2)*sin(x(5) + x(6)*T/2);
    else
        y(1) = x(1) + cos(x(5))*x(4)*T;
        y(2) = x(2) + sin(x(5))*x(4)*T;
    end
    y(3:4) = x(3:4);
    y(5) = x(5) + x(6)*T;
    y(6:end) = x(6:end);
end
