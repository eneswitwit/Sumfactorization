function val = error_L2(f,g,order)
    
    % Initialize gauss weights/points
    [x,w]=int_gauss_weights(order,0,1);
    % Initialize argument of integration
    h=@(x,y) (f(x,y)-g(x,y))^2;
    % Squareroot of integration yields the L2 norm
    val=sqrt(int_gauss(x,w,x,w,h));
    
endfunction