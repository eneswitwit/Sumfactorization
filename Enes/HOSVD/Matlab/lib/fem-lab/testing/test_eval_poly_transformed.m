function test_eval_poly_transformed()
    % This is a rather informal test for correctness
    % The function 'hf_eval_poly_transformed' is supposed to transform a shape function to the local shape function on a mesh cell.
    % To check if our function behaves as desired, we will compare the plots.
    SF=sf_generate(2);
    % x-axis values
    x=0:0.1:1;
    % h is our mesh size
    h=0.5;
    % scale the vector for plotting
    x_h=h*x;
    
    % displacement in x and y direction, given by the coordinates of the first vertex of a cell 
    dis_x=0.5;
    dis_y=0.5;
    % displace vector
    y_h=x_h+dis_y;
    x_h=x_h+dis_x;
    
    [X_h,Y_h]=meshgrid(x_h,y_h);
    
    for i=1:length(SF)
    figure(i)
    mesh(x_h,y_h,hf_eval_poly_transformed(X_h,Y_h,SF(i,:),h,dis_x,dis_y))
    
    [X,Y]=meshgrid(x,x);
    figure(i+length(SF))
    mesh(x,x,hf_eval_poly(X,Y,SF(i,:)))
    

end