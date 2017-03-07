function test_eval_poly_transformed()
    % This is a rather informal test for correctness
    % 
    SF=sf_generate(2);
    x=0:0.1:1;
    h=0.5;
    x_h=h*x;
    dis_x=0.5;
    dis_y=0.5;
    x_h=x_h+dis_x;
    [X_h,Y_h]=meshgrid(x_h,x_h);
    
    for i=1:length(SF)
    figure(i)
    mesh(x_h,x_h,hf_eval_poly_transformed(X_h,Y_h,SF(i,:),h,dis_x,dis_y))
    
    hf_eval_poly_transformed(0.7,0.7,SF,h,dis_x,dis_y)
    
    [X,Y]=meshgrid(x,x);
    figure(i+9)
    mesh(x,x,hf_eval_poly(X,Y,SF(i,:)))
    

end