function m_plot_shape_functions()
    % Iterate over the polynomial degree
    for p=1:5
        polynomial_deg=p;
        % Generate coefficient matrix SF
        SF = sf_generate(polynomial_deg);
        % Initialize sample vector. use scalar to refine mesh for higher order polynomials
        x=0:0.05/p:1;
        y=x;
        % Generate rank one matrices
        [X,Y]=meshgrid(x,y);
        for k=1:rows(SF);
            % Compute shape function values and store in matrix
            Z=hf_eval_poly(X,Y,SF(k,:));
            h=figure(k);
            % actual plot function
            mesh(x,y,Z);
            % Settings for axis, since shape functions can get negative, and can be greater than
            % one for higher order polynomials
            axis ([0 1 0 1 0-(p-1)*0.25 1+(p-1)*0.25])
            xlabel ("x")
            ylabel ("y")
            zlabel ("p(x,y)")
            % Export plot as png with appropriate name
            print(h,'-dpng',num2str(k),[pwd '\Shape Function Plots\Shape Functions of Degree ' num2str(polynomial_deg) '\P' num2str(k)])
         endfor
    endfor
endfunction