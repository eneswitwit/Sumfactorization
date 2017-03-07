function val = hf_eval_poly_transformed(x,y,coefficient_vector,scaling,displacement_x,displacement_y)
    j=1;        % order of basis in the coefficient vector
    level=0;    % the level represents the structure of the ordering of the basisfunctions of the polynomial. It could be described
    val=0;      % as the highest exponent allowed at a specific moment.
    %The following code is almost the exact same code as in hf_eval_poly. We substituted x and y in a way, such that our polynomial is defined on a mesh cell rather than our reference cell. 
    while j < columns(coefficient_vector)
        exponent_x=level;
        exponent_y=0;
        while(exponent_y < level)
            val= val + coefficient_vector(:,j).*(scaling.^(-1).*(x-displacement_x')).^(exponent_x).*(scaling.^(-1).*(y-displacement_y')).^(exponent_y);
            exponent_y++;
            j++;
        endwhile
        exponent_x = 0;
        exponent_y = level;
        while(exponent_x < level)
            val= val+ coefficient_vector(:,j).*(scaling.^(-1).*(x-displacement_x')).^(exponent_x).*(scaling.^(-1).*(y-displacement_y')).^(exponent_y);
            exponent_x++;
            j++;
        endwhile
        val=val+coefficient_vector(:,j).*(scaling.^(-1).*(x-displacement_x')).^(level).*(scaling.^(-1).*(y-displacement_y')).^(level);
        j++;
        level++;
    endwhile
endfunction