function [solution, residual] = ls_solve(A,b)
    % Test whether length of matrix is bigger than 10201. This value is computed by linear_solver_analytics.m.
    % That value tells us that GMRES is better for linear systems of the size 10201 ^2 and bigger. Otherwise MINRES is the better choice.
    if(length(A) > 10201)
        solution = gmres(A,b);
        residual = b-A*solution;
     else
        [solution,residual] = ls_minres(A,b,b,100,10^(-10));
     endif
end