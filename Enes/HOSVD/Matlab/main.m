function main(order)
    % Include necessary libraries and include files
    addpath(genpath([pwd '/includes']));
    addpath(genpath([pwd '/lib']));

    [nodes,weights] = int_gauss_weights(order+1,0,1);
    vertices = equidistant_points(order);

    MASS_TENSOR = mass_tensor(order,weights,nodes,vertices);
    MASS_MATRIX = tensor_to_matrix(MASS_TENSOR,order)
    [U, S, sv] = mlsvd(MASS_TENSOR);
    U_inverse = transpose(U);
    
    % Calculate Inverse of S.
    for i=1:order+1
        for j=1:order+1
            for k=1:order+1
                for l=1:order+1
                    if(S(i,j,k,l)!=0)
                        S_inverse(i,j,k,l)=1/S(i,j,k,l);
                    else 
                        S_inverse(i,j,k,l)=0;
                    endif
                endfor
            endfor
        endfor
    endfor
    
    MASS_TENSOR_INVERSE = lmlragen(U_inverse,S_inverse);
    MASS_MATRIX_PSEUDOINVERSE = tensor_to_matrix(MASS_TENSOR_INVERSE,order)

    % Solve linear system
    maxit = 10000;
    tol = 10^-25;
    restart = 1;
    b= rand((order+1)^2,1)
   
    gmres(MASS_MATRIX,b,restart,tol,maxit,MASS_MATRIX_PSEUDOINVERSE);    
end