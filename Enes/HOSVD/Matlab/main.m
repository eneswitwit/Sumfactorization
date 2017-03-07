function main(order)
    % Include necessary libraries and include files
    addpath(genpath([pwd '/includes']));
    addpath(genpath([pwd '/lib']));

    [nodes,weights] = int_gauss_weights(order+1,0,1);

    MASS_TENSOR = mass_tensor(order,weights,nodes,nodes);
    MASS_MATRIX = tensor_to_matrix(MASS_TENSOR,order)
end