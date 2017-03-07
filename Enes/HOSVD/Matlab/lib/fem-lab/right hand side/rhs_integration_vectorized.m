function rhs = rhs_integration_vectorized(Vertex,Cell,SF,f)
    %Let Cell be the matrix, which stores the vertices for each cell
    %Let Vertex be the matrix, which stores the coordinates for each vertex
    %Let SF be the matrix, containing the coefficients of the shape funtions
    %Let f be the right hand side of the strong form
    
    %This function will give us the right hand side of the linear system
    
    % Initialize Gauss Quadratur
    [sample_points,weights] = int_gauss_weights(5,0,1);
    % Useful computations for later use
    mesh_size=Vertex(2,1)-Vertex(1,1);
    pol_deg = sqrt(length(SF))-1;
    cells_per_row=(1/mesh_size);
    number_of_cells=cells_per_row^2;
    number_of_nodes=(pol_deg*cells_per_row+1)^2;
    
    % Initialize matrix "Nodes" wich contains the global node number for each cell
    Nodes=mesh_nodes(mesh_size,pol_deg);
    
    % Initialize right hand side vector
    rhs=zeros(number_of_nodes,1);
    
    % Iterate over all cells. For one cell, we integrate every shape function over the reference cell at once
    % We use the matrix Nodes to convert the local node numbering into the global numbering
    for k=1:number_of_cells
        % the first vertex of the active cell will give us the displacement of our transformation
        active_vertex=Vertex(Cell(k,1),:);
        % the function g will evaluate the transformed right hand side for every shape function on a given cell
        % we use the octave in-build function "repmat", since "f" shall be multiplied to every shape function
        g =  @(x,y) repmat(f(mesh_size*x+active_vertex(1),mesh_size*y+active_vertex(2)),length(SF),1).*hf_eval_poly(x,y,SF);
        % "Nodes" will select the correct entries of the vector, we only need to add the vector-valued integral
        rhs(int32(Nodes(int32(k),:)))+=int_gauss_vectorized_matrices(sample_points,weights,sample_points,weights,g);
    endfor
    % The transformation formula requires to multiply with the jacobi determinant
    rhs*=(mesh_size^2);
endfunction