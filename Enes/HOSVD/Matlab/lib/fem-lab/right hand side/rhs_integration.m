function rhs = rhs_integration(Vertex,Cell,SF,f)
    % Let cell be the matrix, which stores the vertices for each cell.
    % Let vertex be the matrix, which stores the coordinates for each vertex.
    % Let SF be the matrix, containing the coefficients of the shape funtions.
    
    % This function will give us the right hand side of the linear system.
    % Initialize gauss quadratur.
    [sample_points,weights] = int_gauss_weights(10,0,1);
    
    % Useful computations for later use.
    mesh_size=Vertex(2,1)-Vertex(1,1);
    pol_deg = sqrt(length(SF))-1;
    cells_per_row=int32(1/mesh_size);
    number_of_nodes=int32(((pol_deg*cells_per_row)+1)^2);
    
    % Each node on our mesh represents one basis function. 
    % This means there will be only a small number of cells, which are non-zero for any given integration over the product of the basis functions and the right hand side .
    % We will denote these cells as 'active' cells, the respective shape function on the cell will be denoted as active shape function.
    for i=1:number_of_nodes
        RHS=0;
        % Our function 'node_to_cell' will tell us for a given node, which shape function we have to evaluate on which cell.
        [active_cell,active_sf]=node_to_cells(i,mesh_size,pol_deg);
        for k=1:rows(active_cell)
            active_cell_no = active_cell(k);
            active_vertex = Vertex(Cell(active_cell_no),1),:);
            active_sf_no=active_sf(k);
            % This transformation is further explained in our documentation
            g =  @(x,y) f(mesh_size*x+active_vertex(1),mesh_size*y+active_vertex(2)).*hf_eval_poly(x,y,SF(active_sf_no,:));
            integral=int_gauss(sample_points,weights,sample_points,weights,g);
            RHS=RHS+integral;
        endfor
        rhst(i)=RHS;
    endfor

    rhs=mesh_size^2*rhst';
    
endfunction