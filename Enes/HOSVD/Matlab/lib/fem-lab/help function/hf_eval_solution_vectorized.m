function val=hf_eval_solution_vectorized(x,y,u,Cell,Vertex,SF)
    %Let (x,y) be the point of evaluation,
    %Let u be the solution of the linear system
    %Let Cell be the matrix, wich stores the numbers of the four vertices of the i_th cell in the i_th row
    %Let Vertex be the matrix, in which the coordinates of the i_th vertex are stored in the i_th row 
    %Let pol_deg, be the polynomial degree of the shape functions on the edges of a cell
    %Let SF be the matrix, which stores the coefficients of the shapefunctions on our reference square. Each row corresponds to one shapefunction
    
    %Calculate useful measurements, so we don't have to pass it as a parameter
    pol_deg=sqrt(rows(SF))-1;
    mesh_size=Vertex(2,1)-Vertex(1,1);
    cells_per_row=(1/mesh_size);
    nodes_per_row=(pol_deg*cells_per_row)+1;
    nodes_per_edge=pol_deg+1;
    
    %The following code will check, whether the point (x,y) is located on a node, an edge or the interior of a cell. If (x,y) is on an edge, we also check if it's on the boundary. 
    %If (x,y) is located on a node, we know that only one basis function phi_i is non-zero at this point. Therefore, the value at this point is determined by the coefficient u(i)
    %Let (x,y) be a point which is not located on a node. We know that most of the basis functions are zero at (x,y). Our goal is to find the non-zero basis functions. 
    %Once found, we can evaluate them and their respective coefficient u(i) and sum over all non-zero basisfunctions
    if (floor(x*(nodes_per_row-1))==x*(nodes_per_row-1) && floor(y*(nodes_per_row-1))==y*(nodes_per_row-1))
        node_i=x*(nodes_per_row-1);
        node_j=y*(nodes_per_row-1);
        node_no=node_j*nodes_per_row+node_i+1;
        val=u(node_no);
    else
        if (floor(x*cells_per_row)==x*cells_per_row)       
            cell_j=floor(y*cells_per_row);
            if(x==0)                                                    
                cell_i=0;
                active_cell=(cell_j)*cells_per_row+1;
            else                                                   
                if(x==1)                                                    
                    cell_i=cells_per_row-1;
                    active_cell=(cell_j+1)*cells_per_row;
                else                                                   
                    cell_i=[x*cells_per_row-1;x*cells_per_row];
                    active_cell=[cell_j*cells_per_row+cell_i(1)+1;cell_j*cells_per_row+cell_i(2)+1];
                    cell_j=[cell_j;cell_j];
                endif
            endif
        else                                            
            if (floor(y*cells_per_row)==y*cells_per_row)          
                cell_i=floor(x*cells_per_row);
                if (y==0)                               
                    cell_j=0;
                    active_cell=cell_i+1;
                else                                                    
                    if (y==1)                     
                        cell_j=cells_per_row-1;
                        active_cell=cells_per_row*(cells_per_row-1)+cell_i+1;
                    else                                                      
                        cell_j=[y*cells_per_row-1;y*cells_per_row];
                        active_cell=[cell_j(1)*cells_per_row+cell_i+1;cell_j(2)*cells_per_row+cell_i+1];
                        cell_i=[cell_i;cell_i];
                    endif
                endif
            else                                                    
                cell_i=floor(x*cells_per_row);
                cell_j=floor(y*cells_per_row);
                active_cell=cell_j*cells_per_row+cell_i+1;
            endif
        endif
        %At this point, we computed the cell(s), in which our point (x,y) lies. Now, we want to compute the nodes, that are located in these cells/this cell. Each node correlates to one basis function, which is non-zero at (x,y).
        for i=1:length(active_cell)
            for k=1:nodes_per_edge
                active_node((i-1)*nodes_per_edge^2+nodes_per_edge*(k-1)+1:(i-1)*nodes_per_edge^2+k*nodes_per_edge)=cell_j(i)*(nodes_per_edge-1)*nodes_per_row+cell_i(i)*(nodes_per_edge-1)+1+nodes_per_row*(k-1):cell_j(i)*(nodes_per_edge-1)*nodes_per_row+cell_i(i)*(nodes_per_edge-1)+nodes_per_edge+nodes_per_row*(k-1);
                local_node((i-1)*nodes_per_edge^2+nodes_per_edge*(k-1)+1:(i-1)*nodes_per_edge^2+k*nodes_per_edge)=nodes_per_edge*(k-1)+1:nodes_per_edge*k;
                active_cell_vector((i-1)*nodes_per_edge^2+nodes_per_edge*(k-1)+1:(i-1)*nodes_per_edge^2+k*nodes_per_edge)=active_cell(i);
            endfor
        endfor
        %Note, that the code above only converted the local node numbering to the global node numbering. It is clear, that two adjacent cells share global node points, which lie on the edge bewteen them. 
        %Therefore, we have to remove the duplicate nodes, since exactly one shape function correlates to one node
        [active_node,local_node,active_cell_vector]=hf_remove_duplicates(active_node,local_node,active_cell_vector);
        % We try to avoid any for loops, this is the reason why we went for this construction using the diagonal of a matrix.
        % The rows will store the computations for each cell
        % The columns will store the computations for each active shape function
        % The diagonal of this matrix consists of the correct pairs of shape function and cell and we can use the inner product with the right entries of u
        % If not for this construction, we would have had to loop over each entry of our active_node vector, evaluate a polynomial and multiply with a scalar.
        % Since matrix/vector operations are optimized in octave, this approach is more efficient.
        val=u(int32(active_node))'*diag(hf_eval_poly_transformed(x,y,SF(int32(local_node),:),mesh_size,Vertex(Cell(int32(active_cell_vector),1),1),Vertex(Cell(int32(active_cell_vector),1),2)));
    endif


    
endfunction