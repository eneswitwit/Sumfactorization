function [cell_no,sf_no] = node_to_cells(node_no,mesh_size,pol_deg)
    % Let node_no be the index of a node (global numbering).
    % This function shall help us to construct our basis functions.
    % For a given index i, this function will give us the numbers of mesh-cells [cell_no] and their respective shape function [sf_no]. 
    % When 'glued together', this will yield our i_th basis function.
    
    % Useful computations
    cells_per_row=(1/mesh_size);
    nodes_per_row=(pol_deg*cells_per_row)+1;
    nodes_per_edge=pol_deg+1;
    node_j=0;
    
    % At first, we transform the lexicographic ordering into a x- and y-coordinate.
    while (node_j+1)*nodes_per_row<node_no
        node_j++;
    endwhile
    node_i=node_no-node_j*nodes_per_row-1;
    
    % The following code consists of many tests, to determine where the node is located on the mesh.
    % Generally, we test if the node has the same x- and/or y-coordinate as a vertex from our cell.
    % If both coordinates match, the node is located on a vertex, which means the respective basis function is non-zero on all adjacent cells.
    % If only one coordinate is a match, the node lies on an edge, the adjacent cells are therefore non-zeros.
    % If there aren't any matches, the node lies in the interior of a cell, this cell will be the only non-zero cell.
    % Additionally, in the first two cases we have to check if the node lies on the boundary, since this reduces the number of cells used for computation
    % Given the location of the node, it remains a combinatoric  exercise to determine the active cell/shape function.
    
    % Instead of compairing the actual coordinates, we again use a combinatoric approach since we know how many nodes we have on any edge,
    % Let k be the polynomial degree, then every k-th node in x-direction will match with a vertex' x-coordinate, same for y-direction.
    if ((floor(node_i/pol_deg)==node_i/pol_deg))
        if (floor(node_j/pol_deg)==node_j/pol_deg)  
            if (node_i==0)
                if (node_j==0)
                    cell_no=[1];
                    sf_no=[1];
                else
                    if(node_j==nodes_per_row-1)
                        cell_no=[cells_per_row*(cells_per_row-1)+1];
                        sf_no=[nodes_per_edge^2-pol_deg];
                    else
                        cell_no=[((node_j/pol_deg)-1)*cells_per_row+(node_i/pol_deg)+1;((node_j/pol_deg))*cells_per_row+(node_i/pol_deg)+1];
                        sf_no=[nodes_per_edge^2-pol_deg;1];
                    endif
                endif
            else
                if (node_i==nodes_per_row-1)
                    if (node_j==0)
                        cell_no=[cells_per_row];
                        sf_no=[nodes_per_edge];
                    else
                        if (node_j==nodes_per_row-1)
                            cell_no=[cells_per_row^2];
                            sf_no=[nodes_per_edge^2];
                        else
                            cell_no=[((node_j/pol_deg)-1)*cells_per_row+(node_i/pol_deg);((node_j/pol_deg))*cells_per_row+(node_i/pol_deg)];
                            sf_no=[nodes_per_edge^2;1+pol_deg];
                        endif
                    end    
                else
                    if (node_j==0)
                        cell_no=[((node_j/pol_deg))*cells_per_row+(node_i/pol_deg);((node_j/pol_deg))*cells_per_row+(node_i/pol_deg)+1];
                        sf_no=[1+pol_deg;1];
                    else
                        if (node_j==nodes_per_row-1)
                            cell_no=[((node_j/pol_deg)-1)*cells_per_row+(node_i/pol_deg);((node_j/pol_deg)-1)*cells_per_row+(node_i/pol_deg)+1];
                            sf_no=[nodes_per_edge^2;nodes_per_edge^2-pol_deg];
                        else
                            cell_no(1,1)=((node_j/pol_deg)-1)*cells_per_row+(node_i/pol_deg);
                            cell_no(2,1)=((node_j/pol_deg)-1)*cells_per_row+(node_i/pol_deg)+1;
                            cell_no(3,1)=((node_j/pol_deg))*cells_per_row+(node_i/pol_deg);
                            cell_no(4,1)=((node_j/pol_deg))*cells_per_row+(node_i/pol_deg)+1;
                            sf_no(1,1)=nodes_per_edge^2;
                            sf_no(2,1)=nodes_per_edge^2-pol_deg;
                            sf_no(3,1)=1+pol_deg;
                            sf_no(4,1)=1;
                        endif
                    endif
                endif
            endif    
        else
            if (node_i==0)
                sf_j=(((node_j)/pol_deg)-floor((node_j)/pol_deg))*pol_deg;
                cell_no=[(floor(node_j/pol_deg))*cells_per_row+(node_i/pol_deg)+1];
                sf_no=[1+(nodes_per_edge*sf_j)];
            else
                if (node_i==nodes_per_row-1)
                    sf_j=(((node_j)/pol_deg)-floor((node_j)/pol_deg))*pol_deg;
                    cell_no=[(floor(node_j/pol_deg))*cells_per_row+(node_i/pol_deg)];
                    sf_no=[(nodes_per_edge)*(sf_j+1)];
                else
                    sf_j=(((node_j)/pol_deg)-floor((node_j)/pol_deg))*pol_deg;
                    cell_no=[(floor(node_j/pol_deg))*cells_per_row+(node_i/pol_deg);(floor(node_j/pol_deg))*cells_per_row+(node_i/pol_deg)+1];
                    sf_no=[(nodes_per_edge)*(sf_j+1);1+(nodes_per_edge*sf_j)];
                endif
            endif
        endif
    else
        sf_i=((node_i/pol_deg)-floor((node_i)/pol_deg))*pol_deg;
        if (floor(node_j/pol_deg)==node_j/pol_deg)
            if (node_j/pol_deg==0)
                cell_no=((node_j/pol_deg))*cells_per_row+floor(node_i/pol_deg)+1;
                sf_no=sf_i+1;
            else
                if (node_j==nodes_per_row-1)
                    cell_no=((node_j/pol_deg)-1)*cells_per_row+floor(node_i/pol_deg)+1;
                    sf_no=nodes_per_edge*(nodes_per_edge-1)+sf_i+1  ;
                else
                    cell_no=[((node_j/pol_deg)-1)*cells_per_row+floor(node_i/pol_deg)+1;((node_j/pol_deg))*cells_per_row+floor(node_i/pol_deg)+1];
                    sf_no=[nodes_per_edge*(nodes_per_edge-1)+sf_i+1;sf_i+1]         ;   
                endif
            endif
        else
            sf_j=(((node_j)/pol_deg)-floor((node_j)/pol_deg))*pol_deg;
            cell_no=(floor((node_j)/pol_deg))*cells_per_row+floor((node_i)/pol_deg)+1;
            sf_no=sf_j*nodes_per_edge+sf_i+1;
        endif
    endif
    
    % This was added, because octave would give outputs that weren't integer.
    cell_no=round(cell_no);
    sf_no=round(sf_no);
    
endfunction