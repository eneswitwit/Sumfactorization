function [active_node,local_node,active_cell_vector]=hf_remove_duplicates(active_node,local_node,active_cell_vector)
    % We wil check for duplicate global nodes, which get accessed by more than one cell
    for i=1:length(active_node)-1
        for j=i+1:length(active_node)
            % if two active nodes are the same, one will get ommited
            if (active_node(j)==active_node(i))
                active_node(j)=0;
                local_node(j)=0;
                active_cell_vector(j)=0;
            endif
        end
    endfor
    
    % Every duplicate was swapped with a zero columns
    % Octaves function 'unique' will ommit non-unique columns
    % If the first column is a zero column, it will be skipped
    Aux_matrix=unique([active_node',local_node',active_cell_vector'],"rows");
    if(Aux_matrix(1,1)==0)
        active_node=Aux_matrix(2:rows(Aux_matrix),1)';
        local_node=Aux_matrix(2:rows(Aux_matrix),2)';
        active_cell_vector=Aux_matrix(2:rows(Aux_matrix),3)';
    else
        active_node=Aux_matrix(1:rows(Aux_matrix),1)';
        local_node=Aux_matrix(1:rows(Aux_matrix),2)';
        active_cell_vector=Aux_matrix(1:rows(Aux_matrix),3)';
    endif
end