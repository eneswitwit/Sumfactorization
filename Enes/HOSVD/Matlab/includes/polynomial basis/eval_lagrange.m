function val = eval_lagrange(i,x,vertex)
    val = 1;
    for j=1:size(vertex)
        if i ~= j
            val = val *  ((x - vertex(j)) / (vertex(i) - vertex(j)));
        end
    end
end