function MATRIX = tensor_to_matrix(TENSOR,order)
    for i1 = 1:order+1
        for i2 =1:order+1
            for j1 =1:order+1
                for j2 = 1:order+1
                    k = transform_2d(i1-1, i2-1, j1-1, j2-1, order);
                    MATRIX(k(1)+1,k(2)+1) = TENSOR(i1,i2,j1,j2);
                end
            end
        end
    end
end