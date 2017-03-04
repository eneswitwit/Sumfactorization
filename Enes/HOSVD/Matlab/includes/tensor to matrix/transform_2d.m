function k = transform_2d(i1,i2,j1,j2,order)
    k(1) = transform_1d(i1, i2, order);
    k(2) = transform_1d(j1, j2, order);
end