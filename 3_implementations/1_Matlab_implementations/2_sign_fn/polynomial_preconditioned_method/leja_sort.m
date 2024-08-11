function v = leja_sort(v)
    n = size(v, 1);
    for i = 1:n-1  
        for j = i:n
            if(abs(v(j)) > abs(v(i)))
                temp = v(j);
                v(j) = v(i);
                v(i) = temp;
            end
        end
    end
end