function v = leja_sort(v)
    % Inputs:
    %   v - A vector of numerical values to be sorted.
    %
    % Outputs:
    %   v - The sorted vector, arranged in descending order based on the absolute values of the elements.

    % Get the number of elements in the vector
    n = size(v, 1);

    % Outer loop: iterate over each element except the last one
    for i = 1:n-1  
        % Inner loop: iterate over the elements starting from the current element in the outer loop
        for j = i:n
            % Compare the absolute value of the current element with the absolute value of the element in the outer loop
            if(abs(v(j)) > abs(v(i)))
                % Swap the elements if the current element's absolute value is greater
                temp = v(j);
                v(j) = v(i);
                v(i) = temp;
            end
        end
    end
end