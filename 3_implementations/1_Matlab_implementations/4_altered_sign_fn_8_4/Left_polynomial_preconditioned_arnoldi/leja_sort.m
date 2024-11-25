function v1 = leja_sort(v)
    %% Returns sorted vector based on Leja sort.
    % Inputs:
    %   v - A vector of numerical values to be sorted.
    %
    % Outputs:
    %   v1 - The sorted vector, arranged in descending order based on the absolute values of the elements.

    % Get the number of elements in the vector
    n1 = size(v, 1); % Size of the sorted ritz values.
    v1 = zeros(n1, 1);
    
    n = n1 - 1; % Variying size of the input vector upon sorting cycle

    i = 1;
    [~, idx] = max(abs(v)); % Finding the index and absolute maximum of the 
                            % ritz values
    v1(i) = v(idx); % Saving the absolute maximum theta found to the new 
                    % vector at position i

    v(idx) = []; % Removing the theta that has been found to be absolute 
                 % maximum from the vector

    pdt = ones(n, 1); % Initializing a varying vector for finding the 
                      % product of distance wrt sorted theta's.

    % Starting the for loop for the leja sorting
    for i = 2:n1

        % Evaluating the product of distance between the theta's wrt the
        % sorted theta's
        pdt = pdt .* (abs(v1(i-1) .* ones(n, 1) - v));
        idx = 1; % Setting initial index to 1
        for j = 2:n
            % Checking for a index with larger distance product
            if(pdt(j) > pdt(idx))
                idx = j;
            end
        end

        % Insertion of the found theta and removing it from the orginal
        % vector
        v1(i) = v(idx);
        v(idx) = [];
        pdt(idx) = [];
        n = n - 1;
    end
end