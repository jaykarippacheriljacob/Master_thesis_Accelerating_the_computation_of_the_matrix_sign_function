function matrixData = read_matrix(filename)
    % Input: 
    %      filename - Name of the file in which matrix is present.
    
    % Check if file exists
    if exist(filename, 'file') ~= 2
        error('File does not exist.');
    end
    
    % Load the data from the file
    loadedData = load(filename);
    
    % Extract the matrix from the loaded data
    matrixData = loadedData.D;
end


