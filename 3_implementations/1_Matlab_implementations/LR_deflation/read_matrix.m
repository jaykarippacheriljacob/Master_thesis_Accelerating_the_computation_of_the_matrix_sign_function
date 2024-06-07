function matrixData = read_matrix(filename)
    % Check if file exists
    if exist(filename, 'file') ~= 2
        error('File does not exist.');
    end
    
    % Load the data from the file
    loadedData = load(filename);
    
    % Extract the matrix from the loaded data
    matrixData = loadedData.D;
end


