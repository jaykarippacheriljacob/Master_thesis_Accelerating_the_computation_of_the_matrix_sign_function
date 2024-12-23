function f_Hm = compute_sign_function(Hm)
    % Inputs: 
    %       Hm: Input matrix
    % Outputs:
    %       f_Hm: Resultant sign function applied to Hm  
    
    f_Hm = Hm * (sqrtm(inv(Hm * Hm)));
    % disp(f_Hm);
end
