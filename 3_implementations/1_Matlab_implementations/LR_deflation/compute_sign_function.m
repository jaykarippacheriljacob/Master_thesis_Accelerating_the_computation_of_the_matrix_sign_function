function f_Hm = compute_sign_function(Hm)
    % Inputs: 
    %       Hm: Input matrix
    % Outputs:
    %       f_Hm: Resultant sign function applied to Hm  
    
    f_Hm = Hm * (sqrtm(inv(Hm * Hm)));
    %f_Hm = Hm ./ (Hm .* Hm).^(1/2);
    % disp(f_Hm);
end
