function H = mldivide(H1, H2)
% backslash

if isa(H1,'hss')
    if isscalar(H2)
        % scalar mult
        error('hss scalar mat multiplication is not yet supported')
        %H = 1/H2* H1;
    elseif isa(H2, 'hss')
        error('hss hss backsolver is not yet supported')
        %H = hss_mldivide(H1, H2);
    else
        H = hss_vecsolve(H1,H2);
    end
else
    error('hss inverse is not yet supported')
    % H = inv(H2)/inv(H1);
end
end



