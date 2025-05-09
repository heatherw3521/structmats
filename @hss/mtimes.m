function H = mtimes(H1,H2)

if isa(H1,'hss')
    if isa(H2,'hss')
        error('hss matmat multiplication is not yet supported')
    elseif isscalar(H2)
        error('hss mat scalar multiplication is not yet supported')
    else
        H = hss_matvec(H1,H2);
    end
else
    if isscalar(H1)
        error('hss mat scalar multiplication is not yet supported')
    else
        error('need to code transpose then this will work')
    end
end

end

