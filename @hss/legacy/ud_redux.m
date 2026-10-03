function [Hshrunk,Ss] = ud_redux(H,options)


arguments
    H
    options.incase = 1;
end


[Hshrunk,Ss] = down_recurse(H)


end



function [Hshrunk,Ss] = down_recurse(H)

if H.A11.isleaf

else
    [Hshrunk.A11]
end

end