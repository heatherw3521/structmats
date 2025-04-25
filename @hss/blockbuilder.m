function A = blockbuilder(H,Hparent)
% builds Hblock
% H is the piece we wish to reconstruct
% Hparent is the parent of H

if H.isleaf
    A = leafbuild(H);
    return
end


if H.Ir == Hparent.A12.Ir
    left = leftdrill(Hparent.A11);
    right = rightdrill(Hparent.A22);
else
    left = leftdrill(Hparent.A22);
    right = rightdrill(Hparent.A11);
end


A = left * H.Z * H.lrcomponent * H.Y * right;


end

function left = leftdrill(H)

if H.A11.isleaf
    left = blkdiag(H.A12.Z, H.A21.Z);
else
    left = blkdiag(leftdrill(H.A11), leftdrill(H.A22)) * blkdiag(H.A12.Z, H.A21.Z);
end

end


function right = rightdrill(H)

if H.A11.isleaf
    right = blkdiag(H.A21.Y, H.A12.Y);
else
    right  = blkdiag(H.A21.Y, H.A12.Y)* blkdiag(rightdrill(H.A11), rightdrill(H.A22));
end

end