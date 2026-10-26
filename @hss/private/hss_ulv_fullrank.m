function tf = hss_ulv_fullrank(F)
%HSS_ULV_FULLRANK  True when the ULV factors F (of a square or wide W) show
%   that W has full row rank, i.e. the reduced root block kept every row.
while ~F.isroot
    F = F.next;
end
tf = ~(isfield(F, 'Uperp') && ~isempty(F.Uperp));
end
