function out = leafbuild(H)
%LEAFBUILD  Dense block of a leaf node: D for a diagonal leaf, Z*B*Y otherwise.
if ~H.isleaf
    error('hss:leafbuild:notLeaf', 'leafbuild needs a leaf node.')
elseif H.isdiag
    out = H.D; 
else
    Z = H.Z; 
    Y = H.Y; 
    DD = H.lrcomponent; 
    out= Z*DD*Y;
end
