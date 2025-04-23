function out = leafbuild(H)
% build the leaf component H
%%
if ~H.isleaf
    error('not a leaf')
elseif H.isdiag
    out = H.D; 
else
    Z = H.Z; 
    Y = H.Y; 
    DD = H.lrcomponent; 
    out= Z*DD*Y;
end
