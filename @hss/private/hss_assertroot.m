function hss_assertroot(H, what)
%HSS_ASSERTROOT  Error unless H is a whole HSS matrix (a root), not a piece of one.
%   The pieces H.A11, H.A12, ... of the tree keep the global row/column
%   numbering of the matrix they came from, so used as matrices on their own
%   they give wrong results (audit B05). Methods that treat H as a matrix
%   call this first.
if isempty(H.isroot) || ~H.isroot
    error('hss:notRoot', ['%s needs a whole HSS matrix; this object is a piece of one ' ...
        '(e.g. H.A22), whose indices refer to the full matrix. Use full(...) for ' ...
        'its dense block.'], what);
end
end
