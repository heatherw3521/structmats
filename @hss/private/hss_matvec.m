function y = hss_matvec(H,v)
%HSS_MATVEC  y = H*v for an m x n HSS matrix H and an n x k dense v.
%   Upward pass: compress v through the column bases (ascend); downward
%   pass: expand through the row bases and add the leaf products (descend).
if H.isleaf == false
    dict = dictionary();
    % upward pass: column coefficients of every node
    dict = ascend(H,v,dict);

    % downward pass: row coefficients, then the leaf products
    y = descend(H,v,dict);
    
elseif H.isdiag %leaf diag
    y = H.D * v;
else %leaf off-diag
    y = H.Z*H.lrcomponent*H.Y*v;
end


end

function dict = ascend(H,v,dict)


if H.A11.isleaf == false
    dict = ascend(H.A11,v(1:H.A11.sz(2),:),dict);
    dict = ascend(H.A22,v(H.A11.sz(2)+1:end,:),dict);
    d1 = dict({[H.A11.A21.level,H.A11.A21.coltreeindex]});
    d2 = dict({[H.A11.A12.level,H.A11.A12.coltreeindex]});
    d3 = dict({[H.A22.A21.level,H.A22.A21.coltreeindex]});
    d4 = dict({[H.A22.A12.level,H.A22.A12.coltreeindex]});
    dict({[H.A21.level,H.A21.coltreeindex]}) = {H.A21.Y*[d1{1};d2{1}]};
    dict({[H.A12.level,H.A12.coltreeindex]}) = {H.A12.Y*[d3{1};d4{1}]};
    if H.isroot
        return
    end
else
    dict({[H.A21.level,H.A21.coltreeindex]}) = {H.A21.Y*v(H.A21.Ic(1)-H.Ic(1)+1:H.A21.Ic(2)-H.Ic(1)+1,:)};
    dict({[H.A12.level,H.A12.coltreeindex]}) = {H.A12.Y*v(H.A12.Ic(1)-H.Ic(1)+1:H.A12.Ic(2)-H.Ic(1)+1,:)};
end


end

function [y,dict] = descend(H,v,dict)

if H.isleaf == 0
    d1 = dict({[H.A21.level,H.A21.coltreeindex]});
    d2 = dict({[H.A12.level,H.A12.coltreeindex]});
    if H.isroot
        updateto2 = H.A21.Z*H.A21.lrcomponent*d1{1};
        updateto1 = H.A12.Z*H.A12.lrcomponent*d2{1};
        dict({[H.A21.level,H.A21.coltreeindex]}) = {updateto1};
        dict({[H.A12.level,H.A12.coltreeindex]}) = {updateto2};
    else
        currd = dict({[H.level,H.coltreeindex]});

        updateto1 = H.A12.Z*(currd{1}(1:size(H.A12.Z,2),:)+H.A12.lrcomponent*d2{1});
        updateto2 = H.A21.Z*(currd{1}(size(H.A12.Z,2)+1:end,:)+H.A21.lrcomponent*d1{1});

        dict({[H.A21.level,H.A21.coltreeindex]}) = {updateto1};
        dict({[H.A12.level,H.A12.coltreeindex]}) = {updateto2};
    end
    [y1,dict] = descend(H.A11,v(1:H.A11.sz(2),:),dict);
    [y2,dict] = descend(H.A22,v(H.A11.sz(2)+1:end,:),dict);
    y = [y1;y2];
else
    currd = dict({[H.level,H.rowtreeindex]});
    y = H.D*v + currd{1};
end
   
end