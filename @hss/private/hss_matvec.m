function y = hss_matvec(H,v)
% hss times vector multiplication
% H: mxn hss matrix
% v: nx1 vector
vsize = size(v);
% y = zeros(H.size(1),vsize(2));
if vsize(2) >1
    y = [hss_matvec(H,v(:,1)), hss_matvec(H,v(:,2:end))];
    return
end
if H.isleaf == false
    dict = dictionary();
    % starting from the bottom build dict up
    dict = ascend(H,v,dict);

    % using dict descend 
    y = descend(H,v,dict);
    
elseif H.isdiag %leaf diag
    y = H.D * v;
else %leaf off-diag
    y = H.Z*H.lrcomponent*H.Y*v;
end


end

function dict = ascend(H,v,dict)


if H.A11.isleaf == false
    %dict = ascend(H.A11,v(H.A11.Ic(1)-H.Ic(1)+1:H.A11.Ic(2)-H.Ic(1)+1,1),dict);
    %dict = ascend(H.A22,v(H.A22.Ic(1)-H.Ic(1)+1:H.A22.Ic(2)-H.Ic(1)+1,1),dict);
    dict = ascend(H.A11,v(1:H.A11.size(2),1),dict);
    dict = ascend(H.A22,v(H.A11.size(2)+1:end,1),dict);
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
    dict({[H.A21.level,H.A21.coltreeindex]}) = {H.A21.Y*v(H.A21.Ic(1)-H.Ic(1)+1:H.A21.Ic(2)-H.Ic(1)+1,1)};
    dict({[H.A12.level,H.A12.coltreeindex]}) = {H.A12.Y*v(H.A12.Ic(1)-H.Ic(1)+1:H.A12.Ic(2)-H.Ic(1)+1,1)};
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

        % updateto2 = H.A21.Z*(currd{1}(1:size(H.A21.lrcomponent,1)) + H.A21.lrcomponent*d1{1});
        % updateto1 = H.A12.Z*(currd{1}(size(H.A12.lrcomponent,1)+1:end) + H.A12.lrcomponent*d2{1});
        
        % update = blkdiag(H.A21.Z,H.A12.Z) *(currd{1} + [H.A21.lrcomponent*d1{1};H.A12.lrcomponent*d2{1}]);
        % updateto2 = update(1:size(H.A21.Z,1));
        % updateto1 = update(size(H.A21.Z,1)+1:end);

        % updateto2 = H.A21.Z*(currd{1}(end/2+1:end) + H.A21.lrcomponent*d1{1});
        % updateto1 = H.A12.Z*(currd{1}(1:end/2) + H.A12.lrcomponent*d2{1});

        updateto1 = H.A12.Z*(currd{1}(1:size(H.A12.Z,2))+H.A12.lrcomponent*d2{1});
        updateto2 = H.A21.Z*(currd{1}(size(H.A12.Z,2)+1:end)+H.A21.lrcomponent*d1{1});

        dict({[H.A21.level,H.A21.coltreeindex]}) = {updateto1};
        dict({[H.A12.level,H.A12.coltreeindex]}) = {updateto2};
    end
    [y1,dict] = descend(H.A11,v(1:H.A11.size(2)),dict);
    [y2,dict] = descend(H.A22,v(H.A11.size(2)+1:end),dict);
    y = [y1;y2];
else
    currd = dict({[H.level,H.rowtreeindex]});
    y = H.D*v + currd{1};
end
   
end