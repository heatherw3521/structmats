function spy(H)
%SPY  Picture of the HSS block structure: diagonal leaves shaded, each
%   off-diagonal block labeled with its rank.

m= H.sz(1);
n= H.sz(2);

cutoff = H.levelcount;
set(gca, 'YDir', 'reverse');

hold on

draw_spy(H,[1,1,1],1,cutoff)
xlim([H.Ic(1) H.Ic(2)])
ylim([H.Ir(1) H.Ir(2)])
hold off


end




function draw_spy(H,leafcolor,yestext,cutoff)
if H.isdiag
    leafcolor = leafcolor/2;
end
diagcolor = [0.1, 0.3, 0.95];
textcolor = [0.9290 0.6940 0.1250];

top = H.Ic(1)-0.5;
bottom = H.Ic(2)+0.5;
left = H.Ir(1)-0.5;
right = H.Ir(2)+0.5;
xb = [top,bottom,bottom,top];
yb = [left,left,right,right];

if H.isdiag && H.isleaf
    
    fill(xb,yb,diagcolor)
elseif H.isdiag
    if H.level == cutoff && ~H.isdiag
        yestext = 0;
    elseif H.level == cutoff+1
        yestext = 0;
    end
    draw_spy(H.A11,leafcolor,yestext,cutoff)
    draw_spy(H.A12,leafcolor,yestext,cutoff)
    draw_spy(H.A21,leafcolor,yestext,cutoff)
    draw_spy(H.A22,leafcolor,yestext,cutoff)

    % rank labels in the off-diagonal blocks
    if H.level == cutoff && ~H.isdiag
        rowdims = size(H.Z);
        coldims = size(H.Y);
        text(H.Ic(1)+H.sz(2)/2,H.Ir(2)-H.sz(1)/2,sprintf('%d',min(rowdims(2),coldims(1))),'Color',textcolor,'HorizontalAlignment','center')
    elseif H.level == cutoff+1 && ~H.isdiag && yestext
        rowdims = size(H.Z);
        coldims = size(H.Y);
        text(H.Ic(1)+H.sz(2)/2,H.Ir(2)-H.sz(1)/2,sprintf('%d',min(rowdims(2),coldims(1))),'Color',textcolor,'HorizontalAlignment','center')
    end

else
    fill(xb,yb,leafcolor)
    rowdims = size(H.Z);
    coldims = size(H.Y);
    if yestext
        text(H.Ic(1)+H.sz(2)/2,H.Ir(2)-H.sz(1)/2,sprintf('%d',min(rowdims(2),coldims(1))),'Color',textcolor,'HorizontalAlignment','center')
    end
end

end