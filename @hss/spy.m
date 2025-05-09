function spy(H)
% visualization of rank structure of matrix

m= H.size(1);
n= H.size(2);

cutoff = floor(log2(m/H.blocksize)/2)+2;

set(gca, 'YDir', 'reverse');

hold on

draw_spy(H,[1,1,1],1,cutoff)
xlim([H.Ir(1) H.Ir(2)])
ylim([H.Ic(1) H.Ic(2)])
%axis equal
hold off


end




function draw_spy(H,leafcolor,yestext,cutoff)
if H.isdiag
    leafcolor = leafcolor/2;
end
diagcolor = [0.1, 0.3, 0.95];
% purple
%textcolor = [0.4940 0.1840 0.5560];
% dull yellow
textcolor = [0.9290 0.6940 0.1250];

top = H.Ir(1)-0.5;
bottom = H.Ir(2)+0.5;
left = H.Ic(1)-0.5;
right = H.Ic(2)+0.5;
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

    if H.level == cutoff && ~H.isdiag
        rowdims = size(H.Z);
        coldims = size(H.Y);
        text(H.Ir(1)+H.size(1)/2,H.Ic(2)-H.size(2)/2,sprintf('%d',min(rowdims(2),coldims(1))),'Color',textcolor,'HorizontalAlignment','center')
    elseif H.level == cutoff+1 && ~H.isdiag && yestext
        rowdims = size(H.Z);
        coldims = size(H.Y);
        text(H.Ir(1)+H.size(1)/2,H.Ic(2)-H.size(2)/2,sprintf('%d',min(rowdims(2),coldims(1))),'Color',textcolor,'HorizontalAlignment','center')
    end

else
    fill(xb,yb,leafcolor)
    rowdims = size(H.Z);
    coldims = size(H.Y);
    %text(H.Ir(1)+H.size(1)/2,H.Ic(2) + H.size(2)/2,sprintf('%d',min(rowdims(2),coldims(1))),'Color','red','FontSize',14)
    if yestext
        text(H.Ir(1)+H.size(1)/2,H.Ic(2)-H.size(2)/2,sprintf('%d',min(rowdims(2),coldims(1))),'Color',textcolor,'HorizontalAlignment','center')
    end
end

end