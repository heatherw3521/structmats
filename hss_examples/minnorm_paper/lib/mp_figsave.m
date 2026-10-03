function mp_figsave(fig, figdir, name)
%MP_FIGSAVE  save a figure as PDF and PNG (works in MATLAB and Octave).
if ~exist(figdir, 'dir'), mkdir(figdir); end
set(fig, 'PaperPositionMode', 'auto');
print(fig, fullfile(figdir, [name '.pdf']), '-dpdf');
print(fig, fullfile(figdir, [name '.png']), '-dpng', '-r160');
close(fig);
end
