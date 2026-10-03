function [C, MK] = mp_style()
%MP_STYLE  Validated categorical palette (fixed order) + marker shapes used for
%   secondary encoding in every figure of the suite.
C = [0.165 0.471 0.839;   % #2a78d6 blue
     0.922 0.408 0.204;   % #eb6834 orange
     0.106 0.686 0.478;   % #1baf7a aqua
     0.929 0.631 0.000;   % #eda100 yellow
     0.910 0.482 0.643;   % #e87ba4 magenta
     0.000 0.514 0.000;   % #008300 green
     0.290 0.227 0.655;   % #4a3aa7 violet
     0.890 0.286 0.282];  % #e34948 red
MK = {'o','s','^','d','v','p','x','*'};
set(groot, 'defaultAxesFontSize', 10, 'defaultLineLineWidth', 1.5, ...
    'defaultAxesBox', 'off', 'defaultFigureColor', 'w');
end
