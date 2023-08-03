%Instructions for seting up HM-toolbox and creating add-ons for it. 
% This is required to run the structmat solvers. 
%
% First, install the HM-toolbox via git: 
% (1) Use a terminal to navigate to the directory you wish to install the toolbox. 
% (2) Enter the following command: git clone https://github.com/numpi/hm-toolbox.git
%
% Once the toolbox is added, navigate to the structmats folder and run this file. It will add the
% appropriate path and also add to the toolbox some additional
% functions used by structmat. 
%
% YOU MAY NEED TO ALTER/UPDATE the filepaths in the 'copyfile' commands so that it finds the
% correct location for the hm-toolbox on your machine. 
%
% WARNING: This file invokes a 'clear all, clear classes, clear functions' command. 

%%
addpath('hm-toolbox')
addpath('structmats')

%% set up fADI and URV functions: 
copyfile('toolbox_addons\urv.m', 'hm-toolbox\@hss\urv.m')
copyfile('toolbox_addons\urv_solve.m', 'hm-toolbox\urv_solve.m')
copyfile('toolbox_addons\fADI_col.m', 'hm-toolbox\@hss\private\fADI_col.m')
copyfile('toolbox_addons\fADI_row.m', 'hm-toolbox\@hss\private\fADI_row.m')
copyfile('toolbox_addons\getshifts_adi.m', 'hm-toolbox\@hss\private\getshifts_adi.m')
copyfile('toolbox_addons\hss_build_hss_tree_nudft.m', 'hm-toolbox\@hss\private\hss_build_hss_tree_nudft.m')
copyfile('toolbox_addons\hss_cauchytoeplitz.m', 'hm-toolbox\@hss\private\hss_cauchytoeplitz.m')
copyfile('toolbox_addons\hss_cauchytoeplitz2.m', 'hm-toolbox\@hss\private\hss_cauchytoeplitz2.m')
copyfile('toolbox_addons\hss_nudftv.m', 'hm-toolbox\@hss\private\hss_nudftv.m')
copyfile('toolbox_addons\hss_urv_fact.m', 'hm-toolbox\@hss\private\hss_urv_fact.m')
copyfile('toolbox_addons\vbuildcauchydiags.m', 'hm-toolbox\@hss\private\vbuildcauchydiags.m')
copyfile('toolbox_addons\vfADI_col.m', 'hm-toolbox\@hss\private\vfADI_col.m')
copyfile('toolbox_addons\vfADI_row.m', 'hm-toolbox\@hss\private\vfADI_row.m')
rmpath('toolbox_addons')
%%
%%modify the hss object class file. WARNING: This will delete and replace
% the object class file for hm-toolbox\@hss\. It keeps all class attributes
% from the hm-toolbox as of 8/1/2023, and then adds some additional
% attributes used by structmat. If your hm-toolbox has updates from after
% Aug 1 2023, they may be lost. 
delete hm-toolbox\@hss\hss.m
movefile('structmats\toolbox_addons\hss.m', 'hm-toolbox\@hss\hss.m')
%clear classes/functions to reset with new class str
clear classes
clear functions
clear all

%% does it work? (more thorough tests to come)

nodes = rand(2050,1); modes = 1000; rhs = rand(2050,1);
fail = 0; 
try
    structsolv_nudft2(nodes, modes, rhs, 'tol', 1e-4);
catch 
    fail = 1; 
    error('Installation failed.')
end
if ~fail
    fprintf('Installation appears to be successful. Fair warning, this package is in the early days of construction!')
end










