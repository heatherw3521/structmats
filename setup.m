%Instructions for seting up HM-toolbox and creating add-ons for it. 
% This is required to set up the structmat solvers. It should be run only 
% once. 
%
% Once you have installed the structmat directory, do the following: 
%
% First, install the HM-toolbox via git: 
% (1) Use a terminal to navigate to the directory you wish to install the toolbox. 
% (2) Enter the following command: git clone https://github.com/numpi/hm-toolbox.git
%
% Next, run the setup: 
% (3) In MATLAB, add the paths hm-toolbox/ and structmats/ to the working
% directory (with subfolders)
% (4) Navigate to a parent folder that contains both the structmats directory and the 
% hm-toolbox. From there, run this file. It will add to the toolbox some additional functions used by structmat. 
%
% YOU MAY NEED TO ALTER/UPDATE the filepaths in the 'copyfile' commands so that it finds the
% correct location for the hm-toolbox on your machine. 
%
% WARNING: This file invokes a 'clear all, clear classes, clear functions' command.
%
% Once this setup is complete, it will remove toolbox_addons from the
% structmat file. This directory can then be safely deleted to avoid confusion. 
%
% Check out structmats/Examples_Solvers.m for examples on how to get
% started!


%% set up fADI and URV functions: 
copyfile('structmats\toolbox_addons\urv.m', 'hm-toolbox\@hss\urv.m')
copyfile('structmats\toolbox_addons\urv_solve.m', 'hm-toolbox\urv_solve.m')
copyfile('structmats\toolbox_addons\fADI_col.m', 'hm-toolbox\@hss\private\fADI_col.m')
copyfile('structmats\toolbox_addons\fADI_row.m', 'hm-toolbox\@hss\private\fADI_row.m')
copyfile('structmats\toolbox_addons\getshifts_adi.m', 'hm-toolbox\@hss\private\getshifts_adi.m')
copyfile('structmats\toolbox_addons\hss_build_hss_tree_nudft.m', 'hm-toolbox\@hss\private\hss_build_hss_tree_nudft.m')
copyfile('structmats\toolbox_addons\hss_cauchytoeplitz.m', 'hm-toolbox\@hss\private\hss_cauchytoeplitz.m')
copyfile('structmats\toolbox_addons\hss_cauchytoeplitz2.m', 'hm-toolbox\@hss\private\hss_cauchytoeplitz2.m')
copyfile('structmats\toolbox_addons\hss_nudftv.m', 'hm-toolbox\@hss\private\hss_nudftv.m')
copyfile('structmats\toolbox_addons\hss_urv_fact.m', 'hm-toolbox\@hss\private\hss_urv_fact.m')
copyfile('structmats\toolbox_addons\vbuildcauchydiags.m', 'hm-toolbox\@hss\private\vbuildcauchydiags.m')
copyfile('structmats\toolbox_addons\vfADI_col.m', 'hm-toolbox\@hss\private\vfADI_col.m')
copyfile('structmats\toolbox_addons\vfADI_row.m', 'hm-toolbox\@hss\private\vfADI_row.m')
copyfile('structmats\toolbox_addons\hss_urv_fact_solve.m', 'hm-toolbox\private\hss_urv_fact_solve.m')

%%
%%modify the hss object class file. WARNING: This will delete and replace
% the object class file for hm-toolbox\@hss\. It keeps all class attributes
% from the hm-toolbox as of 8/1/2023, and then adds some additional
% attributes used by structmat. If your hm-toolbox has updates from after
% Aug 1 2023, they may be lost. 
delete hm-toolbox\@hss\hss.m
copyfile('structmats\toolbox_addons\hss.m', 'hm-toolbox\@hss\hss.m')
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
    rmpath('structmats\toolbox_addons\')
end










