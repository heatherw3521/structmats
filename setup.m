%Instructions for seting up HM-toolbox and creating add-ons for it. 
% This is required to run the structmat solvers. 
%
% First, install the HM-toolbox via git: 
% (1) Use a terminal to navigate to the directory you wish to install the toolbox. 
% (2) Enter the following command: git clone https://github.com/numpi/hm-toolbox.git
%
% Once the toolbox is added, run this file. It will add the
% appropriate path and also add to the toolbox some additional
% functions used by structmat. 
%
% WARNING: This file invokes a 'clear all, clear classes, clear functions' command. 

%%
addpath(hm-toolboxroot)

%% set up fADI and URV functions: 
movefile('toolbox_addons\urv.m', 'hm-toolbox\@hss\urv.m')
movefile('toolbox_addons\urv_solve.m', 'hm-toolbox\urv_solve.m')
movefile('toolbox_addons\fADI_col.m', 'hm-toolbox\@hss\private\fADI_col.m')
movefile('toolbox_addons\fADI_row.m', 'hm-toolbox\@hss\private\fADI_row.m')
movefile('toolbox_addons\getshifts_adi.m', 'hm-toolbox\@hss\private\getshifts_adi.m')
movefile('toolbox_addons\hss_build_hss_tree_nudft.m', 'hm-toolbox\@hss\private\hsss_build_hss_tree_nudft.m')
movefile('toolbox_addons\hss_cauchytoeplitz.m', 'hm-toolbox\@hss\private\hss_cauchytoeplitz.m')
movefile('toolbox_addons\hss_cauchytoeplitz2.m', 'hm-toolbox\@hss\private\hss_cauchytoeplitz2.m')
movefile('toolbox_addons\hss_nudftv.m', 'hm-toolbox\@hss\private\hss_nudftv.m')
movefile('toolbox_addons\hss_urv_fact.m', 'hm-toolbox\@hss\private\hss_urv_fact.m')
movefile('toolbox_addons\vbuildcauchydiags.m', 'hm-toolbox\@hss\private\vbuildcauchydiags.m')
movefile('toolbox_addons\vfADI_col.m', 'hm-toolbox\@hss\private\vfADI_col.m')
movefile('toolbox_addons\vfADI_row.m', 'hm-toolbox\@hss\private\vfADI_row.m')
%%
%%modify the hss object class file. WARNING: This will delete and replace
% the object class file for hm-toolbox\@hss\. It keeps all class attributes
% from the hm-toolbox as of 8/1/2023, and then adds some additional
% attributes used by structmat. If your hm-toolbox has updates from after
% Aug 1 2023, they may be lost. 
delete hm-toolbox\@hss\hss.m
movefile('toolbox_addons\hss.m', 'hm-toolbox\@hss\hss.m')
%clear classes/functions to reset with new class str
clear classes
clear functions
clear all










