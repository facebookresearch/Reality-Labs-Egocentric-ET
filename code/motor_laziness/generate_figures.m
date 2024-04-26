% Name: generate_figures.m
%
% Copyright (c) Meta Platforms, Inc. and affiliates.
%
% Author: Charlie S. Burlingham
%
% Usage: First add /matlab_deps/ and /data/ and /motor_laziness/ to path.
%        Run whole script from the /data/ directory containing data.mat.
%
% Date last updated: Feb 7, 2024
%
% dependencies: hline/vline
%       pTestModule
%       fdr_bh
%       dsErrorSurface
%       histogram2Polar
% 
%

d = getData(pwd);


% Body Text Figures
fig1(d);
[e] = fig2(d);
fig3(d,e);
fig4(d);

% Supplementary Figures
figS1(d)
figS2(d)
figS3(d,e)
figS4and5(d,e)
figS6and7(d)
figS8(d,e)
figS9(d)
figS10(d)
figS11(d)

