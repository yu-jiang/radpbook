%% Publish all examples
% Clean up before start
clear; bdclose all; close all; clc
% Save the root level path
cpath = pwd;
% Create options
options = struct('outputDir',[pwd '.\publish\']);
% Clean up the PUBLISH folder
delete([pwd '.\publish\*.*'])

%% Publish Example 2.1
cd('Chapter2_Example1')
publish('Ch2Ex1_main.m', options);
cd(cpath)