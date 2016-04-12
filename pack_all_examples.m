% Create zip files for each example
%% Publish all examples
% Clean up before start
clear; bdclose all; close all; clc
% Clean up the PUBLISH folder
pubdir = [pwd '\publish\'];
if exist('publish', 'dir')
    delete([pubdir, '*.zip'])
else
    mkdir(pubdir);
end

cpath = pwd; % Save the root level path
%% Packing Example 2.1
disp('-->Packing Example 2.1')
cd('Chapter2_Example1')
myfile = fullfile(cpath,'publish','RADP_Ch2Ex1.zip');
zip(myfile, {'*.m', '*.slx'});
cd(cpath)

%% Packing Example 2.2
disp('-->Packing Example 2.2')
cd('Chapter2_Example2')
myfile = fullfile(cpath,'publish','RADP_Ch2Ex2.zip');
zip(myfile, {'*.m'});
cd(cpath)

%% Packing Example 3.1
disp('-->Packing Example 3.1')
cd('Chapter3_Example1')
myfile = fullfile(cpath,'publish','RADP_Ch3Ex1.zip');
zip(myfile, {'*.m'});
cd(cpath)


%% Packing Example 4.4
disp('-->Packing Example 4.4')
cd('Chapter4_Example4')
myfile = fullfile(cpath,'publish','RADP_Ch3Ex1.zip');
zip(myfile, {'*.m'});
cd(cpath)


%% Packing Tools
disp('-->Packing Tools')
cd(cpath)
myfile = fullfile(cpath,'publish','tools.zip');
zip(myfile, 'tools');


