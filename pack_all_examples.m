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

%% Packing Example 4.1
disp('-->Packing Example 4.1')
cd('Chapter4_Example1')
myfile = fullfile(cpath,'publish','RADP_Ch4Ex1.zip');
zip(myfile, {'*.m'});
cd(cpath)

%% Packing Example 4.2
disp('-->Packing Example 4.2')
cd('Chapter4_Example2')
myfile = fullfile(cpath,'publish','RADP_Ch4Ex2.zip');
zip(myfile, {'*.m'});
cd(cpath)

%% Packing Example 4.3
disp('-->Packing Example 4.3')
cd('Chapter4_Example3')
myfile = fullfile(cpath,'publish','RADP_Ch4Ex3.zip');
zip(myfile, {'*.m'});
cd(cpath)

%% Packing Example 4.4
disp('-->Packing Example 4.4')
cd('Chapter4_Example4')
myfile = fullfile(cpath,'publish','RADP_Ch4Ex4.zip');
zip(myfile, {'*.m'});
cd(cpath)

%% Packing Example 5.1
disp('-->Packing Example 5.1')
cd('Chapter5_Example1')
myfile = fullfile(cpath,'publish','RADP_Ch5Ex1.zip');
zip(myfile, {'*.m'});
cd(cpath)

%% Packing Example 5.2
disp('-->Packing Example 5.2')
cd('Chapter5_Example2')
myfile = fullfile(cpath,'publish','RADP_Ch5Ex2.zip');
zip(myfile, {'*.m'});
cd(cpath)

%% Packing Example 6.1
disp('-->Packing Example 6.1')
cd('Chapter6_Example1')
myfile = fullfile(cpath,'publish','RADP_Ch6Ex1.zip');
zip(myfile, {'*.m'});
cd(cpath)

%% Packing Example 7.1
disp('-->Packing Example 7.1')
cd('Chapter7_Example1')
myfile = fullfile(cpath,'publish','RADP_Ch7Ex1.zip');
zip(myfile, {'*.m'});
cd(cpath)

%% Packing Example 7.2
disp('-->Packing Example 7.2')
cd('Chapter7_Example2')
myfile = fullfile(cpath,'publish','RADP_Ch7Ex2.zip');
zip(myfile, {'*.m'});
cd(cpath)

%% Packing Tools
disp('-->Packing Tools')
cd(cpath)
myfile = fullfile(cpath,'publish','tools.zip');
zip(myfile, 'tools');


