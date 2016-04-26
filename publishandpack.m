% Publish and pack all examples

%% Publish all examples
% Clean up before start
clear; bdclose all; close all; clc
% Clean up the PUBLISH folder
pubdir = [pwd '\publish\'];
if exist('publish', 'dir')
    rmdir(pubdir,'s');
end
mkdir(pubdir);

%% Publish Example 2.1
disp('-->Publishing Example 2.1')
cpath = pwd; % Save the root level path
options = struct('outputDir',[cpath '\publish\Ch2Ex1\']);
cd('Chapter2_Example1')
publish('Ch2Ex1_main.m', options);
movefile([cpath '\publish\Ch2Ex1\Ch2Ex1_main.html'], [cpath '\publish\Ch2Ex1\index.html']);
cd(cpath)

%% Publish Example 2.2
disp('-->Publishing Example 2.2')
cpath = pwd; % Save the root level path
options = struct('outputDir',[cpath '\publish\Ch2Ex2\']);
cd('Chapter2_Example2')
publish('Ch2Ex2_main.m', options);
movefile([cpath '\publish\Ch2Ex2\Ch2Ex2_main.html'], [cpath '\publish\Ch2Ex2\index.html']);
cd(cpath)
clear; close all;  % clean up

%% Publish Example 3.1
disp('-->Publishing Example 3.1')
cpath = pwd; % Save the root level path
options = struct('outputDir',[cpath '\publish\Ch3Ex1\']);
cd('Chapter3_Example1')
publish('Ch3Ex1_main.m', options);
movefile([cpath '\publish\Ch3Ex1\Ch3Ex1_main.html'], [cpath '\publish\Ch3Ex1\index.html']);
cd(cpath)
clear; close all;  % clean up

%% Chapter 4 
% Setup CVX
addpath('tools');
run .\tools\setuptools.m
%% Publish Example 4.1
disp('-->Publishing Example 4.1')
cpath = pwd; % Save the root level path
options = struct('outputDir',[cpath '\publish\Ch4Ex1\']);
cd('Chapter4_Example1')
publish('Ch4Ex1_main.m', options);
movefile([cpath '\publish\Ch4Ex1\Ch4Ex1_main.html'], [cpath '\publish\Ch4Ex1\index.html']);
cd(cpath)
clear; close all;  % clean up


%% Publish Example 4.2
disp('-->Publishing Example 4.2')
cpath = pwd; % Save the root level path
options = struct('outputDir',[cpath '\publish\Ch4Ex2\']);
cd('Chapter4_Example2')
publish('Ch4Ex2_main.m', options);
movefile([cpath '\publish\Ch4Ex2\Ch4Ex2_main.html'], [cpath '\publish\Ch4Ex2\index.html']);
cd(cpath)
clear; close all;  % clean up

%% Publish Example 4.3
disp('-->Publishing Example 4.3')
cpath = pwd; % Save the root level path
options = struct('outputDir',[cpath '\publish\Ch4Ex3\']);
cd('Chapter4_Example3')
publish('Ch4Ex3_main.m', options);
movefile([cpath '\publish\Ch4Ex3\Ch4Ex3_main.html'], [cpath '\publish\Ch4Ex3\index.html']);
cd(cpath)
clear; close all;  % clean up

%% Publish Example 4.4
disp('-->Publishing Example 4.4')
cpath = pwd; % Save the root level path
options = struct('outputDir',[cpath '\publish\Ch4Ex4\']);
cd('Chapter4_Example4')
publish('Ch4Ex4_main.m', options);
movefile([cpath '\publish\Ch4Ex4\Ch4Ex4_main.html'], [cpath '\publish\Ch4Ex4\index.html']);
cd(cpath)
clear; close all;  % clean up


%% Publish Example 5.1
disp('-->Publishing Example 5.1')
cpath = pwd; % Save the root level path
options = struct('outputDir',[cpath '\publish\Ch5Ex1\']);
cd('Chapter5_Example1')
publish('Ch5Ex1_main.m', options);
movefile([cpath '\publish\Ch5Ex1\Ch5Ex1_main.html'], [cpath '\publish\Ch5Ex1\index.html']);
cd(cpath)
clear; close all;  % clean up

%% Publish Example 5.2
disp('-->Publishing Example 5.2')
cpath = pwd; % Save the root level path
options = struct('outputDir',[cpath '\publish\Ch5Ex2\']);
cd('Chapter5_Example2')
publish('Ch5Ex2_main.m', options);
movefile([cpath '\publish\Ch5Ex2\Ch5Ex2_main.html'], [cpath '\publish\Ch5Ex2\index.html']);
cd(cpath)
clear; close all;  % clean up


%% Publish Example 6.1
disp('-->Publishing Example 6.1')
cpath = pwd; % Save the root level path
options = struct('outputDir',[cpath '\publish\Ch6Ex1\']);
cd('Chapter6_Example1')
publish('Ch6Ex1_main.m', options);
movefile([cpath '\publish\Ch6Ex1\Ch6Ex1_main.html'], [cpath '\publish\Ch6Ex1\index.html']);
cd(cpath)
clear; close all;  % clean up

%% Publish Example 7.1
disp('-->Publishing Example 7.1')
cpath = pwd; % Save the root level path
options = struct('outputDir',[cpath '\publish\Ch7Ex1\']);
cd('Chapter7_Example1')
publish('Ch7Ex1_main.m', options);
movefile([cpath '\publish\Ch7Ex1\Ch7Ex1_main.html'], [cpath '\publish\Ch7Ex1\index.html']);
cd(cpath)
clear; close all;  % clean up

%% Publish Example 7.2
disp('-->Publishing Example 7.2')
cpath = pwd; % Save the root level path
options = struct('outputDir',[cpath '\publish\Ch7Ex2\']);
cd('Chapter7_Example2')
publish('Ch7Ex2_main.m', options);
movefile([cpath '\publish\Ch7Ex2\Ch7Ex2_main.html'], [cpath '\publish\Ch7Ex2\index.html']);
cd(cpath)
clear; close all;  % clean up

%% Pack all examples
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


