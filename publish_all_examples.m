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