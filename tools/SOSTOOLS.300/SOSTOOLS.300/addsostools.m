% Add SOS Tools to path
%

cm = computer;
if cm(1) == 'M' ||  cm(1)=='G'
    addpath(pwd);
    addpath([pwd '/multipoly']);
    addpath([pwd '/internal']);
    addpath([pwd '/demos']);
    addpath([pwd '/custom']);
elseif cm(1) == 'P'
    addpath(pwd);
    addpath([pwd '\multipoly']);
    addpath([pwd '\internal']);
    addpath([pwd '\demos']);
    addpath([pwd '\custom']);
end

% Alternative single-line syntax to add all subfolders
% addpath(genpath(pwd));    

