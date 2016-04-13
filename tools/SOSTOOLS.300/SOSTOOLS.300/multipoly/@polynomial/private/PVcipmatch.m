function [idx,fstr] = PVcipmatch(clist,str)
% function [Idx,Fstr] = PVcipmatch(Clist,Str)
%
% DESCRIPTION
%   Case insensitive partial match.
%
% INPUTS
%   Clist:
%   Str:
%
% OUTPUTS
%   Idx:
%   Fstr:
%
% SYNTAX
%   [Idx,Fstr] = PVcipmatch(Clist,Str)


% 6/7/2002: Written by Packard

lstr = length(str);
tf = find(strncmpi(clist,str,lstr));
if length(tf)==1
    idx = tf;
    fstr = clist{idx};
elseif isempty(tf)
    idx = [];
    fstr = [];
else
    cnt = 0;
    for i=1:length(tf)
        if length(clist{tf(i)})==lstr
            tfuse = tf(i);
            cnt = cnt + 1;
        end
    end
    if cnt==1
        idx = tfuse;
        fstr = clist{idx};
    else
        tf = find(strncmp(clist,str,lstr));
        if length(tf)==1
            idx = tf;
            fstr = clist{idx};
        elseif isempty(tf)
            idx = [];
            fstr = [];
        else
            cnt = 0;
            for i=1:length(tf)
                if length(clist{tf(i)})==lstr
                    tfuse = tf(i);
                    cnt = cnt + 1;
                end
            end
            if cnt==1
                idx = tfuse;
                fstr = clist{idx};
            else
                error(['Ambiguous property: ''' str '''.']);
            end
        end
    end
end
