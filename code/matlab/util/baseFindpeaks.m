function [Ypk,Xpk,Wpk,Ppk] = baseFindpeaks(Yin,varargin)
% This is necessary since matlab for some reason 
% doesn't have any notion of scope
current_dir = cd(fullfile(matlabroot, 'toolbox/signal/signal'));
[Ypk, Xpk, Wpk, Ppk] = findpeaks(Yin, varargin{:});
cd(current_dir);
end