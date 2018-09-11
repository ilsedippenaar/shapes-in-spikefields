function mt_params = getDefaultMtParams
T = 0.5;
W = 6;
mt_params = [];
mt_params.tapers = [T*W, 2*T*W-1];
mt_params.Fs = 1000;
end