function psd = getPsd(tbl, analysis_idxs, mt_params)
pows = rowfun(@(x) {mtspectrumc(single(x{1}(analysis_idxs,:)), mt_params)}, tbl, 'InputVariables', {'lfps'});
pows = rowfun(@(x) {x{1}/max(x{1})}, pows); % TODO: try w/o normalization
psd = [pows.Var1{:}]';
end