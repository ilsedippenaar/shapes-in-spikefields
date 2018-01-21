function saveToR(x, y, std_err, data_dir, name, varargin)
p = inputParser;
p.KeepUnmatched = true;
p.parse(varargin{:});
args = p.Unmatched;

r_data = [];
r_data.x = x;
r_data.y = y;
r_data.std_err = std_err;

fields = fieldnames(args);
for i=1:numel(fields)
  r_data.(fields{i}) = args.(fields{i});
end

save(fullfile(data_dir, sprintf('../R/jan_grant_app/%s.mat', name)), '-struct', 'r_data');
end