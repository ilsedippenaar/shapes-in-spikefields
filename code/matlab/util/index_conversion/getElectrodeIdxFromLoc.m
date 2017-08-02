function idx = getElectrodeIdxFromLoc(loc)
total_idx = 10*(loc(1)-1) + loc(2);
if any(total_idx == [1 10 91 100])
  idx = [];
elseif total_idx < 10 % first row
  idx = total_idx - 1;
elseif total_idx > 90
  idx = total_idx - 3;
else
  idx = total_idx - 2;
end
end