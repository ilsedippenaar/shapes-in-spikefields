function loc = getLocFromElectrodeIdx(idx)
if idx < 1 || idx > 96
  loc = [];
elseif idx <= 8
  loc = [1 idx+1];
elseif idx >= 89
  loc = [10 idx-87];
else
  loc = [ceil((idx-8)/10)+1 mod(idx-9, 10)+1];
end
end