function trials_out = makeTrialRepresentation(obj, trial_struct, min_reaction_time)
trials_out(obj.num_trials).sections = [];
trials_out(obj.num_trials).saccade = [];
trials_out(obj.num_trials).result = [];
for i=1:obj.num_trials
  tr = trial_struct(i);
  if i == 1
    prev_trial_stop = 1;
  else
    prev_trial_stop = trial_struct(i-1).stopTime;
  end
  if i == obj.num_trials
    % calculate the very end of data that exists
    next_trial_start = size(obj.lfps, 1);
  else
    next_trial_start = trial_struct(i+1).startTime;
  end
  
  fixate = tr.fix.on;
  shape = tr.stim.shapeon;
  saccade = tr.fix.move;
  trials_out(i).sections = { ...
                            double(prev_trial_stop),...
                            double(tr.startTime), ...
                            double(fixate), ...
                            double(tr.stim.on), ...
                            double(shape), ...
                            double(tr.stopTime), ...
                            double(next_trial_start)};
  % saccades cannot be part of "sections" since a saccade can occur
  % anywhere in the trial. Trial section numbers would not be
  % well-defined if they were included.
  trials_out(i).saccade = double(saccade);
   
  if isempty(fixate) % no fixation occurred
    result = 'failed';
  else
    if isempty(shape) % no shape presented
      if isempty(saccade) % no saccade occurred
        result = 'true_negative';
      else
        result = 'false_positive';
      end
    else % shape was presented
      if isempty(saccade)
        result = 'false_negative';
      else
        if saccade > shape + min_reaction_time
          result = 'true_positive';
        else
          result = 'false_positive';
        end
      end
    end
  end
  trials_out(i).result = result;
end
end