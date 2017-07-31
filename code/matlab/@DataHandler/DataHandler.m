classdef DataHandler
  %DATAHANDLER The DATAHANDLER class provides an abstraction from the underlying trial and LFP data.
  properties
    num_trials
    num_lfp_electrodes
    num_units
    num_trial_sections = 6
    
    spike_sample_freq % per second
    lfp_sample_freq
    ts_per_ms % timestamps per millisecond
    electrode_mapping
    unit_snr
    date
    
    section_map = containers.Map(...
      {'pre_trial', 'start_to_fix', 'fix_to_noise', 'noise_to_shape', 'shape_to_stop', 'post_trial'}, ...
      num2cell(1:6));
    
    fixate
    noise
    shape
    saccade
    
    trials
    spikes
    lfps
  end
  methods
    function obj = DataHandler(trial_struct, sorted_trodes_list, lfps, electrode_lfp_mapping, varargin)
      %DATAHANDLER Default constructor for DataHandler
      % This constructor builds a small internal scheme for the trial 
      % structure and makes contiguous spike data.
      % See FROMFILE
      obj.num_trials = numel(trial_struct);
      obj.num_lfp_electrodes = size(lfps, 2);
      
      obj.lfps = int16(lfps);
      obj.electrode_mapping = cell(size(trial_struct(1).spikes, 2), 2);      
      obj.electrode_mapping(electrode_lfp_mapping,1) = num2cell(1:obj.num_lfp_electrodes);
      
      p = inputParser;
      p.addParameter('SNR_cutoff', 2.5);
      p.addParameter('lfp_sample_freq', 1e3);
      p.addParameter('spike_sample_freq', 3e4);
      p.addParameter('save_name', []);
      p.addParameter('date', []);
      p.addParameter('verbose', true);
      p.parse(varargin{:});
      args = p.Results;
      
      obj.spike_sample_freq = args.spike_sample_freq;
      obj.lfp_sample_freq = args.lfp_sample_freq;
      obj.ts_per_ms = obj.spike_sample_freq / 1e3;
      obj.date = args.date;
      
      valid_unit_selec = sorted_trodes_list(:,3) > args.SNR_cutoff;
      obj.unit_snr = sorted_trodes_list(valid_unit_selec, 3);
      valid_spike_electrodes = sorted_trodes_list(valid_unit_selec,1);
      obj.num_units = numel(valid_spike_electrodes);
      for i=1:obj.num_units
        elec_num = valid_spike_electrodes(i);
        obj.electrode_mapping{elec_num,2} = ...
          [obj.electrode_mapping{elec_num,2} i];
      end
      
      obj.trials = obj.makeTrialRepresentation(trial_struct);
      obj.spikes = obj.concatSpikes(trial_struct, valid_unit_selec);
      
      all_events = [obj.trials.sections];
      obj.fixate = int32([all_events{3:7:end}]);
      obj.noise = int32([all_events{4:7:end}]);
      obj.shape = int32([all_events{5:7:end}]);
      obj.saccade = int32([obj.trials.saccade]);
      
      if args.verbose
        printSummary(obj);
      end
      
      if ~isempty(args.save_name)
        fprintf('Saving DataHandler...\n');
        tmp = struct();
        tmp.dh = obj;
        save(args.save_name, '-struct', 'tmp'); % saves with name dh
      end
    end
  end
  methods
    out = select(obj, varargin)
    out = getLfpSlices(obj, spec, varargin)
    out = getDataSlices(obj, data, type, select_range, conditions, varargin)
    obj = clean(obj, varargin)
    function printSummary(obj)
      fprintf('--- Summary for DataHandler %s ---\n', obj.date);
      fprintf('Total time: %.3f minutes\n', size(obj.lfps,1) / obj.lfp_sample_freq / 60);
      fprintf('Number trials: %d\n', obj.num_trials);
      fprintf('Number of LFP electrodes: %d\n', obj.num_lfp_electrodes);
      fprintf('Number of units: %d\n', obj.num_units);
      fprintf('Number fixations: %d\n', numel(obj.fixate));
      fprintf('Number noise stimuli: %d\n', numel(obj.noise));
      fprintf('Number shape stimuli: %d\n', numel(obj.shape));
      fprintf('Number saccades: %d\n', numel(obj.saccade));
    end
  end
  methods (Static)
    function dataHandler = fromFile(trial_struct_filename, lfp_filename, varargin)
      fprintf('Loading trial file...\n');
      trial_file = load(trial_struct_filename);
      fprintf('Loading LFP file...\n');
      lfp_file = load(lfp_filename);
      fprintf('Building DataHandler...\n');
      p = inputParser;
      p.KeepUnmatched = true;
      p.parse(varargin{:});
      args = p.Unmatched;
      args.lfp_sample_freq = lfp_file.NS4.MetaTags.SamplingFreq * lfp_file.NS4.resampleFrac;
      args.spike_sample_freq = lfp_file.NS4.MetaTags.TimeRes;
      dataHandler = DataHandler(trial_file.trials, trial_file.file.sortedtrodes(:,1:3), ...
        lfp_file.NS4.rData, [lfp_file.NS4.ElectrodesInfo.ElectrodeID], args);
    end
    dataHandler = fromDates(dates, trial_dir, lfp_dir, data_handler_dir, varargin)
  end
  methods (Access = private)
    cont_spikes = concatSpikes(obj, trial_struct, valid_unit_selec)
    trials = makeTrialRepresentation(obj, trial_struct)
    out = selectSingle(obj, type, between, buffer, index, trial_num, trial_section)
  end  
end