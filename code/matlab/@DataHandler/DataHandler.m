classdef DataHandler
  %DATAHANDLER The DATAHANDLER class provides an abstraction from the underlying trial, LFP, and spike data.
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
    number_on_date = 1 % given certain selection parameters to avoid noise, there may be multiple data handlers for one day
    
    section_map = containers.Map(...
      {'pre_trial', 'start_to_fix', 'fix_to_noise', 'noise_to_shape', 'shape_to_stop', 'post_trial', ...
      'shape_to_saccade', 'saccade_to_stop'}, ...
      [num2cell(1:6), -1, -2]);
    
    fixate
    noise
    shape
    saccade
    
    trials
    spikes
    lfps
  end
  methods
    function obj = DataHandler(varargin)
      %DATAHANDLER Default constructor for DataHandler
      % This constructor builds a small internal scheme for the trial 
      % structure and makes contiguous spike data.
      % See FROMFILE
      
      p = inputParser;
      % main parameters
      p.addParameter('dh_struct', []);
      p.addParameter('trial_struct', []);
      p.addParameter('sorted_trodes_list', []);
      p.addParameter('lfps', []);
      p.addParameter('electrode_lfp_mapping', []);
      
      % tweaking parameters
      p.addParameter('SNR_cutoff', 2.5);
      p.addParameter('lfp_sample_freq', 1e3);
      p.addParameter('spike_sample_freq', 3e4);
      p.addParameter('save_name', []);
      p.addParameter('date', []);
      p.addParameter('number_on_date', 1);
      p.addParameter('clean', false);
      p.addParameter('verbose', true);
      p.parse(varargin{:});
      args = p.Results;
      
      if ~isempty(args.dh_struct)
        obj = obj.initializeWithStruct(args.dh_struct);
      else
        obj.num_trials = numel(args.trial_struct);
        obj.num_lfp_electrodes = size(args.lfps, 2);

        obj.lfps = int16(args.lfps);
        obj.electrode_mapping = cell(size(args.trial_struct(1).spikes, 2), 2);      
        obj.electrode_mapping(args.electrode_lfp_mapping,1) = num2cell(1:obj.num_lfp_electrodes);

        obj.spike_sample_freq = args.spike_sample_freq;
        obj.lfp_sample_freq = args.lfp_sample_freq;
        obj.ts_per_ms = obj.spike_sample_freq / 1e3;
        obj.date = args.date;

        valid_unit_selec = args.sorted_trodes_list(:,3) > args.SNR_cutoff;
        obj.unit_snr = args.sorted_trodes_list(valid_unit_selec, 3);
        valid_spike_electrodes = args.sorted_trodes_list(valid_unit_selec,1);
        obj.num_units = numel(valid_spike_electrodes);
        for i=1:obj.num_units
          elec_num = valid_spike_electrodes(i);
          obj.electrode_mapping{elec_num,2} = ...
            [obj.electrode_mapping{elec_num,2} i];
        end

        obj.trials = obj.makeTrialRepresentation(args.trial_struct);
        obj.spikes = obj.concatSpikes(args.trial_struct, valid_unit_selec);

        all_events = [obj.trials.sections];
        obj.fixate = int32([all_events{3:obj.num_trial_sections+1:end}]);
        obj.noise = int32([all_events{4:obj.num_trial_sections+1:end}]);
        obj.shape = int32([all_events{5:obj.num_trial_sections+1:end}]);
        obj.saccade = int32([obj.trials.saccade]);
      end
      
      if args.clean
        obj = obj.clean();
      end
      
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
    out = getDataSlices(obj, data, type, select_range, conditions, varargin)
    obj = clean(obj, varargin)
    dh_structs = split(obj, betweens)
    function printSummary(obj, fd)
      if nargin == 1
        fd = 1;
      end
      fprintf(fd,'--- Summary for DataHandler %s (%d) ---\n', obj.date, obj.number_on_date);
      fprintf(fd,'Total time: %.3f minutes\n', size(obj.lfps,1) / obj.lfp_sample_freq / 60);
      fprintf(fd,'Number trials: %d\n', obj.num_trials);
      fprintf(fd,'Number of LFP electrodes: %d\n', obj.num_lfp_electrodes);
      fprintf(fd,'Number of units: %d\n', obj.num_units);
      fprintf(fd,'Number fixations: %d\n', numel(obj.fixate));
      fprintf(fd,'Number noise stimuli: %d\n', numel(obj.noise));
      fprintf(fd,'Number shape stimuli: %d\n', numel(obj.shape));
      fprintf(fd,'Number saccades: %d\n', numel(obj.saccade));
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
      args.trial_struct = trial_file.trials;
      args.sorted_trodes_list = trial_file.file.sortedtrodes(:,1:3);
      args.lfps = lfp_file.NS4.rData;
      args.electrode_lfp_mapping = [lfp_file.NS4.ElectrodesInfo.ElectrodeID];
      dataHandler = DataHandler(args);
    end
    dataHandler = fromDates(dates, trial_dir, lfp_dir, data_handler_dir, varargin)
  end
  methods (Access = private)
    cont_spikes = concatSpikes(obj, trial_struct, valid_unit_selec)
    trials = makeTrialRepresentation(obj, trial_struct)
    out = selectSingle(obj, type, between, buffer, index, trial_num, trial_section)
    function obj = initializeWithStruct(obj, dh_struct)
      names = fieldnames(dh_struct);
      for i=1:numel(names)
        obj.(names{i}) = dh_struct.(names{i});
      end
    end
  end  
end