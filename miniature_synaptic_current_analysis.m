%{
Author: S. Watanabe
Last updated; 7/21/2021

Description
This program is for the detection of miniature synaptic currents in whole-cell voltage clamp data of neurons. 

Inputs
data			[n,1] vector of the current data containing miniature events (n is the number of data points)
data_sd			Value of baseline SD from a data segment containing no events
sampling_freq	Value specifying the sampling frequency
polarity		String, either 'neg' (downward events) or 'pos' (upward events)
params			Structure with the following members	
  search_period				Duration of the segment of data to search a local peak
  baseline_period			Duration of the baseline before the local peak
  baseline_average_period	Duration of the baseline to calculate the baseline value
  decay_search_period		Duration to search the decay point after the peak
  min_separation			Minimum separation between the peaks
  sbtr_length				Window size for event subtraction
  peak_decay_fraction		Level of the trace to detect the decay point, as a fraction of the peak amplitude
  amp_thr_factor			Threshold of event amplitude, in units of baseline SD
  area_thr_factor			Threshold of event area (integral of the amplitude and time to the decay point), in units of baseline SD * millisecond

Outputs
output			Structure containing the following members
  peak_time		Vector of the peak time of detected events
  peak_amp		Vector of the amplitude of detected events
  baseline		Vector of the baseline for detected events
  onset_time	Vector of the onset time (peak time - baseline duration) for detected events
  decay_time 	Vector of the decay time (where the amplitude is peak_decay_fraction times peak_amp) of detected events
  peak_time_pt  Vector of the peak data point of detected events
  onset_time_pt Vector of the onset data point of detected events
  tau_decay_pt  Vector of the decay time constant (in number of data points) of detected events
  data_sbtr		[n,1] vector of the residual current data
  sbtr_wave		Cell array of subtracted data for detected events

Usage
1. Prepare an appropriately filtered data containing miniature events
2. Measure the SD from a segment containing no events
3. Run the program with appropriate parameters
4. Repeat running the program using the residual current data (data_sbtr) to detect double events

Example
output1 = miniature_synaptic_current_analysis(data, data_sd, sampling_freq, polarity, params);
output2 = miniature_synaptic_current_analysis(output1.data_sbtr, data_sd, sampling_freq, polarity, params);
output3 = miniature_synaptic_current_analysis(output2.data_sbtr, data_sd, sampling_freq, polarity, params);

Examples of parameter values
Optimal parameters should be determined by the user.
The following values may be appropriate for AMPA receptor mEPSCs and GABAA receptor mIPSCs in cortical pyramidal neurons.
  search_period = 0.08;
  baseline_period = 0.02; 
  baseline_average_period = 0.01;
  decay_search_period = 0.1;
  min_separation = 0.002;
  sbtr_length = 0.010; %for mEPSC
  sbtr_length = 0.080; %for mIPSC
  peak_decay_fraction = 0.37;
  amp_thr_factor = 2;
  area_thr_factor = 3; %for mEPSC
  area_thr_factor = 6; %for mIPSC
%}


function output = miniature_synaptic_current_analysis(data, data_sd, sampling_freq, polarity, params)

n_pt = length(data);
time = (0:n_pt-1)' / sampling_freq;


search_period = params.search_period;
baseline_period = params.baseline_period;
baseline_average_period = params.baseline_average_period;
decay_search_period = params.decay_search_period;
min_separation = params.min_separation;
sbtr_length = params.sbtr_length;
peak_decay_fraction = params.peak_decay_fraction;
amp_thr_factor = params.amp_thr_factor;
area_thr_factor = params.area_thr_factor;


%
amp_thr = data_sd * amp_thr_factor; 
area_thr = data_sd * area_thr_factor;

if strcmp(polarity, 'neg')
	amp_thr = -amp_thr;
	area_thr = -area_thr;
end

peak_search_period_pt = round(search_period * sampling_freq);
baseline_period_pt = round(baseline_period * sampling_freq);
baseline_average_period_pt = round(baseline_average_period * sampling_freq);
decay_search_period_pt = round(decay_search_period * sampling_freq);


n_segments = floor((time(end)-time(1))/ search_period);

n_peaks = 0;

peak_time = [];
peak_val = [];
baseline_val = []; 
onset_time  = []; 
decay_time  = []; 

 for i=1:n_segments
	
	% Take a segment of data of length peak_search_period and search a local peak
    seg_offset_pt = (i-1)*peak_search_period_pt;
	
    data_seg = data(seg_offset_pt + 1 : seg_offset_pt + peak_search_period_pt);
    if strcmp(polarity, 'neg')
        [peakval, peakpt] = min(data_seg);
    else
        [peakval, peakpt] = max(data_seg);
	end
	
	
	% Calculate the baseline for the data starting from baseline_period before the local peak and duration baseline_average_period
	if seg_offset_pt + peakpt - baseline_period_pt >0 && seg_offset_pt + peakpt + decay_search_period_pt < length(data)
		data_ext =  data( seg_offset_pt + peakpt - baseline_period_pt + 1 :  seg_offset_pt + peakpt + decay_search_period_pt);
		baseline = mean( data_ext(1 : baseline_average_period_pt) );
		
		% Calculate the local peak amplitude
        peakamp = peakval - baseline;

		% Check if the local peak amplitude exceeds amp_thr
		if ( strcmp(polarity, 'neg') && peakamp < amp_thr ) || ( strcmp(polarity, 'pos') && peakamp > amp_thr )
           
			% Detect the onset of the event as the point that reaches 0.5% of the peak amplitude
			if strcmp(polarity, 'neg')
			    onset_pt = find( data_ext(1:baseline_period_pt) > baseline + 0.005 * peakamp, 1, 'last');
			else
			    onset_pt = find( data_ext(1:baseline_period_pt) < baseline + 0.005 * peakamp, 1, 'last');
			end

			% Calculate the decay time point as the point that reaches peak_decay_fraction times the peak amplitude
			if strcmp(polarity, 'neg')
			    decay_pt = baseline_period_pt + find(data_ext(baseline_period_pt+1:end) > baseline + peak_decay_fraction * peakamp, 1, 'first');
			else
			    decay_pt = baseline_period_pt + find(data_ext(baseline_period_pt+1:end) < baseline + peak_decay_fraction * peakamp, 1, 'first');
			end
			
			if isempty(decay_pt)
				decay_pt = peakpt +1;
			end
			

		% Calculate the area under the curve (sum of the amplitude from the peak time point to the decay time point)
			
			if  decay_pt > decay_search_period_pt
				decay_pt = decay_search_period_pt;
			end
			
			if strcmp(polarity, 'neg')
				auc_pt = -sum(data_ext(baseline_period_pt : decay_pt) -baseline);
			else
				auc_pt =  sum(data_ext(baseline_period_pt : decay_pt) -baseline);
			end
			
			auc = auc_pt / sampling_freq *1000; 

		% Compare the area under the curve with area_thr
			if auc > abs(area_thr)
				n_peaks = n_peaks +1;
				peak_time(n_peaks) = (seg_offset_pt + peakpt) / sampling_freq;
				peak_val(n_peaks) = peakval;
				baseline_val(n_peaks) = baseline;
				onset_time(n_peaks) = (seg_offset_pt + peakpt - baseline_period_pt + onset_pt -1) / sampling_freq;
				decay_time(n_peaks) = (seg_offset_pt + peakpt - baseline_period_pt + decay_pt -1) / sampling_freq;
			end
		end
	end
 end

 
% Eliminate redundant peaks

peak_to_eliminate = zeros(1, n_peaks);

for i=1:n_peaks
	for j=i+1:n_peaks
		if abs(peak_time(j)-peak_time(i)) < min_separation
			peak_to_eliminate(j) = 1;
		end
	end
end
n_peaks_s = sum(~peak_to_eliminate);

peak_time_s = peak_time(~peak_to_eliminate);
peak_val_s = peak_val(~peak_to_eliminate);
baseline_val_s = baseline_val(~peak_to_eliminate);
onset_time_s = onset_time(~peak_to_eliminate);
decay_time_s = decay_time(~peak_to_eliminate);

peak_time_pt = round(peak_time_s * sampling_freq);
onset_time_pt = round(onset_time_s * sampling_freq);
tau_decay_pt = round((decay_time_s - peak_time_s) * sampling_freq);

peak_amp_s = peak_val_s - baseline_val_s;


% Subtract detected events from the data
sbtr_size = round(sbtr_length * sampling_freq);

data_sbtr = data;
sbtr_wave = cell(1, n_peaks_s);

for i=1:n_peaks_s

	sbtr_wave{i} = zeros(sbtr_size, 1);
	sbtr_size_pre = peak_time_pt(i) - onset_time_pt(i);
	sbtr_size_post = sbtr_size - sbtr_size_pre;
	sbtr_wave{i}(1 : sbtr_size_pre) =  data(onset_time_pt(i) : peak_time_pt(i)-1 ) - baseline_val_s(i);
	sbtr_wave{i}(sbtr_size_pre+1 : end) =  peak_amp_s(i) * exp(-(0 : sbtr_size_post-1)' / tau_decay_pt(i) );
	
	data_sbtr(onset_time_pt(i) :onset_time_pt(i) + sbtr_size - 1) = data(onset_time_pt(i) :onset_time_pt(i) + sbtr_size - 1) - sbtr_wave{i};
end


output.peak_time = peak_time_s;
output.peak_amp = peak_amp_s;
output.baseline = baseline_val_s;
output.onset_time = onset_time_s;
output.decay_time = decay_time_s;
output.peak_time_pt = peak_time_pt;
output.onset_time_pt = onset_time_pt;
output.tau_decay_pt = tau_decay_pt;

output.data_sbtr = data_sbtr;
output.sbtr_wave = sbtr_wave;




