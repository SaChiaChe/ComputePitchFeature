#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <vector>

using namespace std;

#define MAXDATALEN 100000
#define MAXFEATLEN 10000

// Kaldi's floor behaves weirdly: floor(7) = 6
// Thus we have to implement our own floor
int myfloor(float x){
	int x_ = floor(x);
	if (x_ - x == 0.0)
		x_ -= 1;
	return x_;
}

float FilterFunc(double t, int C, int w){
	float Window, Filter;
	if (fabs(t) < w / (2.0 * C))
		Window = 0.5 * (1 + cos(M_PI * 2 * C / w * t));
	else
		Window = 0.0;  // outside support of window function
	if (t != 0.0)
		Filter = sin(M_PI * 2 * C * t) / (M_PI * t);
	else
		Filter = 2 * C;  // limit of the function at t = 0
	return Filter * Window;
}

float InnerProduct(short *Data, float *Weight, int Len){
	float Sum = 0.0;
	for(int i = 0; i < Len; i++){
		Sum += Data[i] * Weight[i];
	}
	return Sum;
}

double InnerProduct(float *V1, float *V2, int Len){
	double Sum = 0.0;
	for(int i = 0; i < Len; i++)
		Sum += static_cast<double>(V1[i]) * static_cast<double>(V2[i]);
	return Sum;
}

double Sum_(float *V, int Len){
	double Sum = 0.0;
	for(int i = 0; i < Len; i++)
		Sum += static_cast<double>(V[i]);
	return Sum;
}

float Min(float *V, int Len){
	float minimum = V[0];
	for(int i = 1; i < Len; i++){
		if(V[i] < minimum){
			minimum = V[i];
		}
	}
	return minimum;
}

int MinIndex(float *V, int Len){
	int minID = 0;
	for(int i = 1; i < Len; i++){
		if(V[i] < V[minID]){
			minID= i;
		}
	}
	return minID;
}

void Swap(float *V1, float *V2, int Len){
	for(int i = 0; i < Len; i++){
		float temp = V1[i];
		V1[i] = V2[i];
		V2[i] = temp;
	}
}

int Resample(int sample_frequency, int resample_frequency, int C, int w, short *Data, int DataLength, float *Output){
	double window_width = w / (2.0 * C);
	int minIndex = ceil(-window_width * sample_frequency);
	int maxIndex = myfloor(window_width * sample_frequency);
	int window_len = maxIndex - minIndex + 1;

	float weights[window_len];
	for(int i = minIndex; i <= maxIndex; i++){
		double Delta_t = i / static_cast<double>(sample_frequency);
		weights[i - minIndex] = FilterFunc(Delta_t, C, w) / sample_frequency;
	}

	int sample_shift = sample_frequency / resample_frequency;
	int OutputIndex = 0;
	for(int Offset = minIndex; Offset + window_len <= DataLength + maxIndex; Offset += sample_shift, OutputIndex++)
		Output[OutputIndex] = InnerProduct(Data + Offset, weights, window_len);

	return OutputIndex;
}

void Resample(float **input, float **output, int NumofRows, int DimIn, int DimOut, int *FirstIndex, int *WeightDim, float **Weights){
	// Using double pointer is kinda unsafe
	// Could consider modify this to vector
	// input[NumofRows][DimIn]
	// output[NumofRows][DimOut]
	for(int i = 0; i < DimOut; i++){
		float input_part[NumofRows][WeightDim[i]];
		for(int j = 0; j < NumofRows; j++){
			output[j][i] = InnerProduct(input[j] + FirstIndex[i], Weights[i], WeightDim[i]);
		}
	}
}

void SetFirstIndex(int *FirstIndex, float **Weights, double *sample_points, int num_samples_in, int num_samples_out, int upsample_filter_width, int filter_cutoff, int sample_rate_in, int *WeightDim){
	float filter_width = upsample_filter_width / (2.0 * filter_cutoff);
	for(int i = 0; i < num_samples_out; i++){
		float t = sample_points[i];
		float t_min = t - filter_width;
		float t_max = t + filter_width;
		int index_min = ceil(sample_rate_in * t_min);
		int index_max = myfloor(sample_rate_in * t_max);
		if(index_min < 0)
			index_min = 0;
		if(index_max >= num_samples_in)
			index_max = num_samples_in - 1;
		FirstIndex[i] = index_min;
		int weightdim = index_max - index_min + 1;
		Weights[i] = new float[weightdim];
		WeightDim[i] = weightdim;
	}
}

void SetWeights(float **Weights, int *WeightDim, double *sample_points, int *FirstIndex, int num_samples_in, int num_samples_out, int sample_rate_in, int C, int w){
	for(int i = 0; i < num_samples_out; i++){
		for(int j = 0; j < WeightDim[i]; j++){
			double Delta_t = sample_points[i] - (FirstIndex[i] + j) / static_cast<double>(sample_rate_in);
			Weights[i][j] = FilterFunc(Delta_t, C, w) / static_cast<double>(sample_rate_in);
		}
	}
}

int NumFramesAvailable(int SampleLen, int FrameLen, int FrameShift){
	return (SampleLen - FrameLen) / FrameShift + 1;
}

int LagLen(double min_lag, double max_lag, float delta_pitch){
	int Len = 0;
	for(double lag = min_lag; lag <= max_lag; lag *= (1.0 + delta_pitch))
		Len++;
	return Len;
}

void SelectLags(double min_lag, double max_lag, float delta_pitch, double *lags){
	int i;
	double lag;
	for(i = 0, lag = min_lag; lag <= max_lag; lag *= (1.0 + delta_pitch), i++)
		lags[i] = lag;
}

void ExtractFrame(float *DownSampledWave, int DownSampleLen, int Offset, float *Window, int full_frame_length){
	// Intialize window
	for(int i = 0; i < full_frame_length; i++)
		Window[i] = 0.0;

	if(Offset + full_frame_length > DownSampleLen){
		int SubWindow = DownSampleLen - Offset;
		ExtractFrame(DownSampledWave, DownSampleLen, Offset, Window, SubWindow);
		return;
	}

	for(int i = 0; i < full_frame_length; i++){
		Window[i] = DownSampledWave[i + Offset];
	}
	return;
}

void ComputeCorrelation(float *wave, int waveLen, int first_lag, int last_lag, int nccf_window_size, float *inner_prod, float *norm_prod){
	float zero_mean_wave[waveLen];
	memcpy(zero_mean_wave, wave, sizeof(zero_mean_wave));
	float temp = 0.0;
	for(int i = 0; i < nccf_window_size; i++)
		temp += wave[i];
	temp /= nccf_window_size;
	for(int i = 0; i < waveLen; i++)
		zero_mean_wave[i] -= temp;

	float e1, e2, sum;
	float subvec1[nccf_window_size];
	memcpy(subvec1, zero_mean_wave, sizeof(subvec1));
	e1 = InnerProduct(subvec1, subvec1, nccf_window_size);
	for(int lag = first_lag; lag <= last_lag; lag++){
		float subvec2[nccf_window_size];
		float *temp2 = zero_mean_wave + lag;
		memcpy(subvec2, temp2, sizeof(subvec2));
		e2 = InnerProduct(subvec2, subvec2, nccf_window_size);
		sum = InnerProduct(subvec1, subvec2, nccf_window_size);
		inner_prod[lag - first_lag] = sum;
		norm_prod[lag - first_lag] = e1 * e2;
	}
}

void ComputeNccf(float *inner_prod, float *norm_prod, double nccf_ballast, float *nccf_vec, int Len){
	for(int lag = 0; lag < Len; lag++){
		float numerator = inner_prod[lag],
			  denominator = pow(norm_prod[lag] + nccf_ballast, 0.5),
			  nccf;
		if(denominator != 0.0){
			nccf = numerator / denominator;
		}
		else{
			nccf = 0.0;
		}
		nccf_vec[lag] = nccf;
	}
}

void ComputeLocalCost(float *nccf_pitch, double *lags, float *local_cost, int num_states, float soft_min_f0){
	// from the paper, eq. 5, local_cost = 1 - Phi(t,i)(1 - soft_min_f0 L_i)
	// nccf is the nccf on this frame measured at the lags in "lags".
	for(int i = 0; i < num_states; i++)
		local_cost[i] = 1.0;
	for(int i = 0; i < num_states; i++)
		local_cost[i] -= 1.0 * nccf_pitch[i];
	for(int i = 0; i < num_states; i++)
		local_cost[i] += soft_min_f0 * lags[i] * nccf_pitch[i];
}

void ComputeBacktraces(float *nccf_pitch, double *lags, float *prev_forward_cost, vector<pair<int, int> > *index_info, float *this_forward_cost, int *state_info_backpointer, int num_states, float soft_min_f0, float penalty_factor, float delta_pitch){
	float local_cost[num_states];
	ComputeLocalCost(nccf_pitch, lags, local_cost, num_states, soft_min_f0);

	float delta_pitch_sq = pow(log(1.0 + delta_pitch), 2.0);
	float inter_frame_factor = delta_pitch_sq * penalty_factor;

	if(index_info->empty())
		index_info->resize(num_states);

	vector<pair<int, int> > &bounds = *index_info;

	int last_backpointer = 0;
	for(int i = 0; i < num_states; i++){
		int start_j = last_backpointer;
		float best_cost = (start_j - i) * (start_j - i) * inter_frame_factor + prev_forward_cost[start_j];
		int best_j = start_j;
		for(int j = start_j + 1; j < num_states; j++){
			float this_cost = (j - i) * (j - i) * inter_frame_factor + prev_forward_cost[j];
			if(this_cost < best_cost){
				best_cost = this_cost;
				best_j = j;
			}
			else{
				break;
			}
		}
		state_info_backpointer[i] = best_j;
		this_forward_cost[i] = best_cost;
		bounds[i].first = best_j;
		bounds[i].second = num_states - 1;
		last_backpointer = best_j;
	}
	for(int iter = 0; iter < num_states; iter++){
		bool changed = false;
		if(iter % 2 == 0){
			last_backpointer = num_states - 1;
			for(int i = num_states - 1; i >= 0; i--){
				int lower_bound = bounds[i].first;
				int upper_bound = min(last_backpointer, bounds[i].second);
				if(upper_bound == lower_bound){
					last_backpointer = lower_bound;
					continue;
				}
				float best_cost = this_forward_cost[i];
				int best_j = state_info_backpointer[i];
				int initial_best_j = best_j;

				if(best_j == upper_bound){
					last_backpointer = best_j;
					continue;
				}

				for(int j = upper_bound; j > lower_bound + 1; j--){
					float this_cost = (j - i) * (j - i) * inter_frame_factor + prev_forward_cost[j];
					if(this_cost < best_cost){
						best_cost = this_cost;
						best_j = j;
					}
					else{
						if(best_j > j)
							break;
					}
				}

				bounds[i].second = best_j;
				if(best_j != initial_best_j){
					this_forward_cost[i] = best_cost;
					state_info_backpointer[i] = best_j;
					changed = true;
				}
				last_backpointer = best_j;
			}
		}
		else{
			last_backpointer = 0;
			for(int i = 0; i < num_states; i++){
				int lower_bound = max(last_backpointer, bounds[i].first);
				int upper_bound = bounds[i].second;
				if(upper_bound == lower_bound){
					last_backpointer = lower_bound;
					continue;
				}

				float best_cost = this_forward_cost[i];
				int best_j = state_info_backpointer[i];
				int initial_best_j = best_j;

				if(best_j == lower_bound){
					last_backpointer = best_j;
					continue;
				}

				for(int j = lower_bound; j < upper_bound - 1; j++){
					float this_cost = (j - i) * (j - i) * inter_frame_factor + prev_forward_cost[j];
					if(this_cost < best_cost){
						best_cost = this_cost;
						best_j = j;
					}
					else{
						if(best_j < j)
							break;
					}
				}

				bounds[i].first = best_j;
				if(best_j != initial_best_j){
					this_forward_cost[i] = best_cost;
					state_info_backpointer[i] = best_j;
					changed = true;
				}
				last_backpointer = best_j;
			}
		}
		if(!changed)
			break;
	}

	for(int i = 0; i < num_states; i++)
		this_forward_cost[i] += local_cost[i];
}

void SetBestState(int best_state, pair<int, float> *lag_nccf, int NumofFrames, int **state_info_backpointer, float **nccf_pov_resampled){
	for(int i = NumofFrames - 1; i >= 0; i--){
		lag_nccf[i].first = best_state;
		lag_nccf[i].second = nccf_pov_resampled[i][best_state];
		best_state = state_info_backpointer[i][best_state];
	}
}

void GetFeature(float *f0, float *f1, pair<int, float> *lag_nccf, double *lags, int StartFrame, int NumofFrames){
	for(int i = 0; i < NumofFrames; i++){
		f0[StartFrame + i] = lag_nccf[i].second;
		f1[StartFrame + i] = 1.0 / lags[lag_nccf[i].first];
	}
}

/*
Optional parameters:
	min_f0: Minimum possible frequency value (Hz)
	max_f0: Maximum possible frequency value (Hz)
	window_width: Length in seconds of window used for NCCF
	window_shift: Frame_shift, in seconds (should match that used for baseline features e.g. PLP)
	soft_min_f0: Minimum f0, applied in soft way; 
				 must not exceed min_f0.
	nccf_ballast: Increasing this factor reduces NCCF for quiet frames, helping ensure pitch 
				  continuity in unvoiced regions
	penalty_factor: Factor that penalizes frequency change
	delta_pitch: Smallest relative change in pitch that our algorithm measures
	lowpass_cutoff: Low_pass cutoff that we apply to the raw signal
	lowpass_filter_width: Integer that determines filter width of low_pass filter 
						  (more gives wider filter with sharper cutoff)
	resample_frequency: Sample frequency for NCCF;
						must exceed twice lowpass_cutoff.
	upsample_filter_width: Integer that determines filter width when upsampling NCCF
*/

int ComputePitchFeature(int DataLength, short *Data, int sample_frequency, float *f0, float *POV,
						int min_f0=50, int max_f0=400, float penalty_factor=0.1, 
						float window_width=0.025, float window_shift=0.01, float soft_min_f0=10.0, 
						float nccf_ballast=7000, float delta_pitch=0.005, int lowpass_cutoff=1000,
						int lowpass_filter_width=1, int resample_frequency=4000, int upsample_filter_width=5){
	int sample_per_window = (int)(sample_frequency * window_width);
	int sample_per_shift = (int)(sample_frequency * window_shift);
	int head = 0, count = 0;
	while(head + sample_per_shift < DataLength){
		head += sample_per_shift;
		count++;
	}

	// Resample to resample_frequency
	int DownSampleLen = static_cast<int>(DataLength / (sample_frequency / static_cast<float>(resample_frequency))) + 1;
	float DownSampledWave[DownSampleLen];
	DownSampleLen = Resample(sample_frequency, resample_frequency, lowpass_cutoff, lowpass_filter_width, Data, DataLength, DownSampledWave);
	int LastFrameLen = lowpass_filter_width / (2.0 * lowpass_cutoff) * resample_frequency;
	DownSampleLen -= LastFrameLen;
	float *LastFrameWave = (float *)malloc(sizeof(float) * LastFrameLen);
	LastFrameWave = &DownSampledWave[DownSampleLen];

	// For Normalize
	double Sum = Sum_(DownSampledWave, DownSampleLen);
	double SumSquare = InnerProduct(DownSampledWave, DownSampledWave, DownSampleLen);
	double SumLast = Sum_(LastFrameWave, LastFrameLen);
	double SumSquareLast = InnerProduct(LastFrameWave, LastFrameWave, LastFrameLen);

	// Calculate some parameters
	double min_lag = 1.0 / max_f0;
	double max_lag = 1.0 / min_f0;
	int LagsLen = LagLen(min_lag, max_lag, delta_pitch);
	double lags[LagsLen];
	SelectLags(min_lag, max_lag, delta_pitch, lags);
	double outer_min_lag = 1.0 / max_f0 - (upsample_filter_width / (2.0 * resample_frequency));
	double outer_max_lag = 1.0 / min_f0 + (upsample_filter_width / (2.0 * resample_frequency));
	int nccf_first_lag = ceil(resample_frequency * outer_min_lag);
	int nccf_last_lag = myfloor(resample_frequency * outer_max_lag);
	int num_measured_lags = nccf_last_lag - nccf_first_lag + 1;
	int num_resampled_lags = LagsLen;
	int frame_shift = static_cast<int>(resample_frequency * window_shift);
	int basic_frame_length = static_cast<int>(resample_frequency * window_width);
	int full_frame_length = basic_frame_length + nccf_last_lag;

	// Calculate output length
	int NumofFrames = NumFramesAvailable(DownSampleLen, full_frame_length, frame_shift);
	int NumofFramesLast = NumFramesAvailable(DownSampleLen, basic_frame_length, frame_shift) - NumofFrames;
	int StartFrame = 0, EndFrame = NumofFrames;
	int StartFrameLast = NumofFrames, EndFrameLast = NumofFrames + NumofFramesLast;

	if(NumofFrames == 0){
		cout << "No frames output in pitch extraction" << endl;
		return 0;
	}
	
	// Prepare some arrays
	float Window[full_frame_length], 
		  inner_prod[num_measured_lags], 
		  norm_prod[num_measured_lags],
		  nccf_pitch[NumofFrames][num_measured_lags],
		  nccf_pitch_Last[NumofFramesLast][num_measured_lags],
		  nccf_pov[NumofFrames][num_measured_lags],
		  nccf_pov_Last[NumofFramesLast][num_measured_lags],
		  cur_forward_cost[num_resampled_lags];

	// Because the resampling of the NCCF is more efficient when grouped together,
	// we first compute the NCCF for all frames, then resample as a matrix, then
	// do the Viterbi.
	double mean_square = SumSquare / DownSampleLen - pow(Sum / DownSampleLen, 2.0);
	for(int frame = StartFrame; frame < EndFrame; frame++){
		int StartSample = frame * frame_shift;
		ExtractFrame(DownSampledWave, DownSampleLen, StartSample, Window, full_frame_length);
		ComputeCorrelation(Window, full_frame_length, nccf_first_lag, nccf_last_lag, basic_frame_length, inner_prod, norm_prod);
		double nccf_ballast_pov = 0.0,
			   nccf_ballast_pitch = pow(mean_square * basic_frame_length, 2) * nccf_ballast,
			   avg_norm_prod = Sum_(norm_prod, num_measured_lags) / num_measured_lags;
		float *nccf_pitch_row = nccf_pitch[frame - StartFrame];
		ComputeNccf(inner_prod, norm_prod, nccf_ballast_pitch, nccf_pitch_row, num_measured_lags);
		float *nccf_pov_row = nccf_pov[frame - StartFrame];
		ComputeNccf(inner_prod, norm_prod, nccf_ballast_pov, nccf_pov_row, num_measured_lags);	
	}

	// Upsampling
	int FirstIndex[num_resampled_lags];
	float *Weights[num_resampled_lags];
	int WeightDim[num_resampled_lags];
	double lags_[num_resampled_lags];
	double lags_offset = -nccf_first_lag / (1.0 * (resample_frequency));
	for(int i = 0; i < num_resampled_lags; i++)
		lags_[i] = lags[i] + lags_offset;
	SetFirstIndex(FirstIndex, Weights, lags_, num_measured_lags, num_resampled_lags, upsample_filter_width, resample_frequency * 0.5, resample_frequency, WeightDim);
	SetWeights(Weights, WeightDim, lags_, FirstIndex, num_measured_lags, num_resampled_lags, resample_frequency, resample_frequency * 0.5, upsample_filter_width);

	float nccf_pitch_resampled[NumofFrames][num_resampled_lags];
	// Prepare pointer to pointer to pass 2d array
	float *pointer2nccf_pitch[NumofFrames];
	for(int i = 0; i < NumofFrames; i++)
		pointer2nccf_pitch[i] = nccf_pitch[i];
	float *pointer2nccf_pitch_resampled[NumofFrames];
	for(int i = 0; i < NumofFrames; i++)
		pointer2nccf_pitch_resampled[i] = nccf_pitch_resampled[i];
	// Resample
	Resample(pointer2nccf_pitch, pointer2nccf_pitch_resampled, NumofFrames, num_measured_lags, num_resampled_lags, FirstIndex, WeightDim, Weights);
	
	float nccf_pov_resampled[NumofFrames][num_resampled_lags];
	// Prepare pointer to pointer to pass 2d array
	float *pointer2nccf_pov[NumofFrames];
	for(int i = 0; i < NumofFrames; i++)
		pointer2nccf_pov[i] = nccf_pov[i];
	float *pointer2nccf_pov_resampled[NumofFrames];
	for(int i = 0; i < NumofFrames; i++)
		pointer2nccf_pov_resampled[i] = nccf_pov_resampled[i];
	// Resample
	Resample(pointer2nccf_pov, pointer2nccf_pov_resampled, NumofFrames, num_measured_lags, num_resampled_lags, FirstIndex, WeightDim, Weights);

	float forward_cost[num_resampled_lags];
	for(int i = 0; i < num_resampled_lags; i++)
		forward_cost[i] = 0.0;
	
	// Run Viterbi
	vector<pair<int, int> > index_info;
	int backpointer[NumofFrames][num_resampled_lags];
	for(int frame = StartFrame; frame < EndFrame; frame++){
		int frame_index = frame - StartFrame;
		ComputeBacktraces(nccf_pitch_resampled[frame_index], lags, forward_cost, &index_info, cur_forward_cost, backpointer[frame_index], num_resampled_lags, soft_min_f0, penalty_factor, delta_pitch);
		Swap(forward_cost, cur_forward_cost, num_resampled_lags);
		// Renormalize forward_cost so smallest element is zero.
		float remainder = Min(forward_cost, num_resampled_lags);
		for(int i = 0; i < num_resampled_lags; i++)
			forward_cost[i] -= remainder;
	}

	// Prepare pointer to pointer to pass 2d array
	int *pointer2backpointer[NumofFrames];
	for(int i = 0; i < NumofFrames; i++)
		pointer2backpointer[i] = backpointer[i];

	// Trace back the best-path.
	int best_final_state = MinIndex(forward_cost, num_resampled_lags);
	pair<int, float> lag_nccf[NumofFrames];
	SetBestState(best_final_state, lag_nccf, NumofFrames, pointer2backpointer, pointer2nccf_pov_resampled);
	// for(int i = 0; i < NumofFrames; i++)
	// 	cout << lag_nccf[i].second << "\n";
	// cout << endl;

	// Get feature
	GetFeature(f0, POV, lag_nccf, lags, 0, NumofFrames);

	// For the last part
	mean_square = (SumSquare + SumSquareLast) / (DownSampleLen + LastFrameLen) - pow((Sum + SumLast) / (DownSampleLen + LastFrameLen), 2.0);
	for(int frame = StartFrameLast; frame < EndFrameLast; frame++){
		int StartSample = frame * frame_shift;
		ExtractFrame(DownSampledWave, DownSampleLen + LastFrameLen, StartSample, Window, full_frame_length);
		ComputeCorrelation(Window, full_frame_length, nccf_first_lag, nccf_last_lag, basic_frame_length, inner_prod, norm_prod);
		double nccf_ballast_pov = 0.0,
			   nccf_ballast_pitch = pow(mean_square * basic_frame_length, 2) * nccf_ballast,
			   avg_norm_prod = Sum_(norm_prod, num_measured_lags) / num_measured_lags;
		float *nccf_pitch_row_Last = nccf_pitch_Last[frame - StartFrameLast];
		ComputeNccf(inner_prod, norm_prod, nccf_ballast_pitch, nccf_pitch_row_Last, num_measured_lags);
		float *nccf_pov_row_Last = nccf_pov_Last[frame - StartFrameLast];
		ComputeNccf(inner_prod, norm_prod, nccf_ballast_pov, nccf_pov_row_Last, num_measured_lags);
	}


	float nccf_pitch_resampled_Last[NumofFramesLast][num_resampled_lags];
	// Prepare pointer to pointer to pass 2d array
	float *pointer2nccf_pitch_Last[NumofFramesLast];
	for(int i = 0; i < NumofFramesLast; i++)
		pointer2nccf_pitch_Last[i] = nccf_pitch_Last[i];
	float *pointer2nccf_pitch_resampled_Last[NumofFramesLast];
	for(int i = 0; i < NumofFramesLast; i++)
		pointer2nccf_pitch_resampled_Last[i] = nccf_pitch_resampled_Last[i];
	// Resample
	Resample(pointer2nccf_pitch_Last, pointer2nccf_pitch_resampled_Last, NumofFramesLast, num_measured_lags, num_resampled_lags, FirstIndex, WeightDim, Weights);
	
	float nccf_pov_resampled_Last[NumofFramesLast][num_resampled_lags];
	// Prepare pointer to pointer to pass 2d array
	float *pointer2nccf_pov_Last[NumofFramesLast];
	for(int i = 0; i < NumofFramesLast; i++)
		pointer2nccf_pov_Last[i] = nccf_pov_Last[i];
	float *pointer2nccf_pov_resampled_Last[NumofFramesLast];
	for(int i = 0; i < NumofFramesLast; i++)
		pointer2nccf_pov_resampled_Last[i] = nccf_pov_resampled_Last[i];
	// Resample
	Resample(pointer2nccf_pov_Last, pointer2nccf_pov_resampled_Last, NumofFramesLast, num_measured_lags, num_resampled_lags, FirstIndex, WeightDim, Weights);

	// Run Viterbi
	vector<pair<int, int> > index_info_Last;
	int backpointer_Last[NumofFramesLast][num_resampled_lags];
	for(int frame = StartFrameLast; frame < EndFrameLast; frame++){
		int frame_index = frame - StartFrameLast;
		ComputeBacktraces(nccf_pitch_resampled_Last[frame_index], lags, forward_cost, &index_info_Last, cur_forward_cost, backpointer_Last[frame_index], num_resampled_lags, soft_min_f0, penalty_factor, delta_pitch);
		Swap(forward_cost, cur_forward_cost, num_resampled_lags);
		// Renormalize forward_cost so smallest element is zero.
		float remainder = Min(forward_cost, num_resampled_lags);
		for(int i = 0; i < num_resampled_lags; i++)
			forward_cost[i] -= remainder;
	}

	int *pointer2backpointer_Last[NumofFramesLast];
	for(int i = 0; i < NumofFramesLast; i++)
		pointer2backpointer_Last[i] = backpointer_Last[i];

	// Trace back the best-path.
	int best_final_state_Last = MinIndex(forward_cost, num_resampled_lags);
	pair<int, float> lag_nccf_Last[NumofFramesLast];
	SetBestState(best_final_state_Last, lag_nccf_Last, NumofFramesLast, pointer2backpointer_Last, pointer2nccf_pov_resampled_Last);

	// Get feature
	GetFeature(f0, POV, lag_nccf_Last, lags, NumofFrames, NumofFramesLast);

	// return total feature dimension
	return NumofFrames + NumofFramesLast;
}

// For testing
int main(void){
	ifstream fin;
	short Data[MAXDATALEN];
	int DataLength = 0;
	fin.open("Data/M010101_0_1.pcm", ios::in | ios ::binary);
	while(fin.peek() != EOF){
		short s;
		fin.read((char *)&s, 2);
		Data[DataLength] = s;
		DataLength++;
	}

	float f0[MAXFEATLEN], POV[MAXFEATLEN];
	int FeatureLength = ComputePitchFeature(DataLength, Data, 16000, f0, POV, 75, 1000, 0.35);
	// Print out result
	for(int i = 0; i < FeatureLength; i++)
		cout << f0[i] << " " << POV[i] << " " << endl;

	return 0;
}