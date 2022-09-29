/**
 * ------------------------------------------------------------------------------*
 * File: PanTompkinsEmbedded.c                                                           *

 *       ANSI-C implementation of Pan-Tompkins real-time QRS detection algorithm *
 * Author: 
        Rafael de Moura Moreira <rafaelmmoreira@gmail.com>                       *
        Thorold Tronrud <ttronrud@starfishmedical.com>
 * License: MIT License                                                          *
 * ------------------------------------------------------------------------------*
 * ---------------------------------- HISTORY ---------------------------------- *
 *    date   |    author    |                     description                    *
 * ----------| -------------| ---------------------------------------------------*
 * 2022/09/22| Thor J. T.   | - Added embedded implementation with moving average*
 *           |              |   filter system.                                   *
 * 2019/04/11| Rafael M. M. | - Fixed moving-window integral.                    *
 *           |              | - Fixed how to find the correct sample with the    *
 *           |              | last QRS.                                          *
 *           |              | - Replaced constant value in code by its #define.  *
 *           |              | - Added some casting on comparisons to get rid of  *
 *           |              | compiler warnings.                                 *
 * 2019/04/15| Rafael M. M. | - Removed delay added to the output by the filters.*
 *           |              | - Fixed multiple detection of the same peak.       *
 * 2019/04/16| Rafael M. M. | - Added output buffer to correctly output a peak   *
 *           |              | found by back searching using the 2nd thresholds.  *
 * 2019/04/23| Rafael M. M. | - Improved comparison of slopes.                   *
 *           |              | - Fixed formula to obtain the correct sample from  *
 *           |              | the buffer on the back search.                     *
 * ------------------------------------------------------------------------------*
 * MIT License                                                                   *
 *                                                                               *
 * Copyright (c) 2018 Rafael de Moura Moreira                                    *
 *                                                                               *
 * Permission is hereby granted, free of charge, to any person obtaining a copy  *
 * of this software and associated documentation files (the "Software"), to deal *
 * in the Software without restriction, including without limitation the rights  *
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell     *
 * copies of the Software, and to permit persons to whom the Software is         *
 * furnished to do so, subject to the following conditions:                      *
 *                                                                               *
 * The above copyright notice and this permission notice shall be included in all*
 * copies or substantial portions of the Software.                               *
 *                                                                               *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR    *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,      *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE   *
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER        *
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, *
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE *
 * SOFTWARE.                                                                     *
 *-------------------------------------------------------------------------------*
 * Description                                                                   *
 *                                                                               *
 * The main goal of this implementation is to be easy to port to different opera-*
 * ting systems, as well as different processors and microcontrollers, including *
 * embedded systems. It can work both online or offline, depending on whether all*
 * the samples are available or not - it can be adjusted on the input function.  *
 *                                                                               *
 * For a signal peak to be recognized as a fiducial point, its correspondent va- *
 * lue on both the filtered signal and the integrator must be above a certain    *
 * threshold. Additionally, there are time-restraints to prevent a T-wave from   *
 * being mistakenly identified as an R-peak: a hard 200ms restrain (a new peak   *
 * 200ms from the previous one is, necessarily, a T-wave) and a soft 360ms res-  *
 * train (the peak's squared slope must also be very high to be considered as a  *
 * real peak).                                                                   *
 *                                                                               *
 * When a peak candidate is discarded, its value is used to update the noise     *
 * thresholds - which are also used to estimate the peak thresholds.             *
 *                                                                               *
 * Two buffers keep 8 RR-intervals to calculate RR-averages: one of them keeps   *
 * the last 8 RR-intervals, while the other keeps only the RR-intervals that res-*
 * pect certain restrictions. If both averages are equal, the heart pace is con- *
 * sidered normal. If the heart rate isn't normal, the thresholds change to make *
 * it easier to detect possible weaker peaks. If no peak is detected for a long  *
 * period of time, the thresholds also change and the last discarded peak candi- *
 * date is reconsidered.                                                         *
 *-------------------------------------------------------------------------------*/
 
//Use an ifdef to skip compiling this section if the DAQ isn't
//being used for something that needs QRS detection --
//so we can have embedded R detection available without the constant memory usage
#define USING_PAN_TOMPKINS


#include "PanTompkinsEmbedded.h"
#ifdef USING_PAN_TOMPKINS
//These includes may not be available in an embedded context
//depending on how modern your system is
#include <stdlib.h>
//useful for memset
#include <string.h>

#define NOSAMPLE -32000 // An indicator that there are no more samples to read. Use an impossible value for a sample.
//#define FS 250          // Sampling frequency.
#define BUFFSIZE 415    // The size of the buffers (in samples). Must fit more than 1.66 times an RR interval, which
// typically could be around 1 second.

#define DELAY 14		// Delay introduced by the filters. Filter only output samples after this one.
// Set to 0 if you want to keep the delay. Fixing the delay results in DELAY less samples
// in the final end result.
#define  MOVING_AVG_LEN 5

#define WINDOWSIZE 40

#define FS 250

#define OVERWRITE_RS_ON_OVERFLOW

uint32_t sample = 0;

//define these globally, so they're maintained:
// The signal array is where the most recent samples are kept. The other arrays are the outputs of each
// filtering module: DC Block, low pass, high pass, integral etc.
// The output is a buffer where we can change a previous result (using a back search) before outputting.
dataType signal[BUFFSIZE], dcblock[BUFFSIZE], lowpass[BUFFSIZE], highpass[BUFFSIZE], derivative[BUFFSIZE], squared[BUFFSIZE], integral[BUFFSIZE], outputSignal[BUFFSIZE];

// rr1 holds the last 8 RR intervals. rr2 holds the last 8 RR intervals between rrlow and rrhigh.
// rravg1 is the rr1 average, rr2 is the rravg2. rrlow = 0.92*rravg2, rrhigh = 1.08*rravg2 and rrmiss = 1.16*rravg2.
// rrlow is the lowest RR-interval considered normal for the current heart beat, while rrhigh is the highest.
// rrmiss is the longest that it would be expected until a new QRS is detected. If none is detected for such
// a long interval, the thresholds must be adjusted.
int rr1[8], rr2[8], rravg1 = 0, rravg2 = 0, rrlow = 0, rrhigh = 0, rrmiss = 0;


// There are the variables from the original Pan-Tompkins algorithm.
// The ones ending in _i correspond to values from the integrator.
// The ones ending in _f correspond to values from the DC-block/low-pass/high-pass filtered signal.
// The peak variables are peak candidates: signal values above the thresholds.
// The threshold 1 variables are the threshold variables. If a signal sample is higher than this threshold, it's a peak.
// The threshold 2 variables are half the threshold 1 ones. They're used for a back search when no peak is detected for too long.
// The spk and npk variables are, respectively, running estimates of signal and noise peaks.
dataType peak_i = 0, peak_f = 0, threshold_i1 = 0, threshold_i2 = 0, threshold_f1 = 0, threshold_f2 = 0, spk_i = 0, spk_f = 0, npk_i = 0, npk_f = 0;

FilterState saved_filter_state;

// This variable is used as an index to work with the signal buffers. If the buffers still aren't
// completely filled, it shows the last filled position. Once the buffers are full, it'll always
// show the last position, and new samples will make the buffers shift, discarding the oldest
// sample and storing the newest one on the last position.
int current;

//CURRENT STATE OF THE PT RUN
// i and j are iterators for loops.
// sample counts how many samples have been read so far.
// lastQRS stores which was the last sample read when the last R sample was triggered.
// lastSlope stores the value of the squared slope when the last R sample was triggered.
// currentSlope helps calculate the max. square slope for the present sample.
// These are all uint32_t so that very long signals can be read without messing the count.
uint32_t i, j, lastQRS = 0, lastSlope = 0, currentSlope = 0;
dataType last_avg[MOVING_AVG_LEN];
dataType avg = 0;

// qrs tells whether there was a detection or not.
// regular tells whether the heart pace is regular or not.
// prevRegular tells whether the heart beat was regular before the newest RR-interval was calculated.
bool qrs, regular = true, prevRegular;


void PanTompkins_setsaved_filterstate(FilterState *fs)
{
    if(fs == NULL)
        return; //don't try to access a null ptr
    for(int i = 0; i < 8; i++)
	{
		saved_filter_state.rr1[i] = fs->rr1[i];
		saved_filter_state.rr2[i] = fs->rr2[i];
	}
	saved_filter_state.rravg1 = fs->rravg1;
	saved_filter_state.rravg2 = fs->rravg2;
	saved_filter_state.rrlow = fs->rrlow;
	saved_filter_state.rrhigh = fs->rrhigh;
	saved_filter_state.rrmiss = fs->rrmiss;
	saved_filter_state.peak_i = fs->peak_i;
	saved_filter_state.peak_f = fs->peak_f;
	saved_filter_state.threshold_i1 = fs->threshold_i1;
	saved_filter_state.threshold_i2 = fs->threshold_i2;
	saved_filter_state.threshold_f1 = fs->threshold_f1;
	saved_filter_state.threshold_f2 = fs->threshold_f2;
	saved_filter_state.spk_i = fs->spk_i;
	saved_filter_state.spk_f = fs->spk_f;
	saved_filter_state.npk_i = fs->npk_i;
	saved_filter_state.npk_f = fs->npk_f;
}

void PanTompkins_exportsaved_filterstate(FilterState *fs)
{
    if(fs == NULL)
        return; //don't try to access a null ptr
    for(int i = 0; i < 8; i++)
	{
		fs->rr1[i] = saved_filter_state.rr1[i];
		fs->rr2[i] = saved_filter_state.rr2[i];
	}
	fs->rravg1 = saved_filter_state.rravg1;
	fs->rravg2 = saved_filter_state.rravg2;
	fs->rrlow = saved_filter_state.rrlow;
	fs->rrhigh = saved_filter_state.rrhigh;
	fs->rrmiss = saved_filter_state.rrmiss;
	fs->peak_i = saved_filter_state.peak_i;
	fs->peak_f = saved_filter_state.peak_f;
	fs->threshold_i1 = saved_filter_state.threshold_i1;
	fs->threshold_i2 = saved_filter_state.threshold_i2;
	fs->threshold_f1 = saved_filter_state.threshold_f1;
	fs->threshold_f2 = saved_filter_state.threshold_f2;
	fs->spk_i = saved_filter_state.spk_i;
	fs->spk_f = saved_filter_state.spk_f;
	fs->npk_i = saved_filter_state.npk_i;
	fs->npk_f = saved_filter_state.npk_f;
}
void PanTompkins_save_filterstate()
{
	for(int i = 0; i < 8; i++)
	{
		saved_filter_state.rr1[i] = rr1[i];
		saved_filter_state.rr2[i] = rr2[i];
	}
	saved_filter_state.rravg1 = rravg1;
	saved_filter_state.rravg2 = rravg2;
	saved_filter_state.rrlow = rrlow;
	saved_filter_state.rrhigh = rrhigh;
	saved_filter_state.rrmiss = rrmiss;
	saved_filter_state.peak_i = peak_i;
	saved_filter_state.peak_f = peak_f;
	saved_filter_state.threshold_i1 = threshold_i1;
	saved_filter_state.threshold_i2 = threshold_i2;
	saved_filter_state.threshold_f1 = threshold_f1;
	saved_filter_state.threshold_f2 = threshold_f2;
	saved_filter_state.spk_i = spk_i;
	saved_filter_state.spk_f = spk_f;
	saved_filter_state.npk_i = npk_i;
	saved_filter_state.npk_f = npk_f;
}

void PanTompkins_load_filterstate()
{
	for(int i = 0; i < 8; i++)
	{
		rr1[i] = saved_filter_state.rr1[i];
		rr2[i] = saved_filter_state.rr2[i];
	}
	rravg1 = saved_filter_state.rravg1;
	rravg2 = saved_filter_state.rravg2;
	rrlow = saved_filter_state.rrlow;
	rrhigh = saved_filter_state.rrhigh;
	rrmiss = saved_filter_state.rrmiss;
	peak_i = saved_filter_state.peak_i;
	peak_f = saved_filter_state.peak_f;
	threshold_i1 = saved_filter_state.threshold_i1;
	threshold_i2 = saved_filter_state.threshold_i2;
	threshold_f1 = saved_filter_state.threshold_f1;
	threshold_f2 = saved_filter_state.threshold_f2;
	spk_i = saved_filter_state.spk_i;
	spk_f = saved_filter_state.spk_f;
	npk_i = saved_filter_state.npk_i;
	npk_f = saved_filter_state.npk_f;
}

void PanTompkins_init()
{
	//Direct PT to the correct input source
	sample = 0;

	//reset all the detection thresholds
	//and RR data
    rravg1 = 0;
	rravg2 = 0;
	rrlow = 0;
	rrhigh = 0;
	rrmiss = 0;
	//Thresholds for R peak detection
	peak_i = 0;
	peak_f = 0;
	threshold_i1 = 0;
	threshold_i2 = 0;
	threshold_f1 = 0;
	threshold_f2 = 0;
	spk_i = 0;
	spk_f = 0;
	npk_i = 0;
	npk_f = 0;
    // Initializing the RR averages
    for (uint32_t i = 0; i < 8; i++)
    {
        rr1[i] = 0;
        rr2[i] = 0;
    }
    memset(last_avg,0,MOVING_AVG_LEN*sizeof(dataType));
}

bool PanTompkins_SingleStep(dataType inputSample)
{
	if (sample >= BUFFSIZE)
	{
		for (i = 0; i < BUFFSIZE - 1; i++)
		{
			signal[i] = signal[i+1];
			dcblock[i] = dcblock[i+1];
			lowpass[i] = lowpass[i+1];
			highpass[i] = highpass[i+1];
			derivative[i] = derivative[i+1];
			squared[i] = squared[i+1];
			integral[i] = integral[i+1];
			outputSignal[i] = outputSignal[i+1];
		}
		current = BUFFSIZE - 1;
	}
	else
	{
		current = sample;
	}
	signal[current] = inputSample;
	//try an ultra-lightweight moving average
	//set the latest value
	last_avg[MOVING_AVG_LEN-1] = signal[current];
	if(sample > MOVING_AVG_LEN && signal[current] != NOSAMPLE)
	{
		avg = 0;
		for(int q = 0; q < MOVING_AVG_LEN; q++)
		{
			avg += last_avg[q]/MOVING_AVG_LEN;
			if(q < MOVING_AVG_LEN-1)
			{
				last_avg[q] = last_avg[q+1]; //shift average buffer down one for next iteration
			}
		}
		signal[current] -= avg;

	}
	// If no sample was read, stop processing!
	if (signal[current] == NOSAMPLE)
		return false;
	sample++; // Update sample counter

	// DC Block filter
	// This was not proposed on the original paper.
	// It is not necessary and can be removed if your sensor or database has no DC noise.
	if (current >= 1)
		dcblock[current] = signal[current] - signal[current-1] + 0.995*dcblock[current-1];
	else
		dcblock[current] = 0;

	// Low Pass filter
	// Implemented as proposed by the original paper.
	// y(nT) = 2y(nT - T) - y(nT - 2T) + x(nT) - 2x(nT - 6T) + x(nT - 12T)
	// Can be removed if your signal was previously filtered, or replaced by a different filter.
	lowpass[current] = dcblock[current];
	if (current >= 1)
		lowpass[current] += 2*lowpass[current-1];
	if (current >= 2)
		lowpass[current] -= lowpass[current-2];
	if (current >= 6)
		lowpass[current] -= 2*dcblock[current-6];
	if (current >= 12)
		lowpass[current] += dcblock[current-12];

	// High Pass filter
	// Implemented as proposed by the original paper.
	// y(nT) = 32x(nT - 16T) - [y(nT - T) + x(nT) - x(nT - 32T)]
	// Can be removed if your signal was previously filtered, or replaced by a different filter.
	highpass[current] = -lowpass[current];
	if (current >= 1)
		highpass[current] -= highpass[current-1];
	if (current >= 16)
		highpass[current] += 32*lowpass[current-16];
	if (current >= 32)
		highpass[current] += lowpass[current-32];

	// Derivative filter
	// This is an alternative implementation, the central difference method.
	// f'(a) = [f(a+h) - f(a-h)]/2h
	// The original formula used by Pan-Tompkins was:
	// y(nT) = (1/8T)[-x(nT - 2T) - 2x(nT - T) + 2x(nT + T) + x(nT + 2T)]
	derivative[current] = highpass[current];
	if (current > 0)
		derivative[current] -= highpass[current-1];

	// This just squares the derivative, to get rid of negative values and emphasize high frequencies.
	// y(nT) = [x(nT)]^2.
	squared[current] = derivative[current]*derivative[current];

	// Moving-Window Integration
	// Implemented as proposed by the original paper.
	// y(nT) = (1/N)[x(nT - (N - 1)T) + x(nT - (N - 2)T) + ... x(nT)]
	// WINDOWSIZE, in samples, must be defined so that the window is ~150ms.

	integral[current] = 0;
	for (i = 0; i < WINDOWSIZE; i++)
	{
		if (current >= (dataType)i)
			integral[current] += squared[current - i];
		else
			break;
	}
	integral[current] /= (dataType)i;

	qrs = false;

	// If the current signal is above one of the thresholds (integral or filtered signal), it's a peak candidate.
	if (integral[current] >= threshold_i1 || highpass[current] >= threshold_f1)
	{
		peak_i = integral[current];
		peak_f = highpass[current];
	}

	// If both the integral and the signal are above their thresholds, they're probably signal peaks.
	if ((integral[current] >= threshold_i1) && (highpass[current] >= threshold_f1))
	{
		// There's a 200ms latency. If the new peak respects this condition, we can keep testing.
		if (sample > lastQRS + FS/5)
		{
			// If it respects the 200ms latency, but it doesn't respect the 360ms latency, we check the slope.
			if (sample <= lastQRS + (uint32_t)(0.36*FS))
			{
				// The squared slope is "M" shaped. So we have to check nearby samples to make sure we're really looking
				// at its peak value, rather than a low one.
				currentSlope = 0;
				for (j = current - 10; j <= current; j++)
					if (squared[j] > currentSlope)
						currentSlope = squared[j];

				if (currentSlope <= (dataType)(lastSlope/2))
				{
					qrs = false;
				}

				else
				{
					spk_i = 0.125*peak_i + 0.875*spk_i;
					threshold_i1 = npk_i + 0.25*(spk_i - npk_i);
					threshold_i2 = 0.5*threshold_i1;

					spk_f = 0.125*peak_f + 0.875*spk_f;
					threshold_f1 = npk_f + 0.25*(spk_f - npk_f);
					threshold_f2 = 0.5*threshold_f1;

					lastSlope = currentSlope;
					qrs = true;
				}
			}
				// If it was above both thresholds and respects both latency periods, it certainly is a R peak.
			else
			{
				currentSlope = 0;
				for (j = current - 10; j <= current; j++)
					if (squared[j] > currentSlope)
						currentSlope = squared[j];

				spk_i = 0.125*peak_i + 0.875*spk_i;
				threshold_i1 = npk_i + 0.25*(spk_i - npk_i);
				threshold_i2 = 0.5*threshold_i1;

				spk_f = 0.125*peak_f + 0.875*spk_f;
				threshold_f1 = npk_f + 0.25*(spk_f - npk_f);
				threshold_f2 = 0.5*threshold_f1;

				lastSlope = currentSlope;
				qrs = true;
			}
		}
			// If the new peak doesn't respect the 200ms latency, it's noise. Update thresholds and move on to the next sample.
		else
		{
			peak_i = integral[current];
			npk_i = 0.125*peak_i + 0.875*npk_i;
			threshold_i1 = npk_i + 0.25*(spk_i - npk_i);
			threshold_i2 = 0.5*threshold_i1;
			peak_f = highpass[current];
			npk_f = 0.125*peak_f + 0.875*npk_f;
			threshold_f1 = npk_f + 0.25*(spk_f - npk_f);
			threshold_f2 = 0.5*threshold_f1;
			qrs = false;
		}

	}
	// If a R-peak was detected, the RR-averages must be updated.
	if (qrs)
	{
		// Add the newest RR-interval to the buffer and get the new average.
		rravg1 = 0;
		for (i = 0; i < 7; i++)
		{
			rr1[i] = rr1[i+1];
			rravg1 += rr1[i];
		}
		rr1[7] = sample - lastQRS;
		lastQRS = sample;
		rravg1 += rr1[7];
		rravg1 *= 0.125;

		// If the newly-discovered RR-average is normal, add it to the "normal" buffer and get the new "normal" average.
		// Update the "normal" beat parameters.
		if ( (rr1[7] >= rrlow) && (rr1[7] <= rrhigh) )
		{
			rravg2 = 0;
			for (i = 0; i < 7; i++)
			{
				rr2[i] = rr2[i+1];
				rravg2 += rr2[i];
			}
			rr2[7] = rr1[7];
			rravg2 += rr2[7];
			rravg2 *= 0.125;
			rrlow = 0.92*rravg2;
			rrhigh = 1.16*rravg2;
			rrmiss = 1.66*rravg2;
		}

		prevRegular = regular;
		if (rravg1 == rravg2)
		{
			regular = true;
		}
			// If the beat had been normal but turned odd, change the thresholds.
		else
		{
			regular = false;
			if (prevRegular)
			{
				threshold_i1 /= 2;
				threshold_f1 /= 2;
			}
		}
	}
		// If no R-peak was detected, it's important to check how long it's been since the last detection.
	else
	{
		// If no R-peak was detected for too long, use the lighter thresholds and do a back search.
		// However, the back search must respect the 200ms limit and the 360ms one (check the slope).
		if ((sample - lastQRS > (uint32_t)rrmiss) && (sample > lastQRS + FS/5))
		{
			for (i = current - (sample - lastQRS) + FS/5; i < (uint32_t)current; i++)
			{
				if ( (integral[i] > threshold_i2) && (highpass[i] > threshold_f2))
				{
					currentSlope = 0;
					for (j = i - 10; j <= i; j++)
						if (squared[j] > currentSlope)
							currentSlope = squared[j];

					if ((currentSlope < (dataType)(lastSlope/2)) && (i + sample) < lastQRS + 0.36*lastQRS)
					{
						qrs = false;
					}
					else
					{
						peak_i = integral[i];
						peak_f = highpass[i];
						spk_i = 0.25*peak_i+ 0.75*spk_i;
						spk_f = 0.25*peak_f + 0.75*spk_f;
						threshold_i1 = npk_i + 0.25*(spk_i - npk_i);
						threshold_i2 = 0.5*threshold_i1;
						lastSlope = currentSlope;
						threshold_f1 = npk_f + 0.25*(spk_f - npk_f);
						threshold_f2 = 0.5*threshold_f1;
						// If a signal peak was detected on the back search, the RR attributes must be updated.
						// This is the same thing done when a peak is detected on the first try.
						//RR Average 1
						rravg1 = 0;
						for (j = 0; j < 7; j++)
						{
							rr1[j] = rr1[j+1];
							rravg1 += rr1[j];
						}
						rr1[7] = sample - (current - i) - lastQRS;
						qrs = true;
						lastQRS = sample - (current - i);
						rravg1 += rr1[7];
						rravg1 *= 0.125;

						//RR Average 2
						if ( (rr1[7] >= rrlow) && (rr1[7] <= rrhigh) )
						{
							rravg2 = 0;
							for (i = 0; i < 7; i++)
							{
								rr2[i] = rr2[i+1];
								rravg2 += rr2[i];
							}
							rr2[7] = rr1[7];
							rravg2 += rr2[7];
							rravg2 *= 0.125;
							rrlow = 0.92*rravg2;
							rrhigh = 1.16*rravg2;
							rrmiss = 1.66*rravg2;
						}

						prevRegular = regular;
						if (rravg1 == rravg2)
						{
							regular = true;
						}
						else
						{
							regular = false;
							if (prevRegular)
							{
								threshold_i1 /= 2;
								threshold_f1 /= 2;
							}
						}

						break;
					}
				}
			}
		}

		// Definitely no signal peak was detected.
		if (!qrs)
		{
			// If some kind of peak had been detected, then it's certainly a noise peak. Thresholds must be updated accordinly.
			if ((integral[current] >= threshold_i1) || (highpass[current] >= threshold_f1))
			{
				peak_i = integral[current];
				npk_i = 0.125*peak_i + 0.875*npk_i;
				threshold_i1 = npk_i + 0.25*(spk_i - npk_i);
				threshold_i2 = 0.5*threshold_i1;
				peak_f = highpass[current];
				npk_f = 0.125*peak_f + 0.875*npk_f;
				threshold_f1 = npk_f + 0.25*(spk_f - npk_f);
				threshold_f2 = 0.5*threshold_f1;
			}
		}
	}
	return qrs;
}


#else
//If we're not using PT functionality, use as little memory as possible, just allow it to compile
bool PanTompkins_SingleStep(dataType inputSample){ return false; }
void PanTompkins_init(){}
void PanTompkins_setsaved_filterstate(FilterState *fs);
void PanTompkins_exportsaved_filterstate(FilterState *fs);
void PanTompkins_save_filterstate(){}
void PanTompkins_load_filterstate(){}
#endif
