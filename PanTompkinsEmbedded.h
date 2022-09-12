/*
 * PanTompkins.h
 *
 *  Created on: Aug 22, 2022
 *      Author: ttronrud
 */

#ifndef INC_PANTOMPKINS_H_
#define INC_PANTOMPKINS_H_


#include <stdint.h>
typedef int dataType;

/*
 * Execute PT on the initialized signal
 * If a signal has been executed, it will re-run with the
 * final filter state (e.g. re-execute after learning the proper
 * R peak filtering system)
 * This won't reset the next_r variable or anything else,
 * so execution will continue where it left off. If using a new signal,
 * this might be intended. If it isn't, use sigswap.
 */
void PanTompkins();
/*
 * Set the input signal, length
 * This also sets the array of previous RRs and detected R peaks to zero
 * Thresholds are wiped, essentially reverting all learning for a new
 * signal (from a new source). Shouldn't be necessary more than once,
 * but wiping all learning from a poor signal may be useful, who knows.
 */
void PanTompkins_init(dataType *signal, uint32_t sig_len);
/*
 * Swap the active signal source - maintaining buffers and RR averages
 * but directing the flow at a new signal. Not necessary if you are using the
 * original signal input as a ring buffer or some other system that gets
 * overwritten.
 * r_index_start sets the index of the R-peak location array that detected
 * peaks will be written at. If you have used the locations from the last run,
 * set to zero. If you are planning to do an RRI calculation from pieces of ECG,
 * use -1 (or the end of the last R array).
 */
void PanTompkins_sigswap(dataType *new_signal, uint32_t new_sig_len, uint32_t new_index_offset, uint32_t r_index_start);

/*
 * Saving and loading of the filterstate variable.
 * Could be useful in the case of absolutely horrendous data
 * that ruins the state of the thresholds.
 */
void PanTompkins_save_filterstate();
void PanTompkins_load_filterstate();

#endif /* INC_PANTOMPKINS_H_ */
