/**
 * ------------------------------------------------------------------------------*
 * File: panTompkins.c                                                           *

 *       ANSI-C implementation of Pan-Tompkins real-time QRS detection algorithm *
 * Author: Rafael de Moura Moreira <rafaelmmoreira@gmail.com>                    *
 * License: MIT License                                                          *
 * ------------------------------------------------------------------------------*
 * ---------------------------------- HISTORY ---------------------------------- *
 *    date   |    author    |                     description                    *
 * ----------| -------------| ---------------------------------------------------*
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
 * The main function, panTompkings(), calls input() to get the next sample and   *
 * store it in a buffer. Then it runs through a chain of filters: DC block, low  *
 * pass @ 15 Hz and high pass @ 5Hz. The filtered signal goes both through a de- *
 * rivative filter, which output is then squared, and through a windowed-integra-*
 * tor.                                                                          *
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
 *-------------------------------------------------------------------------------*
 * Instructions                                                                  *
 *                                                                               *
 * Here's what you should change to adjust the code to your needs:               *
 *                                                                               *
 * On panTompkins.h:                                                             *
 * - typedef int dataType;                                                       *
 * Change it from 'int' to whatever format your data is (float, unsigned int etc)*
 *                                                                               *
 * On both panTompkins.h and panTompkins.c:                                      *
 * - void init(char file_in[], char file_out[]);                                 *
 * This function is meant to do any kind of initial setup, such as starting a    *
 * serial connection with an ECG sensor. Change its parameters to whatever info  *
 * you need and its content. The test version included here loads 2 plain text   *
 * files: an input file, with the signal as a list of integer numbers in ASCII   *
 * format and an output file where either 0 or 1 will be written for each sample,*
 * whether a peak was detected or not.                                           *
 *                                                                               *
 * On panTompkins.c:                                                             *
 * - #define WINDOWSIZE                                                          *
 * Defines the size of the integration window. The original authors suggest on   *
 * their 1985 paper a 150ms window.                                              *
 *                                                                               *
 * - #define FS                                                                  *
 * Defines the sampling frequency.                                               *
 *                                                                               *
 * - #define NOSAMPLE                                                            *
 * A value to indicate you don't have any more samples to read. Choose a value   *
 * which a sample couldn't possibly have (e.g.: a negative value if your A/D con-*
 * verter only works with positive integers).                                    *
 *                                                                               *
 * - #define BUFFSIZE                                                            *
 * The size of the signal buffers. It should fit at least 1.66 times an RR-inter-*
 * val. Heart beats should be between 60 and 80 BPS for humans. So, considering  *
 * 1.66 times 1 second should be safe.                                           *
 *                                                                               *
 * - #define DELAY 22                                                            *
 * The delay introduced to the output signal. The first DELAY samples will be ig-*
 * nored, as the filters add a delay to the output signal, causing a mismatch    *
 * between the input and output signals. It's easier to compare them this way.   *
 * If you need them both to have the same amount of samples, set this to 0. If   *
 * you're working with different filters and/or sampling rates, you might need to*
 * adjust this value.                                                            *
 *                                                                               *
 * - #include <stdio.h>                                                          *
 * The file, as it is, both gets its inputs and sends its outputs to files. It   *
 * works on both Windows and Linux. If your source isn't a file, and/or your sys-*
 * tem doesn't have the <stdio.h> header, remove it.                             *
 * Include any other headers you might need to make your implementation work,    *
 * such as hardware libraries provided by your microcontroller manufacturer.     *
 *                                                                               *
 * - The input() function                                                        *
 * Change it to get the next sample from your source (a file, a serial device etc*
 * previously set up in your init() function. Return the sample value or NOSAMPLE*
 * if there are no more available samples.                                       *
 *                                                                               *
 * - The output() function                                                       *
 * Change it to output whatever you see fit, however you see fit: an RR-interval *
 * (which can be sent as a parameter to your function using the RR arrays), the  *
 * index of sample or timestamp which caused a R peak, whether a sample was a R  *
 * peak or not etc, and it can be written on a file, displayed on screen, blink a*
 * LED etc.                                                                      *
 *                                                                               *
 * - The panTompkins() function                                                  *
 * This function is almost entirely ANSI C, which means it should work as is on  *
 * most platforms. The only lines you really have to change are the fclose() ones*
 * at the very end, which are only here to allow testing of the code on systems  *
 * such as Windows and Linux for PC. You may wish to create extra variables or   *
 * increase the buffers' size as you see fit to extract different kinds of infor-*
 * mation to output, or fine tune the detection as you see fit - such as adding  *
 * different conditions for verification, initializing self-updating variables   *
 * with empirically-obtained values or changing the filters.                     *
 *-------------------------------------------------------------------------------*
 */
/*
 * Modifications made by Thor Tronrud, StarFish Medical
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * I had to add a moving average filter to the code to make it better-handle noise.
 *
 * I made the buffer arrays global scope, so the filter states don't
 * get wiped with each new run, meaning you can feed it pieces of an ECG and it
 * won't need to re-learn the filtering thresholds.
 */
//Use an ifdef to skip compiling this section if the DAQ isn'tan
//being used for something that needs QRS detection --
//so we can have embedded R detection without the constant memory usage
#define USING_PAN_TOMPKINS


#include "PanTompkinsEmbedded.h"
#ifdef USING_PAN_TOMPKINS
#include <stdlib.h>
#include <stdbool.h>
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


//hardcode for test
//We swap out the signals for the second "iteration" to ensure the state is actually maintained
//if things work correctly, even the first peak of the second iteration should be
//accurately determined!
//Lengths are 2000 each, sample rate 250Hz
//dataType SIGNAL1[] = {-121, -114, -117, -110, -109, -95, -77, -100, -102, -101, -105, -97, -81, -79, -65, -59, -57, -66, -65, -55, -47, -50, -56, -47, -45, -43, -61, -55, -38, -23, -23, -25, -26, -43, -31, -29, -33, -26, -26, -29, -19, -9, -9, -17, -15, -8, -1, 11, 15, 1, -5, -6, 8, -3, -5, -5, -22, -29, -21, -13, -13, -11, -7, -9, -9, -15, -13, -13, -13, -19, -33, -33, -26, -31, -33, -31, -33, -25, -29, -37, -31, -45, -47, -47, -50, -51, -55, -55, -57, -55, -61, -65, -61, -61, -57, -61, -70, -57, -57, -63, -65, -57, -65, -71, -69, -57, -61, -74, -94, -90, -77, -69, -61, -63, -61, -54, -47, -40, -49, -37, -21, 1, 7, 29, 23, 3, -21, -42, -29, -7, -21, -45, -73, -90, -98, -105, -126, -125, -129, -114, -118, -124, -129, -125, -126, -121, -130, -158, -138, -86, 48, 199, 373, 504, 511, 387, 147, -142, -288, -339, -373, -392, -387, -387, -367, -321, -264, -229, -207, -185, -179, -181, -166, -162, -142, -141, -144, -137, -137, -141, -148, -146, -141, -135, -137, -137, -141, -141, -141, -138, -138, -137, -133, -125, -125, -126, -116, -109, -106, -109, -109, -105, -105, -105, -114, -105, -102, -102, -101, -106, -102, -97, -107, -103, -112, -110, -111, -110, -109, -110, -114, -116, -115, -113, -106, -110, -112, -109, -106, -102, -105, -106, -103, -110, -110, -109, -108, -111, -121, -122, -125, -127, -130, -134, -131, -129, -129, -126, -133, -137, -140, -139, -141, -135, -137, -138, -146, -139, -141, -146, -144, -150, -150, -153, -158, -160, -160, -162, -170, -170, -169, -174, -169, -170, -172, -177, -195, -184, -184, -185, -196, -194, -201, -194, -193, -196, -201, -205, -207, -195, -196, -192, -182, -188, -182, -180, -183, -167, -150, -147, -145, -125, -121, -126, -141, -158, -170, -156, -144, -149, -177, -200, -224, -252, -250, -250, -251, -256, -246, -246, -248, -251, -245, -247, -258, -289, -300, -275, -185, -23, 179, 373, 493, 465, 328, 51, -267, -407, -408, -403, -419, -435, -423, -385, -347, -315, -292, -277, -271, -273, -268, -259, -260, -249, -262, -247, -248, -246, -251, -243, -243, -243, -241, -235, -239, -230, -231, -224, -233, -231, -227, -220, -214, -214, -209, -209, -212, -212, -207, -190, -192, -194, -181, -180, -185, -181, -180, -173, -183, -174, -177, -174, -170, -174, -170, -172, -173, -167, -166, -164, -166, -168, -157, -153, -153, -150, -140, -142, -146, -147, -146, -138, -137, -145, -142, -142, -152, -159, -170, -170, -174, -177, -182, -181, -168, -174, -184, -188, -194, -197, -192, -192, -181, -190, -196, -189, -185, -187, -194, -191, -190, -186, -190, -189, -199, -196, -190, -196, -202, -201, -205, -203, -198, -199, -205, -213, -209, -210, -209, -211, -201, -197, -202, -202, -203, -205, -204, -211, -209, -211, -214, -209, -202, -201, -192, -197, -185, -184, -178, -175, -146, -149, -142, -129, -114, -121, -138, -153, -162, -135, -138, -137, -159, -182, -206, -229, -229, -233, -234, -235, -224, -228, -214, -220, -226, -226, -231, -256, -273, -262, -186, -53, 135, 337, 479, 496, 383, 141, -182, -361, -380, -375, -381, -394, -399, -369, -322, -288, -260, -249, -243, -241, -233, -239, -225, -220, -209, -205, -198, -198, -193, -192, -188, -185, -173, -173, -176, -174, -167, -166, -174, -171, -165, -159, -161, -146, -143, -139, -136, -132, -133, -121, -118, -112, -110, -110, -106, -102, -101, -93, -97, -102, -98, -86, -90, -90, -87, -81, -85, -80, -77, -79, -69, -70, -63, -55, -63, -57, -47, -42, -42, -43, -45, -40, -47, -43, -41, -42, -45, -45, -48, -55, -53, -65, -67, -69, -70, -77, -82, -85, -86, -90, -90, -77, -85, -86, -85, -91, -93, -89, -86, -90, -93, -87, -89, -85, -77, -86, -86, -82, -85, -82, -83, -89, -87, -91, -89, -83, -82, -85, -79, -78, -71, -77, -79, -82, -78, -72, -74, -74, -79, -79, -81, -78, -74, -66, -69, -61, -63, -53, -53, -40, -39, -27, -9, 3, 37, 35, 15, -9, -33, -37, -15, -13, -13, -31, -65, -81, -89, -81, -77, -86, -86, -82, -89, -85, -72, -80, -75, -105, -121, -114, -66, 42, 215, 403, 563, 615, 527, 319, 27, -194, -259, -261, -256, -269, -279, -262, -216, -167, -130, -105, -98, -91, -89, -85, -73, -74, -70, -65, -49, -53, -61, -57, -55, -48, -41, -23, -35, -27, -18, -15, -29, -22, -9, -7, -7, -21, -17, 1, 11, 17, 25, 23, 25, 43, 42, 53, 55, 55, 53, 61, 65, 51, 57, 65, 47, 51, 49, 49, 64, 73, 73, 75, 71, 75, 72, 79, 91, 99, 88, 91, 95, 105, 110, 113, 107, 99, 97, 97, 95, 95, 97, 101, 97, 83, 87, 75, 64, 63, 77, 69, 51, 49, 45, 57, 57, 67, 72, 56, 57, 57, 55, 57, 51, 45, 55, 59, 61, 59, 43, 35, 47, 51, 47, 51, 54, 49, 59, 47, 33, 24, 34, 40, 59, 63, 55, 49, 33, 39, 45, 43, 47, 49, 49, 57, 51, 66, 67, 72, 65, 63, 83, 81, 99, 112, 111, 131, 143, 152, 147, 125, 111, 117, 136, 115, 101, 93, 59, 51, 47, 46, 57, 45, 25, 31, 35, 25, 48, 51, 33, 31, 19, 27, 97, 203, 341, 501, 646, 678, 555, 340, 59, -81, -141, -210, -223, -195, -167, -174, -150, -101, -61, -23, 3, 11, 9, 9, 35, 43, 64, 64, 49, 49, 63, 73, 85, 97, 75, 65, 71, 68, 67, 85, 96, 97, 98, 105, 93, 88, 103, 107, 109, 129, 131, 137, 137, 129, 125, 140, 140, 131, 145, 139, 133, 145, 145, 144, 153, 152, 135, 139, 141, 137, 145, 155, 156, 145, 145, 157, 161, 161, 171, 169, 146, 151, 151, 159, 153, 152, 151, 147, 151, 145, 144, 143, 139, 143, 141, 129, 128, 129, 139, 132, 136, 128, 125, 123, 127, 130, 129, 131, 132, 133, 138, 131, 125, 126, 125, 123, 120, 119, 120, 123, 113, 109, 93, 97, 99, 104, 97, 89, 91, 94, 89, 95, 91, 93, 95, 93, 97, 89, 91, 91, 91, 91, 91, 91, 95, 101, 113, 115, 120, 129, 143, 141, 159, 173, 176, 185, 188, 187, 176, 155, 147, 144, 167, 165, 140, 115, 97, 77, 51, 57, 61, 62, 67, 69, 73, 75, 73, 73, 69, 53, 23, 31, 83, 201, 365, 559, 714, 764, 688, 459, 139, -70, -117, -117, -138, -152, -152, -125, -69, -26, 11, 33, 48, 47, 54, 57, 67, 91, 88, 87, 89, 97, 104, 93, 98, 95, 105, 104, 113, 129, 123, 127, 131, 129, 128, 139, 137, 139, 145, 157, 155, 161, 167, 169, 175, 181, 175, 179, 183, 183, 187, 187, 189, 187, 184, 192, 194, 193, 195, 189, 189, 195, 189, 192, 193, 196, 203, 197, 195, 202, 195, 195, 204, 211, 216, 211, 217, 211, 207, 203, 203, 209, 200, 195, 193, 187, 184, 179, 179, 177, 179, 170, 173, 165, 168, 171, 166, 165, 165, 161, 155, 161, 159, 162, 159, 149, 153, 145, 147, 151, 153, 152, 145, 149, 149, 145, 139, 139, 145, 138, 139, 144, 141, 143, 139, 137, 137, 137, 128, 125, 125, 127, 129, 133, 123, 123, 123, 125, 129, 135, 139, 145, 157, 161, 176, 183, 179, 200, 205, 197, 179, 167, 145, 163, 177, 170, 157, 120, 95, 78, 71, 63, 65, 73, 69, 72, 73, 77, 75, 71, 59, 42, 13, 19, 77, 198, 361, 560, 737, 800, 735, 527, 184, -65, -113, -105, -101, -122, -129, -96, -55, -16, 16, 27, 35, 29, 35, 51, 55, 56, 65, 64, 67, 71, 67, 77, 73, 73, 85, 85, 91, 91, 91, 95, 107, 113, 113, 115, 114, 119, 122, 129, 130, 129, 129, 139, 139, 141, 147, 149, 153, 153, 151, 155, 153, 157, 165, 161, 173, 170, 187, 189, 203, 201, 200, 205, 216, 229, 231, 234, 251, 257, 241, 189, 170, 172, 178, 187, 202, 195, 192, 184, 173, 169, 161, 153, 137, 135, 125, 120, 115, 113, 103, 114, 111, 104, 112, 101, 104, 105, 107, 105, 99, 81, 50, 49, 125, 253, 441, 661, 842, 909, 781, 507, 175, -25, -87, -109, -133, -150, -139, -103, -59, -27, -1, 30, 33, 27, 24, 35, 31, 35, 41, 55, 53, 48, 53, 47, 47, 53, 51, 56, 54, 61, 57, 59, 63, 57, 73, 70, 75, 75, 83, 91, 91, 103, 111, 106, 112, 117, 121, 115, 117, 121, 127, 137, 147, 145, 137, 136, 147, 147, 147, 159, 153, 161, 169, 179, 177, 181, 201, 207, 205, 207, 192, 185, 179, 183, 189, 191, 178, 167, 149, 147, 149, 141, 131, 137, 131, 115, 107, 101, 102, 107, 106, 105, 103, 99, 101, 93, 89, 85, 88, 75, 79, 85, 90, 91, 83, 85, 75, 69, 73, 67, 77, 81, 70, 71, 73, 73, 65, 63, 61, 69, 75, 73, 71, 64, 67, 63, 59, 64, 81, 72, 67, 69, 77, 81, 69, 64, 49, 51, 61, 72, 61, 67, 61, 59, 71, 61, 55, 55, 64, 59, 67, 71, 67, 67, 71, 67, 59, 53, 51, 50, 69, 69, 67, 75, 75, 85, 95, 99, 113, 129, 149, 153, 147, 133, 117, 111, 131, 137, 131, 111, 63, 39, 27, 29, 23, 15, 19, 25, 21, 35, 35, 37, 29, 5, -7, 7, 59, 143, 265, 429, 584, 662, 583, 395, 111, -95, -148, -198, -241, -230, -227, -220, -192, -133, -97, -65, -45, -34, -17, -7, 13, 25, 21, 15, 17, 7, 7, 13, 15, 23, 33, 35, 45, 51, 49, 47, 47, 43, 49, 51, 45, 45, 65, 83, 99, 85, 61, 49, 65, 70, 69, 75, 86, 83, 81, 77, 83, 77, 89, 97, 97, 95, 91, 82, 80, 91, 105, 95, 97, 99, 99, 107, 105, 107, 119, 130, 129, 135, 145, 155, 161, 165, 173, 167, 173, 183, 203, 225, 221, 237, 243, 211, 171, 143, 101, 97, 107, 114, 137, 125, 111, 115, 120, 109, 107, 90, 75, 71, 58, 67, 67, 75, 83, 67, 64, 55, 55, 59, 35, 1, -13, 1, 57, 169, 328, 520, 668, 685, 580, 365, 67, -125, -197, -219, -236, -231, -206, -176, -134, -89, -47, -35, -7, 3, -5, -7, -2, 9, 5, 2, -1, 3, 11, 19, 35, 19, 15, 9, 11, 27, 31, 33, 24, 25, 30, 43, 59, 54, 35, 31, 33, 37, 47, 51, 49, 43, 43, 33, 35, 39, 39, 29, 37, 43, 35, 37, 35, 35, 35, 43, 43, 47, 51, 59, 57, 59, 55, 65, 73, 73, 82, 83, 83, 87, 91, 85, 90, 90, 83, 81, 85, 86, 81, 77, 73, 73, 61, 71, 67, 66, 64, 71, 67, 67, 58, 61, 61, 59, 57, 57, 57, 49, 55, 51, 51, 45, 51, 51, 47, 45, 40, 37, 41, 51, 47, 41, 35, 41, 39, 31, 37, 43, 39, 31, 31, 27, 25, 29, 31, 23, 17, 15, 29, 25, 29, 35, 35, 35, 27, 35, 33, 33, 32, 29, 24, 33, 33, 31, 43, 27, 35, 35, 31, 37, 43, 40, 41, 51, 47, 55, 67, 66, 81, 91, 99, 107, 117, 130, 131, 113, 107, 90, 83, 99, 97, 95, 66, 49, 24, 3, 7, 3, 11, 11, 17, 19, 13, 22, 5, -15, -41, -59, -25, 61, 213, 402, 607, 732, 725, 613, 326, -5, -150, -141, -141, -162, -185, -164, -124, -75, -47, -37, -21, -17, -23, -11, -15, -3, 3, 10, 9, 17, 11, 11, 18, 19, 23, 23, 23, 23, 22, 31, 37, 41, 46, 47, 51, 44, 49, 55, 48, 64, 70, 67, 65, 73, 75, 83, 87, 85, 81, 87, 89, 88, 91, 99, 97, 93, 95, 100, 108, 105, 106, 104, 107, 103, 105, 105, 109, 107, 109, 122, 126, 123, 125, 135, 142, 137, 139, 144, 144, 147, 141, 139, 135, 128, 131, 123, 117, 115, 107, 109, 104, 103, 89, 91, 89, 87, 90, 91, 85, 85, 80, 85, 75, 57, 53, 56, 58, 51, 51, 49};
//dataType SIGNAL2[] = {49, 43, 43, 43, 37, 41, 39, 37, 39, 34, 37, 31, 31, 30, 25, 23, 25, 25, 19, 19, 21, 26, 17, 19, 22, 19, 19, 11, 11, 17, 18, 19, 15, 15, 7, 15, 19, 23, 23, 27, 39, 41, 53, 65, 72, 81, 79, 93, 97, 87, 67, 51, 53, 72, 75, 64, 34, 13, -5, -25, -37, -39, -35, -37, -31, -33, -29, -35, -40, -53, -77, -92, -97, -21, 83, 249, 441, 613, 678, 609, 413, 82, -166, -210, -193, -193, -218, -221, -195, -154, -121, -84, -65, -61, -61, -59, -54, -42, -37, -43, -30, -21, -29, -25, -29, -18, -19, -15, -13, -3, -3, 1, 10, 11, 11, 5, 9, 17, 17, 22, 21, 23, 29, 35, 35, 41, 31, 35, 33, 39, 40, 37, 43, 43, 43, 35, 37, 41, 51, 51, 41, 45, 41, 49, 43, 31, 39, 43, 47, 47, 55, 59, 59, 63, 61, 69, 75, 79, 81, 72, 78, 75, 73, 73, 67, 69, 65, 63, 54, 51, 45, 53, 43, 43, 41, 55, 47, 47, 47, 37, 27, 19, 22, 23, 21, 21, 17, 19, 17, 11, 18, 9, 7, 5, -1, 1, 2, -1, 1, 3, 1, -1, -2, -10, -5, -13, -7, -9, -13, -11, -9, -13, -7, -7, -13, -13, -15, -15, -13, -15, -15, -21, -15, -23, -13, -14, -9, -11, -1, 19, 17, 25, 35, 56, 65, 81, 83, 67, 61, 47, 40, 43, 61, 45, 25, -5, -33, -49, -55, -57, -63, -53, -45, -48, -48, -47, -53, -56, -69, -90, -108, -86, 1, 151, 340, 526, 633, 612, 465, 179, -111, -237, -232, -217, -238, -248, -248, -204, -165, -133, -102, -87, -77, -85, -81, -77, -71, -71, -77, -63, -59, -53, -50, -43, -57, -57, -49, -47, -45, -25, -25, -29, -29, -24, -7, -1, -3, -1, -9, -18, -17, -7, -5, 1, 9, 13, 9, 11, 9, 19, 31, 33, 31, 17, 11, 23, 33, 31, 27, 33, 31, 35, 35, 27, 35, 40, 29, 45, 47, 53, 59, 51, 56, 57, 67, 65, 65, 64, 57, 55, 54, 53, 59, 61, 59, 57, 51, 43, 41, 40, 45, 37, 31, 17, 3, 7, 19, 11, 15, 19, 9, 11, 1, -6, 1, 3, 5, -1, 15, -1, -19, -17, -13, -15, -21, -25, -29, -23, -23, -22, -13, -13, -23, -23, -25, -25, -37, -47, -47, -31, -29, -30, -43, -52, -50, -43, -27, -21, -30, -43, -55, -57, -50, -29, -29, -10, -24, -37, -33, -22, -10, 3, 16, 27, 39, 57, 55, 11, -5, 19, 23, 15, -15, -33, -55, -85, -101, -97, -90, -109, -101, -85, -81, -82, -73, -82, -123, -129, -118, -69, 23, 166, 330, 473, 531, 443, 249, -7, -206, -274, -313, -376, -384, -373, -342, -313, -254, -214, -193, -187, -186, -177, -158, -146, -138, -130, -130, -150, -153, -153, -142, -148, -122, -97, -98, -134, -159, -138, -130, -126, -117, -106, -130, -117, -117, -125, -113, -113, -134, -125, -123, -110, -89, -99, -117, -129, -140, -137, -119, -105, -98, -120, -141, -138, -118, -127, -137, -138, -150, -143, -114, -101, -121, -147, -150, -136, -116, -129, -130, -133, -133, -121, -117, -128, -122, -106, -116, -138, -149, -158, -153, -141, -137, -141, -118, -130, -155, -172, -174, -174, -178, -173, -166, -161, -144, -149, -175, -190, -190, -187, -192, -214, -219, -215, -214, -216, -230, -232, -230, -234, -241, -233, -241, -242, -241, -244, -250, -246, -254, -254, -263, -258, -256, -259, -264, -255, -269, -268, -276, -270, -273, -274, -273, -273, -278, -270, -274, -282, -273, -274, -265, -257, -249, -252, -247, -237, -218, -214, -203, -201, -198, -223, -240, -251, -237, -231, -237, -247, -276, -310, -327, -340, -340, -345, -340, -345, -337, -343, -331, -345, -346, -369, -386, -387, -363, -279, -136, 43, 213, 313, 275, 123, -149, -420, -531, -559, -565, -583, -589, -583, -545, -501, -455, -435, -419, -412, -402, -399, -394, -386, -379, -381, -371, -376, -373, -369, -369, -356, -363, -361, -367, -359, -352, -349, -349, -359, -359, -363, -361, -357, -355, -347, -340, -346, -340, -329, -330, -336, -332, -325, -325, -328, -323, -329, -331, -325, -329, -320, -324, -325, -325, -326, -326, -324, -327, -326, -334, -332, -324, -321, -319, -318, -317, -318, -310, -312, -313, -311, -308, -313, -316, -317, -312, -323, -320, -322, -327, -328, -330, -333, -339, -338, -344, -345, -347, -340, -347, -352, -351, -351, -345, -351, -360, -358, -358, -360, -363, -369, -377, -375, -373, -379, -377, -382, -383, -385, -381, -383, -387, -385, -385, -386, -382, -383, -382, -385, -384, -383, -381, -391, -397, -393, -385, -390, -392, -392, -391, -392, -394, -392, -383, -375, -367, -375, -363, -347, -347, -340, -332, -311, -316, -301, -301, -305, -330, -340, -347, -322, -315, -324, -353, -375, -399, -413, -412, -419, -411, -415, -402, -405, -403, -401, -403, -397, -431, -455, -467, -431, -337, -194, 9, 200, 304, 269, 115, -176, -469, -585, -573, -561, -581, -595, -573, -535, -487, -467, -448, -433, -439, -431, -427, -413, -409, -407, -397, -403, -401, -397, -393, -385, -385, -385, -388, -375, -378, -371, -375, -363, -361, -361, -352, -361, -360, -356, -352, -347, -347, -342, -337, -337, -331, -330, -333, -337, -321, -320, -324, -314, -313, -310, -308, -307, -314, -311, -306, -299, -304, -302, -299, -293, -290, -288, -288, -281, -270, -269, -268, -258, -259, -256, -260, -256, -259, -251, -254, -253, -248, -254, -261, -266, -269, -273, -270, -273, -279, -281, -285, -290, -287, -279, -282, -284, -280, -269, -277, -268, -265, -268, -266, -270, -260, -269, -273, -270, -273, -274, -264, -265, -274, -270, -270, -274, -269, -270, -261, -261, -258, -266, -267, -266, -266, -262, -264, -265, -267, -260, -269, -266, -261, -266, -264, -258, -256, -250, -258, -258, -256, -247, -239, -241, -230, -223, -210, -202, -185, -185, -165, -155, -160, -177, -197, -205, -189, -173, -184, -209, -232, -249, -268, -265, -273, -262, -265, -261, -262, -254, -252, -256, -250, -262, -288, -286, -257, -162, -9, 181, 367, 455, 406, 227, -65, -330, -425, -409, -413, -415, -427, -422, -371, -333, -296, -281, -263, -264, -263, -251, -241, -239, -239, -236, -222, -227, -220, -220, -218, -206, -196, -197, -191, -196, -194, -189, -182, -179, -168, -157, -153, -149, -146, -150, -134, -130, -129, -133, -133, -121, -119, -115, -128, -129, -114, -107, -117, -109, -99, -95, -94, -103, -97, -93, -94, -96, -90, -91, -87, -86, -78, -75, -77, -69, -63, -45, -57, -55, -49, -37, -37, -31, -39, -42, -41, -40, -45, -40, -39, -45, -51, -49, -61, -57, -57, -62, -53, -53, -57, -61, -61, -64, -66, -55, -71, -59, -69, -67, -69, -69, -64, -75, -87, -85, -79, -69, -61, -73, -73, -79, -79, -81, -77, -77, -78, -77, -77, -79, -69, -78, -75, -69, -89, -90, -86, -85, -89, -86, -80, -77, -85, -66, -50, -66, -82, -88, -80, -76, -73, -53, -37, -33, -37, -15, 9, 8, 19, 31, 25, 7, -7, 1, 14, 41, 17, -19, -43, -77, -82, -83, -85, -94, -88, -79, -69, -73, -69, -83, -73, -80, -90, -63, -21, 59, 175, 337, 490, 580, 528, 344, 73, -153, -243, -287, -306, -315, -286, -282, -251, -212, -169, -147, -110, -99, -101, -106, -85, -85, -70, -61, -69, -73, -47, -41, -57, -41, -37, -33, -15, -19, -35, -41, -47, -41, -29, -17, -5, 8, 5, -5, -21, -13, 1, 19, 31, 45, 19, -8, -7, -5, -1, 9, 25, 21, 11, 3, 5, 19, 19, 15, 11, 11, 11, 25, 51, 45, 31, 25, 39, 43, 41, 35, 40, 42, 51, 47, 59, 49, 37, 38, 51, 50, 53, 59, 45, 41, 29, 37, 39, 31, 27, 19, 25, 29, 27, 31, 27, 23, 25, 27, 27, 19, 18, 19, 19, 15, 11, 11, 9, 8, 5, 2, 7, 3, -5, -1, -3, -8, -9, -8, -1, -2, -6, -7, -9, -11, -7, -5, -13, -11, -8, -9, -5, -7, -3, -5, -3, 0, 1, 1, 3, 10, 15, 5, 9, 23, 37, 41, 54, 64, 75, 85, 91, 93, 88, 72, 59, 45, 61, 79, 67, 57, 27, 7, -15, -37, -22, -25, -21, -25, -17, -18, -13, -13, -14, -15, -45, -53, -35, 45, 165, 328, 523, 665, 696, 589, 367, 49, -147, -194, -210, -210, -223, -224, -205, -145, -101, -65, -53, -37, -31, -23, -21, -15, -5, -5, 1, 21, 17, 9, 10, 24, 19, 29, 31, 35, 39, 42, 48, 49, 49, 41, 53, 57, 57, 63, 71, 67, 75, 88, 89, 87, 95, 99, 95, 99, 109, 107, 104, 111, 115, 121, 121, 120, 121, 125, 120, 123, 125, 123, 127, 125, 127, 136, 133, 131, 134, 143, 147, 147, 147, 143, 143, 139, 141, 143, 147, 145, 147, 141, 137, 141, 138, 141, 133, 133, 129, 122, 121, 123, 122, 121, 126, 119, 125, 126, 128, 117, 120, 115, 115, 115, 114, 115, 115, 104, 110, 107, 105, 106, 103, 97, 107, 95, 99, 99, 95, 91, 95, 99, 95, 99, 95, 95, 97, 89, 91, 93, 93, 89, 91, 89, 93, 93, 95, 99, 109, 111, 113, 123, 131, 136, 143, 149, 164, 169, 185, 193, 186, 171, 141, 145, 161, 179, 169, 141, 114, 97, 70, 59, 64, 69, 65, 75, 75, 70, 69, 75, 75, 63, 35, 15, 39, 120, 275, 465, 667, 792, 798, 650, 363, 41, -109, -98, -91, -108, -127, -121, -69, -23, 11, 29, 41, 45, 53, 65, 75, 83, 89, 95, 89, 97, 95, 101, 105, 105, 113, 119, 127, 127, 133, 133, 136, 139, 145, 151, 155, 155, 167, 157, 158, 173, 173, 181, 187, 191, 197, 197, 207, 204, 200, 209, 208, 211, 207, 216, 222, 231, 233, 243, 257, 265, 267, 271, 265, 287, 296, 299, 313, 317, 305, 262, 236, 227, 234, 249, 253, 253, 253, 239, 225, 229, 223, 209, 197, 192, 181, 185, 182, 183, 187, 185, 179, 171, 174, 179, 177, 180, 179, 175, 176, 165, 145, 136, 189, 304, 488, 709, 905, 1014, 929, 675, 331, 83, -5, -13, -37, -45, -49, -17, 31, 67, 96, 105, 121, 129, 128, 131, 137, 134, 135, 136, 143, 149, 151, 153, 153, 157, 159, 163, 171, 171, 167, 175, 173, 187, 189, 185, 185, 185, 185, 197, 193, 209, 210, 211, 210, 213, 221, 226, 229, 225, 227, 229, 241, 242, 236, 235, 245, 256, 272, 281, 285, 277, 278, 289, 296, 299, 306, 310, 308, 307, 305, 306, 318, 313, 303, 285, 287, 287, 285, 271, 267, 256, 255, 247, 233, 231, 221, 222, 219, 227, 219, 217, 209, 201, 207, 207, 223, 223, 221, 208, 205, 198, 205, 215, 195, 190, 185, 192, 197, 197, 195, 193, 194, 193, 187, 181, 171, 185, 193, 189, 195, 203, 196, 177, 171, 163, 171, 166, 171, 173, 171, 163, 171, 183, 179, 173, 168, 179, 165, 166, 167, 163, 187, 191, 163, 153, 153, 165, 171, 177, 184, 185, 179, 177, 157, 163, 171, 177, 175, 192, 195, 195, 199, 189, 199, 211, 219, 235, 265, 289, 277, 241, 207, 200, 205, 223, 241, 217, 189, 155, 135, 137, 123, 107, 125, 141, 127, 121, 125, 141, 133, 127, 85, 89, 143, 223, 322, 461, 628, 743, 722, 580, 315, 67, -27, -77, -137, -168, -148, -125, -93, -56, -25, -3, 11, 59, 99, 97, 65, 59, 74, 95, 94, 88, 83, 81, 91, 81, 89, 99, 101, 127, 155, 130, 95, 90, 89, 97, 109, 129, 152, 152, 136, 117, 104, 115, 127, 152, 163, 151, 123, 113, 123, 127, 134, 139, 125, 118, 139, 149, 139, 121, 119, 123, 131, 152, 157, 143, 139, 145, 159, 179, 169, 145, 139, 145, 159, 167, 163, 163, 157, 155, 165, 161, 159, 157, 147, 149, 152, 154, 143, 131, 139, 136, 135, 131, 133, 131, 133, 129, 125, 127, 115, 107, 109, 109, 103, 99, 99, 95, 97, 90, 83, 81, 79, 78, 75, 67, 64, 56, 55, 59, 51};
dataType *SIGNAL;
uint32_t index_offset = 0;

int next_out = 0;
int *Rs;//[RS_LEN] = {0};
int RS_LEN = 32;
int next_r = 0;
uint32_t SIGNAL_LEN = 0;
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

dataType input(uint32_t index)
{
    int num = NOSAMPLE;
    if (index < SIGNAL_LEN)
    {
        num = SIGNAL[index];
    }
    return num;
}

void output(int out)
{
	if(out)
	{
		Rs[next_r++] = next_out + index_offset;
#ifdef OVERWRITE_RS_ON_OVERFLOW //If we are fine with Rs being overwritten as the array overflows
		next_r = next_r % RS_LEN; //prevent bad memory access
#endif
	}
	next_out++;
	//SIGNAL_OUT[next_out++] = out;
}
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

void PanTompkins_init(dataType *signal, uint32_t sig_len, int *Rs_out, int in_RS_LEN)
{
	//Direct PT to the correct input source
	SIGNAL = signal;
	SIGNAL_LEN = sig_len;
    Rs = Rs_out;
    RS_LEN = in_RS_LEN;
	next_out = 0;
	sample = 0;
	next_r = 0;
	index_offset = 0;
	for(int i = 0; i < RS_LEN; i++)
	{
        //Set the input to an invalid flag
        //Useful to determine whether the Rs
        //have been written up to this point
        //or not
		Rs[i] = -1; 
	}
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
}
void PanTompkins_sigswap(dataType *new_signal, uint32_t new_sig_len, uint32_t new_index_offset, uint32_t r_index_start, int *new_Rs, int new_RS_LEN)
{
	SIGNAL = new_signal;
	SIGNAL_LEN = new_sig_len;
    Rs = new_Rs;
    RS_LEN = new_RS_LEN;
	index_offset = new_index_offset;
	if(r_index_start >= 0)
	{
		next_r = r_index_start;
	}
	else
	{
		next_r = 0; //start at beginning
	}
}


void PanTompkins()
{

    dataType last_avg[MOVING_AVG_LEN];
    dataType avg = 0;
    memset(last_avg,0,MOVING_AVG_LEN*sizeof(dataType));
	next_out = 0;
	sample = 0;
	//next_r = 0;


	// i and j are iterators for loops.
	// sample counts how many samples have been read so far.
	// lastQRS stores which was the last sample read when the last R sample was triggered.
	// lastSlope stores the value of the squared slope when the last R sample was triggered.
	// currentSlope helps calculate the max. square slope for the present sample.
	// These are all uint32_t so that very long signals can be read without messing the count.
	uint32_t i, j, lastQRS = 0, lastSlope = 0, currentSlope = 0;

    // qrs tells whether there was a detection or not.
    // regular tells whether the heart pace is regular or not.
    // prevRegular tells whether the heart beat was regular before the newest RR-interval was calculated.
    bool qrs, regular = true, prevRegular;


    // The main loop where everything proposed in the paper happens. Ends when there are no more signal samples.
    do{
        // Test if the buffers are full.
        // If they are, shift them, discarding the oldest sample and adding the new one at the end.
        // Else, just put the newest sample in the next free position.
        // Update 'current' so that the program knows where's the newest sample.
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
        signal[current] = input(sample);
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
            break;
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
                outputSignal[current] = qrs;
                if (sample > DELAY + BUFFSIZE)
                    output(outputSignal[0]);
                continue;
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

                if (qrs)
                {
                    outputSignal[current] = false;
                    outputSignal[i] = true;
                    if (sample > DELAY + BUFFSIZE)
                        output(outputSignal[0]);
                    continue;
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
        // The current implementation outputs '0' for every sample where no peak was detected,
        // and '1' for every sample where a peak was detected. It should be changed to fit
        // the desired application.
        // The 'if' accounts for the delay introduced by the filters: we only start outputting after the delay.
        // However, it updates a few samples back from the buffer. The reason is that if we update the detection
        // for the current sample, we might miss a peak that could've been found later by backsearching using
        // lighter thresholds. The final waveform output does match the original signal, though.
        outputSignal[current] = qrs;
        if (sample > DELAY + BUFFSIZE)
            output(outputSignal[0]);
    } while (signal[current] != NOSAMPLE);

    // Output the last remaining samples on the buffer
    for (i = 1; i < BUFFSIZE; i++)
        output(outputSignal[i]);

}
#else
//If we're not using PT functionality, use as little memory as possible, just allow it to compile
void panTompkins(){}
void PanTompkins_init(dataType *signal, uint32_t sig_len, int *Rs_out, int in_RS_LEN){}
void PanTompkins_sigswap(dataType *new_signal, uint32_t new_sig_len, uint32_t new_index_offset, uint32_t r_index_start, int *new_Rs, int new_RS_LEN){}
void PanTompkins_setsaved_filterstate(FilterState *fs);
void PanTompkins_exportsaved_filterstate(FilterState *fs);
void PanTompkins_save_filterstate(){}
void PanTompkins_load_filterstate(){}
#endif
