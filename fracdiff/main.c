//
//  main.c
//  fracdiff

// code for exploring fractional differencing and fractional integration

// diffent.com

// inspired by part of this python notebook:
// https://github.com/Dhiraj96/blog_posts/blob/master/post_1%2C2%2C3/Post1_Fractional_Differencing.ipynb

// ... 3 fractional differencing blog posts from the same author:
// http://kidquant.blogspot.com

// ... and this nice article, which is also accessible to the less-than-PhD level reader,
// but you do need some math background to understand
// https://towardsdatascience.com/preserving-memory-in-stationary-time-series-6842f7581800

// also see gretl econometrics package fracdiff function for reference
// http://gretl.sourceforge.net/gretl-help/funcref.html#fracdiff

// the gretl fracdiff seems to return 1 less value in its output series versus its input series,
// unlike this routine which returns the same number of points in its output as in its input

// For a slightly different take on fractional time series modeling methods, see the Fractional Empirical Motion
// feature of our MCarloRisk3D iPad app, where the Hurst exponent can be specified for projecting
// forward price / probability envelopes.  FEM is a generalization of FBM (Fractional Brownian Motion)
// which does not assume a Gaussian for the underlying distribution, but rather allows
// an empirical distribution to be used:

// https://apps.apple.com/us/app/mcarlorisk3d/id641208540

// In particular, see slide set 9 in this training guide for examples, and how to access the Hurst parameter:
// https://diffent.com/mcrtrain/mcrtrain.html

// Notes:

// Ordinary differencing of a time series can cause an analyst to lose much long-term
// memory of what has gone on in the past in that time series, since differencing
// is a type of high pass filter operation, with lower frequencies (slow changes over time) in the series
// being attenuated.  In fact, consulting a Fourier transform table, we are reminded
// that simple differentiation truncates low frequencies linearly w/ increasing frequency:

// see row 6:  d^n(f(t))/dt^n of this somewhat easy to read (but not necessarily to understand) Fourier transform table:
// https://ethz.ch/content/dam/ethz/special-interest/baug/ibk/structural-mechanics-dam/education/identmeth/fourier.pdf

// note:  lower case Greek w (omega) is commonly used as "frequency" in these transform tables.

// e.g. 0 frequency (w = 0) ("direct current," DC, or constant value in the case of financial data: e.g. long term avg price level)
// is fully filtered out by the derivative operation (or difference operation for discrete time),
// very low frequencies are mostly attenuated, and, as the frequency increases,
// more of the signal is kept.  Higher frequency data (possibly noise?) is amplified (which is not necessarily
// a good thing).

// The presence of j (the imaginary unit) in the frequency domain (right side of the table) suggests that a phase shift occurs (which
// may be thought of as a shift in where in time the peaks are on the filtered time series versus
// the time-positioned peaks on the original series)
// when we apply (fractional) differencing, which shift is common with backwards-looking time series operations.  Think
// 20 day moving average:  For an ordinary backward looking 20 day moving average on an asset price series,
// a peak in the actual data occurs, then this
// peak is gradually included in the height of the average chart as you move forward in time, with the peak
// of the moving average curve being much later than when the peak actually occured in the raw price data.

// Moreover, we see from line 8 of the above transform table that definite integration (the inverse of our differencing)
// suggests a filter of type 1/w (e.g. 1/frequency), ignoring for the moment the complex j phase shift.

// This makes sense, since we know from elementary calculus that integration and differentiation
// (or summing and differencing in computer-land discrete time) are inverse operations in the time
// domain, so it is natural that w and 1/w are (multiplicative inverses) in the frequency domain.

// Without proof here, and setting aside the complex unit j for a moment:
// for fractional differencing, the equivalent filter transform in frequency domain
// seems to be:  w^d where d is the differencing parameter; a simple power function of frequency
// (the same as integer differencing but allowing a fractional power).
// Note that this works for both fractional differencing and fractional integration, since w^(-d) = 1/(w^d),
// and note again that w^(-d) * w^d = 1 (e.g. fractional integration and fractional differencing to
// the same fraction "d" are inverse operations as we would expect).

// In the frequency domain, we can multiply the transform of differencing (w) by
// the transform of integration (1/w) and note that this yields 1, the "pass thru" filter.
// Happily (and not accidentally), the j is in the numerator for differencing, and in the denominator
// for integrating, and so the imaginary j units cancel out, as we hope they should for a true inverse operation.
// E.g. the phase shifts of differencing ane re-integrating cancel each other out when those 2 inverse operations
// are applied sequentially to the same data.

// Various authors have proposed fractionally differencing (e.g. price) time series with a non-integer differencing parameter
// typically less than 1 and set to a value such that the resulting fractionally differenced series is just about stationary
// (see above links)
// rather than doing the harsh integer differencing (ordinary P(t) - P(t-1) differencing) which causes
// more loss or filtering of long term behavior of the series.

// Once a time series is made stationary, the analyst's toolbox opens up as to what kind of analysis
// can be properly performed, as more model types work well on stationary time series
// than on non-stationary time series.  By "work well," we mean "give better predictions."
// The hazard is that if you apply a model that is expecting stationary data to a non stationary data set,
// you may get a forecast... the model may run.  But such forecast will likely be inaccurate.
// In other words, you will get a misleading forecast, which may be worse than no forecast at all.

// The inverse fractional (difference/integral) operation can then be performed on any synthetic
// data generated from fractionally differenced data to put the values back into the original
// problem units (e.g. price).  In other words, the fractional differencing operation is reversable (up
// to a constant in theory, but fully reversable in practice [disregarding a bit of numeric noise
// that may occur due to limited precision operations as in float], in this implementation.

// To invert the fracdiff operation and recover the original series, just reverse the sign on the fractional param d.
// An example of this is included in the main() routine below.

// While some things in C code are obscure, this can also be said for Python,
// so I tried to write this fractional difference code using very basic operations
// for clarity.

// This C code is easily translatable to a variety of other languages as well, since
// it is low level.

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// uncomment to turn on printf statements for testing
//#define printf(...)

// orig built under xcode 10.1 on MacOS, but MacOS-specific routines
// are not used, so this main.c should build on a variety of systems
// w/ a standard C compiler

// on a linux type system or mac, you can also compile and run this routine like this:

// gcc main.c
// ./a.out (default executable name is a.out from gcc compiler)

// On Mac, the gcc command didnt need the standard c library to be added to the gcc command,
// but you may have to specify a -L and/or -l argument to add a system library or two
// if you are on some type of linux or windows system.  However, if you regularly
// compile C code using a command line, you probably know this already.

// ----
// This function allocs the weight array (returned value),
// caller must free() the returned pointer.

// written using single precision floats for speed

// This weights routine is almost verbatim from the python,
// but array operators don't exist in C directly,
// so these are done manually here.

// see notes below for argument descriptions

float * findWeights_ffd(float d, int length, float threshold, int useNWeights) {
    
    float * w = calloc(length, sizeof(float)); // allocate & zero out
    
    w[0] = 1;
    int wcount = 1;
    int k = 1;
    float w_curr = 1;

    // while we still have more weights to process, do the following
    while (k < length) {
        printf("in loop k = %d\n", k);
        // note that wcount-1 is the end of the array's current valid data
        w_curr = (-w[wcount-1]*(d-k+1))/k; // [A] the main recurrence relation for weights
        printf("w_curr = %f\n", w_curr);
        if (fabsf(w_curr) <= threshold) break;
        if (useNWeights > 0 && k >= useNWeights) break;
        
        w[wcount] = w_curr; // add an item to the weights array and
        wcount++;           // and increment where the end of valid data is
        k++;
    }
    return w;

}

// Fractional differencing:
// Each item of the output series is merely a weighted sum of all prior values in the series,
// with weights determined by the recurrence relation at line [A] above.

// series is just an array of floats e.g. prices in the case of econometric / stock data
// len is the length of the input series

// d is the differencing to apply:
// 1 = resolves to ordinary differencing
// e.g. 0.5 = fractional half difference
// e.g. -0.5 = fractional half integrate (undo the related fractional difference of 0.5)

// threshold = weights less than this will be set to zero
// if threhshold is zero, use all weights

// this routine has mainly been tested with threshold set to 0 (keep all weights)

// This is a special case of a more general digital filtering operation in discrete time domain
// e.g. if we use different weights than above, we can filter the input differently.
// However, we are no longer fractionally differencing if we use different weights,
// we are applying more general filter operations, and the frequency domain properties of the (fractional)
// difference such as noted in the aforementioned Fourier transform table no longer apply.

// if useNWeights = 0 use all weights
// if useNWeights > 0 use that number of weights

// If we use fewer than the total number of weights as computed by the weight generation function, it is not a true
// fractional difference, and it is not fully invertible by reversing the difference parameter.
// Similarly, such windowed or truncated fractional differences have different frequency responses
// than the full fractional difference as we might expect:  a general
// characteristic of digital filters:  if you change the weights (e.g. set some to 0), you change the frequency
// characteristics.

float * fracDiff(float * series, int len, float d, float threshold, int useNWeights) {
    
    float * weights = findWeights_ffd(d, len, threshold, useNWeights); // generate the weights
    
    float * df_temp = calloc(len, sizeof(float)); // for output
    
    // for every value in the original series
    for (int i = 0; i < len; i++)
    {
        // dot product of (orig array from i to the end) DOT (full weights array)
        // taking care not to roll past the end of either warray
        
         float sum = 0;

        // go from the given item in the orig series to the end, multiply each
        // of these values by a corresponding weight,
        // the sum up the result
        
        for (int j = i; j < len; j++) {
            sum += series[j] * weights[j-i];  // note that since j starts at i for this loop, we always start with weights[0]
        }
        df_temp[i] = sum;
    }
    
    // Theoretical papers often leave out important points such as this:
    
    // Much like ordinary differencing, once we hit the last element in the input data array,
    // there is no 'next' element to subtract from it.  so what do we do?  what can we do?
    // we are out of data.
    // however, the problem is worse with fractional differencing,
    // because as we index forward in the input array (e.g. backwards in time in this implementation), we continually run out of
    // more and more of the later terms in the fractional difference weighted sum.
    // An interesting thing is that if we just leave this as-is, when re-(fractionally)integrate (with the same
    // limitation of running out of data at the end of the array), the original series is re-generated,
    // perhaps with a bit of numerical noise).  So the run-out-of-data problem is taken into account in
    // the inverse operation too.  A happy accident or purposeful design?  You decide.
    
    // TLDR:  The last value in the input array comes back as itself since there is no data left
    // beyond it to difference with.
    
    free(weights);
                   
    return df_temp;

}

// main program to test the algorithm w/ some default data

int main(int argc, const char * argv[]) {
    
        // test the weight generation routine
        
        int len = 10;

        // the "series" array is the time series we will difference

        // if this is considered a price time series, the most recent price is first
        // and the oldest price is last
    
        // so for this example, an ordinary difference of the first value - the value
        // immediately previous in time is:  series[0]-series[1] = 2-1 = 1
    
        // if you add or delete values from this array, don't forget to change the variable len above correspondingly
    
        float series[] = {2, 1, 3, 5, 6,   0, -1, 2, 2, 5};

        float tolerance = 0; // just use all weights for now, don't set a cutoff tolerance
    
        int useNWeights = 0; // if useNWeights != 0, the operation isn't revers
    
        // difflevel notes:
    
        // if we set this to 1, fd (fractional difference) series becomes the ordinary day to day difference
    
        // if we set this to -1, the result is an ordinary definite integral (cumulative sum, assuming time deltas of 1)
        // for the case of -1 (ordinary definite integral):
        // first value of the resulting array should match the full sum of the above series (computed below for reference).
        // e.g. array element 0 to the end of the array (index len-1)
    
        // the second value of the resulting array should match the sum of the series from array element [1] (2nd element since we are
        // in C code with 0-indexed arrays) to the end of the series
    
        // and so on
    
        // interesting experiment to try:
        // set difflevel close to 0 (e.g. 0.01), you get a slightly modified version of your orig series
        // a very mild high-pass filter
    
        float difflevel = 0.5;
    
        // for ref:  sum of orig series (integral with delta t = 1)
    
        float sumseries = 0;
        for (int i = 0; i < len; i++) sumseries += series[i];
    
        printf("sum of orig series = %f\n", sumseries);
    
        float * w = findWeights_ffd(difflevel, len, tolerance, useNWeights);
    
        // sum the weights to see what they look like in aggregate.
        // this is just for exploration purposes.  do they sum to something
        // interesting like 0, 1, or -1?
    
        float sum = 0;
        for (int i = 0; i < len; i++) {
            printf("w[%d] = %f\n", i, w[i]);
            sum += w[i];
        }
        
        printf("sum of wts = %f\n", sum);
        
        free(w);
        
        // test the fractional differencing on above defined series
    
        // fd = fractional difference
    
        float * fd = fracDiff(series, len, difflevel, tolerance, useNWeights);
        
        // dump out results, compare to orig
        
        for (int i = 0; i < len; i++)
            printf("orig = %f fd = %f\n", series[i], fd[i]);
        
        // re-integrate the differenced series by using -difflevel,
        // should get orig result back maybe w/ a little numeric noise
    
        // fi = fractional integral
    
        float * fi = fracDiff(fd, len, -difflevel, tolerance, useNWeights);
        
        for (int i = 0; i < len; i++)
            printf("redo orig = %f fi = %f\n", series[i], fi[i]);
        
        free(fd);
        free(fi);
    
        return 0;
}
