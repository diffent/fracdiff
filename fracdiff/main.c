//
//  main.c
//  fracdiff

// code for exploring fractional differencing and fractional integration

// diffent.com

// inspired by part of this python notebook:
// https://github.com/Dhiraj96/blog_posts/blob/master/post_1%2C2%2C3/Post1_Fractional_Differencing.ipynb

// ... 3 fractional differencing blog posts from the same author:
// http://kidquant.blogspot.com

// ... and this nice article, which is also accessable to the less-than-PhD level reader
// https://towardsdatascience.com/preserving-memory-in-stationary-time-series-6842f7581800

// ordinary differencing of a time series can cause an analyst to lose some, or, more likely, all, long-term
// memory of what has gone on in the past in that time series, so various authors have
// proposed fractionally differencing (e.g. price) time series with a non-integer differencing parameter
// set to a value such that the resulting fractionally differenced series is jussst about stationary

// once a time series is made stationary, the analyst's toolbox opens up as to what kind of analysis
// can be properly performed

// the inverse fractional (difference/integral) operation can then be performed on any synthetic
// data generated from fractionally differenced data to put the values back into the original
// problem units (e.g. price)

// while some things in C code are obscure, this can also be said for python,
// so i tried to make this fractional difference code using very basic operations
// for clarity

// C code is easily translatable to a variety of other languages as well, since
// it is low level

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// orig built under xcode 10.1 on mac, but mac-specific routines
// are not used, so this main.c should build on a variety of systems
// w/ a standard C compiler

// on a linux type system or mac, you can also compile and run this routine like this:
// gcc main.c
// ./a.out (default executable name is a.out from gcc compiler)

// on mac, the gcc command didnt need the standard c library to be added to the gcc command,
// but you may have to specify a -L and/or -l argument to add a system library or two
// if you are on some type of linux or windows system.  however, if you regularly
// compile C code using a command line, you probably know this already.

// this function allocs weight array (returned value),
// caller must free() the returned value

// written using single precision floats for speed

float * findWeights_ffd(float d, int length, float threshold) {
    
    float * w = calloc(length, sizeof(float)); // zero'd
    
    w[0] = 1;
    int wcount = 1;
    int k = 1;
    float w_curr = 1;
    
    // while we still have more weights to process, do the following
    while (k < length) {
        printf("in loop k = %d\n", k);
        w_curr = (-w[wcount-1]*(d-k+1))/k;
        printf("w_curr = %f\n", w_curr);
        if (fabsf(w_curr) <= threshold) break;
        w[wcount] = w_curr;
        wcount++;
        k++;
    }
    return w;

}

// series is just an array of floats e.g. prices in the case of econometric / stock data
// len is length of series
// d is the differencing to apply:
// 1 = resolves to ordinary differencing
// e.g. 0.5 = fractional difference
// e.g. -0.5 = fractional integrate (undo the related difference of 0.5)
// threshold = weights less than this will be set to zero

float * fracDiff(float * series, int len, float d, float threshold) {
    
    float * weights = findWeights_ffd(d, len, threshold);
    
    float * df_temp = calloc(len, sizeof(float));
    
    // for every value in the original series
    for (int i = 0; i < len; i++)
    {
        // dot product of orig array and a particular subset of the weight array
        
         float sum = 0;

        // go from the given item in the orig series to the end, multiply each
        // of these values from the orig series to the end by a corresponding weight,
        // the sum up the result
        
        for (int j = i; j < len; j++) {
            sum += series[j] * weights[j-i];  // note that since j starts at i for this loop, we always start with weights[0]
        }
        df_temp[i] = sum;
    }
    
    // much like ordinary differencing, once we hit the last element in the input data array,
    // there is no 'next' element to subtract from it.  so what do we do?  what can we do?
    // we are out of data.
    // however, the problem is worse with fractional differencing,
    // because as we move forward in the input array, we continually run out of the later terms
    // in the fractional difference weighted sum.
    // the interesting thing is that if we just leave this as-is, when re-(fractionally)integrate (with the same
    // limitation of running out of data at the end of the array), the original series is re-generated
    // (perhaps with a bit of numerical noise).  so the run-out-of-data problem is taken into account in
    // the inverse operation.
    
    free(weights);
                   
    return df_temp;

}
                    
int main(int argc, const char * argv[]) {
    
        // test the weight generation routine
        
        int len = 20;

       // the "series" array is the time series we will difference

        // if this is considered a price time series, the most recent price is first
        // and the oldest price is last
    
        // so for this example, an ordinary difference of the first value - the value
        // immediately previous is:  series[0]-series[1] = 2-1 = 1
    
        // if you add or delete values from this array, don't forget to change len above correspondingly
    
        float series[] = {2, 1, 3, 5, 6,   0, -1, 2, 2, 5,
                          2, 1, 3, 5, 6,   0, -1, 2, 2, 5};

        float tolerance = 0; // just use all weights for now, don't set a cutoff tolerance
    
        // notes:
        // if we set this to 1, fd (fractional difference) series becomes the ordinary day to day difference
    
        // if we set this to -1, the result is an ordinary definite integral (sum, assuming time deltas of 1)
    
        // first value of the resulting array should match the full sum of the above series (computed below for reference).
        // e.g. array element 0 to the end of the array
    
        // the second value of the resulting array should match the sum of the series from array element [1] (2nd element since we are
        // in C code) to the end of series
    
        // and so on
    
        float difflevel = 0.5;
    
        // for ref:  sum of orig series (integral with delta t = 1)
    
        float sumseries = 0;
        for (int i = 0; i < len; i++) sumseries += series[i];
    
        printf("sum of orig series = %f\n", sumseries);
    
        float * w = findWeights_ffd(difflevel, len, tolerance);
    
        // sum the weights to see what they look like in aggregate.
        // this is just for exploration purposes.  do they sum to something
        // interesting like 0, 1, or -1?  not for small len in the general case.
    
        float sum = 0;
        for (int i = 0; i < len; i++) {
            printf("w[%d] = %f\n", i, w[i]);
            sum += w[i];
        }
        
        printf("sum of wts = %f\n", sum);
        
        free(w);
        
        // test the fractional differencing on above defined series
    
        // fd = fractional difference
    
        float * fd = fracDiff(series, len, difflevel, tolerance);
        
        // dump out results, compare to orig
        
        for (int i = 0; i < len; i++)
            printf("orig = %f fd = %f\n", series[i], fd[i]);
        
        // re-integrate the differenced series by using -difflevel,
        // should get orig result back maybe w/ a little numeric noise
    
        // fi = fractional integral
    
        float * fi = fracDiff(fd, len, -difflevel, tolerance);
        
        for (int i = 0; i < len; i++)
            printf("orig = %f fi = %f\n", series[i], fi[i]);
        
        free(fd);
        free(fi);
    
        return 0;
}
