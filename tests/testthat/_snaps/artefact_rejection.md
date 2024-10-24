# Calculate epoch stats

    Code
      demo_epo
    Output
           measure epoch       A5       A13       A21      A29       A31        B5
            <char> <num>    <num>     <num>     <num>    <num>     <num>     <num>
        1:     max     1 11.56975  5.843069  6.950029 16.05196  5.835857 16.037969
        2:     max     2 11.47310  3.076121 17.336609 25.07866  5.495951 13.935791
        3:     max     3 10.00425  5.781016 17.194586 17.19246  9.782940 12.822427
        4:     max     4  9.31763 10.440170 15.260184 17.20443 13.233307  4.500334
        5:     max     5 13.83158  7.853081 13.628233 21.79226 11.007382 13.648832
       ---                                                                        
      396:  minmax    76 46.45077 15.341318 64.726750 34.25538 45.579252 56.265198
      397:  minmax    77 30.11000 20.520123 29.243788 49.73817 33.340733 40.902041
      398:  minmax    78 18.74668 10.368268 19.936261 31.79093 23.085611 23.651198
      399:  minmax    79 25.99021 18.232140 19.778950 36.60640 24.057504 36.741414
      400:  minmax    80 33.72606 29.993020 25.637490 31.08364 28.710451 37.933588
                  B6        B8        B16       B18       B26
               <num>     <num>      <num>     <num>     <num>
        1: 15.222137  5.365861  0.8970604  2.939259 15.516402
        2: 10.753399 10.561023  1.4575628  8.013282  7.306200
        3:  9.603629  5.616830  8.1576152 16.020897 14.130052
        4:  8.287086 12.999384  7.1457154  6.542394  6.952915
        5: 11.821007 10.215768 12.8578863  8.408418 13.708609
       ---                                                   
      396: 42.861184 41.866046 21.9149196 15.964518 41.049195
      397: 39.543097 34.720623 47.8854240 48.952999 35.531254
      398: 18.209071 16.666973 20.3557738 18.072725 19.773843
      399: 20.970510 21.315809 15.1661056 23.930244 27.508591
      400: 30.095758 27.686008 24.1916472 24.933046 37.598064

# Calculating channel stats

    Code
      chan_stats
    Output
          electrode      means      sds variance  kurtosis   minmax
      A5         A5 -2.0201400 5.872054 34.48102 0.8110859 55.99172
      A13       A13 -0.7884497 3.922905 15.38918 0.8728392 39.18374
      A21       A21  3.2864957 6.088564 37.07061 1.4236771 68.28431
      A29       A29  4.1693878 7.781989 60.55936 0.1267196 57.26694
      A31       A31  1.0943630 5.607857 31.44806 0.8183168 47.41699
      B5         B5 -2.1856578 7.477209 55.90865 1.0617711 73.42224
      B6         B6 -2.8632122 5.963742 35.56622 0.3695797 46.77692
      B8         B8 -1.6502435 6.782543 46.00289 0.9660805 66.80369
      B16       B16 -1.4560545 5.200351 27.04365 2.0987794 51.42989
      B18       B18 -0.4696431 5.460927 29.82173 2.3864875 63.04594
      B26       B26  2.8831543 6.146965 37.78518 0.7989692 56.77661

# ar_thresh runs correctly.

    Code
      ar_thresh(demo_epochs, 20)
    Message
      250 (0.34%) samples above 20 uV threshold.
      205 (0.28%) samples below -20 uV threshold.
      57 epochs contain samples above threshold.
    Output
      Epoched EEG data
      
      Number of channels	: 11 
      Number of epochs	: 80 
      Epoch limits		: -0.197 - 0.451 seconds
      Electrode names		: A5 A13 A21 A29 A31 B5 B6 B8 B16 B18 B26 
      Sampling rate		: 128  Hz
      Reference		: average 

---

    Code
      ar_thresh(test_data, 30)
    Message
      10867 (4.42%) samples above 30 uV threshold.
      230400 (93.75%) samples below -30 uV threshold.
    Output
      EEG data
      
      Number of channels	: 16 
      Electrode names		: A1 A2 A3 A4 A5 A6 A7 A8 A9 A10 A11 A12 A13 A14 A15 A16 
      Sampling rate		: 256 Hz
      Reference		: 
      Signal length: 0 59.996 seconds

