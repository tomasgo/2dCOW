2DCOW
=====

1. Please refer details to

D. Zhang, X. Huang, F.E. Regnier, and M. Zhang (2008). Two-dimensional correlation optimized warping algorithm for aligning GCXGC-MS data. Analytical Chemistry, 80 (8): 2664-2671.

2. It calls the function "cow" in the WarpingTB package, which is available from <http://www.models.kvl.dk/source/DTW_COW/ directly.

3. Use "Batch2DAlign" (in Batch2DAlign.m) or "AlignGC2Data" (in AlignGC2Data.m) to align a batch of 2-dimensional images.

4. I would suggest to correct shifts between images, using "OptimShift" (in OptimShift.m) and "ShiftImage" (in ShiftImage.m), before alignment.

5. 2DCOW is also available at <http://ccehub.org/tools/gc2msclass> for free on-site use (i.e., you can register and load your data to use the pacakge through a friendly GUI).