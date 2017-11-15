# FACETS suite

# 1. Initial Parameters

Specify default setting for generating SNP counts and running FACETS

1.1. Min.het
------------

Use min.het=15 for purity and high-sensitivity runs

1.2. SNP set
------------

Stick with the one that we have\...

# 2. Pre-fit

Criteria for removing samples from consideration based on the SNP counts
only.

2.1. Low number of het SNPs
---------------------------

If there are fewer than 1500 het SNPs in an IMPACT sample, check sample.
If there is neither compelling CNA nor compelling diploidy, then reject
the sample.

2.2. Potential mis-match
------------------------

Use concordance per segment?

2.3. Contamination & Coverage
-----------------------------

These can be found in the pairing file.

# 3. Refit

Criteria for rejecting a FACETS fit and procedure for selecting
parameters for a subsequent run.

3.1. Segmentation
-----------------

### 3.1.1. hypersegmentation

If there are an implausibly large number of segments, try refitting with
a lower cval.

total number of segments, minimum number of segments per chr

### 3.1.2. hyposegmentation

If a single segment appears to be an \"average\" of multiple compelling
CNAs, then try refitting with a lower cval.

Use concordance per segment?

3.2. Diploid LogR
-----------------

### 3.2.1. Too low

If there are less than 2 balanced diploid segments, check if the diploid
logR has been set too high. If the diploid segments are likely artifacts
or if there are no compelling triploid segments, choose the next highest
logR corresponding to a balanced segment.

### 3.2.2. Too high

If more than 2% of the genome is covered by \"clonal\" (cf \> 0.5)
homdels or if more than 50% of the genome is covered by LOH, check if
the diploid logR has been set too high.

# 4. Post-fit

Final criteria for accepting a FACETS fit as well as procedures for
improving the purity estimate, gaps in the segments etc.
