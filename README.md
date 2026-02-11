#raman_spec_tools

---

Version 1.1.

2/10/2026





# Motivation

---

 The purpose of this repo is have only the rs_tools.py code. We (Doughty and Hill) also believe that the code and data use to do science should be made easily available after publication of a paper. This way others can easily use it for their own research, and also potentially detect bugs or problems in the software.

# Description 

---

There are two main steps, for spectral analysis: preprocessing and analysis. As preprocessing can take a long time, this step is done separately from the 'analysis' step, and data is saved in .rdta, .rmta, and .rflu files, corresponding to data, meta, and fluorescence files, respectively. Preprocessing steps can include removal of bad pixels, removal of cosmic rays, calculation of background, removal of background, removal of fluorescence, removal of saturated spectra, and interpolation across bad pixels, and combination of replicates into an average, median, or summed spectrum. Preprocessing is done in python. The file [rs_tools](rs_tools.py) includes routines to accomplish all of these tasks, as well as many other routines for plotting, and general processing of Raman spectral data. At the end of the preprocessing step, data are saved as two dimensional files, with each row (other than the header) containing a Raman spectrum (.rdta), a fluorescence spectrum (.rflu), or metadata associated with the associated spectra (.rmta).


# Public Domain Dedication

The text of the CC0 license can be found [here](https://creativecommons.org/publicdomain/zero/1.0/legalcode.txt). A local copy of the license is in the file license.txt.
