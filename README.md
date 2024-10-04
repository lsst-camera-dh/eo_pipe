This package contains the LSST Camera team code for electro-optical test code used for sensor testing.
The algorithms are based on the original LSST CCD testing specification, [LCA-128](https://docushare.lsst.org/docushare/dsweb/Get/LCA-128),
and the CCD Electro-Optical Testing Methods document [LCA-10103](https://docushare.lsst.org/docushare/dsweb/Get/LCA-10103).

The code in this package analyzes outputs from the LSST calibration pipelines, [cp_pipe](https://github.com/lsst/cp_pipe), and implements some custom pipelines using standard LSST code.
In particular, the instrument signature removal code is from the [ip_isr](https://github.com/lsst/ip_isr) package.  The cp_pipe and ip_isr repos are useful references for the default config settings.
