# [Target Acquisition Report](@id report_ta)

## Requirements
- target list (multiple formats supported), e.g., `TargetWizard.all.TW.target.csv`.
- traditional and targeted mass spectrometry data, e.g., `DDA.raw` and `TMS.raw`.
  - The raw data should be converted into an open-source format such as MS1/MS2. [ThermoRawRead](http://thermorawread.ctarn.io) is recommended.
- (filtered) identification results of targeted mass spectrometry data, e.g., `TMS_fdr.pfind.csv`.
- optional: precursor list detected by [`PepPre`](http://peppre.ctarn.io)

## Parameters
#### Max. MS1 Mass Error
mass error used to match targets, PSMs, and MS scans.

#### FDR Threshold
used to filter PSM list.

#### MS Sim. Thres.
used to match traditional and targeted MS scans.

## Output Results
Once finished, TargetWizard will save two reports to `Output Directroy`.
- `csv` report of all targets, e.g., `TMS.target.TargetAcquisitionReport.csv`.
- `csv` report of all PSMs, e.g., `TMS.psm.TargetAcquisitionReport.csv`.

## Usage
![Target Acquisition Report](../../assets/report/TargetAcquisitionReport.png)
