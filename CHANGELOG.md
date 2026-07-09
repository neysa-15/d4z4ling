# Changelog

## [0.2.0] - 2026-07-09
### Added
- File restructuring for organisation
- Changed `fshd1_script.sh` to `d4z4ling.sh` to emphasize main script
- Changed `fshd1_annotate_fasta_script.sh` to `d4z4ling_no_methylation.sh` to emphasize main script
- Added option to use conda environment
- Added `environment.yml` for conda env setup
- Conda documentation for different environment:
  - NCI gadi HPC: `README.gadi.md` to `setup/RUN_ON_NCI_GADI.md`
  - Other HPC: `setup/RUN_ON_OTHER_HPC.md`
- Minor changes in `helper/read_classification` to a type safe initial value for "percent_identity"
- Added a `demo/check_summary_tsv.py` for sanity check of script installation
- Updated `README.md` docs
- Added `.gitignore`
- Added `software_version.md` 

## [0.1.1] - 2026-07-03
### Added
- Empty file checks to exit without error for:
    - No D4Z4 region sequenced (`fshd1_script.sh`)
    - No/minimal methylation data (`fshd1_script.sh`)
    - No 4A161L haplotype (`distal_haplotype_blast.sh`)
    
## [0.1.0] - 2026-04-21
### Added
- Initial release