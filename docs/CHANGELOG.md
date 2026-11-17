# Changelog

## [0.2.2] - 2026-11-17
### Added
- Added `--threads` as variable
- Added `genome/` for user installed input instead of being in `inputs/`
  - Added `genome/meryl.sh` for meryl DB setup
  - Added `genome/download_ref_chm13.sh` for reference genome installation and indexing

### Changed
- Moved user installed inputs (e.g. reference) to a separate directory
- Meryl setup no longer in `d4z4ling.sh` but in `genome/meryl.sh` for multiple runs purpose
- Adjusted default variable based on genome file moves in `d4z4ling.sh`
- Set tools used in script as variable
- Adjusted SSLP column to exist by default

### Fixed
- Fixed relative path issues
- Remove extra spacing in the file
- Remove NCI gadi user specific module calls in main scripts
- Updated setup documentation according to the changes in this version

## [0.2.1] - 2026-10-18
### Fixed
- Adjusted `helper/minimap_features.sh` to fix duplicate reads on rerun into the same output directory
- Change `--PLAM` to `--plam` in both `setup/RUN_ON_NCI_GADI.md` and `setup/RUN_ON_OTHER_HPC.md`

## [0.2.0] - 2026-07-09
### Added
- Added option to use conda environment
- Added `environment.yml` for conda env setup
- Conda documentation for different environment:
  - NCI gadi HPC: `README.gadi.md` to `setup/RUN_ON_NCI_GADI.md`
  - Other HPC: `setup/RUN_ON_OTHER_HPC.md`
- Added a `demo/check_summary_tsv.py` for sanity check of script installation
- Added `.gitignore`
- Added `software_version.md` 

### Changed
- File restructuring for organisation
- Changed `fshd1_script.sh` to `d4z4ling.sh` to emphasize main script
- Changed `fshd1_annotate_fasta_script.sh` to `d4z4ling_no_methylation.sh` to emphasize main script
- Updated `README.md` docs
- Moved `CHANGELOG.md` to `docs/`

### Fixed
- Minor changes in `helper/read_classification` to a type safe initial value for `percent_identity`

## [0.1.1] - 2026-07-03
### Added
- Empty file checks to exit without error for:
    - No D4Z4 region sequenced (`fshd1_script.sh`)
    - No/minimal methylation data (`fshd1_script.sh`)
    - No 4A161L haplotype (`distal_haplotype_blast.sh`)
    
## [0.1.0] - 2026-04-21
### Added
- Initial release