# LazypipeX

Customizable metagenomics pipeline for fast and sensitive virus discovery from NGS data, with chained annotation strategies and interactive IGV reporting.

For prior versions of the pipeline (Lazypipe 1–3) see [bitbucket.org/plyusnin/lazypipe](https://bitbucket.org/plyusnin/lazypipe/).

## Documentation

Full documentation is available in [docs/UserGuide.v3.1.md](docs/UserGuide.v3.1.md).

## Features

- Quality preprocessing and background filtering (host/contaminant removal)
- De novo assembly and modular annotation strategies (one-round, two-round, chained)
- Search engines: Minimap2, BLASTN/P/X, DIAMOND, SANSparallel, HMMscan
- Databases: NCBI NT, RefSeq, UniRef100, Pfam, NeoRdRp
- Taxonomic binning with ICTV VMR metadata integration
- Interactive Krona graphs and IGV reports

## Citing LazypipeX

Weinstein I, Vapalahti O, Kant R, Smura T. *LazypipeX: Customizable Virome Analysis Pipeline Enabling Fast and Sensitive Virus Discovery from NGS data.* npj Viruses (2026, in press).

Earlier versions:

- Plyusnin I et al. Enhanced Viral Metagenomics with Lazypipe 2. *Viruses* 15(2):431 (2023). https://doi.org/10.3390/v15020431
- Plyusnin I et al. Novel NGS Pipeline for Virus Discovery from a Wide Spectrum of Hosts and Sample Types. *Virus Evolution*, veaa091 (2020). https://doi.org/10.1093/ve/veaa091

## Quick access on CSC

Lazypipe is available as a [preinstalled module](https://docs.csc.fi/apps/lazypipe/) at the [Finnish Centre of Scientific Computing](https://research.csc.fi/).

## License

GNU GPLv3 — see LICENSE and COPYRIGHT.

## Contact

Project website: https://www.helsinki.fi/en/projects/lazypipe  
Email: grp-lazypipe@helsinki.fi

