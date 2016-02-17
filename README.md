# LMD Paper

This repository contains the code and text for the analysis and manuscript from the following article:

T.W.R. Harrop, I. Ud Din, V. Gregis, M. Osnato, S. Jouannic, H. Adam and M. Kater (2016). Gene expression profiling of reproductive meristem types in early rice inflorescences by laser microdissection. Accepted article, doi: **doi**

## Data availability

Sequence data from this article have been deposited with the National Center for Biotechnology Information Sequence Read Archive (SRA) under accession number [SRP067488](http://www.ncbi.nlm.nih.gov/sra/SRP067488).

## Repeatability

The analysis can be repeated by cloning the repository and downloading the data files from the SRA to the `data` directory. The steps in the analysis are shown in the pipeline PDF file in the `ruffus` directory. The steps may be run individually by executing the scripts in `src` in the order shown in the pipeline PDF. Provided the required software is installed, the whole pipeline may also be run by executing `src/py/pipeline.py`. `src/py/pipeline.py` will have to be modified to run on your system; it is currently set up to run jobs on the `SLURM` workload manager.

After running the analysis, the manuscript can be reproduced from the `ms.Rmd` file. Several files will be needed, including the style for for the Plant Journal (`csl/the-plant-journal.csl`), a template document for the MS Word document (`csl/wordStyle.docx`), a bibtex bibliography (`bib/papersLib.bibtex`) and two figures (`data/lmdFigure` and `data/isFigure.pdf`).

See the scripts, `SessionInfo.txt` and `METADATA.csv` files for the required software and the versions used in the analysis.

## Contact

For questions, contact the corresponding authors or open an issue in the GitHub repository.