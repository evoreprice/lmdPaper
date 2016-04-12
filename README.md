## LMD Paper

* This repository contains the code to generate the manuscript, figures and data for the following article:

    > Harrop, T. W.R., Ud Din, I., Gregis, V., Osnato, M., Jouannic, S., Adam, H. and Kater, M. (2016), Gene expression profiling of reproductive meristem types in early rice inflorescences by laser microdissection. The Plant Journal, 86: 75–88. [doi:10.1111/tpj.13147](http://dx.doi.org/10.1111/tpj.13147)

### Data availability

* Sequence data from this article have been deposited with the National Center for Biotechnology Information Sequence Read Archive (SRA) under accession number [SRP067488](http://www.ncbi.nlm.nih.gov/sra/SRP067488).

### Reproducibility

* The analysis can be repeated by cloning the repository and downloading the data files from the SRA to the `data` directory. The steps in the analysis are shown in the pipeline PDF file in the `ruffus` directory. The steps may be run individually by executing the scripts in `src` in the order shown in the pipeline PDF. Provided the required software is installed, the whole pipeline may also be run by executing `src/py/pipeline.py`. `src/py/pipeline.py` is set up to run jobs using the `SLURM` workload manager and will require modification to run on other systems.

* See the scripts, `SessionInfo.txt` and `METADATA.csv` files for the required software and the versions used in the analysis.

* After running the analysis, the manuscript can be reproduced from the `ms.Rmd` file. Several external files will be needed, including the style for for the Plant Journal (`csl/the-plant-journal.csl`), a template document for the MS Word document (`csl/wordStyle.docx`), a bibtex bibliography (`bib/papersLib.bibtex`) and two figures (`data/lmdFigure` and `data/isFigure.pdf`).

### Contact

* For questions, contact the corresponding authors or open an issue in the GitHub repository.