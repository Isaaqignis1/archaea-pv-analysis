# archaeal-phasome

Scripts for identifying putative phase-variable (PPV) loci in archaeal
genomes, with three bacterial control genera as reference points.

The repo contains SLURM jobs that run Prokka, PhasomeIt and eggNOG-mapper,
the python scripts that parse PhasomeIt HTML output, the scripts that
count and normalise PPV loci, the script that joins PhasomeIt PV genes
to eggNOG annotations, and three Jupyter notebooks that produce the
report figures.

## Layout

```
repo/
├── README.md
├── LICENSE
├── config.sh.example              copy to config.sh and edit
├── .gitignore
│
├── standards/
│   ├── colors.py                  per-genus colour map and plot ordering
│   ├── SpeciesCallsArchaea.csv
│   └── SpeciesCallsBacteria.tsv
│
├── annotation/
│   ├── prokka_archaea.job         prokka, kingdom Archaea
│   └── prokka_controls.job        prokka, kingdom Bacteria
│
├── phasomeit/
│   ├── setup_archaea.py           stage per-genus folders from <Genus>_gffs.txt
│   ├── setup_controls.py          stage per-genus folders from a sample/species tsv
│   └── phasomeit_run.job          phasomeit array runner, cutoffs -c 7 6 0 5 5
│
├── parsing/
│   ├── parse_groups.py            phasomeit group html -> 4 csvs
│   ├── parse_tracts.py            phasomeit strain html -> tract rows
│   └── summarise_runs.py          per-genus run completeness
│
├── eggnog/
│   ├── eggnog_input_merger.job    concat prokka .faa files with sample-id prefixes
│   ├── run_eggnog_archaea.job     emapper, tax_scope Archaea
│   ├── run_eggnog_bacteria.job    emapper, default tax scope
│   ├── parse_eggnog.py            emapper.annotations -> tidy tsvs with broad roles
│   └── link_eggnog.py             join phasomeit pv members to eggnog by locus tag
│
├── quantification/
│   ├── pv_per_fasta.py            pv count per genome
│   ├── genome_lengths.py          genome size in bp and mb
│   └── pv_per_mb.py               per-genome and per-genus pv/mb
│
└── figures/
    ├── build_master_tables.ipynb
    ├── plot_pv_burden.ipynb
    └── plot_cog_functions.ipynb
```



