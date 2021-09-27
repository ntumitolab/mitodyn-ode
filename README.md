# Mitochondrial dynamics model

## Model description

Mathematical descriptions are in [model_desc.md](model_desc.md).

You can export it to MS Word `.docx` file using [typora](https://typora.io) GUI or [pandoc](https://pandoc.org) in the CLI:


```bash
pandoc -s model_desc.md -o  model_desc.docx
```

## How to run the code

### Installing Julia

Download and install Julia (*version 1.6+*) from [the official website](https://julialang.org/downloads/). Choose "add Julia to system PATH" during installation.

### Download the code

Open the terminal (either Windows terminal, powershell, or bash) and enter the following command:

```bash
git clone https://github.com/NTUMitoLab/MitochondrialDynamics.git
cd MitochondrialDynamics
```

If that does not work, download this repositroy via the big green `Code` button => `download zip`. You will see the project folder `MitochondrialDynamics` in the extracted zip archive.

### Running the code

In the terminal, type the following commands to run the code to generate figures and data. Make sure the working directory (output of `pwd()`) is the project folder `MitochondrialDynamics`.

```bash
pwd # Make sure it's the project directory, MitochondrialDynamics
julia --color=yes allfigs.jl
```

### Trouble-shooting


## Appendix: Convert pdf figures to tiff ones

*Might not work in Windows systems*

`poppler-utils` is needed for `pdftoppm`. One can install it in the software repos or via `conda` package manager.

```bash
conda install poppler -c conda-forge
```

The following command converts a `pdf` file to a 300 DPI `tiff` file using `lzw` compression.

```bash
pdftoppm -tiff -tiffcompression lzw -r 300 Fig2.pdf
```
