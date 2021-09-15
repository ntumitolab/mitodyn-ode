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

### Environment Setup

1. Download this repositroy via the green button.
2. Extract the downloaded archive, the `zip` file. You will see the project folder `MitochondrialDynamics`.
3. Choose `Open with terminal` in the folder's right-click context menu.
4. In the terminal, type the following to install dependencies

```bash
julia --project=@. 'import Pkg; Pkg.instantiate()'
```

### Running the code

In the terminal, type the following to generate the figures. For example, to generate Figure 2:

```bash
julia --project=@. fig2.jl
```

### Trouble-shooting


## How to convert pdf figures to tiff

`poppler-utils` is needed for `pdftoppm`. In Windows, you can install it via the [conda package manager](https://www.anaconda.com/products/individual):

```bash
conda install poppler -c conda-forge
```

To convert a `pdf` file to a 300 DPI `tiff` file using `lzw` compression.

```bash
pdftoppm -tiff -tiffcompression lzw -r 300 in.pdf outname
```
