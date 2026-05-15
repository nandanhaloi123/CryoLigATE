![CryoLigate Pipeline](docs/Intro_for_webpage.png)

## Introduction
The full potential of cryo-EM in drug discovery remains limited by poor density resolvability
at ligand-binding interfaces. Although recent advances in deep learning have transformed cryo-
EM map enhancement, existing approaches largely focus on protein regions and often neglect
ligand-containing sites. Here, we present CryoLigate, an AI framework specifically designed
to enhance the density resolvability of protein–ligand interfaces. We trained and evaluated CryoLigate across a structurally diverse dataset including pharmaceutical drugs, lipids, steroids, and carbohydrates.

CryoLigate features a streamlined, single-command interface. It automatically isolates the target sub-volume using a preliminary atomic model as a spatial reference, requiring no manual box curation. The pipeline is computationally efficient, with localized refinement completed in seconds on standard desktop-grade GPU hardware.

## Installation

> **Note:** We strongly recommend installing CryoLigate in a fresh Python or Conda environment to avoid dependency conflicts.

### Set up the Conda Environment First
Create and activate the environment using the provided `environment.yml` file:
```bash
conda env create -f environment.yml
conda activate CryoLigate
```

Install CryoLigate via PyPI (Recommended):
```bash
pip install cryoligate -U
```

Or install directly from GitHub for the latest development updates:

```Bash
git clone [https://github.com/nandanhaloi123/CryoLigate.git](https://github.com/nandanhaloi123/CryoLigate.git)
cd CryoLigate
pip install -e .
```

If you are installing on CPU-only or non-CUDA GPU hardware, the pipeline will automatically fall back to CPU processing. Note that the CPU version is significantly slower than the GPU version for 3D volumetric refinement.


## Inference

Before running inference or fine-tuning, download the pre-trained weights:

```bash
mkdir weights
wget -O weights/cryoligate_v1.0.0.pth [https://github.com/nandanhaloi123/CryoLigate/releases/download/v1.0.0/cryoligate_v1.0.0.pth](https://github.com/nandanhaloi123/CryoLigate/releases/download/v1.0.0/cryoligate_v1.0.0.pth)
```


You can run inference using CryoLigate with:
```bash
CryoLigate-infer --weights cryoligate_v1.0.0.pth --map example/8ioe/emd_35617.map --pdb example/8ioe/8ioe.cif --resname TPP --chain A --resid 801
```
