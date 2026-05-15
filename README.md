![CryoLigate Pipeline](docs/assets/Intro.png)

# Introduction
The full potential of cryo-EM in drug discovery remains limited by poor density resolvability
at ligand-binding interfaces. Although recent advances in deep learning have transformed cryo-
EM map enhancement, existing approaches largely focus on protein regions and often neglect
ligand-containing sites. Here, we present CryoLigate, an AI framework specifically designed
to enhance the density resolvability of protein–ligand interfaces. We trained and evaluated CryoLigate across a structurally diverse dataset including pharmaceutical drugs, lipids, steroids, and carbohydrates.

CryoLigate features a streamlined, single-command interface. It automatically isolates the target sub-volume using a preliminary atomic model as a spatial reference, requiring no manual box curation. The pipeline is computationally efficient, with localized refinement completed in seconds on standard desktop-grade GPU hardware.

# CryoLigate
An AI framework specifically designed to enhance the density resolvability of protein–ligand interface.


## 📥 Download Weights

Before running inference or fine-tuning, download the pre-trained weights:

```bash
mkdir weights
wget -O weights/best_model.pth [https://github.com/nandanhaloi123/CryoLigate/releases/download/v1.0.0/best_model.pth](https://github.com/nandanhaloi123/CryoLigate/releases/download/v1.0.0/best_model.pth)
