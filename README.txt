# 2D Morphometric Analysis of Hominin Crania

This repository contains a suite of MATLAB scripts for the two-dimensional geometric morphometric analysis of hominin crania. The software uses Fourier analysis to parameterize skull outlines from 2D images and calculates various metrics related to skull shape.

This code is associated with the following publication:
> Mishol N., Herzlinger G., Rak Y., Smilansky U., Carmel L., Gokhman D. (2025). Candidate Denisovan fossils identified through gene regulatory phenotyping. PNAS

---

## 🧭 Overview

The core methodology involves extracting the 2D outline of a skull from an image, representing this outline using Fourier coefficients, and then performing a series of automated and semi-automated analyses to quantify shape. The toolbox is divided into three main analytical modules based on the cranial view.

---

## ✨ Features

* **Lateral View Analysis (`SklMorph2D.m`):** Computes skull top flatness by fitting a parametric function to the vault and measures forehead height and projection based on curve tangents.
* **Inferior View Analysis (`SkulInfView.m`):** Analyzes the geometry of the zygomatic bones from an inferior view aligned to the Frankfurt plane.
* **Superior View Analysis (`SkulSupView.m`):** Quantifies the globularity of the frontal bone in the glabella region from a superior view.
* **Batch Processing:** Includes wrapper scripts to analyze multiple images at once for each view.
* **Data Export:** Automatically saves results to `.xlsx` files and generates high-quality SVG figures for publication.

---

## 🔧 Requirements

* **MATLAB** (R2021a or newer recommended)
* **Statistics and Machine Learning Toolbox** (for the `pdist2` function)
* **Computer Vision Toolbox** (for the `ransac` function)

---

## ⚙️ Installation

1.  Clone this repository or download and extract the ZIP file to a local directory.
    ```sh
    git clone <repository-url>
    ```
2.  Open MATLAB.
3.  Add the project directory to your MATLAB path using `addpath` or the "Set Path" dialog.
    ```matlab
    addpath('path/to/your/project/folder');
    ```

---

## 🚀 Usage

The primary way to use this software is through the batch processing scripts, which allow for the analysis of multiple specimens at once.

1.  Organize your skull images (e.g., `.jpg`, `.png`) into a single folder. The images should be of skulls on a plain white background.
2.  Run the appropriate batch script from the MATLAB command window for the desired analysis:
    * For lateral view (flatness/forehead): `Multi_Read()` 
    * For inferior view (zygomatic): `Multi_Read_InfView()`
    * For superior view (globularity): `Multi_Read_SupView()` 
3.  A file selection dialog will open. Navigate to your image folder and select all the images you wish to analyze.
4.  For analyses requiring landmarking (`SklMorph2D`, `SkulInfView`), an interactive window will appear for each image. Follow the on-screen prompts to select the required anatomical points.
5.  The scripts will process all selected files, and the results (`.xlsx` data sheets and `.svg` figures) will be automatically saved in the same folder as your images.

---

## 📜 Citation

If you use this software in your research, please cite the associated publication:

> Mishol N., Herzlinger G., Rak Y., Smilansky U., Carmel L., Gokhman D. (2025). Candidate Denisovan fossils identified through gene regulatory phenotyping. PNAS

---

## ⚖️ License

This project is licensed under the MIT License. See the [LICENSE]file for details.