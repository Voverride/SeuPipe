# SeuPipe: An Interactive Cloud-Based Collaborative Platform for Three-Dimensional Reconstruction of Spatial Transcriptomics

## Overview

**SeuPipe** is an interactive cloud-based collaborative analysis platform designed for three-dimensional reconstruction of spatial transcriptomics data. It provides an end-to-end workflow from raw tissue section images and gene expression matrices to single-cell-resolution expression matrices, cell annotation, spatial alignment, and 3D reconstruction outputs.

SeuPipe adopts a browser/server (B/S) architecture and supports private deployment in local laboratory environments. By keeping data storage and computation on local servers, the platform enables in situ analysis of large-scale spatial transcriptomics datasets while reducing cross-network data transfer overhead and improving data security.

The platform integrates interactive visualization, task progress monitoring, multi-user state synchronization, project-based data management, and fine-grained permission control. It is designed to help researchers perform spatial transcriptomics analysis in a more accessible, traceable, collaborative, and reproducible manner.

## Key Features

- **End-to-End Spatial Transcriptomics Workflow**  
  Supports the complete analysis chain from tissue image and gene expression matrix input to single-cell-resolution data generation and multi-slice spatial alignment.

- **Region Clipping and Sample Splitting**  
  Provides interactive region-of-interest selection for large tissue sections and multi-sample slices.

- **Cell Segmentation**  
  Integrates image-based cell segmentation methods to generate initial cell masks from tissue staining images.

- **Cell Boundary Expansion**  
  Supports cell boundary refinement and expansion to improve transcript assignment and single-cell expression matrix construction.

- **Segmentation Result Visualization and Export**  
  Enables visual comparison of segmentation results, mask overlay, contour display, and standardized data export.

- **Single-Cell Filtering**  
  Provides interactive scatterplot-based filtering, including manual selection and refined cell subset extraction.

- **Cell Annotation**  
  Supports annotation of spatially resolved cells by integrating reference single-cell transcriptomics data and spatial transcriptomics data.

- **Multi-Slice Spatial Alignment**  
  Supports automatic spatial alignment and manual fine-tuning for continuous tissue sections.

- **Interactive Web Interface**  
  Built with Dash, providing browser-based access and interactive visualization without requiring users to operate through command-line scripts.

- **Real-Time Collaboration and Task Monitoring**  
  Uses WebSocket-based communication to support task progress updates, multi-user state synchronization, and collaborative analysis.

- **Project-Based Data Management and Traceability**  
  Organizes analysis outputs by project and preserves intermediate results, parameters, and workflow states for reproducible analysis.

- **Permission Management**  
  Provides user authentication and role-based access control for collaborative laboratory use.

## Core Modules

SeuPipe consists of the following major functional modules:

1. **RegionClip**
2. **Segmentation**
3. **Expansion**
4. **MaskViewer**
5. **CellSelector**
6. **Annotation**
7. **Alignment**
8. **Passcode**

Together, these modules form a complete workflow for spatial transcriptomics data processing, single-cell-resolution expression matrix generation, annotation, and three-dimensional reconstruction.

## Installation

### 1. Clone the Repository

```bash
git clone https://github.com/Voverride/SeuPipe.git
cd SeuPipe
````

### 2. Create the Conda Environment

SeuPipe provides an `environment.yml` file for dependency installation. Create the environment by running:

```bash
conda env create -f environment.yml
```

### 3. Activate the Environment

```bash
conda activate SeuPipe
```

If your Conda version does not support `conda activate`, use:

```bash
source activate SeuPipe
```

## Running the Application

After installing the dependencies and activating the Conda environment, start the application from the project directory:

```bash
python main.py
```

By default, SeuPipe runs on port `8088` and automatically detects the available network IP address.

You can also specify a custom port:

```bash
python main.py --port <port_number>
```

For example:

```bash
python main.py --port 8088
```

After the application starts, open a web browser and visit:

```text
http://<auto-detected-network-IP>:8088
```

For example:

```text
http://10.86.x.x:8088
```

You should then see the SeuPipe web interface.

## Input Data

SeuPipe is designed to work with spatial transcriptomics data consisting of:

* **Tissue staining images**, such as TIFF images
* **Gene expression matrices**, such as CSV or TSV files
* **AnnData files**, such as `.h5ad` files, for standardized downstream analysis

The platform uses AnnData as a unified data container to store expression matrices, spatial coordinates, cell annotations, intermediate results, and downstream analysis outputs.

## Recommended Use Case

SeuPipe is suitable for research scenarios involving:

* Spatial transcriptomics preprocessing
* Single-cell-resolution expression matrix construction
* Tissue section cell segmentation
* Cell boundary refinement
* Spatial cell annotation
* Continuous tissue section alignment
* Three-dimensional reconstruction of spatial transcriptomics data
* Multi-user collaborative analysis in laboratory environments

## Citation

If you use SeuPipe in your research, please cite the corresponding thesis, paper, or project repository once available.

## Contact

For questions, issues, or contributions, please contact:

* **Email:** [1917926203@qq.com](mailto:1917926203@qq.com)
* **Repository:** [https://github.com/Voverride/SeuPipe](https://github.com/Voverride/SeuPipe)