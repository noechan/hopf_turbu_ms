import os
import nibabel as nib
import numpy as np
import pandas as pd
from nilearn import datasets, image
import abagen

# ---------------------------------------------------------
# 1. Paths and output directory
# ---------------------------------------------------------

# Output folder: ./gene_niftis/schaefer1000_2mm/ (relative to current dir)
output_dir = os.path.join(os.getcwd(), "gene_niftis", "schaefer1000_2mm")
os.makedirs(output_dir, exist_ok=True)
print("Saving gene NIfTI maps to:", output_dir)

# ---------------------------------------------------------
# 2. Load Schaefer-1000 (7-networks, 2mm) atlas
# ---------------------------------------------------------

schaefer = datasets.fetch_atlas_schaefer_2018(
    n_rois=1000,
    yeo_networks=7,
    resolution_mm=2
)
atlas_img = nib.load(schaefer.maps)   # NIfTI image for the parcellation
atlas_data = atlas_img.get_fdata()    # 3D array with ROI indices
n_rois = int(atlas_data.max())
print(f"Schaefer atlas loaded with {n_rois} ROIs")

# ---------------------------------------------------------
# 3. Get AHBA expression data parcellated to Schaefer-1000
# ---------------------------------------------------------

expression_data = abagen.get_expression_data(
    atlas=atlas_img,      # pass atlas as NIfTI
    ibf_threshold=0,      # standard thresholding
    missing="centroids",  # handle missing data via centroid assignment
    norm_matched=False    # recommended when using missing="centroids"
)

print("Expression data shape:", expression_data.shape)
print("First few index labels:", expression_data.index[:5])
print("First few gene names:", expression_data.columns[:5])

# ---------------------------------------------------------
# 4. Define genes of interest and subset expression data
# ---------------------------------------------------------

genes_of_interest = [
    "MAPT", "APOE", "PICALM", "BIN1", "CLU", "CR1", "ABCA7", "SORL1", "PLEKHA1",
    "CD2AP", "CD33", "APP", "PSEN1", "PSEN2", "CASS4", "EPHA1", "PTK2B",
    "INPP5D", "MEF2C", "CELF1", "MADD"
]

# Keep only genes that are present in the expression matrix
available_genes = [g for g in genes_of_interest if g in expression_data.columns]
missing_genes = [g for g in genes_of_interest if g not in expression_data.columns]

print("Available genes:", available_genes)
if missing_genes:
    print("Warning: these genes were not found in the expression data and will be skipped:", missing_genes)

gene_expression_data = expression_data[available_genes]

# ---------------------------------------------------------
# 5. Function to create and save a gene-expression NIfTI
# ---------------------------------------------------------

def create_gene_nifti(gene_name, atlas_data, atlas_img, gene_expression_data, output_dir):
    """
    Generates and saves a NIfTI image for a given gene expression profile.

    Parameters
    ----------
    gene_name : str
        Name of the gene (must be a column in gene_expression_data).
    atlas_data : ndarray
        3D array with integer ROI labels (Schaefer-1000).
    atlas_img : Nifti1Image
        NIfTI image providing affine and header.
    gene_expression_data : DataFrame
        DataFrame with ROIs as index and genes as columns.
    output_dir : str
        Directory where the NIfTI will be saved.
    """
    nifti_data = np.zeros_like(atlas_data, dtype=float)

    # Determine how ROIs are indexed in the DataFrame (int vs str)
    index_as_str = isinstance(gene_expression_data.index[0], str)

    # Fill volume with gene expression values by ROI
    for roi_idx in range(1, int(atlas_data.max()) + 1):
        mask = (atlas_data == roi_idx)

        if index_as_str:
            roi_key = str(roi_idx)
        else:
            roi_key = roi_idx

        if roi_key in gene_expression_data.index:
            value = gene_expression_data.loc[roi_key, gene_name]
            nifti_data[mask] = value

    # Create NIfTI image and save
    nifti_img = nib.Nifti1Image(nifti_data, affine=atlas_img.affine)
    out_name = f"{gene_name}_expression_sch1000_2mm.nii.gz"
    out_path = os.path.join(output_dir, out_name)
    nib.save(nifti_img, out_path)

    print(f"Saved: {out_path}")
    return nifti_img

# ---------------------------------------------------------
# 6. Generate NIfTI images for all available genes
# ---------------------------------------------------------

for gene in available_genes:
    print(f"Processing gene: {gene}")
    create_gene_nifti(
        gene_name=gene,
        atlas_data=atlas_data,
        atlas_img=atlas_img,
        gene_expression_data=gene_expression_data,
        output_dir=output_dir
    )

print("All available gene NIfTI images have been generated.")
