	This data collection contains spatially resolved single-cell transcriptomics datasets acquired using MERFISH on the adult whole mouse brain collected by the Xiaowei Zhuang Laboratory at Harvard University and Howard Hughes Medical Institute.

	This dataset contains MERFISH images of 92 experiments, which include 217 coronal (one hemisphere) slices and 28 sagittal slices of the mouse brain (10 um thick slices) collected from four male or female C57BL/6J mice aged 57–63 days: 
	Mouse1_coronal: male mouse sliced along the anterior-posterior axis at 200 um intervals, containing 67 coronal slices that passed QC.
	Mouse2_coronal: female mouse sliced along the anterior-posterior axis at 100 um intervals, containing 150 coronal slices that passed QC.
	Mouse3_sagittal: male mouse sliced along the medial-lateral axis at 200 um intervals, containing 25 sagittal slices that passed QC.
	Mouse4_sagittal: male mouse sliced along the medial-lateral axis at 200 um intervals, containing 3 sagittal slices that passed QC.

Raw images:
	Each folder named by the animal id contains the raw images of multiple experiments collected from one animal. Each experiment has a unique name (experiment_id) and contains one or multiple tissue slices. The raw images of each experiment were archived by tar and compressed by zstandard to create one single file named {experiment_id}.tar.zst. After decompression and unarchiving, each experiment folder contains image files of multiple fields of view (FOVs) from all imaging rounds in the "data" subfolder. Each dax file corresponds to the image of one FOV of one single round. PolyT and DAPI images are also included for cell segmentation. See corresponding data_organization.csv file for detailed channel and imaging round information for each experiment. In addition to the dax files, the inf files contain the data type, frame dimensions and number of frames information corresponding to the FOV, xml files contain filming and stage positions, and off files contain the z-plane positions corresponding to the FOV. 

Processed data:
	The processed_data folder contains cell segmentation files, decoded spots files, and cell by gene matrices.

	Cell boundaries files:
	The segmented cell boundary coordinates for each experiment are included in the {experiment _id}.csv files and are grouped by animal ids into separate folders. The index contains unique cell identifiers, and the columns contain the x, y coordinates of each z-plane for each segmented cell. All coordinates are in the unit of microns.
 
	Decoded spots files:
	The decoded spots for each experiment are included in the spots_{experiment_id}.csv files and are grouped by animal ids into separate folders. The files contain decoded spot location (global_x, global_y, global_z) in the unit of microns and their target gene identity for each experiment. The same coordinate system was used for the cell segmentation file and decoded spots file for each experiment, and hence the spots can be parsed into the segmented cells by comparison of their coordinates with the cell boundary locations. 

	Cell by gene matrices: 
	The counts_{animal_id}.h5ad files contain the cell by gene matrix (adata.X) for each animal as well as cell metadata information (adata.obs). The matrices contain RNA counts that are normalized by the cell volume of each cell and batch corrected based on the median count of each experiment.
	The raw_counts_{animal_id}.h5ad files contain the cell by gene matrix (adata.X) with raw RNA counts for each animal as well as cell metadata information (adata.obs). The matrices contain raw RNA counts that are measured, without normalizing by the cell volume of each cell or any further processing.


Additional files:
	We also provide these files and folders that are associated to this dataset:
•	experiment_metadata.csv: Provides the associating codebook and data organization files for each experiment.
•	codebook.csv: Provides the barcodes that encode individual genes measured in the combinatorial imaging rounds. Two different versions of codebooks containing a total of 1124 or 1147 genes were used in this dataset. Both encode individual genes with error-robust barcodes with 32-bit, Hamming distance 4, Hamming weight 4 binary code.
•	probes.fasta: Provides sequences of all encoding probes used for hybridization for the genes measured in the combinatorial imaging rounds. Two different versions of probes containing a total of 1124 or 1147 genes were used, corresponding to the two versions of codebooks described above.
•	probes_sequential.csv: Provides sequences of all encoding probes used for hybridization for the genes measured in the sequential imaging rounds.
•	probes_readout.csv: Provides sequences, fluorophore information, and purpose (indicates whether each probe was used in the combinatorial imaging rounds as a specific bit associated with the MERFISH barcodes or was used in the sequential rounds to measure individual genes) of all readout probes used in this dataset.
•	dataorganization folder: Contains the data organization files that associated with each experiment, each named as dataorganization_{experiment _id}.csv and were grouped by animal ids into separate subfolders. The data organization file provides information on how individual channels and z-planes are ordered for each experiment in the multi-frame dax file for each field of view of raw images.

CHANGE LOG:
March 15 2023 - The following files were added to the dataset in the "additional_files" directory:
probes_readout.xlsx  probes_sequential.xlsx  probes_v1.fasta  probes_v2.fasta
April 28 2023 - The file /bil/data/23/bil/data/29/3c/293cc39ceea87f6d/processed_data/raw_counts_mouse1_coronal.h5ad was updated to a new version. The README was updated 
accurately reflect the number of slices that passed QC for mouse 3 and mouse 4. Old versions may be available upon request.
