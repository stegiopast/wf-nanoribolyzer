from pathlib import Path
import pod5
import remora
from remora import io, refine_signal_map, util
import numpy as np
import argparse
import pysam
from tqdm import tqdm
import math
from itertools import repeat
import json
import re
import pandas as pd
import time
import polars as pl
import h5py
import os


# Bascalling and alignment should have been executed with the following flags in order for this script to function
#
# dorado basecaller ~/dorado-0.4.3-linux-x64/bin/rna002_70bps_hac@v3 extraction_final.pod5 --emit-moves |\
#       ~/samtools-1.18/bin/samtools fastq -T "*" | minimap2 -y --MD -ax map-ont Gm_randomers.fasta - |\
#       ~/samtools-1.18/bin/samtools view -b -F 4 > extraction_final.bam
#
# Move table and MD tags must be integrated in the final bam file


class RawCurrentReadDataset:
    """
    A class to initialize and process signal data from POD5 and BAM files, supporting both basecall and reference-anchored analysis.

    Attributes:
    ----------
    id : str
        The read ID used to identify and load the data.
    pod5_read : optional
        Single read from a POD5 dataset (default is None).
    pod5_dr : pod5.dataset.DatasetReader
        DatasetReader object to access POD5 data (default is None).
    pod5_filename : Path
        Path to the POD5 file.
    bam_read : optional
        Single read alignment from a BAM file (default is None).
    bam_fh : remora.io.ReadIndexedBam
        BAM file handler used to read the alignment.
    bam_filename : Path
        Path to the BAM file.
    kmer_table : Path
        Path to the k-mer model file used for refining signal mapping.
    sig_map_refiner : optional
        Signal mapping refiner to map signals to k-mers (default is None).

    Methods:
    -------
    __init__()
        Initializes the dataset by loading POD5 and BAM data, processing both basecall-anchored and reference-anchored signals.
    """
    def __init__(
        self,
        id: str = "",
        pod5_read=None,
        pod5_dr: pod5.dataset.DatasetReader = None,
        pod5_filename: Path = "",
        bam_read=None,
        bam_fh: remora.io.ReadIndexedBam = None,
        bam_filename: Path = "",
        kmer_table: Path = "",
        sig_map_refiner=None,
    ):
        """
        Initializes a RawCurrentReadDataset object by loading POD5 and BAM data,
        as well as basecall-anchored and reference-anchored signal information.

        Parameters:
        ----------
        id : str
            The read ID for identifying the read in the POD5 and BAM files.
        pod5_read : optional
            Single read from POD5 dataset, if available.
        pod5_dr : pod5.dataset.DatasetReader
            DatasetReader object to access POD5 data.
        pod5_filename : Path
            Path to the POD5 file.
        bam_read : optional
            Single read alignment from BAM file, if available.
        bam_fh : remora.io.ReadIndexedBam
            BAM file handler to read alignment data.
        bam_filename : Path
            Path to the BAM file.
        kmer_table : Path
            Path to the k-mer model file used for refining signal mapping.
        sig_map_refiner : optional
            Signal mapping refiner to map signals to k-mers.

        Raises:
        ------
        AttributeError
            If an attribute is missing or fails to initialize.
        remora.RemoraError
            If Remora-specific errors occur.
        """
        try:
            self.id = id
            self.kmer_table = kmer_table
            self.sig_map_refiner = sig_map_refiner
            # Check if SigMapRefiner has been loaded before the class construction
            if self.sig_map_refiner is None:
                self.sig_map_refiner = refine_signal_map.SigMapRefiner(
                    kmer_model_filename=self.kmer_table,
                    do_rough_rescale=True,
                    scale_iters=1,
                    do_fix_guage=True,
                )
            # Check if single pod5 read has been loaded before the class construction
            if pod5_read != None:
                self.pod5_read = pod5_read
            # Check if pod5 DataSetReader has been loaded before the class construction
            elif pod5_dr is None:
                self.pod5_dr = pod5.DatasetReader(pod5_filename)
                self.pod5_read = self.pod5_dr.get_read(id)
            else:
                self.pod5_dr = pod5_dr
                self.pod5_read = self.pod5_dr.get_read(id)
            # Check if single bam has been loaded before class construction
            if bam_read != None:
                self.bam_read = bam_read
            # Check if bam dataset has been loaded before class construction
            elif bam_fh is None:
                self.bam_fh = io.ReadIndexedBam(bam_filename)
                self.bam_read = self.bam_fh.get_first_alignment(id)
            else:
                self.bam_fh = bam_fh
                self.bam_read = self.bam_fh.get_first_alignment(id)
            # Determine read direction
            if self.bam_read.is_forward:
                flip = True
            else:
                flip = False

            # Initialize read for basecall anchored analysis
            self.io_read_basecall = io.Read.from_pod5_and_alignment(
                self.pod5_read, self.bam_read, reverse_signal=flip
            )
            self.io_read_basecall.set_refine_signal_mapping(
                self.sig_map_refiner, ref_mapping=False
            )
            # Extract basecall anchored information
            self.basecall_anchored_data = self.io_read_basecall.extract_basecall_region(
                end_base=self.io_read_basecall.ref_seq_len
            )
            self.basecall_sequence = str(self.basecall_anchored_data.seq)
            self.basecall_sequence_to_signal_map = (
                self.basecall_anchored_data.seq_to_sig_map
            )
            self.basecall_norm_signal = self.basecall_anchored_data.norm_signal
            basecall_dwell_mean_sd = self.io_read_basecall.compute_per_base_metric(
                "dwell_mean_sd", ref_mapping=False
            )
            basecall_dwell_trimmean_sd = self.io_read_basecall.compute_per_base_metric(
                "dwell_trimmean_trimsd", ref_mapping=False
            )
            self.basecall_mean = basecall_dwell_mean_sd["mean"]
            self.basecall_sd = basecall_dwell_mean_sd["sd"]
            self.basecall_trimmean = basecall_dwell_trimmean_sd["trimmean"]
            self.basecall_trimsd = basecall_dwell_trimmean_sd["trimsd"]
            self.basecall_dwell = basecall_dwell_trimmean_sd["dwell"]
            self.basecall_resquiggled = []
            for index, value in enumerate(self.basecall_sequence):
                associated_signal = self.basecall_anchored_data.norm_signal[
                    self.basecall_anchored_data.seq_to_sig_map[
                        index
                    ] : self.basecall_anchored_data.seq_to_sig_map[index + 1]
                ]
                self.basecall_resquiggled.append(list(associated_signal))

            # Initialize read for reference anchored analysis
            self.io_read_reference = io.Read.from_pod5_and_alignment(
                self.pod5_read, self.bam_read, reverse_signal=True
            )
            self.io_read_reference.set_refine_signal_mapping(
                self.sig_map_refiner, ref_mapping=False
            )
            self.io_read_reference.set_refine_signal_mapping(
                self.sig_map_refiner, ref_mapping=True
            )
            # Extract reference anchored information
            self.reference_anchored_data = self.io_read_reference.extract_ref_reg(
                self.io_read_reference.ref_reg
            )
            self.reference_sequence = str(self.reference_anchored_data.seq)
            self.contig = self.reference_anchored_data.ref_reg.ctg
            self.reference_norm_signal = self.reference_anchored_data.norm_signal
            reference_dwell_mean_sd = self.io_read_reference.compute_per_base_metric(
                "dwell_mean_sd", ref_mapping=True
            )
            reference_dwell_trimmean_trimsd = (
                self.io_read_reference.compute_per_base_metric(
                    "dwell_trimmean_trimsd", ref_mapping=True
                )
            )
            self.reference_mean = reference_dwell_mean_sd["mean"]
            self.reference_sd = reference_dwell_mean_sd["sd"]
            self.reference_trimmean = reference_dwell_trimmean_trimsd["trimmean"]
            self.reference_trimsd = reference_dwell_trimmean_trimsd["trimsd"]
            self.reference_dwell = reference_dwell_mean_sd["dwell"]
            self.reference_resquiggled = []
            for index, value in enumerate(self.reference_sequence):
                associated_signal = self.reference_norm_signal[
                    self.reference_anchored_data.seq_to_sig_map[
                        index
                    ] : self.reference_anchored_data.seq_to_sig_map[index + 1]
                ]
                self.reference_resquiggled.append(list(associated_signal))
            # Index read to reference alignment
            index_reference = [
                i
                for i in range(
                    self.reference_anchored_data.ref_reg.start,
                    self.reference_anchored_data.ref_reg.end,
                )
            ]
            self.index_ref_dict = {}
            for index, value in enumerate(self.reference_sequence):
                self.index_ref_dict[index_reference[index]] = index

        except AttributeError:
            print(f"An Attribute error occured with read {id}")
        except remora.RemoraError:
            pass


    def extract_signal_reference_coordinates(
        self, coordinate: int = 0, bases_upstream: int = 4, bases_downstream: int = 4
    ):
        """
        Extracts signal data from a specified coordinate within a reference sequence, 
        including a defined number of bases upstream and downstream from the coordinate.

        Parameters:
        ----------
        coordinate : int
            The specific coordinate within the reference sequence from which to extract data (default is 0).
        bases_upstream : int
            Number of bases upstream from the coordinate to include in the extracted region (default is 4).
        bases_downstream : int
            Number of bases downstream from the coordinate to include in the extracted region (default is 4).

        Returns:
        -------
        tuple
            A tuple containing:
                - extracted_motif: Sequence segment surrounding the coordinate, including upstream and downstream bases.
                - extracted_signal: List of measurements from the resquiggled signal in the extracted region.
                - extracted_signal_means: Mean signal values in the extracted region.
                - extracted_signal_trimmeans: Trimmed mean signal values in the extracted region.
                - extracted_signal_dwell: Dwell time data in the extracted region.

            If an error occurs during extraction, returns a tuple of None values.

        Raises:
        ------
        KeyError, AttributeError, IndexError
            In case of missing data attributes or out-of-bounds coordinates.
        """
        try:
            extracted_motif = self.reference_sequence[
                self.index_ref_dict[coordinate]
                - bases_upstream : self.index_ref_dict[coordinate]
                + (bases_downstream + 1)
            ]
            extracted_signal = []
            for squiggle in self.reference_resquiggled[
                self.index_ref_dict[coordinate]
                - bases_upstream : self.index_ref_dict[coordinate]
                + (bases_downstream + 1)
            ]:
                for measurement in squiggle:
                    extracted_signal.append(measurement)
            extracted_signal_means = self.reference_mean[
                self.index_ref_dict[coordinate]
                - bases_upstream : self.index_ref_dict[coordinate]
                + (bases_downstream + 1)
            ]
            extracted_signal_trimmeans = self.reference_trimmean[
                self.index_ref_dict[coordinate]
                - bases_upstream : self.index_ref_dict[coordinate]
                + (bases_downstream + 1)
            ]
            extracted_signal_dwell = self.reference_dwell[
                self.index_ref_dict[coordinate]
                - bases_upstream : self.index_ref_dict[coordinate]
                + (bases_downstream + 1)
            ]
            return (
                extracted_motif,
                extracted_signal,
                extracted_signal_means,
                extracted_signal_trimmeans,
                extracted_signal_dwell,
            )
        except:
            return None, None, None, None, None

    def extract_signal_reference_motif(
        self,
        motif: str = "AG",
        all_or_first: str = "all",
        bases_upstream: int = 0,
        bases_downstream: int = 0,
    ):
        """
        Extracts signal data from a specified motif within a reference sequence.

        Parameters:
        ----------
        motif : str
            The motif sequence to search for within the reference sequence (default is "AG").
        all_or_first : str
            Determines whether to extract data from the "first" instance of the motif or from "all" instances (default is "all").
        bases_upstream : int
            Number of bases upstream from the motif to include in the extracted region (default is 0).
        bases_downstream : int
            Number of bases downstream from the motif to include in the extracted region (default is 0).

        Returns:
        -------
        tuple
            If `all_or_first` is "first":
                - extracted_motif: The sequence containing the motif plus upstream and downstream bases.
                - extracted_signal: List of measurements from the resquiggled signal within the extracted region.
                - extracted_signal_means: Mean signal values in the extracted region.
                - extracted_signal_trimmeans: Trimmed mean signal values in the extracted region.
                - extracted_signal_dwells: Dwell time data in the extracted region.

            If `all_or_first` is "all":
                - extracted_motifs: List of sequences each containing one motif with surrounding upstream and downstream bases.
                - extracted_signals: List of lists of resquiggled signal measurements for each extracted region.
                - extracted_signals_means: List of mean signal values for each extracted region.
                - extracted_signals_trimmeans: List of trimmed mean signal values for each extracted region.
                - extracted_signals_dwells: List of dwell times for each extracted region.

            If an error occurs or if lengths are inconsistent, returns a tuple of None values.

        Raises:
        ------
        KeyError
            If extraction from the reference motif fails due to a missing key in the data.
        AttributeError
            If extraction fails due to missing data attributes.
        """
        try:
            if all_or_first == "first":
                coordinate = self.reference_sequence.find(motif)
                length = len(motif)
                extracted_motif = self.reference_sequence[
                    coordinate
                    - bases_upstream : coordinate
                    + bases_downstream
                    + length
                    + 1
                ]
                extracted_signal = []
                extracted_signal_means = self.reference_mean[
                    coordinate
                    - bases_upstream : coordinate
                    + bases_downstream
                    + length
                    + 1
                ]
                extracted_signal_trimmeans = self.reference_trimmean[
                    coordinate
                    - bases_upstream : coordinate
                    + bases_downstream
                    + length
                    + 1
                ]
                extracted_signal_dwells = self.reference_dwell[
                    coordinate
                    - bases_upstream : coordinate
                    + bases_downstream
                    + length
                    + 1
                ]
                for squiggle in self.reference_resquiggled[
                    coordinate
                    - bases_upstream : coordinate
                    + bases_downstream
                    + length
                    + 1
                ]:
                    for measurement in squiggle:
                        extracted_signal.append(measurement)
                if len(extracted_signal_means) == (
                    len(motif) + bases_upstream + bases_downstream
                ):
                    return (
                        extracted_motif,
                        extracted_signal,
                        extracted_signal_means,
                        extracted_signal_trimmeans,
                        extracted_signal_dwells,
                    )
                else:
                    return None, None, None, None, None
            elif all_or_first == "all":
                coordinates = [
                    (m.start() - bases_upstream, m.end() + bases_downstream)
                    for m in re.finditer(motif, self.reference_sequence)
                ]
                extracted_motifs = []
                extracted_signals = []
                extracted_signals_means = []
                extracted_signals_trimmeans = []
                extracted_signals_dwells = []
                for coordinate in coordinates:
                    current_extracted_signal = []
                    current_extracted_signal_means = self.reference_mean[
                        coordinate[0] : coordinate[1]
                    ]
                    current_extracted_signal_trimmeans = self.reference_trimmean[
                        coordinate[0] : coordinate[1]
                    ]
                    current_extracted_signal_dwells = self.reference_dwell[
                        coordinate[0] : coordinate[1]
                    ]
                    current_extracted_motif = self.reference_sequence[
                        coordinate[0] : coordinate[1]
                    ]

                    for squiggle in self.reference_resquiggled[
                        coordinate[0] : coordinate[1]
                    ]:
                        [
                            current_extracted_signal.append(measurement)
                            for measurement in squiggle
                        ]

                    if len(current_extracted_signal_means) == (
                        len(motif) + bases_upstream + bases_downstream
                    ):
                        extracted_motifs.append(current_extracted_motif)
                        extracted_signals.append(current_extracted_signal)
                        extracted_signals_means.append(current_extracted_signal_means)
                        extracted_signals_trimmeans.append(
                            current_extracted_signal_trimmeans
                        )
                        extracted_signals_dwells.append(current_extracted_signal_dwells)
                return (
                    extracted_motifs,
                    extracted_signals,
                    extracted_signals_means,
                    extracted_signals_trimmeans,
                    extracted_signals_dwells,
                )
        except KeyError:
            # print("Extraction of signal from reference motif failed.\n")
            return None, None, None, None, None
        except AttributeError:
            # print("Extraction of signal from reference motif failed.\n")
            return None, None, None, None, None

    def extract_signal_basecall_motif(
        self,
        motif: str = "AG",
        all_or_first: str = "all",
        bases_upstream: int = 0,
        bases_downstream: int = 0,
    ):
        """
        Extracts signal data from a specified motif within a basecalled sequence.

        Parameters:
        ----------
        motif : str
            The motif sequence to search for within the basecall sequence (default is "AG").
        all_or_first : str
            Specifies whether to extract data from the "first" instance of the motif or from "all" instances (default is "all").
        bases_upstream : int
            Number of bases upstream from the motif to include in the extracted region (default is 0).
        bases_downstream : int
            Number of bases downstream from the motif to include in the extracted region (default is 0).

        Returns:
        -------
        tuple
            If `all_or_first` is "first":
                - extracted_motif: The sequence containing the motif with upstream and downstream bases.
                - extracted_signal: List of measurements from the resquiggled signal within the extracted region.
                - extracted_signal_means: Mean signal values in the extracted region.
                - extracted_signal_trimmeans: Trimmed mean signal values in the extracted region.
                - extracted_signal_dwells: Dwell time data in the extracted region.

            If `all_or_first` is "all":
                - extracted_motifs: List of sequences each containing one motif with surrounding upstream and downstream bases.
                - extracted_signals: List of lists of resquiggled signal measurements for each extracted region.
                - extracted_signals_means: List of mean signal values for each extracted region.
                - extracted_signals_trimmeans: List of trimmed mean signal values for each extracted region.
                - extracted_signals_dwells: List of dwell times for each extracted region.

            If an error occurs or if lengths are inconsistent, returns a tuple of None values.

        Raises:
        ------
        KeyError, AttributeError
            In case of missing data attributes or if extraction from the basecall motif fails.
        """
        try:
            if all_or_first == "first":
                coordinate = self.basecall_sequence.find(motif)
                length = len(motif)
                extracted_motif = self.basecall_sequence[
                    coordinate
                    - bases_upstream : coordinate
                    + bases_downstream
                    + length
                    + 1
                ]
                extracted_signal = []
                extracted_signal_means = self.basecall_mean[
                    coordinate
                    - bases_upstream : coordinate
                    + bases_downstream
                    + length
                    + 1
                ]
                extracted_signal_trimmeans = self.basecall_trimmean[
                    coordinate
                    - bases_upstream : coordinate
                    + bases_downstream
                    + length
                    + 1
                ]
                extracted_signal_dwells = self.basecall_dwell[
                    coordinate
                    - bases_upstream : coordinate
                    + bases_downstream
                    + length
                    + 1
                ]

                for squiggle in self.basecall_resquiggled[
                    coordinate
                    - bases_upstream : coordinate
                    + bases_downstream
                    + length
                    + 1
                ]:
                    for measurement in squiggle:
                        extracted_signal.append(measurement)
                if len(extracted_signal_means) == (
                    len(motif) + bases_upstream + bases_downstream
                ):
                    return (
                        extracted_motif,
                        extracted_signal,
                        extracted_signal_means,
                        extracted_signal_trimmeans,
                        extracted_signal_dwells,
                    )
                else:
                    return None, None, None, None, None
            elif all_or_first == "all":
                coordinates = [
                    (m.start() - bases_upstream, m.end() + bases_downstream)
                    for m in re.finditer(motif, self.basecall_sequence)
                ]
                extracted_motifs = []
                extracted_signals = []
                extracted_signals_means = []
                extracted_signals_trimmeans = []
                extracted_signals_dwells = []
                for coordinate in coordinates:
                    current_extracted_signal = []
                    current_extracted_signal_means = self.basecall_mean[
                        coordinate[0] : coordinate[1]
                    ]
                    current_extracted_signal_trimmeans = self.basecall_trimmean[
                        coordinate[0] : coordinate[1]
                    ]
                    current_extracted_signal_dwells = self.basecall_dwell[
                        coordinate[0] : coordinate[1]
                    ]
                    current_extracted_motif = self.basecall_sequence[
                        coordinate[0] : coordinate[1]
                    ]

                    for squiggle in self.basecall_resquiggled[
                        coordinate[0] : coordinate[1]
                    ]:
                        [
                            current_extracted_signal.append(measurement)
                            for measurement in squiggle
                        ]

                    if len(current_extracted_signal_means) == (
                        len(motif) + bases_upstream + bases_downstream
                    ):
                        extracted_motifs.append(current_extracted_motif)
                        extracted_signals.append(current_extracted_signal)
                        extracted_signals_means.append(current_extracted_signal_means)
                        extracted_signals_trimmeans.append(
                            current_extracted_signal_trimmeans
                        )
                        extracted_signals_dwells.append(current_extracted_signal_dwells)
                return (
                    extracted_motifs,
                    extracted_signals,
                    extracted_signals_means,
                    extracted_signals_trimmeans,
                    extracted_signals_dwells,
                )
        except KeyError:
            # print("Extraction of signal from basecall motif failed.\n")
            return None, None, None, None, None
        except AttributeError:
            return None, None, None, None, None

    def write_dataset(self, output_path: str):
        """
        Writes the dataset information to an HDF5 file.

        Parameters:
        ----------
        output_path : str
            The file path where the dataset will be written. If the file does not exist, it will be created.

        Writes:
        ------
        HDF5 File
            The dataset is stored in an HDF5 file with the following structure:
                - Each read ID is a group containing the following datasets:
                    - "read_id": The read ID.
                    - "basecall_sequence": Sequence from the basecall data.
                    - "basecall_mean": List of basecall mean values.
                    - "basecall_sd": List of basecall standard deviation values.
                    - "basecall_trimmean": List of trimmed mean values from the basecall.
                    - "basecall_trimsd": List of trimmed standard deviation values from the basecall.
                    - "basecall_dwell": List of basecall dwell times.
                    - "reference_sequence": Sequence from the reference data.
                    - "reference_mean": List of reference mean values.
                    - "reference_sd": List of reference standard deviation values.
                    - "reference_trimmean": List of trimmed mean values from the reference.
                    - "reference_trimsd": List of trimmed standard deviation values from the reference.
                    - "reference_dwell": List of reference dwell times.
                    - "index_ref_dict": Dictionary mapping index positions in the reference sequence to index positions in the read.

        Raises:
        ------
        AttributeError
            If an attribute is missing or cannot be accessed, the method will skip writing and pass silently.
        """
        try:
            dataset_dict = {}
            dataset_dict[self.id] = {
                "read_id": str(self.id),
                "basecall_sequence": str(self.basecall_sequence),
                "basecall_mean": list(self.basecall_mean),
                "basecall_sd": list(self.basecall_sd),
                "basecall_trimmean": list(self.basecall_trimmean),
                "basecall_trimsd": list(self.basecall_trimsd),
                "basecall_dwell": list(self.basecall_dwell),
                "reference_sequence": str(self.reference_sequence),
                "reference_mean": list(self.reference_mean),
                "reference_sd": list(self.reference_sd),
                "reference_trimmean": list(self.reference_trimmean),
                "reference_trimsd": list(self.reference_trimsd),
                "reference_dwell": list(self.reference_dwell),
                "index_ref_dict": [
                    (int(i), int(j)) for i, j in self.index_ref_dict.items()
                ],
            }
            if dataset_dict != None:
                if not os.path.exists(output_path):
                    with open(output_path, "w") as hf:
                        hf.write("")
                with h5py.File(output_path, "a") as hf:
                    for grp_name in dataset_dict:
                        grp = hf.create_group(grp_name)
                        for dset_name in dataset_dict[grp_name]:
                            dset = grp.create_dataset(
                                dset_name, data=dataset_dict[grp_name][dset_name]
                            )

        except AttributeError:
            pass

    def get_dataset(self):
        """Returns itself

        Returns:
            RawCurrentReadDataset
        """
        return self

    def get_dictionairy(self):
        """
        Generates a dictionary containing the dataset information for the current read.

        Returns:
        -------
        dict
            A dictionary with the dataset information structured as follows:
                - "read_id": The read ID.
                - "basecall_sequence": Sequence from the basecall data.
                - "basecall_mean": List of basecall mean values.
                - "basecall_sd": List of basecall standard deviation values.
                - "basecall_trimmean": List of trimmed mean values from the basecall.
                - "basecall_trimsd": List of trimmed standard deviation values from the basecall.
                - "basecall_dwell": List of basecall dwell times.
                - "reference_sequence": Sequence from the reference data.
                - "reference_mean": List of reference mean values.
                - "reference_sd": List of reference standard deviation values.
                - "reference_trimmean": List of trimmed mean values from the reference.
                - "reference_trimsd": List of trimmed standard deviation values from the reference.
                - "reference_dwell": List of reference dwell times.
                - "index_ref_dict": List of tuples mapping index positions in the reference sequence to index positions in the read.

        Raises:
        ------
        AttributeError
            If an attribute is missing or cannot be accessed, the method returns None.
        """
        try:
            dataset_dict = {}
            dataset_dict[self.id] = {
                "read_id": str(self.id),
                "basecall_sequence": str(self.basecall_sequence),
                "basecall_mean": list(self.basecall_mean),
                "basecall_sd": list(self.basecall_sd),
                "basecall_trimmean": list(self.basecall_trimmean),
                "basecall_trimsd": list(self.basecall_trimsd),
                "basecall_dwell": list(self.basecall_dwell),
                "reference_sequence": str(self.reference_sequence),
                "reference_mean": list(self.reference_mean),
                "reference_sd": list(self.reference_sd),
                "reference_trimmean": list(self.reference_trimmean),
                "reference_trimsd": list(self.reference_trimsd),
                "reference_dwell": list(self.reference_dwell),
                "index_ref_dict": [
                    (int(i), int(j)) for i, j in self.index_ref_dict.items()
                ],
            }
            return dataset_dict
        except AttributeError:
            return None
