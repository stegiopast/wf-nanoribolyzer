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

"""
Bascalling and alignment should have been executed with the following flags in order for this script to function

dorado basecaller ~/dorado-0.4.3-linux-x64/bin/rna002_70bps_hac@v3 extraction_final.pod5 --emit-moves |\
      ~/samtools-1.18/bin/samtools fastq -T "*" | minimap2 -y --MD -ax map-ont Gm_randomers.fasta - |\
      ~/samtools-1.18/bin/samtools view -b -F 4 > extraction_final.bam

Move table and MD tags must be integrated in the final bam file
"""

"""
Creates an object of a single read in a pod5 dataset, providing direct access to resquiggleing data performed with ONT's remora tool. 
"""


class RawCurrentReadDataset:
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

    """
    Iterate over resquiggled data and
    reconstruct the raw_current at a given coordinate +/-
    upstream and downstream bases
    Also extract the mean signal for each squiggle
    """

    def extract_signal_reference_coordinates(
        self, coordinate: int = 0, bases_upstream: int = 4, bases_downstream: int = 4
    ):
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

    """
    In 'first' only the first occurence of the motif in the reference is extracted
    In 'all' all the occurences of the motif in the reference are extracted
    Signal is extracted from position of substring in the reference +/- up- and downstream bases buffer
    The extracted motif will be returned
    The signal will be returned as a list of raw current and as a list with means for each squiggle
    """

    def extract_signal_reference_motif(
        self,
        motif: str = "AG",
        all_or_first: str = "all",
        bases_upstream: int = 0,
        bases_downstream: int = 0,
    ):
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

    """
    In 'first' only the first occurence of the motif in the basecall sequence is extracted
    In 'all' all the occurences of the motif in the basecall sequence are extracted
    Signal is extracted from position of substring in the basecall +/- up- and downstream bases buffer
    The extracted motif will be returned
    The signal will be returned as a list of raw current and as a list with means for each squiggle
    """

    def extract_signal_basecall_motif(
        self,
        motif: str = "AG",
        all_or_first: str = "all",
        bases_upstream: int = 0,
        bases_downstream: int = 0,
    ):
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
        return self

    def get_dictionairy(self):
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
