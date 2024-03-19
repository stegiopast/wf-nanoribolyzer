from pathlib import Path
import pod5
import remora
from remora import io,refine_signal_map, util
import numpy as np
import argparse
import pysam 
from tqdm import tqdm
import math
from itertools import repeat
import json
import re
import pandas as pd

"""
Bascalling and alignment should have been executed with the following flags in order for this script to function

dorado basecaller ~/dorado-0.4.3-linux-x64/bin/rna002_70bps_hac@v3 extraction_final.pod5 --emit-moves |\
      ~/samtools-1.18/bin/samtools fastq -T "*" | minimap2 -y --MD -ax map-ont Gm_randomers.fasta - |\
      ~/samtools-1.18/bin/samtools view -b -F 4 > extraction_final.bam

Move table and MD tags must be integrated in the final bam file
"""

class RawCurrentReadDataset():
    def __init__(self,id: str="" ,pod5_filename: Path="" ,bam_filename: Path="", pod5_dr:pod5.dataset.DatasetReader=None, bam_fh:remora.io.ReadIndexedBam=None, kmer_table:Path=""):
        try:
            self.id = id
            self.kmer_table = kmer_table
            self.sig_map_refiner = refine_signal_map.SigMapRefiner(kmer_model_filename=self.kmer_table,do_rough_rescale=True,scale_iters=1,do_fix_guage=True)       
            if pod5_dr is None:
                self.pod5_dr = pod5.DatasetReader(pod5_filename)
            else:
                self.pod5_dr = pod5_dr
            if bam_fh is None:
                self.bam_fh = io.ReadIndexedBam(bam_filename)
            else:
                self.bam_fh = bam_fh
            #print(type(self.bam_fh)
            self.bam_read = self.bam_fh.get_first_alignment(id)
            if self.bam_read.is_forward:
                flip = True
            else:
                flip = False
            self.pod5_read = self.pod5_dr.get_read(id)
            self.io_read_basecall = io.Read.from_pod5_and_alignment(self.pod5_read, self.bam_read,reverse_signal=flip)
            self.io_read_basecall.set_refine_signal_mapping(self.sig_map_refiner, ref_mapping=False)
            #self.basecall_metrics = self.io_read_basecall.compute_per_base_metric
            #print(self.basecall_metrics)
            
            #Extract basecall anchored information
            self.basecall_anchored_data = self.io_read_basecall.extract_basecall_region(end_base=self.io_read_basecall.seq_len)
            self.basecall_sequence = str(self.basecall_anchored_data.seq)
            self.basecall_sequence_to_signal_map = self.basecall_anchored_data.seq_to_sig_map
            self.basecall_norm_signal = self.basecall_anchored_data.norm_signal
            self.basecall_resquiggled = []
            for index,value in enumerate(self.basecall_sequence):
                associated_signal = self.basecall_anchored_data.norm_signal[self.basecall_anchored_data.seq_to_sig_map[index]:self.basecall_anchored_data.seq_to_sig_map[index + 1]]
                self.basecall_resquiggled.append(list(associated_signal))
            
            #Extract reference anchored information
            self.io_read_reference = io.Read.from_pod5_and_alignment(self.pod5_read, self.bam_read,reverse_signal=True)
            self.io_read_reference.set_refine_signal_mapping(self.sig_map_refiner, ref_mapping=False)
            self.io_read_reference.set_refine_signal_mapping(self.sig_map_refiner, ref_mapping=True)
            #self.reference_metrics = self.io_read_reference.compute_per_base_metric
            self.reference_anchored_data = self.io_read_reference.extract_ref_reg(self.io_read_reference.ref_reg)
            self.reference_sequence = str(self.reference_anchored_data.seq)
            self.contig = self.reference_anchored_data.ref_reg.ctg
            #Maps read to reference 
            index_reference = [i for i in range(self.reference_anchored_data.ref_reg.start,self.reference_anchored_data.ref_reg.end)]
            self.index_ref_dict = {}
            for index,value in enumerate(self.reference_sequence):
                self.index_ref_dict[index_reference[index]] = index
            self.reference_resquiggled = []
            self.reference_norm_signal = self.reference_anchored_data.norm_signal
            for index,value in enumerate(self.reference_sequence):
                associated_signal = self.reference_norm_signal[self.reference_anchored_data.seq_to_sig_map[index]:self.reference_anchored_data.seq_to_sig_map[index + 1]]
                self.reference_resquiggled.append(list(associated_signal))
        except:
            #print("Error in extracting reads")
            return
    
    def extract_signal_reference_coordinates(self,coordinate:int = 0,bases_upstream:int  = 4,bases_downstream:int = 4):
        #Iterate over resquiggled data and 
        #reconstruct the raw_current at a given coordinate +/- 
        #upstream and downstream bases
        #Also extract the mean signal for each squiggle
        try:
            extracted_motif = self.reference_sequence[self.index_ref_dict[coordinate]-bases_upstream: self.index_ref_dict[coordinate]+(bases_downstream+1)]
            extracted_signal = []
            extracted_signal_means = []
            for squiggle in self.reference_resquiggled[self.index_ref_dict[coordinate]-bases_upstream: self.index_ref_dict[coordinate]+(bases_downstream+1)]:
                for measurement in squiggle:
                    extracted_signal.append(measurement)
                #print(squiggle)
                if squiggle is None or squiggle == []:
                    extracted_signal_means.append(100)
                elif np.array(squiggle).size >= 2:
                    #print(np.array(squiggle))
                    extracted_signal_means.append(np.array(squiggle).mean())
                elif np.array(squiggle).size == 1:
                    extracted_signal_means.append(squiggle[0])
                #print(extracted_signal_means)  
            #print(extracted_signal_means) 
            if 100 in extracted_signal_means:
                return None,None,None
            #print(extracted_signal_means)
            if len(extracted_signal_means) == (1 + bases_upstream + bases_downstream): 
                return extracted_motif,extracted_signal, extracted_signal_means
            else:
                return None,None,None

        except KeyError:
            #print("Extraction of signal from reference coordinate failed.\n")
            return None,None,None
        except AttributeError:
            return None,None,None

    def extract_signal_reference_motif(self,motif:str="AG",all_or_first:str="all",bases_upstream:int=0,bases_downstream:int=0):
        #In 'first' only the first occurence of the motif in the reference is extracted
        #In 'all' all the occurences of the motif in the reference are extracted 
        #Signal is extracted from position of substring in the reference +/- up- and downstream bases buffer
        #The extracted motif will be returned
        #The signal will be returned as a list of raw current and as a list with means for each squiggle  
        try:
            if all_or_first == "first":
                coordinate = self.reference_sequence.find(motif)
                length = len(motif)
                extracted_motif = self.reference_sequence[coordinate-bases_upstream:coordinate+bases_downstream+length+1]
                extracted_signal = []
                extracted_signal_means = []
                for squiggle in self.reference_resquiggled[coordinate-bases_upstream:coordinate+bases_downstream+length+1]:
                    for measurement in squiggle:
                        extracted_signal.append(measurement)
                    if squiggle is None or squiggle == []:
                        extracted_signal_means.append(100)
                    elif np.array(squiggle).size >= 2:
                        #print(np.array(squiggle))
                        extracted_signal_means.append(np.array(squiggle).mean())
                    elif np.array(squiggle).size == 1:
                        extracted_signal_means.append(squiggle[0])
                if 100 in extracted_signal_means:
                    return None,None,None
                    #print(extracted_signal_means)
                if len(extracted_signal_means) == (len(motif) + bases_upstream + bases_downstream): 
                    return extracted_motif,extracted_signal, extracted_signal_means
                else:
                    return None,None,None
            elif all_or_first == "all":
                coordinates = [(m.start()-bases_upstream,m.end()+bases_downstream) for m in re.finditer(motif, self.reference_sequence)]
                extracted_motifs = []
                extracted_signals = []
                extracted_signals_means = []
                for coordinate in coordinates:
                    current_extracted_signal = []
                    current_extracted_signal_means = []
                    current_extracted_motif = self.reference_sequence[coordinate[0]:coordinate[1]]
                    for squiggle in self.reference_resquiggled[coordinate[0]:coordinate[1]]:
                        [current_extracted_signal.append(measurement) for measurement in squiggle]   
                        if squiggle is None or squiggle == []:
                            current_extracted_signal_means.append(100)
                        elif np.array(squiggle).size >= 2:
                            current_extracted_signal_means.append(np.array(squiggle).mean())
                        elif np.array(squiggle).size == 1:
                             current_extracted_signal_means.append(squiggle[0])
                    if 100 in current_extracted_signal_means:
                        continue
                    if len(current_extracted_signal_means) == (len(motif) + bases_upstream + bases_downstream): 
                        extracted_motifs.append(current_extracted_motif)
                        extracted_signals.append(current_extracted_signal)
                        extracted_signals_means.append(current_extracted_signal_means)
                return extracted_motifs,extracted_signals,extracted_signals_means
        except KeyError:
            #print("Extraction of signal from reference motif failed.\n")
            return None,None,None
        except AttributeError:
            #print("Extraction of signal from reference motif failed.\n")
            return None,None,None
        
    def extract_signal_basecall_motif(self,motif:str="AG",all_or_first:str="all",bases_upstream:int=0,bases_downstream:int=0):
        #In 'first' only the first occurence of the motif in the basecall sequence is extracted
        #In 'all' all the occurences of the motif in the basecall sequence are extracted 
        #Signal is extracted from position of substring in the basecall +/- up- and downstream bases buffer
        #The extracted motif will be returned
        #The signal will be returned as a list of raw current and as a list with means for each squiggle  
        #try:
        if all_or_first == "first":
            coordinate = self.basecall_sequence.find(motif)
            length = len(motif)
            extracted_motif = self.basecall_sequence[coordinate-bases_upstream:coordinate+bases_downstream+length+1]
            extracted_signal = []
            extracted_signal_means = []
            for squiggle in self.basecall_resquiggled[coordinate-bases_upstream:coordinate+bases_downstream+length+1]:
                for measurement in squiggle:
                    extracted_signal.append(measurement)
                if squiggle is None or squiggle == []:
                    extracted_signal_means.append(100)
                elif np.array(squiggle).size >= 2:
                    #print(np.array(squiggle))
                    extracted_signal_means.append(np.array(squiggle).mean())
                elif np.array(squiggle).size == 1:
                    extracted_signal_means.append(squiggle[0])
            if 100 in extracted_signal_means:
                return None,None,None
                #print(extracted_signal_means)
            if len(extracted_signal_means) == (len(motif) + bases_upstream + bases_downstream): 
                return extracted_motif,extracted_signal, extracted_signal_means
            else:
                return None,None,None
        elif all_or_first == "all":
            coordinates = [(m.start()-bases_upstream,m.end()+bases_downstream) for m in re.finditer(motif, self.basecall_sequence)]
            extracted_motifs = []
            extracted_signals = []
            extracted_signals_means = []
            for coordinate in coordinates:
                current_extracted_signal = []
                current_extracted_signal_means = []
                current_extracted_motif = self.basecall_sequence[coordinate[0]:coordinate[1]]
                for squiggle in self.basecall_resquiggled[coordinate[0]:coordinate[1]]:
                    [current_extracted_signal.append(measurement) for measurement in squiggle]   
                    if squiggle is None or squiggle == []:
                        current_extracted_signal_means.append(100)
                    elif np.array(squiggle).size >= 2:
                        current_extracted_signal_means.append(np.array(squiggle).mean())
                    elif np.array(squiggle).size == 1:
                            current_extracted_signal_means.append(squiggle[0])
                if 100 in current_extracted_signal_means:
                    continue
                if len(current_extracted_signal_means) == (len(motif) + bases_upstream + bases_downstream): 
                    extracted_motifs.append(current_extracted_motif)
                    extracted_signals.append(current_extracted_signal)
                    extracted_signals_means.append(current_extracted_signal_means)
            return extracted_motifs,extracted_signals,extracted_signals_means
        #except:
        #    print("Extraction of signal from basecall motif failed.\n")
        #    return None,None,None


        


if __name__ =='__main__':
    import matplotlib.pyplot as plt
    pod5_dr = pod5.DatasetReader("/home/stefan/Data_nano_ribolyzer/GM_randomers/randomer_extraction_final.pod5")
    read_ids = [str(read_record.read_id) for read_record in pod5_dr.reads()]

    bam_fh = io.ReadIndexedBam("/home/stefan/Data_nano_ribolyzer/GM_randomers/randomer_extraction_final.bam")
    kmer_table = "/home/stefan/kmer_models/rna_r9.4_180mv_70bps/5mer_levels_v1.txt"
    mean_signals = []
    fig,(ax,ax2,ax3,ax4,ax5,ax6,ax7) = plt.subplots(nrows=7,ncols=1, figsize=(14,14))
    mean_signals1 = []
    mean_signals2 = []
    for id in tqdm(read_ids[0:1000]):#['00677d43-1ae3-474f-9d03-0c7f97d3b495', '0060c229-4802-43c4-8524-e64a3179ceac', '00b0ff84-2604-4ab3-b455-0b8cb74ff8a2', '00c2dc01-e935-4400-a08e-d2e4749d14d8', '00e3792d-1516-4bf2-b427-2cc500f8eb95', '01300f5b-7b38-4d20-b4f1-6161ea104f5e', '01d8b539-f49b-400a-90e7-696f5d4baf13', '05598ed6-3dfb-465d-a12c-e735e3cba960']:#, '057055c5-586d-469f-9c99-5d1df569d45b', '0c04c630-bb5e-479c-8ae8-cdcd57b46f45', '0d87dd21-af19-446d-b8d4-2b91a74d5ff6', '00bf2f25-e186-44ab-ae39-444bd0f69e42', '00c5d4ed-f619-409a-b27a-0b11fed0c574', '174f5127-f4a0-4a8f-9d0e-d9af04cce4c3', '0003c7cf-fb1e-422d-858d-15ca20fca2ca', '0051edf5-8cff-4aab-90e9-9f8e5cf7c202', '0059a274-fc67-42cf-aeab-6de6a6a382a4', '0052a214-742b-4f75-9196-bf25af1ee4ca', '0075832b-a9e7-4901-ae40-c4dd32ebdadd', '0094b3f5-8d3b-4ec9-8b48-a191ccf28448', '00bdeccb-e13f-4a54-9944-b5656fb44647', '00df1698-0cdb-4071-9db4-bf3b93a96b23', '00eca5aa-fd58-4866-940a-eaf86d91dd47', '00cba7e1-161e-429e-b3af-d7f02d44bdbe', '01000291-082f-4841-ad63-a23bd7035321', '01328f13-252d-43b2-ad22-fb7e333c0137']:
    #for id in tqdm(['00b0ff84-2604-4ab3-b455-0b8cb74ff8a2' for i in range(100)]):
        #print(f"\n\n\n{id}")
        dataset = RawCurrentReadDataset(id = id,
                                        #pod5_filename="/home/stefan/Data_nano_ribolyzer/GM_randomers/randomer_extraction_final.pod5",
                                        #bam_filename="/home/stefan/Data_nano_ribolyzer/GM_randomers/randomer_extraction_final.bam",
                                        pod5_dr = pod5_dr,
                                        bam_fh = bam_fh,
                                        kmer_table=kmer_table)
        #print(dataset.basecall_norm_signal)
        #print("\n")
        #print(dataset.reference_norm_signal)
        #print(dataset.reference_anchored_data.seq_to_sig_map)
        #print("\n")
        #plt.plot([i for i in range(len(dataset.basecall_norm_signal))],dataset.basecall_norm_signal, color = "blue", marker = "o")
        #plt.scatter([i for i in range(len(dataset.reference_norm_signal))],dataset.reference_norm_signal, color = "red", marker = "o")


        motif, signal ,mean_signal = dataset.extract_signal_reference_coordinates(637,20,20)
        if motif is None or signal is None or mean_signal is None:
            continue

        if mean_signal is not None or mean_signal != []:
            ax.plot([j for j in range(len(mean_signal))],mean_signal,color = "blue",marker="o",markerfacecolor='green', markersize=6, linewidth=0.22,alpha = 0.2)

        
        if mean_signal is not None or mean_signal != []:
            mean_signals1.append(mean_signal)
            mean_of_means = np.mean(np.array(mean_signals1),axis=0)
            std_of_means = np.std(np.array(mean_signals1),axis=0)
            

        motif, signal ,mean_signal = dataset.extract_signal_reference_coordinates(700,20,20)
        if motif is None or signal is None or mean_signal is None:
             continue

        if mean_signal is not None or mean_signal != []:
            mean_signals2.append(mean_signal)
            mean_of_means2 = np.mean(np.array(mean_signals2),axis=0)
            std_of_means2 = np.std(np.array(mean_signals2),axis=0)

        if mean_signal is not None or mean_signal != []:
            ax2.plot([j for j in range(len(mean_signal))],mean_signal,color = "red",marker="o",markerfacecolor='green', markersize=6, linewidth=0.2,alpha=0.2)
        

        ref_motifs, ref_signals, ref_mean_signals = dataset.extract_signal_reference_motif(motif = "AATATGGCCANNNNNNNNNNGNNNNNNNNNNGCGAGGAGCT", all_or_first="all",bases_upstream=0,bases_downstream=0)
        if ref_motifs is None or ref_signals is None or ref_mean_signals is None:
            continue
        # print("\nExtraction Motif in reference")
        # [print(i,j) for i,j in zip(ref_motifs,ref_mean_signals)]
        for i in ref_mean_signals:
              if i is not None:
                ax5.plot([j for j in range(len(i))],i,color="blue",marker="o",markerfacecolor='green', markersize=6, linewidth=0.2, alpha = 0.2)

                
        ref_motifs, ref_signals, ref_mean_signals = dataset.extract_signal_reference_motif(motif = "CTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAG", all_or_first="all",bases_upstream=0,bases_downstream=0)
        if ref_motifs is None or ref_signals is None or ref_mean_signals is None:
            continue

        for i in ref_mean_signals:
              if i is not None:
                ax6.plot([j for j in range(len(i))],i,color="red",marker="o",markerfacecolor='green', markersize=6, linewidth=0.2, alpha= 0.2)

        bcall_motifs, bcall_signals, bcall_mean_signals = dataset.extract_signal_basecall_motif(motif = "CTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCA", all_or_first="all",bases_upstream=0,bases_downstream=0)
        if bcall_motifs is None or bcall_signals is None or bcall_mean_signals is None:
            continue
        #print("\nExtraction Motif in basecalls")
        #[print(i,j) for i,j in zip(bcall_motifs,bcall_mean_signals)]
        for i in bcall_mean_signals:
            if i is not None:
                ax7.plot([j for j in range(len(i))],i,color="red",marker="o",markerfacecolor='green', markersize=6, linewidth=0.2, alpha = 0.2)
    
    ax3.plot([j for j in range(len(mean_signals1[0]))],mean_of_means,color = "blue")
    ax3.fill_between([j for j in range(len(mean_signals1[0]))],mean_of_means - std_of_means,mean_of_means + std_of_means, color = "blue", alpha=0.2)
    ax4.plot([j for j in range(len(mean_signals1[0]))],mean_of_means2,color = "red")
    ax4.fill_between([j for j in range(len(mean_signals1[0]))],mean_of_means2 - std_of_means2,mean_of_means2 + std_of_means2, color = "red", alpha=0.2)
    
    ax.set_ylim([-2,2])
    ax2.set_ylim([-2,2])
    ax3.set_ylim([-2,2])
    ax4.set_ylim([-2,2])
    ax5.set_ylim([-2,2])
    ax6.set_ylim([-2,2])
    plt.show()




    


        

        


