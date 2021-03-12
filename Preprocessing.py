#!/usr/bin/env python
__author__ = "Yadollah Shahryary Dizaji"
__description__ = "Deep Learning Project"
__version__ = "1.0.0"
__email__ = "y.shahryary@tum.de"

# loading libraries

from gffutils.iterators import DataIterator
import os
import subprocess
import csv
import numpy as np
import re
import random
import glob
import pandas as pd
import re
import time
import h5py
import string


class Preprocess:

    ##############################################
    # General Functions
    ##############################################
    def replace_chr_text_in_file(self, filename):
        # Read in the file
        with open(filename, 'r') as file:
            file_content = file.read()
        # Replace the target string
        pattern = re.compile("chr", re.IGNORECASE)
        file_content = pattern.sub("", file_content)

        # Write the file out again
        with open(filename, 'w') as file:
            file.write(file_content)

    # Calculate file length
    def reg_file_len(self, filename):
        with open(filename) as f:
            for i, l in enumerate(f):
                pass
        return i + 1
    ##############################################
    # End General Functions
    ##############################################

    def create_output_genes_file(self):
        if len(glob.glob("*.gff")) > 0:
            input_file_gff = glob.glob("*.gff")[0]
        elif len(glob.glob("*.gff3")) > 0:
            input_file_gff = glob.glob("*.gff3")[0]
        else:
            print("no GFF file found")

        if not os.path.exists("output_genes.txt"):
            with open(str("output_genes.txt"), 'w') as fout:
                for feature in DataIterator(str(input_file_gff)):
                    if feature.featuretype == "gene":
                        fout.write(str(feature[0]) + '\t' + str(feature[3]) + '\t' + str(feature[4]) + '\t' +
                                   (feature[8]["gene_id"])[0] + '\n')
            self.replace_chr_text_in_file(str("output_genes.txt"))

        output_file_gff = "output_genes.txt"
        return input_file_gff, output_file_gff

    ##############################################
    # Calculate the lenght of the genome using the GFF file
    ##############################################
    #safe cast
    def safe_cast(self, val, to_type, default=None):
        try:
            return to_type(val)
        except (ValueError, TypeError):
            return default

    def genome_lenght(self, input_file_gff):
        seq_regions_maize_tmp = []
        list_order = []
        with open(input_file_gff) as f_original_gff:
            content_original_gff = f_original_gff.readlines()

        for hashtag_row in range(self.reg_file_len(input_file_gff)):
            if content_original_gff[hashtag_row].startswith('##'):
                if content_original_gff[hashtag_row].startswith('##sequence-region'):
                    temp_region = content_original_gff[hashtag_row].split()
                    if type(self.safe_cast(temp_region[1], int)) == int or str(temp_region[1]).startswith("chr"):
                        seq_regions_maize_tmp.append(self.safe_cast(temp_region[3], int))
                        list_order.append(self.safe_cast(temp_region[1], int) - 1)
            else:
                break

        seq_regions_maize = [None] * len(seq_regions_maize_tmp)
        for c, elem in enumerate(list_order):
            seq_regions_maize[elem] = seq_regions_maize_tmp[c]

        return sum(seq_regions_maize), seq_regions_maize



    ##############################################
    # process BAM (single file)
    ##############################################
    def process_bam_single(self, num_chr, bam_path):
        # index the RNA BAM file and create BED file for the BAM file
        BAM_list = glob.glob(bam_path)
        for bam_file in BAM_list:

            idxFile = str(bam_file) + ".bai"

            if not os.path.isfile(idxFile):
                subprocess.call(["samtools", "index", bam_file])

            dot_idx_reverse = (len(str(bam_file)) - str(bam_file).find('.')) - 1
            file_without_bam = str(bam_file)[:-dot_idx_reverse]
            bed_file_name = file_without_bam + "bed"

            if not os.path.isfile(bed_file_name):
                with open('temp_bash.sh', 'w') as the_file:
                    the_file.write("bamToBed < " + str(bam_file) + " > " + str(bed_file_name))
                os.system("sh './temp_bash.sh'")
                self.replace_chr_text_in_file(bed_file_name)

            dir_name = self.seperate_chromosomes(bed_file_name, num_chr)

        return dir_name

    ##############################################
    # seperate_chromosomes for RNA and GFF
    ##############################################
    def seperate_chromosomes(self, file_name, num_chr):
        file_ext = str(file_name)[-4:]
        filename_without_ext = str(file_name)[:-4]
        chr_dir_name = filename_without_ext + "_chromosomes"

        # create directory if not exists
        if not os.path.exists(chr_dir_name):
            try:
                os.mkdir(chr_dir_name)
            except OSError:
                print("Creation of the directory failed")

        for x in range(0, num_chr):
            if not os.path.isfile(chr_dir_name + "/" + str(x + 1) + file_ext):
                # print "test1", x
                with open('temp_bash.sh', 'w') as the_file:
                    the_file.write(
                        "awk -F, '$1 + 0 == " + str(x + 1) + "' " + file_name + " > " + chr_dir_name + "/" + str(
                            x + 1) + file_ext)
                os.system("sh './temp_bash.sh'")
        return chr_dir_name

    ##############################################
    # allc info:
    # https://github.com/yupenghe/methylpy#output-format
    # https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1541962
    ##############################################
    def process_methylated_data(self):
        allc_fname = ""
        for fname in os.listdir('.'):  # change directory as needed
            if 'allc' in fname:
                allc_fname = fname

        # create directory if not exists
        output_methylated_chr_dirName = "output_methylated_chr"
        if not os.path.exists(output_methylated_chr_dirName):
            try:
                os.mkdir(output_methylated_chr_dirName)
            except OSError:
                print("Creation of the directory failed -- output_methylated_chr_dirName")

        if len(allc_fname) > 0:
            if os.path.exists(allc_fname):
                for a in range(1, 11):
                    if os.path.exists(output_methylated_chr_dirName + "/" + str(a)):
                        pass
                    else:
                        with open('temp_bash.sh', 'w') as the_file:
                            the_file.write("awk -F\\\\t \'{print>\"output_methylated_chr/\"$1}\' " + allc_fname)
                        os.system("sh './temp_bash.sh'")
                        break

        return output_methylated_chr_dirName

    ##############################################
    # process excel file
    ##############################################
    def process_excel(self, HiChip_marks):
        list_filename_data = []
        list_filename_headers = []
        for h in HiChip_marks:
            filename_data = str(h + "_data.txt")
            filename_headers = str(h + "_headers.txt")
            if not os.path.isfile(filename_data) and not os.path.isfile(filename_headers):
                with open(h+".csv") as f:
                    reader = csv.reader(f, delimiter="\t")
                    d = list(reader)

                head = d[0:1]
                d = d[1:]
                d = np.array(d, dtype='float')

                np.savetxt(filename_data, d, fmt='%i\t%i\t%i\t%i\t%i')

                df_headers = np.array(list(head))
                np.savetxt(filename_headers, df_headers, fmt='%s\t%s\t%s\t%s\t%s')

            list_filename_data.append(filename_data)
            list_filename_headers.append(filename_headers)

        return list_filename_data, list_filename_headers
    ##############################################
    # HiChip_marks_data
    ##############################################
    def HiChip_marks_data(self, filename_data):

        final_matrix = []
        for fd in filename_data:
            with open(fd) as f:
                reader = csv.reader(f, delimiter="\t")
                d = list(reader)
            d = np.array(d, dtype='float')
            final_matrix.extend(d)
        final_matrix = np.array(final_matrix)
        # remove rows where regions have exactly same start and stop coordinates for both s1/e1 and s2/e2
        final_matrix = final_matrix[np.unique(final_matrix[:, [0, 1, 2, 3, 4]], return_index=True, axis=0)[1]]

        final_matrix = np.array(final_matrix, dtype='object')
        final_matrix[:, 0:6] = final_matrix[:, 0:6].astype(int)
        #final_matrix[:, 6] = final_matrix[:, 6].astype(float)
        
        final_matrix=final_matrix[np.argsort(final_matrix[:, 0])]
        # label with all positives
        pos_labels = np.full(len(final_matrix), 1)
        pos_labels = pos_labels.reshape(len(pos_labels), -1)
        final_matrix = np.append(final_matrix, pos_labels, axis=1)
        print(final_matrix[1:4])
        # difference in base pairs between region 1 and region 2
        diff_r1_r2 = final_matrix[:, 3] - final_matrix[:, 2]

        print(np.mean(diff_r1_r2), "mean difference in base pairs between region 1 and region 2 - interactions data")
        print(min(diff_r1_r2), "min difference in base pairs between region 1 and region 2 - interactions data")
        print(max(diff_r1_r2), "max difference in base pairs between region 1 and region 2 - interactions data")

        return final_matrix

    ##############################################


    ##############################################
    # round functions
    ##############################################
    def round_down(self, num, divisor=5000):
        return num - (num % divisor)

    def round_nearest(self, x, num=10000):
        return int(round(float(x) / num) * num)

    ##############################################
    # random_gen_negative_samples_matrix
    ##############################################
    def random_gen_negative_samples_matrix(self, HiChip_matrix, seq_regions_maize):
        negative_samples = []
        randomly_selected_chr_list = []
        toSample_start_idx_region_1 = []
        if not os.path.exists("random_gen_negative_samples_matrix_2500.txt"):
            for neg_sample in range(len(HiChip_matrix)):
                # for neg_sample in range(25):
                if neg_sample % 50 == 0:
                    print(neg_sample, "neg_sample number -- ln451")
                # randomly_selected_chr = random.randint(0,len(seq_regions_maize)-1)
                # randomly_selected_chr_list.append(randomly_selected_chr+1)
                while True:
                    randomly_selected_chr = random.randint(0, len(seq_regions_maize) - 1)

                    toSample = int(subprocess.check_output(
                        ["shuf", "-i", "1-" + str(seq_regions_maize[randomly_selected_chr] - 45001), "-n", "1"]))
                    toSample = self.round_nearest(toSample)

                    check1 = False
                    for rand_chr, t in zip(randomly_selected_chr_list, toSample_start_idx_region_1):
                        if randomly_selected_chr + 1 == rand_chr and t == toSample:
                            check1 = True
                            break

                    check2 = False
                    if not check1:
                        for rand_chr, t in zip(HiChip_matrix[:, 0], HiChip_matrix[:, 1]):
                            if randomly_selected_chr + 1 == rand_chr and t == toSample:
                                check2 = True
                                break

                    if not check1 and not check2:
                        randomly_selected_chr_list.append(randomly_selected_chr + 1)
                        toSample_start_idx_region_1.append(toSample)
                        break
            np.savetxt('random_gen_negative_samples_matrix.txt',
                       np.c_[randomly_selected_chr_list, toSample_start_idx_region_1], fmt='%i', delimiter='\t')

        print("negative samples generated")
        #neg_samples_matrix = np.loadtxt('random_gen_negative_samples_matrix.txt')
        #randomly_selected_chr_list = neg_samples_matrix[:, 0]
        #toSample_start_idx_region_1 = neg_samples_matrix[:, 1]
        neg_samples_matrix = np.loadtxt('random_gen_negative_samples_matrix_2500.txt',dtype=int)
        negative_samples=neg_samples_matrix
        #print(negative_samples[1:5])
        
        '''
        for rand_chr, t in zip(randomly_selected_chr_list, toSample_start_idx_region_1):
            rand_chr = int(rand_chr)
            t = int(t)
            temp_neg_sample= [rand_chr, t, t + 1000, t + 29000, t + 30000, 0]
            negative_samples.append(temp_neg_sample)
            # print rand_chr, t
        # create negatives set
        negative_samples = np.array(negative_samples)
        '''
        print(negative_samples[5:10])
        return negative_samples

    ##############################################
    # Combine negative samples with positives and shuffle
    # Train & Test Samples
    ##############################################
    def combine_neg_pos(self, HiChip_marks_data_matrix, negative_samples):
        HiChip_marks_data_matrix_with_neg = np.concatenate((HiChip_marks_data_matrix, negative_samples))
        # shuffle the negatives and positives
        np.random.shuffle(HiChip_marks_data_matrix_with_neg)
        print("len Shuffle:", len(HiChip_marks_data_matrix_with_neg))
        return HiChip_marks_data_matrix_with_neg

    ##############################################
    # process_FASTA_vector
    ##############################################
    def process_FASTA_vector(self, region_s, region_e, chr_num_in_line, retrieve_FASTA_matrix):

        FASTA_sequences_filename = str(region_s) + "_" + str(region_e) + ".txt"

        if not os.path.exists("FASTA_sequences/" + chr_num_in_line + "/" + FASTA_sequences_filename):
            # create FASTA per_base_cov directory if not exists
            if not os.path.exists("FASTA_sequences"):
                try:
                    os.mkdir("FASTA_sequences")
                except OSError:
                    print("Creation of the directory failed -- FASTA_sequences")

            if not os.path.exists("FASTA_sequences/" + chr_num_in_line):
                try:
                    os.mkdir("FASTA_sequences/" + chr_num_in_line)
                except OSError:
                    print("Creation of the directory failed -- " + "FASTA_sequences/" + chr_num_in_line)

            FASTA_sequences_filename = str(region_s) + "_" + str(region_e) + ".txt"

            if not os.path.exists("FASTA_sequences/" + chr_num_in_line + "/" + FASTA_sequences_filename):
                with open('temp_bash.sh', 'w') as the_file:
                    the_file.write(
                        "bedtools getfasta -fi Zea_mays_B73_v4.fasta -bed temp_region.txt > FASTA_sequences/" + chr_num_in_line + "/" + FASTA_sequences_filename)
                subprocess.call(["sh", "./temp_bash.sh"])

        if retrieve_FASTA_matrix and os.path.exists(
                "FASTA_sequences/" + chr_num_in_line + "/" + FASTA_sequences_filename):
            FASTA_matrix = np.zeros((4, 2500))  # ACGT
            fp = open("FASTA_sequences/" + chr_num_in_line + "/" + FASTA_sequences_filename)
            for i, line in enumerate(fp):
                if i == 1:
                    line_character_count = 0
                    for c in line:
                        # print c
                        if c == "A":
                            FASTA_matrix[0][line_character_count] = 1
                        if c == "C":
                            FASTA_matrix[1][line_character_count] = 1
                        if c == "G":
                            FASTA_matrix[2][line_character_count] = 1
                        if c == "T":
                            FASTA_matrix[3][line_character_count] = 1
                        line_character_count += 1
            fp.close()

            return FASTA_matrix





    ##############################################
    # create Matrix
    ##############################################
    def create_matrix(self, chr_num, region_s, region_e, using_fasta_seq):

        with open('temp_region.txt', 'w') as the_file:
            the_file.write(str(chr_num) + "\t" + str(region_s) + "\t" + str(region_e))

        print("subprocess call will start")

        # empty the temp_bash.sh file
        open('temp_bash.sh', 'w').close()

        subprocess_bedtools_job_tracker_counter = 0

        if not os.stat("temp_bash.sh").st_size == 0:
            for subprocess_count in range(subprocess_bedtools_job_tracker_counter):
                with open('temp_bash.sh', 'a') as the_file:
                    letter = string.ascii_lowercase[subprocess_count]
                    the_file.write("wait $" + letter + "\n" "echo \"job " + letter + " returned $?\" \n")
                    the_file.close()

        # subprocess call
        subprocess.call(["sh", "./temp_bash.sh"])
        print("subprocess call ended")

        # generate sequence
        if using_fasta_seq:
            self.process_FASTA_vector(region_s, region_e, chr_num, retrieve_FASTA_matrix=False)

        final_matrix = []

        final_matrix = np.array(final_matrix, dtype=int)

        # get ACGT sequence
        if using_fasta_seq:
            FASTA_matrix = self.process_FASTA_vector(region_s, region_e, chr_num, retrieve_FASTA_matrix=True)
            if len(final_matrix) == 0:
                final_matrix = FASTA_matrix
            else:
                final_matrix = np.array(final_matrix, dtype=int)
                final_matrix = np.vstack((final_matrix, FASTA_matrix))

        final_matrix = np.array(final_matrix, dtype=int)
        print(final_matrix)
        print(final_matrix.shape, "final_matrix shape")
        return final_matrix


    ##############################################
    # Create hdf5 file of the data (sequence information and labels)
    ##############################################
    def data_to_h5(self, HiChip_marks_data_matrix_with_neg,
                                           total_num_seq_signals, HDF5_samples_num,
                                           using_fasta_seq):

        h5_filename = "regions_labels_rand_neg_" + str(HDF5_samples_num) + ".h5"

        if not os.path.isfile(h5_filename):
            hf = h5py.File(h5_filename, 'w')
            hf.create_dataset('re_region_matrices', shape=(HDF5_samples_num, total_num_seq_signals, 2500), dtype="int",
                              data=None)
            hf.create_dataset('interactions_region_matrices', shape=(HDF5_samples_num, total_num_seq_signals, 2500),
                              dtype="int", data=None)

            # get labels
            print("retrieving labels")
            labels = HiChip_marks_data_matrix_with_neg[:, -1]
            labels = labels.astype(int)
            hf.create_dataset('labels', data=labels[:HDF5_samples_num])

            # get all regulatory element region matrices
            print("retrieving regulatory element region matrices")
            # re_region_matrices = []
            for data_count, data_sample in enumerate(HiChip_marks_data_matrix_with_neg[:HDF5_samples_num]):
                print(data_count, " retrieving regulatory element region matrices")
                chr_num = data_sample[0]
                re_region_s = data_sample[1]
                re_region_e = data_sample[2]

                re_region_matrix = self.create_matrix(str(chr_num), re_region_s, re_region_e,using_fasta_seq)

                with h5py.File(h5_filename, 'a') as hf:
                    hf["re_region_matrices"][data_count] = re_region_matrix

            # get all interactions region matrices
            print("retrieving interactions region matrices")
            # interactions_region_matrices = []
            for data_count, data_sample in enumerate(HiChip_marks_data_matrix_with_neg[:HDF5_samples_num]):
                print(data_count, " retrieving interactions region matrices")
                chr_num = data_sample[0]
                in_region_s = data_sample[3]
                in_region_e = data_sample[4]

                interactions_region_matrix = self.create_matrix(str(chr_num), in_region_s, in_region_e,using_fasta_seq)
                with h5py.File(h5_filename, 'a') as hf:
                    hf["interactions_region_matrices"][data_count] = interactions_region_matrix

            hf.close()

            print("done creating HDF5 file")

        return h5_filename


    ############################################################################################
    # Main run
    ############################################################################################
    def __init__(self, data_dir, HiChip_marks_list,
                 num_seq_signals_excluding_histones, HDF5_samples_num,
                 using_fasta_seq):
        os.chdir(data_dir)
        # data size
        # Create gene regions file (includes indices of genes corresponding chromosome number)
        self.input_file_gff, self.output_file_gff = self.create_output_genes_file()
        print("Reference Genome: ", self.input_file_gff, "Output file: ", self.output_file_gff)
        self.genome_lenght, self.seq_regions_maize = self.genome_lenght(self.input_file_gff)


        # get genome length and Chromosome lengths
        print(self.genome_lenght, 'is the genome lenght')
        print(self.seq_regions_maize, 'are the lenghts of each chromosome by order')
        print(len(self.seq_regions_maize), 'are the total number of chromosomes')


        # process excel
        self.list_filename_data, self.list_filename_headers = self.process_excel(HiChip_marks_list)

        # get HiChip_marks_data
        self.HiChip_marks_data_matrix = self.HiChip_marks_data(self.list_filename_data)
        
        
        # get negative samples
        self.negative_samples = self.random_gen_negative_samples_matrix(self.HiChip_marks_data_matrix, self.seq_regions_maize)

        # combine negative samples with positives and shuffle
        self.HiChip_marks_data_matrix_with_neg = self.combine_neg_pos(self.HiChip_marks_data_matrix, self.negative_samples)

        # create Histone samples of the data

        self.total_num_seq_signals = num_seq_signals_excluding_histones

        # saving to hdf5 file
        self.h5_filename = self.data_to_h5(self.HiChip_marks_data_matrix_with_neg,
                                           self.total_num_seq_signals, HDF5_samples_num,
                                           using_fasta_seq )

