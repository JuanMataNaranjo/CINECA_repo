from abc import ABCMeta, abstractmethod
import re
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from collections import defaultdict
import subprocess
import glob


class LogMain(metaclass=ABCMeta):
    """
    Abstract object
    All classes that will inherit from this class will need to include the following methods
    """

    @abstractmethod
    def read_log(self, *args, **kwargs):
        """
        Method to read log and provide as a file type which is most according depending on the .log file
        """
        raise NotImplementedError()

    @abstractmethod
    def check_log(self, *args, **kwargs):
        """
        Perform checks on the log. Ideally each class will have other methods which will help this method out
        """
        raise NotImplementedError()

    # TODO: Look at the format of the initial champiaons (name consistency)


class Overall:
    """
    Different to the next classes, this class will look at the baseline samples to check that things are in order,
    i.e. this class will not use the logs
    Due to the large size of files, we will directly call the functions using the terminal, without reading data into
    python
    """

    def __init__(self, sample, path, table_path='data/fastq.csv'):
        self.sample = sample
        self.path = path
        self.paired = self.single_paired(table_path)

    def single_paired(self, table_path='data/fastq.csv'):
        """
        Check whether the sample is paired (R1 and R2) or single (only R1). This method is relevant
        for cases in which two different log files are generated
        :param table_path: Path in which we can find the fastq.csv (file containing this information)
        :return: Boolean value which will be stored as part of the class variables
        """

        try:
            df = pd.read_csv(table_path)
            bool_ = len(df[df.Sample == self.sample]) == 2
            return bool_
        except IsADirectoryError as e:
            path = table_path + '/*.gz'
            files = glob.glob(path)
            bool_ = len(files) == 2
            return bool_

    def check_sample_length(self):
        """
        Check that the sample has the correct length:

        - All have to be multiple of 4
        - If the sample is paired, _R1 and _R2 need to be of the same length
        """

        if self.paired:
            path_R1 = self.path + self.sample + '_R1.gastq.gz'
            path_R2 = self.path + self.sample + '_R2.gastq.gz'
            res_R1 = subprocess.run(['zcat', path_R1, '|', 'wc', '-l'])
            res_R2 = subprocess.run(['zcat', path_R2, '|', 'wc', '-l'])
            if res_R1 != res_R2:
                raise Exception('check_sample_length:', self.sample, ' is paired but does not have the same length')
        else:
            path_R1 = self.path + self.sample + '.gastq.gz'
            res = subprocess.run(['zcat', path_R1, '|', 'wc', '-l'])
            if res % 4 != 0:
                raise Exception('check_sample_length:', self.sample, ' is single but does not have the correct '
                                                                     'sample length')

class Parent(LogMain):
    """
    This class will be a place holder for multiple methods used for BaseRecalibrator, ApplyBQSR and HaploType
    """

    def __init__(self, path, sample):
        self.path = path
        self.sample = sample
        self.chr_count = defaultdict(int)
        self.chr_time = defaultdict(float)
        self.chr_reads = defaultdict(int)

    @staticmethod
    def entropy(dict_):
        """
        Entropy of the dictionary output
        :param dict_: Dict over which we will perform the evaluation
        :return: Value of entropy
        """
        vals = list(dict_.values())
        res = 0
        for val in vals:
            res += val * np.log(val)

        return 5 if 1 / abs(res) == np.inf else 1 / abs(res)

    @staticmethod
    def distance(dict_, correct_dist_dict_):
        """
        KL Divergence between a dictionary and a dict which we consider to be "correct"
        :param dict_: Dict over which we will do the evaluation
        :param correct_dist_dict_: Supposidly correct dictionary
        :return: Value with the distance
        """
        for keys, values in correct_dist_dict_.items():
            if keys in dict_.keys():
                continue
            else:
                dict_[keys] = 0

        p = list(dict_.values())
        q = list(correct_dist_dict_.values())

        return abs(np.nansum([p[i] * np.log2(p[i] / q[i]) for i in range(len(p))]))

    @staticmethod
    def std(dict_):
        """
        Estimate variance of the dictionary
        :param dict_: Dict over which estimate will be performed
        :return: Std Value
        """
        std = np.std(np.array(list(dict_.values())))
        return std

    @staticmethod
    def normalize(dict_, target=1.0):
        """
        Method to normalize a dictionary based on its maximum value

        :param dict_: Dictionary to normzalie
        :return: Normalized dict
        """
        raw = max(dict_.values())
        factor = target / raw
        return {key: value * factor for key, value in dict_.items()}

    def single_paired(self, table_path='data/fastq.csv'):
        """
        Check whether the sample is paired (R1 and R2) or single (only R1). This method is relevant
        for cases in which two different log files are generated
        :param table_path: Path in which we can find the fastq.csv (file containing this information)
        :return: Boolean value which will be stored as part of the class variables
        """

        try:
            df = pd.read_csv(table_path)
            bool_ = len(df[df.Sample == self.sample]) == 2
            return bool_
        except FileNotFoundError as e:
            path = table_path + '/*.gz'
            files = glob.glob(path)
            bool_ = len(files) == 2
            return bool_

    def read_log(self, end_part, haplo_prefix=None):
        """
        Method to store the log file as part of the class variables
        """
        if haplo_prefix:
            with open(self.path + haplo_prefix + self.sample + end_part) as f:
                self.log_file = f.readlines()
        else:
            with open(self.path + self.sample + end_part) as f:
                self.log_file = f.readlines()

    def _read_template(self, path='template_recaldat.log'):
        """
        We read the global flag template so that we can compare this part more easily
        """
        with open(path) as f:
            self.log_template = f.readlines()

    def check_output_exists(self, file='data/OUTPUT/something.sam'):
        """
        Method to check that the output has been generated correctly
        """
        file_exist = not os.path.exists(file)
        if file_exist:
            raise Exception('check_output_exists: ' + self.sample + ' did not generate the output file')

    def check_running(self):
        """
        Check that the log file contains a row written running:

        - ``Running:``
        """
        if self.log_file[1][:-1] != 'Running:':
            raise Exception('check_running: ' + self.sample + ' did not start running...')

    def check_correct_sample(self):
        """
        Make sure that the file that is being processed is the actual sample

        - If we are processing sample HSRR062650 we should only have this value in the command line string
        """
        if len(re.findall(self.sample, self.log_file[2][:-1])) == 0:
            raise Exception('check_correct_sample: ' + self.sample + ' should be processed however another sample has '
                                                                     'been processed instead')

    def check_global_flags_start(self):
        """
        Check we have a section with global flags:

        - ``[Global flags]``
        """
        if self.log_file[3][:-1] != '[Global flags]':
            raise Exception('check_global_flags: ' + self.sample + ' defining the global flags did not seem to work')

    def check_final_section_success(self):
        """
        Makes sure the Recalibration has been success
        """
        if self.final_section[1][:-1] != 'SUCCESS':
            raise Exception('check_final_section_success: ' + self.sample + ' has not been successful during '
                                                                            'BaseRecalibration')

    def check_final_section_others(self):
        """
        Check that the other rows also have the expected output, in particular we check that the following strings are
        present:

        - ``PSYoungGen``
        - ``ParOldGen``
        - ``Metaspace``
        """
        s = ''.join(self.final_section[3:])
        if len(re.findall(r"\bPSYoungGen\b|\bParOldGen\b|\bMetaspace\b", s)) != 3:
            raise Exception('check_final_section_others: ' + self.sample + ' not all the expected fields are present '
                                                                           'in the last part of the log')

    def check_global_flags_length(self):
        """
        Method to check the length of the global flags section

        - The length should be 722
        """
        if len(self.global_flags) != 722:
            raise Exception('check_global_flags_length: ' + self.sample + ' does not have the expected number of rows')

    def check_global_flags_variables(self):
        """
        Method to check that all the global flags are set as expected

        - We compare the global flags with that of the template (check template_recaldat.log file for more details)
        - There are only 5 rows which are not supposed to be equal (they change based on the computer session used)
            - ``uintx CICompilerCount``
            - ``uintx InitialHeapSize``
            - ``uintx MaxHeapSize``
            - ``uintx MaxNewSize``
            - ``uintx NewSize``
            - ``uintx OldSize``
        - All other variables should be the same
        """
        for row in range(len(self.global_flags)):
            original = self.global_flags[row]
            template = self.log_template[row]

            if (original != template) & (row not in [55, 257, 304, 316, 349, 360]):
                raise Exception('check_global_flags_variables: ' + self.sample + ' does not have the right global '
                                                                                 'flags')

    # Removed CHRY
    def check_progressmeter_chromosomes(self):
        """
        Check all chromosomes are present in the progressmeter section
        - We also take the chance to check all chromsome ids are positive integers
        - We also check that the regions are also positive integers
        """
        # Split the data for further analysis
        chromosome = []
        chromosome_id = []
        regions = []

        chr_temp = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12',
                    'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22',
                    'chrX']

        for row in self.progressmeter[2:-1]:
            _, _, _, _, chrom, _, region, _ = row.split()
            chrom, chrom_id = chrom.split(':')
            chromosome.append(chrom)
            chromosome_id.append(int(chrom_id))
            regions.append(int(region))

        t1 = list(dict.fromkeys(chromosome)) != chr_temp
        t2 = any(i < 0 for i in chromosome_id)
        t3 = any(i < 0 for i in regions)

        if t1 | t2 | t3:
            error = ','.join(filter(None, [t1 * 't1', t2 * 't2', t3 * 't3']))
            self.progressmeter_analysis(title='ApplyBQSR')
            if t1:
                error = error + ' --> ' + list(set(chr_temp).difference(set(list(dict.fromkeys(chromosome)))))[0]
            raise Exception('check_progressmeter_chromosomes: ' + self.sample + ' has some strange chromosome values '
                                                                                'or is missing some chromosomes to be '
                                                                                'inspected. Issue in '
                                                                                'condition/s: ', error)

    def check_progressmeter_len(self):
        """
        Check correct length of this section:
        - For paired it should be 19 rows
        - For single it should be 18 rows
        """
        if self.paired:
            if len(self.progressmeter) != 19:
                raise Exception('check_progressmeter_len: ' + self.sample + ' does not have the right number of rows')
        else:
            if len(self.progressmeter) != 18:
                raise Exception('check_progressmeter_len: ' + self.sample + ' does not have the right number of rows')

    def check_progressmeter_start_end(self):
        """
        Check correct start and end statements
        - ``Starting traversal``
        - ``Traversal complete. Processed 180757 total reads in 2.6 minutes.``
        - We check that the strings are correct and also that the numeric part is larger than 0
        """
        start_text = self.progressmeter[0][-19:-1]
        end_text = re.findall('Traversal complete. Processed', self.progressmeter[-1][15:])
        end_nums = list(map(float, re.findall(r"[-+]?\d*\.\d+|\d+", self.progressmeter[-1][15:])))

        t1 = start_text != 'Starting traversal'
        t2 = len(end_text) != 1
        t3 = any(i <= 0 for i in end_nums)

        if t1 | t2 | t3:
            error = ','.join(filter(None, [t1*'t1', t2*'t2', t3*'t3']))
            raise Exception('check_progressmeter_start_end: ' + self.sample + ' does not have the correct start and '
                                                                              'end statements in the ProgressMeter '
                                                                              'section. Issue in condition/s: ', error)

    def progressmeter_analysis(self, title=False):
        """
        Visual test to see whether the output is in line with our expectations
        """

        for row in self.progressmeter[2:-1]:
            row_split = row.split()
            self.chr_count[re.findall(r'chr.*:', row_split[4])[0][3:-1]] += 1
            self.chr_time[re.findall(r'chr.*:', row_split[4])[0][3:-1]] += float(row_split[5])
            self.chr_reads[re.findall(r'chr.*:', row_split[4])[0][3:-1]] += int(row_split[6])


def df_func(list_, header=['Sample', 'Sample_Score']):
    """
    Aux function to filter and return output
    """

    df = pd.DataFrame(list_, columns=header)
    df = df.sort_values(by=[header[1]]).reset_index(drop=True)
    return df
