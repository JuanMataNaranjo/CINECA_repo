from abc import ABCMeta, abstractmethod
import re
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from collections import defaultdict
import subprocess


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

        df = pd.read_csv(table_path)

        bool_ = len(df[df.Sample == self.sample]) == 2

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
        #TODO: Put this outside and read instead
        self.true_apply_chr_count = {'1': 1.0, '2': 0.8, '3': 0.8, '4': 0.8, '5': 0.8, '6': 0.8, '7': 0.6000000000000001,
                                '8': 0.4, '9': 0.6000000000000001, '10': 0.4, '11': 0.6000000000000001,
                                '12': 0.6000000000000001, '13': 0.2, '14': 0.4, '15': 0.4, '16': 0.2, '17': 0.4,
                                '18': 0.4, '19': 0.2, '20': 0.2, '21': 0.2, '22': 0.2, 'X': 0.2, 'Y': 0.4}
        self.true_apply_chr_time = {'1': 0.125, '2': 0.25, '3': 0.385, '4': 0.515, '5': 0.65, '6': 0.785, '7': 0.675,
                               '8': 0.49000000000000005, '9': 0.8, '10': 0.5750000000000001, '11': 0.925, '12': 1.0,
                               '13': 0.35000000000000003, '14': 0.7250000000000001, '15': 0.76, '16': 0.39, '17': 0.81,
                               '18': 0.8400000000000001, '19': 0.435, '20': 0.44000000000000006, '21': 0.45,
                               '22': 0.45999999999999996, 'X': 0.4650000000000001, 'Y': 0.9550000000000001}
        self.true_apply_chr_reads = {'1': 0.11632613625962351, '2': 0.24521240351073462, '3': 0.3809838599774223,
                                '4': 0.5150772422630611, '5': 0.6435770439452033, '6': 0.7685579749204187,
                                '7': 0.6631138954712338, '8': 0.4850041189093535, '9': 0.7918882911102749,
                                '10': 0.5711045796169923, '11': 0.9216390208182901, '12': 1.0,
                                '13': 0.35080903515819667, '14': 0.727846878273516, '15': 0.7625169078686425,
                                '16': 0.3941745400551222, '17': 0.8142829538173646, '18': 0.8488207715073175,
                                '19': 0.4373874927537706, '20': 0.4459812665900516, '21': 0.4544835091073663,
                                '22': 0.46197890711605155, 'X': 0.46945396483163326, 'Y': 0.9574905875293663}

        self.true_base_chr_count = {'1': 0.96, '2': 1.0, '3': 0.8, '4': 0.8, '5': 0.72, '6': 0.72, '7': 0.68, '8': 0.64,
                               '9': 0.56, '10': 0.6, '11': 0.6, '12': 0.56, '13': 0.44, '14': 0.4, '15': 0.36,
                               '16': 0.36, '17': 0.36, '18': 0.32, '19': 0.28, '20': 0.28, '21': 0.16, '22': 0.16,
                               'X': 0.4, 'Y': 0.08}
        self.true_base_chr_time = {'1': 0.10115722266560254, '2': 0.30925778132482046, '3': 0.3978451715881883,
                              '4': 0.5313248204309656, '5': 0.5919792498004787, '6': 0.6997206703910613,
                              '7': 0.7609736632083002, '8': 0.8038707102952913, '9': 0.773343974461293,
                              '10': 0.9014365522745412, '11': 0.9766560255387073, '12': 0.979050279329609,
                              '13': 0.8150438946528332, '14': 0.7757382282521949, '15': 0.7272545889864326,
                              '16': 0.7541899441340781, '17': 0.7811252992817238, '18': 0.716879489225858,
                              '19': 0.644852354349561, '20': 0.6614126097366321, '21': 0.3852753391859537,
                              '22': 0.39066241021548287, 'X': 1.0, 'Y': 0.20371109337589782}
        self.true_base_chr_reads = {'1': 0.09953048645616466, '2': 0.31896821035301304, '3': 0.41314227468936626,
                               '4': 0.554865131721466, '5': 0.6201158847901782, '6': 0.7334718234866225,
                               '7': 0.7933752756311296, '8': 0.8296525508563661, '9': 0.7921097097407528,
                               '10': 0.9187119047564754, '11': 0.9917019922971503, '12': 0.9917225149872645,
                               '13': 0.8250714303630919, '14': 0.7847261018974366, '15': 0.7342790492977819,
                               '16': 0.7603474263406447, '17': 0.786744166425335, '18': 0.7222299499018331,
                               '19': 0.648519287908259, '20': 0.6643673652742401, '21': 0.38703285226627504,
                               '22': 0.39215440359922377, 'X': 0.9999999999999999, 'Y': 0.20331589892803148}

        self.true_haplo_chr_count = {'1': 0.9509803921568627, '2': 1.0, '3': 0.8529411764705882, '4': 0.8333333333333334,
                                '5': 0.7549019607843137, '6': 0.7549019607843137, '7': 0.6568627450980392,
                                '8': 0.5882352941176471, '9': 0.5098039215686274, '10': 0.5980392156862745,
                                '11': 0.5686274509803921, '12': 0.5588235294117647, '13': 0.4411764705882353,
                                '14': 0.37254901960784315, '15': 0.3431372549019608, '16': 0.39215686274509803,
                                '17': 0.3529411764705882, '18': 0.3333333333333333, '19': 0.22549019607843138,
                                '20': 0.3333333333333333, '21': 0.22549019607843138, '22': 0.18627450980392157,
                                'X': 0.3333333333333333, 'Y': 0.08823529411764705}
        self.true_haplo_chr_time = {'1': 0.09779541337539592, '2': 0.31184611024165443, '3': 0.4358648905100495,
                               '4': 0.5770496247643231, '5': 0.6511602114628654, '6': 0.773256602053014,
                               '7': 0.7719750089341832, '8': 0.7695720218363755, '9': 0.7268974355814614,
                               '10': 0.9235973333004719, '11': 0.949044350515718, '12': 1.0, '13': 0.8371514128331843,
                               '14': 0.7393313534362722, '15': 0.7071929413794384, '16': 0.8402568115441971,
                               '17': 0.7854317366819069, '18': 0.7662448089317189, '19': 0.5318118522717474,
                               '20': 0.8061960098091168, '21': 0.5589964140038698, '22': 0.4699873072989195,
                               'X': 0.8595423233804487, 'Y': 0.23152472612108604}
        self.true_haplo_chr_reads = {'1': 0.10026755156243426, '2': 0.3170263732140641, '3': 0.43986302238764197,
                                '4': 0.5779748238479174, '5': 0.6493011479754309, '6': 0.7680823428086322,
                                '7': 0.7687299296240473, '8': 0.7679727665566674, '9': 0.7274631021751309,
                                '10': 0.9229204425520959, '11': 0.948005031443449, '12': 1.0, '13': 0.8364298290370424,
                                '14': 0.7386816863074641, '15': 0.707876279923501, '16': 0.8379892646383053,
                                '17': 0.7806048958615608, '18': 0.7624005472129715, '19': 0.5296087956892019,
                                '20': 0.8008151274835964, '21': 0.5515464567877524, '22': 0.4621923427323935,
                                'X': 0.853342230274628, 'Y': 0.23219735892811305}

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

        df = pd.read_csv(table_path)

        bool_ = len(df[df.Sample == self.sample]) == 2

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

    def progressmeter_analysis(self, title='BaseRecalibrator'):
        """
        Visual test to see whether the output is in line with our expectations
        """

        for row in self.progressmeter[2:-1]:
            row_split = row.split()
            self.chr_count[re.findall(r'chr.*:', row_split[4])[0][3:-1]] += 1
            self.chr_time[re.findall(r'chr.*:', row_split[4])[0][3:-1]] += float(row_split[5])
            self.chr_reads[re.findall(r'chr.*:', row_split[4])[0][3:-1]] += int(row_split[6])

        if title:

            fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, sharex=True, figsize=(15, 15))

            ax1.bar(range(len(self.chr_count)), list(self.chr_count.values()), align='center')
            ax1.set_xticks(range(len(self.chr_count)))
            ax1.set_xticklabels(list(self.chr_count.keys()))
            ax1.set_ylabel('Number of Chromosomes \n Processed')

            ax2.bar(range(len(self.chr_time)), list(self.chr_time.values()), align='center')
            ax2.set_xticks(range(len(self.chr_count)))
            ax2.set_xticklabels(list(self.chr_count.keys()))
            ax2.set_ylabel('Time Required for \n Processing (min)')

            ax3.bar(range(len(self.chr_reads)), list(self.chr_reads.values()), align='center')
            ax3.set_xticks(range(len(self.chr_count)))
            ax3.set_xticklabels(list(self.chr_count.keys()))
            ax3.set_ylabel('Number of Reads \n Performed')
            ax3.set_xlabel('Chromosomes affected')

            fig.suptitle('Overview of ' + title + ' Process \n ' + self.sample, fontsize=20)

            plt.show()

    def compute_score(self, true_dict):
        """
        In order to estimate which samples are more likely to be corrupt we will first compute a score which will
        later be used to perform a filter
        :return: Score value
        """

        score = 0

        score += self.entropy(self.normalize(self.chr_count))
        score += self.entropy(self.normalize(self.chr_time))
        score += self.entropy(self.normalize(self.chr_reads))
        score += self.std(self.normalize(self.chr_count))
        score += self.std(self.normalize(self.chr_time))
        score += self.std(self.normalize(self.chr_reads))
        score += self.distance(self.normalize(self.chr_count), true_dict[0])
        score += self.distance(self.normalize(self.chr_time), true_dict[1])
        score += self.distance(self.normalize(self.chr_reads), true_dict[2])

        return score

def df_func(list_, header=['Sample', 'Sample_Score']):
    """
    Aux function to filter and return output
    """

    df = pd.DataFrame(list_, columns=header)
    df = df.sort_values(by=[header[1]]).reset_index(drop=True)
    return df