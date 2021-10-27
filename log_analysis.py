from abc import ABCMeta, abstractmethod
import re
import numpy as np
import pandas as pd


class LogMain(metaclass=ABCMeta):
    """
    Abstract object
    All classes that will inherit from this class will need to include the following methods
    """

    @abstractmethod
    def read_log(self):
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


class Bwa(LogMain):
    """
    This class will check the bwa log
    """
    # TODO: Can we somehow include the numbers that the log provides us for anything? (e.g. percentiles, mean, std, ...)
    # TODO: What is the threshold for "enough pairs"?

    def __init__(self, path, sample):
        self.log_file = None
        self.path = path
        self.sample = sample
        self.paired = self.single_paired()
        self.dict_ = None

    def single_paired(self, table_path='data/fastq.csv'):
        """
        Method which will check whether the sample is paired (R1 and R2) or single (only R1). This method is relevant
        for cases in which this difference produces different log files
        :param table_path: Path in which we can find this information
        :return:
        """

        df = pd.read_csv(table_path)

        bool_ = len(df[df.Sample == self.sample]) == 2

        return bool_

    def read_log(self):
        with open(self.path + self.sample + '.log') as f:
            self.log_file = f.readlines()

    def check_log(self):

        self._check_lines()
        self._check_num_sequence()
        self._check_consistency()
        self._check_start_statement()
        self._check_correct_sample()
        self._check_candidate_pairs()
        self._check_not_enough_pairs(threshold=100)
        self._check_enough_pairs(threshold=100)
        print('Everything fine :)')

    def _check_lines(self):
        """
        Check the number of lines of the log
        """
        if self.paired:
            if len(self.log_file) != 16:
                raise Exception('check_lines: ' + self.sample + ' has the wrong number of log lines')
        else:
            if len(self.log_file) != 6:
                raise Exception('check_lines: ' + self.sample + ' has the wrong number of log lines')

    def _check_num_sequence(self):
        """
        Make sure we are reading a positive number of sequences
        """
        seq = re.findall(r'\d+', self.log_file[1])
        if any(int(i) < 0 for i in seq):
            raise Exception('check_num_sequence: ' + self.sample + ' did not read a positive number of sequences')

    def _check_consistency(self):
        """
        Check consistency in terms of single-end and paired-end sequences
        """
        if self.paired:
            seq1 = int(re.findall(r'\d+', self.log_file[1])[0])
            seq2 = list(map(int, re.findall(r'\d+', self.log_file[2])))
            if seq1 != sum(seq2):
                raise Exception('check_consistency: ' + self.sample + ' has inconsistency in terms of paired-end and '
                                                                      'single-end sequences')
        else:
            seq1 = int(re.findall(r'\d+', self.log_file[1])[0])
            seq2 = int(re.findall(r'\d+', self.log_file[2])[0])
            if seq1 != seq2:
                raise Exception('check_consistency: ' + self.sample + ' has inconsistency in terms of read sequences '
                                                                      'and processed sequences')

    def _check_start_statement(self):
        """
        Make sure that the initial log statement is the expected one
        """
        if self.log_file[0][:-1] != '[M::bwa_idx_load_from_disk] read 3171 ALT contigs':
            raise Exception('check_start_statement: ' + self.sample + ' does not have the correct log starting '
                                                                      'statement')

    def _check_finish_statement(self):
        """
        The last log statement needs to provide us with positive runtime values
        """
        text = self.log_file[-1][:17]
        nums = list(map(int, re.findall(r"[-+]?\d*\.\d+|\d+", self.log_file[-1])))
        if any(nums <= 0 for i in nums) | text != '[main] Real time:':
            raise Exception('check_finish_statement: ' + self.sample + ' does not have the final statement we expected')

    def _check_correct_sample(self):
        """
        Make sure that the file that is being processed is the actual sample
        """

        if len(re.findall(self.sample, self.log_file[-2])) == 0:
            raise Exception('check_correct_sample: ' + self.sample + ' should be processed however another sample has '
                                                                     'been processed instead')

    def _check_candidate_pairs(self):
        """
        Check the candidate pairs are being evaluated
        """
        if self.paired:
            if '[M::mem_pestat] # candidate unique pairs for (FF, FR, RF, RR):' != self.log_file[3][:62]:
                raise Exception(
                    'check_candidate_pairs: ' + self.sample + ' does not have the expected candidate pairs output')

    def _check_not_enough_pairs(self, threshold=100):
        """
        Whenever there are not enough pairs the command skips those pairs
        """
        if self.paired:
            FF, FR, RF, RR = list(map(int, self.log_file[3][64:-2].strip().split(',')))
            self.dict_ = {'FF':FF,
                     'FR': FR,
                     'RF': RF,
                     'RR': RR}
            for i in self.log_file[4:-4]:
                if i[:32] == '[M::mem_pestat] skip orientation':
                    if self.dict_[i[33:35]] > threshold:
                        raise Exception('check_not_enough_pairs ' + self.sample + ' has skipped an orientation due to low '
                                                                              'number of pairs where this is not the '
                                                                              'case')

    def _check_enough_pairs(self, threshold=100):
        """
        When we have enough pairs, there are several things that need to be checked
        :param threshold: Threshold which will define when we have enough pairs or not
        """
        if self.paired:
            for num, i in enumerate(self.log_file[4:-4]):
                if (i[-4:-1] == '...'):
                    if self.dict_[self.log_file[4:-4][1][-6:-4]] > threshold:
                        if ((self.log_file[4 + num + 1][:40] != '[M::mem_pestat] (25, 50, 75) percentile:') |
                                (self.log_file[4 + num + 2][
                                 :71] != '[M::mem_pestat] low and high boundaries for computing mean and std.dev:') |
                                (self.log_file[4 + num + 3][:33] != '[M::mem_pestat] mean and std.dev:') |
                                (self.log_file[4 + num + 4][
                                 :57] != '[M::mem_pestat] low and high boundaries for proper pairs:')):
                            raise Exception('check_enough_pairs ' + self.sample + ' not all the steps of BWA have been '
                                                                                  'executed')
                    else:
                        raise Exception('check_enough_pairs ' + self.sample + ' not enough pairs, however code has '
                                                                              'considered to have enough pairs anyway')






class Fastqc(LogMain):
    """
    Probaably this is not very usefull since the log of fastq is very basic...
    """

    def __init__(self, path, sample):
        self.log_file_1 = None
        self.log_file_2 = None
        self.path = path
        self.sample = sample

    def read_log(self):
        try:
            with open(self.path + self.sample + '_R1_fastqc.log') as f:
                self.log_file_1 = f.readlines()
        except:
            self.log_file_1 = False

        try:
            with open(self.path + self.sample + '_R2_fastqc.log') as f:
                self.log_file_2 = f.readlines()
        except:
            self.log_file_2 = False

    def check_log(self):

        self._check_lines()
        self._check_start_end()
        print('Everything fine :)')

    def _check_lines(self):
        """
        Check correct number of lines in log
        """
        if self.log_file_1:
            if len(self.log_file_1) != 22:
                raise Exception('check_lines: ' + self.sample + '_R1 does not have the correct number of lines')

        if self.log_file_2:
            if len(self.log_file_2) != 22:
                raise Exception('check_lines: ' + self.sample + '_R2 does not have the correct number of lines')

    def _check_start_end(self):
        """
        Check start and end statements are the expected ones
        """
        if self.log_file_1:
            if (self.log_file_1[0][:16] != 'Started analysis') | (self.log_file_1[-1][:17] != 'Analysis complete'):
                raise Exception('check_lines: ' + self.sample + '_R1 does not seem to have been processed properly')

        if self.log_file_2:
            if (self.log_file_2[0][:16] != 'Started analysis') | (self.log_file_2[-1][:17] != 'Analysis complete'):
                raise Exception('check_lines: ' + self.sample + '_R2 does not seem to have been processed properly')
