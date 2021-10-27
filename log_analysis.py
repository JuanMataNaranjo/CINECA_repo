from abc import ABCMeta, abstractmethod
import re
import numpy as np


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


class Seqtk(LogMain):
    """
    This class will check the seqtk log
    """

    def __init__(self, path, sample):
        self.log_file = None
        self.path = path
        self.sample = sample

    def read_log(self):
        with open(self.path + self.sample + '.log') as f:
            self.log_file = f.readlines()

    def check_log(self):

        self._check_lines()
        self._check_num_sequence()
        self._check_consistency()
        print('Everything fine :)')

    def _check_lines(self):
        """
        Check the number of lines of the log
        """
        if len(self.log_file) != 16:
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
        seq1 = int(re.findall(r'\d+', self.log_file[1])[0])
        seq2 = list(map(int, re.findall(r'\d+', self.log_file[2])))
        if seq1 != sum(seq2):
            raise Exception('check_consistency: ' + self.sample + ' has inconsistency in terms of paired-end and '
                                                                  'single-end sequences')


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
