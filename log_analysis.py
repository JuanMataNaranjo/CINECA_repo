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
    # TODO: Check if the output expected actually exists

    def __init__(self, path, sample):
        self.log_file = None
        self.path = path
        self.sample = sample
        self.dict_ = None
        self.paired = self.single_paired()
        self.read_log()

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

    def read_log(self):
        """
        Method to store the log file as part of the class variables
        """
        with open(self.path + self.sample + '.log') as f:
            self.log_file = f.readlines()

    def check_log(self):
        """
        This method will run all the methods implemented for this class
        """
        self.check_lines()
        self.check_num_sequence()
        self.check_consistency()
        self.check_start_statement()
        self.check_correct_sample()
        self.check_candidate_pairs()
        self.check_not_enough_pairs(threshold=100)
        self.check_enough_pairs(threshold=100)

    def check_lines(self):
        """
        Check the number of lines of the log:

        - In case of paired samples (R1 + R2) the log should contain 16 lines
        - In case of single samples (R1) the log should contain 6 lines
        """
        if self.paired:
            if len(self.log_file) != 16:
                raise Exception('check_lines: ' + self.sample + ' has the wrong number of log lines')
        else:
            if len(self.log_file) != 6:
                raise Exception('check_lines: ' + self.sample + ' has the wrong number of log lines')

    def check_num_sequence(self):
        """
        Make sure we are reading a positive number of sequences

        - ``[M::process] read 400000 sequences (40000000 bp)...``
        - 400000 must be a positive number
        """
        seq = re.findall(r'\d+', self.log_file[1])
        if any(int(i) < 0 for i in seq):
            raise Exception('check_num_sequence: ' + self.sample + ' did not read a positive number of sequences')

    def check_consistency(self):
        """
        Check consistency in terms of single-end and paired-end sequences. This check is only done for paired samples

        - ``[M::process] read 400000 sequences (40000000 bp)...``
        - ``[M::process] 0 single-end sequences; 400000 paired-end sequences``
        - 40000 (in the first line) has to be equal to 0 + 400000 (second line)
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

    def check_start_statement(self):
        """
        Make sure that the initial log statement is the expected one

        - ``[M::bwa_idx_load_from_disk] read 3171 ALT contigs``
        """
        if self.log_file[0][:-1] != '[M::bwa_idx_load_from_disk] read 3171 ALT contigs':
            raise Exception('check_start_statement: ' + self.sample + ' does not have the correct log starting '
                                                                      'statement')

    def check_finish_statement(self):
        """
        Make sure that the final log statement is the expected one. We check that the text is the expected one and also
        that the runtime is a positive value

        - ``[main] Real time: 23.460 sec; CPU: 133.799 sec``
        - Check for [main] Real time:
        - Check real time and CPU times are both positive
        """
        text = self.log_file[-1][:17]
        nums = list(map(int, re.findall(r"[-+]?\d*\.\d+|\d+", self.log_file[-1])))
        if any(i <= 0 for i in nums) | text != '[main] Real time:':
            raise Exception('check_finish_statement: ' + self.sample + ' does not have the final statement we expected')

    def check_correct_sample(self):
        """
        Make sure that the file that is being processed is the actual sample

        - ``[main] CMD: bwa mem -p -t 8 -R @RG\tID:HSRR062625\tLB:HSRR062625\tSM:HSRR062625\tPU:unknown\ ...``
        - If we are processing sample HSRR062625 we should only have this value in this string
        """
        if len(re.findall(self.sample, self.log_file[-2])) == 0:
            raise Exception('check_correct_sample: ' + self.sample + ' should be processed however another sample has '
                                                                     'been processed instead')

    def check_candidate_pairs(self):
        """
        Check the candidate pairs are being evaluated. This check is only done for paired samples

        - ``[M::mem_pestat] # candidate unique pairs for (FF, FR, RF, RR):``
        - We make sure that this string is present in the log file
        """
        if self.paired:
            if '[M::mem_pestat] # candidate unique pairs for (FF, FR, RF, RR):' != self.log_file[3][:62]:
                raise Exception(
                    'check_candidate_pairs: ' + self.sample + ' does not have the expected candidate pairs output')

    def check_not_enough_pairs(self, threshold=100):
        """
        Whenever there are not enough pairs the command skips those pairs. This check is only done for paired samples

        - ``[M::mem_pestat] skip orientation FF as there are not enough pairs``
        - We make sure that the orientation FF actually does not have enough pairs

        :param threshold: Threshold used to define whether we actually have enough pairs
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

    def check_enough_pairs(self, threshold=100):
        """
        When we have enough pairs, there are several things that need to be checked. This check is only done for paired
        samples

        - ``[M::mem_pestat] analyzing insert size distribution for orientation FR...``
        - This log line is only shown for orientations with sufficient pairs (FR in this case). Once there are enough pairs 4 steps are done on them. We check that these steps are really done

        :param threshold: Threshold used to define whether we actually have enough pairs
        """
        if self.paired:
            for num, i in enumerate(self.log_file[4:-4]):
                if i[-4:-1] == '...':
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
    This class will check the fastqc log
    """

    def __init__(self, path, sample):
        self.log_file_1 = None
        self.log_file_2 = None
        self.path = path
        self.sample = sample
        self.read_log()

    def read_log(self):
        """
        Method to store the log file as part of the class variables
        """
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
        """
        This method will run all the methods implemented for this class
        """
        self.check_lines()
        self.check_start_end()

    def check_lines(self):
        """
        Check correct number of lines in log

        - We expect 22 lines in the log
        """
        if self.log_file_1:
            if len(self.log_file_1) != 22:
                raise Exception('check_lines: ' + self.sample + '_R1 does not have the correct number of lines')

        if self.log_file_2:
            if len(self.log_file_2) != 22:
                raise Exception('check_lines: ' + self.sample + '_R2 does not have the correct number of lines')

    def check_start_end(self):
        """
        Check start and end statements are the expected ones

        - ``Started analysis`` is the start log line
        - ``Analysis complete`` is the end log line
        """
        if self.log_file_1:
            if (self.log_file_1[0][:16] != 'Started analysis') | (self.log_file_1[-1][:17] != 'Analysis complete'):
                raise Exception('check_lines: ' + self.sample + '_R1 does not seem to have been processed properly')

        if self.log_file_2:
            if (self.log_file_2[0][:16] != 'Started analysis') | (self.log_file_2[-1][:17] != 'Analysis complete'):
                raise Exception('check_lines: ' + self.sample + '_R2 does not seem to have been processed properly')


class SamSort(LogMain):
    """
    This class will check the samsort log
    """
    # TODO: Make sure check_unmated makes sense (I think the mate them later on so we should expect them to be 0 anyway
    #  but not sure)

    def __init__(self, path, sample):
        self.log_file = None
        self.path = path
        self.sample = sample
        self.paired = self.single_paired()
        self.read_log()
        self.dups = None

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

    def read_log(self):
        """
        Method to store the log file as part of the class variables
        """
        with open(self.path + self.sample + '_sort_nodup.sam.log') as f:
            self.log_file = f.readlines()

    def check_log(self):
        """
        This method will run all the methods implemented for this class
        """
        self.check_lines()
        self.check_start_statement()
        self.check_finish_statement()
        self.check_correct_sample()
        self.check_third_line()
        self.check_unmated()
        self.check_header()
        self.check_rows()
        self.check_table_sums()
        self.check_removals()

    def check_lines(self):
        """
        Check the number of lines of the log:

        - In case of paired samples (R1 + R2) the log should contain 15 lines
        - In case of single samples (R1) the log should contain 14 lines
        """
        if self.paired:
            if len(self.log_file) != 15:
                raise Exception('check_lines: ' + self.sample + ' has the wrong number of log lines')
        else:
            if len(self.log_file) != 14:
                raise Exception('check_lines: ' + self.sample + ' has the wrong number of log lines')

    def check_start_statement(self):
        """
        Make sure that the initial log statement is the expected one

        - ``samblaster: Version 0.1.26``
        """
        if self.log_file[0][:-1] != 'samblaster: Version 0.1.26':
            raise Exception('check_start_statement: ' + self.sample + ' does not have the correct log starting '
                                                                      'statement')

    def check_finish_statement(self):
        """
        Make sure that the final log statement is the expected one.

        - ``[bam_sort_core] merging from files and in-memory blocks...``
        - We remove the numbers from the string we check since they are variable
        """
        if re.sub(" \d+", "", self.log_file[-1][:-1]) != '[bam_sort_core] merging from files and in-memory blocks...':
            raise Exception('check_finish_statement: ' + self.sample + ' does not have the final statement we expected')

    def check_correct_sample(self):
        """
        Make sure that the file that is being processed is the actual sample

        - ``samblaster: Opening OUTPUT_NEW/HSRR062650/bwa/HSRR062650.sam for read``
        - If we are processing sample HSRR062650 we should only have this value in this string
        """
        if len(re.findall(self.sample, self.log_file[1][:-1])) == 0:
            raise Exception('check_correct_sample: ' + self.sample + ' should be processed however another sample has '
                                                                     'been processed instead')

    def check_third_line(self):
        """
        The third line of these logs are all constant:

        - ``samblaster: Outputting to stdout``
        """
        if self.log_file[2][:-1] != 'samblaster: Outputting to stdout':
            raise Exception('check_third_line: ' + self.sample + ' does not have the expected output in line 3')

    def check_unmated(self):
        """
        Make sure that we have 0 mated pairs

        - ``samblaster: Found 0 of 200000 (0.000%) total read ids are marked paired yet are unmated.``
        - We check 0 and 0.000
        """
        if ((int(re.findall(r'\d+', self.log_file[4])[0]) != 0) |
            (float(re.findall(r"[-+]?\d*\.\d+|\d+", self.log_file[4])[-1]) != float(0))):
            raise Exception('check_unmated: ' + self.sample + ' has mated pairs when it should not have them yet')

    def check_header(self):
        """
        Check that the header of the table is correct
        """
        if self.log_file[6].replace(" ", "")[:-1] != 'samblaster:PairTypeType_ID_Count%Type/All_IDsDup_ID_Count%Dups/Type_ID_Count%Dups/All_Dups%Dups/All_IDs':
            raise Exception('check_header: ' + self.sample + ' does not have the expected table header')


    @staticmethod
    def _check_digit(x):
        """
        Aux function, not relevant
        """
        dig = False
        try:
            float(x)
            dig = True
        except:
            pass
        return dig

    def check_rows(self):
        """
        Check that the row names are all present

        - In the case of paired it should be ['Both Unmapped', 'Orphan/Singleton', 'Both', 'Total']
        - In the case of singles it should be ['Unmapped Orphan/Singleton', 'Mapped Orphan/Singleton', 'Total']
        """
        # Extract Data
        list_ = []
        for i in self.log_file[8:-3]:
            x = [j for j in i[12:-1].split() if not self._check_digit(j)]
            list_.append(' '.join(x))

        # Check values
        if self.paired:
            if list_ != ['Both Unmapped', 'Orphan/Singleton', 'Both Mapped', 'Total']:
                raise Exception('check_rows: ' + self.sample + ' does not have the expected row names')
        else:
            if list_ != ['Unmapped Orphan/Singleton', 'Mapped Orphan/Singleton', 'Total']:
                raise Exception('check_rows: ' + self.sample + ' does not have the expected row names')

    def check_table_sums(self):
        """
        Check that the table provided by the log file is correct:

        - We are checking that the values and the sums add up
        """
        # Extract data from table
        list_ = []
        for i in self.log_file[8:-3]:
            x = [float(j) for j in i[12:-1].split() if self._check_digit(j)]
            list_.append(x)

        self.dups = int(list_[-1][2])

        # Perform the summation accordingly
        sum_ = lambda x: np.isclose(sum(x[:-1]), x[-1], atol=0.001)
        if any(list(~sum_(np.array(list_)))):
            raise Exception('check_table_sums: ' + self.sample + ' has a mismatch in the total sums')

    def check_removals(self):
        """
        Check that the number of removed duplicates matches the number in the table and that the text of the second last
        line of the log are as expected
        """
        nums = list(map(float, re.findall(r"[-+]?\d*\.\d+|\d+", self.log_file[-2])))
        text = re.sub(r"\d+", "", self.log_file[-2])[:-1]
        # noinspection PyTypeChecker
        if ((text != 'samblaster: Removed         of      (.%) total read ids as duplicates using k memory in .S CPU seconds and S wall time.') |
            (int(nums[0]) != self.dups) |
            (any(int(i) < 0 for i in nums))):
            raise Exception('check_removals: ' + self.sample + ' has some issue (text or numeric related)')

