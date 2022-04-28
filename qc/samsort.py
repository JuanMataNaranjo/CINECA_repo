from .log_analysis_new import LogMain
import pandas as pd
import re
import numpy as np

class SamSort(LogMain):
    """
    This class will check the samsort log
    """

    def __init__(self, path, sample, table_path='data/fastq.csv'):
        self.log_file = None
        self.path = path
        self.sample = sample
        self.paired = self.single_paired(table_path)
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
        with open(self.path + self.sample + '_samblaster.log') as f:
            self.log_file = f.readlines()
        with open(self.path + self.sample + '_sort_nodup.sam.log') as f:
            self.log_file_2 = f.readlines()

    def check_log(self):
        """
        This method will run all the methods implemented for this class

        - ``check_lines()``
        - ``check_start_statement()``
        - ``check_finish_statement()``
        - ``check_correct_sample()``
        - ``check_third_line()``
        - ``check_unmated()``
        - ``check_header()``
        - ``check_rows()``
        - ``check_table_sums()``
        - ``check_removals()``
        - ``check_output_exists()``
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

        - In case of paired samples (R1 + R2) the log should contain 14 lines
        - In case of single samples (R1) the log should contain 13 lines
        """
        if self.paired:
            if len(self.log_file) != 14:
                raise Exception('check_lines: ' + self.sample + ' which is paired has the wrong number of log lines')
        else:
            if len(self.log_file) != 13:
                raise Exception('check_lines: ' + self.sample + ' which is single has the wrong number of log lines')

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
        if re.sub(" \d+", "", self.log_file_2[-1][:-1]) != '[bam_sort_core] merging from files and in-memory blocks...':
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
        for i in self.log_file[8:-2]:
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
        for i in self.log_file[8:-2]:
            x = [float(j) for j in i[12:-1].split() if self._check_digit(j)]
            list_.append(x)

        self.dups = int(list_[-1][2])

        # Perform the summation accordingly
        sum_ = lambda x: np.isclose(sum(x[:-1]), x[-1], atol=0.001)
        if any(list(~sum_(np.array(list_))[:3])+list(~sum_(np.array(list_))[4:])):
            raise Exception('check_table_sums: ' + self.sample + ' has a mismatch in the total sums')

    def check_removals(self):
        """
        Check that the number of removed duplicates matches the number in the table and that the text of the second last
        line of the log are as expected
        """
        nums = list(map(float, re.findall(r"[-+]?\d*\.\d+|\d+", self.log_file[-1])))
        text = re.sub(r"\d+", "", self.log_file[-2])[:-1]

        t1 = int(nums[0]) != self.dups
        t2 = any(int(i) < 0 for i in nums)

        # noinspection PyTypeChecker
        if t1 | t2:
            error = ','.join(filter(None, [t1 * 't1', t2 * 't2']))
            raise Exception('check_removals: ' + self.sample + ' has some issue (text or numeric related). Issue in '
                                                               'condition/s: ' + error)