from .log_analysis_new import LogMain


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

        - ``check_lines()``
        - ``check_start_end()``
        - ``check_output_exists()``
        """
        self.check_lines()
        self.check_start_end()

    def check_lines(self):
        """
        Check correct number of lines in log

        - We expect 21 lines in the log
        """
        if self.log_file_1:
            if (len(self.log_file_1) != 21) | (len(self.log_file_1) != 22):
                raise Exception('check_lines: ' + self.sample + '_R1 does not have the correct number of lines')

        if self.log_file_2:
            if (len(self.log_file_2) != 21) | (len(self.log_file_2) != 22):
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