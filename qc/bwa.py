from .log_analysis_new import LogMain
import pandas as pd
import re
import glob
import os 


class Bwa(LogMain):
    """
    This class will check the bwa log
    """

    def __init__(self, path, sample, table_path='data/fastq.csv'):
        self.log_file = None
        self.path = path
        self.sample = sample
        self.dict_ = None
        self.paired = self.single_paired(table_path)
        self.read_log()
        self.process = []
        self.mem_pestat = []
        self.mem_process_seqs = []
        self.split_log()

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

    def read_log(self):
        """
        Method to store the log file as part of the class variables
        """
        with open(self.path + self.sample + '.log') as f:
            self.log_file = f.readlines()

    def split_log(self):
        """
        This method is used to split the log into smaller and more manageable files (due to the large size of the file
        its better to split before to simplify our life). The outputs are stored as part of the class variables

        - ``self.process``:
        - ``self.mem_pestat``:
        - ``self.mem_process_seqs``:
        """

        for row in self.log_file:
            if row[:12] == '[M::process]':
                self.process.append(row)
            elif row[:15] == '[M::mem_pestat]':
                self.mem_pestat.append(row)
            elif row[:21] == '[M::mem_process_seqs]':
                self.mem_process_seqs.append(row)

    def check_log(self, check_process=True, check_mem_process_seqs=True, check_mem_pestat=True, check_start_statement=True, 
                  check_correct_sample=True, check_tmp=True):
        """
        This method will run all the methods implemented for this class

        - ``check_process()``
        - ``check_start_statement()``
        - ``check_correct_sample()``
        - ``check_candidate_pairs()``
        - ``check_not_enough_pairs(threshold=100)``
        - ``check_enough_pairs(threshold=100)``
        - ``check_output_exists()``
        """

        if check_process:
            self.check_process()
        if check_mem_process_seqs:
            self.check_mem_process_seqs()
        if check_mem_pestat:
            self.check_mem_pestat()
        if check_start_statement:
            self.check_start_statement()
        if check_correct_sample:
            self.check_correct_sample()
        if check_tmp:
            self.check_tmp_files()

    def _batch(self, iterable, n=1):
        l = len(iterable)
        for ndx in range(0, l, n):
            yield iterable[ndx:min(ndx + n, l)]

    def check_process(self):
        """
        Check the log file defined as process (which starts with [M::process]. The checks done are:

        ``check_num_sequence``
        ``check_consistency``
        """
        for batch in self._batch(self.process, 2):
            self.check_num_sequence(batch[0])
            # self.check_consistency(batch)

    def check_mem_process_seqs(self):
        """
        Check the log file defined as mem_process_seqs (which starts with [M::mem_process_seqs]. The checks done are:

        """
        for batch in self.mem_process_seqs:
            self.check_positive_nums(batch)

    def check_mem_pestat(self):
        """
        Check the log file defined as mem_pestat (which starts with [M::mem_pestat]. The checks done are:

        ``check_not_enough_pairs``
        ``check_enough_pairs``
        """
        if self.paired:

            batch_nums = []
            for num, i in enumerate(self.mem_pestat):
                if i[:27] == '[M::mem_pestat] # candidate':
                    batch_nums.append(num)

            for i in range(len(batch_nums)-1):
                batch = self.mem_pestat[batch_nums[i]: batch_nums[i + 1]]
                self.check_not_enough_pairs(batch)
                self.check_enough_pairs(batch)

    def check_positive_nums(self, row):
        """
        Make sure we are processing a positive number of sequences and the run times are all positive

        - ``[M::process] read 400000 sequences (40000000 bp)...``
        - 400000 must be a positive number
        """
        nums = list(map(float, re.findall(r"[-+]?\d*\.\d+|\d+", row)))
        if any(i <= 0 for i in nums):
            raise Exception('check_positive_nums: ' + self.sample + ' did not process a positive number of sequences')

    def check_num_sequence(self, row):
        """
        Make sure we are reading a positive number of sequences

        - ``[M::process] read 400000 sequences (40000000 bp)...``
        - 400000 must be a positive number
        """
        seq = re.findall(r'\d+', row)
        if any(int(i) < 0 for i in seq):
            raise Exception('check_num_sequence: ' + self.sample + ' did not read a positive number of sequences')

    def check_consistency(self, rows):
        """
        Check consistency in terms of single-end and paired-end sequences. This check is only done for paired samples

        - ``[M::process] read 400000 sequences (40000000 bp)...``
        - ``[M::process] 0 single-end sequences; 400000 paired-end sequences``
        - 40000 (in the first line) has to be equal to 0 + 400000 (second line)
        """
        if self.paired:
            seq1 = int(re.findall(r'\d+', rows[0])[0])
            seq2 = list(map(int, re.findall(r'\d+', rows[1])))
            if seq1 != sum(seq2):
                raise Exception('check_consistency: ' + self.sample + ' has inconsistency in terms of paired-end and '
                                                                      'single-end sequences')
        else:
            seq1 = int(re.findall(r'\d+', rows[0])[0])
            seq2 = int(re.findall(r'\d+', rows[0])[0])
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

        - ``[main] CMD: bwa mem -p -t 8 -R @RG\tID:HSRR062625\tLB:HSRR062625\tSM:HSRR062625\tPU:unknown ...``
        - If we are processing sample HSRR062625 we should only have this value in this string
        """
        if len(re.findall(self.sample, self.log_file[-2])) == 0:
            raise Exception('check_correct_sample: ' + self.sample + ' should be processed however another sample has '
                                                                     'been processed instead')

    def check_not_enough_pairs(self, batch, threshold=10):
        """
        Whenever there are not enough pairs the command skips those pairs. This check is only done for paired samples

        - ``[M::mem_pestat] skip orientation FF as there are not enough pairs``
        - We make sure that the orientation FF actually does not have enough pairs

        :param threshold: Threshold used to define whether we actually have enough pairs
        """
        if self.paired:
            FF, FR, RF, RR = list(map(int, batch[0][64:-2].strip().split(',')))
            self.dict_ = {'FF':FF,
                     'FR': FR,
                     'RF': RF,
                     'RR': RR}
            for i in batch[1:]:
                if (i[:32] == '[M::mem_pestat] skip orientation') & (i[36:-1] == 'as there are not enough pairs'):
                    if self.dict_[i[33:35]] > threshold:
                        raise Exception('check_not_enough_pairs ' + self.sample + ' has skipped an orientation due to low '
                                                                              'number of pairs where this is not the '
                                                                              'case')

    def check_enough_pairs(self, batch, threshold=10):
        """
        When we have enough pairs, there are several things that need to be checked. This check is only done for paired
        samples

        - ``[M::mem_pestat] analyzing insert size distribution for orientation FR...``
        - This log line is only shown for orientations with sufficient pairs (FR in this case). Once there are enough pairs 4 steps are done on them. We check that these steps are really done

        :param threshold: Threshold used to define whether we actually have enough pairs
        """
        if self.paired:
            for num, i in enumerate(batch):
                if i[-4:-1] == '...':
                    if self.dict_[i[-6:-4]] >= threshold:
                        if ((batch[num + 1][:40] != '[M::mem_pestat] (25, 50, 75) percentile:') |
                                (batch[num + 2][
                                 :71] != '[M::mem_pestat] low and high boundaries for computing mean and std.dev:') |
                                (batch[num + 3][:33] != '[M::mem_pestat] mean and std.dev:') |
                                (batch[num + 4][
                                 :57] != '[M::mem_pestat] low and high boundaries for proper pairs:')):
                            raise Exception('check_enough_pairs ' + self.sample + ' not all the steps of BWA have been '
                                                                                  'executed')
                    else:
                        raise Exception('check_enough_pairs ' + self.sample + ' not enough pairs, however code has '
                                                                              'considered to have enough pairs anyway')

    def check_tmp_files(self):
        """
        Ensure that no tmp files are still laying around which would imply that the pipeline has not been completed correctly
        """

        files = [f for f in os.listdir(self.path) if os.path.isfile(f)]
        tmp_files = [i for i in files if 'tmp' in i]
        if len(tmp_files) > 0:
            raise Exception('check_tmp_files ' + self.sample + ' still has some tmp files laying around: \n', tmp_files)
