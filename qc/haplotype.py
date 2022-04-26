from .log_analysis_new import Parent
import re
import matplotlib.pyplot as plt
from collections import defaultdict

class HaploType(Parent):
    """
    This class will check the Haplotype log
    """

    def __init__(self, path, sample, table_path='data/fastq.csv'):
        super().__init__(path, sample)
        self.log_file = None
        self.path = path
        self.sample = sample
        self.paired = self.single_paired(table_path)
        self.read_log(end_part='_sort_nodup.g.vcf.log')
        self.haplotype = []
        self.featuremanager = []
        self.progressmeter = []
        self.warning = []
        self.split_log()

    def split_log(self):
        """
        This method is used to split the log into smaller and more manageable files (due to the large size of the file
        its better to split before to simplify our life). The outputs are stored as part of the class variables

        - ``self.haplotype``: Information on Haplotype
        - ``self.featuremanager``: Information on FeatureManager
        - ``self.progressmeter``: Information on the ProgressMeter
        """

        haplotype_bool = False
        featuremanager_bool = False
        progressmeter_bool = False
        warn_bool = False

        for row in self.log_file:

            if bool(re.search(r'INFO  HaplotypeCaller', row)):
                haplotype_bool = True

            if bool(re.search(r'read\(s\) filtered by:', row)):
                haplotype_bool = True

            if bool(re.search(r'INFO  ProgressMeter', row)):
                progressmeter_bool = True

            if bool(re.search(r'INFO  FeatureManager', row)):
                featuremanager_bool = True

            if bool(re.search(r' WARN ', row)):
                warn_bool = True

            if haplotype_bool:
                self.haplotype.append(row)
                haplotype_bool = False
            elif featuremanager_bool:
                self.featuremanager.append(row)
                featuremanager_bool = False
            elif progressmeter_bool:
                self.progressmeter.append(row)
                progressmeter_bool = False
            elif warn_bool:
                self.warning.append(row)
                warn_bool = False

    def check_log(self, warning_plot=True, title='HaploTypeCaller', score=True):
        """
        This method will run all the methods implemented for this class
        - ``check_running()``
        - ``check_correct_sample()``
        - ``check_haplotype()``
        - ``check_featuremanager()``
        - ``check_warnings()``
        - ``check_progressmeter()``
        - ``check_output_exists()``
        """

        self.check_len()
        self.check_running()
        self.check_correct_sample()
        self.check_haplotype()
        self.check_featuremanager()
        if warning_plot:
            self.check_warnings()
        self.check_progressmeter()
        self.progressmeter_analysis(title=title)
        if score:
            return self.compute_score([self.true_base_chr_count, self.true_base_chr_time, self.true_base_chr_reads])

    def check_len(self):
        """
        Checks if the log file has the number of expected lines
        """
        rows = len(self.haplotype)
        if rows < 32:
            raise Exception('check_len: ' + self.sample + ' does not have the expected log file length')


    def check_featuremanager(self):
        """
        Method to run checks on the feature manager part of the log, including the following methods:

        - ``check_featuremanager_files``
        """

        self.check_featuremanager_files()

    def check_haplotype(self):
        """
        Method to check everything regarding the INFO Haplotype part. Methods included are:

        - ``check_haplotype_engine``
        - ``check_haplotype_filters``
        """

        self.check_haplotype_engine()
        self.check_haplotype_filters()

    def check_warnings(self):
        # TODO: How shall we interpret these warnings?
        """
        There are multiple warnings in the log file. We will make sure that these warnings are serious or not

        """
        self.warning_stats()

    def check_progressmeter(self):
        """
        Method to run checks on the progress meter part of the log, including the following methods:

        - ``check_progressmeter_chromosomes``
        - ``check_progressmeter_start_end``
        """

        self.check_progressmeter_chromosomes()
        self.check_progressmeter_start_end()

    def warning_stats(self):
        """
        Method to plot some insights on the warning outputs
        """

        DepthPerSampleHC = []
        summary_depth = defaultdict(int)
        StrandBiasBySample = []
        summary_strand = defaultdict(int)
        InbreedingCoeff = []
        summary_inbreed = defaultdict(int)
        for row in self.warning:
            if bool(re.search(r'DepthPerSampleHC', row)):
                DepthPerSampleHC.append(row)
                summary_depth[re.findall(r'chr.*:', row)[0][3:-1]] += 1
            elif bool(re.search(r'StrandBiasBySample', row)):
                StrandBiasBySample.append(row)
                summary_strand[re.findall(r'chr.*:', row)[0][3:-1]] += 1
            elif bool(re.search(r'InbreedingCoeff', row)):
                InbreedingCoeff.append(row)
                summary_inbreed[re.findall(r'chr.*:', row)[0][3:-1]] += 1

        plt.figure(figsize=(12, 8))
        plt.bar(range(len(summary_depth)), list(summary_depth.values()), align='center')
        plt.xticks(range(len(summary_depth)), list(summary_depth.keys()))
        plt.ylabel('Number of Errors')
        plt.xlabel('Chromosomes affected')
        plt.title(self.sample + '\n' + ' WARN  DepthPerSampleHC')
        plt.show()

        plt.figure(figsize=(12, 8))
        plt.bar(range(len(summary_strand)), list(summary_strand.values()), align='center')
        plt.xticks(range(len(summary_strand)), list(summary_strand.keys()))
        plt.ylabel('Number of Errors')
        plt.xlabel('Chromosomes affected')
        plt.title(self.sample + '\n' + ' WARN  StrandBiasBySample')
        plt.show()

        plt.figure(figsize=(12, 8))
        plt.bar(range(len(summary_inbreed)), list(summary_inbreed.values()), align='center')
        plt.xticks(range(len(summary_inbreed)), list(summary_inbreed.keys()))
        plt.ylabel('Number of Errors')
        plt.xlabel('Chromosomes affected')
        plt.title(self.sample + '\n' + ' WARN  InbreedingCoeff')
        plt.show()

    def check_haplotype_engine(self):
        """
        Method to check that the base haplotype engine starts and ends correctly. We make sure that the following
        are present:

        - ``Initializing engine``
        - ``Done initializing engine``
        - ``Shutting down engine``
        """
        t1 = self.haplotype[19][-20:-1] != 'Initializing engine'
        t2 = self.haplotype[20][-25:-1] != 'Done initializing engine'
        t3 = self.haplotype[-1][-21:-1] != 'Shutting down engine'

        if t1 | t2 | t3:
            error = ','.join(filter(None, [t1 * 't1', t2 * 't2', t3 * 't3']))
            raise Exception('check_haplotype_engine: ' + self.sample + ' haplotype engine did not work properly. '
                                                                       'Issue in condition/s: ', error)

    def check_haplotype_filters(self):
        """
        Method to check that the haplotype engine has execited the correct filters:

        - ``MappingQualityAvailableReadFilter``
        - ``MappedReadFilter``
        - ``NotSecondaryAlignmentReadFilter``
        - ``NotDuplicateReadFilter``
        - ``PassesVendorQualityCheckReadFilter``
        - ``NonZeroReferenceLengthAlignmentReadFilter``
        - ``GoodCigarReadFilter``
        - ``WellformedReadFilter``
        """
        t1 = self.haplotype[25][-35:-2] != 'MappingQualityAvailableReadFilter'
        t2 = self.haplotype[26][-18:-2] != 'MappedReadFilter'
        t3 = self.haplotype[27][-33:-2] != 'NotSecondaryAlignmentReadFilter'
        t4 = self.haplotype[28][-24:-2] != 'NotDuplicateReadFilter'
        t5 = self.haplotype[29][-36:-2] != 'PassesVendorQualityCheckReadFilter'
        t6 = self.haplotype[30][-43:-2] != 'NonZeroReferenceLengthAlignmentReadFilter'
        t7 = self.haplotype[31][-21:-2] != 'GoodCigarReadFilter'
        t8 = self.haplotype[32][-22:-2] != 'WellformedReadFilter'

        if t1 | t2 | t3 | t4 | t5 | t6 | t7 | t8:
            error = ','.join(filter(None, [t1 * 't1', t2 * 't2', t3 * 't3', t4 * 't4', t5 * 't5',
                                           t6 * 't6', t7 * 't7', t8 * 't8']))
            raise Exception('check_haplotype_filters: ' + self.sample + ' not all the correct filters have been '
                                                                               'used on the haplotype. '
                                                                        'Issue in condition/s: ', error)

    def check_featuremanager_files(self):
        """
        Check that the files used to extract the features are the correct ones. The files which should be used are:

        - ``hg38_resources/wgs_calling_regions.hg38.interval_list``
        """

        list_ = ['hg38_resources/wgs_calling_regions.hg38.interval_list']

        s_check = ''.join(['\\' + 'b' + i + '\\' + 'b|' for i in list_])[:-1]
        s = ''.join(self.featuremanager)
        if len(re.findall(s_check, s)) != 1:
            raise Exception('check_featuremanager_files: ' + self.sample + ' did not get the features from the correct '
                                                                           'input files')