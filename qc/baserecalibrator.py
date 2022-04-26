from .log_analysis_new import Parent
import re

class BaseRecalibrator(Parent):
    """
    This class will check the baserecalibrator log
    """

    def __init__(self, path, sample, table_path='data/fastq.csv', 
                 input_parameter_path='template_recaldat.log'):
        super().__init__(path, sample)
        self.log_file = None
        self.log_template = None
        self.path = path
        self.sample = sample
        self.paired = self.single_paired(table_path)
        self.read_log(end_part='_sort_nodup.recaldat.log')
        self._read_template(path=input_parameter_path)
        self.global_flags = []
        self.baserecalibrator = []  # Contains info on BaseRecalibrationEngine and BaseRecalibrator
        self.featuremanager = []
        self.progressmeter = []
        self.final_section = None
        self.split_log()

    def split_log(self):
        """
        This method is used to split the log into smaller and more manageable files (due to the large size of the file
        its better to split before to simplify our life). The outputs are stored as part of the class variables

        - ``self.global_flags``: Everything after [Global Flags] until that part ends
        - ``self.baserecalibrator``: Information on BaseRecalibrationEngine and BaseRecalibrator
        - ``self.featuremanager``: Information on FeatureManager
        - ``self.progressmeter``: Information on the ProgressMeter
        - ``self.final_section``: Last part of the log with final results
        """

        global_flags_bool = False
        baserecalibrator_bool = False
        featuremanager_bool = False
        progressmeter_bool = False

        for row in self.log_file:
            if row[:-1] == '[Global flags]':
                global_flags_bool = True
                continue

            if bool(re.search(r'(\d+-\d+-\d+)', row[:10])):
                global_flags_bool = False

            if bool(re.search(r'INFO  BaseRecalibrat', row)):
                baserecalibrator_bool = True

            if bool(re.search(r'read\(s\) filtered by:', row)):
                baserecalibrator_bool = True

            if bool(re.search(r'INFO  ProgressMeter', row)):
                progressmeter_bool = True

            if bool(re.search(r'INFO  FeatureManager', row)):
                featuremanager_bool = True

            if global_flags_bool:
                self.global_flags.append(row)
            elif baserecalibrator_bool:
                self.baserecalibrator.append(row)
                baserecalibrator_bool = False
            elif featuremanager_bool:
                self.featuremanager.append(row)
                featuremanager_bool = False
            elif progressmeter_bool:
                self.progressmeter.append(row)
                progressmeter_bool = False

        self.final_section = self.log_file[-11:]

    def check_log(self, title='BaseRecalibrator', score=True):
        """
        This method will run all the methods implemented for this class

        - ``check_running()``
        - ``check_correct_sample()``
        - ``check_global_flags_start()``
        - ``check_final_section()``
        - ``check_global_flags()``
        - ``check_baserecalibrator()``
        - ``check_featuremanager()``
        - ``check_progressmeter()``
        - ``check_output_exists()``
        """

        self.check_running()
        self.check_correct_sample()
        self.check_global_flags_start()
        self.check_final_section()
        self.check_global_flags()
        self.check_baserecalibrator()
        self.check_featuremanager()
        self.check_progressmeter()
        self.progressmeter_analysis(title=title)
        if score:
            return self.compute_score([self.true_base_chr_count, self.true_base_chr_time, self.true_base_chr_reads])

    def check_final_section(self):
        """
        Method to run checks on the final section part of the log, including the following methods:

        - ``check_final_section_length``
        - ``check_final_section_success``
        - ``check_final_section_others``

        """

        self.check_final_section_success()
        self.check_final_section_others()

    def check_global_flags(self):
        """
        Method to run checks on the global flags part of the log, including the following methods:

        - ``check_global_flags_length``
        - ``check_global_flags_variables``
        """

        self.check_global_flags_length()
        self.check_global_flags_variables()

    def check_baserecalibrator(self):
        """
        Method to run checks on the base recalibrator part of the log, including the following methods:

        - ``check_baserecalibrator_engine``
        - ``check_baserecalibrator_covariates``
        - ``check_baserecalibrator_filters``
        - ``check_baserecalibrator_quantization``
        """

        self.check_baserecalibrator_len()
        self.check_baserecalibrator_engine()
        self.check_baserecalibrator_covariates()
        self.check_baserecalibrator_filters()
        self.check_baserecalibrator_quantization()

    def check_featuremanager(self):
        """
        Method to run checks on the feature manager part of the log, including the following methods:

        - ``check_featuremanager_files``
        """

        self.check_featuremanager_files()

    def check_progressmeter(self):
        """
        Method to run checks on the progress meter part of the log, including the following methods:

        - ``check_progressmeter_chromosomes``
        - ``check_progressmeter_start_end``
        """

        self.check_progressmeter_chromosomes()
        self.check_progressmeter_start_end()

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

    def check_featuremanager_files(self):
        """
        Check that the files used to extract the features are the correct ones. The files which should be used are:

        - ``hg38_resources/dbsnp_reannotated.vcf``
        - ``hg38_resources/Mills_and_1000G_gold_standard.indels.hg38.vcf``
        - ``hg38_resources/1000G_omni2.5.hg38.vcf``
        - ``hg38_resources/wgs_calling_regions.hg38.interval_list``
        """

        list_ = ['hg38_resources/dbsnp_reannotated.vcf',
                 'hg38_resources/Mills_and_1000G_gold_standard.indels.hg38.vcf',
                 'hg38_resources/1000G_omni2.5.hg38.vcf',
                 'hg38_resources/wgs_calling_regions.hg38.interval_list']

        s_check = ''.join(['\\' + 'b' + i + '\\' + 'b|' for i in list_])[:-1]
        s = ''.join(self.featuremanager)
        if len(re.findall(s_check, s)) != 4:
            raise Exception('check_featuremanager_files: ' + self.sample + ' did not get the features from the correct '
                                                                           'input files')

    
    def check_baserecalibrator_len(self):
        """
        Method to check that the base recalibrator engine starts and ends correctly. We make sure that the following
        are present:

        - ``Initializing engine``
        - ``Done initializing engine``
        - ``Shutting down engine``
        """

        l = len(self.baserecalibrator)

        if l < 21:
            raise Exception('check_baserecalibrator_len: ' + self.sample + ' baserecalibrator engine is not as long as expected')
    
    
    def check_baserecalibrator_engine(self):
        """
        Method to check that the base recalibrator engine starts and ends correctly. We make sure that the following
        are present:

        - ``Initializing engine``
        - ``Done initializing engine``
        - ``Shutting down engine``
        """
        t1 = self.baserecalibrator[19][-20:-1] != 'Initializing engine'
        t2 = self.baserecalibrator[20][-25:-1] != 'Done initializing engine'
        t3 = self.baserecalibrator[-1][-21:-1] != 'Shutting down engine'

        if t1 | t2 | t3:
            error = ','.join(filter(None, [t1 * 't1', t2 * 't2', t3 * 't3']))
            raise Exception('check_baserecalibrator_engine: ' + self.sample + ' baserecalibrator engine did not work '
                                                                              'properly. Issue in condition/s: ', error)

    def check_baserecalibrator_covariates(self):
        """
        Method to check that the base recalibrator engine considers the correct covariate values:

        - ``ReadGroupCovariate``
        - ``QualityScoreCovariate``
        - ``ContextCovariate``
        - ``CycleCovariate``
        """
        t1 = self.baserecalibrator[22][-19:-1] != 'ReadGroupCovariate'
        t2 = self.baserecalibrator[23][-22:-1] != 'QualityScoreCovariate'
        t3 = self.baserecalibrator[24][-17:-1] != 'ContextCovariate'
        t4 = self.baserecalibrator[25][-15:-1] != 'CycleCovariate'

        if t1 | t2 | t3 | t4:
            error = ','.join(filter(None, [t1 * 't1', t2 * 't2', t3 * 't3', t4 * 't4']))
            raise Exception('check_baserecalibrator_engine: ' + self.sample + ' not all the covariates have been used. '
                                                                              'Issue in condition/s: ', error)

    def check_baserecalibrator_filters(self):
        """
        Method to check that the base recalibrator engine has execited the correct filters:

        - ``MappingQualityAvailableReadFilter``
        - ``MappedReadFilter``
        - ``NotSecondaryAlignmentReadFilter``
        - ``NotDuplicateReadFilter``
        - ``PassesVendorQualityCheckReadFilter``
        - ``WellformedReadFilter``
        """
        t1 = self.baserecalibrator[27][-35:-2] != 'MappingQualityAvailableReadFilter'
        t2 = self.baserecalibrator[28][-18:-2] != 'MappedReadFilter'
        t3 = self.baserecalibrator[29][-33:-2] != 'NotSecondaryAlignmentReadFilter'
        t4 = self.baserecalibrator[30][-24:-2] != 'NotDuplicateReadFilter'
        t5 = self.baserecalibrator[31][-36:-2] != 'PassesVendorQualityCheckReadFilter'
        t6 = self.baserecalibrator[32][-22:-2] != 'WellformedReadFilter'

        if t1 | t2 | t3 | t4 | t5 | t6:
            error = ','.join(filter(None, [t1 * 't1', t2 * 't2', t3 * 't3', t4 * 't4', t5 * 't5', t6 * 't6']))
            raise Exception('check_baserecalibrator_filters: ' + self.sample + ' not all the correct filters have been '
                                                                               'used on the base recalibrator. '
                                                                               'Issue in condition/s: ', error)

    def check_baserecalibrator_quantization(self):
        """
        Method to check that the base recalibrator engine generates the report well and that the quantization passes

        - ``Calculating quantized quality scores...``
        - ``Writing recalibration report...``
        - ``...done!``
        - ``BaseRecalibrator was able to recalibrate 371510 reads``
        - We also check that the recalibration is done over a number greater than 0
        """
        t1 = self.baserecalibrator[33][-40:-1] != 'Calculating quantized quality scores...'
        t2 = self.baserecalibrator[34][-32:-1] != 'Writing recalibration report...'
        t3 = self.baserecalibrator[35][-9:-1] != '...done!'
        t4 = re.sub(r'[0-9]', '', self.baserecalibrator[36][38:-1]) != 'BaseRecalibrator was able to recalibrate  reads'
        t5 = int(re.findall(r'\d+', self.baserecalibrator[36][38:-1])[0]) < 0

        if t1 | t2 | t3 | t4 | t5:
            error = ','.join(filter(None, [t1 * 't1', t2 * 't2', t3 * 't3', t4 * 't4', t5 * 't5']))
            raise Exception('check_baserecalibrator_quantization: ' + self.sample + ' the quantization part did not '
                                                                                    'work properly. Issue in '
                                                                                    'condition/s: ', error)