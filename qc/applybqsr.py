from .log_analysis_new import Parent
import re

class ApplyBQSR(Parent):
    """
    This class will check the Apply BQSR log
    """

    def __init__(self, path, sample, input_parameter_path='template_recaldat.log'):
        super().__init__(path, sample)
        self.path = path
        self.sample = sample
        self.log_file = None
        self.log_template = None
        self.read_log(end_part='_sort_nodup.bqsr.log')
        self._read_template(path=input_parameter_path)
        self.global_flags = []
        self.applybqsr = []  # Contains info on BaseRecalibrationEngine and BaseRecalibrator
        self.featuremanager = []
        self.progressmeter = []
        self.final_section = None
        self.split_log()


    def split_log(self):
        """
        This method is used to split the log into smaller and more manageable files (due to the large size of the file
        its better to split before to simplify our life). The outputs are stored as part of the class variables

        - ``self.global_flags``: Everything after [Global Flags] until that part ends
        - ``self.applybqsr``: Information on ApplyBQSR and BaseRecalibrator
        - ``self.featuremanager``: Information on FeatureManager
        - ``self.progressmeter``: Information on the ProgressMeter
        - ``self.final_section``: Last part of the log with final results
        """

        global_flags_bool = False
        applybqsr_bool = False
        featuremanager_bool = False
        progressmeter_bool = False

        for row in self.log_file:
            if row[:-1] == '[Global flags]':
                global_flags_bool = True
                continue

            if bool(re.search(r'(\d+-\d+-\d+)', row[:10])):
                global_flags_bool = False

            if bool(re.search(r'INFO  ApplyBQSR', row)):
                applybqsr_bool = True

            if bool(re.search(r'INFO  ProgressMeter', row)):
                progressmeter_bool = True

            if bool(re.search(r'INFO  FeatureManager', row)):
                featuremanager_bool = True

            if global_flags_bool:
                self.global_flags.append(row)
            elif applybqsr_bool:
                self.applybqsr.append(row)
                applybqsr_bool = False
            elif featuremanager_bool:
                self.featuremanager.append(row)
                featuremanager_bool = False
            elif progressmeter_bool:
                self.progressmeter.append(row)
                progressmeter_bool = False

        self.final_section = self.log_file[-9:]

    def check_log(self, title='ApplyBQSR', check_running=True, check_correct_sample=True, check_global_flags_start=True,
                  check_final_section=True, check_global_flags=True, check_applybqsr=True, check_featuremanager=True, check_progressmeter=True,
                  progressmeter_analysis=True):
        """
        This method will run all the methods implemented for this class

        - ``check_running()``
        - ``check_correct_sample()``
        - ``check_global_flags_start()``
        - ``check_final_section()``
        - ``check_global_flags()``
        - ``check_applybqsr()``
        - ``check_featuremanager()``
        - ``check_progressmeter()``
        - ``check_output_exists()``
        """

        if check_running:
            self.check_running()
        if check_correct_sample:
            self.check_correct_sample()
        if check_global_flags_start:
            self.check_global_flags_start()
        if check_final_section:
            self.check_final_section()
        if check_global_flags:
            self.check_global_flags()
        if check_applybqsr:
            self.check_applybqsr()
        if check_featuremanager:
            self.check_featuremanager()
        if check_progressmeter:
            self.check_progressmeter()
        #if progressmeter_analysis:
        #    self.progressmeter_analysis(title=title)

    def check_final_section(self):
        """
        Method to run checks on the final section part of the log, including the following methods:

        - ``check_final_section_others``

        """
        self.check_final_section_others()

    def check_global_flags(self):
        """
        Method to run checks on the global flags part of the log, including the following methods:

        - ``check_global_flags_length``
        - ``check_global_flags_variables``
        """
        self.check_global_flags_length()
        self.check_global_flags_variables()

    def check_applybqsr(self):
        """
        Method to run checks on the ApplyBQSR part of the log, including the following methods:

        - ``check_applybqsr_engine``
        - ``check_applybqsr_quantization``
        """
        self.check_applybqsr_engine()
        self.check_applybqsr_quantization()

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

        #self.check_progressmeter_chromosomes()
        self.check_progressmeter_start_end()

    def check_final_section_others(self):
        """
        Check that the other rows also have the expected output, in particular we check that the following strings are
        present:

        - ``PSYoungGen``
        - ``ParOldGen``
        - ``Metaspace``
        """
        s = ''.join(self.final_section[1:])
        if len(re.findall(r"\bPSYoungGen\b|\bParOldGen\b|\bMetaspace\b", s)) != 3:
            raise Exception('check_final_section_others: ' + self.sample + ' not all the expected fields are present '
                                                                           'in the last part of the log')

    def check_featuremanager_files(self):
        """
        Check that the files used to extract the features are the correct ones. The files which should be used are:

        - ``hg38_resources/wgs_calling_regions.hg38.interval_list``
        """

        list_ = ['hg38_resources/wgs_calling_regions.hg38.interval_list']

        s_check = ''.join(['\\' + 'b' + i + '\\' + 'b|' for i in list_])[:-1]
        s = ''.join(self.featuremanager)
        if len(re.findall(s_check, s)) != len(list_):
            raise Exception('check_featuremanager_files: ' + self.sample + ' did not get the features from the correct '
                                                                           'input files')

    def check_applybqsr_engine(self):
        """
        Method to check that the applybqsr engine starts and ends correctly. We make sure that the following
        are present:

        - ``Initializing engine``
        - ``Done initializing engine``
        - ``Shutting down engine``
        """
        t1 = self.applybqsr[19][-20:-1] != 'Initializing engine'
        t2 = self.applybqsr[20][-25:-1] != 'Done initializing engine'
        t3 = self.applybqsr[-1][-21:-1] != 'Shutting down engine'

        if t1 | t2 | t3:
            error = ','.join(filter(None, [t1 * 't1', t2 * 't2', t3 * 't3']))
            raise Exception('check_applybqsr_engine: ' + self.sample + ' applybqsr engine did not work '
                                                                              'properly. Issue in condition/s: ', error)

    def check_applybqsr_quantization(self):
        """
        Method to check that the applybqsr engine generates the report well and that the quantization passes
        """
        t1 = self.applybqsr[-2][-22:-2] != 'WellformedReadFilter'

        if t1:
            raise Exception('check_applybqsr_quantization: ' + self.sample + ' the quantization part did not '
                                                                                    'work properly')
