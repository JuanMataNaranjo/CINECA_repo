import pandas as pd
import os
import sys
import argparse
import glob

from qc.log_analysis_new import *
from qc.applybqsr import ApplyBQSR
from qc.baserecalibrator import BaseRecalibrator
from qc.bwa import Bwa
from qc.fastqc import Fastqc
from qc.haplotype import HaploType
from qc.samsort import SamSort

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Quality Checks")

    parser.add_argument("--checks", type=str, required=True)
    parser.add_argument("--path", type=str, default='/data/output')
    args = parser.parse_args()

    checks = args.checks.split(',')
    path = args.path

    print(checks)
    print(path)

    centers = ['1a77728b-011a-4407-b60e-ec90b6430b99', 'ed615b8d-fa86-48bb-b0f8-c841e1aeb0eb', 'd7c41726-13a8-4abd-b185-68198fec12f4', 
               'fb203f69-94b9-42aa-9f34-c9ee8219a22e', 
               'eaabd6f1-1b34-4196-882c-198465045d71']

    for center in centers:
        try:
            for study in os.listdir(path + '/' + center + '/'):
                for sample in os.listdir(path + '/' + center + '/' + study + '/'):
                    path_parent = path + '/' + center + '/' + study + '/' + sample + '/'
                    try:
                        samples = [os.path.basename(x) for x in glob.glob(path_parent + '/*.gz')][0]
                        samples = samples[:-12]
                    except IndexError as e:
                        print('Cannot find .gz files for this sample')
                        continue
                    if 'bwa' in checks:
                        bwa_class = Bwa(path=path_parent + 'bwa/', sample=samples, table_path=path_parent + 'bwa/')
                        try:
                            bwa_class.check_log()
                        except Exception as e:
                            print(center + '/' + study + '/' + sample)
                            print('Issue with Bwa: ', samples)
                            print(e)
                            print('===================')
                    
                    if 'fastqc' in checks:
                        fastqc_class = Fastqc(path=path_parent +'fastqc/', sample=samples)
                        try:
                            fastqc_class.check_log()
                        except Exception as e:
                            print(center + '/' + study + '/' + sample)
                            print('Issue with Fastqc: ', samples)
                            print(e)
                            print('===================')

                    if 'samsort' in checks:
                        try:
                            samsort_class = SamSort(path=path_parent+'bwa/', sample=samples, table_path='data/fastq.csv')
                        except FileNotFoundError as e:
                            print('SAMSORTBLASTER file not found')
                        try:
                            samsort_class.check_log()
                        except Exception as e:
                            print(center + '/' + study + '/' + sample)
                            print('Issue with SamSort: ', samples)
                            print(e)
                            print('===================')

                    if 'baserecalibrator' in checks:
                        try:
                            baserecalibrator_class = BaseRecalibrator(path=path_parent + 'gatk_bsr/', sample=samples, 
                                                                table_path='data/fastq.csv')
                        except FileNotFoundError as e:
                            print('GATK_BSR files not found')
                        try:
                            baserecalibrator_class.check_log(title=None, progressmeter_analysis=False, check_running=True, check_correct_sample=True, check_global_flags_start=True,
                  check_final_section=True, check_global_flags=True, check_baserecalibrator=True, check_featuremanager=True,
                  check_progressmeter=True)
                                    
                        except Exception as e:
                            print(center + '/' + study + '/' + sample)
                            print('Issue with BaseRecalibrator: ', samples)
                            print(e)
                            print('=================')

                    if 'applybqsr' in checks:
                        try:
                            applybqsr_class = ApplyBQSR(path=path_parent+'gatk_bsr/', sample=samples)
                        except FileNotFoundError as e:
                            print('APPLYBQSR File not Found')
                        try:
                            applybqsr_class.check_log(title=None, check_running=True, check_correct_sample=False, check_global_flags_start=False,
                  check_final_section=False, check_global_flags=False, check_applybqsr=False, check_featuremanager=False, check_progressmeter=False,
                  progressmeter_analysis=True)
                                    
                        except Exception as e:
                            print(center + '/' + study + '/' + sample)
                            print('Issue with ApplyBQSR: ', samples)
                            print(e)
                            print('================')

                    if 'haplotype' in checks:
                        try:
                            haplo_class = HaploType(path=path_parent+'gatk_gvcf/', sample=samples)
                        except FileNotFoundError as e:
                            print('HAPLOTYPE File not Found')
                        try:
                            haplo_class.check_log(warning_plot=False, title=None)
                                    
                        except Exception as e:
                            print(center + '/' + study + '/' + sample)
                            print('Issue with HaploType: ', sample)
                            print(e)
                            print('================')


                    # break
        except PermissionError as e:
            print('No permissions for ', center)
            continue

