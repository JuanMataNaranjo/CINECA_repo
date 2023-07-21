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
               'eaabd6f1-1b34-4196-882c-198465045d71', 'gatk_db']

    for center in centers:
        print('=============================')
        print('CENTER:')
        print(center)
        try:
            for study in os.listdir(path + '/' + center + '/'):
                for sample in os.listdir(path + '/' + center + '/' + study + '/'):
                    print(center + '||' + study + '||' + sample)
                    path_parent = path + '/' + center + '/' + study + '/' + sample + '/'
                    try:
                        samples = [os.path.basename(x) for x in glob.glob(path_parent + '/*.gz')][0]
                        samples = samples[:-12]
                    except IndexError as e:
                        print('Cannot find .gz files for this sample')
                        continue
                    if 'bwa' in checks:
                        print('BWA')
                        bwa_class = Bwa(path=path_parent + 'bwa/', sample=samples, table_path=path_parent + 'bwa/')
                        try:
                            bwa_class.check_log()
                        except Exception as e:
                            print('Issue with Bwa: ', sample)
                            print(e)
                    
                    if 'fastqc' in checks:
                        print('FASTQC')
                        fastqc_class = Fastqc(path=path_parent +'fastqc/', sample=samples)
                        try:
                            fastqc_class.check_log()
                        except Exception as e:
                            print('Issue with Fastqc: ', sample)
                            print(e)

                    # if 'samsort' in checks:
                    #     print('SAMSORT')
                    #     samsort_class = SamSort(path=parent_path+'samsort/', sample=sample, table_path='data/fastq.csv')
                    #     try:
                    #         samsort_class.check_log()
                    #     except Exception as e:
                    #         print('Issue with SamSort: ', sample)
                    #         print(e)

                    # if 'baserecalibrator' in checks:
                    #     print('BaseRecalibrator')
                    #     baserecalibrator_class = BaseRecalibrator(path=parent_path + 'baserecalibrator/', sample=sample, 
                    #                                             table_path='data/fastq.csv')
                    #     try:
                    #         baserecalibrator_class.check_log(title=None, score=False)
                                    
                    #     except Exception as e:
                    #         print('Issue with BaseRecalibrator: ', sample)
                    #         print(e)

                    # if 'applybqsr' in checks:
                    #     print('ApplyBQSR')
                    #     applybqsr_class = ApplyBQSR(path=parent_path+'applybqsr/', sample=sample)
                    #     try:
                    #         applybqsr_class.check_log(title=None, score=False)
                                    
                    #     except Exception as e:
                    #         print('Issue with ApplyBQSR: ', sample)
                    #         print(e)

                    # if 'haplotype' in checks:
                    #     print('HaploType')
                    #     haplo_class = HaploType(path=parent_path+'haplotypecaller/', sample=sample, table_path='data/fastq.csv')
                    #     try:
                    #         haplo_class.check_log(warning_plot=False, title=None, score=False)
                                    
                    #     except Exception as e:
                    #         print('Issue with HaploType: ', sample)
                    #         print(e)


                    # break
        except PermissionError as e:
            print('No permissions for ', center)
            continue

        
    #     print('ApplyBQSR')
    #     applybqsr_class = ApplyBQSR(path=parent_path+'applybqsr/', sample=sample)
    #     try:
    #         score += applybqsr_class.check_log(title=None, score=True)
                    
    #     except Exception as e:
    #         print('Issue with ApplyBQSR: ', sample)
    #         print(e)
        
    #     print('HaploType')
    #     haplo_class = HaploType(path=parent_path+'haplotypecaller/', sample=sample, table_path='data/fastq.csv')
    #     try:
    #         score += haplo_class.check_log(warning_plot=False, title=None, score=True)
                    
    #     except Exception as e:
    #         print('Issue with HaploType: ', sample)
    #         print(e)
            
    #     list_.append([sample, round(score, 2)])

    #     print('===============')
    #     print('')






