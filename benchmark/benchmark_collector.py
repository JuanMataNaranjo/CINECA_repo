## Quick code to collect all benchmark information

import os
import pandas as pd
import re

def extract_information(sample, path=''):
    """
    Function to extract information for a given sample
    :params path: parent directory from which all the remaining directories can be found
    :params sample: name of the sample so that we can extract info regarding all those samples
    :return: The output is a dictionary for each rule with the associated **minutes** that each rule required to execute
    """
    
    benchmark_dict = {}

    list_files_fastqc = os.listdir(path+'fastqc_out/')
    list_benchmark = [i for i in list_files_fastqc if re.search('.benchmark$', i)]

    value = 0
    for i in list_benchmark:
            
        df = pd.read_csv(path+'fastqc_out/'+i, delimiter=' ')
        value += df.s.values[0]

    value_seqtk = 0
    for i in os.listdir(path+'OUTPUT3/'+sample+'/seqtk/'):

        df = pd.read_csv(path+'OUTPUT3/'+sample+'/seqtk/'+i, delimiter=' ')
        value_seqtk += df.s.values[0]

    ## Fastqc
    benchmark_dict['fastqc'] = round(value/60, 2)

    ## Seqtk
    benchmark_dict['Seqtk'] = round(value_seqtk/60, 2) 

    ## Bwa
    benchmark_dict['Bwa'] = round(pd.read_csv(path+'OUTPUT3/'+sample+'/bwa/'+sample+'.benchmark', delimiter=' ').s.values[0]/60, 2)

    ## SamBlaster
    benchmark_dict['SamBlaster'] = round(pd.read_csv(path+'OUTPUT3/'+sample+'/bwa/'+sample+'_samblaster.benchmark', delimiter=' ').s.values[0]/60, 2)

    ## SamSort
    benchmark_dict['SamSort'] = round(pd.read_csv(path+'OUTPUT3/'+sample+'/bwa/'+sample+'_sort_nodup.sam.benchmark', delimiter=' ').s.values[0]/60, 2)

    ## SamIndex
    benchmark_dict['SamIndex'] = round(pd.read_csv(path+'OUTPUT3/'+sample+'/bwa/'+sample+'_sort_nodup.benchmark', delimiter=' ').s.values[0]/60, 2)

    ## BaseRecalibrator
    benchmark_dict['BaseRecalibrator'] = round(pd.read_csv(path+'OUTPUT3/'+sample+'/gatk_bsr/'+sample+'_sort_nodup.recaldat.benchmark', delimiter=' ').s.values[0]/60, 2)

    ## ApplyBQSR
    benchmark_dict['ApplyBQSR'] = round(pd.read_csv(path+'OUTPUT3/'+sample+'/gatk_bsr/'+sample+'_sort_nodup.bqsr.benchmark', delimiter=' ').s.values[0]/60, 2)

    ## HaplotypeCaller
    benchmark_dict['HaploType'] = round(pd.read_csv(path+'OUTPUT3/'+sample+'/gatk_gvcf/'+sample+'_sort_nodup.g.vcf.benchmark', delimiter=' ').s.values[0]/60, 2)
            
    return benchmark_dict


if __name__ == "__main__":
    
    print('Time is reported in MINUTES')
    df = pd.DataFrame.from_dict(extract_information(sample='HSRR062634'), orient='index', columns=['HSRR062634'])
    print(df)

