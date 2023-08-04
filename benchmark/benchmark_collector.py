## Quick code to collect all benchmark information

import glob
import os
import pandas as pd
import re
import argparse


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


def extract_info_list(id_, sample, path):

    l = []
    deli = '\t'

    list_files_fastqc = os.listdir(path+'fastqc/')
    list_benchmark = [i for i in list_files_fastqc if re.search('.benchmark$', i)]

    value = 0
    for i in list_benchmark:

        df = pd.read_csv(path+'fastqc/'+i, delimiter=deli)
        value += df.s.values[0]

    # ID
    l.append(id_)

    ## Fastqc
    l.append(round(value/60, 2))

    ## Bwa
    try:
        l.append(round(pd.read_csv(path+'bwa/'+sample+'.benchmark', delimiter=deli).s.values[0]/60, 2))
    except FileNotFoundError as e:
        l.append(0)

    ## SamBlaster
    try:
        l.append(round(pd.read_csv(path+'bwa/'+sample+'_samblaster.benchmark', delimiter=deli).s.values[0]/60, 2))
    except FileNotFoundError as e:
        l.append(0)

    ## SamSort
    try:
        l.append(round(pd.read_csv(path+'bwa/'+sample+'_sort_nodup.sam.benchmark', delimiter=deli).s.values[0]/60, 2))
    except FileNotFoundError as e:
        l.append(0)

    ## BaseRecalibrator
    try:
        l.append(round(pd.read_csv(path+'gatk_bsr/'+sample+'_sort_nodup.recaldat.benchmark', delimiter=deli).s.values[0]/60, 2))
    except FileNotFoundError as e:
        l.append(0)

    ## ApplyBQSR
    try:
        l.append(round(pd.read_csv(path+'gatk_bsr/'+sample+'_sort_nodup.bqsr.benchmark', delimiter=deli).s.values[0]/60, 2))
    except FileNotFoundError as e:
        l.append(0)

    ## Haplotype
    try:
        l.append(round(pd.read_csv(path+'gatk_gvcf/tmp_'+sample+'_sort_nodup.g.vcf.benchmark', delimiter=deli).s.values[0]/60, 2))
    except FileNotFoundError as e:
        l.append(0)

    ## File size
    size = round(os.path.getsize(path + sample + '_R1.fastq.gz')*1e-9, 2)
    size += round(os.path.getsize(path + sample + '_R1.fastq.gz')*1e-9, 2)
    l.append(size)

    return l


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='Collect Benchmark Information')
    parser.add_argument("--path", type=str)

    args = parser.parse_args()
    path = args.path

    centers = ['1a77728b-011a-4407-b60e-ec90b6430b99', 'ed615b8d-fa86-48bb-b0f8-c841e1aeb0eb', 'd7c41726-13a8-4abd-b185-68198fec12f4', 
               'fb203f69-94b9-42aa-9f34-c9ee8219a22e', 
               'eaabd6f1-1b34-4196-882c-198465045d71']


    data = [['id', 'fastqc', 'bwa', 'samblaster', 'samsort', 'base', 'apply', 'haplo', 'size (GB)', 'total time (min)']]

    for center in centers:
        for study in os.listdir(path + '/' + center + '/'):
            for sample in os.listdir(path + '/' + center + '/' + study + '/'):
                path_parent = path + '/' + center + '/' + study + '/' + sample + '/'
                try:
                    samples = [os.path.basename(x) for x in glob.glob(path_parent + '/*.gz')][0]
                    samples = samples[:-12]
                except IndexError as e:
                    print('Cannot find .gz files for this sample')
                    continue

                tmp = extract_info_list(id_=center + '||' + study + '||' + sample, sample = samples, path = path_parent)
                tmp.append(sum(tmp[1:8]))
                data.append(tmp)


    print('FINISHED')
    print('==================================')
    print('==================================')
    print('==================================')
    df = pd.DataFrame(data[1:], columns=data[0])
    bins = [0, 5, 10, 15, 20, 25]
    bins1 = [i for i in range(25)]
    df['size_bin'] = pd.cut(df['size (GB)'], bins)
    df['size_bin_small'] = pd.cut(df['size (GB)'], bins1)   
    
    x = df.groupby(['size_bin'])['total time (min)'].mean()
    plot = x.plot(kind='bar', xlabel='Size Bucket', ylabel='Processing Time (min)', title='Computational Processing \n Benchmark')
    fig = plot.get_figure()
    fig.savefig("output.png", bbox_inches='tight')
    # df.to_csv('res.csv')
    df1 = df[(df.size_bin_small == pd.Interval(5, 6, closed='right')) | 
             (df.size_bin_small == pd.Interval(10, 11, closed='right')) | 
             (df.size_bin_small == pd.Interval(22, 23, closed='right'))]
    print(df1.head())
    x = df1.groupby(['size_bin_small'])[['fastqc', 'bwa', 'samsort', 'base', 'apply', 'haplo']].mean().dropna()
    print(x)
    plot = x.T.plot(kind='bar')
    fig = plot.get_figure()
    fig.savefig("output1.png", bbox_inches='tight')


    

