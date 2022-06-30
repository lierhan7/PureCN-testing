#!/usr/bin/env python


import os
import os.path
from os import listdir
from os.path import isfile, join
import subprocess
import argparse
from joblib import Parallel, delayed
import pyranges as pr   # header needs to match as Chromosome Start End
import pandas as pd
import glob

__author__ = 'Erhan Li'
__contact__ = 'leh@maijinggene.com'
__version__ = '1.0.0'
__doc__ = '''PureCN wrapper to run in parallel to evaluate germline/somatic mutations and CNA.'''


def parallel_run(output_dir, sn, rscript, purecn_exe, bam_files_dir, sn_bam_dict, intervals_output_file):
    """
    creates files:
    _coverage_loess.png
    _coverage_loess.txt
    _coverage_loess_qc.txt
    _coverage.txt

    """
    coverage_command = [rscript, purecn_exe + '/Coverage.R',
                       '--force',
                       '--out-dir', output_dir + '/' + sn,
                       '--bam', '/'.join([bam_files_dir, sn_bam_dict[sn]]),
                       '--intervals', intervals_output_file]
    print(' '.join(coverage_command))

    out, error = subprocess.Popen(" ".join(coverage_command), stdout=subprocess.PIPE, shell=True).communicate()
    if error:
        print("An error occured: {error}".format(error=error))

def parallel_run_pureCN_SGZ(output_dir, sn, rscript, purecn_exe, sgz_dir, input_bed_file, tumorcov, normaldb_file, mapping_bias, genome, tumorvcf, tumorstats, intervals_output_file, snp_blacklist):

    pureCN_command = [rscript, purecn_exe + '/PureCN.R',
                       '--out', output_dir + '/' + sn + '/',
                       '--tumor', tumorcov,
                       '--sampleid', sn,
                       '--vcf', tumorvcf,
                       '--stats-file', tumorstats,
                       '--fun-segmentation PSCBS',
                       '--normaldb', normaldb_file,
                       '--mapping-bias-file', mapping_bias,
                       '--intervals', intervals_output_file,
                       '--snp-blacklist', snp_blacklist,
                       '--max-copy-number 8',
                       '--genome', genome,
                       '--model betabin',
                       '--force --post-optimize --seed 123']

    print(' '.join(pureCN_command))

    out, error = subprocess.Popen(" ".join(pureCN_command), stdout=subprocess.PIPE, shell=True).communicate()
    if error:
        print("An error occured: {error}".format(error=error))


    sgz_command = [rscript, purecn_exe + '/TumorOnlySGZ_v01.R',
                       '--rds', output_dir + '/' + sn + '/' + sn +'.rds',
                       '--out', output_dir + '/' + sn + '/' + sn,
                       '--vcf', tumorvcf,
                       '--sgzdir', sgz_dir,
                       '--callable', input_bed_file]
    print("sgz_cmd: " + " ".join(sgz_command))
    out, error = subprocess.Popen(" ".join(sgz_command), stdout=subprocess.PIPE, shell=True).communicate()
    if error:
        print("An error occured: {error}".format(error=error))


def __main__():

    parser = argparse.ArgumentParser(description='Parse Carlsbad Archer VCs and output into single file')
    parser.add_argument('-t', '--tumor_bam_dir',  help="Directory containing tumor bam files", required=True)
    parser.add_argument('-n', '--normal_bam_dir',  help="Directory containing normal bam files", required=False)
    parser.add_argument('-b', '--bed',  help="bed file", required=False)
    parser.add_argument('-v', '--tumor_vcf_dir',  help="Directory containing tumor vcf files", required=True)
    parser.add_argument('-p', '--normal_panel',  help="normal vcf file", required=False)
    parser.add_argument('-o', '--output_dir',  help="Directory pipe outputs", required=True)

    args = parser.parse_args()

    rscript = '/home/lierhan/Program/miniconda3/bin/Rscript'
    purecn_exe = '/home/lierhan/Program/miniconda3/lib/R/library/PureCN/extdata'
    sgz_dir = '/home/lierhan/Program/SGZ-master'

    ref_fa = '/extend2/Share/Databases/MY_hg19/hg19.fa'
    genome = 'hg19'
    mappability_file = '/home/lierhan/Project/PureCN_20220418/20220428/DB/wgEncodeCrgMapabilityAlign100mer.bigWig'
    reptiming_file = '/home/lierhan/Project/PureCN_20220418/20220428/DB/wgEncodeUwRepliSeqK562WaveSignalRep1.bigWig'
    intervals_output_file =args.output_dir +'/baits_' + genome + '_intervals.txt'
    snp_blacklist = '/home/lierhan/Project/PureCN_20220418/20220428/DB/snp_blacklist/hg19_simpleRepeats.bed'

    if not args.bed:
        input_bed_file = '/isilon/RnD/tools/custom_script/cnvkit/bed/QIAseq323.CDHS-24104Z-14204.refseq-anno.roi.exons.nolowcovrois_genes_only.bed' #Note, had to remove the header line to make work
    else:
        input_bed_file = args.bed

    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)

    print('------------------------------\n| ** generate an interval ** |\n------------------------------\n')
    prepare_command = [rscript, purecn_exe + '/IntervalFile.R',
                        '--force',
                        '--in-file', input_bed_file,
                        '--fasta', ref_fa,
                        '--out-file', intervals_output_file,
                        '--off-target',
                        '--genome', genome,
                        '--export', args.output_dir + '/baits_optimized' + genome + '.bed',
                        '--mappability', mappability_file,
                        '--reptiming', reptiming_file]

    cmd1 = " ".join(prepare_command)
    print('cmd: {cmd}'.format(cmd=cmd1))
    out, error = subprocess.Popen(cmd1, stdout=subprocess.PIPE, shell=True).communicate()
    if error:
        print("An error occured: {error}".format(error=error))

    print('** finished  ** |\n------------------\n')

    EXT = (".vcf.gz",".vcf")
    if args.tumor_vcf_dir:
        tumorvcffiles = [f for f in sorted(listdir(args.tumor_vcf_dir)) if isfile(join(args.tumor_vcf_dir, f)) and f.endswith(EXT)]
        tumorstatsfiles = [ f for f in sorted(listdir(args.tumor_vcf_dir)) if isfile(join(args.tumor_vcf_dir, f)) and f.endswith(".stats")]

    print('--------------------------------\n| ** generate tumor coverage ** |\n--------------------------------\n')
    # 4. PureCN Coverage
    if args.tumor_bam_dir:
        tumorbamfiles = [f for f in sorted(listdir(args.tumor_bam_dir)) if isfile(join(args.tumor_bam_dir, f)) and f.endswith('.bam')]
        tumor_sample_list = []
        tumor_sn_bam_dict = {}
        tumor_coverage_list = []

        for bam in tumorbamfiles:
            sn = os.path.basename(bam).split('.')[0]
            tumor_sample_list.append(sn)
            tumor_sn_bam_dict[sn] = bam
            if not os.path.exists(args.output_dir + '/' + sn):
                os.makedirs(args.output_dir + '/' + sn)

            tumor_coverage_list.append(args.output_dir + '/' + sn + '/' + sn + '_coverage_loess.txt.gz')

        # determine the number of CPUs
        if len(tumorbamfiles) > 20:
            NCPU = 20
        else:
            NCPU=int(len(tumorbamfiles))

        Parallel(n_jobs=NCPU)(delayed(parallel_run)(output_dir=args.output_dir, \
                              sn=tumor_sample_list[i], rscript=rscript, purecn_exe=purecn_exe, bam_files_dir=args.tumor_bam_dir, \
                              sn_bam_dict=tumor_sn_bam_dict, intervals_output_file=intervals_output_file) for i in range(0,len(tumor_sample_list)))

    print('--------------------------------\n| ** finished generating tumor coverage ** |\n--------------------------------\n')

    if args.normal_bam_dir:
        print('--------------------------------\n| ** generate normal coverage ** |\n--------------------------------\n')
        normalbamfiles = [f for f in sorted(listdir(args.normal_bam_dir)) if isfile(join(args.normal_bam_dir, f)) and f.endswith('.bam')]
        normal_sample_list = []
        normal_sn_bam_dict = {}
        normal_coverage_list = []
        for bam in normalbamfiles:
            sn = os.path.basename(bam).split('.')[0]
            normal_sample_list.append(sn)
            normal_sn_bam_dict[sn] = bam
            if not os.path.exists(args.output_dir + '/' + sn):
                os.makedirs(args.output_dir + '/' + sn)

            normal_coverage_list.append(args.output_dir + '/' + sn + '/' + sn + '_coverage_loess.txt.gz')

        if normal_coverage_list:
            normal_coverage_list_file_name = args.output_dir + '/normal_coverage_list_file.txt'
            normal_coverage_list_file = open(normal_coverage_list_file_name, 'w')
            for i in normal_coverage_list:
                normal_coverage_list_file.write(i + '\n')
            normal_coverage_list_file.flush()
            normal_coverage_list_file.close()

        Parallel(n_jobs=NCPU)(delayed(parallel_run)(output_dir=args.output_dir, \
                              sn=normal_sample_list[i], rscript=rscript, purecn_exe=purecn_exe, bam_files_dir=args.normal_bam_dir, \
                              sn_bam_dict=normal_sn_bam_dict, intervals_output_file=intervals_output_file) for i in range(0,len(normal_sample_list)))
        print('--------------------------------\n| ** finished generating normal coverage ** |\n--------------------------------\n')

        print('--------------------------------\n| ** generate normalDB ** |\n--------------------------------\n')
        normal_panel = args.normal_panel    # normal vcf

        cmd2 = '{Rscript} {PURECN}/NormalDB.R --out-dir {outdir} --coverage-files {outdir}/normal_coverage_list_file.txt --force --genome hg19 --normal-panel {normalpanel}'.format(Rscript=rscript, PURECN = purecn_exe, outdir=args.output_dir, normalpanel=normal_panel)
        print('cmd:{cmd2}'.format(cmd2=cmd2))
        out, error = subprocess.Popen(cmd2, stdout=subprocess.PIPE, shell=True).communicate()
        if error:
            print("An error occured: {error}".format(error=error))

        normaldb_file = args.output_dir + '/' + 'normalDB_hg19.rds'
        mapping_bias =  args.output_dir + '/' + 'mapping_bias_hg19.rds'
    else:
        #recycle files
        normaldb_file = '/isilon/RnD/tools/custom_script/purecn/mapping_bias/qiaseq323/pon_min1_qiaseq323_14normals_normalDB_hg19.rds'
        mapping_bias =  '/isilon/RnD/tools/custom_script/purecn/mapping_bias/qiaseq323/pon_min1_qiaseq323_14normals_mapping_bias_hg19.rds'

    # # 5. Run PureCN for Tumor Files
    tumorvcfs = []
    tumorstats_list = []
    for v in tumorvcffiles:
        tumorvcfs.append(args.tumor_vcf_dir+'/'+ v)
    for s in tumorstatsfiles:
        tumorstats_list.append(args.tumor_vcf_dir+'/'+s)
    print('--------------------------------\n| ** run PureCN followed by SGZ ** |\n--------------------------------\n')
    if len(tumorbamfiles) == len(tumorvcffiles):
        Parallel(n_jobs=NCPU)(delayed(parallel_run_pureCN_SGZ)(output_dir=args.output_dir, \
                          sn=os.path.basename(tumor_coverage_list[i]).split('_coverage')[0], rscript=rscript, purecn_exe=purecn_exe, sgz_dir=sgz_dir, \
                          input_bed_file=input_bed_file, tumorcov=tumor_coverage_list[i], normaldb_file=normaldb_file, mapping_bias=mapping_bias, genome=genome, tumorvcf=tumorvcfs[j], \
                          tumorstats=tumorstats_list[k], intervals_output_file=intervals_output_file, snp_blacklist=snp_blacklist) for i,j,k in zip(range(0,len(tumor_coverage_list)),range(0,len(tumorvcfs)),range(0,len(tumorstats_list))))

        print('--------------------------------\n| ** finished runing PureCN and SGZ ** |\n--------------------------------\n')


        '''
        create a final report
        '''
        print('--------------------------------\n| ** generate a final report for PureCN ** |\n--------------------------------\n')

        fields1 = ['Sampleid','chr', 'start','end','arm','C','M','type','seg.mean','M.flagged']
        fields2 = ['Sampleid','chr', 'start','end','ID','REF','ALT','ML.SOMATIC','POSTERIOR.SOMATIC','ML.LOH','log.ratio','depth','prior.somatic','gene.symbol']
        outputF = '{out}/{out}_purecn_final_report.xlsx'.format(out=args.output_dir)
        x = []

        flag=0
        print('--------------------------------\n| generate a final report |\n--------------------------------\n')
        with pd.ExcelWriter(outputF, mode='w') as writer:
            lohFiles=sorted(glob.glob('{}/[Acc,NTP,MOL]*/*_loh.csv'.format(args.output_dir)))
            genesFiles=sorted(glob.glob('{}/[Acc,NTP,MOL]*/*_variants.csv'.format(args.output_dir)))

            for lohF, geneF in zip(lohFiles,genesFiles):
                sn1 = '-'.join(os.path.basename(geneF).split('-')[0:3])
                try:
                    df1 = pd.read_csv(lohF, delimiter=',', comment='#', usecols=fields1, skip_blank_lines=True,encoding='utf-8')
                    df2 = pd.read_csv(geneF, delimiter=',', comment='#', usecols=fields2, skip_blank_lines=True,encoding='utf-8')

                    if df1.shape[0]!=0:
                        df1 = df1[(df1['C']!=2) | (df1['M']==0)]  # exclude CN=2 and yet keep whoe neutral arm
                        df1.sort_values(by=['chr'])
                        #df2 = df2[(df2['prior.somatic']>0.1)]   # dbsnp includes some somatic so hold off applying now
                        x.append(df1)
                        if flag==0:
                            df1.to_excel(writer,sheet_name='master',index=False)
                            flag=flag+1

                        # create PyRanges-objects from the dfs
                        # pr requires a certain format with column names Chromosome Start End
                        df1.rename(columns={'chr':'Chromosome','start':'Start','end':'End'},inplace=True)
                        df2.rename(columns={'chr':'Chromosome','start':'Start','end':'End'},inplace=True)

                        gr1, gr2 = pr.PyRanges(df1),pr.PyRanges(df2)

                        # intersect the two
                        gr = gr2.intersect(gr1)

                        gr.df.to_excel(writer, sheet_name=sn1, index=False)
                except:
                    print('sample {sn} is not available'.format(sn=sn1))
            x1 = pd.concat(x)
            x1.to_excel(writer, sheet_name='master',index=False)

if __name__=="__main__": __main__()
