import hailtop.batch as hb
import pandas as pd;
import hail as hl
import argparse
import json

hl.init(log = '/dev/null')

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', required = True)
parser.add_argument('-m', '--masks', nargs = '+', required = True)
parser.add_argument('-o', '--output', required = True)
args = parser.parse_args()

core_img = '12arnavg/repository@sha256:e42aa3bd202f827b840579e8b0318eb51958088516f58b6af47b68b7a14f2632'
nirvana_img, working_bkt, billing_project = 'annotation/nirvana:3.14', 'ttn-neb-analysis-files', 'arnav-trial'

'''backend = hb.ServiceBackend(billing_project = billing_project, bucket = working_bkt)
batch = hb.Batch(name = 'batch', backend = backend)'''

regions = ['chr2:151483334-151736474', 'chr2:178523989-178809423'] # NEB, TTN
repetitive_regions = ['chr2:151579109-151588808', 'chr2:151588872-151599459', 'chr2:151599478-151609050',
    'chr2:178653204-178655062', 'chr2:178656913-178659722', 'chr2:178658723-178663615'] # 3 NEB, 3 TTN
lofs = ['splice_acceptor_variant', 'splice_donor_variant', 'stop_gained', 'frameshift_variant']

in_path = 'gs://' + working_bkt + '/input/' + args.input
in_csv = hl.hadoop_open(in_path)
in_df = pd.read_csv(in_csv)
'''
for index, row in in_df.iterrows():

    cram_path = row['cram_path']
    cram_full_name = cram_path[5:][:-5]
    cram_arr = cram_full_name.split('/')
    cram_bkt, cram_name = cram_arr[0], cram_arr[-1]

    for mask in args.masks:

        vcf_path = 'gs://' + working_bkt + '/output/' + mask + '/' + cram_full_name + '.vcf.gz'

        if not hl.hadoop_is_file(vcf_path):

            vcf_job = batch.new_bash_job(name = 'vcf_job')
            vcf_job.image(core_img)

            vcf_job.gcsfuse(bucket = cram_bkt, mount_point = cram_bkt)
            vcf_job.gcsfuse(bucket = working_bkt, mount_point = working_bkt)

            ref_path = working_bkt + '/input/Homo_sapiens_assembly38.fasta'
            regions_str = ' '.join(regions)
            vcf_job.command(f'samtools view -b {cram_full_name}.cram -T {ref_path} {regions_str} > {cram_name}.bam')

            vcf_job.command(f'samtools index {cram_name}.bam')
            vcf_job.command(f'samtools sort -n {cram_name}.bam > {cram_name}.sorted1.bam')
            vcf_job.command(f'bedtools bamtofastq -i {cram_name}.sorted1.bam -fq {cram_name}.r1.fastq -fq2 {cram_name}.r2.fastq')

            mask_path = working_bkt + '/input/Homo_sapiens_assembly38.chr2.' + mask + '.fasta'

            vcf_job.command(f'bwa mem {mask_path} {cram_name}.r1.fastq {cram_name}.r2.fastq > {cram_name}.sam')

            vcf_job.command(f'samtools view -b {cram_name}.sam > {cram_name}.bam')
            vcf_job.command(f'samtools sort {cram_name}.bam > {cram_name}.sorted2.bam')
            vcf_job.command(f'samtools index {cram_name}.sorted2.bam')

            vcf_job.command(f'java -jar gatk-package-4.2.0.0-local.jar AddOrReplaceReadGroups -I {cram_name}.sorted2.bam -LB 1 -PL ILLUMINA -PU 1 -SM out --VALIDATION_STRINGENCY LENIENT -O {cram_name}.sorted.rg.bam')
            vcf_job.command(f'java -jar gatk-package-4.2.0.0-local.jar MarkDuplicates --CREATE_INDEX -I {cram_name}.sorted.rg.bam --METRICS_FILE /dev/null -O {cram_name}.sorted.deduped.bam')
            vcf_job.command(f'java -jar gatk-package-4.2.0.0-local.jar HaplotypeCaller -I {cram_name}.sorted.deduped.bam -O {cram_name}.vcf.gz -R {mask_path}')

            vcf_job.command(f'cp {cram_name}.vcf.gz {vcf_job.vcf}')
            batch.write_output(vcf_job.vcf, vcf_path)
            vcf_job.command(f'cp {cram_name}.vcf.gz.tbi {vcf_job.vcf_index}') 
            batch.write_output(vcf_job.vcf_index, vcf_path + '.tbi')
    
        json_path = 'gs://' + working_bkt + '/output/' + mask + '/' + cram_full_name + '.json.gz'

        if not hl.hadoop_is_file(json_path):

            json_job = batch.new_bash_job(name = 'json_job')
            json_job.image(nirvana_img)

            json_job.gcsfuse(bucket = working_bkt, mount_point = working_bkt)

            json_job.command(f'dotnet /opt/nirvana/Nirvana.dll -c {working_bkt}/Nirvana/Data/Cache/GRCh38/Both -r {working_bkt}/Nirvana/Data/References/Homo_sapiens.GRCh38.Nirvana.dat --sd {working_bkt}/Nirvana/Data/SupplementaryAnnotation/GRCh38 -i {vcf_path[5:]} -o {cram_name}')

            json_job.command(f'cp {cram_name}.json.gz {json_job.json}')
            batch.write_output(json_job.json, json_path)
            json_job.command(f'cp {cram_name}.json.gz.jsi {json_job.json_index}')
            batch.write_output(json_job.json_index, json_path + '.jsi')

            json_job.depends_on(vcf_job)

batch.run()'''

out_arr = []

for index, row in in_df.iterrows():

    for mask in args.masks:

        cram_path = row['cram_path']
        cram_full_name = cram_path[5:][:-5]

        json_path = 'gs://' + working_bkt + '/output/' + mask + '/' + cram_full_name + '.json.gz'
        json_dict = json.loads(hl.hadoop_open(json_path).read())

        for position in json_dict['positions']:

            variants = position['variants']

            for i in range(len(variants)):

                variant = variants[i]

                consequences = []

                if 'transcripts' in variant:

                    for transcript in variant['transcripts']:

                        consequences.extend(transcript['consequence'])
                
                consequences = list(set(consequences))

                is_lof = False
                lof_consequences = []

                for consequence in consequences:

                    if consequence in lofs:

                        is_lof = True
                        lof_consequences.append(consequence)

                out_line = [cram_path, position['chromosome'], position['position'], position['samples'][0]['genotype'], position['refAllele'], 
                    position['altAlleles'][i], variant['variantType'], consequences, lof_consequences, is_lof]

                neb_arr = regions[0].split(':')[1].split('-')

                if out_line[2] >= int(neb_arr[0]) and out_line[2] <= int(neb_arr[1]):
                    gene = 'NEB'
                else:
                    gene = 'TTN'
                
                in_repetitive_regions = False

                for repetitive_region in repetitive_regions:

                    region_arr = repetitive_region.split(':')[1].split('-')
    
                    if out_line[2] >= int(region_arr[0]) and out_line[2] <= int(region_arr[1]):

                        in_repetitive_regions = True
                
                clinvar_significances = []

                if 'clinvar' in variant:

                    for clinvar in variant['clinvar']:

                        clinvar_significances.extend(clinvar['significance'])
                
                clinvar_significances = list(set(clinvar_significances))

                gnomad_coverage = 0
                gnomad_allAf = 0
    
                if 'gnomad' in variant:

                    gnomad_coverage = variant['gnomad']['coverage']
                    gnomad_allAf = variant['gnomad']['allAf']
                
                out_line.extend([mask, gene, in_repetitive_regions, clinvar_significances, gnomad_coverage, gnomad_allAf])
                out_arr.append(out_line)

out_df = pd.DataFrame(out_arr, columns = ['cram_path', 'chromosome', 'position', 'genotype', 'ref_allele', 'alt_allele', 'variant_type', 'consequences',
    'lof_consequences', 'is_lof', 'mask', 'gene', 'in_repetitive_regions', 'clinvar_significances', 'gnomad_coverage', 'gnomad_allAf'])

out_df.to_csv(args.output, index = False)