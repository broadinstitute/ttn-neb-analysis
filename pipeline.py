import hailtop.batch as hb
import pandas as pd
import argparse
import hail as hl

#
parser = argparse.ArgumentParser()
parser.add_argument('-r', '--region')
args = parser.parse_args()

#
billing = 'arnav-trial' 
bucket = 'ttn-neb-analysis-files'
digest = 'sha256:e42aa3bd202f827b840579e8b0318eb51958088516f58b6af47b68b7a14f2632'
image = '12arnavg/repository@' + digest

#
backend = hb.ServiceBackend(billing_project = billing, bucket = bucket)
batch = hb.Batch(name = 'batch', backend = backend)

#
paths_path = 'cram_paths.tsv'
df = pd.read_table(paths_path)

for i in df.index:
    

    #
    cram_path = df['cram_path'][i][5:]
    cram_arr = cram_path.split('/')
    cram_bucket = cram_arr[0]
    cram_name = cram_arr[-1][:-5]

    vcf_path = 'gs://' + bucket + '/output/' + cram_path[:-5] + '.vcf.gz'

    if not hl.hadoop_is_file(vcf_path):

        #
        job = batch.new_bash_job(name = 'job')
        job.image(image)

        #
        job.gcsfuse(bucket = bucket, mount_point = bucket)
        job.gcsfuse(bucket = cram_bucket, mount_point = cram_bucket)

        #
        job.command(f'samtools view -b {cram_path} {args.region} > {cram_name}.bam')
        job.command(f'samtools index {cram_name}.bam')
        job.command(f'samtools sort -n {cram_name}.bam > {cram_name}.sorted1.bam')
        job.command(f'bedtools bamtofastq -i {cram_name}.sorted1.bam -fq {cram_name}.r1.fastq -fq2 {cram_name}.r2.fastq')

        #
        fasta_path = bucket + '/input/Homo_sapiens_assembly38.chr2.fasta'
        job.command(f'bwa mem {fasta_path} {cram_name}.r1.fastq {cram_name}.r2.fastq > {cram_name}.sam')

        #
        job.command(f'samtools view -b {cram_name}.sam > {cram_name}.bam')
        job.command(f'samtools sort {cram_name}.bam > {cram_name}.sorted2.bam')
        job.command(f'samtools index {cram_name}.sorted2.bam')

        #
        job.command(f'java -jar gatk-package-4.2.0.0-local.jar AddOrReplaceReadGroups -I {cram_name}.sorted2.bam -LB 1 -PL ILLUMINA -PU 1 -SM out --VALIDATION_STRINGENCY LENIENT -O {cram_name}.sorted.rg.bam')
        job.command(f'java -jar gatk-package-4.2.0.0-local.jar MarkDuplicates --CREATE_INDEX -I {cram_name}.sorted.rg.bam --METRICS_FILE /dev/null -O {cram_name}.sorted.deduped.bam')
        job.command(f'java -jar gatk-package-4.2.0.0-local.jar HaplotypeCaller -I {cram_name}.sorted.deduped.bam -O {cram_name}.vcf.gz -R {fasta_path}')
    
        #
        job.command(f'cp {cram_name}.vcf.gz {job.vcf_gz}')
        batch.write_output(job.vcf_gz, vcf_path)
        job.command(f'cp {cram_name}.vcf.gz.tbi {job.vcf_gz_tbi}') 
        batch.write_output(job.vcf_gz_tbi, vcf_path + '.tbi')

batch.run()
