import hailtop.batch as hb

billing_project = 'arnav-trial'
bucket = 'ttn-neb-analysis-files'
digest = 'sha256:d8a1a842a93adcb4dcf1e784d1743fd21fb84d2b404b5350c1b34bdbfd208264'
image = '12arnavg/repository@' + digest

backend = hb.ServiceBackend(billing_project = billing_project, bucket = bucket)
batch = hb.Batch(name = 'batch', backend = backend)
job = batch.new_bash_job(name = 'job')
job.image(image)
job.gcsfuse(bucket = bucket, mount_point = bucket)

in_path = 'ttn-neb-analysis-files/pipeline/'

in_bam = in_path + 'in.bam'

#region = 'chr2:()-()'
#job.command('samtools view {in_bam} region > {in_bam}')
job.command(f'samtools sort -n {in_bam} > in.sorted.bam')
job.command('bedtools bamtofastq -i in.sorted.bam -fq in.r1.fastq -fq2 in.r2.fastq')

in_fasta = in_path + 'in.fasta'

job.command(f'bwa mem {in_fasta} in.r1.fastq in.r2.fastq > out.sam')

job.command('samtools view -b out.sam > out.bam')
job.command('samtools sort out.bam > out.sorted.bam')
job.command('samtools index out.sorted.bam')

job.command('java -jar gatk-package-4.2.0.0-local.jar AddOrReplaceReadGroups -I out.sorted.bam -LB 1 -PL ILLUMINA -PU 1 -SM out --VALIDATION_STRINGENCY LENIENT -O out.sorted.rg.bam')
job.command('java -jar gatk-package-4.2.0.0-local.jar MarkDuplicates --CREATE_INDEX -I out.sorted.rg.bam --METRICS_FILE /dev/null -O out.sorted.deduped.bam')
job.command(f'java -jar gatk-package-4.2.0.0-local.jar HaplotypeCaller -I out.sorted.deduped.bam -O out.vcf.gz -R {in_fasta}')

out_path = 'gs://ttn-neb-analysis-files/pipeline/'

job.command(f'cp out.vcf.gz {job.out_vcf_gz}')
batch.write_output(job.out_vcf_gz, out_path + 'out.vcf.gz')
job.command(f'cp out.vcf.gz.tbi {job.out_vcf_gz_tbi}')
batch.write_output(job.out_vcf_gz_tbi, out_path + 'out.vcf.gz.tbi')

batch.run()