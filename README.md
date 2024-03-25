# nf-CHIPseeker

```
PROFILE=raven
nextflow run nf-CHIPseeker -params-file ./nf-CHIPseeker/params.slurm.json -entry images -profile ${PROFILE}  && \
nextflow run nf-CHIPseeker -params-file ./nf-CHIPseeker/params.slurm.json -profile ${PROFILE}
```