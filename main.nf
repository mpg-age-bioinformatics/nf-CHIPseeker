#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process get_images {
  stageInMode 'symlink'
  stageOutMode 'move'

  script:
    """

    if [[ "${params.run_type}" == "r2d2" ]] || [[ "${params.run_type}" == "raven" ]] || [[ "${params.run_type}" == "studio" ]]; 
      then
        cd ${params.image_folder}
        if [[ ! -f chipseeker-3.18.0.sif ]] ;
          then
            singularity pull chipseeker-3.18.0.sif docker://index.docker.io/mpgagebioinformatics/chipseeker:3.18.0
        fi
    fi


    if [[ "${params.run_type}" == "local" ]] ; 
      then
        docker pull mpgagebioinformatics/chipseeker:3.18.0
    fi

    """

}

process chipseeker_R {
  stageInMode 'symlink'
  stageOutMode 'move'

  when:
    ( ! file("${params.diffbind_out}/tmp/peakAnno.Rdata").exists() ) 
  
  script:
    """
#!/usr/bin/Rscript

setwd("${params.diffbind_out}/")
print(Sys.getenv("R_LIBS_USER")) # update
print("loading libraries")
library(ChIPseeker)
library(GenomicRanges)
library(ReactomePA)
library(openxlsx)

if (!dir.exists('${params.diffbind_out}/tmp')) {
  dir.create('${params.diffbind_out}/tmp')
}

stat_results = list()
if("${params.TxDb}" != ""){

library(${params.TxDb})
txdb <- ${params.TxDb}

all_files = list.files('.', pattern = 'all_peaks')
for(f in all_files){
  if(grepl('.tsv', f)){
    print(paste('Annotate peaks in', f))
    comp_name = gsub('.tsv', '', gsub('all_peaks.', '', f))
    
    D = read.delim(f, as.is =T)
    D.gr = GRanges(seqnames =  D[,'seqnames'],
               ranges = IRanges(D[,'start'], D[,'end']),
               strand = D[,'strand'],
               peak_width  = D[,'width'],
               Conc = D[,'Conc'],
               Fold = D[,'Fold'],
               FDR = D[,'FDR']
               )
    
    if ("${params.organism}" != "c_albicans_sc5314") {
        seqlevelsStyle(D.gr) <- "UCSC"
        }
    stat_results[[comp_name]] = D.gr
    promoter <- getPromoters(TxDb=txdb, upstream=2000, downstream=2000)
    tagMatrix <- getTagMatrix(D.gr, windows=promoter, weightCol = 'Conc')
    
    if (!is.null(dim(tagMatrix))){
      pdf(paste0('tagMatrix.', comp_name, '.pdf'))
    
        tagHeatmap(tagMatrix, xlim=c(-2000, 2000), color="red")
        plotAvgProf(tagMatrix, xlim=c(-2000, 2000),
              xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
        plotAvgProf(tagMatrix, xlim=c(-2000, 2000), conf = 0.95, resample = 1000)
      dev.off()
    }
    
    peakAnno <- annotatePeak(D.gr, tssRegion=c(-2000, 200),
                         TxDb=txdb, annoDb="${params.annoDb}")
    peakAnno
    
    pdf(paste0('annotated.', comp_name, '.pdf'), width = 10, height = 7)
      plotAnnoPie(peakAnno)
      plotDistToTSS(peakAnno,
              title="Distribution of transcription factor-binding loci\nrelative to TSS")
    dev.off()
    
    write.table(as.data.frame(peakAnno), paste0('annotated_peaks.', comp_name, '.tsv'), sep = '\t', quote = F, row.names = F) 
    write.xlsx(as.data.frame(peakAnno), paste0('annotated_peaks.', comp_name, '.xlsx'), row.names = F)
    library(ggplot2)
    peaks = as.data.frame(peakAnno)
    peaks[,'peak.group'] = ifelse(grepl('Exon', peaks[, 'annotation']), 'Exon', as.character(peaks[,'annotation']))
    peaks[,'peak.group'] = ifelse(grepl('Intron', peaks[,'peak.group']), 'Intron', as.character(peaks[, 'peak.group']))
    pdf(paste0('feature_foldChange_violin.', comp_name, '.pdf'), height = 7, width = 14)
    p = ggplot(peaks, aes(x = peak.group, y = Fold, fill = peak.group)) +
      geom_violin() +
      ggtitle(comp_name) + geom_hline(yintercept = 0, lwd = 0.7, lty = 2) + 
      stat_summary(geom = 'crossbar', size = 0.1, fun = 'median') 
    print(p)
    
    dev.off()
    
    if ("${params.organism}" != "c_albicans_sc5314") {
        # functional enrichment
        pathway1 <- enrichPathway(as.data.frame(peakAnno)[, 'geneId'], organism ='${params.peakAnno_organism}')
        if(nrow(as.data.frame(pathway1)) >= 1){
          pdf(paste0('all_peak_functional_enrichment.', comp_name, '.pdf'), width = 12, height = 6)
            dp = dotplot(pathway1)
            print(dp)
          dev.off()
        }
        write.table(as.data.frame(pathway1), paste0('enriched_reactome_using_all_peaks', comp_name, '.tsv'), sep = '\t', quote = F, row.names = F)
        write.xlsx(as.data.frame(pathway1), paste0('enriched_reactome_using_all_peaks', comp_name, '.xlsx'), row.names = F)
        # functional enrichment
        pathway1 <- enrichPathway(subset(as.data.frame(peakAnno), FDR <= 0.05)[, 'geneId'], organism = '${params.peakAnno_organism}')
        if(nrow(as.data.frame(pathway1)) >= 1){
          pdf(paste0('significant_peak_functional_enrichment.', comp_name, '.pdf'), width = 12, height = 6)
            dp = dotplot(pathway1)
            print(dp)
          dev.off()
        }
        write.table(as.data.frame(pathway1), paste0('enriched_reactome_using_significant_peaks', comp_name, '.tsv'), sep = '\t', quote = F, row.names = F)
        write.xlsx(as.data.frame(pathway1), paste0('enriched_reactome_using_significant_peaks', comp_name, '.xlsx'), row.names = F)
    }
  }
}
# annotate consensus peaks and add stats for master table
D = read.delim('consensus_peaks.tsv', as.is =T)
D.gr = makeGRangesFromDataFrame(D, keep.extra.columns=TRUE, seqnames.field = 'chrom_names', start.field = 'START', end.field = 'END')
if ("${params.organism}" != "c_albicans_sc5314") {
    seqlevelsStyle(D.gr) <- "UCSC"
    }
peakAnno <- annotatePeak(D.gr, tssRegion=c(-2000, 200),
                         TxDb=txdb, annoDb="${params.annoDb}")
DAnno = as.data.frame(peakAnno)
DAnno[, "merge_id"] = paste(DAnno[, 'seqnames'], DAnno[,'start'], DAnno[, 'end'])
if(length(stat_results) > 0){
  for(i in names(stat_results)){
    tmp = as.data.frame(stat_results[[i]])
    tmp[, "merge_id"] = paste(tmp[, 'seqnames'], tmp[,'start'], tmp[, 'end'])
    tmp = tmp[, c('merge_id', 'Conc', 'Fold', 'FDR')]
    names(tmp)[2:4] = paste(names(tmp)[2:4], i, sep = "|")
    DAnno = merge(DAnno, tmp, by = "merge_id", all.x = TRUE)
  }
}
DAnno = DAnno[,-1]
write.xlsx(DAnno, "annotated_master_table.xlsx", row.names = FALSE)
 
} else {
  print("THERE IS NO ANNOTATION DATABASE FOR ${params.organism}")
}
save.image("${params.diffbind_out}/tmp/peakAnno.Rdata")
sessionInfo()


    """
}




workflow images {
  main:
    get_images()
}

workflow {
  chipseeker_R()
}

