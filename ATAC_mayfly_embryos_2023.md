# ATACseq analysis in embryos 2023
# Sample download
We download the samples raw data from the page provided by the sequencer
Alignment to the genome and peakcalling are done using a pipeline available at https://github.com/apposada/ptychodera_cisreg_development/blob/main/code/code_markdowns/04_chromatin/atac_mapping.md
First the files need to be renamed 
# Renaming the files
Some samples are in two different fastq.gz files and need to be merged
```sh
#Rename files in a new directory

for files in /data/mayfly/ATAC_mayfly/ATAC_2023/ALMUDIISA_02_2/e*.fastq.gz; do
 file=$(basename $files)
 ext=$(echo $file | cut -d "." -f 2-)
 filename=$(echo $file | cut -d "." -f -1)
 new_file=$(echo $filename | awk -F '_' -v OFS='_' '{print $1, $2, $7, $6}')
 cp "$files" proves/"$new_file"."$ext"
 echo "$filename is now $new_file"
done

#Merge the data that comes from the different reads, in one sigle fastq.gz 

for file in proves/*.fastq.gz; do
    filename=$(basename "$file")
    
    zcat "$file" "raw_data/$filename" | gzip > "merge/$filename"
done
```
# Quality control
The fasta files for each read for each sample come with a link that contains the quality control results for each file
We merge all the quality control in a multiqc
# Alignment to genome
```sh
mkdir bt_aligment

for files in ./raw_data/*R1*; do
    
    #get filename
    file=$(basename "$files")
    filename=$(echo "$file" |  cut -d "." -f 1 | cut -d "_" -f -3)

    #define the two pair-end reads
    fq1="$filename"_R1.fastq.gz
    fq2="$filename"_R2.fastq.gz  

    echo "Treating sample ${filename} ..." 
 
     perl ATAC_pipe.pl \
        -f1 $fq1 -f2 $fq2 \
        -o bt_aligment/$"filename" \
        -s /data/mayfly/genome/genome_CADEPI01_scaffoldslen.tbl \
        -i cdip_CADEPI01 \
        -p 12 -bp /data/mayfly/genome/bowtie2_index/ \
        -tmp tmp/ \
        -ov_th \
    2> bt_aligment/"$filename"_ATAC_pipe.log 1>&2
done 
```

# With the aligned files lets create the Consensus peaks using the MACS2, IDR, and bedtools
The reason why we do this is to ensure that the peaks are free of backgroung noise
```sh
for files in bt_aligment/*_1_*nucfree.bed; do

    #get filename
    file=$(basename "$files")
    filename=$(echo "$file" |  cut -d "." -f 1 | cut -d "_" -f 1)

    rep1="$filename"_1_001_nucfree
    rep2="$filename"_2_001_nucfree
   
    mkdir ./peaks/"$filename"
    
    bash idr_mod_job_customc.sh \
        $filename \
        180287000 \
        $rep1 \
        $rep2 \
        0.05

done

#idr_mod_job_customc.sh is the name of the shell script with the code to run everything
#the $filename, 18..., rep1, rep2 and 0.05 correspond to the name of the file, GENOME_size, Nucfree1, Nucfree2 and pvalue_macs2 necessaries to run the idr script
```
Now, we make an unique consensus peak file, from all the IDR analysis, after this we generate a Counts files that will be used DESeq2 to analyse differential chromatin accessibility.

```sh
#generating the consensus peak file

    cat *ConsPeaks.bed | \
    sort -k1,1 -k2,2n | \
    cut -f1,2,3,4,5,6 | \
    mergeBed | sortBed | \
    awk 'BEGIN {OFS="\t"} {print $1,$2,$3,"peak"NR}' \
    > pfla_all_peaks.bed

#Generating the Counts files

for i in ../bt_aligment/*_nucfree.bed ; do
  
  x=${i##*/}
  
  z=${x%_nucfree.bed}
  
  echo Starting intersect ... sample $x
  
  mkdir ${z}
  
  intersectBed -c -a pfla_all_peaks.bed -b $i -nonamecheck > ${z}/${z}_counts.txt  #intersectBed -c -a pfla_all_peaks.bed -b ../bt_aligment/*e10_1_001_nucfree.bed  -nonamecheck > e10_1_001_prova_counts.txt
  
  echo done. Starting parsing to bed... sample $x
  
  awk 'BEGIN {OFS="\t"} {print $1,$2,$3,"peak"NR,$5,"+"}' ${z}/${z}_counts.txt > ${z}/${z}.bed  #awk 'BEGIN {OFS="\t"} {print $1,$2,$3,"peak"NR,$5,"+"}' e10_1_001_prova_counts.txt > e10_1_001_prova_counts.bed
  
  echo done. Starting parsing to counts... sample $x
  
  cut -f 4,5 ${z}/${z}.bed > ${z}/${z}.counts   

  rm ${z}/${z}_counts.txt
  
  echo done with sample $x .

done

echo Done. 
  
```
bedToBigBed to visualize the peaks in IGV
```sh
for files in ./bedtools/renamed_embryo_idrConsPeaks/e*_idrConsPeaks.bed; do
    #get filename
    file=$(basename "$files")
    echo "$file"

    sort -k1,1 -k2,2n ./bedtools/renamed_embryo_idrConsPeaks/"$file" > temp_"$file"_sorted.bed
    echo "temp_"$file"_sorted.bed"
    echo "starting BB creation"
    bedToBigBed temp_"$file"_sorted.bed /data/mayfly/genome/genome_CADEPI01_scaffoldslen.tbl > "$file"_myBigBed.bb
    #rm temp_"$file"_sorted.bed
done    

```
# Classification of the peaks in regulatory regions ("zones")
We got the peaks free of background noise using IDR, now we classify them in six different categories: 
1- Promotor zone (Transcription start site -1000bases + 500bases)
2- Proximal zone (up to 4000bases away from the promotor zone)
3- Gene body zone (from the end of the promotor zone to the end of the gene)
4- Distal zone (all the peaks that cant be classified in the previous categories)
5- Basal zone (Transcription start site -1000bases + 5000bases if there is another TTS, the zone stops there)
6- Great zone (1megabase away from TTS upstream and downstream if there is a basal zone, the zone stops there)
In order to do these zones files we use bedtools to properly sort and make sure that the selected regions are inside the scaffold lenght.
We will use bedtools slop, closest and bedtools intersect

Creating all the regions("zones") to classify the peaks, the python script used is available at: https://github.com/m-rossello/GeneRegLocator
```sh
# extractig the TSS bed file. $3 is TSS+1bp to ensure the bedtools can run properly
awk '{print $1, $2, $2+1, $4, $5, $6}' clodip_v2_annotation_onlygenes_corrected.bed > TSS.bed
# promotor_zone, instead of 500bp we substract 499 to compensate the extra bp added when extracting the TSS bed file
# -a imput file
# -s number of bases to add on the start coordinate (to substract put a negative value)
# -e number of bases to add on the end coordinate (to substract put a negative value)
# -f (optional) used to filter 
# -o ouput
python3 /home/mayfly/maria_r/my_scripts/define_reg_zones.py -a TSS.bed -s -1000 -e 499 -l /data/mayfly/genome/genome_CADEPI01_scaffoldslen.tbl -o promoter_zone.bed

# proximal_zone, 1001 is substracted from the end coord to compensate the extra bp from the TSS
python3 /home/mayfly/maria_r/my_scripts/define_reg_zones.py -a TSS.bed -s -5000 -e -1001 -l /data/mayfly/genome/genome_CADEPI01_scaffoldslen.tbl -o zones/promoter.bed -o proximal_zone.bed

# genebody_zone, we first mmerge promotor and proximal zone to filter out possible overlaps
cat promoter_zone.bed proximal_zone.bed | sort -k1,1 -k2,2n | groupBy -g 1,4 -c 4,2,3 -o count,min,max | awk -v OFS='\t' '{print $1, $4, $5, $2, ".", "."}' > promoter_and_proximal.bed

python3 /home/mayfly/maria_r/my_scripts/define_reg_zones.py -a clodip_v2_annotation_onlygenes_corrected.bed -s 500 -e 0 -l /data/mayfly/genome/genome_CADEPI01_scaffoldslen.tbl -f promoter_and_proximal.bed -o genebody_zone.bed

# basal_zone, 4999 is added in the end coord to compensate the extra bp from the TSS
python3 /home/mayfly/maria_r/my_scripts/define_reg_zones.py -a TSS.bed -s -1000 -e 4999 -l /data/mayfly/genome/genome_CADEPI01_scaffoldslen.tbl -o basal_zone.bed

#great zone, 999.999 is added in the end coord to compensate the extra bp from the TSS
python3 /home/mayfly/maria_r/my_scripts/define_reg_zones.py -a TSS.bed -s -1000000 -e 999999 -l /data/mayfly/genome/genome_CADEPI01_scaffoldslen.tbl -f basal_zone.bed -o great_zone.bed

```
Peak classification
```sh
# with the defined zones, lets classify each peak into ints corresponding region
for files in ../counts/e*_idrConsPeaks.bed; do
    #get filename
    file=$(basename "$files")
    filename=$(echo "$file" |  cut -d "." -f 1)
    name=$(echo "$filename" | cut -d "_" -f 1);
    echo "$name"

    #sort files, create a temporary file for each
    sort -k1,1 -k2,2n -k3,3n ../counts/"$filename".bed > sorted_${filename}_temp.bed

    bedtools intersect \
        -a  promoter_zone.bed \
        -b  sorted_${filename}_temp.bed  \
        -wb \
        -F 0.501 \
        > peaks_classification/"$name"_promoter_zone.bed
    echo "promoter zone created"

    bedtools subtract -a ../counts/"$filename".bed -b "peaks_classification/"$name"_promoter_zone.bed" -A > temp_${name}_all_but_promoter.bed

    bedtools intersect \
        -a  proximal_zone.bed \
        -b  temp_${name}_all_but_promoter.bed  \
        -wb \
        -F 0.501 \
        > peaks_classification/"$name"_proximal_zone.bed
    echo "proximal zone created"    

    bedtools subtract -a temp_${name}_all_but_promoter.bed -b "peaks_classification/${name}_proximal_zone.bed" -A > temp_${name}_all_but_proximal.bed      
    
    bedtools intersect \
        -a  genebody_zone.bed \
        -b  temp_${name}_all_but_proximal.bed  \
        -wb \
        -F 0.501 \
        > peaks_classification/"$name"_genebody_zone.bed
    echo "gemebody zone created"      

     bedtools subtract -a temp_${name}_all_but_proximal.bed -b "peaks_classification/${name}_genebody_zone.bed" -A > ${name}_distal_zone.bed

    rm temp_${name}_all_but_proximal.bed  
    rm temp_${name}_all_but_promoter.bed
    sorted_${filename}_temp.bed

done

```
We make sure that the files are different, and dont have similar lines (OPTIONAL, making sure that the distal_zone is done properly)
```sh
compare_bed_files() {
    file1=$1
    file2=$2

    if [ ! -f "$file1" ] || [ ! -f "$file2" ]; then
        echo "Error: File not found."
    fi

    if [ "$(wc -l < "$file1")" -ne "$(wc -l < "$file2")" ]; then
        echo "The BED files are different."
    fi

    if ! cmp -s "$file1" "$file2"; then
        echo "The BED files are different."
    fi

    echo "The BED files are identical."
}
# Usage
file1=distal/"e10_distal_zone.bed"
file2="e10_proximal_zone.bed"
compare_bed_files "$file1" "$file2"


# no similar lines
awk 'FNR==NR{arr[$4]=1; next} $4 in arr' distal/"e10_distal_zone.bed" "e10_proximal_zone.bed"
```
# Peaks to genes association
```sh
#the previously defined zones are used to associate peaks.
#Othologous pairs file, comes from an orthology association with broccoli
bedtools intersect -a zones/great_zone.bed -b counts/pfla_all_peaks.bed -wa -wb -F 0.5 | awk  -v FS='\t' -v OFS='\t' '{print $10, $4}' > peaks_to_genes.tbl

awk -v FS='\t' -v OFS='\t' 'FNR==NR {a[$2]=$0; next} $1 in a {print a[$1], $0}' peaks_to_genes.tbl orthogrups/orthologous_pairs_cdip_annotated.tbl | awk -v FS='\t' -v OFS='\t' '{print $1,$2,$4,$5,$6,$7}' > peaks_to_genes_annotation.tbl
```

# Merge counts files for each sample in a table, and extract peak width to normalize
Lets start with merging all the counts files into a big file, before normalizing by peak width and library size

```python
#import the following packages
import pandas as pd
import glob
import os

#Get the counts files for each sample
#Getting all samples counts files, it has to match the pattern "e*.counts"
embryo_counts_files = 'C:/Users/usuario/Desktop/TFM/ATAC-seq/results/counts/counts_files/e*.counts'
delimiter = '\t'  

#find file paths matching the pattern
embryo_sample_file_paths = glob.glob(embryo_counts_files)

#print to see if its getting all samples
print (embryo_sample_file_paths)


if not embryo_sample_file_paths:
    print("No files found matching the pattern:", embryo_counts_files)
else:
    # Read files and store them in separate data frames
    dfs = []
    for embryo_sample_file_path in embryo_sample_file_paths:
        file_name = os.path.basename(embryo_sample_file_path)
        column_name = os.path.splitext(file_name)[0]
        df = pd.read_csv(embryo_sample_file_path, delimiter="\t", header=None, names=["Peaks", column_name])
        dfs.append(df)
        print("Read file:", embryo_sample_file_path)

    # Merge data frames based on the common column
    merged_df = dfs[0]  # Initialize with the first data frame
    for df in dfs[1:]:
        merged_df = pd.merge(merged_df, df, on="Peaks")

    # Display the merged data frame
    print(merged_df)

    # Save merged output as .counts file with tab-separated values
    output_file_path = 'C:/Users/usuario/Desktop/TFM/ATAC-seq/results/counts/counts_files/merged_embryos.counts'
    merged_df.to_csv(output_file_path, sep="\t", index=False)

    #Print confirmation message 
    print("Merged output saved to:", output_file_path)

#The merged_counts files gives us a table with all the counts for each sample (each sample is a column), this can be used to normalize the data
```

Get peak width to normalize
```sh
#Just print the width of each peak, this will be used to normalize the data
awk '{print $4,$3-$2}' atac_embryo_2023/counts/pfla_all_peaks.bed > pfla_all_peaks_width.txt
```
Now we move to Rstudio to normalize the data with limma analysis and proceed with analysis.