# 1st batch ovary snRNAseq
cd ~/ovary_SC

cellranger count --id=ov-cor-ct \
                   --transcriptome=/home/buckcenter.org/fwu/cellranger/refdata-gex-GRCh38-2020-A \
                   --fastqs=/bigrock/FurmanLab/DataRepository/U54_Ovary_demultipx/ov-cor-ct \
                   --localcores=32

cellranger count --id=ov-cor-dox \
                   --transcriptome=/home/buckcenter.org/fwu/cellranger/refdata-gex-GRCh38-2020-A \
                   --fastqs=/bigrock/FurmanLab/DataRepository/U54_Ovary_demultipx/ov-cor-dox \
                   --localcores=32

cellranger count --id=ov-med-ct \
                   --transcriptome=/home/buckcenter.org/fwu/cellranger/refdata-gex-GRCh38-2020-A \
                   --fastqs=/bigrock/FurmanLab/DataRepository/U54_Ovary_demultipx/ov-med-ct \
                   --localcores=16

cellranger count --id=ov-med-dox \
                   --transcriptome=/home/buckcenter.org/fwu/cellranger/refdata-gex-GRCh38-2020-A \
                   --fastqs=/bigrock/FurmanLab/DataRepository/U54_Ovary_demultipx/ov-med-dox \
                   --localcores=16

# 2nd batch ovary snRNAseq
mkdir /data/array2/fwu/ovary_SC_2nd
cd /data/array2/fwu/ovary_SC_2nd

cellranger count --id=RTL_856_ov-cor-ct-1 \
                   --transcriptome=/home/buckcenter.org/fwu/cellranger/refdata-gex-GRCh38-2020-A \
                   --fastqs=/bigrock/FurmanLab/DataRepository/U54/U54_Data/Bone_10X/BI012-10X/HCCCLDRX3_1/outs/fastq_path/HCCCLDRX3_1 \
                   --sample=RTL_856_CTX_CTRL-1 \
                   --localcores=32

cellranger count --id=RTL_856_ov-cor-dox-1 \
                   --transcriptome=/home/buckcenter.org/fwu/cellranger/refdata-gex-GRCh38-2020-A \
                   --fastqs=/bigrock/FurmanLab/DataRepository/U54/U54_Data/Bone_10X/BI012-10X/HCCCLDRX3_1/outs/fastq_path/HCCCLDRX3_1 \
                   --sample=RTL_856_CTX_DOXO-1 \
                   --localcores=32

cellranger count --id=RTL_856_ov-cor-dox-2 \
                   --transcriptome=/home/buckcenter.org/fwu/cellranger/refdata-gex-GRCh38-2020-A \
                   --fastqs=/bigrock/FurmanLab/DataRepository/U54/U54_Data/Bone_10X/BI012-10X/HCCCLDRX3_1/outs/fastq_path/HCCCLDRX3_1 \
                   --sample=RTL_856_Cortex_DOXO-2 \
                   --localcores=32

cellranger count --id=RTL_856_ov-med-ct-1 \
                   --transcriptome=/home/buckcenter.org/fwu/cellranger/refdata-gex-GRCh38-2020-A \
                   --fastqs=/bigrock/FurmanLab/DataRepository/U54/U54_Data/Bone_10X/BI012-10X/HCCCLDRX3_1/outs/fastq_path/HCCCLDRX3_1 \
                   --sample=RTL_856_Medulla_CTRL-1 \
                   --localcores=32

cellranger count --id=RTL_856_ov-med-dox-1 \
                   --transcriptome=/home/buckcenter.org/fwu/cellranger/refdata-gex-GRCh38-2020-A \
                   --fastqs=/bigrock/FurmanLab/DataRepository/U54/U54_Data/Bone_10X/BI012-10X/HCCCLDRX3_1/outs/fastq_path/HCCCLDRX3_1 \
                   --sample=RTL_856_Medulla_DOXO-1 \
                   --localcores=32

cellranger count --id=RTL_856_ov-med-dox-2 \
                   --transcriptome=/home/buckcenter.org/fwu/cellranger/refdata-gex-GRCh38-2020-A \
                   --fastqs=/bigrock/FurmanLab/DataRepository/U54/U54_Data/Bone_10X/BI012-10X/HCCCLDRX3_1/outs/fastq_path/HCCCLDRX3_1 \
                   --sample=RTL_856_Medulla_DOXO-2 \
                   --localcores=32

cellranger count --id=RTL_857_ov-cor-ct-1 \
                   --transcriptome=/home/buckcenter.org/fwu/cellranger/refdata-gex-GRCh38-2020-A \
                   --fastqs=/bigrock/FurmanLab/DataRepository/U54/U54_Data/Bone_10X/BI012-10X/HCCCLDRX3_1/outs/fastq_path/HCCCLDRX3_1 \
                   --sample=RTL_857_CTX_CTRL-1 \
                   --localcores=32

cellranger count --id=RTL_857_ov-cor-dox-1 \
                   --transcriptome=/home/buckcenter.org/fwu/cellranger/refdata-gex-GRCh38-2020-A \
                   --fastqs=/bigrock/FurmanLab/DataRepository/U54/U54_Data/Bone_10X/BI012-10X/HCCCLDRX3_1/outs/fastq_path/HCCCLDRX3_1 \
                   --sample=RTL_857_CTX_DOXO-1 \
                   --localcores=32