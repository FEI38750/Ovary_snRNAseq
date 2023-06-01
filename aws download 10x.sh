# Confirm that the S3 location is accessible using the ls command
aws s3 ls --profile wasabi --endpoint-url=https://s3.wasabisys.com s3://seqmatic-data-releases/tommy_tran/
# Transfer data to local directory using the sync command
aws s3 sync --profile wasabi --endpoint-url=https://s3.wasabisys.com s3://seqmatic-data-releases/tommy_tran/230328_BI012_Rapa_Bone_SEQ23-SQ01188_10X /bigrock/FurmanLab/DataRepository/U54/U54_Data/Bone_10X