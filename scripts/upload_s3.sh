#!/bin/bash
# Usage: ./upload.sh /path/to/local_folder

LOCAL_DIR=$1
PROFILE="default"
BUCKET_BASE="s3://omics-data-002313225286/clients/gtria/prod/inputs"

# 1. Error handling: check if user provided a folder
if [ -z "$LOCAL_DIR" ]; then
    echo "Error: Please provide the folder path. Example: ./upload.sh ./Batch_01"
    exit 1
fi

echo "Checking AWS session..."
if ! aws sts get-caller-identity --profile $PROFILE > /dev/null 2>&1; then
    echo "Session expired. Opening browser for SSO login..."
    aws sso login --profile $PROFILE
fi

# 2. Extract the folder name (This becomes your BATCH_NAME)
BATCH_NAME=$(basename "$LOCAL_DIR")

echo "--- Preparing to upload Batch: $BATCH_NAME ---"

# 3. Sync to the specific Batch Folder (Archive)
aws s3 sync "$LOCAL_DIR" "$BUCKET_BASE/$BATCH_NAME/" \
    --profile $PROFILE \
    --exclude "*" \
    --include "*.bam"

# 4. Sync to 'latest' (The Trigger for the EC2)
# We still use a 'latest' folder so the EC2 knows what to download TODAY
echo "Updating the 'latest' workspace for the EC2..."
aws s3 sync "$BUCKET_BASE/$BATCH_NAME/" "s3://omics-data-002313225286/latest/" \
    --profile $PROFILE \
    --delete

echo "--- Done! Batch $BATCH_NAME is now on S3. ---"
