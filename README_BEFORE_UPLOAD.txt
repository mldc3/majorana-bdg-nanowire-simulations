MAJORANA GITHUB FLAT-UPLOAD PACKAGE

This package is designed for a browser workflow where files are uploaded loose and folders are not preserved.

1. Run the safe-delete prompt from PROMPTS_DO_NOT_UPLOAD.
2. Upload every loose file from 01_UPLOAD_BATCH_CORE_FLAT to the repository root.
3. Upload every loose file from 02_UPLOAD_BATCH_PARAMETERS_01_10_FLAT to the same root.
4. Upload every loose file from 03_UPLOAD_BATCH_PARAMETERS_11_20_FLAT to the same root.
5. Do not upload the PROMPTS_DO_NOT_UPLOAD directory.
6. After all three batches are present, run the organization-and-validation prompt.

All temporary upload filenames are globally unique. The Copilot prompt uses 000_UPLOAD_MAPPING.json to place and rename every file exactly.
