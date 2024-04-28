import os
import glob

dir_name = 'fastq_files'
for filename in glob.glob(f'{dir_name}/*'):
    sra_code = filename.split('/')[1].split('_')[0]
    suffix = filename.split('_')[-1]
    new_name = f'{dir_name}/{sra_code}_{suffix}'
    os.rename(filename, new_name)
    print(new_name)
