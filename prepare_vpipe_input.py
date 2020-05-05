import os
import glob
from pathlib import Path

import pandas as pd

from tqdm import tqdm


def main(root_dir):
    df = pd.read_csv(root_dir / 'results/selected_samples.csv')
    accession_list = df['accession'].tolist()

    df_final = pd.read_csv(root_dir / 'results/final_dataframe.csv', index_col=0)
    readlen_dict = df_final.to_dict()['avg_read_length']

    name_template = 'samples_{type_}'
    dummy_date = '19700101'

    for accession in tqdm(accession_list):
        available_files = glob.glob(str(root_dir / f'data/{accession}*.fastq'))

        if len(available_files) == 1:
            # SE
            type_ = 'SE'
        elif len(available_files) == 2:
            # PE
            type_ = 'PE'
        elif len(available_files) == 3:
            # PE
            type_ = 'PE'
        else:
            raise RuntimeError(available_files)

        name = Path(name_template.format(type_=type_))
        target = name / accession / dummy_date / 'raw_data'

        target.mkdir(parents=True, exist_ok=True)

        for path in available_files:
            basename = os.path.basename(path)

            if len(available_files) == 3:
                # if there is a varying number of reads per spot,
                # we only consider PE reads
                if '_1' not in basename and '_2' not in basename:
                    print('skipping', path)
                    continue

            # make V-pipe recognize PE files
            basename = basename.replace('_1', '_R1')
            basename = basename.replace('_2', '_R2')

            # create hard link
            dest = target / basename
            # print(accession, path, dest)

            os.link(path, dest)

        with open(f'{name}.tsv', 'a') as fd:
            fd.write(f'{accession}\t{dummy_date}\t{readlen_dict[accession]}\n')


if __name__ == '__main__':
    main(Path('pipeline_run/'))
