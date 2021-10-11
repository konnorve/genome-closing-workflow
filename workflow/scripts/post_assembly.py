import json
from pathlib import Path
import pandas as pd

genome_report_paths = snakemake.input.sample_reports # Path('/Users/kve/Desktop/Chisholm_lab_non-DB/genome_closing')

strain_ID_df = pd.read_csv(Path(snakemake.input.s_ids_tsv), sep="\t")

strain_ID_df['forward barcode'] = strain_ID_df['reverse barcode'] = strain_ID_df['Barcode']

genome_df_list = []

for report_path_str in genome_report_paths:
    report_path = Path(report_path_str)

    with open(report_path) as f:
        data = json.load(f)

    columns_dict_list = data.get('tables')[0].get('columns')

    barcode_df = pd.DataFrame()

    barcode_str = report_path.parent.parent.name
    fwd_bc, rev_bc = (barcode_str, barcode_str)

    for c, vals in [(c.get('header'), c.get('values')) for c in columns_dict_list]:
        barcode_df[c] = vals

    sort_order = ['forward barcode', 'reverse barcode'] + barcode_df.columns.tolist()

    barcode_df['forward barcode'] = fwd_bc
    barcode_df['reverse barcode'] = rev_bc

    barcode_df = barcode_df[sort_order]

    genome_df_list.append(barcode_df)

genome_df = pd.concat(genome_df_list)

genome_df = strain_ID_df.merge(genome_df, how='right',on=['forward barcode', 'reverse barcode'])

df_out_path = Path(snakemake.output.dataframe_path)

df_out_path.parent.mkdir(parents=True, exist_ok=True)

genome_df.to_csv(df_out_path, sep='\t')
