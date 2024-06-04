# main.py
import argparse
import os
import csv
from lib.db_utils import *
import pkg_resources


# Utilisation de pkg_resources pour obtenir les chemins des fichiers de données
def get_resource_path(relative_path):
    return pkg_resources.resource_filename(__name__, relative_path)


# Chemins des fichiers de données
DATABASE_PATH = get_resource_path('data/pankegg.db')
PATHWAY_FILE = get_resource_path('data/kegg_map_orthologs.tsv')
KO_FILE = get_resource_path('data/ko.txt')
CSV_FILE_PATH = get_resource_path('data/sample_path.csv')

TABLES_TO_DROP = [
    'taxonomy', 'bin', 'map', 'kegg', 'bin_map_kegg',
    'bin_map', 'map_kegg', 'bin_extra', 'bin_extra_kegg',
    'sample'
]


def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Process KEGG data and store it in a database.',
        epilog='Format du fichier CSV attendu: Sample name, Annotation_dir, classification_dir, Checkm2_dir'
    )
    parser.add_argument('-i', '--input', required=True, help='Path to the input CSV file.')
    parser.add_argument('-o', '--output', default='pankegg', help='Name for the output database file (without extension).')
    parser.add_argument('--output_dir', default='./db_output', help='Directory path for the output database file.')
    return parser.parse_args()


def main():
    print("Starting...")
    args = parse_arguments()
    csv_file_path = args.input
    database_name = args.output
    directory_output = args.output_dir

    if not directory_output.endswith('/'):
        directory_output += '/'
    if not os.path.exists(directory_output):
        os.makedirs(directory_output)

    database_path = os.path.join(directory_output, f'{database_name}.db')

    conn = connect_db(database_path)
    cur = conn.cursor()

    drop_all_tables(cur, TABLES_TO_DROP)
    create_tables(cur)
    ko_translate_table = load_ko_descriptions(KO_FILE)
    map_translate_table = load_pathways(PATHWAY_FILE)

    with open(csv_file_path, mode='r', newline='') as csvfile:
        csvreader = csv.DictReader(csvfile)

        for row in csvreader:
            sample = row['Sample name']
            annotation_path = row['Annotation_dir']
            sourmash_txt_dir = row['classification_dir']
            checkm2_dir = row['Checkm2_dir']

            print("Processing : ", sample)
            sample_id = add_sample_in_db(cur, sample)
            sample_preliminary_process(cur, annotation_path, sample_id)
            process_taxonomy_files(cur, sourmash_txt_dir)
            add_bin_quality_values(cur, checkm2_dir)
            map_ko_dict = process_annotation_file(cur, annotation_path, map_translate_table, sample_id)
            set_ko(cur, ko_translate_table, map_ko_dict)
            link_maps_to_kos(cur, map_ko_dict, map_translate_table)
            set_full_annot_table(cur, annotation_path, sample_id)
            link_bin_to_map_keg(cur, annotation_path, sample_id)

    conn.commit()
    cur.close()
    conn.close()

    print("Programme successfully completed.")
    print("Database: ", database_name, " Created at directory: ", database_path)


if __name__ == "__main__":
    main()
