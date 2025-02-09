# db_utils.py

import sqlite3
import csv
import glob
import os
from .sql_commands import *


def connect_db(path):
    return sqlite3.connect(path)


def create_tables(cur):
    for table, sql in CREATE_TABLES.items():
        cur.execute(sql)
    # conn.commit()


def drop_all_tables(cur, tables):
    for table in tables:
        cur.execute(DROP_TABLE.format(table_name=table))
    # conn.commit()


def drop_table(conn, table_name):
    cur = conn.cursor()
    cur.execute(DROP_TABLE.format(table_name=table_name))
    conn.commit()


def insert_taxonomy(cur, tax_data):
    # Vérifier si l'entrée existe déjà
    cur.execute(SELECT_TAXONOMY_ID, tax_data[2:9])
    result = cur.fetchone()
    if result:
        return result[0]
    else:
        # Insérer la nouvelle entrée taxonomique
        cur.execute(INSERT_TAXONOMY, tax_data[2:9])
        return cur.lastrowid


def insert_bin(cur, bin_name, taxonomic_id, sample_id):
    cur.execute(SELECT_BIN_ID, (bin_name, sample_id))
    bin_id = cur.fetchone()
    if not bin_id:  # Insert only if the bin does not exist for this sample
        cur.execute(INSERT_BIN, (bin_name, taxonomic_id, sample_id))

def load_pathways(pathway_file):
    """Charge les données des pathways depuis un fichier et les stocke dans un dictionnaire."""
    pathways_dict = {}
    with open(pathway_file, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            map_id = row[0].strip()
            pathway_description = row[1].strip() if len(row) > 1 else None
            # pathways_orthologs_number = int(row[2].strip())
            pathways_orthologs_number = row[2].split(",")
            #  print(pathways_orthologs_number)
            #  pathways_dict[map_id] = pathway_description
            pathways_dict[map_id] = [pathway_description, pathways_orthologs_number]
    return pathways_dict


def load_ko_descriptions(ko_file):
    """Charge les descriptions des KO depuis un fichier et les stocke dans un dictionnaire."""
    ko_dict = {}
    with open(ko_file, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            ko_id = row[0].strip()
            description = row[1].split(";") if len(row) > 1 else ["Unknown"]
            ko_dict[ko_id] = description
    return ko_dict


# def process_pathways_and_kos(cur, pathway, maps_dict, pathways_dict, kos, bin_map_list):
#     if pathway not in bin_map_list:
#         bin_map_list.append(pathway)
#
#     if pathway not in maps_dict:
#         # Tente d'insérer le pathway dans la base de données si ce n'est pas déjà fait
#         description = pathways_dict.get(pathway, "Description not available")  # Utilise une description par défaut si non trouvée
#         cur.execute(INSERT_MAP, (description, pathway, pathway))
#         maps_dict[pathway] = kos
#     else:
#         map_ko_list = maps_dict.get(pathway)
#         original_length = len(map_ko_list)
#         map_ko_list.extend([ko for ko in kos if ko not in map_ko_list])
#         if len(map_ko_list) > original_length:
#             maps_dict[pathway] = map_ko_list


def process_pathways_and_kos(cur, pathway, maps_dict, pathways_dict, kos, bin_map_list):
    """
    Process the pathways and associated KEGG Orthology identifiers (KOs), and inserts them into the database.

    Args:
        cur (sqlite3.Cursor): Database cursor for executing SQL commands.
        pathway (str): The pathway map number (e.g., "map00010").
        maps_dict (dict): Dictionary tracking the current state of pathways and associated KOs in the database.
        pathways_dict (dict): Dictionary containing pathway descriptions and number of orthologs.
        kos (list): List of KEGG Orthology identifiers associated with the pathway.
        bin_map_list (list): List tracking which pathways have been processed.
    """
    if pathway not in bin_map_list:
        bin_map_list.append(pathway)

    if pathway not in maps_dict:
        # Retrieve the pathway information from the pathways_dict, using a default if not found
        pathway_info = pathways_dict.get(pathway, ["Description not available", []])
        description, total_orthologs = pathway_info  # Extract both description and total orthologs

        # Attempt to insert the pathway into the database if it's not already done
        cur.execute(INSERT_MAP, (description, pathway, len(total_orthologs), pathway))
        maps_dict[pathway] = kos

    else:
        # If the pathway is already in maps_dict, update it with any new KOs
        map_ko_list = maps_dict.get(pathway, []).copy()
        original_length = len(map_ko_list)
        # Extend the list only with new, unique KOs
        map_ko_list.extend([ko for ko in kos if ko not in map_ko_list])
        if len(map_ko_list) > original_length:
            maps_dict[pathway] = map_ko_list


def link_bins_to_pathways(cur, bin_name, bin_map_list, sample_id):
    for map_number in bin_map_list:
        cur.execute(SELECT_BIN_ID, (bin_name, sample_id))  # Include sample_id
        bin_id = cur.fetchone()
        if bin_id:
            bin_id = bin_id[0]
            cur.execute(SELECT_MAP_ID, (map_number,))
            map_id = cur.fetchone()
            if map_id:
                map_id = map_id[0]
                cur.execute(INSERT_BIN_MAP, (bin_id, map_id, bin_id, map_id))


def link_maps_to_kos(cur, maps_dict, map_translate_table):
    for map_name in maps_dict:
        cur.execute(SELECT_MAP_ID, (map_name,))
        map_id = cur.fetchone()
        if map_id:
            map_id = map_id[0]
            for ko in maps_dict.get(map_name):
                if not ko == "":
                    ko_entry = ko.split(":")[1]
                    map_kegg_list = map_translate_table.get(map_name, [None, []])[1]
                    # try:
                    #     map_kegg_list = map_translate_table.get(map_name)[1]
                    # except:
                    #     map_kegg_list = []
                    ko_is_in_map = 0
                    if ko_entry in map_kegg_list:
                        ko_is_in_map = 1

                    cur.execute(SELECT_KEGG_ID, (ko_entry,))
                    ko_id = cur.fetchone()
                    if ko_id:
                        ko_id = ko_id[0]
                        # Insérer la jointure dans bin_map si elle n'existe pas déjà
                        cur.execute(INSERT_MAP_KEGG, (map_id, ko_id, ko_is_in_map, map_id, ko_id, ko_is_in_map))


def add_bin_quality_values(cur, quality_report_path, sample_id):  # Include sample_id
    with open(quality_report_path, 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            bin_name = row['Name']
            if bin_name.endswith(".fa"):
                bin_name = bin_name[:-3]
            completeness = float(row['Completeness'])
            contamination = float(row['Contamination'])

            # Update with sample_id to ensure correct association
            cur.execute(UPDATE_BIN_QUALITY, (completeness, contamination, bin_name, sample_id))


def set_ko(cur, kegg_translate_table, maps_dict):
    unique_kos = set()
    for kos in maps_dict.values():
        unique_kos.update(kos)
    unique_ko_list = list(unique_kos)
    for ko in unique_ko_list:
        if not ko == "":
            ko_id = ko.split(":")[1]
            kegg_info = kegg_translate_table.get(ko_id)
            if kegg_info:
                kegg_name = kegg_info[0]
                kegg_full_name = kegg_info[1]
                cur.execute(INSERT_KEGG, (ko_id, kegg_name, kegg_full_name, ko_id))

            else:
                cur.execute(INSERT_KEGG, (ko_id, "null", "null", ko_id))


def process_taxonomy_files(cur, file_path, sample_id):
    files = glob.glob(file_path)
    for file_path in files:
        with open(file_path, 'r') as file:
            next(file)  # Skip header
            line = next(file).strip()
            data = line.split(',')
            assigned_status = data[1]
            bin_name = data[0]
            if bin_name.endswith(".fa"):
                bin_name = bin_name[:-3]
            if assigned_status != "nomatch":
                tax_id = insert_taxonomy(cur, data)
                cur.execute(UPDATE_BIN_TAXONOMY, (tax_id, bin_name, sample_id))


def sample_preliminary_process(cur, annotation_files_path, sample_id):
    for filename in glob.glob(annotation_files_path):
        base_name = os.path.basename(filename)
        name_part = base_name.split('.')[0] + '.' + base_name.split('.')[1]

        if name_part.endswith(".fa"):
            name_part = name_part[:-3]

        cur.execute(SELECT_BIN_ID, (name_part, sample_id))  # Ensure unique bins per sample
        bin_id = cur.fetchone()
        if not bin_id:
            insert_bin(cur, name_part, "null", sample_id)


def process_annotation_file(cur, annotation_files_path, pathways_dict, sample_id):
    maps_dict = {}

    for filename in glob.glob(annotation_files_path):
        with open(filename, 'r', newline='', encoding='utf-8') as file:
            base_name = os.path.basename(filename)
            name_part = base_name.split('.')[0] + '.' + base_name.split('.')[1]

            # Read the file, skipping the first line (##)
            lines = file.readlines()
            lines = lines[1:]  # Skip first line (##) to use the correct header

            # Create a DictReader with the corrected header
            reader = csv.DictReader(lines, delimiter='\t')

            # Print detected headers for debugging
            #print(f"Processing file: {filename}")
            #print(f"Detected headers: {reader.fieldnames}")

            # Ensure 'KEGG_Pathway' exists
            if 'KEGG_Pathway' not in reader.fieldnames:
                print(f"Skipping {filename}: 'KEGG_Pathway' column not found!")
                continue

            bin_map_list = []
            for row in reader:
                #print(f"Row data: {row}")  # Debugging
                kegg_pathways = row.get('KEGG_Pathway', '')  # Use get() to avoid KeyError
                kegg_kos = row.get('KEGG_ko', '')

                kos = kegg_kos.split(',')
                if kegg_pathways:
                    for pathway in kegg_pathways.split(','):
                        pathway = pathway.strip()
                        if pathway.startswith('map'):
                            process_pathways_and_kos(cur, pathway, maps_dict, pathways_dict, kos, bin_map_list)
                else:
                    pathway = name_part + "_not_mapped"
                    process_pathways_and_kos(cur, pathway, maps_dict, pathways_dict, kos, bin_map_list)

            bin_name = name_part
            if bin_name.endswith(".fa"):
                bin_name = bin_name[:-3]
            link_bins_to_pathways(cur, bin_name, bin_map_list, sample_id)

    return maps_dict


def link_full_line_with_kos(cur, bin_id, kegg_gos, kegg_kos, kegg_free_desc):
    # Insert the full annotation line
    cur.execute(INSERT_BIN_EXTRA, (bin_id, kegg_gos, kegg_kos, kegg_free_desc))
    kegg_extra_line_id = cur.lastrowid

    if kegg_extra_line_id and kegg_kos:  # Ensure kegg_kos is not None or empty
        kos = kegg_kos.split(',')
        for ko in kos:
            if ":" not in ko:  # Skip invalid KO entries
                #print(f"Skipping invalid KEGG_ko entry: {ko}")
                continue

            ko_entry = ko.split(":")[1]  # Extract KO ID safely

            # Check if KO entry exists in the database
            cur.execute(SELECT_KEGG_ID, (ko_entry,))
            ko_id = cur.fetchone()

            if ko_id:
                ko_id = ko_id[0]
                cur.execute(INSERT_BIN_EXTRA_KEGG, (kegg_extra_line_id, ko_id, kegg_extra_line_id, ko_id))
            else:
                #print(f"Warning: KEGG ID {ko_entry} not found in the database")
                continue


def set_full_annot_table(cur, annotation_files_path, sample_id):
    for filename in glob.glob(annotation_files_path):
        with open(filename, 'r', newline='', encoding='utf-8') as file:
            base = os.path.basename(filename)
            name_part = base.split('.')[0] + '.' + base.split('.')[1]
            bin_name = name_part
            if bin_name.endswith(".fa"):
                bin_name = bin_name[:-3]

            cur.execute(SELECT_BIN_ID, (bin_name, sample_id))
            bin_id = cur.fetchone()
            if bin_id:
                bin_id = bin_id[0]

                # Read file and skip first line
                lines = file.readlines()
                lines = lines[1:]  # Skip first line (##)

                reader = csv.DictReader(lines, delimiter='\t')

                # Debugging: Print the actual column names detected
                #print(f"Processing file: {filename}")
                #print(f"Detected headers: {reader.fieldnames}")

                if 'KEGG_ko' not in reader.fieldnames:
                    print(f"Skipping {filename}: 'KEGG_ko' column not found!")
                    continue  # Skip this file instead of crashing

                for row in reader:
                    #print(f"Row data: {row}")  # Debugging

                    # Use get() to avoid crashing
                    kegg_kos = row.get('KEGG_ko', None)
                    kegg_gos = row.get('GOs', None)
                    kegg_free_desc = row.get('eggNOG free text desc.', None)

                    if kegg_kos:
                        link_full_line_with_kos(cur, bin_id, kegg_gos, kegg_kos, kegg_free_desc)


# def link_bin_to_map_keg(cur, annotation_files_path):
#     for filename in glob.glob(annotation_files_path):
#         with open(filename, 'r', newline='') as file:
#             base_name = os.path.basename(filename)
#             name_part = base_name.split('.')[0] + '.' + base_name.split('.')[1]
#             bin_name = name_part + ".fa"
#
#             cur.execute(SELECT_BIN_ID, (bin_name,))
#             bin_id = cur.fetchone()[0]  # Supposons qu'il y a toujours un bin_id valide
#
#             reader = csv.DictReader(file, delimiter='\t')
#             for row in reader:
#                 kegg_pathways = row['KEGG_Pathway']
#                 kegg_kos = row['KEGG_ko']
#                 kos = kegg_kos.split(',')
#
#                 for pathway in kegg_pathways.split(',') if kegg_pathways else [name_part + "_not_mapped"]:
#                     pathway = pathway.strip()
#                     if pathway.startswith('map'):
#                         cur.execute(SELECT_MAP_ID, (pathway,))
#                         map_id = cur.fetchone()[0]
#
#                         for ko_entry in kos:
#                             ko_id = ko_entry.split(":")[1]
#                             cur.execute(SELECT_KEGG_ID, (ko_id,))
#                             kegg_id = cur.fetchone()[0]  # Supposons que kegg_id est toujours valide
#
#                             cur.execute(SELECT_MAP_KEGG_ID, (map_id, kegg_id,))
#                             map_kegg_id = cur.fetchone()[0]  # Supposons que map_kegg_id est toujours valide
#
#                             cur.execute(INSERT_BIN_MAP_KEGG, (bin_id, map_kegg_id, bin_id, map_kegg_id))


def count_map_kegg_ids(cur):
    query = """
    SELECT map_kegg_id, COUNT(bin_id) as bin_count
    FROM bin_map_kegg
    GROUP BY map_kegg_id
    ORDER BY bin_count DESC;
    """
    cur.execute(query)
    results = cur.fetchall()
    return results


def total_shared_associations(cur):
    # Cette requête calcule la somme totale des associations partagées,
    # où chaque bin_count est le nombre de bins partageant un même map_kegg_id.
    query = """
    SELECT SUM(bin_count)
    FROM (
        SELECT COUNT(bin_id) as bin_count
        FROM bin_map_kegg
        GROUP BY map_kegg_id
    ) as counts;
    """
    cur.execute(query)
    result = cur.fetchone()  # Fetch the result which is the total sum of shared associations
    return result[0] if result else 0  # Return the sum or 0 if the result is None


def link_kegg_to_map(cur, kos, map_id, bin_id):
    for ko_entry in kos:
        if ko_entry:
            ko_id = ko_entry.split(":")[1]
            cur.execute(SELECT_KEGG_ID, (ko_id,))
            kegg_id = cur.fetchone()[0]  # Supposons que kegg_id est toujours valide

            cur.execute(SELECT_MAP_KEGG_ID, (map_id, kegg_id,))
            map_kegg_id = cur.fetchone()[0]  # Supposons que map_kegg_id est toujours valide

            cur.execute(INSERT_BIN_MAP_KEGG, (bin_id, map_kegg_id, bin_id, map_kegg_id))


def link_bin_to_map_keg(cur, annotation_files_path, sample_id):
    for filename in glob.glob(annotation_files_path):
        with open(filename, 'r', newline='', encoding='utf-8') as file:
            base_name = os.path.basename(filename)
            name_part = base_name.split('.')[0] + '.' + base_name.split('.')[1]
            bin_name = name_part
            if bin_name.endswith(".fa"):
                bin_name = bin_name[:-3]

            cur.execute(SELECT_BIN_ID, (bin_name, sample_id))
            bin_id = cur.fetchone()
            if not bin_id:
                print(f"Warning: No bin_id found for {bin_name}")
                continue
            bin_id = bin_id[0]

            # Skip first line (##) and parse as CSV
            lines = file.readlines()
            lines = lines[1:]  # Skip header comment
            reader = csv.DictReader(lines, delimiter='\t')

            print(f"Processing file: {filename}")
            print(f"Detected headers: {reader.fieldnames}")

            for row in reader:
                if 'KEGG_Pathway' not in row:
                    print(f"Skipping row due to missing KEGG_Pathway: {row}")
                    continue  # Skip row if the column is missing

                kegg_pathways = row.get('KEGG_Pathway', '')  # Use .get() to avoid KeyError
                kegg_kos = row.get('KEGG_ko', '')

                if not kegg_pathways or kegg_pathways == '-':
                    print(f"Skipping row: No valid KEGG_Pathway found")
                    continue  # Skip empty pathway rows

                kos = kegg_kos.split(',')

                for pathway in kegg_pathways.split(','):
                    pathway = pathway.strip()
                    if pathway.startswith('map'):
                        cur.execute(SELECT_MAP_ID, (pathway,))
                        map_id = cur.fetchone()
                        if not map_id:
                            print(f"Warning: No map_id found for pathway {pathway}")
                            continue
                        map_id = map_id[0]

                        link_kegg_to_map(cur, kos, map_id, bin_id)


def add_sample_in_db(cur, name):
    cur.execute(INSERT_SAMPLE, (name,))
    return cur.lastrowid
