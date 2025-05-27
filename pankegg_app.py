import os
import sqlite3
from flask import Flask, render_template, request, session, redirect, url_for, Response, jsonify, abort
import csv
import io
import json
import click
import sys
from importlib.resources import files
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist
# Import your robust DB print function
from lib.db_utils import print_db_status

# set database path
if "--d" in sys.argv or "--db" in sys.argv:
    try:
        if "--d" in sys.argv:
            idx = sys.argv.index("--d")
        else:
            idx = sys.argv.index("--db")
        db_path_arg = sys.argv[idx + 1]
        os.environ["PANKEGG_DB_PATH"] = db_path_arg
        print(f"[CLI] Using database: {db_path_arg}")
    except (IndexError, ValueError):
        print("ERROR: No database path provided after --d/--db flag.")
        sys.exit(1)

# Helper: Get DB path (env > default)
def resolve_db_path():
    # Try to get from environment variable first
    db_path = os.getenv('PANKEGG_DB_PATH')
    if db_path:
        return db_path
    # Fallback: use importlib.resources to locate within package
    try:
        with importlib.resources.path('data', 'pankegg.db') as p:
            return str(p)
    except Exception:
        # Fallback to relative path if not installed as package
        return os.path.join('data', 'pankegg.db')


app = Flask(__name__)
app.secret_key = 'local'  # Change this for production!


def get_db_connection():
    db_path = resolve_db_path()
    print(f"[Flask DEBUG] Opening DB at: {os.path.abspath(db_path)}")
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    print_db_status(conn)
    return conn

def handle_sql_error(e):
    err = str(e)
    if "no such table" in err:
        return (f"<h2>[ERROR] Table missing in database: {e}</h2>"
                f"<p>Check if the database file is correct, up-to-date, and contains all tables. "
                f"Rebuild the DB if necessary.</p>"), 500
    else:
        return f"<h2>[ERROR] Database operational error: {e}</h2>", 500

def taxon_unclassified(val):
    v = val or ''
    if v in ['', 'G__', 'g__', 'S__', 's__', 'D__', 'd__', 'P__', 'p__', 'C__', 'c__', 'O__', 'o__', 'F__', 'f__']:
        return '<span class="text-muted fst-italic">Unclassified</span>'
    return v
app.jinja_env.filters['taxon_unclassified'] = taxon_unclassified


@app.errorhandler(404)
def page_not_found(error):
    return "<h2>404 Not Found</h2><p>The item you requested was not found.</p>", 404


@app.route('/get_taxonomy_data', methods=['POST'])
def get_taxonomy_data():
    try:
        rank = request.form.get('rank')

        if rank not in ['_kingdom_', '_phylum_', '_class_', '_order_', '_family_', '_genus_', '_species_']:
            return jsonify({'error': 'Invalid taxonomic rank selected.'})

        conn = get_db_connection()
        cur = conn.cursor()

        # Query to get the count of bins with and without classification
        bins_count_query = """
        SELECT s.sample_name, 
            COUNT(CASE WHEN t.id IS NOT NULL THEN 1 END) AS classified_bins,
            COUNT(*) AS total_bins
        FROM bin b
        JOIN sample s ON b.sample_id = s.id
        LEFT JOIN taxonomy t ON b.taxonomic_id = t.id
        GROUP BY s.sample_name
        """
        cur.execute(bins_count_query)
        bins_counts = cur.fetchall()

        bins_count_data = {row[0]: {'classified_bins': row[1], 'total_bins': row[2]} for row in bins_counts}

        # Query to get the taxonomic composition data
        query = f"""
        SELECT s.sample_name, t.{rank}, COUNT(*)
        FROM bin b
        JOIN sample s ON b.sample_id = s.id
        JOIN taxonomy t ON b.taxonomic_id = t.id
        GROUP BY s.sample_name, t.{rank}
        """
        cur.execute(query)
        rows = cur.fetchall()
        cur.close()
        conn.close()

        data = {}
        sample_totals = {}

        for row in rows:
            sample = row[0]
            taxon = row[1] if row[1] else 'Unknown'
            count = row[2]
            if sample not in data:
                data[sample] = {}
                sample_totals[sample] = 0
            data[sample][taxon] = count
            sample_totals[sample] += count

        for sample in data:
            for taxon in data[sample]:
                data[sample][taxon] = (data[sample][taxon] / sample_totals[sample]) * 100

        # Include bins count data in the response
        for sample in data:
            data[sample]['bins_info'] = bins_count_data[sample]

        return jsonify(data)
    except sqlite3.OperationalError as e:
        return handle_sql_error(e)


@app.route('/get_heatmap_data_for_taxonomy', methods=['POST'])
def get_heatmap_data_for_taxonomy():
    try:
        sample = request.form.get('sample')
        rank = request.form.get('rank')

        if rank not in ['_kingdom_', '_phylum_', '_class_', '_order_', '_family_', '_genus_', '_species_']:
            return jsonify({'error': 'Invalid taxonomic rank selected.'})

        conn = get_db_connection()
        cur = conn.cursor()

        query = f"""
        SELECT m.map_number, t.{rank}, COUNT(*)
        FROM bin b
        JOIN sample s ON b.sample_id = s.id
        JOIN taxonomy t ON b.taxonomic_id = t.id
        JOIN bin_map bm ON b.id = bm.bin_id
        JOIN map m ON bm.map_id = m.id
        WHERE s.sample_name = ? AND m.map_number NOT LIKE '%not_mapped%'
        GROUP BY m.map_number, t.{rank}
        """
        cur.execute(query, (sample,))
        rows = cur.fetchall()
        cur.close()
        conn.close()

        if not rows:
            return jsonify({'error': 'No data found for the selected sample and taxonomic rank.'})

        df = pd.DataFrame(rows, columns=['map_number', 'taxon', 'count'])

        # Filter out empty and unassigned taxon values
        df = df[(df['taxon'] != '') & (df['taxon'] != 'unassigned')]

        # Create a pivot table for the heatmap data
        heatmap_data = df.pivot_table(index='map_number', columns='taxon', values='count', fill_value=0)

        # Standardize the data before clustering
        scaler = StandardScaler()
        heatmap_data_scaled = scaler.fit_transform(heatmap_data)

        # Perform hierarchical clustering
        row_clusters = linkage(pdist(heatmap_data_scaled, metric='euclidean'), method='ward')
        col_clusters = linkage(pdist(heatmap_data_scaled.T, metric='euclidean'), method='ward')

        # Create a dendrogram and get the order of rows and columns
        row_dendrogram = dendrogram(row_clusters, no_plot=True)
        col_dendrogram = dendrogram(col_clusters, no_plot=True)

        ordered_rows = [heatmap_data.index[i] for i in row_dendrogram['leaves']]
        ordered_cols = [heatmap_data.columns[i] for i in col_dendrogram['leaves']]

        heatmap_data = heatmap_data.loc[ordered_rows, ordered_cols]

        z = heatmap_data.values.tolist()
        x = heatmap_data.columns.tolist()
        y = heatmap_data.index.tolist()

        return jsonify({'heatmap_data': {'z': z, 'x': x, 'y': y}})
    except sqlite3.OperationalError as e:
        return handle_sql_error(e)


@app.route('/taxonomy_comparison', methods=['GET'])
def taxonomy_comparison():
    return render_template('./taxonomy_comparison.html')


#  All for the bin vs bin page ------------------------------------


@app.route('/get_heatmap_data_for_bins', methods=['POST'])
def get_heatmap_data_for_bins():
    try:
        bin1 = request.form.get('bin1')
        bin2 = request.form.get('bin2')

        conn = get_db_connection()
        cur = conn.cursor()

        query = """
        SELECT m.pathway_name,
            COUNT(CASE WHEN b.bin_name = ? THEN 1 END) AS bin1_count,
            COUNT(CASE WHEN b.bin_name = ? THEN 1 END) AS bin2_count
        FROM map m
        JOIN bin_map bm ON m.id = bm.map_id
        JOIN bin b ON bm.bin_id = b.id
        GROUP BY m.pathway_name
        """
        cur.execute(query, (bin1, bin2))
        rows = cur.fetchall()
        cur.close()
        conn.close()

        data = {
            "pathway_names": [],
            "bin1_counts": [],
            "bin2_counts": []
        }

        for row in rows:
            data["pathway_names"].append(row['pathway_name'])
            data["bin1_counts"].append(row['bin1_count'])
            data["bin2_counts"].append(row['bin2_count'])

        return jsonify(data)
    except sqlite3.OperationalError as e:
        return handle_sql_error(e)


@app.route('/bin_vs_bin', methods=['GET'])
def bin_vs_bin():
    return render_template('bin_vs_bin.html', pathway_groups=pathway_groups)


@app.route('/get_bins', methods=['POST'])
def get_bins():
    try:
        sample_name = request.form.get('sample')

        conn = get_db_connection()
        cur = conn.cursor()
        cur.execute("SELECT bin_name FROM bin JOIN sample ON bin.sample_id = sample.id WHERE sample.sample_name = ?",
                    (sample_name,))
        bins = [row[0] for row in cur.fetchall()]
        cur.close()
        conn.close()

        return jsonify(bins=bins)
    except sqlite3.OperationalError as e:
        return handle_sql_error(e)


@app.route('/get_common_pathways_for_bins', methods=['POST'])
def get_common_pathways_for_bins():
    try:
        bin1 = request.form.get('bin1')
        bin2 = request.form.get('bin2')

        conn = get_db_connection()
        cur = conn.cursor()

        query = """
        SELECT m.pathway_name,
            COUNT(CASE WHEN b.bin_name = ? THEN 1 END) AS bin1_count,
            COUNT(CASE WHEN b.bin_name = ? THEN 1 END) AS bin2_count,
            COUNT(CASE WHEN b.bin_name IN (?, ?) THEN 1 END) AS both_count
        FROM map m
        JOIN bin_map bm ON m.id = bm.map_id
        JOIN bin b ON bm.bin_id = b.id
        WHERE m.pathway_name != 'Description not available'
        GROUP BY m.pathway_name
        """
        cur.execute(query, (bin1, bin2, bin1, bin2))
        rows = cur.fetchall()
        cur.close()
        conn.close()

        data = []
        for row in rows:
            data.append({
                'pathway_name': row['pathway_name'],
                'both': row['both_count'],
                'bin1': row['bin1_count'],
                'bin2': row['bin2_count']
            })

        return jsonify(data)
    except sqlite3.OperationalError as e:
        return handle_sql_error(e)


#  All for the sample vs sample page ------------------------------------


def parse_pathway_groups(file_path):
    try:
        pathway_groups = {}
        current_group = None

        with open(file_path, 'r') as file:
            for line in file:
                line = line.strip()
                if not line:
                    continue
                if line[0].isdigit() and '.' in line:
                    # This is a group header
                    current_group = line
                    pathway_groups[current_group] = []
                else:
                    # This is a map number
                    if current_group is not None:
                        pathway_groups[current_group].append(f"map{line}")
                    else:
                        print(f"Warning: Found map number '{line}' without a group header.")

        return pathway_groups
    except sqlite3.OperationalError as e:
        return handle_sql_error(e)


# Initialize pathway_groups from the file
try:
    resource = files('data').joinpath('pathway_groups.txt')
    pathway_groups = parse_pathway_groups(str(resource))
except Exception:
    rel_path = os.path.join('data', 'pathway_groups.txt')
    pathway_groups = parse_pathway_groups(rel_path)


# Route to get list of samples


@app.route('/get_samples', methods=['GET'])
def get_samples():
    try:
        conn = get_db_connection()
        cur = conn.cursor()
        cur.execute("SELECT sample_name FROM sample")
        samples = [row[0] for row in cur.fetchall()]
        cur.close()
        conn.close()
        return jsonify(samples=samples)
    except sqlite3.OperationalError as e:
        return handle_sql_error(e)


@app.route('/get_bins_list', methods=['GET'])
def get_bins_list():
    try:
        filter_mode = request.args.get('mode', 'sample')  # 'sample' or 'bin'
        conn = get_db_connection()
        cur = conn.cursor()
        if filter_mode == 'sample':
            cur.execute("SELECT sample_name FROM sample")
            items = [row[0] for row in cur.fetchall()]
        else:  # 'bin'
            cur.execute("""
                SELECT bin.bin_name || ' (Sample: ' || sample.sample_name || ')' AS bin_sample
                FROM bin
                JOIN sample ON bin.sample_id = sample.id
            """)
            items = [row[0] for row in cur.fetchall()]
        cur.close()
        conn.close()
        return jsonify(items=items)
    except Exception as e:
        return jsonify(items=[])


# Route pour obtenir les données de la heatmap
@app.route('/get_heatmap_data', methods=['POST'])
def get_heatmap_data():
    try:
        sample = request.form.get('sample')
        selected_groups = json.loads(request.form.get('groups', '[]'))

        conn = get_db_connection()
        cur = conn.cursor()
        query = """
        SELECT bin.bin_name, m.map_number, m.pathway_name, COUNT(DISTINCT k.ko_id) as ko_count, m.pathway_total_orthologs
        FROM map m
        LEFT JOIN map_kegg mk ON m.id = mk.map_id
        LEFT JOIN kegg k ON mk.kegg_id = k.id
        JOIN bin_map_kegg bmk ON mk.id = bmk.map_kegg_id
        JOIN bin ON bmk.bin_id = bin.id
        JOIN sample ON bin.sample_id = sample.id
        WHERE sample.sample_name = ? AND mk.real_pathway_id = 1
        GROUP BY bin.bin_name, m.map_number, m.pathway_name
        """
        cur.execute(query, (sample,))
        rows = cur.fetchall()
        cur.close()
        conn.close()

        if not rows:
            return jsonify({'error': 'No data found for the selected sample.'})

        df = pd.DataFrame(rows, columns=['bin_name', 'map_number', 'pathway_name', 'ko_count', 'pathway_total_orthologs'])
        df['completion_percentage'] = (df['ko_count'] / df['pathway_total_orthologs']) * 100

        heatmaps_data = {}
        for group, maps in pathway_groups.items():
            if group not in selected_groups:
                continue
            df_group = df[df['map_number'].isin(maps)]
            if not df_group.empty:
                df_group.set_index('map_number', inplace=True)
                heatmap_data = df_group.pivot_table(index='map_number', columns='bin_name', values='completion_percentage', fill_value=0)
                z = heatmap_data.values.tolist()
                x = heatmap_data.columns.tolist()
                y = heatmap_data.index.tolist()
                heatmaps_data[group] = {'z': z, 'x': x, 'y': y}

        return jsonify(heatmaps_data=heatmaps_data, context=f'Heatmap for Sample: {sample}')
    except sqlite3.OperationalError as e:
        return handle_sql_error(e)


@app.route('/get_map_details', methods=['POST'])
def get_map_details():
    try:
        data = request.get_json()
        map_numbers = data.get('map_numbers', [])
        print("mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm")
        print(map_numbers)
        if not map_numbers:
            return jsonify({'error': 'No map numbers provided.'})

        conn = get_db_connection()
        cur = conn.cursor()
        query = """
        SELECT map_number, pathway_name
        FROM map
        WHERE map_number IN ({})
        """.format(','.join('?' * len(map_numbers)))

        cur.execute(query, map_numbers)
        rows = cur.fetchall()
        cur.close()
        conn.close()

        map_details = {row['map_number']: row['pathway_name'] for row in rows}
        print(map_details)
        return jsonify(map_details=map_details)
    except sqlite3.OperationalError as e:
        return handle_sql_error(e)


@app.route('/get_scatterplot_data', methods=['POST'])
def get_scatterplot_data():
    try:
        sample = request.form.get('sample')

        conn = get_db_connection()
        cur = conn.cursor()
        query = """
        SELECT bin.bin_name, bin.completeness, bin.contamination
        FROM bin
        JOIN sample ON bin.sample_id = sample.id
        WHERE sample.sample_name = ?
        """
        cur.execute(query, (sample,))
        rows = cur.fetchall()
        cur.close()
        conn.close()

        if not rows:
            return jsonify({'error': 'No data found for the selected sample.'})

        df = pd.DataFrame(rows, columns=['bin_name', 'completeness', 'contamination'])

        in_filter = df[df['completeness'] - 5 * df['contamination'] > 50]
        out_filter = df[df['completeness'] - 5 * df['contamination'] <= 50]

        data = {
            'inFilter': {
                'x': in_filter['completeness'].tolist(),
                'y': in_filter['contamination'].tolist(),
                'text': in_filter['bin_name'].tolist()
            },
            'outFilter': {
                'x': out_filter['completeness'].tolist(),
                'y': out_filter['contamination'].tolist(),
                'text': out_filter['bin_name'].tolist()
            }
        }

        return jsonify(data)
    except sqlite3.OperationalError as e:
        return handle_sql_error(e)


@app.route('/get_common_pathways_data', methods=['POST'])
def get_common_pathways_data():
    try:
        sample1 = request.form.get('sample1')
        sample2 = request.form.get('sample2')

        conn = get_db_connection()
        cur = conn.cursor()

        query = """
        SELECT m.pathway_name,
            COUNT(CASE WHEN s.sample_name = ? THEN 1 END) AS sample1_count,
            COUNT(CASE WHEN s.sample_name = ? THEN 1 END) AS sample2_count,
            COUNT(CASE WHEN s.sample_name IN (?, ?) THEN 1 END) AS both_count
        FROM map m
        JOIN bin_map bm ON m.id = bm.map_id
        JOIN bin b ON bm.bin_id = b.id
        JOIN sample s ON b.sample_id = s.id
        WHERE m.pathway_name != 'Description not available'
        GROUP BY m.pathway_name
        """
        cur.execute(query, (sample1, sample2, sample1, sample2))
        rows = cur.fetchall()
        cur.close()
        conn.close()

        data = []
        for row in rows:
            data.append({
                'pathway_name': row['pathway_name'],
                'both': row['both_count'],
                'sample1': row['sample1_count'],
                'sample2': row['sample2_count']
            })

        return jsonify(data)
    except sqlite3.OperationalError as e:
        return handle_sql_error(e)


@app.route('/get_pca_data', methods=['POST'])
def get_pca_data():
    try:
        sample = request.form.get('sample')

        conn = get_db_connection()
        cur = conn.cursor()

        # Récupérer les bins et les maps associés pour le sample sélectionné
        query = """
        SELECT bin.bin_name, map.map_number
        FROM bin
        JOIN bin_map ON bin.id = bin_map.bin_id
        JOIN map ON bin_map.map_id = map.id
        JOIN sample ON bin.sample_id = sample.id
        WHERE sample.sample_name = ?
        """
        cur.execute(query, (sample,))
        rows = cur.fetchall()
        cur.close()
        conn.close()

        if not rows:
            return jsonify({'error': 'No data found for the selected sample.'})

        # Créer une matrice binaire
        bin_names = sorted(set(row[0] for row in rows))
        map_numbers = sorted(set(row[1] for row in rows))
        bin_index = {bin_name: i for i, bin_name in enumerate(bin_names)}
        map_index = {map_number: i for i, map_number in enumerate(map_numbers)}

        matrix = np.zeros((len(bin_names), len(map_numbers)))

        for row in rows:
            bin_name, map_number = row
            matrix[bin_index[bin_name], map_index[map_number]] = 1

        # Calculer les clusters K-means
        n_clusters = 3  # Par exemple, on choisit de regrouper les bins en 3 clusters
        kmeans = KMeans(n_clusters=n_clusters, random_state=42)
        clusters = kmeans.fit_predict(matrix)

        # Calculer la PCA
        pca = PCA(n_components=2)
        pca_result = pca.fit_transform(matrix)

        pca_data = {
            'x': pca_result[:, 0].tolist(),
            'y': pca_result[:, 1].tolist(),
            'bins': bin_names,
            'clusters': clusters.tolist()
        }

        return jsonify(pca_data)
    except sqlite3.OperationalError as e:
        return handle_sql_error(e)


@app.route('/sample_vs_sample', methods=['GET'])
def sample_vs_sample():
    return render_template('sample_vs_sample.html', pathway_groups=pathway_groups)


#  sep--------------------------------------------------------------------------------------------------------


@app.route('/compar_full', methods=['GET', 'POST'])
def compar_full():
    try:
        if request.method == 'POST':
            pca_type = request.form.get('pca_type')
            taxonomy_level = request.form.get('taxonomy_level', '_species_')

            conn = get_db_connection()
            cur = conn.cursor()

            if pca_type == 'maps':
                query = """
                SELECT sample.sample_name, map.map_number, COUNT(bin.id) AS count_bins
                FROM sample
                JOIN bin ON sample.id = bin.sample_id
                JOIN bin_map ON bin.id = bin_map.bin_id
                JOIN map ON bin_map.map_id = map.id
                GROUP BY sample.sample_name, map.map_number
                """
                context = "PCA based on Maps"
            elif pca_type == 'kos':
                query = """
                SELECT sample.sample_name, k.ko_id, COUNT(DISTINCT bin.id) AS count_bins
                FROM sample
                JOIN bin ON sample.id = bin.sample_id
                JOIN bin_map_kegg bmk ON bin.id = bmk.bin_id
                JOIN map_kegg mk ON bmk.map_kegg_id = mk.id
                JOIN kegg k ON mk.kegg_id = k.id
                GROUP BY sample.sample_name, k.ko_id
                """
                context = "PCA based on KOs"
            elif pca_type == 'taxonomy':
                query = f"""
                SELECT sample.sample_name, t.{taxonomy_level} AS taxon, COUNT(bin.id) AS count_bins
                FROM sample
                JOIN bin ON sample.id = bin.sample_id
                JOIN taxonomy t ON bin.taxonomic_id = t.id
                GROUP BY sample.sample_name, t.{taxonomy_level}
                """
                context = "PCA based on Taxonomy level : " + taxonomy_level
            else:
                return jsonify({
                    'context': 'Error',
                    'pca_data': [],
                    'error': 'Invalid PCA type selected.'
                })

            cur.execute(query)
            rows = cur.fetchall()
            cur.close()
            conn.close()

            # Convert data to DataFrame
            if pca_type == 'maps':
                df = pd.DataFrame(rows, columns=['sample_name', 'map_number', 'count_bins'])
            elif pca_type == 'kos':
                df = pd.DataFrame(rows, columns=['sample_name', 'ko_id', 'count_bins'])
            elif pca_type == 'taxonomy':
                df = pd.DataFrame(rows, columns=['sample_name', 'taxon', 'count_bins'])

            # Debug: Show the DataFrame before pivot
            print("DataFrame before pivot:\n", df.head())

            if pca_type == 'maps':
                df_pivot = df.pivot_table(index='sample_name', columns='map_number', values='count_bins', fill_value=0)
            elif pca_type == 'kos':
                df_pivot = df.pivot_table(index='sample_name', columns='ko_id', values='count_bins', fill_value=0)
            elif pca_type == 'taxonomy':
                df_pivot = df.pivot_table(index='sample_name', columns='taxon', values='count_bins', fill_value=0)

            # Debug: Show the pivoted DataFrame
            print("Pivoted DataFrame:\n", df_pivot.head())

            # Remove constant columns
            df_pivot = df_pivot.loc[:, (df_pivot != df_pivot.iloc[0]).any()]

            # Debug: Show the DataFrame after removing constant columns
            print("DataFrame after removing constant columns:\n", df_pivot.head())

            # Check if dataframe is empty after removing constant columns
            if df_pivot.empty:
                return jsonify({
                    'context': context,
                    'pca_data': [],
                    'error': 'No variability in data after removing constant columns.'
                })

            # Normalization
            scaler = StandardScaler()
            df_pivot_normalized = scaler.fit_transform(df_pivot)

            # Debug: Show the normalized DataFrame
            print("Normalized DataFrame:\n", df_pivot_normalized)

            # Perform PCA
            pca = PCA(n_components=2)
            try:
                pca_results = pca.fit_transform(df_pivot_normalized)
                explained_variance = pca.explained_variance_ratio_
            except Exception as e:
                return jsonify({
                    'context': context,
                    'pca_data': [],
                    'error': str(e)
                })

            # Debug: Show the PCA results and explained variance
            print("PCA Results:\n", pca_results)
            print("Explained Variance Ratio:\n", explained_variance)

            # Convert PCA results to DataFrame for visualization
            pca_df = pd.DataFrame(pca_results, columns=['PC1', 'PC2'])
            pca_df['sample_name'] = df_pivot.index

            # Debug: Show the PCA DataFrame
            print("PCA DataFrame:\n", pca_df)

            # Return JSON response for AJAX
            return jsonify({
                'context': context,
                'pca_data': pca_df.to_dict(orient='records'),
                'explained_variance': explained_variance.tolist()
            })

        return render_template('compar_full.html')
    except sqlite3.OperationalError as e:
        return handle_sql_error(e)


@app.route('/comparison')
def main_compare():
    return render_template('comp_menu.html', content="Comp")


#  Export functions ===============================================================================================================


@app.route('/export_bins', methods=['POST'])
def export_bins():
    try:
        map_number = request.values.get('map_number')
        kegg_id = request.values.get('kegg_id')
        taxon = request.values.get('taxon')
        search_query = request.values.get('search_query')
        gtdb_filter = session.get('gtdb_filter', False)
        sort_filter = session.get('selected_sort_option', False)
        bin_name_sort_sql_command = " ORDER BY bin_number ASC"
        completeness_sort_sql_command = " ORDER BY bin.completeness DESC"
        contamination_sort_sql_command = " ORDER BY bin.contamination DESC"

        context = "Display of all bins"  # Initialisation par défaut du contexte

        conn = get_db_connection()
        cur = conn.cursor()

        # Préparer les colonnes à récupérer de chaque table
        bin_columns = ['bin.id as bin_id', 'bin.bin_name', 'bin.completeness', 'bin.contamination', "sample.sample_name",
                    "CAST(SUBSTR(bin_name, INSTR(bin_name, '.') + 1) AS INTEGER) AS bin_number"]
        taxonomy_columns = [
            'taxonomy."_kingdom_" as kingdom',
            'taxonomy."_phylum_" as phylum',
            'taxonomy."_class_" as class',
            'taxonomy."_order_" as "order"',
            'taxonomy."_family_" as family',
            'taxonomy."_genus_" as genus',
            'taxonomy."_species_" as species'
        ]

        # Joindre les tables bin et taxonomy
        join_query = "LEFT JOIN taxonomy ON bin.taxonomic_id = taxonomy.id"

        # Construire la requête de base avec DISTINCT pour éviter les duplications
        query = f"SELECT DISTINCT {', '.join(bin_columns + taxonomy_columns)} FROM bin {join_query}"
        query += " JOIN sample ON sample.id = bin.sample_id "

        conditions = []
        params = []

        # Ajouter des conditions basées sur map_number ou kegg_id
        if map_number:
            query += " JOIN bin_map_kegg bmk ON bin.id = bmk.bin_id"
            query += " JOIN map_kegg mk ON bmk.map_kegg_id = mk.id"
            query += " JOIN map m ON mk.map_id = m.id"
            conditions.append("m.map_number = ?")
            params.append(map_number)
            # Fetch pathway name for nicer context
            cur.execute("SELECT pathway_name FROM map WHERE map_number = ?", (map_number,))
            map_row = cur.fetchone()
            pathway_name = map_row[0] if map_row else None
            if pathway_name:
                context = f"Display of bins for Pathway: {pathway_name}"
            else:
                context = f"Display of bins for Map number: {map_number}"
        elif kegg_id:
            # Modification ici pour utiliser ko_id
            query += " JOIN bin_map_kegg bmk ON bin.id = bmk.bin_id"
            query += " JOIN map_kegg mk ON bmk.map_kegg_id = mk.id"
            query += " JOIN kegg k ON mk.kegg_id = k.id"
            conditions.append("k.ko_id = ?")
            params.append(kegg_id)
            context = f"Display of bins for KEGG ID: {kegg_id}"
        elif taxon:
            # Modification ici pour utiliser taxonomy_id
            conditions.append("? IN (kingdom, phylum, class, \"order\", family, genus, species)")
            params.append(taxon if taxon != "none" else "")
            context = f"Display of bins for taxonomy entry: {taxon}"

        if gtdb_filter:
            conditions.append("(completeness - 5 * contamination > 50)")
        if search_query:
            search_condition = "(sample.sample_name LIKE ? OR bin.bin_name LIKE ?)"
            conditions.append(search_condition)
            search_pattern = f"%{search_query}%"
            params.extend([search_pattern, search_pattern])
            context += f" with search pattern: {search_query}"

        if conditions:
            query += " WHERE " + " AND ".join(conditions)

        if sort_filter == "option1":
            query += bin_name_sort_sql_command
        elif sort_filter == "option2":
            query += completeness_sort_sql_command
        else:
            query += contamination_sort_sql_command

        cur.execute(query, params)

        rows = cur.fetchall()
        cur.close()
        conn.close()

        # Préparer les noms des colonnes pour l'affichage
        display_column_labels = [col.split(' as ')[1] if ' as ' in col else col.split('.')[1] for col in bin_columns]

        # Organiser les données pour le CSV
        bins = []
        sample_bins = {}
        for row in rows:
            bin_data = dict(
                zip(display_column_labels + ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'], row))
            bins.append(bin_data)
            sample_name = bin_data['sample_name']
            if sample_name not in sample_bins:
                sample_bins[sample_name] = []
            sample_bins[sample_name].append(bin_data['bin_name'])

        # Générer le CSV
        output = io.StringIO()
        writer = csv.writer(output)
        writer.writerow(["Pankegg", "bin"])
        writer.writerow(["Filter", context])
        for sample_name, bin_names in sample_bins.items():
            writer.writerow([sample_name, ', '.join(bin_names)])

        response = Response(output.getvalue(), mimetype='text/csv')
        response.headers.set("Content-Disposition", "attachment", filename="bins_export.csv")
        return response
    except sqlite3.OperationalError as e:
        return handle_sql_error(e)


@app.route('/export_maps', methods=['POST'])
def export_maps():
    try:
        bin_id = request.values.get('bin_id')
        ko_id = request.values.get('ko_id')
        taxon = request.values.get('taxon')
        search_query = request.values.get('search_query')

        context = "Display of all maps"  # Initialisation par défaut du contexte

        conn = get_db_connection()
        cur = conn.cursor()

        base_query = """
            SELECT m.map_number, m.pathway_name, k.ko_id, k.kegg_name, m.pathway_total_orthologs, mk.real_pathway_id
            FROM map m
            LEFT JOIN map_kegg mk ON m.id = mk.map_id
            LEFT JOIN kegg k ON mk.kegg_id = k.id
        """
        conditions = []
        params = []

        if bin_id:
            # Fetch both bin_name and sample_name
            cur.execute("""
                SELECT bin.bin_name, sample.sample_name
                FROM bin
                JOIN sample ON bin.sample_id = sample.id
                WHERE bin.id = ?
            """, (bin_id,))
            bin_info = cur.fetchone()
            if bin_info:
                bin_name, sample_name = bin_info
            else:
                bin_name, sample_name = "Unknown bin", "Unknown sample"
            context = f"Maps associated with bin: <strong>{bin_name}</strong> (Sample: <strong>{sample_name}</strong>)"
            base_query += """
                JOIN bin_map_kegg bmk ON mk.id = bmk.map_kegg_id
            """
            conditions.append("bmk.bin_id = ?")
            params.append(bin_id)
        elif ko_id:
            context = f"Maps containing the KEGG identifier {ko_id}"
            base_query = """
                    SELECT m.map_number, m.pathway_name, k2.ko_id, k2.kegg_name, m.pathway_total_orthologs, mk2.real_pathway_id
                    FROM map m
                    LEFT JOIN map_kegg mk ON m.id = mk.map_id
                    LEFT JOIN kegg k ON mk.kegg_id = k.id
                    LEFT JOIN map_kegg mk2 ON m.id = mk2.map_id
                    LEFT JOIN kegg k2 ON mk2.kegg_id = k2.id
                    """
            conditions.append("k.ko_id = ?")
            params.append(ko_id)
        elif taxon:
            context = f"Maps associated with taxonomy: {taxon}"
            base_query += """
                JOIN bin_map_kegg bmk ON mk.id = bmk.map_kegg_id
                JOIN bin bi ON bi.id = bmk.bin_id
                JOIN taxonomy t ON bi.taxonomic_id = t.id
            """
            conditions.append("? IN (t._kingdom_, t._phylum_, t._class_, t._order_, t._family_, t._genus_, t._species_)")
            if taxon == "none":
                params.append("")
            else:
                params.append(taxon)
        else:
            context = "All Maps"

        if search_query:
            search_condition = "(m.map_number LIKE ? OR m.pathway_name LIKE ?)"
            conditions.append(search_condition)
            search_pattern = f"%{search_query}%"
            params.extend([search_pattern, search_pattern])
            context += f" with search pattern: {search_query}"

        if conditions:
            base_query += " WHERE " + " AND ".join(conditions)

        cur.execute(base_query, params)

        rows = cur.fetchall()
        cur.close()
        conn.close()

        # Organiser les données pour le CSV
        maps = []
        map_kos = {}
        for row in rows:
            map_key = (row[0], row[1])  # map_number, pathway_name
            if map_key not in map_kos:
                map_kos[map_key] = []
            map_kos[map_key].append(row[2])  # ko_id

        # Générer le CSV
        output = io.StringIO()
        writer = csv.writer(output)
        writer.writerow(["Pankegg", "maps"])
        writer.writerow(["Filter", context])
        for map_key, ko_ids in map_kos.items():
            # Filtrer les None avant de joindre les ko_ids
            ko_ids = [ko_id for ko_id in ko_ids if ko_id is not None]
            writer.writerow([map_key[0], ', '.join(ko_ids)])

        response = Response(output.getvalue(), mimetype='text/csv')
        response.headers.set("Content-Disposition", "attachment", filename="maps_export.csv")
        return response
    except sqlite3.OperationalError as e:
        return handle_sql_error(e)


@app.route('/export_kegg', methods=['POST'])
def export_kegg():
    try:
        ko_id = request.values.get('ko_id')
        bin_id = request.values.get('bin_id')
        taxon = request.values.get('taxon')
        search_query = request.values.get('search_query')

        context = "Display of all KEGG IDs"  # Initialisation par défaut du contexte

        conn = get_db_connection()
        cur = conn.cursor()

        query = """
        SELECT k.ko_id, k.kegg_name, k.kegg_full_name
        FROM kegg k
        LEFT JOIN bin_extra_kegg bek ON k.id = bek.kegg_id
        LEFT JOIN bin_extra be ON bek.extra_id = be.id
        LEFT JOIN bin b ON be.bin_id = b.id
        """

        conditions = []
        params = []

        if ko_id:
            context = f"Display for KEGG ID: {ko_id}"
            conditions.append("k.ko_id = ?")
            params.append(ko_id)
        elif bin_id:
            cur.execute("SELECT bin_name FROM bin WHERE id = ?", (bin_id,))
            bin_name_result = cur.fetchone()
            bin_name = bin_name_result[0] if bin_name_result else "Unknown bin"
            context = f"KEGG inputs associated with {bin_name}"
            conditions.append("b.id = ?")
            params.append(bin_id)
        elif taxon:
            context = f"KEGG inputs associated with taxonomy: {taxon}"
            query += " JOIN taxonomy t ON t.id = b.taxonomic_id"
            conditions.append("? IN (t._kingdom_, t._phylum_, t._class_, t._order_, t._family_, t._genus_, t._species_)")
            if taxon == "none":
                params.append("")
            else:
                params.append(taxon)

        if search_query:
            search_condition = "(k.ko_id LIKE ? OR k.kegg_name LIKE ? OR k.kegg_full_name LIKE ?)"
            conditions.append(search_condition)
            search_pattern = f"%{search_query}%"
            params.extend([search_pattern, search_pattern, search_pattern])
            context += f" with search pattern: {search_query}"

        if conditions:
            query += " WHERE " + " AND ".join(conditions)

        cur.execute(query, params)

        rows = cur.fetchall()
        cur.close()
        conn.close()

        # Organiser les données pour le CSV
        kegg_entries = []
        for row in rows:
            kegg_entries.append(row)  # ko_id, kegg_name, kegg_full_name

        # Générer le CSV
        output = io.StringIO()
        writer = csv.writer(output)
        writer.writerow(["Pankegg", "kegg"])
        writer.writerow(["Filter", context])
        for entry in kegg_entries:
            writer.writerow(entry)

        response = Response(output.getvalue(), mimetype='text/csv')
        response.headers.set("Content-Disposition", "attachment", filename="kegg_export.csv")
        return response
    except sqlite3.OperationalError as e:
        return handle_sql_error(e)

#  working======================================================================================================


@app.route('/bin_query')
def show_bins2():
    try:
        # Obtenir les colonnes demandées, sinon utiliser un ensemble par défaut
        requested_columns = request.args.getlist('columns')
        # print(requested_columns)
        default_columns = ['id', 'bin_name', 'completeness', 'contamination', 'taxonomic_id']
        columns = requested_columns if requested_columns else default_columns

        # Construire une requête SQL sécurisée en vérifiant que chaque colonne demandée est valide
        safe_columns = [col for col in columns if col in default_columns]
        if not safe_columns:
            safe_columns = default_columns  # Utiliser les colonnes par défaut si aucune colonne demandée n'est valide

        conn = get_db_connection()
        cur = conn.cursor()
        query = f"SELECT {', '.join(safe_columns)} FROM bin"
        # print(query)
        cur.execute(query)
        bins = cur.fetchall()
        conn.close()
        return render_template('bin.html', bins=bins, columns=safe_columns)
    except sqlite3.OperationalError as e:
        return handle_sql_error(e)

@app.route('/taxonomy')
def taxonomy():
    try:
        level = request.args.get('level')
        conn = get_db_connection()
        cur = conn.cursor()

        if level:
            # Mise à jour de la requête pour joindre les tables et compter les bins associés
            query = f"""
            SELECT TRIM(LOWER(t._{level}_)) AS taxon, COUNT(b.id) AS bins_associated
            FROM taxonomy t
            LEFT JOIN bin b ON t.id = b.taxonomic_id
            GROUP BY taxon
            ORDER BY taxon
            """
            cur.execute(query)
            results = cur.fetchall()
            taxons = [{'name': row['taxon'].capitalize() if row['taxon'] else 'none',
                    'count': row['bins_associated']} for row in results]
        else:
            taxons = []

        cur.close()
        conn.close()

        return render_template('taxonomy.html', taxons=taxons, level=level or "none")
    except sqlite3.OperationalError as e:
        return handle_sql_error(e)

@app.route('/kegg', methods=['GET', 'POST'])
def kegg():
    def normalize_taxon_case(taxon):
        # Handles cases like O__christensenellales → o__Christensenellales
        if not taxon or len(taxon) < 4:
            return taxon
        prefix = taxon[:3].lower()
        rest = taxon[3:]
        if len(rest) == 0:
            return prefix
        return prefix + rest[0].upper() + rest[1:]

    try:
        ko_id = request.values.get('ko_id')
        bin_id = request.values.get('bin_id')
        taxon = request.values.get('taxon')
        search_query = request.form.get('search_query', '')

        conn = get_db_connection()
        cur = conn.cursor()
        kegg_entries = {}
        bin_name = None

        query = """
        SELECT k.ko_id, k.kegg_name, k.kegg_full_name, b.bin_name, be.go, be.ko, be.eggnog_desc
        FROM kegg k
        LEFT JOIN bin_extra_kegg bek ON k.id = bek.kegg_id
        LEFT JOIN bin_extra be ON bek.extra_id = be.id
        LEFT JOIN bin b ON be.bin_id = b.id
        """
        conditions = []
        params = []

        if ko_id:
            context = f"Display for KEGG ID: <strong>{ko_id}</strong>"
            conditions.append("k.ko_id = ?")
            params.append(ko_id)
        elif bin_id:
            cur.execute("SELECT bin_name FROM bin WHERE id = ?", (bin_id,))
            bin_name_result = cur.fetchone()
            bin_name = bin_name_result[0] if bin_name_result else "Unknown bin"
            context = f"KEGG inputs associated with <strong>{bin_name}</strong>"
            conditions.append("b.id = ?")
            params.append(bin_id)
        elif taxon:
            # --- Begin Taxon Filtering ---
            normalized_taxon = normalize_taxon_case(taxon)
            prefix = normalized_taxon[:3].lower()
            rank_map = {
                'd__': '_kingdom_',
                'p__': '_phylum_',
                'c__': '_class_',
                'o__': '_order_',
                'f__': '_family_',
                'g__': '_genus_',
                's__': '_species_'
            }
            rank = rank_map.get(prefix)
            query += " JOIN taxonomy t ON t.id = b.taxonomic_id "
            if rank in ['_genus_', '_species_']:
                cleaned_taxon = normalized_taxon.strip()
                if cleaned_taxon.lower() in ['g__', 's__']:
                    # Unclassified: match g__/s__, empty, or NULL
                    conditions.append(
                        f"(LOWER(t.{rank}) = ? OR t.{rank} = '' OR t.{rank} IS NULL)"
                    )
                    params.append(cleaned_taxon.lower())
                else:
                    # Normal classified genus/species
                    conditions.append(f"t.{rank} LIKE ?")
                    params.append(cleaned_taxon + '%')
            elif rank:
                conditions.append(f"t.{rank} = ?")
                params.append(normalized_taxon)
            else:
                # fallback: match on any rank as before
                conditions.append(
                    "(t._kingdom_ = ? OR t._phylum_ = ? OR t._class_ = ? "
                    "OR t._order_ = ? OR t._family_ = ? OR t._genus_ = ? OR t._species_ = ?)"
                )
                params.extend([normalized_taxon]*7)
            context = f"KEGG inputs associated with taxonomy: <strong>{taxon}</strong>"
            # --- End Taxon Filtering ---
        else:
            context = "Display of all KEGG IDs"

        if search_query:
            search_condition = "(k.ko_id LIKE ? OR k.kegg_name LIKE ? OR k.kegg_full_name LIKE ?)"
            conditions.append(search_condition)
            search_pattern = f"%{search_query}%"
            params.extend([search_pattern, search_pattern, search_pattern])
            context += f" with search pattern: <strong>{search_query}</strong>"

        if conditions:
            query += " WHERE " + " AND ".join(conditions)

        cur.execute(query, params)

        rows = cur.fetchall()
        for row in rows:
            ko_key = (row[0], row[1], row[2])
            if ko_key not in kegg_entries:
                kegg_entries[ko_key] = []
            go_terms = row[4].split(',') if row[4] else []
            kegg_entries[ko_key].append((row[3], go_terms, row[5], row[6]))

        cur.close()
        conn.close()
        return render_template('kegg.html', context=context, kegg_entries=kegg_entries.items())
    except sqlite3.OperationalError as e:
        return handle_sql_error(e)



@app.route('/map', methods=['GET', 'POST'])
def show_maps():
    try:
        # 1. Gather all filtering params (from POST preferred, fallback GET)
        filter_mode = request.form.get('filter_mode', request.args.get('filter_mode', 'sample'))  # 'sample' or 'bin'
        active_filters = request.form.getlist('active_filters') or request.args.getlist('active_filters')
        bin_id = request.values.get('bin_id')
        ko_id = request.values.get('ko_id')
        taxon = request.values.get('taxon')
        search_query = request.form.get('search_query', '') or request.values.get('search_query', '')

        conn = get_db_connection()
        cur = conn.cursor()

        maps = {}
        context_tags = []
        filters_clause = []
        params = []

        base_query = """
            SELECT m.map_number, m.pathway_name, k.ko_id, k.kegg_name, m.pathway_total_orthologs, mk.real_pathway_id
            FROM map m
            LEFT JOIN map_kegg mk ON m.id = mk.map_id
            LEFT JOIN kegg k ON mk.kegg_id = k.id
            JOIN bin_map_kegg bmk ON mk.id = bmk.map_kegg_id
            JOIN bin b ON bmk.bin_id = b.id
            JOIN sample s ON b.sample_id = s.id
        """

        # --- Taxon-based filtering (if present, apply as a first filter)
        if taxon:
            filters_clause.append(
                "(b.taxonomic_id IN (SELECT id FROM taxonomy WHERE "
                "(LOWER(_kingdom_) = LOWER(?) OR LOWER(_phylum_) = LOWER(?) OR LOWER(_class_) = LOWER(?) "
                "OR LOWER(_order_) = LOWER(?) OR LOWER(_family_) = LOWER(?) OR LOWER(_genus_) = LOWER(?) OR LOWER(_species_) = LOWER(?))))"
            )
            params.extend([taxon]*7)
            context_tags.append(f"<span class='tag-chip'>{taxon}</span>")

        # --- Tag filtering (sample/bin, additive) ---
        if active_filters:
            if filter_mode == "sample":
                filters_clause.append('s.sample_name IN ({})'.format(','.join('?' for _ in active_filters)))
                params.extend(active_filters)
                context_tags += [f"<span class='tag-chip'>{x}</span>" for x in active_filters]
            elif filter_mode == "bin":
                parsed_bins = []
                for entry in active_filters:
                    if ' (Sample: ' in entry:
                        bin_name, rest = entry.split(' (Sample: ')
                        sample_name = rest.replace(')', '')
                    elif ' [' in entry and entry.endswith(']'):
                        bin_name, sample_name = entry[:-1].split(' [')
                    else:
                        bin_name, sample_name = entry, ""
                    parsed_bins.append( (bin_name, sample_name) )
                if parsed_bins:
                    sub_clauses = []
                    for bin_name, sample_name in parsed_bins:
                        sub_clauses.append('(b.bin_name = ? AND s.sample_name = ?)')
                        params.extend([bin_name, sample_name])
                    filters_clause.append('(' + ' OR '.join(sub_clauses) + ')')
                    context_tags += [f"<span class='tag-chip'>{bin_name} [{sample_name}]</span>" for bin_name, sample_name in parsed_bins]
        elif bin_id:
            cur.execute("SELECT bin_name, sample_name FROM bin JOIN sample ON bin.sample_id = sample.id WHERE bin.id = ?", (bin_id,))
            bin_row = cur.fetchone()
            if bin_row:
                bin_name, sample_name = bin_row
                context_tags.append(f"<span class='tag-chip'>{bin_name} [{sample_name}]</span>")
            else:
                context_tags.append(f"<span class='tag-chip'>Bin id: {bin_id}</span>")
            filters_clause.append("b.id = ?")
            params.append(bin_id)

        # --- Pathway search (additive, always works with tag filters) ---
        if search_query:
            filters_clause.append("(m.map_number LIKE ? OR m.pathway_name LIKE ?)")
            params.extend([f"%{search_query}%", f"%{search_query}%"])
            context_tags.append(f"<span class='tag-chip'>{search_query}</span>")

        # --- Final context string for display ---
        context = "Maps associated with search: " + ' '.join(context_tags) if context_tags else "All Maps"

        if filters_clause:
            base_query += " WHERE " + " AND ".join(filters_clause)

        cur.execute(base_query, params)

        # --- Pathway completion calculation as before ---
        map_completions = {}
        for row in cur.fetchall():
            map_key = (row[0], row[1])  # map_number, pathway_name
            if map_key not in maps:
                maps[map_key] = {
                    'kegg_ids': [],
                    'completion': '0.00%'
                }
            if row[2] and row[3]:  # Ensure ko_id and kegg_name are not None
                maps[map_key]['kegg_ids'].append((row[2], row[3], row[5]))  # Include real_pathway_id
            pathway_total_orthologs = row[4]
            if map_key[0] not in map_completions:
                map_completions[map_key[0]] = pathway_total_orthologs

        for map_key in maps.keys():
            kegg_ids = maps[map_key]['kegg_ids']
            filtered_kegg_ids = [kegg_id[0] for kegg_id in kegg_ids if kegg_id[2] == 1]
            count_id = len(set(filtered_kegg_ids))
            pathway_total = map_completions[map_key[0]]
            if pathway_total > 0:
                completion_ratio = count_id / pathway_total * 100
                maps[map_key]['completion'] = completion_ratio
            else:
                maps[map_key]['completion'] = None

        cur.close()
        conn.close()

        return render_template(
            'maps.html',
            maps=maps.items(),
            context=context,
            filter_mode=filter_mode,
            active_filters=active_filters or [],
            search_query=search_query,
            taxon=taxon,
        )
    except sqlite3.OperationalError as e:
        return handle_sql_error(e)








@app.route('/bin', methods=['GET', 'POST'])
def show_bins():
    def normalize_taxon_case(taxon):
        # Handles cases like O__christensenellales → o__Christensenellales
        if not taxon or len(taxon) < 4:
            return taxon
        prefix = taxon[:3].lower()
        rest = taxon[3:]
        if len(rest) == 0:
            return prefix
        return prefix + rest[0].upper() + rest[1:]

    try:
        map_number = request.values.get('map_number')
        kegg_id = request.values.get('kegg_id')
        taxon = request.values.get('taxon')
        search_query = request.form.get('search_query', '') or request.values.get('search_query', '')
        gtdb_filter = session.get('gtdb_filter', False)
        sort_filter = session.get('selected_sort_option', False)
        bin_name_sort_sql_command = " ORDER BY bin_number ASC"
        completeness_sort_sql_command = " ORDER BY bin.completeness DESC"
        contamination_sort_sql_command = " ORDER BY bin.contamination DESC"
        search_fields = request.form.getlist('search_fields') or ['sample_name', 'bin_name']

        context = "Display of all bins"

        conn = get_db_connection()
        cur = conn.cursor()

        bin_columns = [
            'bin.id as bin_id',
            'bin.bin_name',
            'bin.completeness',
            'bin.contamination',
            'sample.sample_name',
            'CAST(SUBSTR(bin.bin_name, INSTR(bin.bin_name, \'.\') + 1) AS INTEGER) AS bin_number'
        ]
        taxonomy_columns = [
            'taxonomy."_kingdom_" as kingdom',
            'taxonomy."_phylum_" as phylum',
            'taxonomy."_class_" as class',
            'taxonomy."_order_" as "order"',
            'taxonomy."_family_" as family',
            'taxonomy."_genus_" as genus',
            'taxonomy."_species_" as species'
        ]
        join_query = "LEFT JOIN taxonomy ON bin.taxonomic_id = taxonomy.id"
        query = f"SELECT DISTINCT {', '.join(bin_columns + taxonomy_columns)} FROM bin {join_query}"
        query += " JOIN sample ON sample.id = bin.sample_id "

        conditions = []
        params = []

        # Map, KEGG, Taxon filters
        if map_number:
            query += " JOIN bin_map_kegg bmk ON bin.id = bmk.bin_id"
            query += " JOIN map_kegg mk ON bmk.map_kegg_id = mk.id"
            query += " JOIN map m ON mk.map_id = m.id"
            conditions.append("m.map_number = ?")
            params.append(map_number)
            cur.execute("SELECT pathway_name FROM map WHERE map_number = ?", (map_number,))
            map_row = cur.fetchone()
            pathway_name = map_row[0] if map_row else None
            if pathway_name:
                context = f"Display of bins for Pathway: {pathway_name}"
            else:
                context = f"Display of bins for Map number: {map_number}"
        elif kegg_id:
            query += " JOIN bin_map_kegg bmk ON bin.id = bmk.bin_id"
            query += " JOIN map_kegg mk ON bmk.map_kegg_id = mk.id"
            query += " JOIN kegg k ON mk.kegg_id = k.id"
            conditions.append("k.ko_id = ?")
            params.append(kegg_id)
            context = f"Display of bins for KEGG ID: <strong>{kegg_id}</strong>"
        elif taxon:
            normalized_taxon = normalize_taxon_case(taxon)
            prefix = normalized_taxon[:3].lower()
            rank_map = {
                'd__': 'kingdom',
                'p__': 'phylum',
                'c__': 'class',
                'o__': 'order',
                'f__': 'family',
                'g__': 'genus',
                's__': 'species'
            }
            rank = rank_map.get(prefix)
            cleaned_taxon = normalized_taxon.strip()
            if rank in ['genus', 'species']:
                if cleaned_taxon.lower() in ['g__', 's__']:
                    # Match unclassified: g__/s__, empty, or NULL
                    conditions.append(
                        f"(LOWER(taxonomy.\"_{rank}_\") = ? OR taxonomy.\"_{rank}_\" = '' OR taxonomy.\"_{rank}_\" IS NULL)"
                    )
                    params.append(cleaned_taxon.lower())
                else:
                    conditions.append(f'taxonomy."_{rank}_" LIKE ?')
                    params.append(cleaned_taxon + '%')
            elif rank:
                conditions.append(f'taxonomy."_{rank}_" = ?')
                params.append(normalized_taxon)
            else:
                # fallback: match on any rank as before (just in case)
                conditions.append(
                    "(taxonomy.\"_kingdom_\" = ? OR taxonomy.\"_phylum_\" = ? OR taxonomy.\"_class_\" = ? "
                    "OR taxonomy.\"_order_\" = ? OR taxonomy.\"_family_\" = ? OR taxonomy.\"_genus_\" = ? OR taxonomy.\"_species_\" = ?)"
                )
                params.extend([normalized_taxon]*7)
            context = f"Display of bins for taxonomy entry: <strong>{taxon}</strong>"

        if gtdb_filter:
            conditions.append("(bin.completeness - 5 * bin.contamination > 50)")

        # -- Main Search Logic --
        if search_query:
            search_column_map = {
                'sample_name': 'sample.sample_name',
                'bin_name': 'bin.bin_name',
                'kingdom': 'taxonomy.\"_kingdom_\"',
                'phylum': 'taxonomy.\"_phylum_\"',
                'class': 'taxonomy.\"_class_\"',
                'order': 'taxonomy.\"_order_\"',
                'family': 'taxonomy.\"_family_\"',
                'genus': 'taxonomy.\"_genus_\"',
                'species': 'taxonomy.\"_species_\"'
            }
            search_subclauses = []
            search_pattern = f"%{search_query}%"
            # Check for "unclassified" (case-insensitive, any fragment)
            is_unclassified = "unclassified" in search_query.lower()
            for field in search_fields:
                dbcol = search_column_map.get(field)
                if dbcol:
                    # Taxonomic search for unclassified
                    if is_unclassified and field in ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']:
                        prefix = field[0].lower() + '__'  # e.g., 'g__' for genus
                        search_subclauses.append(
                            f"({dbcol} IS NULL OR {dbcol} = '' OR LOWER({dbcol}) = ?)"
                        )
                        params.append(prefix)
                    else:
                        search_subclauses.append(f"{dbcol} LIKE ?")
                        params.append(search_pattern)
            if search_subclauses:
                search_condition = '(' + ' OR '.join(search_subclauses) + ')'
                conditions.append(search_condition)
            context += f" with search pattern: <strong>{search_query}</strong>"


        if conditions:
            query += " WHERE " + " AND ".join(conditions)

        # Sorting
        if sort_filter == "option1":
            query += bin_name_sort_sql_command
        elif sort_filter == "option2":
            query += completeness_sort_sql_command
        else:
            query += contamination_sort_sql_command

        cur.execute(query, params)

        rows = cur.fetchall()
        cur.close()
        conn.close()

        display_column_labels = [col.split(' as ')[1] if ' as ' in col else col.split('.')[1] for col in bin_columns]
        bins = []
        sample_names = set()
        for row in rows:
            bin_data = dict(zip(display_column_labels + ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'], row))
            bins.append(bin_data)
            sample_names.add(bin_data['sample_name'])
        sample_names = sorted(sample_names)

        return render_template(
            'bin.html',
            bins=bins,
            columns=display_column_labels,
            context=context,
            sample_names=sample_names,
            search_fields=search_fields,
            search_query=search_query,
            taxon=taxon  # persist in form
        )
    except sqlite3.OperationalError as e:
        return handle_sql_error(e)




@app.route('/toggle_gtdb_filter', methods=['POST'])
def toggle_gtdb_filter():
    try:
        session['gtdb_filter'] = request.form.get('gtdb_filter') == 'on'

        # Rediriger en conservant les paramètres actuels
        map_number = request.form.get('map_number')
        kegg_id = request.form.get('kegg_id')
        taxon = request.form.get('taxon')  # Utilisez request.form au lieu de request.args

        if map_number:
            return redirect(url_for('show_bins', map_number=map_number))
        elif kegg_id:
            return redirect(url_for('show_bins', kegg_id=kegg_id))
        elif taxon:
            return redirect(url_for('show_bins', taxon=taxon))
        else:
            return redirect(url_for('show_bins'))
    except sqlite3.OperationalError as e:
        return handle_sql_error(e)

@app.route('/set_sort_option', methods=['POST'])
def set_sort_option():
    try:
        selected_option = request.form.get('sort_option')
        session['selected_sort_option'] = selected_option
        return redirect(url_for('show_bins'))
    except sqlite3.OperationalError as e:
        return handle_sql_error(e)

@app.route('/')
def home():
    try:
        if 'gtdb_filter' not in session:
            session['gtdb_filter'] = False  # Valeur par défaut
        if 'selected_sort_option' not in session:
            session['selected_sort_option'] = 'option1'  # ou toute autre valeur par défaut
        return render_template('index.html', content="Testing")
    except Exception as e:
        return handle_sql_error(e)


def get_default_db_path():
    try:
        resource = files('data').joinpath('pankegg.db')
        return str(resource)
    except Exception:
        return os.path.join('data', 'pankegg.db')

@click.command()
@click.option('--database', "--d", default=get_default_db_path(), help='Path to the SQLite database file.')
def start_server(database):
    global db_path
    db_path = database
    print(db_path)
    app.run(host='0.0.0.0', port=5000, debug=True)


if __name__ == '__main__':
    # app.run(debug=True)
    start_server()
