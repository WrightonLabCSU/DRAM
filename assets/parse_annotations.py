import sqlite3
import sys

def extract_query_ids(db_path, ko_list):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    ko_terms = ko_list.split(';')
    query = f"SELECT DISTINCT sample, query_id FROM annotations WHERE gene_id IN ({', '.join('?' for _ in ko_terms)})"
    cursor.execute(query, ko_terms)
    results = cursor.fetchall()
    conn.close()
    return results

def main():
    if len(sys.argv) != 4:
        print("Usage: python parse_annotations.py <db_path> <ko_list> <output_file>")
        sys.exit(1)

    db_path, ko_list, output_file = sys.argv[1:]
    results = extract_query_ids(db_path, ko_list)

    with open(output_file, 'w') as file:
        for sample, query_id in results:
            file.write(f"{sample}\t{query_id}\n")

    print(f"Extracted {len(results)} entries written to {output_file}")

if __name__ == '__main__':
    main()