import os
import subprocess
from os import path, mkdir
import argparse
import shutil
from shutil import move, rmtree
from glob import glob
import requests

DEFAULT_THREADS = 1

# Define database URLs
database_urls = {
    "dbcan": {
        "hmm_url": "https://bcb.unl.edu/dbCAN2/download/Databases/V12/dbCAN-HMMdb-V12.txt",
        "family_url": "https://bcb.unl.edu/dbCAN2/download/Databases/V12/CAZyDB.08062022.fam.subfam.ec.txt",
        "subfamily_url": "https://bcb.unl.edu/dbCAN2/download/Databases/V12/CAZyDB.08062022.fam-activities.txt"
    },
    "cant_hyd": {
        "base_url": "https://api.github.com/repos/dgittins/CANT-HYD-HydrocarbonBiodegradation/contents/HMMs/concatenated%20HMMs",
    },
    "camper": {
        "hmm_url": "https://raw.githubusercontent.com/WrightonLabCSU/CAMPER/main/CAMPER.hmm",
        "faa_url": "https://raw.githubusercontent.com/WrightonLabCSU/CAMPER/main/CAMPER_blast.faa",
        "mmseq_scores_url": "https://raw.githubusercontent.com/WrightonLabCSU/CAMPER/main/CAMPER_blast_scores.tsv",
        "hmm_scores_url": "https://raw.githubusercontent.com/WrightonLabCSU/CAMPER/main/CAMPER_hmm_scores.tsv"
    },
    "vogdb": {
        "hmm_url": "https://fileshare.lisc.univie.ac.at/vog/latest/vog.hmm.tar.gz",
        "annotations_url": "https://fileshare.lisc.univie.ac.at/vog/latest/vog.annotations.tsv.gz"
    },
    "pfam": {
        "hmm_url": "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz",
        "mmseq_url": "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.full.gz"
    },
    "merops": {
        "mmseqs_url": "ftp://ftp.ebi.ac.uk/pub/databases/merops/current_release/pepunit.lib"
    },
    "kofam": {
        "hmm_url": "ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz",
        "ko_list_url": "ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz"
    },
    "fegenie": {
        "base_url": "https://api.github.com/repos/Arkadiy-Garber/FeGenie/contents/hmms/iron/",
        "directories": [
            "iron_aquisition-heme_oxygenase",
            "iron_aquisition-heme_transport",
            "iron_aquisition-iron_transport",
            "iron_aquisition-siderophore_synthesis",
            "iron_aquisition-siderophore_transport",
            "iron_aquisition-siderophore_transport_potential",
            "iron_gene_regulation",
            "iron_oxidation",
            "iron_storage",
            "iron_reduction",
            "magnetosome_formation"
        ],
        "files": [
            "HMM-bitcutoffs.txt"
        ]
    },
    "sulfur": {
        "hmm_path": "/home/reedrich/Wrighton-Lab/Projects/DRAM-Main-Main-Project/DRAM2/DRAM-2-Nextflow/databases/sulfur/sulfur.hmm"
    },
    "methyl": {
        "faa_path": "/home/reedrich/Wrighton-Lab/Projects/DRAM-Main-Main-Project/DRAM2/DRAM-2-Nextflow/databases/methyl/methylotrophy.faa"
    }
}

# Define rename_mapping dictionaries for each database
rename_mapping_dbcan = {
    "hmm_url": "dbcan.hmm",
    "family_url": "dbcan.fam-activities.tsv",
    "subfamily_url": "dbcan.fam.subfam.ec.tsv"
}

rename_mapping_cant_hyd = {
    "hmm_file": "cant_hyd.hmm"
}

rename_mapping_camper = {
    "CAMPER.hmm": "hmm/camper.hmm",
    "CAMPER_blast.faa": "mmseqs/camper_blast.faa",
    "CAMPER_blast_scores.tsv": "mmseqs/camper_blast_scores.tsv",
    "CAMPER_hmm_scores.tsv": "hmm/camper_hmm_scores.tsv"
}

rename_mapping_vogdb = {
    "hmm_url": "vog_latest_hmms.hmm",
    "annotations_url": "vog_annotations_latest.tsv.gz"
}

rename_mapping_pfam = {
    "Pfam-A.hmm.dat.gz": "hmm/Pfam-A.hmm.dat",
    "Pfam-A.full.gz": "mmseqs/Pfam-A.full.gz"
}

rename_mapping_merops = {
    "pepunit.lib": "mmseqs/merops_peptidases_nr.faa",
}

rename_mapping_kofam = {
    "hmm_url": "kofam_profiles.hmm",
    "ko_list": "kofam_ko_list.tsv"
}

rename_mapping_fegenie = {
    "hmm_file": "fegenie.hmm",
    "HMM-bitcutoffs.txt": "fegenie_iron_cut_offs.txt",
}

rename_mapping_sulfur = {
    "hmm_file": "sulfur.hmm"
}

rename_mapping_methyl = {
    "hmm_file": "methyl.faa"
}

def download_cant_hyd_files(output_dir):
    base_url = database_urls["cant_hyd"]["base_url"]
    response = requests.get(base_url)
    
    if response.status_code == 200:
        files = response.json()
        for file in files:
            if file["type"] == "file" and file["name"] == "CANT-HYD.hmm":
                download_url = file["download_url"]
                file_name = "cant_hyd.hmm"  # Rename the file to cant_hyd.hmm
                download_file_path = path.join(output_dir, file_name)
                
                file_response = requests.get(download_url)
                if file_response.status_code == 200:
                    with open(download_file_path, "wb") as f:
                        f.write(file_response.content)
                    print(f"Downloaded file: {file_name}")
                else:
                    print(f"Failed to download file: {file_name}. Status code: {file_response.status_code}")
    else:
        print(f"Failed to fetch file list from GitHub. Status code: {response.status_code}")

    print("cant_hyd files downloaded successfully.")

def download_fegenie_files(output_dir):
    repo_url = "https://github.com/Arkadiy-Garber/FeGenie.git"
    fegenie_dir = path.join(output_dir, "FeGenie")

    # Remove the directory and its contents if it exists
    if path.exists(fegenie_dir):
        print(f"Removing contents of the directory {fegenie_dir}.")
        shutil.rmtree(fegenie_dir)

    # Clone the FeGenie repository
    git_clone_repo(repo_url, fegenie_dir)

    # Move the contents of the cloned repository to the output directory
    for item in os.listdir(fegenie_dir):
        item_path = path.join(fegenie_dir, item)
        destination_path = path.join(output_dir, item)

        # Move the item to the destination directory
        if path.exists(destination_path):
            if path.isdir(destination_path):
                shutil.rmtree(destination_path)
            else:
                os.remove(destination_path)
        move(item_path, output_dir)
        print(f"Moved {item} to {output_dir}")

    # Move and rename HMM-bitcutoffs.txt
    bitcutoffs_src = path.join(output_dir, "iron", "HMM-bitcutoffs.txt")
    bitcutoffs_dst = path.join(output_dir, "fegenie_iron_cut_offs.txt")
    if path.exists(bitcutoffs_src):
        move(bitcutoffs_src, bitcutoffs_dst)
        print(f"Moved and renamed HMM-bitcutoffs.txt to {bitcutoffs_dst}")

    # Remove the cloned repository directory
    rmtree(fegenie_dir)

    print("FeGenie files downloaded successfully.")


def download_fegenie_file(base_url, file_name, output_dir):
    """
    Download a file from the FeGenie repository.
    
    Args:
        base_url (str): The base URL of the FeGenie repository.
        file_name (str): The name of the file to download.
        output_dir (str): The directory where the downloaded file will be saved.
    """
    download_url = f"{base_url}/{file_name}"
    download_file_path = path.join(output_dir, file_name)

    # Send a GET request to the download URL
    response = requests.get(download_url)

    if response.status_code == 200:
        # Write the content of the response to the file
        with open(download_file_path, "wb") as file:
            file.write(response.content)
        print(f"Downloaded file: {file_name}")
    else:
        print(f"Failed to download file: {file_name}. Status code: {response.status_code}")

def download_databases(output_dir, databases):
    for database in databases:
        if database in database_urls:
            urls = database_urls[database]
            database_output_dir = path.join(output_dir, database)
            if not path.exists(database_output_dir):
                mkdir(database_output_dir)

            if database == "fegenie":
                download_fegenie_files(database_output_dir)
            elif database == "sulfur" or database == "methyl":
                continue  # Skip downloading for sulfur and methyl
            elif database == "dbcan":
                for url_key, url_value in urls.items():
                    download_url = url_value  # Get the download URL from the dictionary
                    download_file_name = path.basename(download_url)  # Extract the file name from the URL
                    download_file_path = path.join(database_output_dir, download_file_name)

                    # Download the file using curl
                    curl_command = f"curl -L -o {download_file_path} {download_url}"
                    subprocess.run(curl_command, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

                    print(f"Downloaded file: {download_file_name}")
            elif database == "cant_hyd":
                download_cant_hyd_files(database_output_dir)
            elif database == "kofam":
                for url_key, url_value in urls.items():
                    download_url = url_value  # Get the download URL from the dictionary
                    download_file_name = path.basename(download_url)  # Extract the file name from the URL
                    download_file_path = path.join(database_output_dir, download_file_name)

                    # Download the file using curl
                    curl_command = f"curl -L -o {download_file_path} {download_url}"
                    subprocess.run(curl_command, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

                    print(f"Downloaded file: {download_file_name}")
            else:
                for url_key, url_value in urls.items():
                    download_url = urls[url_key]  # Get the download URL from the dictionary
                    download_file_name = path.basename(download_url)  # Extract the file name from the URL
                    download_file_path = path.join(database_output_dir, download_file_name)

                    # Download the file using curl
                    curl_command = f"curl -L -o {download_file_path} {download_url}"
                    subprocess.run(curl_command, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

                    print(f"Downloaded file: {download_file_name}")

            print(f"All files downloaded for {database} database.")
        else:
            print(f"URLs for {database} database are not defined. Skipping.")

def process_vogdb_database(database_file, output_dir, num_threads=1):
    print("Processing VogDB database...")
    if path.exists(database_file):
        # Extract the downloaded tar.gz file
        tar_command = f"tar -xzf {database_file} -C {output_dir}"
        subprocess.run(tar_command, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        # Locate the .gz file in the extracted directory (if any)
        hmm_files = glob(path.join(output_dir, "hmm", "*.hmm.gz"))
        print("Extracted .hmm.gz files:", hmm_files)

        # Decompress each .hmm.gz file using pigz
        for hmm_file in hmm_files:
            decompress_command = f"pigz -d -p {num_threads} {hmm_file}"
            subprocess.run(decompress_command, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            print(f"Decompressed file: {hmm_file}")

        print("VogDB processed successfully.")
    else:
        print("VogDB file not found. Skipping processing.")

def process_camper_database(output_dir):
    print("Processing Camper database...")

    for original_file, new_file in rename_mapping_camper.items():
        original_path = path.join(output_dir, original_file)
        new_path = path.join(output_dir, new_file)
        new_dir = path.dirname(new_path)

        # Create target directory if it doesn't exist
        if not path.exists(new_dir):
            os.makedirs(new_dir)

        if path.exists(original_path):
            move(original_path, new_path)
            print(f"File {original_file} renamed and moved to {new_path}")
        else:
            print(f"File {original_file} not found. Skipping.")

    # Index the HMM file
    hmm_file = path.join(output_dir, "hmm", "camper.hmm")
    try:
        index_hmm(path.join(output_dir, "hmm"), hmm_file)
    except subprocess.CalledProcessError as e:
        print(f"Error running hmmpress on {hmm_file}: {e}")

    # Rename MMseqs files to follow the desired naming scheme
    mmseqs_dir = path.join(output_dir, "mmseqs")
    mmseqs_files = {
        "camper_blast_scores.tsv": "camper_scores.tsv"
    }
    for old_name, new_name in mmseqs_files.items():
        old_path = path.join(mmseqs_dir, old_name)
        new_path = path.join(mmseqs_dir, new_name)
        if path.exists(old_path):
            move(old_path, new_path)
            print(f"File {old_name} renamed and moved to {new_path}")
        else:
            print(f"File {old_name} not found. Skipping.")

    print("Camper processed successfully.")

def process_pfam_database(output_dir, num_threads=1):
    print("Processing Pfam database...")
    
    for original_file, new_file in rename_mapping_pfam.items():
        original_path = path.join(output_dir, original_file)
        new_path = path.join(output_dir, new_file)
        new_dir = path.dirname(new_path)

        # Create target directory if it doesn't exist
        if not path.exists(new_dir):
            os.makedirs(new_dir)

        if path.exists(original_path):
            if original_file.endswith(".gz"):
                decompress_command = f"pigz -d -p {num_threads} {original_path}"
                subprocess.run(decompress_command, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                original_path = original_path[:-3]  # Remove .gz extension after decompression
            move(original_path, new_path)
            print(f"File {original_file} decompressed (if needed), renamed and moved to {new_path}")
        else:
            print(f"File {original_file} not found. Skipping.")

    print("Pfam processed successfully.")

def process_merops_database(output_dir, num_threads=1):
    print("Processing Merops database...")
    
    for original_file, new_file in rename_mapping_merops.items():
        original_path = path.join(output_dir, original_file)
        new_path = path.join(output_dir, new_file)
        new_dir = path.dirname(new_path)

        # Create target directory if it doesn't exist
        if not path.exists(new_dir):
            os.makedirs(new_dir)

        if path.exists(original_path):
            move(original_path, new_path)
            print(f"File {original_file} renamed and moved to {new_path}")
        else:
            print(f"File {original_file} not found. Skipping.")

    print("Merops processed successfully.")

    # Create the MMseqs2 database
    create_mmseqs(new_path, new_path + "_db", num_threads)

def process_fegenie_database(output_dir):
    print("Processing FeGenie database...")
    
    # Download and prepare FeGenie files
    download_fegenie_files(output_dir)

    # Concatenate all .hmm files into a single file
    hmm_output_file = path.join(output_dir, "fegenie.hmm")
    concatenate_hmm_files(output_dir, hmm_output_file)

    # Sanitize the concatenated HMM file
    sanitized_hmm_file = path.join(output_dir, "sanitized_fegenie.hmm")
    sanitize_hmm_file(hmm_output_file, sanitized_hmm_file)

    # Move the sanitized HMM file to replace the original one
    final_hmm_file = path.join(output_dir, "fegenie.hmm")
    os.remove(hmm_output_file)
    move(sanitized_hmm_file, final_hmm_file)

    # Index the sanitized HMM file
    try:
        index_hmm(output_dir, final_hmm_file)
    except subprocess.CalledProcessError as e:
        print(f"Error running hmmpress: {e}")
        print("Please check the format of the sanitized concatenated HMM file.")

    # Remove the directories after concatenation
    directories = database_urls["fegenie"]["directories"]
    for directory in directories:
        dir_path = path.join(output_dir, directory)
        if path.exists(dir_path):
            rmtree(dir_path)
            print(f"Removed directory: {dir_path}")

    # Delete all files except fegenie.hmm, fegenie.hmm index files, and fegenie_iron_cut_offs.txt
    for file in os.listdir(output_dir):
        if not (file.startswith("fegenie.hmm") or file == "fegenie_iron_cut_offs.txt"):
            file_path = path.join(output_dir, file)
            if path.isdir(file_path):
                rmtree(file_path)
            else:
                os.remove(file_path)
            print(f"Removed file: {file_path}")

    print("FeGenie processed successfully.")
 


def process_sulfur_database(output_dir):
    print("Processing Sulfur database...")
    sulfur_hmm_path = database_urls["sulfur"]["hmm_path"]
    
    # Create the sulfur output directory if it doesn't exist
    if not path.exists(output_dir):
        mkdir(output_dir)
    
    # Move the sulfur HMM file to the sulfur output directory
    destination_path = path.join(output_dir, path.basename(sulfur_hmm_path))
    shutil.copyfile(sulfur_hmm_path, destination_path)
    print(f"Copied sulfur HMM file to: {destination_path}")

    # Index the HMM file
    try:
        index_hmm(output_dir, output_dir)
    except subprocess.CalledProcessError as e:
        print(f"Error running hmmpress: {e}")
        print("Please check the format of the sulfur HMM file.")
    
    print("Sulfur processed successfully.")

def process_methyl_database(output_dir, num_threads=1):
    print("Processing Methyl database...")
    
    # Paths for methyl database files
    methyl_faa_path = database_urls["methyl"]["faa_path"]
    methyl_faa_dest = path.join(output_dir, "methyl.faa")

    # Create the output directory if it doesn't exist
    if not path.exists(output_dir):
        mkdir(output_dir)

    # Copy the methyl .faa file to the output directory
    shutil.copyfile(methyl_faa_path, methyl_faa_dest)
    print(f"Copied methyl .faa file to: {methyl_faa_dest}")

    print("Methyl processed successfully.")

def process_cant_hyd_database(output_dir):
    print("Processing cant_hyd database...")
    repo_url = "https://github.com/dgittins/CANT-HYD-HydrocarbonBiodegradation.git"
    cant_hyd_dir = path.join(output_dir, "cant_hyd_repo")
    hmm_dir = path.join(output_dir, "hmm")

    # Remove the directory and its contents if it exists
    if path.exists(cant_hyd_dir):
        print(f"Removing contents of the directory {cant_hyd_dir}.")
        shutil.rmtree(cant_hyd_dir)

    # Clone the cant_hyd repository
    git_clone_repo(repo_url, cant_hyd_dir)

    # Create the hmm subdirectory if it doesn't exist
    if not path.exists(hmm_dir):
        os.makedirs(hmm_dir)

    # Move the concatenated HMM file to the hmm subdirectory
    concatenated_hmm_path = path.join(cant_hyd_dir, "HMMs", "concatenated HMMs", "CANT-HYD.hmm")
    final_hmm_path = path.join(hmm_dir, "cant_hyd.hmm")
    if path.exists(concatenated_hmm_path):
        move(concatenated_hmm_path, final_hmm_path)
        print(f"Moved CANT-HYD.hmm to {final_hmm_path}")
    else:
        print(f"Concatenated HMM file not found in the repository.")

    # Remove the cloned repository directory
    rmtree(cant_hyd_dir)

    # Index the HMM file
    try:
        index_hmm(hmm_dir, final_hmm_path)
    except subprocess.CalledProcessError as e:
        print(f"Error running hmmpress: {e}")
        print("Please check the format of the cant_hyd HMM file.")

    # Ensure no extra cant_hyd.hmm files persist outside the hmm subdirectory
    if path.exists(path.join(output_dir, "cant_hyd.hmm")):
        os.remove(path.join(output_dir, "cant_hyd.hmm"))
        print("Removed extraneous cant_hyd.hmm from the output directory.")

    print("cant_hyd processed successfully.")

def process_dbcan_database(output_dir):
    print("Processing dbcan database...")
    rename_mapping = {
        "dbCAN-HMMdb-V12.txt": "dbcan.hmm",
        "CAZyDB.08062022.fam.subfam.ec.txt": "dbcan.fam.subfam.ec.tsv",
        "CAZyDB.08062022.fam-activities.txt": "dbcan.fam-activities.tsv"
    }

    for original_file, new_file in rename_mapping.items():
        original_path = path.join(output_dir, original_file)
        new_path = path.join(output_dir, new_file)
        new_dir = path.dirname(new_path)

        # Create target directory if it doesn't exist
        if not path.exists(new_dir):
            os.makedirs(new_dir)

        if path.exists(original_path):
            move(original_path, new_path)
            print(f"File {original_file} renamed and moved to {new_path}")
        else:
            print(f"File {original_file} not found. Skipping.")

    # Index the HMM file
    try:
        index_hmm(output_dir, path.join(output_dir, "dbcan.hmm"))
    except subprocess.CalledProcessError as e:
        print(f"Error running hmmpress: {e}")
        print("Please check the format of the dbcan HMM file.")

    print("dbcan processed successfully.")

def process_kofam_database(output_dir, num_threads=1):
    print("Processing kofam database...")
    hmm_dir = path.join(output_dir, "hmm")

    # Create the hmm subdirectory if it doesn't exist
    if not path.exists(hmm_dir):
        os.makedirs(hmm_dir)

    # Extract the profiles.tar.gz file
    profiles_tar_path = path.join(output_dir, "profiles.tar.gz")
    if path.exists(profiles_tar_path):
        tar_command = f"tar -xzf {profiles_tar_path} -C {hmm_dir}"
        subprocess.run(tar_command, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        print(f"Extracted {profiles_tar_path} to {hmm_dir}")
    else:
        print("Profiles tar file not found. Skipping extraction.")

    # Decompress ko_list.gz and rename
    ko_list_path = path.join(output_dir, "ko_list.gz")
    if path.exists(ko_list_path):
        decompress_command = f"pigz -d -p {num_threads} {ko_list_path}"
        subprocess.run(decompress_command, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        ko_list_path = ko_list_path[:-3]  # Remove .gz extension after decompression
        move(ko_list_path, path.join(output_dir, "kofam_ko_list.tsv"))
        print(f"Decompressed and renamed ko_list to kofam_ko_list.tsv")
    else:
        print("ko_list.gz file not found. Skipping decompression.")

    # Concatenate all HMM files into a single file
    kofam_hmm_file = path.join(output_dir, "kofam_profiles.hmm")
    concatenate_hmm_files(hmm_dir, kofam_hmm_file)

    # Index the concatenated HMM file
    try:
        index_hmm(output_dir, kofam_hmm_file)
    except subprocess.CalledProcessError as e:
        print(f"Error running hmmpress on {kofam_hmm_file}: {e}")
        print("Please check the format of the kofam HMM files.")

    # Remove the hmm directory and profiles.tar.gz file
    if path.exists(hmm_dir):
        rmtree(hmm_dir)
        print(f"Removed directory: {hmm_dir}")
    if path.exists(profiles_tar_path):
        os.remove(profiles_tar_path)
        print(f"Removed file: {profiles_tar_path}")

    print("kofam processed successfully.")

def generate_readme(output_dir, databases_info):
    def human_readable_size(size, decimal_places=1):
        for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
            if size < 1024.0:
                return f"{size:.{decimal_places}f} {unit}"
            size /= 1024.0

    readme_content = "Database Information:\n\n"
    for database_name, database_info in databases_info.items():
        readme_content += f"Database: {database_name}\n"
        readme_content += "URLs used for download:\n"
        for url_key, url_value in database_info["urls"].items():
            readme_content += f"- {url_key}: {url_value}\n"
        readme_content += f"\nDate and Time Downloaded: {database_info['download_date']}\n\n"
        
        # List contents of the directory with sizes, including subdirectories
        database_dir = path.join(output_dir, database_name)
        if path.exists(database_dir) and path.isdir(database_dir):
            readme_content += f"{database_name}/\n"
            for root, dirs, files in os.walk(database_dir):
                for file_name in files:
                    file_path = path.join(root, file_name)
                    size = os.path.getsize(file_path)
                    relative_path = path.relpath(file_path, database_dir)
                    readme_content += f"{relative_path} ({human_readable_size(size)})\n"
            readme_content += "\n"

    readme_file_path = path.join(output_dir, "README.txt")
    with open(readme_file_path, "w") as readme_file:
        readme_file.write(readme_content)

    print(f"README.txt generated successfully at {readme_file_path}")

def index_hmm(hmm_dir, hmm_file):
    try:
        # Run HMMER to index the HMM file
        index_command = f"hmmpress -f {hmm_file}"
        result = subprocess.run(index_command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print(result.stdout.decode('utf-8'))
        print(f"Indexed HMM file: {path.basename(hmm_file)}")
    except subprocess.CalledProcessError as e:
        print(f"Error indexing HMM file {path.basename(hmm_file)}: {e}")
        print(e.stderr.decode('utf-8'))

def create_mmseqs(fasta_loc, output_loc, threads):
    """Takes a fasta file and makes a mmseqs2 database for use in blast searching and hmm searching with mmseqs2."""
    print(f"Creating MMseqs2 database for {fasta_loc}...")
    subprocess.run(["mmseqs", "createdb", fasta_loc, output_loc], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    tmp_dir = path.join(path.dirname(output_loc), "tmp")
    subprocess.run(["mmseqs", "createindex", output_loc, tmp_dir, "--threads", str(threads)], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    print(f"MMseqs2 database created at {output_loc}")
    
    # Remove the temporary directory
    if path.exists(tmp_dir):
        rmtree(tmp_dir)
        print(f"Removed temporary directory: {tmp_dir}")

def sanitize_hmm_file(input_file, output_file):
    with open(input_file, 'r') as in_file:
        lines = in_file.readlines()

    sanitized_lines = []
    current_entry = []
    entry_started = False
    entry_names = {}
    name_suffix_counter = {}

    for line in lines:
        line = line.strip()
        if line.startswith("HMMER3/f"):
            if entry_started:
                # Process current entry before starting a new one
                sanitized_entry = process_entry(current_entry, entry_names, name_suffix_counter)
                sanitized_lines.extend(sanitized_entry)
                current_entry = []
            entry_started = True
            current_entry.append(line)
        elif entry_started:
            current_entry.append(line)
            if line == "//":
                entry_started = False
                # Process the last entry
                sanitized_entry = process_entry(current_entry, entry_names, name_suffix_counter)
                sanitized_lines.extend(sanitized_entry)
                current_entry = []

    if entry_started:
        # Process the final entry if it wasn't already processed
        sanitized_entry = process_entry(current_entry, entry_names, name_suffix_counter)
        sanitized_lines.extend(sanitized_entry)

    with open(output_file, 'w') as out_file:
        out_file.write('\n'.join(sanitized_lines))

    print(f"Sanitized HMM file written to: {output_file}")

def process_entry(entry_lines, entry_names, name_suffix_counter):
    name_line_index = None
    acc_line_index = None
    entry_name = None

    for i, line in enumerate(entry_lines):
        if line.startswith("NAME"):
            name_line_index = i
            entry_name = line.split()[1]
        elif line.startswith("ACC"):
            acc_line_index = i

    if entry_name:
        if entry_name in entry_names:
            name_suffix_counter[entry_name] += 1
            unique_suffix = f"_{name_suffix_counter[entry_name]}"
            entry_name += unique_suffix
            entry_lines[name_line_index] = f"NAME  {entry_name}"
            if acc_line_index is not None:
                acc_line = entry_lines[acc_line_index]
                acc_line_parts = acc_line.split()
                if len(acc_line_parts) > 1:
                    acc_line_parts[1] += unique_suffix
                    entry_lines[acc_line_index] = " ".join(acc_line_parts)
        else:
            entry_names[entry_name] = True
            name_suffix_counter[entry_name] = 0

    return entry_lines

def concatenate_hmm_files(hmm_dir, output_file):
    hmm_files = []
    for root, _, files in os.walk(hmm_dir):
        for file in files:
            if file.endswith(".hmm"):
                hmm_files.append(path.join(root, file))

    hmm_files.sort()

    with open(output_file, 'w') as out_file:
        for hmm_file in hmm_files:
            with open(hmm_file, 'r') as in_file:
                lines = in_file.readlines()
                out_file.writelines(lines)
                if not lines[-1].strip() == "//":
                    out_file.write("\n//\n")
                else:
                    out_file.write("\n")

    print(f"Concatenated HMM files into: {output_file}")

def git_clone_repo(repo_url, output_dir):
    clone_command = f"git clone {repo_url} {output_dir}"
    subprocess.run(clone_command, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

def main():
    parser = argparse.ArgumentParser(description="Prepare databases for annotation")
    parser.add_argument("--output_dir", help="Output directory path", required=True)
    parser.add_argument("--databases", "-db", nargs='+', help="List of databases to update", required=True)
    parser.add_argument("--threads", type=int, default=DEFAULT_THREADS, help="Number of threads for processing")
    parser.add_argument("--verbose", action="store_true", help="Verbose mode")
    parser.add_argument("--download_date", help="Download date for databases")
    args = parser.parse_args()

    # Download databases
    download_databases(args.output_dir, args.databases)

    # Process databases
    for database in args.databases:
        database_output_dir = path.join(args.output_dir, database)
        if database == "vogdb":
            process_vogdb_database(path.join(database_output_dir, "vog.hmm.tar.gz"), database_output_dir, args.threads)
            concatenate_hmm_files(path.join(database_output_dir, "hmm"), path.join(database_output_dir, "vog_latest_hmms.hmm"))
            rmtree(path.join(database_output_dir, "hmm"))
            os.remove(path.join(database_output_dir, "vog.hmm.tar.gz"))
            index_hmm(database_output_dir, path.join(database_output_dir, "vog_latest_hmms.hmm"))
        elif database == "camper":
            process_camper_database(database_output_dir)
            create_mmseqs(path.join(database_output_dir, "mmseqs", "camper_blast.faa"), path.join(database_output_dir, "mmseqs", "camper.mmsdb"), args.threads)
        elif database == "pfam":
            process_pfam_database(database_output_dir, args.threads)
        elif database == "merops":
            process_merops_database(database_output_dir, args.threads)
        elif database == "fegenie":
            process_fegenie_database(database_output_dir)
        elif database == "sulfur":
            process_sulfur_database(database_output_dir)
        elif database == "methyl":
            process_methyl_database(database_output_dir)
            create_mmseqs(path.join(database_output_dir, "methyl.faa"), path.join(database_output_dir, "methyl"), args.threads)
        elif database == "dbcan":
            process_dbcan_database(database_output_dir)
        elif database == "cant_hyd":
            process_cant_hyd_database(database_output_dir)
        elif database == "kofam":
            process_kofam_database(database_output_dir, args.threads)

    # Generate README
    databases_info = {}
    for database in args.databases:
        if database in database_urls:
            databases_info[database] = {
                "urls": database_urls[database],
                "download_date": args.download_date
            }

    generate_readme(args.output_dir, databases_info)

if __name__ == "__main__":
    main()