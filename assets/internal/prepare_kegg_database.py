from os import path, mkdir
from datetime import datetime
from shutil import move, rmtree
from glob import glob
import logging
import subprocess
from skbio import write as write_sequence, read as read_sequence
from collections import defaultdict
import gzip
import argparse
from pathlib import Path


LOGGER = logging.getLogger("database_processing.log")


def prepare_databases(
    output_dir,
    kegg_loc=None,
    gene_ko_link_loc=None,
    kegg_download_date=None,
    threads=10,
    verbose=True,
):
    # setup temp, logging, and db_handler
    output_dir = Path(output_dir)
    if not output_dir.exists():
        output_dir.mkdir(parents=True)
    temporary = output_dir / "database_files"
    temporary.mkdir()

    LOGGER.info("Database preparation started")
    LOGGER.info("Processing KEGG database")

    processed_locs = process_kegg(
        kegg_loc=kegg_loc,
        output_dir=output_dir,
        gene_ko_link_loc=gene_ko_link_loc,
        logger=LOGGER,
        threads=threads,
        verbose=verbose,
        download_date=kegg_download_date,
    )
    # Process databases

    for k, v in processed_locs.items():
        final_dest = output_dir / path.basename(v)
        if v != final_dest:
            for db_file in glob("%s*" % v):
                move(db_file, path.join(output_dir, path.basename(db_file)))
            v = path.join(output_dir, path.basename(v))
        # update_dram_forms the settings per OUTPUT fill, including the process_settings
        #  and database_settings, which are per input file.
        LOGGER.info(f"Moved {k} to final destination")


def process_kegg(
    kegg_loc,
    output_dir,
    logger,
    gene_ko_link_loc=None,
    download_date=None,
    threads=10,
    verbose=True,
):
    threads = threads or 10  # make sure cli option is >=1 and not None
    if download_date is None:
        download_date = get_iso_date()
    if gene_ko_link_loc is not None:
        # add KOs to end of header where KO is not already there
        kegg_mod_loc = path.join(output_dir, "kegg.mod.fa")
        write_sequence(
            generate_modified_kegg_fasta(kegg_loc, gene_ko_link_loc),
            format="fasta",
            into=kegg_mod_loc,
        )
    else:
        kegg_mod_loc = kegg_loc
    # make mmseqsdb from modified kegg fasta
    kegg_mmseqs_db = path.join(output_dir, "kegg.%s.mmsdb" % download_date)
    create_mmseqs(
        kegg_mod_loc,
        kegg_mmseqs_db,
        #logger,
        #create_index=True,
        threads=threads,
        #verbose=verbose,
    )
    LOGGER.info("KEGG database processed")
    return {"kegg": kegg_mmseqs_db}


def create_mmseqs(fasta_loc, output_loc, threads):
    """Takes a fasta file and makes a mmseqs2 database for use in blast searching and hmm searching with mmseqs2."""
    print(f"Creating MMseqs2 database for {fasta_loc}...")
    subprocess.run(
        ["mmseqs", "createdb", fasta_loc, output_loc],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    tmp_dir = path.join(path.dirname(output_loc), "tmp")
    print(f"Created MMseqs2 database for {fasta_loc}... Now creating index")
    subprocess.run(
        ["mmseqs", "createindex", output_loc, tmp_dir, "--threads", str(threads)],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    print(f"MMseqs2 database created at {output_loc}")

    # Remove the temporary directory
    if path.exists(tmp_dir):
        rmtree(tmp_dir)
        print(f"Removed temporary directory: {tmp_dir}")

# def create_mmseqs(fasta_loc, output_loc, threads):
#     """Takes a fasta file and makes a mmseqs2 database for use in blast searching and hmm searching with mmseqs2."""
#     print(f"Creating MMseqs2 database for {fasta_loc}...")
#     print( ["mmseqs", "createdb", fasta_loc, output_loc])
#     # subprocess.run(
#     #     ["mmseqs", "createdb", fasta_loc, output_loc],
#     #     check=True,
#     #     stdout=subprocess.DEVNULL,
#     #     stderr=subprocess.DEVNULL,
#     # )
#     tmp_dir = path.join(path.dirname(output_loc), "tmp")
#     print(["mmseqs", "createindex", output_loc, tmp_dir, "--threads", str(threads)])
#     # subprocess.run(
#     #     ["mmseqs", "createindex", output_loc, tmp_dir, "--threads", str(threads)],
#     #     check=True,
#     #     stdout=subprocess.DEVNULL,
#     #     stderr=subprocess.DEVNULL,
#     # )
#     print(f"MMseqs2 database created at {output_loc}")
#
#     # Remove the temporary directory
#     # if path.exists(tmp_dir):
#     #     rmtree(tmp_dir)
#     #     print(f"Removed temporary directory: {tmp_dir}")

def generate_modified_kegg_fasta(kegg_fasta, gene_ko_link_loc=None):
    """
    Takes kegg fasta file and gene ko link file, adds kos not already in headers to headers
    Whish I knew about this, oh well I may split this out.
    """
    genes_ko_dict = defaultdict(list)
    if gene_ko_link_loc is not None:
        if gene_ko_link_loc.endswith(".gz"):
            gene_ko_link_fh = gzip.open(gene_ko_link_loc, "rt")
        else:
            gene_ko_link_fh = open(gene_ko_link_loc)
        for line in gene_ko_link_fh:
            gene, ko = line.strip().split()
            genes_ko_dict[gene].append(remove_prefix(ko, "ko:"))
    for seq in read_sequence(kegg_fasta, format="fasta"):
        new_description = seq.metadata["description"]
        for ko in genes_ko_dict[seq.metadata["id"]]:
            if ko not in new_description:
                new_description += "; %s" % ko
        seq.metadata["description"] = new_description
        yield seq


def remove_prefix(text, prefix):
    if text.startswith(prefix):
        return text[len(prefix) :]
    return text  # or whatever


def get_iso_date():
    return datetime.today().strftime("%Y%m%d")


def main():
    parser = argparse.ArgumentParser(description="Prepare KEGG database")
    parser.add_argument("--kegg_loc", type=str, help="Path to the KEGG fasta file")
    parser.add_argument("--output_dir", type=str, help="Path to the output directory")
    parser.add_argument(
        "--gene_ko_link_loc", type=str, help="Path to the gene KO link file"
    )
    parser.add_argument("--download_date", type=str, help="Date of the KEGG download")
    parser.add_argument("--threads", type=int, help="Number of threads to use", default=10)
    args = parser.parse_args()

    prepare_databases(
        kegg_loc=args.kegg_loc,
        output_dir=args.output_dir,
        gene_ko_link_loc=args.gene_ko_link_loc,
        kegg_download_date=args.download_date,
        threads=args.threads,
        verbose=True,
    )


if __name__ == "__main__":
    main()
