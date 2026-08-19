#!/usr/bin/env python

import argparse
import os
from typing import Optional

import pandas as pd

from database_utils import (
    get_tools_info,
    insert_many_annotations_to_toolset_db,
    generate_toolset_hash,
    _build_dsn,
)

def _parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Load annotated variants into Postgres (partitioned)")

    # Input
    parser.add_argument(
        "-m", "--maf", type=str, required=True, help="Input MAF/TSV file with annotations"
    )
    parser.add_argument(
        "-st", "--sample_type", type=str, default="somatic",
        choices=["somatic", "germline"], help="Sample type"
    )

    # Tooling/config
    parser.add_argument(
        "-c", "--config", type=str, required=True,
        help="Path to CancerVar/InterVar (Annovar) config file"
    )
    parser.add_argument(
        "-f", "--funcotator", type=str, required=True,
        help="Path to Funcotator database folder"
    )
    parser.add_argument(
        "-t", "--technology", type=str, default="illumina",
        help="Sequencing technology label"
    )
    parser.add_argument(
        "-av", "--annovar_version", type=str, default="v0",
        help="Annovar version"
    )
    parser.add_argument(
        "-gb", "--genome_build", type=str, default="GRCh38",
        help="Genome build (e.g., GRCh38)"
    )
    parser.add_argument(
        "--tissue_type",
        type=str,
        default=None,
        help="Tumor tissue type (optional, used for database schema annotation)",
    )

    # Postgres connection
    parser.add_argument(
        "--dsn",
        type=str,
        default=None,
        help="PostgreSQL DSN connection string (e.g., postgresql://user:pass@host:5432/dbname)"
    )
    parser.add_argument("--pg-host", type=str, default=None, help="Postgres host (fallback if --dsn missing)")
    parser.add_argument("--pg-port", type=int, default=None, help="Postgres port (fallback if --dsn missing)")
    parser.add_argument("--pg-user", type=str, default=None, help="Postgres user (fallback if --dsn missing)")
    parser.add_argument("--pg-password", type=str, default=None, help="Postgres password (fallback if --dsn missing)")
    parser.add_argument("--pg-database", type=str, default=None, help="Postgres database name (fallback if --dsn missing)")

    # Target schema & batching
    parser.add_argument(
        "--schema", type=str, default=None,
        help="Override schema name (otherwise derived from toolset hash)"
    )
    parser.add_argument(
        "--batch-size", type=int, default=10_000,
        help="Batch size for batched UPSERTs"
    )

    return parser.parse_args()

def main():
    args = _parse_args()

    # Check if PostgreSQL credentials are provided
    if args.dsn or (args.pg_host and args.pg_user):

        if not os.path.exists(args.maf):
            raise ValueError(f"Maf file {args.maf} does not exist!")

        # Load MAF/TSV
        maf_df = pd.read_csv(args.maf, sep="\t", header=0)

        # Build tools_info → toolset hash → schema (if not overridden)
        tools_info = get_tools_info(
            annovar_config_file=args.config,
            funcotator_db_path=args.funcotator,
            technology=args.technology,
            annovar_version=args.annovar_version,
            germline=(args.sample_type == "germline"),
            genome_build=args.genome_build,
        )
        
        tools_info.append({
                "tool_info": {"name": "tissue_type", "version": "tissue_type"},
                "db_info": {"name": "tissue_type", "version": args.tissue_type if args.tissue_type else "unknown"},
            })
        # PostgreSQL credentials provided - insert annotations
        # Priority: 1) pg_dsn, 2) pg_host+port+user+password
        dsn_list = []
        
        # Try pg_dsn first if provided
        if args.dsn:
            dsn_list.append(("DSN", args.dsn))
        
        # Add host+port+user+password as fallback (or primary if no DSN)
        if args.pg_host and args.pg_user:
            try:
                fallback_dsn = _build_dsn(
                    explicit_dsn=None,
                    host=args.pg_host,
                    port=args.pg_port,
                    db=args.pg_database,
                    user=args.pg_user,
                    password=args.pg_password,
                )
                dsn_list.append(("host+port+user+password", fallback_dsn))
            except ValueError:
                pass
        
        # Try each DSN in order
        success = False
        for method, dsn in dsn_list:
            try:
                # Insert in batches with UPSERT into partitioned table
                toolset_id, schema_name = insert_many_annotations_to_toolset_db(
                    maf=maf_df,
                    tools_info=tools_info,
                    dsn=dsn,
                    schema_name=args.schema,      # if None → derived from hash
                    batch_size=args.batch_size,
                )

                print(f"[OK] Loaded {len(maf_df)} rows into schema '{schema_name}' using {method}.")
                print(f"   Toolset ID: {toolset_id[:12]}...")
                success = True
                break
            except Exception as e:
                if method == dsn_list[-1][0]:  # Last attempt
                    print(f"[WARNING] Failed to insert annotations into PostgreSQL: {e}")
                    print(f"   Annotations were not saved to database.")
                else:
                    print(f"[WARNING] Insertion failed using {method}, trying next method...")
                    continue
    else:
        print(f"[INFO] No PostgreSQL credentials provided (--dsn or --pg-host/--pg-user).")
        print(f"   Skipping database insertion. Annotations will not be cached for future use.")


if __name__ == "__main__":
    main()
