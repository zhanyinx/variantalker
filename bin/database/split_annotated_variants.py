#!/usr/bin/env python
import argparse
import os

import pandas as pd
import psycopg

from database_utils import (
    get_tools_info,
    _build_dsn,
    get_annotation_from_toolset_db,
    get_or_create_schema_for_tools,
    get_column_order,
)


def _parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-v",
        "--vcf",
        type=str,
        default=None,
        required=True,
        help="vcf file",
    )
    parser.add_argument(
        "-f",
        "--funcotator",
        type=str,
        default=None,
        required=True,
        help="Path to funcotator database folder.",
    )
    parser.add_argument(
        "-c",
        "--config",
        type=str,
        default=None,
        required=True,
        help="Path to cancervar or intervar config file.",
    )
    parser.add_argument(
        "-d",
        "--db_dir",
        type=str,
        default=None,
        help="Path to sqlites3 databases folder.",
    )
    parser.add_argument(
        "-st", "--sample_type", type=str, default="somatic",
        choices=["somatic", "germline"], help="Sample type"
    )
    parser.add_argument(
        "-gb", "--genome_build", type=str, default="GRCh38", help="Genome build"
    )
    parser.add_argument(
        "-av", "--annovar_version", type=str, default="v0", help="Annovar version"
    )
    parser.add_argument(
        "-t",
        "--technology",
        type=str,
        default="Illumina",
        help="Technology used for sequencing.",
    )
    parser.add_argument(
        "--dsn",
        type=str,
        default=None,
        help="PostgreSQL DSN connection string (e.g., postgresql://user:pass@host:5432/dbname)",
    )
    parser.add_argument(
        "--pg-host",
        type=str,
        default=None,
        help="PostgreSQL host (fallback if --dsn not provided)",
    )
    parser.add_argument(
        "--pg-port",
        type=int,
        default=None,
        help="PostgreSQL port (fallback if --dsn not provided)",
    )
    parser.add_argument(
        "--pg-user",
        type=str,
        default=None,
        help="PostgreSQL user (fallback if --dsn not provided)",
    )
    parser.add_argument(
        "--pg-password",
        type=str,
        default=None,
        help="PostgreSQL password (fallback if --dsn not provided)",
    )
    parser.add_argument(
        "--pg-database",
        type=str,
        default=None,
        help="PostgreSQL database name (fallback if --dsn not provided)",
    )
    parser.add_argument(
        "--tissue_type",
        type=str,
        default=None,
        help="Tumor tissue type (optional, used for database schema annotation)",
    )

    return parser.parse_args()


def _safe_divide(num, den):
    den = float(den)
    if den == 0:
        return 0.0
    return float(num) / den


def _parse_info_field(info_str):
    info_dict = {}
    for item in info_str.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            info_dict[k] = v
    return info_dict


def _extract_ad_counts(format_str, sample_str):
    """
    Return (ref_count, alt_count) from FORMAT/sample columns if AD is present,
    otherwise return (None, None).
    """
    fmt_keys = format_str.split(":")
    sample_vals = sample_str.split(":")

    # Protect against malformed rows with mismatched FORMAT/sample lengths
    sample_map = dict(zip(fmt_keys, sample_vals))
    ad = sample_map.get("AD")
    if not ad:
        return None, None

    parts = ad.split(",")
    if len(parts) < 2:
        return None, None

    try:
        ref_count = int(parts[0])
        alt_count = int(parts[1])
    except ValueError:
        return None, None

    return ref_count, alt_count


def main():
    """Generate raw report in csv format."""
    args = _parse_args()

    if not os.path.exists(args.vcf):
        raise ValueError(f"VCF file {args.vcf} does not exist!")

    if not os.path.exists(args.funcotator):
        raise ValueError(
            f"Funcotator database folder {args.funcotator} does not exist!"
        )

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

    # Database setup
    use_database = False
    conn = None
    cur = None
    schema_name = None

    if args.dsn or (args.pg_host and args.pg_user):
        dsn_list = []

        if args.dsn:
            dsn_list.append(("DSN", args.dsn))

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

        for i, (method, dsn) in enumerate(dsn_list):
            try:
                print(f"[INFO] Attempting to connect to PostgreSQL database using {method}...")
                print(f"[DEBUG] DSN: {dsn}")
                print(
                    f"[DEBUG] Host: {args.pg_host}, Port: {args.pg_port}, "
                    f"User: {args.pg_user}, Database: {args.pg_database}"
                )
                print(method, dsn)

                toolset_id, schema_name = get_or_create_schema_for_tools(dsn, tools_info)
                print(f"[OK] Using schema '{schema_name}' for toolset ID {toolset_id}.")

                conn = psycopg.connect(dsn)
                cur = conn.cursor()
                use_database = True

                print(f"[OK] Connected to PostgreSQL database using {method}.")
                print("   Variants will be split based on database lookup.")
                print(f"   Schema: {schema_name}")
                break

            except Exception as e:
                print(e)
                is_last_attempt = (i == len(dsn_list) - 1)
                if is_last_attempt:
                    print(f"[WARNING] Failed to connect to PostgreSQL database: {e}")
                    print(
                        "   Proceeding without database lookup - all variants will be annotated."
                    )
                    use_database = False
                else:
                    print(f"[WARNING] Connection failed using {method}, trying next method...")
                    continue
    else:
        print("[INFO] No PostgreSQL credentials provided (--dsn or --pg-host/--pg-user).")
        print("   Skipping database lookup - all variants will be annotated.")

    skip_rows = []
    normal_sample_name = None
    tumor_sample_name = None
    sample_name = None
    tumor_id = None
    normal_id = None

    output_full_file = args.vcf.replace(".vcf", "_noskip.vcf")
    output_skip_file = args.vcf.replace(".vcf", "_skip.maf")

    try:
        with open(args.vcf, "r", encoding="utf-8") as f, open(
            output_full_file, "w", encoding="utf-8"
        ) as output_full:

            maf_header = []

            for line in f:
                if line.startswith("##"):
                    output_full.write(line)
                    maf_header.append(
                        line if line.startswith("## ")
                        else line.replace("##", "## ", 1)
                    )

                    if line.startswith("##normal_sample"):
                        normal_sample_name = line.split("=", 1)[1].rstrip("\n")
                    if line.startswith("##tumor_sample"):
                        tumor_sample_name = line.split("=", 1)[1].rstrip("\n")
                    continue

                if line.startswith("#CHROM"):
                    output_full.write(line)
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) == 11:
                        tumor_id = fields.index(tumor_sample_name)
                        normal_id = fields.index(normal_sample_name)
                    else:
                        sample_name = fields[9]
                    continue

                fields = line.rstrip("\n").split("\t")

                if len(fields[3]) > len(fields[4]) and fields[3].startswith(fields[4]):
                    position = str(int(fields[1]) + 1)
                    ref = fields[3][len(fields[4]) :]
                    alt = "-"
                elif len(fields[4]) > len(fields[3]) and fields[4].startswith(fields[3]):
                    position = fields[1]
                    ref = "-"
                    alt = fields[4][len(fields[3]) :]
                else:
                    position = fields[1]
                    ref = fields[3]
                    alt = fields[4]

                variant = {
                    "chromosome": fields[0],
                    "position": position,
                    "ref": ref,
                    "alt": alt,
                }

                raw_data = None
                if use_database and cur and schema_name:
                    raw_data = get_annotation_from_toolset_db(
                        variant=variant,
                        cur=cur,
                        schema_name=schema_name,
                    )

                if raw_data:
                    # Turn DB result into a plain row dict instead of a one-row DataFrame
                    row = dict(raw_data)

                    info_dict = _parse_info_field(fields[7])

                    # Preserve previous behavior:
                    # overwrite only keys already present in DB result
                    for key, value in info_dict.items():
                        if key in row:
                            row[key] = value

                    if (args.sample_type == "germline") or len(fields) == 10:
                        sample_ref, sample_alt = _extract_ad_counts(fields[8], fields[9])
                        if sample_ref is not None and sample_alt is not None:
                            row["t_alt_count"] = sample_alt
                            row["t_ref_count"] = sample_ref
                            row["tumor_f"] = _safe_divide(
                                sample_alt, sample_alt + sample_ref
                            )

                        if (args.sample_type == "germline"):
                            row["Matched_Norm_Sample_Barcode"] = sample_name
                        else:
                            row["Tumor_Sample_Barcode"] = sample_name

                    else:
                        tumor_ref, tumor_alt = _extract_ad_counts(fields[8], fields[tumor_id])
                        if tumor_ref is not None and tumor_alt is not None:
                            row["t_alt_count"] = tumor_alt
                            row["t_ref_count"] = tumor_ref
                            row["tumor_f"] = _safe_divide(
                                tumor_alt, tumor_alt + tumor_ref
                            )

                        normal_ref, normal_alt = _extract_ad_counts(
                            fields[8], fields[normal_id]
                        )
                        if normal_ref is not None and normal_alt is not None:
                            row["n_alt_count"] = normal_alt
                            row["n_ref_count"] = normal_ref

                        if normal_sample_name:
                            row["Matched_Norm_Sample_Barcode"] = normal_sample_name
                        if tumor_sample_name:
                            row["Tumor_Sample_Barcode"] = tumor_sample_name

                    skip_rows.append(row)
                else:
                    output_full.write(line)

        if len(skip_rows) == 0:
            open(output_skip_file, "a").close()
            print(
                f"[OK] Split complete: {len(skip_rows)} variants from database (skip), rest need annotation (noskip)"
            )
            return

        # NOTE: when every variant is served from the cache (empty_full is True),
        # the noskip.vcf is intentionally left header-only. The workflow detects
        # this and routes the chunk past the annotation chain straight to the
        # clinical summary / filtering steps, using this skip.maf as the source.

        skip_df = pd.DataFrame(skip_rows)

        # Render skip.maf in the schema's canonical column order so it stays
        # aligned with freshly-annotated MAFs (and across chunks). Canonical
        # columns first (in canonical order), then any remaining columns.
        if use_database and cur and schema_name:
            try:
                canonical = get_column_order(cur, schema_name)
            except Exception as e:
                print(f"[WARNING] Could not fetch canonical column order: {e}")
                canonical = None
            if canonical:
                ordered = [c for c in canonical if c in skip_df.columns]
                ordered += [c for c in skip_df.columns if c not in ordered]
                skip_df = skip_df[ordered]

        with open(output_skip_file, "w", encoding="utf-8") as f:
            f.writelines(maf_header)

        skip_df.to_csv(
            output_skip_file,
            sep="\t",
            index=False,
            header=True,
            encoding="utf-8",
            mode="a",
        )

        print(
            f"[OK] Split complete: {len(skip_df)} variants from database (skip), rest need annotation (noskip)"
        )

    finally:
        if cur is not None:
            cur.close()
        if conn is not None:
            conn.close()


if __name__ == "__main__":
    main()