import configparser
import glob
import hashlib
import json
import os
import re
import sqlite3  # remove if you've fully migrated to Postgres
import time
from typing import Dict, List, Optional, Tuple

import pandas as pd

import psycopg
from psycopg.rows import dict_row

# support for psycopg3 execute_values
def execute_values(cur, sql, argslist, template=None, page_size=1000):
        cur.executemany(sql if template is None else sql, argslist)

# ---------- Registry (public) ----------

DDL_VERSION = "2.0-toolset"  # bump if your per-toolset schema changes

def generate_toolset_hash(tools_info: List[Dict]) -> Tuple[str, str]:
    """
    Generate a stable SHA256 hash for a tool+db combo list.
    Returns (toolset_id, tools_json_string).
    """
    # normalize and sort for stability
    normalized = sorted(
        [
            {
                "tool": {
                    "name": t["tool_info"]["name"].strip().lower(),
                    "version": str(t["tool_info"]["version"]).strip(),
                },
                "db": {
                    "name": t["db_info"]["name"].strip().lower(),
                    "version": str(t["db_info"]["version"]).strip(),
                },
            }
            for t in tools_info
        ],
        key=lambda x: (x["tool"]["name"], x["tool"]["version"], x["db"]["name"], x["db"]["version"]),
    )
    json_str = json.dumps(normalized, sort_keys=True, separators=(",", ":"))
    toolset_id = hashlib.sha256(json_str.encode("utf-8")).hexdigest()
    return toolset_id, json_str

def _build_dsn(
    *,
    explicit_dsn: Optional[str],
    host: Optional[str],
    port: Optional[str],
    db: Optional[str],
    user: Optional[str],
    password: Optional[str],
) -> str:
    """Resolve a Postgres DSN from --dsn or PG* env vars / flags."""
    if explicit_dsn:
        return explicit_dsn

    # Fallback to env vars then flags
    pg_host = host or os.getenv("PGHOST", "localhost")
    pg_port = port or os.getenv("PGPORT", "5432")
    pg_db   = db   or os.getenv("PGDATABASE")
    pg_user = user or os.getenv("PGUSER")
    pg_pwd  = password or os.getenv("PGPASSWORD", "")

    if not pg_db or not pg_user:
        raise ValueError(
            "Missing database/user. Provide --dsn or --pg-db and --pg-user "
            "(or set PGDATABASE and PGUSER)."
        )

    auth = f"{pg_user}:{pg_pwd}@" if pg_pwd else f"{pg_user}@"
    return f"postgresql://{auth}{pg_host}:{pg_port}/{pg_db}"

def init_registry_db(dsn: str):
    """
    Creates the central registry in the 'public' schema, if needed.
    """
    with psycopg.connect(dsn, autocommit=True) as conn:
        with conn.cursor() as cur:
            cur.execute("""
            CREATE TABLE IF NOT EXISTS public."Toolset_Registry" (
                toolset_id   TEXT PRIMARY KEY,                -- SHA256 of tools_info
                tools_json   JSONB NOT NULL,                  -- JSONB for easy querying
                schema_name  TEXT NOT NULL,                   -- per-toolset schema
                ddl_version  TEXT NOT NULL,                   -- to detect schema drift
                created_at   TIMESTAMPTZ DEFAULT now()
            );
            """)



# ---------- Toolset schema (per toolset) ----------

def init_toolset_db(dsn: str, schema_name: str):
    """
    Initialize the per-toolset database schema.

    This creates a dedicated PostgreSQL schema for the toolset and defines
    a single table, <schema_name>.variant, used to store variant records.
    If the schema or table already exist, they are left unchanged.

    The following objects are created:

        SCHEMA:
            <schema_name>

        TABLE:
            <schema_name>.variant
                - variant_id      BIGSERIAL PRIMARY KEY
                - contig          TEXT NOT NULL
                - position        INTEGER NOT NULL
                - ref_allele      TEXT NOT NULL
                - alt_allele      TEXT NOT NULL
                - raw_annotation  JSONB NOT NULL
                - annotation_hash CHAR(64) NOT NULL
                - gene_symbol     TEXT
                - cancervar       TEXT
                - intervar        TEXT
                - clinvar         TEXT
                - cosmic          TEXT
                - UNIQUE (contig, position, ref_allele, alt_allele)

        INDEX:
            idx_variant_gene ON <schema_name>.variant(gene_symbol)
    """

    with psycopg.connect(dsn, autocommit=True) as conn:
        with conn.cursor() as cur:
            # Create schema
            cur.execute(f'CREATE SCHEMA IF NOT EXISTS "{schema_name}";')
            # Create table
            cur.execute(f"""
            CREATE TABLE IF NOT EXISTS "{schema_name}".variant (
                variant_id      BIGSERIAL PRIMARY KEY,
                contig          TEXT NOT NULL,
                position        INTEGER NOT NULL,
                ref_allele      TEXT NOT NULL,
                alt_allele      TEXT NOT NULL,
                ref_allele_hash BYTEA NOT NULL,
                alt_allele_hash BYTEA NOT NULL,
                raw_annotation  JSONB NOT NULL,
                annotation_hash CHAR(64) NOT NULL,
                gene_symbol     TEXT,
                cancervar       TEXT,
                intervar        TEXT,
                clinvar         TEXT,
                cosmic          TEXT
            );
            """)

            # Unique index for variant identity
            cur.execute(f"""
            CREATE UNIQUE INDEX IF NOT EXISTS variant_unique_idx
            ON "{schema_name}".variant
            (contig, position, ref_allele_hash, alt_allele_hash);
            """)

            # Optional gene index
            cur.execute(f"""
            CREATE INDEX IF NOT EXISTS idx_variant_gene
            ON "{schema_name}".variant (gene_symbol);
            """)

            # Canonical MAF column order for this toolset (single-row table).
            # Stores the union of all column orders ever seen, merged with an
            # order-preserving rule (see merge_column_order). Used to render the
            # *_skip.maf deterministically when no fresh Funcotator MAF exists.
            cur.execute(f"""
            CREATE TABLE IF NOT EXISTS "{schema_name}".column_order (
                id         INT PRIMARY KEY DEFAULT 1 CHECK (id = 1),
                columns    JSONB NOT NULL,
                updated_at TIMESTAMPTZ DEFAULT now()
            );
            """)

def _schema_from_toolset(toolset_id: str) -> str:
    """
    Derive a compact, valid schema name from the toolset hash.
    """
    short = toolset_id[:12]
    return f"toolset_{short}"  # e.g., toolset_a1b2c3d4


def get_or_create_schema_for_tools(
    dsn: str,
    tools_info: List[Dict],
    *,
    ddl_version: str = DDL_VERSION,
) -> Tuple[str, str]:
    """
    Ensure a per-toolset schema exists and is registered, returning (toolset_id, schema_name).

    - Inserts a row in public.Toolset_Registry if missing.
    - Creates/initializes the per-toolset schema (partitioned tables, indexes).
    """
    toolset_id, tools_json = generate_toolset_hash(tools_info)
    schema_name = _schema_from_toolset(toolset_id)

    init_registry_db(dsn)  # idempotent

    with psycopg.connect(dsn, autocommit=True) as conn:
        with conn.cursor(row_factory=dict_row) as cur:
            # Try to insert the mapping (race-safe)
            cur.execute("""
                INSERT INTO public."Toolset_Registry" (toolset_id, tools_json, schema_name, ddl_version)
                VALUES (%s, %s::jsonb, %s, %s)
                ON CONFLICT (toolset_id) DO NOTHING;
            """, (toolset_id, tools_json, schema_name, ddl_version))

            # Fetch current mapping
            cur.execute("""
                SELECT toolset_id, tools_json, schema_name, ddl_version
                FROM public."Toolset_Registry"
                WHERE toolset_id = %s;
            """, (toolset_id,))
            row = cur.fetchone()

            # If an older ddl_version is stored, you can trigger a migration here
            if row and row["ddl_version"] != ddl_version:
                # Optional: run migrations for existing schema, then bump version in registry
                cur.execute('UPDATE public."Toolset_Registry" SET ddl_version = %s WHERE toolset_id = %s;',
                            (ddl_version, toolset_id))

    # Ensure the per-toolset schema and partitioned tables exist
    # (this is your previously provided init_toolset_db(dsn, schema_name))
    init_toolset_db(dsn, schema_name)

    return toolset_id, schema_name


def merge_column_order(canonical: List[str], new: List[str]) -> List[str]:
    """
    Merge two column orderings into a union that preserves relative order.

    The canonical ordering is authoritative for columns shared by both lists
    (canonical wins on conflicts). Columns present only in ``new`` are inserted
    immediately after the column they follow in ``new`` (their most recent
    "anchor", i.e. a column also present in ``canonical``). New-only columns
    that appear before any anchor are placed at the front.

    Example:
        canonical = [A, B], new = [A, C]   (C comes right after A)
        -> [A, C, B]
    """
    canonical = list(canonical)
    new = list(new)
    canon_set = set(canonical)

    # Bucket new-only columns by the anchor they follow (None -> front).
    front: List[str] = []
    insert_after: Dict[str, List[str]] = {}
    last_anchor: Optional[str] = None
    for col in new:
        if col in canon_set:
            last_anchor = col
        else:
            if last_anchor is None:
                front.append(col)
            else:
                insert_after.setdefault(last_anchor, []).append(col)

    result: List[str] = []
    seen = set()

    def _emit(col: str) -> None:
        if col not in seen:
            result.append(col)
            seen.add(col)

    for col in front:
        _emit(col)
    for col in canonical:
        _emit(col)
        for nc in insert_after.get(col, []):
            _emit(nc)

    return result


def get_column_order(cur, schema_name: str) -> Optional[List[str]]:
    """Return the stored canonical column order for a schema, or None if unset."""
    cur.execute(f'SELECT columns FROM "{schema_name}".column_order WHERE id = 1;')
    row = cur.fetchone()
    if not row:
        return None
    cols = row[0] if not isinstance(row, dict) else row.get("columns")
    if isinstance(cols, str):
        try:
            cols = json.loads(cols)
        except Exception:
            return None
    return list(cols) if cols else None


def update_column_order(dsn: str, schema_name: str, new_columns: List[str]) -> List[str]:
    """
    Merge ``new_columns`` into the schema's stored canonical column order and
    persist the result. Returns the merged order.
    """
    new_columns = list(new_columns)
    with psycopg.connect(dsn, autocommit=True) as conn:
        with conn.cursor() as cur:
            existing = get_column_order(cur, schema_name)
            merged = (
                merge_column_order(existing, new_columns)
                if existing
                else new_columns
            )
            cur.execute(
                f"""
                INSERT INTO "{schema_name}".column_order (id, columns, updated_at)
                VALUES (1, %s::jsonb, now())
                ON CONFLICT (id) DO UPDATE
                SET columns = EXCLUDED.columns,
                    updated_at = now();
                """,
                (json.dumps(merged),),
            )
    return merged


def hash_annotation(annotation_dict, decimal_places=8):
    """Generate a hash of the annotation dict after normalizing floats."""

    def normalize_value(val):
        if isinstance(val, float):
            return round(val, decimal_places)
        elif isinstance(val, str):
            try:
                return round(float(val), decimal_places)
            except ValueError:
                return val
        return val

    normalized = {k: normalize_value(v) for k, v in annotation_dict.items()}
    json_str = json.dumps(normalized)
    return hashlib.sha256(json_str.encode("utf-8")).hexdigest(), json_str


def get_annotation_from_toolset_db(
    variant: dict,
    cur,          # psycopg 3 cursor
    schema_name: str,
):
    """
    Retrieve the annotation for a given variant from the per-toolset Postgres table.

    Parameters
    ----------
    variant : dict
        Must contain: "chromosome" (or "contig"), "position", "ref", "alt".
    cur : psycopg.Cursor
        psycopg 3 cursor.
    schema_name : str
        Per-toolset schema name (e.g., 'toolset_ab12cd34').

    Returns
    -------
    dict | None
        The parsed annotation (JSON) if found, else None.
    """
    # Accept 'chromosome' or 'contig'
    contig = variant.get("contig", variant.get("chromosome"))
    if contig is None:
        raise ValueError("variant must include 'chromosome' or 'contig'")

    position   = int(variant["position"])
    ref_allele = variant["ref"]
    alt_allele = variant["alt"]

    table = f'"{schema_name}".variant'

    ref_allele_hash = hashlib.sha256(ref_allele.encode()).digest()
    alt_allele_hash = hashlib.sha256(alt_allele.encode()).digest()

    # Direct lookup on the unique key (contig, position, ref_allele_hash, alt_allele_hash)
    cur.execute(
        f"""
        SELECT raw_annotation
        FROM {table}
        WHERE contig     = %s
          AND position   = %s
          AND ref_allele_hash = %s
          AND alt_allele_hash = %s
        LIMIT 1
        """,
        (contig, position, ref_allele_hash, alt_allele_hash),
    )
    row = cur.fetchone()
    if not row and contig == "chrM":
        contig = "MT"
        # Direct lookup on the unique key (contig, position, ref_allele_hash, alt_allele_hash)
        cur.execute(
            f"""
            SELECT raw_annotation
            FROM {table}
            WHERE contig     = %s
            AND position   = %s
            AND ref_allele_hash = %s
            AND alt_allele_hash = %s
            LIMIT 1
            """,
            (contig, position, ref_allele_hash, alt_allele_hash),
        )
        row = cur.fetchone()
    if not row:
        print(f"[ERROR] Variant {contig}:{position}:{ref_allele}:{alt_allele} not found.")
        return None

    # psycopg3 may return JSONB as a Python object or as a str, depending on settings
    raw = row[0] if not isinstance(row, dict) else row.get("raw_annotation")
    if isinstance(raw, (dict, list)):
        return raw
    try:
        return json.loads(raw)
    except Exception:
        print("[WARNING] Failed to decode annotation JSON.")
        return None


def insert_many_annotations_to_toolset_db(
    maf: pd.DataFrame,
    tools_info: List[Dict],
    *,
    dsn: str,
    schema_name: Optional[str] = None,
    batch_size: int = 10_000,
) -> Tuple[str, str]:
    """
    Insert many variant annotations into the per-toolset Postgres table
    with batched UPSERTs.

    Parameters
    ----------
    maf : pd.DataFrame
        MAF-like DataFrame. Must provide:
        - Chromosome (or 'contig')
        - Start_Position
        - Reference_Allele
        - Tumor_Seq_Allele2
        Other columns are stored in raw_annotation (JSONB).
    tools_info : list[dict]
        Tool/db descriptors used to derive the toolset hash -> schema.
    dsn : str
        Postgres DSN (e.g., 'postgresql://user:password@host:5432/dbname').
    schema_name : str | None
        If None, it's created/resolved via get_or_create_schema_for_tools().
    batch_size : int
        Number of rows per batch for execute_values.

    Returns
    -------
    (toolset_id, schema_name)
    """

    toolset_id: str
    if schema_name is None:
        toolset_id, schema_name = get_or_create_schema_for_tools(dsn, tools_info)
    else:
        # still compute the id for the caller's awareness
        toolset_id, _ = generate_toolset_hash(tools_info)

    table = f'"{schema_name}".variant'

    upsert_sql = f"""
        INSERT INTO {table} (
            contig, position, ref_allele, alt_allele, ref_allele_hash, alt_allele_hash,
            raw_annotation, annotation_hash,
            gene_symbol, cancervar, intervar, clinvar, cosmic
        )
        VALUES (
            %s, %s, %s, %s, %s, %s,
            %s::jsonb, %s,
            %s, %s, %s, %s, %s
        )
        ON CONFLICT (contig, position, ref_allele_hash, alt_allele_hash)
        DO UPDATE SET
            raw_annotation = EXCLUDED.raw_annotation,
            annotation_hash = EXCLUDED.annotation_hash,
            gene_symbol = EXCLUDED.gene_symbol,
            cancervar = EXCLUDED.cancervar,
            intervar = EXCLUDED.intervar,
            clinvar = EXCLUDED.clinvar,
            cosmic = EXCLUDED.cosmic
        WHERE {table}.annotation_hash IS DISTINCT FROM EXCLUDED.annotation_hash
    """

    # Template: cast raw_annotation to jsonb on the server
    values_template = "(%s,%s,%s,%s,%s,%s,%s::jsonb,%s,%s,%s,%s,%s,%s)"

    def build_batch_params(df: pd.DataFrame) -> List[tuple]:
        params: List[tuple] = []
        for rec in df.to_dict(orient="records"):
            contig = rec.get("contig", rec.get("Chromosome"))
            if contig is None:
                raise ValueError("Missing 'Chromosome'/'contig' column")

            try:
                position = int(rec["Start_Position"])
            except KeyError:
                raise ValueError("Missing 'Start_Position' column")

            ref_allele = rec.get("Reference_Allele")
            alt_allele = rec.get("Tumor_Seq_Allele2")
            if ref_allele is None or alt_allele is None:
                raise ValueError("Missing 'Reference_Allele' or 'Tumor_Seq_Allele2'")

            print(f"[INSERT] {contig}:{position}:{ref_allele}:{alt_allele}")

            ref_allele_hash = hashlib.sha256(ref_allele.encode()).digest()
            alt_allele_hash = hashlib.sha256(alt_allele.encode()).digest()

            # Treat nan
            cleaned = {}
            for k, v in rec.items():
                if pd.isna(v):
                    cleaned[k] = None
                else:
                    cleaned[k] = v
            annotation_dict = cleaned
            ann_hash, raw_json = hash_annotation(annotation_dict)


            params.append((
                contig,
                position,
                ref_allele,
                alt_allele,
                ref_allele_hash,
                alt_allele_hash,
                raw_json,   # -> %s::jsonb
                ann_hash,
                annotation_dict.get("Hugo_Symbol"),
                annotation_dict.get("CancerVar"),
                annotation_dict.get("InterVar"),
                annotation_dict.get("ClinVar_VCF_CLNSIG"),
                annotation_dict.get("cosmic"),
            ))
        return params

    with psycopg.connect(dsn) as conn:
        with conn.cursor() as cur:
            n = len(maf)
            if n == 0:
                return toolset_id, schema_name

            for start in range(0, n, batch_size):
                end = min(start + batch_size, n)
                chunk = maf.iloc[start:end]
                params = build_batch_params(chunk)            

                execute_values(
                    cur,
                    upsert_sql,
                    params,
                    template=values_template,
                    page_size=len(params),
                )
            conn.commit()

    # Persist the canonical column order (union, order-preserving) for this toolset
    # so cached variants can be rendered deterministically without re-annotation.
    try:
        update_column_order(dsn, schema_name, list(maf.columns))
    except Exception as e:
        print(f"[WARNING] Failed to update canonical column order: {e}")

    print(f"[OK] Inserted/updated {len(maf)} rows into schema '{schema_name}'.")
    return toolset_id, schema_name


def get_version(file: str):
    """
    file: configuration file used to extract database version.
    """
    with open(file) as f:
        for line in f:
            if line.startswith("version = "):
                v = line.split("=")[1]
                return v.strip()


def get_name(file: str):
    """
    file: configuration file used to extract database version.
    """
    with open(file) as f:
        for line in f:
            if line.startswith("name = "):
                v = line.split("=")[1]
                return v.strip().lower()


def get_funcotator_db(funcotator_db_path: str, tools_info: list):
    configs = glob.glob(f"{funcotator_db_path}/*/*/*config")

    # Funcotator
    tool_name = "funcotator"
    tool_version = "4.5.0.0"

    for config in configs:
        version = get_version(config)
        name = get_name(config)
        current_dict = {
            "tool_info": {"name": tool_name, "version": tool_version},
            "db_info": {"name": name, "version": version},
        }
        tools_info.append(current_dict)
    return tools_info


def get_annovar_db(config_file: str, tools_info: list, annovar_version="v0"):
    # Preprocess to strip inline comments (lines with '#')
    clean_lines = []
    with open(config_file, "r") as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                # Remove trailing comments
                if "#" in line:
                    line = line.split("#")[0].strip()
                clean_lines.append(line)

    # Use configparser to parse cleaned config
    config = configparser.ConfigParser()
    config.read_string("\n".join(clean_lines))

    # Extract Annovar database names
    annovar_dbs = config.get("Annovar", "database_names").split()

    # Extract Intervar version
    try:
        intervar, version = config.get("Other", "current_version").split("_")

    except:
        intervar, version = "CancerVar", "v0"
    tools_info.append(
        {
            "tool_info": {"name": "annovar", "version": annovar_version},
            "db_info": {"name": intervar, "version": version},
        }
    )
    for name in annovar_dbs:
        version = annovar_version
        current_dict = {
            "tool_info": {"name": "annovar", "version": annovar_version},
            "db_info": {"name": name, "version": version},
        }
        tools_info.append(current_dict)

    return tools_info


def get_tools_info(
    funcotator_db_path: str,
    annovar_config_file: str,
    technology: str,
    annovar_version="v0",
    genome_build="GRCh38",
    germline=False,
):
    if germline:
        germline = "germline"
    else:
        germline = "somatic"

    tools_info = [
        {
            "tool_info": {"name": technology, "version": "v0"},
            "db_info": {"name": "none", "version": "none"},
        },
        {
            "tool_info": {"name": genome_build, "version": "v0"},
            "db_info": {"name": genome_build, "version": "v0"},
        },
        {
            "tool_info": {"name": germline, "version": "v0"},
            "db_info": {"name": germline, "version": "v0"},
        },
    ]
    tools_info = get_funcotator_db(funcotator_db_path, tools_info)
    tools_info = get_annovar_db(annovar_config_file, tools_info, annovar_version)
    return tools_info
