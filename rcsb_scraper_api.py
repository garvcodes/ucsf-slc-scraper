import pandas as pd
import requests
import time

# --- Configuration ---
# filename = 'all_slc_genes.csv'       # use this to run on all genes
# filename = 'new_slc_genes.csv'       # use this to run on the updated gene list
filename = 'new_slc_genes.csv'           # use this to run on the updated gene list
output_filename = 'raw_output_web_scraping.csv'

SEARCH_API = 'https://search.rcsb.org/rcsbsearch/v2/query'
DATA_API = 'https://data.rcsb.org/rest/v1/core'


# --- Prokaryotic homolog mapping ---
#
# The RCSB gene-symbol search misses prokaryotic structural homologs of human SLCs
# because those structures are deposited under organism-specific gene names
# (e.g. GltPh, LeuT, XylE) rather than the human symbol. The well-characterized
# prokaryotic SLC homologs in the structural biology literature are listed below
# and pulled in via a separate phase after the main per-SLC scrape.
#
# Each entry attributes the homolog to ONE representative human SLC gene in its
# family (to avoid inflating per-gene counts across SLC1A1/A2/A3/A4/A5 etc.).
#
# Fields:
#   slc            : representative SLC gene in the gene list (must be present, or row is skipped)
#   label          : human-readable description for logging
#   search         : RCSB gene-name terms to OR-search (case-insensitive)
#   organism_filter: source organism substring required to keep a hit (protects against
#                    gene-name homonyms; e.g. `xylE` is also a Pseudomonas dioxygenase).
#                    Empty string = accept any prokaryotic source.
#   pdbs           : explicit PDB IDs to include unconditionally (catches structures with
#                    empty rcsb_gene_name fields that the search misses). If such a PDB is
#                    deposited as "synthetic construct" (no organism in the record), it is
#                    kept and classified as 'Synthetic' rather than guessing its kingdom.
SLC_PROKARYOTIC_HOMOLOGS = [
    # SLC1 family — glutamate / neutral-amino-acid transporters
    {"slc": "SLC1A1", "label": "GltPh (Pyrococcus horikoshii)",
     "search": ["gltPh", "PH1295", "O59010"],
     "organism_filter": "Pyrococcus",
     "pdbs": ["2NWL", "2NWW", "2NWX", "4OYF", "5DWX"]},
    {"slc": "SLC1A1", "label": "GltTk (Thermococcus kodakarensis)",
     "search": ["TK0986"],
     "organism_filter": "Thermococcus",
     "pdbs": []},

    # SLC2 family — facilitative sugar transporters (MFS fold)
    {"slc": "SLC2A1", "label": "XylE (E. coli)",
     "search": ["xylE"],
     "organism_filter": "Escherichia coli",
     "pdbs": []},

    # SLC5 family — Na+/sugar cotransporters (Faham 2008 Science)
    {"slc": "SLC5A1", "label": "vSGLT (Vibrio parahaemolyticus)",
     "search": [],
     "organism_filter": "Vibrio parahaemolyticus",
     "pdbs": ["2XQ2", "3DH4"]},

    # SLC6 family — Na+/Cl− neurotransmitter transporters (LeuT-fold)
    {"slc": "SLC6A4", "label": "LeuT (Aquifex aeolicus)",
     "search": ["leuT", "snf", "aq_2077"],
     "organism_filter": "Aquifex aeolicus",
     "pdbs": ["2A65"]},

    # SLC7 family — amino acid antiporters (APC superfamily)
    {"slc": "SLC7A11", "label": "AdiC (E. coli)",
     "search": ["adiC"],
     "organism_filter": "Escherichia coli",
     "pdbs": []},
    {"slc": "SLC7A5", "label": "ApcT (bacterial/archaeal homologs)",
     "search": ["apcT", "MJ0609"],
     "organism_filter": "",
     "pdbs": ["3GIA"]},

    # SLC11 family — divalent metal transporters (Nramp)
    {"slc": "SLC11A1", "label": "Bacterial Nramp / MntH",
     "search": ["mntH"],
     "organism_filter": "",
     "pdbs": []},

    # SLC13 family — Na+/dicarboxylate cotransporters
    {"slc": "SLC13A5", "label": "VcINDY (Vibrio cholerae)",
     "search": ["VC_A0025"],
     "organism_filter": "Vibrio cholerae",
     "pdbs": []},

    # SLC15 family — proton-coupled oligopeptide transporters (POT/PTR)
    {"slc": "SLC15A1", "label": "Bacterial POT / PepT homologs",
     "search": ["yjdL", "SO_1277"],
     "organism_filter": "",
     "pdbs": ["2XUT", "4UVM"]},

    # SLC28 family — concentrative nucleoside transporters
    {"slc": "SLC28A3", "label": "vcCNT (Vibrio cholerae)",
     "search": ["VC_2352"],
     "organism_filter": "Vibrio cholerae",
     "pdbs": []},

    # SLC30 family — Zn2+ exporters (CDF fold)
    {"slc": "SLC30A8", "label": "YiiP (E. coli)",
     "search": ["yiiP"],
     "organism_filter": "Escherichia coli",
     "pdbs": []},

    # SLC40 family — ferroportin
    {"slc": "SLC40A1", "label": "BbFpn (Bdellovibrio bacteriovorus)",
     "search": ["Bd2019"],
     "organism_filter": "Bdellovibrio",
     "pdbs": []},
]


# --- Helper functions ---

def search_pdb_ids(gene_name):
    """Search RCSB for all PDB IDs matching an exact gene name."""
    query = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_entity_source_organism.rcsb_gene_name.value",
                "operator": "exact_match",
                "value": gene_name
            }
        },
        "return_type": "entry",
        "request_options": {
            "paginate": {"start": 0, "rows": 10000}
        }
    }
    response = requests.post(SEARCH_API, json=query)
    if response.status_code == 204:
        return []  # no results
    response.raise_for_status()
    data = response.json()
    return [hit['identifier'] for hit in data.get('result_set', [])]


def get_entry_data(pdb_id):
    """Fetch entry-level metadata (method, resolution, authors, etc.)."""
    resp = requests.get(f'{DATA_API}/entry/{pdb_id}')
    resp.raise_for_status()
    return resp.json()


def get_assembly_data(pdb_id):
    """Fetch assembly data for oligomeric state."""
    resp = requests.get(f'{DATA_API}/assembly/{pdb_id}/1')
    if resp.status_code == 404:
        return {}
    resp.raise_for_status()
    return resp.json()


def get_all_entities(pdb_id, entry):
    """Fetch every polymer entity for a PDB. Returns {entity_id: entity_dict}.

    Heteromer structures (SLC + fusion partner / toxin / viral protein) have multiple
    polymer entities, and the SLC is not always entity 1 — we have to walk them all.
    """
    eids = entry.get('rcsb_entry_container_identifiers', {}).get('polymer_entity_ids', ['1'])
    entities = {}
    for eid in eids:
        resp = requests.get(f'{DATA_API}/polymer_entity/{pdb_id}/{eid}')
        if resp.status_code == 404:
            continue
        resp.raise_for_status()
        entities[eid] = resp.json()
        time.sleep(0.1)
    return entities


def resolve_homolog_source(entities, organism_filter, is_explicit):
    """Resolve a source_organism dict for a prokaryotic homolog PDB.

    1. Walk all entities + source_organism records, picking the first one that is
       (a) prokaryotic (Bacteria/Archaea), AND (b) matches the organism_filter
       substring if one was provided.
    2. If none match and the PDB came from an explicit ID list, keep it anyway and
       return whatever the record actually carries (typically 'synthetic construct'
       with no taxonomy) — extract_metadata classifies that as 'Synthetic' rather
       than fabricating an organism.
    3. Otherwise return None — the PDB came from a gene-name search but doesn't
       actually match the homolog's expected organism (a false positive).
    """
    target = organism_filter.lower() if organism_filter else ""
    for ent in entities.values():
        for src in (ent.get('rcsb_entity_source_organism') or []):
            lin = src.get('taxonomy_lineage', []) or []
            sk = next((t['name'] for t in lin if t.get('name') in {'Bacteria', 'Archaea'}), '')
            if not sk:
                continue
            if target:
                org = (src.get('ncbi_scientific_name') or '').lower()
                if target not in org:
                    continue
            return src
    if is_explicit:
        for ent in entities.values():
            srcs = ent.get('rcsb_entity_source_organism') or []
            if srcs:
                return srcs[0]
        return {'ncbi_scientific_name': 'synthetic construct', 'taxonomy_lineage': [], 'rcsb_gene_name': []}
    return None


def pick_source_organism(entities, gene):
    """Find the source_organism dict whose rcsb_gene_name matches `gene`.

    Walks every entity AND every source_organism within each entity, since a single
    entity can carry multiple source organisms when it's a fusion construct (e.g. an
    E. coli cybC tag fused to a human SLC). Falls back to the first available source
    organism if no gene-name match is found.
    """
    target = gene.upper()
    for ent in entities.values():
        for src in (ent.get('rcsb_entity_source_organism') or []):
            for g in (src.get('rcsb_gene_name') or []):
                if (g.get('value') or '').upper() == target:
                    return src
    for ent in entities.values():
        srcs = ent.get('rcsb_entity_source_organism') or []
        if srcs:
            return srcs[0]
    return {}


def extract_metadata(gene, pdb_id, entry, assembly, source_organism, partner_genes):
    """Extract all needed fields from the API responses."""
    row = {'GENE': gene, 'PDB ID': pdb_id}
    row['PARTNER GENES'] = ', '.join(sorted(partner_genes))

    # Method (EM / X-RAY DIFFRACTION / etc.)
    exptl = entry.get('exptl', [{}])
    row['METHOD'] = exptl[0].get('method', '') if exptl else ''

    # Resolution — EM or X-ray
    em3d = entry.get('em_3d_reconstruction', [])
    refine = entry.get('refine', [])
    if em3d and em3d[0].get('resolution'):
        row['RESOLUTION'] = em3d[0]['resolution']
    elif refine and refine[0].get('ls_d_res_high'):
        row['RESOLUTION'] = refine[0]['ls_d_res_high']
    else:
        row['RESOLUTION'] = ''

    # Release date
    row['RELEASE DATE'] = entry.get('rcsb_accession_info', {}).get('initial_release_date', '')

    # Source organism (the one whose rcsb_gene_name matches the SLC, picked upstream)
    row['SOURCE ORGANISM'] = source_organism.get('ncbi_scientific_name', '')

    # Kingdom: pull the superkingdom node from the NCBI taxonomy lineage
    lineage = source_organism.get('taxonomy_lineage', []) or []
    row['KINGDOM'] = next(
        (t['name'] for t in lineage
         if t.get('name') in {'Bacteria', 'Archaea', 'Eukaryota', 'Viruses'}),
        ''
    )

    # Prokaryotic / Eukaryotic category derived from the superkingdom above,
    # materialized as its own column so downstream analysis can group on it directly.
    # Structures deposited as "synthetic construct" carry no superkingdom, so they get
    # their own honest 'Synthetic' bucket rather than being guessed or left blank.
    row['PROK_EUK'] = {
        'Bacteria': 'Prokaryotic',
        'Archaea': 'Prokaryotic',
        'Eukaryota': 'Eukaryotic',
        'Viruses': 'Virus',
    }.get(row['KINGDOM'], '')
    if not row['PROK_EUK'] and 'synthetic construct' in (row['SOURCE ORGANISM'] or '').lower():
        row['PROK_EUK'] = 'Synthetic'

    # Oligomeric state
    pdbx_assembly = assembly.get('pdbx_struct_assembly', {})
    row['OLIGOMERIC STATE'] = pdbx_assembly.get('oligomeric_details', '')

    # Molecular weight
    row['MOLECULAR WEIGHT'] = entry.get('rcsb_entry_info', {}).get('molecular_weight', '')

    # Last (senior) author
    authors = entry.get('audit_author', [])
    row['STRUCTURE AUTHOR'] = authors[-1].get('name', '') if authors else ''

    # Publication year
    citations = entry.get('citation', [])
    primary = [c for c in citations if c.get('id') == 'primary']
    row['PUBLICATION YEAR'] = primary[0].get('year', '') if primary else ''

    return row


# --- Main script ---

columns = [
    'GENE', 'PDB ID', 'PARTNER GENES', 'METHOD', 'RESOLUTION', 'RELEASE DATE',
    'SOURCE ORGANISM', 'KINGDOM', 'PROK_EUK', 'OLIGOMERIC STATE', 'MOLECULAR WEIGHT',
    'STRUCTURE AUTHOR', 'PUBLICATION YEAR'
]


if __name__ == '__main__':
    print('Program starting')
    genes_df = pd.read_csv(filename)
    genes = [r.iloc[0].strip() for _, r in genes_df.iterrows()]

    # Phase 1: search every gene, build the PDB -> {genes} map
    print('\n--- Searching ---')
    gene_to_pdbs = {}
    pdb_to_genes = {}
    for gene in genes:
        pdb_ids = search_pdb_ids(gene)
        gene_to_pdbs[gene] = pdb_ids
        print(f'  {gene}: {len(pdb_ids)} structures')
        for pdb_id in pdb_ids:
            pdb_to_genes.setdefault(pdb_id, set()).add(gene)

    # Phase 2: fetch metadata once per unique PDB ID, emit one row per (gene, pdb_id)
    print('\n--- Fetching metadata ---')
    entry_cache, assembly_cache, entities_cache = {}, {}, {}
    all_rows = []

    for gene in genes:
        pdb_ids = gene_to_pdbs[gene]
        if not pdb_ids:
            continue
        print(f'\n--- {gene} ---')
        for pdb_id in pdb_ids:
            try:
                if pdb_id not in entry_cache:
                    entry_cache[pdb_id] = get_entry_data(pdb_id)
                    assembly_cache[pdb_id] = get_assembly_data(pdb_id)
                    entities_cache[pdb_id] = get_all_entities(pdb_id, entry_cache[pdb_id])
                    time.sleep(0.1)  # be polite to the API

                partners = pdb_to_genes[pdb_id] - {gene}
                source_org = pick_source_organism(entities_cache[pdb_id], gene)
                metadata = extract_metadata(
                    gene, pdb_id,
                    entry_cache[pdb_id], assembly_cache[pdb_id], source_org,
                    partners,
                )
                if metadata['PUBLICATION YEAR'] == '':
                    print(f'  {pdb_id}: skipped (no publication year)')
                    continue
                all_rows.append(metadata)
                tag = f' [partners: {metadata["PARTNER GENES"]}]' if partners else ''
                print(f'  {pdb_id}: {metadata["METHOD"]}, res={metadata["RESOLUTION"]}{tag}')
            except Exception as e:
                print(f'  {pdb_id}: ERROR - {e}')

    # Phase 3: prokaryotic homologs missed by the gene-symbol search
    print('\n--- Adding prokaryotic homologs ---')
    existing_rows = {(r['GENE'], r['PDB ID']) for r in all_rows}
    homolog_added = 0
    for homolog in SLC_PROKARYOTIC_HOMOLOGS:
        slc = homolog['slc']
        if slc not in gene_to_pdbs:
            print(f"\n  {homolog['label']} -> {slc}: skipped (representative SLC not in gene list)")
            continue

        explicit_pdbs = set(homolog['pdbs'])
        search_pdbs = set()
        for term in homolog['search']:
            try:
                search_pdbs.update(search_pdb_ids(term))
                time.sleep(0.1)
            except Exception as e:
                print(f"  search {term!r} failed: {e}")

        candidates = sorted(explicit_pdbs | search_pdbs)
        added, skipped_org, skipped_dup = 0, 0, 0
        print(f"\n--- {homolog['label']} -> {slc}: {len(candidates)} candidates "
              f"({len(explicit_pdbs)} explicit, {len(search_pdbs)} from search) ---")

        for pdb_id in candidates:
            if (slc, pdb_id) in existing_rows:
                skipped_dup += 1
                continue
            try:
                if pdb_id not in entry_cache:
                    entry_cache[pdb_id] = get_entry_data(pdb_id)
                    assembly_cache[pdb_id] = get_assembly_data(pdb_id)
                    entities_cache[pdb_id] = get_all_entities(pdb_id, entry_cache[pdb_id])
                    time.sleep(0.1)

                source_org = resolve_homolog_source(
                    entities_cache[pdb_id],
                    homolog['organism_filter'],
                    is_explicit=(pdb_id in explicit_pdbs),
                )
                if source_org is None:
                    print(f"  {pdb_id}: skipped (no prokaryotic source matching "
                          f"{homolog['organism_filter']!r})")
                    skipped_org += 1
                    continue

                metadata = extract_metadata(
                    slc, pdb_id,
                    entry_cache[pdb_id], assembly_cache[pdb_id], source_org, set(),
                )
                if metadata['PUBLICATION YEAR'] == '':
                    print(f"  {pdb_id}: skipped (no publication year)")
                    skipped_org += 1
                    continue
                all_rows.append(metadata)
                existing_rows.add((slc, pdb_id))
                added += 1
                homolog_added += 1
                print(f"  {pdb_id}: {metadata['SOURCE ORGANISM']} "
                      f"({metadata['KINGDOM']}) res={metadata['RESOLUTION']}")
            except Exception as e:
                print(f"  {pdb_id}: ERROR - {e}")

        print(f"  --> {added} added, {skipped_dup} duplicates, {skipped_org} skipped")

    print(f"\n--- Homolog phase done: {homolog_added} prokaryotic rows added ---")

    output_df = pd.DataFrame(all_rows, columns=columns)
    output_df.to_csv(output_filename, index=False)
    xlsx_filename = output_filename.replace('.csv', '.xlsx')
    output_df.to_excel(xlsx_filename, index=False)
    print(f'\nDone! Saved {len(all_rows)} rows to {output_filename} and {xlsx_filename}')
