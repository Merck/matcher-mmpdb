import psycopg2
from psycopg2.extras import execute_values
import time
from rdkit import Chem

# Objective: setup everything needed in the DB beyond what's provided in the original mmpdb open source package
def extend_postgres_build(connection=None, cursor=None):
    if connection is None or cursor is None:
        raise Exception("Both connection and cursor objects must be passed to this function")
    conn = connection
    c = cursor

    DDL_statements = [
        # Create molecule columns for rdkit cartridge structure searching
        "ALTER TABLE compound ADD clean_smiles_mol mol",    
        "ALTER TABLE rule_smiles ADD smiles_mol mol",    
        "ALTER TABLE constant_smiles ADD smiles_mol mol",
        # Create "construct" tables
        # A construct is a way to represent a compound as a combination of a constant and variable part
        # One compound can be represented by multiple constructs depending on how the molecule is fragmented into constant/variable parts
        # We create two construct tables to map to the directionality of the rule table: from_construct maps to from_rule_smiles, and to_construct maps to to_rule_smiles
        # These construct tables establish a direct link between the constant/variable fragments and the compounds they belong to, which is the backbone of matcher's open-ended querying
        """
CREATE TABLE from_construct (
    id SERIAL PRIMARY KEY,
    constant_id INTEGER NOT NULL REFERENCES constant_smiles(id),
    rule_smiles_id INTEGER NOT NULL REFERENCES rule_smiles(id),
    compound_id INTEGER NOT NULL REFERENCES compound (id),
    num_frags INTEGER,
    rule_smiles_num_heavies INTEGER,
    compound_num_heavies INTEGER
)""",
        """
CREATE TABLE to_construct (
    id SERIAL PRIMARY KEY,
    constant_id INTEGER NOT NULL REFERENCES constant_smiles(id),
    rule_smiles_id INTEGER NOT NULL REFERENCES rule_smiles(id),
    compound_id INTEGER NOT NULL REFERENCES compound (id),
    num_frags INTEGER,
    rule_smiles_num_heavies INTEGER,
    compound_num_heavies INTEGER
)""",
        # Optional indexing to improve speed of construct table writing
        "CREATE INDEX pair_c_re_c1_c2 on pair (constant_id, rule_environment_id, compound1_id, compound2_id)",
    ]

    for statement in DDL_statements:
        c.execute(statement)
    conn.commit()

    # Add molecule objects with explicit H to molecule columns
    s = """
SELECT id, clean_smiles
FROM compound
"""
    c.execute(s)
    compound_id_smiles_list = c.fetchall()
    compound_id_pkl_list = []
    for compound_id, smiles in compound_id_smiles_list:
        # The whole reason we're doing all this is to add explicit H to the skeleton of the cartridge molecules,
        # so that the cartridge sss functions work intuitively (i.e. in a similar fashion to reaxys, scifinder, etc.),
        # in all test cases, based on what users would draw in a sketcher
        pkl = Chem.AddHs(Chem.MolFromSmiles(smiles)).ToBinary()
        compound_id_pkl_list.append((compound_id, pkl))

    update_statement = """
UPDATE compound
SET clean_smiles_mol = mol_from_pkl(update_batch.pkl)
FROM (VALUES %s) AS update_batch (id, pkl)
WHERE compound.id = update_batch.id
"""
    execute_values(c, update_statement, compound_id_pkl_list)
    conn.commit()

    s = """
SELECT id, smiles
FROM constant_smiles
"""
    c.execute(s)
    compound_id_smiles_list = c.fetchall()
    compound_id_pkl_list = []
    for compound_id, smiles in compound_id_smiles_list:
        # For constant_smiles and rule_smiles, we use He, Ne, and Ar as dummy atoms to represent attachment points in a way that is respected by the cartridge substructure search
        pkl = Chem.AddHs(Chem.MolFromSmiles(smiles.replace('*:1', 'He').replace('*:2', 'Ne').replace('*:3', 'Ar'))).ToBinary()
        compound_id_pkl_list.append((compound_id, pkl))
    update_statement = """
UPDATE constant_smiles
SET smiles_mol = mol_from_pkl(update_batch.pkl)
FROM (VALUES %s) AS update_batch (id, pkl)
WHERE constant_smiles.id = update_batch.id
"""
    execute_values(c, update_statement, compound_id_pkl_list)
    conn.commit()

    s = """
SELECT id, smiles
FROM rule_smiles
"""
    c.execute(s)
    compound_id_smiles_list = c.fetchall()
    compound_id_pkl_list = []
    for compound_id, smiles in compound_id_smiles_list:
        # For constant_smiles and rule_smiles, we use He, Ne, and Ar as dummy atoms to represent attachment points in a way that is respected by the cartridge substructure search
        pkl = Chem.AddHs(Chem.MolFromSmiles(smiles.replace('*:1', 'He').replace('*:2', 'Ne').replace('*:3', 'Ar'))).ToBinary()
        compound_id_pkl_list.append((compound_id, pkl))

    update_statement = """
UPDATE rule_smiles
SET smiles_mol = mol_from_pkl(update_batch.pkl)
FROM (VALUES %s) AS update_batch (id, pkl)
WHERE rule_smiles.id = update_batch.id
"""
    execute_values(c, update_statement, compound_id_pkl_list)
    conn.commit()

    # FROM constructs
    s = """
INSERT INTO from_construct (constant_id, rule_smiles_id, compound_id)
SELECT DISTINCT pair.constant_id, rule.from_smiles_id, pair.compound1_id
FROM rule, rule_environment, pair
WHERE pair.rule_environment_id = rule_environment.id
AND rule_environment.rule_id = rule.id
"""
    c.execute(s)
    conn.commit()

    # TO constructs
    s = """
INSERT INTO to_construct (constant_id, rule_smiles_id, compound_id)
SELECT DISTINCT pair.constant_id, rule.to_smiles_id, pair.compound2_id
FROM rule, rule_environment, pair
WHERE pair.rule_environment_id = rule_environment.id
AND rule_environment.rule_id = rule.id
"""
    c.execute(s)
    conn.commit()

    # Transfering data to construct tables - redundant, but for potentially faster querying

    construct_update_statements = [
        """
UPDATE from_construct
SET num_frags = (
    SELECT (CHAR_LENGTH(rule_smiles.smiles) - CHAR_LENGTH(REPLACE(rule_smiles.smiles, '*', ''))) / CHAR_LENGTH('*')
    FROM rule_smiles
    WHERE rule_smiles.id = from_construct.rule_smiles_id
)
WHERE EXISTS (
    SELECT 1
    FROM rule_smiles
    WHERE rule_smiles.id = from_construct.rule_smiles_id
)
""",
        """
UPDATE from_construct
SET rule_smiles_num_heavies = (
    SELECT rule_smiles.num_heavies
    FROM rule_smiles
    WHERE rule_smiles.id = from_construct.rule_smiles_id
)
WHERE EXISTS (
    SELECT 1
    FROM rule_smiles
    WHERE rule_smiles.id = from_construct.rule_smiles_id
)
""",
        """
UPDATE from_construct
SET compound_num_heavies = (
    SELECT compound.clean_num_heavies
    FROM compound
    WHERE compound.id = from_construct.compound_id
)
WHERE EXISTS (
    SELECT 1
    FROM compound
    WHERE compound.id = from_construct.compound_id
)
""",
        """
UPDATE to_construct
SET num_frags = (
    SELECT (CHAR_LENGTH(rule_smiles.smiles) - CHAR_LENGTH(REPLACE(rule_smiles.smiles, '*', ''))) / CHAR_LENGTH('*')
    FROM rule_smiles
    WHERE rule_smiles.id = to_construct.rule_smiles_id
)
WHERE EXISTS (
    SELECT 1
    FROM rule_smiles
    WHERE rule_smiles.id = to_construct.rule_smiles_id
)
""",
        """
UPDATE to_construct
SET rule_smiles_num_heavies = (
    SELECT rule_smiles.num_heavies
    FROM rule_smiles
    WHERE rule_smiles.id = to_construct.rule_smiles_id
)
WHERE EXISTS (
    SELECT 1
    FROM rule_smiles
    WHERE rule_smiles.id = to_construct.rule_smiles_id
)
""",
        """
UPDATE to_construct
SET compound_num_heavies = (
    SELECT compound.clean_num_heavies
    FROM compound
    WHERE compound.id = to_construct.compound_id
)
WHERE EXISTS (
    SELECT 1
    FROM compound
    WHERE compound.id = to_construct.compound_id
)
"""
    ]

    for s in construct_update_statements:
        c.execute(s)
        conn.commit()

    # Add functions for structure filtering of rule_smiles and constant_smiles
    # We need to use a plpgsql procedure (below), to force the query planner to expect a certain number of rows to result from the substructure search
    # Currently the row estimate is set to 100,000 for constant_smiles (see last row of the procedure)
    # Otherwise, the planner tends to vastly underestimate how many rows will be returned, causing bad plans such as nested loop joins between large relations
    #  which can, for example, cause a 15 second query (with below code) to take 10+ minutes
    for s in (
        """
CREATE OR REPLACE FUNCTION filter_rule_smiles ( pkl IN BYTEA, nf IN INTEGER, min_heavies IN INTEGER DEFAULT 0, max_heavies IN INTEGER DEFAULT 0) RETURNS TABLE(id INTEGER)
AS $$
DECLARE
BEGIN
    IF (min_heavies > 0 AND max_heavies = 0) THEN
        RETURN QUERY
            SELECT rule_smiles.id AS id
            FROM rule_smiles
            WHERE rule_smiles.smiles_mol@>mol_from_pkl(pkl)
            AND rule_smiles.num_frags = nf
            AND rule_smiles.num_heavies >= min_heavies;
    ELSIF (min_heavies > 0 AND max_heavies > 0) THEN
        RETURN QUERY
            SELECT rule_smiles.id AS id
            FROM rule_smiles
            WHERE rule_smiles.smiles_mol@>mol_from_pkl(pkl)
            AND rule_smiles.num_frags = nf
            AND rule_smiles.num_heavies >= min_heavies
            AND rule_smiles.num_heavies <= max_heavies;
    ELSIF (min_heavies = 0 AND max_heavies > 0) THEN
        RETURN QUERY
            SELECT rule_smiles.id AS id
            FROM rule_smiles
            WHERE rule_smiles.smiles_mol@>mol_from_pkl(pkl)
            AND rule_smiles.num_frags = nf
            AND rule_smiles.num_heavies <= max_heavies;
    ELSIF (min_heavies = 0 AND max_heavies = 0) THEN
        RETURN QUERY
            SELECT rule_smiles.id AS id
            FROM rule_smiles
            WHERE rule_smiles.smiles_mol@>mol_from_pkl(pkl)
            AND rule_smiles.num_frags = nf;
    END IF;
END;
$$ LANGUAGE plpgsql
   ROWS 100
""",
        """
CREATE OR REPLACE FUNCTION filter_constant_smiles ( pkl IN BYTEA, nf IN INTEGER ) RETURNS TABLE(id integer)
AS $$
DECLARE
BEGIN
  RETURN QUERY
  SELECT constant_smiles.id AS id
  FROM constant_smiles
  WHERE constant_smiles.smiles_mol@>mol_from_pkl(pkl)
  AND constant_smiles.num_frags = nf;
END;
$$ LANGUAGE plpgsql
   ROWS 100000
"""):
        c.execute(s)
        conn.commit()

    DDL_statements = [
        """
CREATE TABLE query (
    id SERIAL PRIMARY KEY, 
    inserted_at TIMESTAMP
)""",
        """
CREATE TABLE query_result (
    id SERIAL PRIMARY KEY,
    query_id INTEGER REFERENCES query(id) ON UPDATE CASCADE ON DELETE CASCADE,
    rule_id INTEGER REFERENCES rule(id),
    use_original_direction BOOLEAN,
    from_construct_id INTEGER REFERENCES from_construct(id),
    to_construct_id INTEGER REFERENCES to_construct(id),
    environment_smarts VARCHAR(4000),
    from_smiles_env VARCHAR(4000),
    to_smiles_env VARCHAR(4000)
)""",
        """
CREATE TABLE sketched_content (
    id SERIAL PRIMARY KEY,
    molfile TEXT,
    variable_atoms VARCHAR(4000),
    variable_bonds VARCHAR(4000),
    environment_atoms VARCHAR(4000),
    environment_bonds VARCHAR(4000)
)""",
    # Snapquery captures an input state, to allow for reproducible and shareable queries
        """
CREATE TABLE snapquery (
    id SERIAL PRIMARY KEY,
    version_id INTEGER,
    query_id INTEGER REFERENCES query(id) ON UPDATE CASCADE ON DELETE SET NULL,
    query_type VARCHAR(4000),
    transform_order VARCHAR(4000),
    mol1_sketched_content_id INTEGER REFERENCES sketched_content(id),
    mol2_sketched_content_id INTEGER REFERENCES sketched_content(id),
    REQUIRED_properties VARCHAR(4000),
    OPTIONAL_properties VARCHAR(4000),
    variable_min_heavies INTEGER,
    variable_max_heavies INTEGER,
    compound_min_heavies INTEGER,
    compound_max_heavies INTEGER,
    aggregation_type VARCHAR(4000)
)""",
    # Snapfilter captures the filter state of an output, to allow for sharing curated analyses
        """
CREATE TABLE snapfilter (
    id SERIAL PRIMARY KEY,
    version_id INTEGER,
    snapfilter_string TEXT
)""",
        """
CREATE TABLE snapshot (
    id SERIAL PRIMARY KEY,
    snapquery_id INTEGER REFERENCES snapquery(id),
    snapfilter_id INTEGER REFERENCES snapfilter(id)
)""",

        # We use num_frags in rule_smiles / constant_smiles in partial indexes for substructure searches, making them faster
        "ALTER TABLE rule_smiles ADD num_frags INTEGER",
        "ALTER TABLE constant_smiles ADD num_frags INTEGER",
        "UPDATE constant_smiles SET num_frags = (CHAR_LENGTH(constant_smiles.smiles) - CHAR_LENGTH(REPLACE(constant_smiles.smiles, '*', ''))) / CHAR_LENGTH('*')",
        "UPDATE rule_smiles SET num_frags = (CHAR_LENGTH(rule_smiles.smiles) - CHAR_LENGTH(REPLACE(rule_smiles.smiles, '*', ''))) / CHAR_LENGTH('*')",

        # We store base (i.e. log, negative log, or raw), units (e.g. M), and displayed names of properties here
        # API endpoints have been written to fetch values in certain formats, e.g. get the hERG IC50 in uM, but this property may have been loaded in as negative log of molarity
        # By definining and populating these columns, we will know what to do with respect to transforming the data to the desired format
        "ALTER TABLE property_name ADD base VARCHAR(4000)",
        "ALTER TABLE property_name ADD unit VARCHAR(4000)",
        "ALTER TABLE property_name ADD display_name VARCHAR(4000)",
        "ALTER TABLE property_name ADD display_base VARCHAR(4000)",
        "ALTER TABLE property_name ADD display_unit VARCHAR(4000)",
        "ALTER TABLE property_name ADD change_displayed VARCHAR(4000)",
		"""ALTER TABLE rule_environment_statistics ADD COLUMN query_id INT 
CONSTRAINT restats_query_fk_id REFERENCES query(id)
ON UPDATE CASCADE ON DELETE SET NULL""",

        "CREATE INDEX compound_clean_smiles_mol_idx ON compound USING gist(clean_smiles_mol)",
        # Use partial indexes for faster searching, because we know the num_frags (number of bond cuts in our query) ahead of the search
        "CREATE INDEX constant_smiles_smiles_mol_idx_nf1 ON constant_smiles USING gist (smiles_mol) WHERE (num_frags = 1)",
        "CREATE INDEX constant_smiles_smiles_mol_idx_nf2 ON constant_smiles USING gist (smiles_mol) WHERE (num_frags = 2)",
        "CREATE INDEX constant_smiles_smiles_mol_idx_nf3 ON constant_smiles USING gist (smiles_mol) WHERE (num_frags = 3)",
        "CREATE INDEX rule_smiles_smiles_mol_idx_nf1 ON rule_smiles USING gist (smiles_mol) WHERE (num_frags = 1)",
        "CREATE INDEX rule_smiles_smiles_mol_idx_nf2 ON rule_smiles USING gist (smiles_mol) WHERE (num_frags = 2)",
        "CREATE INDEX rule_smiles_smiles_mol_idx_nf3 ON rule_smiles USING gist (smiles_mol) WHERE (num_frags = 3)",

        "CREATE INDEX rule_smiles_smiles_mol_idx ON rule_smiles USING gist(smiles_mol)",
        # We always start the construct indices with num_frags
        # Then we allow filtering starting either with constant_id or rule_smiles_id
        "CREATE INDEX from_const_nf_c_rs_cmpd on from_construct (num_frags, constant_id, rule_smiles_id, compound_id)",
        "CREATE INDEX from_const_nf_rs_c_cmpd on from_construct (num_frags, rule_smiles_id, constant_id, compound_id)",
        "CREATE INDEX to_const_nf_c_rs_cmpd on to_construct (num_frags, constant_id, rule_smiles_id, compound_id)",
        "CREATE INDEX to_const_nf_rs_c_cmpd on to_construct (num_frags, rule_smiles_id, constant_id, compound_id)",
        "CREATE INDEX cprop_pnid_cid_v on compound_property (property_name_id, compound_id, value)",
        "CREATE INDEX rule_from_to_id on rule (from_smiles_id, to_smiles_id, id)",
        "CREATE INDEX rule_to_from_id on rule (to_smiles_id, from_smiles_id, id)",
        "CREATE INDEX property_name_name on property_name (name)",

        "CREATE INDEX query_result_q_r ON query_result USING btree (query_id, rule_id)",
        "CREATE INDEX query_result_queryid_toenv_fromenv ON query_result(query_id, to_smiles_env, from_smiles_env)",
        "CREATE INDEX query_result_queryid_fromenv_toenv ON query_result(query_id, from_smiles_env, to_smiles_env)",

        # These indexes on from/to constructs are experimental, may not be needed if they are not used by the optimizer / not faster even if used:
        "CREATE INDEX from_const_nf_rsNH_c_rs_cmpd on from_construct (num_frags, rule_smiles_num_heavies, constant_id, rule_smiles_id, compound_id)",
        "CREATE INDEX from_const_nf_rsNH_rs_c_cmpd on from_construct (num_frags, rule_smiles_num_heavies, rule_smiles_id, constant_id, compound_id)",
        "CREATE INDEX from_const_nf_cmpdNH_c_rs_cmpd on from_construct (num_frags, compound_num_heavies, constant_id, rule_smiles_id, compound_id)",
        "CREATE INDEX from_const_nf_cmpdNH_rs_c_cmpd on from_construct (num_frags, compound_num_heavies, rule_smiles_id, constant_id, compound_id)",
        "CREATE INDEX from_const_nf_rsNH_cmpdNH_c_rs_cmpd on from_construct (num_frags, rule_smiles_num_heavies, compound_num_heavies, constant_id, rule_smiles_id, compound_id)",
        "CREATE INDEX from_const_nf_rsNH_cmpdNH_rs_c_cmpd on from_construct (num_frags, rule_smiles_num_heavies, compound_num_heavies, rule_smiles_id, constant_id, compound_id)",
    
        "CREATE INDEX to_const_nf_rsNH_c_rs_cmpd on to_construct (num_frags, rule_smiles_num_heavies, constant_id, rule_smiles_id, compound_id)",
        "CREATE INDEX to_const_nf_rsNH_rs_c_cmpd on to_construct (num_frags, rule_smiles_num_heavies, rule_smiles_id, constant_id, compound_id)",
        "CREATE INDEX to_const_nf_cmpdNH_c_rs_cmpd on to_construct (num_frags, compound_num_heavies, constant_id, rule_smiles_id, compound_id)",
        "CREATE INDEX to_const_nf_cmpdNH_rs_c_cmpd on to_construct (num_frags, compound_num_heavies, rule_smiles_id, constant_id, compound_id)",
        "CREATE INDEX to_const_nf_rsNH_cmpdNH_c_rs_cmpd on to_construct (num_frags, rule_smiles_num_heavies, compound_num_heavies, constant_id, rule_smiles_id, compound_id)",
        "CREATE INDEX to_const_nf_rsNH_cmpdNH_rs_c_cmpd on to_construct (num_frags, rule_smiles_num_heavies, compound_num_heavies, rule_smiles_id, constant_id, compound_id)",

        "ANALYZE",
    ]

    for statement in DDL_statements:
        c.execute(statement)
        conn.commit()

    return