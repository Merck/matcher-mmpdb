# mmpdb - matched molecular pair database generation and analysis
#
# Copyright (c) 2015-2017, F. Hoffmann-La Roche Ltd.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#    * Redistributions of source code must retain the above copyright
#      notice, this list of conditions and the following disclaimer.
#    * Redistributions in binary form must reproduce the above
#      copyright notice, this list of conditions and the following
#      disclaimer in the documentation and/or other materials provided
#      with the distribution.
#    * Neither the name of F. Hoffmann-La Roche Ltd. nor the names of
#      its contributors may be used to endorse or promote products
#      derived from this software without specific prior written
#      permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#

from __future__ import print_function, absolute_import, division

import sys
import os.path
import datetime
import json
import itertools
from collections import OrderedDict, defaultdict

from . import reporters
from . import config
from ._compat import basestring

from . import peewee
from psycopg2.extras import execute_values

SCHEMA_FILENAME = os.path.join(os.path.dirname(__file__), "schema.sql")
CREATE_INDEX_FILENAME = os.path.join(os.path.dirname(__file__), "create_index.sql")
DROP_INDEX_FILENAME = os.path.join(os.path.dirname(__file__), "drop_index.sql")


_schema_template = None    
def get_schema_template():
    global _schema_template
    if _schema_template is None:
        with open(SCHEMA_FILENAME) as infile:
            _schema_template = infile.read()
    return _schema_template


_create_index_sql = None    
def get_create_index_sql():
    global _create_index_sql
    if _create_index_sql is None:
        with open(CREATE_INDEX_FILENAME) as infile:
            _create_index_sql = infile.read()
    return _create_index_sql


_create_oracle_index_sql = None    
def get_create_oracle_index_sql():
    global _create_oracle_index_sql
    if _create_oracle_index_sql is None:
        with open(CREATE_INDEX_FILENAME) as infile:
            _create_oracle_index_sql = infile.read()
            _create_oracle_index_sql += "-- OracleDB"
    return _create_oracle_index_sql

_drop_index_sql = None    
def get_drop_index_sql():
    global _drop_index_sql
    if _drop_index_sql is None:
        with open(DROP_INDEX_FILENAME) as infile:
            _drop_index_sql = infile.read()
    return _drop_index_sql


class SQLiteConfig:
    PRIMARY_KEY = "INTEGER PRIMARY KEY"
    COLLATE = ""     # default collation is binary so this isn't needed
    DATE = "DATETIME"

class MySQLConfig:
    PRIMARY_KEY = "INTEGER AUTO_INCREMENT"
    COLLATE = "COLLATE latin_bin"  # default is case-insensitive; force binary
    DATE = "DATETIME"

class PostgresConfig:
    PRIMARY_KEY = "SERIAL PRIMARY KEY"
    COLLATE = ""
    DATE = "TIMESTAMP"

class OracleConfig:
    PRIMARY_KEY = "INTEGER PRIMARY KEY"
    COLLATE = ""
    DATE = "TIMESTAMP"
    VARCHAR = "VARCHAR2(4000)"
    
    # Unlike SQLite, Oracle 11.2 does not automatically generate/insert
    # auto-incremented primary key values,
    # when values are inserted into other columns of a table.
    # Therefore, the below sequences and triggers are required,
    # in order for the loadprops_command function in do_database module to work.
    s1 = """
    CREATE SEQUENCE cmpd_prop_insert
    """

    s2 = """
    CREATE OR REPLACE TRIGGER cmpd_prop_on_insert
        BEFORE INSERT ON compound_property
        FOR EACH ROW
    BEGIN
        SELECT cmpd_prop_insert.nextval
        INTO :new.id
        FROM dual;
    END;
    """

    s3 = """
    CREATE SEQUENCE re_stats_insert
    """

    s4 = """
    CREATE OR REPLACE TRIGGER re_stats_on_insert
        BEFORE INSERT ON rule_environment_statistics
        FOR EACH ROW
    BEGIN
        SELECT re_stats_insert.nextval
        INTO :new.id
        FROM dual;
    END;
    """
    
    AUTO_INC_statements = [s1,s2,s3,s4]

def get_schema_for_database(db_config):
    template = get_schema_template()
    
    # Handle some non-portable SQL
    text = template.replace("$PRIMARY_KEY$", db_config.PRIMARY_KEY)
    schema = text.replace("$COLLATE$", db_config.COLLATE)
    schema2 = schema.replace("$DATE$", db_config.DATE)
    
    if db_config == OracleConfig:
        schema3 = schema2.replace("VARCHAR(255)", db_config.VARCHAR)
        schema4 = schema3.replace("VARCHAR(2000)", db_config.VARCHAR)
        schema5 = schema4.replace("VARCHAR(43)", db_config.VARCHAR)
        # add "-- OracleDB" tag so that the DB is recognized as an oracle DB by other functions
        #schema5 += db_config.AUTO_INC_PK + "-- OracleDB"
        schema5 += "-- OracleDB"
        return schema5

    return schema2


def _get_sql_statements(bulk_sql):
    # Remove the license preamble, if it exists
    i = bulk_sql.find("END OF LICENSE")
    if i != -1:
        newline_i = bulk_sql.index("\n",i)
        bulk_sql = bulk_sql[newline_i+1:]

    # Need to modify character appearing at end of SQL statements to be compatible with cx_Oracle
    separator = ""
    j = bulk_sql.find("-- OracleDB")
    if j == -1:
        separator = ";"
    # Need to delete -- OracleDB tag to prevent Oracle SQL error
    else:
        bulk_sql = bulk_sql[0:-11]

    blocks = bulk_sql.split(";")
    # Can only execute one statement at a time,

    for statement in blocks:
        statement = statement.strip()
        if not statement:
            continue
        statement += separator
        #print("execute", repr(statement))
        yield statement


def _execute_sql(c, bulk_sql):
    for statement in _get_sql_statements(bulk_sql):
        try:            
            c.execute(statement)
        except Exception:
            sys.stderr.write("Failed to execute the following SQL:\n")
            sys.stderr.write(statement)
            raise


def create_schema_for_sqlite(sqlite_db):
    c = sqlite_db.cursor()
    _execute_sql(c,
                 get_schema_for_database(SQLiteConfig))

def create_schema_for_oracle(oracle_db):
    c = oracle_db.cursor()
    _execute_sql(c,
                 get_schema_for_database(OracleConfig))
    for statement in OracleConfig.AUTO_INC_statements:
        c.execute(statement)

def create_schema_for_postgres(postgres_db):
    c = postgres_db.cursor()
    _execute_sql(c,
                 get_schema_for_database(PostgresConfig))

def create_index(c):
    _execute_sql(c, get_create_index_sql())

def create_oracle_index(c):
    _execute_sql(c, get_create_oracle_index_sql())

def drop_index(c):
    _execute_sql(c, get_drop_index_sql())


INSERT_RULE_PROPERTY_SQL = (
    "insert into rule_property (rule_id, property_name_id, count, avg, std, kurtosis, "
    " skewness, min, q1, median, q3, max, paired_t, p_value) values (?,?,?,?,?,?,?,?,?,?,?,?,?,?)")
    
UPDATE_RULE_PROPERTY_SQL = (
    "update rule_property set count=?, avg=?, std=?, kurtosis=?, "
    " skewness=?, min=?, q1=?, median=?, q3=?, max=?, paired_t=?, p_value=? where id=?")

########## Work with an existing dataset via an ORM-like API

# Note: an earlier system stored multiple datasets in the same
# database.  This turned out to be complicated because all of the
# tables needed an extra dataset identifier. It also made it slow to
# load a dataset because I couldn't drop the indices, load, and build
# the indices without taking all of the datasets offline.

# We decided to get rid of MySQL/multi-dataset support and instead
# store each dataset in its own SQLite file, which can be indexed
# independently. I could therefore merge the MMPDatabase and
# MMPDataset into one class, but part of me thinks it's nicer to keep
# them distinct.
    

class PropertyNameRow(object):
    def __init__(self, id, name):
        self.id = id
        self.name = name


class TableSizes(object):
    def __init__(self, num_compounds, num_rules, num_pairs,
                 num_rule_environments, num_rule_environment_stats):
        self.num_compounds = num_compounds
        self.num_rules = num_rules
        self.num_pairs = num_pairs
        self.num_rule_environments = num_rule_environments
        self.num_rule_environment_stats = num_rule_environment_stats
    
        
class MMPDatabase(object):
    def __init__(self, db):
        self.db = db
        self._dataset = None

    def SQL(self, sql):
        return sql.replace("?", self.db.interpolation)
    
    def execute(self, sql, args=(), cursor=None):
        if cursor is None:
            cursor = self.db.get_cursor()
        sql = sql.replace("?", self.db.interpolation)
        
        # for Oracle DB with cx_Oracle
        # and Postgres DB with psycopg2
        try:
            if self.db.as_statement is not None:
                sql = sql.replace(" as ", self.db.as_statement)
                # to accelerate large select statements with oracle
                if type(self.db) == peewee.CustomPostgresqlDatabase:
                    cursor.itersize = 100000
                else:
                    try:
                        cursor.prefetchrows = 100000
                        cursor.arraysize = 100000
                    except:
                        pass
        except:
            pass

        if 0:
            import time
            print("EXECUTE")
            print(sql)
            print(repr(args))
            t1 = time.time()
        cursor.execute(sql, args)
        if 0:
            t2 = time.time()
            print("Elapsed1:", t2-t1)
            ## if "SELECT" in sql:
            ##     n = len(list(cursor))
            ##     t3 = time.time()
            ##     print("Elapsed2:", t3-t1, "for", n, "rows")
            ##     cursor.execute(sql, args)
        return cursor

    def atomic(self):
        return self.db.atomic()

    def get_cursor(self, cursor=None):
        if cursor is None:
            return self.db.get_cursor()
        return cursor

    def get_dataset(self, cursor=None):
        if self._dataset is None:
            cursor = self.get_cursor(cursor)
            cursor.execute("SELECT id, mmpdb_version, title, creation_date, "
                           "fragment_options, index_options, is_symmetric, "
                           "num_compounds, num_rules, num_pairs, num_rule_environments, "
                           "num_rule_environment_stats FROM dataset WHERE id = 1")
            for (id, mmpdb_version, title, creation_date,
                 fragment_options, index_options, is_symmetric,
                 num_compounds, num_rules, num_pairs, num_rule_environments,
                 num_rule_environment_stats) in cursor:

                # The MySQL adapter returns a datetime object. The SQLite adapter
                # returns a string. This bit of hack converts the datetime string
                # a real datetime object.
                if isinstance(creation_date, basestring):
                    creation_date = creation_date.rstrip("Z")
                    creation_date = datetime.datetime.strptime(creation_date, "%Y-%m-%d %H:%M:%S.%f")

                table_sizes = TableSizes(num_compounds, num_rules, num_pairs,
                                         num_rule_environments, num_rule_environment_stats)
                                         
                self._dataset = MMPDataset(self, title, creation_date,
                                           fragment_options, index_options, is_symmetric, table_sizes)
                break
            else:
                raise AssertionError("Database missing dataset with id 1")
        return self._dataset


class Pair(object):
    def __init__(self, pair_id, rule_environment_id, compound1_id, compound2_id, constant_id):
        self.pair_id = pair_id
        self.rule_environment_id = rule_environment_id
        self.compound1_id = compound1_id
        self.compound2_id = compound2_id
        self.constant_id = constant_id


def _get_one(cursor):
    for row in cursor:
        return row[0]
    raise AssertionError("missing result")


def _get_one_or_none(cursor):
    for row in cursor:
        return row[0]
    return None


class MMPDataset(object):
    def __init__(self, mmpa_db, title, creation_date,
                 fragment_options_str, index_options_str,
                 is_symmetric, table_sizes):
        self.mmpa_db = mmpa_db
        self.title = title
        self.creation_date = creation_date
        self.fragment_options_str = fragment_options_str
        self.index_options_str = index_options_str
        self.is_symmetric = is_symmetric
        self.table_sizes = table_sizes
        self._rule_ids = None
        self._clean_smiles_cache = {}
        self._property_name_table = None
        
        if type(self.mmpa_db.db) in [peewee.OracleDatabase, peewee.CustomPostgresqlDatabase]:
            self.re_stats_updates = []
            self.re_stats_inserts = []
            self.add_rule_environment_statistics = self.add_rule_environment_statistics_toList
            self.update_rule_environment_statistics = self.update_rule_environment_statistics_toList
            if type(self.mmpa_db.db) == peewee.OracleDatabase:
                self.batch_update_re_stats = self.batch_update_re_stats_to_oracle
                self.batch_insert_re_stats = self.batch_insert_re_stats_to_oracle
                self.validate_statistics = self.validate_oracle_statistics
            elif type(self.mmpa_db.db) == peewee.CustomPostgresqlDatabase:
                self.batch_update_re_stats = self.batch_update_re_stats_to_postgres
                self.batch_insert_re_stats = self.batch_insert_re_stats_to_postgres
                self.validate_statistics = self.validate_postgres_statistics
        else:
            self.add_rule_environment_statistics = self.add_rule_environment_statistics_toDB
            self.update_rule_environment_statistics = self.update_rule_environment_statistics_toDB
            
    def get_cursor(self, cursor=None):
        return self.mmpa_db.get_cursor(cursor)

    def get_table_sizes(self, recount=False, cursor=None):
        if not recount:
            return self.table_sizes
        cursor = self.get_cursor(cursor)
        return TableSizes(
            num_compounds = self.get_num_compounds(cursor),
            num_rules = self.get_num_rules(cursor),
            num_pairs = self.get_num_pairs(cursor),
            num_rule_environments = self.get_num_rule_environments(cursor),
            num_rule_environment_stats = self.get_num_rule_environment_stats(cursor),
            )

    def get_property_name_rows(self, cursor=None):
        c = self.mmpa_db.execute("SELECT id, name from property_name ORDER BY id",
                                 cursor=cursor)
        return [PropertyNameRow(id, name) for (id, name) in c]
    
    def get_property_names(self, cursor=None):
        c = self.mmpa_db.execute("SELECT name from property_name ORDER BY id",
                                 cursor=cursor)
        return [name for (name,) in c]

    def get_property_name_id(self, name, default=None):
        return self.get_property_names_table().get(name, default)

    def get_property_names_table(self, cursor=None):
        table = self._property_name_table
        if table is None:
            table = OrderedDict()
            for propname_row in self.get_property_name_rows(cursor):
                table[propname_row.name] = propname_row.id
            self._property_name_table = table
        return table

    def get_num_compounds(self, cursor=None):
        return _get_one(self.mmpa_db.execute("SELECT count(*) from compound", cursor=cursor))

    def get_num_compound_properties(self, property_name_id, cursor=None):
        return _get_one(self.mmpa_db.execute(
            "SELECT count(*) from compound_property "
            " WHERE property_name_id = ?",
            (property_name_id,),
            cursor=cursor))

    def get_num_rules(self, cursor=None):
        return _get_one(self.mmpa_db.execute("SELECT count(*) from rule", cursor=cursor))
    
    def get_num_pairs(self, property_name_id=None, cursor=None):
        return _get_one(self.mmpa_db.execute("SELECT count(*) FROM pair", cursor=cursor))

    def get_num_rule_environments(self, cursor=None):
        return _get_one(self.mmpa_db.execute("SELECT count(*) from rule_environment", cursor=cursor))
    
    def get_num_rule_environment_stats(self, cursor=None):
        return _get_one(self.mmpa_db.execute("SELECT count(*) from rule_environment_statistics", cursor=cursor))
    
    def get_num_rule_smiles(self, cursor=None):
        return _get_one(self.mmpa_db.execute("SELECT count(*) from rule_smiles", cursor=cursor))
    
    def get_num_constant_smiles(self, cursor=None):
        return _get_one(self.mmpa_db.execute("SELECT count(*) from constant_smiles", cursor=cursor))
        
    def get_property_names_and_counts(self, cursor=None):
        c = self.mmpa_db.execute("""
SELECT property_name.name, count(property_name_id)
  FROM compound_property, property_name
 WHERE compound_property.property_name_id = property_name.id
 GROUP BY property_name_id
 ORDER BY property_name.name""", cursor=cursor)
        return list(c)

    def get_public_id_to_id_table(self, cursor=None):
        c = self.mmpa_db.execute("SELECT public_id, id FROM compound", cursor=cursor)
        return dict(c)

    def get_or_add_property_name(self, property_name, cursor=None):
        cursor = self.get_cursor(cursor)
        table = self.get_property_names_table(cursor=cursor)
        property_name_id = table.get(property_name, None)
        if property_name_id is None:
            # Get the new id for the property name
            if table:
                property_name_id = max(table.values())+1
            else:
                property_name_id = 0
            # Insert
            self.mmpa_db.execute(
                "INSERT INTO property_name (id, name) VALUES (?, ?)",
                (property_name_id, property_name), cursor=cursor)
            # Add to my local cache
            self._property_name_table[property_name] = property_name_id
        return property_name_id
    
    def get_fragment_options(self, cursor=None):
        from . import fragment_io
        cursor = self.mmpa_db.execute("SELECT fragment_options FROM dataset WHERE id = 1", (),
                                      cursor)
        for (fragment_options,) in cursor:
            import json
            d = json.loads(fragment_options)
            # {"max_heavies": 100, "max_rotatable_bonds": 10,
            #  "rotatable_smarts": "[!$([NH]!@C(=O))&!D1&!$(*#*)]-&!@[!$([NH]!@C(=O))&!D1&!$(*#*)]",
            #  "cut_smarts": "[#6+0;!$(*=,#[!#6])]!@!=!#[!#0;!#1;!$([CH2]);!$([CH3][CH2])]", "num_cuts": 3,
            #   "method": "dalke", "salt_remover": "<default>"
            return config.FragmentOptions(
                max_heavies = d["max_heavies"],
                max_rotatable_bonds = d["max_rotatable_bonds"],
                rotatable_smarts = str(d["rotatable_smarts"]),  # otherwise it's unicode
                cut_smarts = str(d["cut_smarts"]),   # otherwise it's unicode
                num_cuts = d["num_cuts"],
                salt_remover = d["salt_remover"],
                method = d["method"],
                min_heavies_per_const_frag = d["min_heavies_per_const_frag"]
                )
        raise AssertionError("dataset 1 is supposed to exist")

    def get_rule_smiles_id(self, smiles, cursor=None):
        c = self.mmpa_db.execute(
            "SELECT id FROM rule_smiles WHERE smiles = ?",
            (smiles,), cursor)
        return _get_one_or_none(c)

    def find_rule_environments_for_transform(self, smiles_id, possible_env_fps, max_variable_size=9999,
                                             cursor=None):
        assert len(possible_env_fps) > 0, possible_env_fps
        cursor = self.mmpa_db.get_cursor(cursor)

        test_fp_in = " OR ".join(("environment_fingerprint.fingerprint = ?",)*len(possible_env_fps))
        
        execute_args = (smiles_id, ) + tuple(possible_env_fps) + (max_variable_size,)

        # Find the rule environments which use this SMILES on the LHS
        sql = ("SELECT rule_environment.rule_id, rule_environment.id, rule_smiles.smiles, 0\n"  # 0 = forward
               "  FROM rule_environment, rule, environment_fingerprint, rule_smiles\n"
               " WHERE rule_environment.rule_id = rule.id\n"
               "   AND rule.from_smiles_id = ?\n"  # NOTE: *from*_smiles_id
               "   AND rule_environment.environment_fingerprint_id = environment_fingerprint.id\n"
               "   AND (" + test_fp_in + ")\n"
               "   AND rule.to_smiles_id = rule_smiles.id"  # NOTE: *to*_smiles_id
               "   AND rule_smiles.num_heavies <= ?\n"
               )
        if not self.is_symmetric:
            sql += (
               "UNION\n"
               "SELECT rule_environment.rule_id, rule_environment.id, rule_smiles.smiles, 1\n"  # 1 = reverse
               "  FROM rule_environment, rule, environment_fingerprint, rule_smiles\n"
               " WHERE rule_environment.rule_id = rule.id\n"
               "   AND rule.to_smiles_id = ?\n"   # NOTE: *to*_smiles_id
               "   AND rule_environment.environment_fingerprint_id = environment_fingerprint.id\n"
               "   AND (" + test_fp_in + ")\n"
               "   AND rule.from_smiles_id = rule_smiles.id\n"  # NOTE: *from*_smiles_id
               "   AND rule_smiles.num_heavies <= ?"
               )
            execute_args = execute_args + execute_args  # Double the args, one for each direction

        #print("SQL", sql, execute_args)
        for x in self.mmpa_db.execute(sql, execute_args, cursor):
            yield x
        
    def get_fingerprint_ids(self, fingerprints, cursor=None):
        fpids = set()
        if fingerprints:
            in_term = "', '".join(fingerprints)
            c = self.mmpa_db.execute("SELECT id FROM environment_fingerprint "
                                     "WHERE fingerprint IN ('""" + in_term + "')",
                                     cursor=cursor)
            fpids.update(fpid for (fpid,) in c)
            
        return fpids
            
    def iter_selected_property_rules(self, from_smiles, to_smiles, property_id, cursor=None):
        assert from_smiles is not None
        assert to_smiles is not None
        sql = """
SELECT rule.id, 0, from_smiles.smiles, from_smiles.num_heavies, to_smiles.smiles, to_smiles.num_heavies,
        rule_environment.id, rule_environment.radius,
        rule_environment.environment_fingerprint_id, environment_fingerprint.fingerprint,
        rule_environment_statistics.id, count, avg, std, kurtosis, skewness, min, q1, median, q3, max, paired_t, p_value
  FROM rule, rule_environment, environment_fingerprint, rule_environment_statistics,
          rule_smiles as from_smiles, rule_smiles as to_smiles
 WHERE from_smiles.smiles = ?
   AND to_smiles.smiles = ?
   AND rule.from_smiles_id = from_smiles.id
   AND rule.to_smiles_id = to_smiles.id
   AND rule_environment.rule_id = rule.id
   AND rule_environment.environment_fingerprint_id = environment_fingerprint.id
   AND rule_environment_statistics.rule_environment_id = rule_environment.id
   AND rule_environment_statistics.property_name_id = ?
        """
        args = (from_smiles, to_smiles, property_id)
        
        if not self.is_symmetric:
            sql += """\
UNION
SELECT rule.id, 1, to_smiles.smiles, to_smiles.num_heavies, from_smiles.smiles, from_smiles.num_heavies,
        rule_environment.id, rule_environment.radius,
        rule_environment.environment_fingerprint_id, environment_fingerprint.fingerprint,
        rule_environment_statistics.id, count, -avg, std, kurtosis, skewness, -max, -q3, -median, -q1, -min, paired_t, p_value
  FROM rule, rule_environment, environment_fingerprint, rule_environment_statistics, 
        rule_smiles as from_smiles, rule_smiles as to_smiles
 WHERE from_smiles.smiles = ?
   AND to_smiles.smiles = ?
   AND rule.from_smiles_id = from_smiles.id
   AND rule.to_smiles_id = to_smiles.id
   AND rule_environment.rule_id = rule.id
   AND rule_environment.environment_fingerprint_id = environment_fingerprint.id
   AND rule_environment_statistics.rule_environment_id = rule_environment.id
   AND rule_environment_statistics.property_name_id = ?
"""
            args = args + (to_smiles, from_smiles, property_id)

        ## print(sql)
        ## print((from_smiles, to_smiles, property_id,
        ##                                to_smiles, from_smiles, property_id))
        # If this query is slow, and the indices are present, then you might
        # try running "analyze" in the SQLite terminal.
        c = self.mmpa_db.execute(sql, args, cursor=cursor)
        for (rule_id, is_reversed, from_smiles, from_num_heavies, to_smiles, to_num_heavies,
                 rule_environment_id, radius, fingerprint_id, fingerprint,
                 rule_environment_statistics_id, count, avg, std, kurtosis, skewness, min, q1, median, q3, max, paired_t, p_value) in c:
            yield PropertyRule(
                rule_id, is_reversed, from_smiles, from_num_heavies, to_smiles, to_num_heavies,
                rule_environment_id, radius, fingerprint_id, fingerprint,
                rule_environment_statistics_id, count, avg, std, kurtosis, skewness, min, q1, median, q3, max, paired_t, p_value)
        
    def get_property_values(self, property_name_id, cursor=None):
        c = self.mmpa_db.execute(
            "SELECT compound_id, value FROM compound_property "
            "   WHERE property_name_id = ?",
            (property_name_id,), cursor=cursor)
        return dict(c)

    def get_property_rule(self, property_name_id, rule_environment_id, is_reversed, cursor=None):
        c = self.mmpa_db.execute(
            "SELECT rule_environment.rule_id, from_smiles.smiles, from_smiles.num_heavies, to_smiles.smiles, to_smiles.num_heavies, "
            "          rule_environment.radius, rule_environment.environment_fingerprint_id, "
            "          environment_fingerprint.fingerprint, "
            "          rule_environment_statistics.id, count, avg, std, kurtosis, skewness, min, q1, median, q3, max, paired_t, p_value "
            "  FROM rule, rule_environment, environment_fingerprint, rule_environment_statistics, "
            "          rule_smiles as from_smiles, rule_smiles as to_smiles "
            " WHERE rule_environment.id = ? "
            "   AND rule_environment_statistics.rule_environment_id = ? "
            "   AND rule_environment_statistics.property_name_id = ? "
            "   AND rule_environment.rule_id = rule.id "
            "   AND rule_environment.environment_fingerprint_id = environment_fingerprint.id "
            "   AND rule.from_smiles_id = from_smiles.id "
            "   AND rule.to_smiles_id = to_smiles.id ",
            (rule_environment_id, rule_environment_id, property_name_id), cursor=cursor)
        rows = list(c)
        if not rows:
            # This can happen if the rule environment exists but there are no pairs for
            # the given property
            return None
        assert len(rows) == 1, (property_name_id, rule_environment_id, is_reversed, rows)
        
        rule_id, from_smiles, from_num_heavies, to_smiles, to_num_heavies, radius, fingerprint_id, fingerprint, \
             rule_environment_statistics_id, count, avg, std, kurtosis, skewness, min, q1, median, q3, max, paired_t, p_value = rows[0]
        if is_reversed:
            from_smiles, to_smiles = to_smiles, from_smiles
            from_num_heavies, to_num_heavies = to_num_heavies, from_num_heavies
            avg = -avg
            min, q1, median, q3, max = -max, -q3, -median, -q1, -min

        return PropertyRule(
            rule_id, is_reversed, from_smiles, from_num_heavies, to_smiles, to_num_heavies,
            rule_environment_id, radius, fingerprint_id, fingerprint,
            rule_environment_statistics_id, count, avg, std, kurtosis, skewness, min, q1, median, q3, max, paired_t, p_value)

    def get_rule_environment_statistics_mapping(self, property_name_ids=None, cursor=None):
        if property_name_ids is None:
            # Get all of them
            c = self.mmpa_db.execute(
                "SELECT rule_environment_id, property_name_id, id "
                "  FROM rule_environment_statistics",
                cursor=cursor)
        elif not property_name_ids:
            # Get none of them
            return {}
        else:
            # Get the selected rule ids
            ids_str = ",".join("%d" % (id,) for id in property_name_ids)
            c = self.mmpa_db.execute(
                "SELECT rule_environment_id, property_name_id, id "
                "  FROM rule_environment_statistics "
                " WHERE property_name_id IN (%s)" % (ids_str,),
                cursor=cursor)
            
        d = {}
        for rule_environment_id, property_name_id, rule_environment_statistics_id in c:
            d[rule_environment_id, property_name_id] = rule_environment_statistics_id
        return d

    # Adding executemany functionality for Oracle DB
    # When it's an Oracle DB, instead of loading everything row by row to the DB,
    # instead accumulate the data into lists, then write the lists (in batches) to DB

    # According to oracle 11.2 documentation:
    # The REAL data type is a floating-point number with a binary precision of 63, or 18 decimal.
    # Based on empirical testing, the below error is thrown when a number is > the large number below, or < the small number below:
    # error cx_Oracle.DatabaseError: DPI-1044: value cannot be represented as an Oracle number
    # tested upper limit (for integers): 10000000000000000000000000000000000000000, where integer is > this number, and ends not with 0
    # tested lower limit: 1.0e-130
    # Therefore the numbers are adjusted as below, if they would otherwise trigger a DPI-1044 error

    def validate_oracle_statistics(self, statistics):
        for i, stat in enumerate(statistics):
            if stat:
                if abs(stat) > 10000000000000000000000000000000000000000:
                    statistics[i] = 10000000000000000000000000000000000000000
                elif abs(stat) < 1.0e-130:
                    statistics[i] = 1.0e-130
        return statistics

    def validate_postgres_statistics(self, statistics):
        for i, stat in enumerate(statistics):
            if stat:
                # Note: haven't tested the overflow limit
                if abs(stat) > 10000000000000000000000000000000000000000:
                    statistics[i] = 10000000000000000000000000000000000000000
                # According to an error message, the underflow limit is 6 decimal places for reals
                elif abs(stat) < 0.000001:
                    if stat > 0:
                        statistics[i] = 0.000001
                    elif stat == 0:
                        statistics[i] = 0
                    else:
                        statistics[i] = -0.000001
        return statistics


    def add_rule_environment_statistics_toList(
            self, rule_environment_id, property_name_id, statistics, cursor=None):

        count, avg, std, kurtosis, skewness, min, q1, median, q3, max, paired_t, p_value = statistics
        self.re_stats_inserts.append((rule_environment_id, property_name_id, count, avg, std, kurtosis,
                skewness, min, q1, median, q3, max, paired_t, p_value))

    def add_rule_environment_statistics_toDB(
            self, rule_environment_id, property_name_id, statistics, cursor=None):
        count, avg, std, kurtosis, skewness, min, q1, median, q3, max, paired_t, p_value = statistics
        self.mmpa_db.execute(
                "INSERT INTO rule_environment_statistics "
                "  (rule_environment_id, property_name_id, count, avg, std, kurtosis, "
                "       skewness, min, q1, median, q3, max, paired_t, p_value) "
                "  VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                (rule_environment_id, property_name_id, count, avg, std, kurtosis,
                skewness, min, q1, median, q3, max, paired_t, p_value),
                cursor=cursor)

    def update_rule_environment_statistics_toList(
            self, rule_environment_statistics_id, statistics, cursor=None):
        
        count, avg, std, kurtosis, skewness, min, q1, median, q3, max, paired_t, p_value = statistics
        self.re_stats_updates.append((count, avg, std, kurtosis, skewness,
             min, q1, median, q3, max, paired_t, p_value,
             rule_environment_statistics_id))

    def update_rule_environment_statistics_toDB(
            self, rule_environment_statistics_id, statistics, cursor=None):
        count, avg, std, kurtosis, skewness, min, q1, median, q3, max, paired_t, p_value = statistics
        self.mmpa_db.execute(
                "UPDATE rule_environment_statistics "
                " SET count=?, avg=?, std=?, kurtosis=?, skewness=?, "
                "     min=?, q1=?, median=?, q3=?, max=?, paired_t=?, p_value=? "
                " WHERE id = ?",
                (count, avg, std, kurtosis, skewness,
                min, q1, median, q3, max, paired_t, p_value,
                rule_environment_statistics_id),
                cursor=cursor)

    def batch_update_re_stats_to_postgres(self):

        c = self.mmpa_db.get_cursor()
        update_statement_begin = b"""UPDATE rule_environment_statistics
SET count = update_batch.count, avg = update_batch.avg, std = update_batch.std, kurtosis = update_batch.kurtosis, skewness = update_batch.skewness, min = update_batch.min,
    q1 = update_batch.q1, median = update_batch.median, q3 = update_batch.q3, max = update_batch.max, paired_t = update_batch.paired_t, p_value = update_batch.p_value
FROM (VALUES """
        update_values = b','.join(c.mogrify("(%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)", x) for x in self.re_stats_updates)
        update_values = update_values.replace(b"NULL", b"NULL::REAL")

        update_statement_end = b""") AS update_batch (count, avg, std, kurtosis, skewness, min, q1, median, q3, max, paired_t, p_value, id)
WHERE rule_environment_statistics.id = update_batch.id"""

        update_statement = update_statement_begin + update_values + update_statement_end
        c.execute(update_statement)
        # Clear the buffer, conserve memory
        self.re_stats_updates = []

    def batch_insert_re_stats_to_postgres(self):
        insert_statement = b"""INSERT INTO rule_environment_statistics 
(rule_environment_id, property_name_id, count, avg, std, kurtosis, 
skewness, min, q1, median, q3, max, paired_t, p_value) 
VALUES """
        c = self.mmpa_db.get_cursor()
        insert_values = b','.join(c.mogrify("(%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)", x) for x in self.re_stats_inserts)
        c.execute(insert_statement + insert_values)
        # Clear the buffer, conserve memory
        self.re_stats_inserts = []

    def batch_update_re_stats_to_oracle(self):
        update_statement = """
        UPDATE rule_environment_statistics 
        SET count=?, avg=?, std=?, kurtosis=?, skewness=?,
            min=?, q1=?, median=?, q3=?, max=?, paired_t=?, p_value=?
        WHERE id = ?
        """.replace('?', ':a')
        c = self.mmpa_db.get_cursor()
        # Pass up to 10 million rows at a time, into an executemany statement.
        # Adding this due to previous error when writing 98,208,426 pair rows:
        # cx_Oracle.DatabaseError: DPI-1015: array size of 98208426 is too large
        c.executemany(update_statement, self.re_stats_updates)
        # Clear the buffer, conserve memory
        self.re_stats_updates = []

    def batch_insert_re_stats_to_oracle(self):
        insert_statement = """
        INSERT INTO rule_environment_statistics 
                    (rule_environment_id, property_name_id, count, avg, std, kurtosis, 
                    skewness, min, q1, median, q3, max, paired_t, p_value) 
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """.replace('?', ':a')
        c = self.mmpa_db.get_cursor()
        c.executemany(insert_statement, self.re_stats_inserts)
        # Clear the buffer, conserve memory
        self.re_stats_inserts = []

    def delete_property_name_id(self, property_name_id, cursor=None):
        c = self.mmpa_db.get_cursor(cursor)
        self.mmpa_db.execute("DELETE FROM property_name "
                             " WHERE id = ?",
                             (property_name_id,), cursor=c)
        self.mmpa_db.execute("DELETE FROM compound_property "
                             " WHERE property_name_id = ?",
                             (property_name_id,), cursor=c)
        self.mmpa_db.execute("DELETE FROM rule_environment_statistics "
                             " WHERE property_name_id = ?",
                             (property_name_id,), cursor=c)
        
    def delete_compound_properties(self, compound_properties_to_delete, cursor=None):
        to_delete = sorted(compound_properties_to_delete)
        cursor = self.mmpa_db.get_cursor(cursor)
        cursor.executemany(self.mmpa_db.SQL("DELETE FROM compound_property "
                           " WHERE compound_id = ? AND property_name_id = ?"),
                           to_delete)
        
    def delete_rule_environment_statistics(self, ids_to_delete, cursor=None):
        ids_to_delete = [(id,) for id in sorted(ids_to_delete)]
        cursor = self.mmpa_db.get_cursor(cursor)

        # With cx_Oracle 8.1.0, talking to Oracle 11.2 DB,
        # an executemany on inserting 100 million pairs failed.
        # Will set to 10 million here
        #from .peewee import OracleDatabase
        #if isinstance(self.mmpa_db.db, OracleDatabase):
        if type(self.mmpa_db.db) in [peewee.OracleDatabase, peewee.CustomPostgresqlDatabase]:
            batch_size = 10000000
        else:
            batch_size = 10000

        while ids_to_delete:
            # Process 10,000 at a time because I don't know how many
            # I can really pass to it. (XXX Is this really a problem?)
            cursor.executemany(self.mmpa_db.SQL("DELETE FROM rule_environment_statistics "
                               " WHERE id = ?"),
                               ids_to_delete[:batch_size])
            del ids_to_delete[:batch_size]

    def iter_pairs(self, cursor=None):
        c = self.mmpa_db.execute("""
SELECT id, rule_environment_id, compound1_id, compound2_id, constant_id
  FROM pair
 ORDER BY rule_environment_id""",
                                cursor=cursor)
        for id, rule_environment_id, compound1_id, compound2_id, constant_id in c:
            yield Pair(id, rule_environment_id, compound1_id, compound2_id, constant_id)

    # tupleized_string contains compound_ids where a property value was inserted or updated
    def iter_pairs_with_property(self, tupleized_string):
        statement = """
        SELECT pair.id, pair.rule_environment_id, pair.compound1_id, pair.compound2_id, pair.constant_id
        FROM pair 
        WHERE (pair.compound1_id, 0) IN {} 
        AND (pair.compound2_id, 0) IN {} 
        ORDER BY pair.rule_environment_id
        """.format(tupleized_string, tupleized_string)
        pair_results = self.mmpa_db.execute(statement)
        for id, rule_environment_id, compound1_id, compound2_id, constant_id in pair_results:
            yield Pair(id, rule_environment_id, compound1_id, compound2_id, constant_id)

    # tupleized_string contains compound_ids where a property value was inserted or updated
    def get_pairs_with_property(self, tupleized_string):
        statement = """
        SELECT pair.id, pair.rule_environment_id, pair.compound1_id, pair.compound2_id, pair.constant_id
        FROM pair 
        WHERE (pair.compound1_id, 0) IN {} 
        AND (pair.compound2_id, 0) IN {} 
        ORDER BY pair.rule_environment_id
        """.format(tupleized_string, tupleized_string)
        pair_results = self.mmpa_db.execute(statement)
        pair_list = []
        i = 1
        for id, rule_environment_id, compound1_id, compound2_id, constant_id in pair_results:
            pair_list.append(Pair(id, rule_environment_id, compound1_id, compound2_id, constant_id))
    
    def get_property_rule_pairs(self, rule, property_name_id, cursor=None):
        sql = """
SELECT pair.id,
       lhs_compound.clean_smiles, rhs_compound.clean_smiles,
       lhs_compound.public_id, rhs_compound.public_id,
       lhs_property.value, rhs_property.value,
       rhs_property.value-lhs_property.value
  FROM pair,
       compound as lhs_compound, compound as rhs_compound,
       compound_property as lhs_property, compound_property as rhs_property
 WHERE pair.rule_environment_id = ?
   AND pair.compound1_id = lhs_compound.id
   AND pair.compound2_id = rhs_compound.id
   AND lhs_property.compound_id = lhs_compound.id
   AND lhs_property.property_name_id = ?
   AND rhs_property.compound_id = rhs_compound.id
   AND rhs_property.property_name_id = ?
        """
        c = self.mmpa_db.execute(
            sql, (rule.rule_environment_id, property_name_id, property_name_id),
            cursor=cursor)

        for (pair_id, lhs_smiles, rhs_smiles, lhs_public_id, rhs_public_id,
             lhs_value, rhs_value, delta) in c:
            if rule.is_reversed:
                yield PropertyRulePair(
                    rule, pair_id, rhs_smiles, lhs_smiles, rhs_public_id, lhs_public_id,
                    rhs_value, lhs_value, -delta)
            else:
                yield PropertyRulePair(
                    rule, pair_id, lhs_smiles, rhs_smiles, lhs_public_id, rhs_public_id,
                    lhs_value, rhs_value, delta)

    def iter_compounds(self, cursor=None):
        c = self.mmpa_db.execute(
            "SELECT id, public_id, input_smiles, clean_smiles, clean_num_heavies "
            "  FROM compound "
            " ORDER BY public_id")
        for row in c:
            yield CompoundRow(*row)


class PropertyRulePair(object):
    def __init__(self, rule, pair_id, lhs_smiles, rhs_smiles, lhs_public_id, rhs_public_id,
                 lhs_value, rhs_value, delta):
        self.__dict__.update(rule.to_dict())  # Incorporate all of the rule fields
        self.rule = rule
        
        self.pair_id = pair_id
        self.lhs_smiles = lhs_smiles
        self.rhs_smiles = rhs_smiles
        self.lhs_public_id = lhs_public_id
        self.rhs_public_id = rhs_public_id
        self.lhs_value = lhs_value
        self.rhs_value = rhs_value
        self.delta = delta
    

class PropertyRule(object):
    # Note: analysis_algorithms.check_eval_names() uses the following to
    # check if the passed-in eval'ed code contains unexpected names.
    __slots__ = ("rule_id", "is_reversed", "from_smiles", "from_num_heavies", "to_smiles", "to_num_heavies", "smirks",
                 "rule_environment_id", "radius", "fingerprint_id", "fingerprint",
                 "rule_environment_statistics_id", "count", "avg", "std", "kurtosis", "skewness", "min", "q1", "median", "q3", "max",
                 "paired_t", "p_value",
                 "is_bidirectional")
    
    def __init__(self, rule_id, is_reversed, from_smiles, from_num_heavies, to_smiles, to_num_heavies,
                 rule_environment_id, radius, fingerprint_id, fingerprint,
                 rule_environment_statistics_id, count, avg, std, kurtosis, skewness, min, q1, median, q3, max, paired_t, p_value):
        self.rule_id = rule_id
        self.is_reversed = is_reversed
        
        self.from_smiles = from_smiles
        self.from_num_heavies = from_num_heavies
        self.to_smiles = to_smiles
        self.to_num_heavies = to_num_heavies
        self.smirks = from_smiles + ">>" + to_smiles

        self.rule_environment_id = rule_environment_id
        self.radius = radius
        self.fingerprint_id = fingerprint_id
        self.fingerprint = fingerprint
        
        self.rule_environment_statistics_id = rule_environment_statistics_id
        self.count = count
        self.avg = avg
        self.std = std
        self.kurtosis = kurtosis
        self.skewness = skewness
        self.min = min
        self.q1 = q1
        self.median = median
        self.q3 = q3
        self.max = max
        self.paired_t = paired_t
        self.p_value = p_value
        self.is_bidirectional = None

    def to_dict(self):
        return {
            "rule_id": self.rule_id,
            "is_reversed": self.is_reversed,

            "from_smiles": self.from_smiles,
            "from_num_heavies": self.from_num_heavies,
            "to_smiles": self.to_smiles,
            "to_num_heavies": self.to_num_heavies,
            "smirks": self.smirks,

            "rule_environment_id": self.rule_environment_id,
            "radius": self.radius,
            "fingerprint": self.fingerprint,

            "rule_environment_statistics_id": self.rule_environment_statistics_id, 
            "count": self.count,
            "avg": self.avg,
            "std": self.std,
            "kurtosis": self.kurtosis,
            "skewness": self.skewness,
            "min": self.min,
            "q1": self.q1,
            "median": self.median,
            "q3": self.q3,
            "max": self.max,
            "paired_t": self.paired_t,
            "p_value": self.p_value,

            "is_bidirectional": self.is_bidirectional,
            }


class CompoundRow(object):
    __slots__ = ("id", "public_id", "input_smiles", "clean_smiles", "clean_num_heavies")

    def __init__(self, id, public_id, input_smiles, clean_smiles, clean_num_heavies):
        self.id = id
        self.public_id = public_id
        self.input_smiles = input_smiles
        self.clean_smiles = clean_smiles
        self.clean_num_heavies = clean_num_heavies
