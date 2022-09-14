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

from __future__ import absolute_import, print_function

# Some low-level/backend writers for MMPWriter.
# The options are:
#   - save as a simple line-oriented format
#   - save the data to a SQLite database using Python's sqlite3 library
#   - save the data to a SQLite database using the third-party apsw library.

import os
import sqlite3
import datetime
import json
import collections
import itertools
import tempfile
import atexit
import shutil
from . import dbutils
from . import matcher_setup

try:
    import cx_Oracle
except ImportError:
    cx_Oracle = None

# To install apsw, do:
# pip install --user https://github.com/rogerbinns/apsw/releases/download/3.16.2-r1/apsw-3.16.2-r1.zip \
# --global-option=fetch --global-option=--version --global-option=3.16.2 --global-option=--all \
# --global-option=build --global-option=--enable-all-extensions

try:
    import apsw
except ImportError:
    apsw = None
    ## The speedup is only about 10%. Not enough to make a firm suggestion.
    ## import sys
    ## sys.stderr.write("You should install apsw.\n")
    import sqlite3

from . import schema
from .fragment_algorithm import get_num_heavies_from_smiles

nan = float("nan")


class TableIndexWriter(object):
    def __init__(self, outfile):
        self.outfile = outfile
        self._W = outfile.write

    def close(self):
        self.outfile.close()
        
    def rollback(self):
        self._W("ROLLBACK")
        self.close()
        
    def commit(self):
        self._W("COMMIT")
        self.close()
        
    def start(self, fragment_options, index_options):
        self._W("VERSION\tmmpa/3\n")
        self._W("FRAGMENT_OPTIONS\t%s\n" % (json.dumps(list(fragment_options.to_dict().items())),))
        self._W("INDEX_OPTIONS\t%s\n" % (json.dumps(list(index_options.to_dict().items())),))

    def add_property_name(self, property_name_idx, property_name):
        self._W("PROPNAME\t%d\t%s\n" % (property_name_idx, property_name))
        
    def add_rule_smiles(self, smiles_idx, smiles):
        self._W("RULE_SMILES\t%d\t%s\n" % (smiles_idx, smiles))

    def add_rule(self, rule_idx, from_smiles_idx, to_smiles_idx):
        self._W("RULE\t%d\t%d\t%d\n" % (rule_idx, from_smiles_idx, to_smiles_idx))

    def add_environment_fingerprint(self, fp_idx, environment_fingerprint):
        self._W("FINGERPRINT\t%d\t%s\n" % (fp_idx, environment_fingerprint))

    # Added to account for parents in indexing (only has an effect on SQLite Tables)
    def add_environment_fingerprint_parent(self, fp_idx, environment_fingerprint, parent_idx):
        self._W("FINGERPRINT\t%d\t%s\n" % (fp_idx, environment_fingerprint, parent_idx))

    def add_rule_environment(self, rule_env_idx, rule_idx, env_fp_idx, radius):
        self._W("RULEENV\t%d\t%d\t%d\t%d\n" % (rule_env_idx, rule_idx, env_fp_idx, radius))

    def add_compound(self, compound_idx, compound_id, input_smiles,
                     normalized_smiles, num_normalized_heavies):
        self._W("COMPOUND\t%d\t%s\t%s\t%s\t%d\n" % (
            compound_idx, compound_id, input_smiles,
            normalized_smiles, num_normalized_heavies))
        
    def add_constant_smiles(self, smiles_idx, constant_smiles):
        self._W("CONSTANT_SMILES\t%d\t%s\n" % (smiles_idx, constant_smiles))

    def add_rule_environment_pair(self, pair_idx, env_idx, compound1_idx, compound2_idx, constant_idx):
        self._W("PAIR%d\t\t%d\t%d\t%d\t%d\n" % (pair_idx, env_idx, compound1_idx, compound2_idx, constant_idx))

    def add_compound_property(self, compound_idx, property_name_idx, value):
        self._W("PROP\t%d\t%d\t%s\n" % (compound_idx, property_name_idx, value))

    def add_rule_environment_statistics(self, rule_env_idx, property_name_idx, values):
        self._W("RULEENV_STATS\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %
                    ((rule_env_idx, property_name_idx) + tuple(values)))

    def end(self, reporter):
        pass


def open_table_index_writer(outfile):
    return TableIndexWriter(outfile)


class BaseSqliteIndexWriter(object):
    def __init__(self, db, conn, title):
        self.db = db
        self.conn = conn
        self.title = title

    def start(self, fragment_options, index_options):
        creation_date = datetime.datetime.now().isoformat(sep=" ")
        fragment_options_str = json.dumps(fragment_options.to_dict())
        index_options_str = json.dumps(index_options.to_dict())
    
        self.conn.execute("INSERT INTO dataset (id, mmpdb_version, title, creation_date, "
                          "    fragment_options, index_options, is_symmetric) "
                          "    VALUES (1, 2, ?, ?, ?, ?, ?)",
                        (self.title, creation_date, fragment_options_str,
                         index_options_str, index_options.symmetric))

    def add_property_name(self, property_name_idx, property_name):
        self.conn.execute("INSERT INTO property_name (id, name) VALUES (?, ?)",
                          (property_name_idx, property_name))
        
    def add_rule_smiles(self, smiles_idx, smiles):
        self.conn.execute("INSERT INTO rule_smiles (id, smiles, num_heavies) VALUES (?, ?, ?)",
                          (smiles_idx, smiles, get_num_heavies_from_smiles(smiles)))

    def add_rule(self, rule_idx, from_smiles_idx, to_smiles_idx):
        self.conn.execute("INSERT INTO rule (id, from_smiles_id, to_smiles_id) "
                          "  VALUES (?, ?, ?)",
                          (rule_idx, from_smiles_idx, to_smiles_idx))

    def add_environment_fingerprint(self, fp_idx, environment_fingerprint):
        self.conn.execute("INSERT INTO environment_fingerprint (id, fingerprint) "
                          " VALUES (?, ?)",
                          (fp_idx, environment_fingerprint))

    # Added to include a parent idx in the EnvironmentFingerprint table
    def add_environment_fingerprint_parent(self, fp_idx, environment_fingerprint, parent_idx):
        self.conn.execute("INSERT INTO environment_fingerprint (id, fingerprint, parent_id) "
                          " VALUES (?, ?, ?)",
                          (fp_idx, environment_fingerprint, parent_idx))

    def add_rule_environment(self, rule_env_idx, rule_idx, env_fp_idx, radius):
        self.conn.execute("INSERT INTO rule_environment (id, rule_id, environment_fingerprint_id,  radius) "
                          "  VALUES (?, ?, ?, ?)",
                          (rule_env_idx, rule_idx, env_fp_idx, radius))

    def add_compound(self, compound_idx, compound_id, input_smiles,
                     normalized_smiles, num_normalized_heavies):
        self.conn.execute("INSERT INTO compound (id, public_id, input_smiles, clean_smiles, clean_num_heavies) "
                          "   VALUES (?, ?, ?, ?, ?)",
                          (compound_idx, compound_id, input_smiles, normalized_smiles, num_normalized_heavies))
        
    def add_constant_smiles(self, smiles_idx, constant_smiles):
        self.conn.execute("INSERT INTO constant_smiles (id, smiles) VALUES (?, ?)",
                          (smiles_idx, constant_smiles))

    def add_rule_environment_pair(self, pair_idx, env_idx, compound1_idx, compound2_idx, constant_idx):
        self.conn.execute("INSERT INTO pair (id, rule_environment_id, compound1_id, compound2_id, constant_id) "
                          "  VALUES (?, ?, ?, ?, ?)",
                          (pair_idx, env_idx, compound1_idx, compound2_idx, constant_idx))

    def add_compound_property(self, compound_idx, property_name_idx, value):
        self.conn.execute("INSERT INTO compound_property (compound_id, property_name_id, value) VALUES (?, ?, ?)",
                          (compound_idx, property_name_idx, value))

    def add_rule_environment_statistics(self, rule_env_idx, property_name_idx, values):
        count, avg, std, kurtosis, skewness, min, q1, median, q3, max, paired_t, p_value = values
        assert rule_env_idx is not None
        assert property_name_idx is not None
        assert count is not None
        assert avg is not None
        assert min is not None
        assert q1 is not None
        assert median is not None
        assert q3 is not None
        assert max is not None
        # XXX check for too-large/infinite values?

        self.conn.execute("INSERT INTO rule_environment_statistics "
                          "  (rule_environment_id, property_name_id, count, avg, std, kurtosis, "
                          "       skewness, min, q1, median, q3, max, paired_t, p_value) "
                          "  VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                          (rule_env_idx, property_name_idx, count, avg, std,
                           kurtosis, skewness, min, q1, median, q3, max, paired_t, p_value))
        
    def end(self, reporter):
        reporter.update("Building index ...")
        schema.create_index(self.conn)
        
        # Improve SQLite query planning
        reporter.update("Analyzing database ...")
        self.conn.execute("ANALYZE")
        
        reporter.update("Computing sizes ...")
        num_compounds = schema._get_one(self.conn.execute("SELECT count(*) from compound"))
        num_rules = schema._get_one(self.conn.execute("SELECT count(*) from rule"))
        num_pairs = schema._get_one(self.conn.execute("SELECT count(*) from pair"))
        num_envs = schema._get_one(self.conn.execute("SELECT count(*) from rule_environment"))
        num_stats = schema._get_one(self.conn.execute("SELECT count(*) from rule_environment_statistics"))
        self.conn.execute("UPDATE dataset set num_compounds=?, num_rules=?, num_pairs=?, "
                          "num_rule_environments=?, num_rule_environment_stats=? WHERE id = 1",
                          (num_compounds, num_rules, num_pairs, num_envs, num_stats))
        
        reporter.update("")


class SQLiteIndexWriter(BaseSqliteIndexWriter):
    def close(self):
        self.conn.close()
        self.db.commit()
        self.db.close()

    def commit(self):
        self.conn.close()
        self.db.commit()
        self.db.close()

    def rollback(self):
        self.conn.close()
        self.db.close()


class APSWIndexWriter(BaseSqliteIndexWriter):
    def start(self, fragment_options, index_options):
        self.conn.execute("BEGIN TRANSACTION")
        super(APSWIndexWriter, self).start(fragment_options, index_options)
    
    def close(self):
        self.conn.close()
        self.db.execute("COMMIT")
        self.db.close()

    def commit(self):
        self.conn.execute("COMMIT")
        self.conn.close()
        self.db.close()

    def rollback(self):
        #self.conn.execute("ROLLBACK")
        self.conn.close()
        self.db.close()

import time

# Copied and pasted BaseSqliteIndexWriter + SQLiteIndexWriter, but changed syntax in executed statements, to make syntax compatible with cx_Oracle    
class OracleIndexWriter(object):
    def __init__(self, db, conn, title):
        self.db = db
        self.conn = conn
        self.title = title
        # Adding lists for batch statements with Oracle DB
        # Need to do batch statements with large datasets, otherwise latency to DB will kill indexing speed
        # Currently intending to do single batch statement for writing below lists to DB,
        # at end of index_algorithm.write_matched_molecule_pairs function
        # If memory runs out while indexing large datasets, may need to break the process up,
        # executing a batch statement before memory runs out, deleting the below lists, then refilling them
        self.rule_smiles_listForDB = []
        self.rule_listForDB = []
        self.environment_fingerprint_listForDB = []
        self.environment_fingerprint_parent_listForDB = []
        self.rule_environment_listForDB = []
        self.compound_listForDB = []
        self.constant_smiles_listForDB = []
        self.rule_environment_pair_listForDB = []
        self.compound_property_listForDB = []

    def start(self, fragment_options, index_options):
        # Oracle TIMESTAMP variable requires date format such as "15-Mar-2021 10:29:00.123456"
        creation_date = datetime.datetime.now().strftime("%d-%b-%Y %I:%M:%S.%f %p")
        fragment_options_str = json.dumps(fragment_options.to_dict())
        index_options_str = json.dumps(index_options.to_dict())
    
        self.conn.execute("INSERT INTO dataset (id, mmpdb_version, title, creation_date, "
                          "    fragment_options, index_options, is_symmetric) "
                          "    VALUES (1, 2, :a, :b, :c, :d, :e)",
                        [self.title, creation_date, fragment_options_str,
                         index_options_str, index_options.symmetric])

    def add_property_name(self, property_name_idx, property_name):
        self.conn.execute("INSERT INTO property_name (id, name) VALUES (:a, :b)",
                          [property_name_idx, property_name])
        
    def add_rule_smiles(self, smiles_idx, smiles):
        self.rule_smiles_listForDB.append((smiles_idx, smiles, get_num_heavies_from_smiles(smiles)))
        
        #self.conn.execute("INSERT INTO rule_smiles (id, smiles, num_heavies) VALUES (:a, :b, :c)",
                          #[smiles_idx, smiles, get_num_heavies_from_smiles(smiles)])

    def add_rule(self, rule_idx, from_smiles_idx, to_smiles_idx):
        self.rule_listForDB.append((rule_idx, from_smiles_idx, to_smiles_idx))
        
        #self.conn.execute("INSERT INTO rule (id, from_smiles_id, to_smiles_id) "
                          #"  VALUES (:a, :b, :c)",
                          #[rule_idx, from_smiles_idx, to_smiles_idx])

    def add_environment_fingerprint(self, fp_idx, environment_fingerprint):
        self.environment_fingerprint_listForDB.append((fp_idx, environment_fingerprint))
        
        #self.conn.execute("INSERT INTO environment_fingerprint (id, fingerprint) "
                          #" VALUES (:a, :b)",
                          #[fp_idx, environment_fingerprint])

    # Added to include a parent idx in the EnvironmentFingerprint table
    def add_environment_fingerprint_parent(self, fp_idx, environment_fingerprint, parent_idx):
        self.environment_fingerprint_parent_listForDB.append((fp_idx, environment_fingerprint, parent_idx))
        
        #self.conn.execute("INSERT INTO environment_fingerprint (id, fingerprint, parent_id) "
                          #" VALUES (:a, :b, :c)",
                          #[fp_idx, environment_fingerprint, parent_idx])

    def add_rule_environment(self, rule_env_idx, rule_idx, env_fp_idx, radius):
        self.rule_environment_listForDB.append((rule_env_idx, rule_idx, env_fp_idx, radius))
        
        #self.conn.execute("INSERT INTO rule_environment (id, rule_id, environment_fingerprint_id,  radius) "
                          #"  VALUES (:a, :b, :c, :d)",
                          #[rule_env_idx, rule_idx, env_fp_idx, radius])

    def add_compound(self, compound_idx, compound_id, input_smiles,
                     normalized_smiles, num_normalized_heavies):
        self.compound_listForDB.append((compound_idx, compound_id, input_smiles, normalized_smiles, num_normalized_heavies))
        
        #self.conn.execute("INSERT INTO compound (id, public_id, input_smiles, clean_smiles, clean_num_heavies) "
                          #"   VALUES (:a, :b, :c, :d, :e)",
                          #[compound_idx, compound_id, input_smiles, normalized_smiles, num_normalized_heavies])
        
    def add_constant_smiles(self, smiles_idx, constant_smiles):
        self.constant_smiles_listForDB.append((smiles_idx, constant_smiles))
        
        #self.conn.execute("INSERT INTO constant_smiles (id, smiles) VALUES (:a, :b)",
                          #[smiles_idx, constant_smiles])

    def add_rule_environment_pair(self, pair_idx, env_idx, compound1_idx, compound2_idx, constant_idx):
        self.rule_environment_pair_listForDB.append((pair_idx, env_idx, compound1_idx, compound2_idx, constant_idx))
        
        #self.conn.execute("INSERT INTO pair (id, rule_environment_id, compound1_id, compound2_id, constant_id) "
                          #"  VALUES (:a, :b, :c, :d, :e)",
                          #[pair_idx, env_idx, compound1_idx, compound2_idx, constant_idx])

    def add_compound_property(self, compound_idx, property_name_idx, value):
        self.compound_property_listForDB.append((compound_idx, property_name_idx, value))
        
        #self.conn.execute("INSERT INTO compound_property (compound_id, property_name_id, value) VALUES (:a, :b, :c)",
                          #[compound_idx, property_name_idx, value])

    def batch_write_to_oracle(self):
        
        def to_batches(original_list, batch_size):
            return(original_list[a:a+batch_size] for a in range(0,len(original_list),batch_size))

        def write_batches(statement, original_list, batch_size):
            for batch in to_batches(original_list, batch_size):
                self.conn.executemany(statement, batch)

        start_batch_write_time = time.time()
        print("Starting oracle writing at " + str(start_batch_write_time))

        # Pass up to 10 million rows at a time, into an executemany statement.
        # Adding this due to previous error when writing 98,208,426 pair rows:
        # cx_Oracle.DatabaseError: DPI-1015: array size of 98208426 is too large
        batch_size = 6100000

        statement = "INSERT INTO rule_smiles (id, smiles, num_heavies) VALUES (:a, :b, :c)"
        write_batches(statement, self.rule_smiles_listForDB, batch_size)
        t_1 = time.time()
        print(str(len(self.rule_smiles_listForDB)) + " rule_smiles took " + str(t_1 - start_batch_write_time) + " seconds to write")

        statement = """
        INSERT INTO rule (id, from_smiles_id, to_smiles_id)  
        VALUES (:a, :b, :c)
        """
        write_batches(statement, self.rule_listForDB, batch_size)
        t_2 = time.time()
        print(str(len(self.rule_listForDB)) + " rules took " + str(t_2 - t_1) + " seconds to write")
        
        num_fps = 0

        if self.environment_fingerprint_listForDB:
            statement = """
            INSERT INTO environment_fingerprint (id, fingerprint) 
            VALUES (:a, :b)
            """
            write_batches(statement, self.environment_fingerprint_listForDB, batch_size)
            num_fps = len(self.environment_fingerprint_listForDB)
            
        if self.environment_fingerprint_parent_listForDB:
            statement = """
            INSERT INTO environment_fingerprint (id, fingerprint, parent_id) 
            VALUES (:a, :b, :c)
            """
            write_batches(statement, self.environment_fingerprint_parent_listForDB, batch_size)
            num_fps = len(self.environment_fingerprint_parent_listForDB)
        
        t_3 = time.time()
        print(str(num_fps) + " efps took " + str(t_3 - t_2) + " seconds to write")
        
        statement = """
        INSERT INTO rule_environment (id, rule_id, environment_fingerprint_id,  radius) 
        VALUES (:a, :b, :c, :d)
        """
        write_batches(statement, self.rule_environment_listForDB, batch_size)
        t_4 = time.time()
        print(str(len(self.rule_environment_listForDB)) + " rule_envs took " + str(t_4 - t_3) + " seconds to write")
        
        statement = """
        INSERT INTO compound (id, public_id, input_smiles, clean_smiles, clean_num_heavies) 
        VALUES (:a, :b, :c, :d, :e)
        """
        write_batches(statement, self.compound_listForDB, batch_size)
        t_5 = time.time()
        print(str(len(self.compound_listForDB)) + " compounds took " + str(t_5 - t_4) + " seconds to write")
        
        statement = """
        INSERT INTO constant_smiles (id, smiles) VALUES (:a, :b) 
        """
        write_batches(statement, self.constant_smiles_listForDB, batch_size)
        t_6 = time.time()
        print(str(len(self.constant_smiles_listForDB)) + " constant_smiles took " + str(t_6 - t_5) + " seconds to write")
        
        statement = """
        INSERT INTO pair (id, rule_environment_id, compound1_id, compound2_id, constant_id) 
        VALUES (:a, :b, :c, :d, :e)
        """
        write_batches(statement, self.rule_environment_pair_listForDB, batch_size)
        t_7 = time.time()
        print(str(len(self.rule_environment_pair_listForDB)) + " pairs took " + str(t_7 - t_6) + " seconds to write")
        
        if self.compound_property_listForDB:
            statement = """
            INSERT INTO compound_property (compound_id, property_name_id, value) 
            VALUES (:a, :b, :c)
            """
            write_batches(statement, self.compound_property_listForDB, batch_size)

        self.db.commit()

        # Empty out lists after batch writing, to conserve memory
        self.rule_smiles_listForDB = []
        self.rule_listForDB = []
        self.environment_fingerprint_listForDB = []
        self.environment_fingerprint_parent_listForDB = []
        self.rule_environment_listForDB = []
        self.compound_listForDB = []
        self.constant_smiles_listForDB = []
        self.rule_environment_pair_listForDB = []
        self.compound_property_listForDB = []

    def add_rule_environment_statistics(self, rule_env_idx, property_name_idx, values):
        count, avg, std, kurtosis, skewness, min, q1, median, q3, max, paired_t, p_value = values
        assert rule_env_idx is not None
        assert property_name_idx is not None
        assert count is not None
        assert avg is not None
        assert min is not None
        assert q1 is not None
        assert median is not None
        assert q3 is not None
        assert max is not None
        # XXX check for too-large/infinite values?

        self.conn.execute("INSERT INTO rule_environment_statistics "
                          "  (rule_environment_id, property_name_id, count, avg, std, kurtosis, "
                          "       skewness, min, q1, median, q3, max, paired_t, p_value) "
                          "  VALUES (:a, :b, :c, :d, :e, :f, :g, :h, :i, :j, :k, :l, :m, :n)",
                          [rule_env_idx, property_name_idx, count, avg, std,
                           kurtosis, skewness, min, q1, median, q3, max, paired_t, p_value])
        
    def end(self, reporter):
        reporter.update("Building index ...")
        # Added a new oracle function pathway, due to need to remove ; in create_index.sql for oracle statements
        schema.create_oracle_index(self.conn)
        
        # Improve SQLite query planning
        reporter.update("Analyzing database ...")

        # Oracle does not tolerate below analyze syntax
        # self.conn.execute("ANALYZE")
        
        reporter.update("Computing sizes ...")
        num_compounds = schema._get_one(self.conn.execute("SELECT count(*) from compound"))
        num_rules = schema._get_one(self.conn.execute("SELECT count(*) from rule"))
        num_pairs = schema._get_one(self.conn.execute("SELECT count(*) from pair"))
        num_envs = schema._get_one(self.conn.execute("SELECT count(*) from rule_environment"))
        num_stats = schema._get_one(self.conn.execute("SELECT count(*) from rule_environment_statistics"))
        self.conn.execute("UPDATE dataset set num_compounds=:a, num_rules=:b, num_pairs=:c, "
                          "num_rule_environments=:d, num_rule_environment_stats=:e WHERE id = 1",
                          [num_compounds, num_rules, num_pairs, num_envs, num_stats])
        
        reporter.update("")
    
    def close(self):
        self.conn.close()
        self.db.commit()
        self.db.close()

    def commit(self):
        self.conn.close()
        self.db.commit()
        self.db.close()

    def rollback(self):
        self.conn.close()
        self.db.close()

class PostgresIndexWriter(OracleIndexWriter):
    def start(self, fragment_options, index_options):
        # Example of valid Postgres TIMESTAMP is "2021-09-24 16:20:00"
        creation_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        fragment_options_str = json.dumps(fragment_options.to_dict())
        index_options_str = json.dumps(index_options.to_dict())
        # Need to convert boolean to integer, otherwise psycopg2 throws an error when trying to insert into an integer column is_symmetric
        if index_options.symmetric == False:
            is_symmetric_int = 0
        else:
            is_symmetric_int = 1
        self.conn.execute("INSERT INTO dataset (id, mmpdb_version, title, creation_date, "
                          "    fragment_options, index_options, is_symmetric) "
                          "    VALUES (1, 2, %s, %s, %s, %s, %s)",
                        [self.title, creation_date, fragment_options_str,
                         index_options_str, is_symmetric_int])

    def add_property_name(self, property_name_idx, property_name):
        self.conn.execute("INSERT INTO property_name (id, name) VALUES (%s, %s)",
                          [property_name_idx, property_name])

    def batch_write_to_postgres(self):
        
        def to_batches(original_list, batch_size):
            return(original_list[a:a+batch_size] for a in range(0,len(original_list),batch_size))

        def write_batches(statement_begin, statement_end, original_list, batch_size):
            for batch in to_batches(original_list, batch_size):
                mogrified_batch = b','.join(self.conn.mogrify(statement_end, x) for x in batch)
                #mogrified_batch = ','.join(self.conn.mogrify("(%s,%s,%s,%s,%s,%s,%s,%s,%s)", x) for x in batch)
                #self.conn.execute("INSERT INTO table VALUES " + mogrified_batch)
                self.conn.execute(statement_begin + mogrified_batch)

        start_batch_write_time = time.time()
        print("Starting postgres writing at " + str(start_batch_write_time))

        # Pass up to 10 million rows at a time, into an executemany statement.
        # Adding this due to previous error when writing 98,208,426 pair rows:
        # cx_Oracle.DatabaseError: DPI-1015: array size of 98208426 is too large
        batch_size = 610000

        statement_begin = b"INSERT INTO rule_smiles (id, smiles, num_heavies) VALUES "
        statement_end = "(%s,%s,%s)"
        write_batches(statement_begin, statement_end, self.rule_smiles_listForDB, batch_size)
        t_1 = time.time()
        print(str(len(self.rule_smiles_listForDB)) + " rule_smiles took " + str(t_1 - start_batch_write_time) + " seconds to write")

        statement_begin = b"INSERT INTO rule (id, from_smiles_id, to_smiles_id) VALUES "
        statement_end = "(%s,%s,%s)"
        write_batches(statement_begin, statement_end, self.rule_listForDB, batch_size)
        t_2 = time.time()
        print(str(len(self.rule_listForDB)) + " rules took " + str(t_2 - t_1) + " seconds to write")
        
        num_fps = 0

        if self.environment_fingerprint_listForDB:
            statement_begin = b"INSERT INTO environment_fingerprint (id, fingerprint) VALUES "
            statement_end = "(%s,%s)"
            write_batches(statement_begin, statement_end, self.environment_fingerprint_listForDB, batch_size)
            num_fps = len(self.environment_fingerprint_listForDB)
            
        if self.environment_fingerprint_parent_listForDB:
            statement_begin = b"INSERT INTO environment_fingerprint (id, fingerprint, parent_id) VALUES "
            statement_end = "(%s,%s,%s)"
            write_batches(statement_begin, statement_end, self.environment_fingerprint_parent_listForDB, batch_size)
            num_fps = len(self.environment_fingerprint_parent_listForDB)
        
        t_3 = time.time()
        print(str(num_fps) + " efps took " + str(t_3 - t_2) + " seconds to write")
        
        statement_begin = b"INSERT INTO rule_environment (id, rule_id, environment_fingerprint_id, radius) VALUES "
        statement_end = "(%s,%s,%s,%s)"
        write_batches(statement_begin, statement_end, self.rule_environment_listForDB, batch_size)
        t_4 = time.time()
        print(str(len(self.rule_environment_listForDB)) + " rule_envs took " + str(t_4 - t_3) + " seconds to write")
        
        statement_begin = b"INSERT INTO compound (id, public_id, input_smiles, clean_smiles, clean_num_heavies) VALUES "
        statement_end = "(%s,%s,%s,%s,%s)"
        write_batches(statement_begin, statement_end, self.compound_listForDB, batch_size)
        t_5 = time.time()
        print(str(len(self.compound_listForDB)) + " compounds took " + str(t_5 - t_4) + " seconds to write")
        
        statement_begin = b"INSERT INTO constant_smiles (id, smiles) VALUES "
        statement_end = "(%s,%s)"
        write_batches(statement_begin, statement_end, self.constant_smiles_listForDB, batch_size)
        t_6 = time.time()
        print(str(len(self.constant_smiles_listForDB)) + " constant_smiles took " + str(t_6 - t_5) + " seconds to write")
        
        statement_begin = b"INSERT INTO pair (id, rule_environment_id, compound1_id, compound2_id, constant_id) VALUES "
        statement_end = "(%s,%s,%s,%s,%s)"
        write_batches(statement_begin, statement_end, self.rule_environment_pair_listForDB, batch_size)
        t_7 = time.time()
        print(str(len(self.rule_environment_pair_listForDB)) + " pairs took " + str(t_7 - t_6) + " seconds to write")
        
        if self.compound_property_listForDB:
            statement_begin = b"INSERT INTO compound_property (compound_id, property_name_id, value) VALUES "
            statement_end = "(%s,%s,%s)"
            write_batches(statement_begin, statement_end, self.compound_property_listForDB, batch_size)

        self.db.commit()

        # Empty out lists after batch writing, to conserve memory
        self.rule_smiles_listForDB = []
        self.rule_listForDB = []
        self.environment_fingerprint_listForDB = []
        self.environment_fingerprint_parent_listForDB = []
        self.rule_environment_listForDB = []
        self.compound_listForDB = []
        self.constant_smiles_listForDB = []
        self.rule_environment_pair_listForDB = []
        self.compound_property_listForDB = []

    def add_rule_environment_statistics(self, rule_env_idx, property_name_idx, values):
        count, avg, std, kurtosis, skewness, min, q1, median, q3, max, paired_t, p_value = values
        assert rule_env_idx is not None
        assert property_name_idx is not None
        assert count is not None
        assert avg is not None
        assert min is not None
        assert q1 is not None
        assert median is not None
        assert q3 is not None
        assert max is not None
        # XXX check for too-large/infinite values?

        self.conn.execute("INSERT INTO rule_environment_statistics "
                          "  (rule_environment_id, property_name_id, count, avg, std, kurtosis, "
                          "       skewness, min, q1, median, q3, max, paired_t, p_value) "
                          "  VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)",
                          [rule_env_idx, property_name_idx, count, avg, std,
                           kurtosis, skewness, min, q1, median, q3, max, paired_t, p_value])
        
    def end(self, reporter):
        reporter.update("Building index ...")
        # Added a new oracle function pathway, due to need to remove ; in create_index.sql for oracle statements
        schema.create_oracle_index(self.conn)
        
        reporter.update("Computing sizes ...")
        # With psycopg2, passing the self.conn.execute(statement) as an argument causes it to evaluate to None
        # It works if we execute the statement first, then pass the cursor as an argument
        self.conn.execute("SELECT count(*) from compound")
        num_compounds = schema._get_one(self.conn)
        self.conn.execute("SELECT count(*) from rule")
        num_rules = schema._get_one(self.conn)
        self.conn.execute("SELECT count(*) from pair")
        num_pairs = schema._get_one(self.conn)
        self.conn.execute("SELECT count(*) from rule_environment")
        num_envs = schema._get_one(self.conn)
        self.conn.execute("SELECT count(*) from rule_environment_statistics")
        num_stats = schema._get_one(self.conn)
        self.conn.execute("UPDATE dataset set num_compounds=%s, num_rules=%s, num_pairs=%s, "
                          "num_rule_environments=%s, num_rule_environment_stats=%s WHERE id = 1",
                          [num_compounds, num_rules, num_pairs, num_envs, num_stats])

        reporter.update("Adding Matcher-specific data (e.g. cartridge molecules, constructs)")
        t1 = time.time()
        matcher_setup.extend_postgres_build(connection=self.db, cursor=self.conn)
        t2 = time.time()
        print("Matcher-specific postgres extension took (seconds): " + str(round(t2-t1, 3)))
        
        ## Skip analyze for now, let's analyze after all the props are added later
        # Improve query planning
        #reporter.update("Analyzing database ...")
        #self.conn.execute("ANALYZE")

        reporter.update("")

def open_sqlite_index_writer(filename, title):
    if filename != ":memory:":
        if os.path.exists(filename):
            os.unlink(filename)
    if apsw is None:
        db = sqlite3.connect(filename)
        klass = SQLiteIndexWriter
    else:
        db = apsw.Connection(filename)
        klass = APSWIndexWriter
    
    schema.create_schema_for_sqlite(db)
    conn = db.cursor()

    return klass(db, conn, title)

def open_oracle_index_writer(filename, title):
    db = dbutils.open_oracle_database(filename)
    schema.create_schema_for_oracle(db)
    conn = db.cursor()
   
    return OracleIndexWriter(db, conn, title)

def open_postgres_index_writer(filename, title):
    db = dbutils.open_postgres_database(filename)
    schema.create_schema_for_postgres(db)
    conn = db.cursor()
   
    return PostgresIndexWriter(db, conn, title)