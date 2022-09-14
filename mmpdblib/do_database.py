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

from __future__ import print_function, absolute_import

# list [database*]

import sys
import json
import itertools
import time
import datetime

from . import schema
from . import command_support
from . import dbutils
from . import index_algorithm
from . import fileio    

from . import peewee
# mmpdb list [--all] [--quiet] [--recount] [filename]*

def list_command(parser, args):
    reporter = command_support.get_reporter(args.quiet)
    databases = args.databases

    name_list = []
    num_compounds_list = []
    num_rules_list = []
    num_pairs_list = []
    num_envs_list = []
    num_stats_list = []
    titles = []
    property_names = []

    name_width = 5
    num_compounds_width = 6
    num_rules_width = 6
    num_pairs_width = 6
    num_envs_width = 6
    num_stats_width = 6
    title_width = 7

    all_fragment_options = []
    all_index_options = []
    for dbinfo in dbutils.iter_dbinfo(databases, reporter):
        reporter.update("Opening %r ... " % (dbinfo.get_human_name(),))
        database = None
        try:
            database = dbinfo.open_database()
            dataset = database.get_dataset()
        except dbutils.DBError as err:
            reporter.update("")
            reporter.report("Skipping %s: %s"
                            % (dbinfo.get_human_name(), err))
            if database is not None:
                database.close()
            continue
        reporter.update("")

        name = dbinfo.name
        name_width = max(name_width, len(name))
        name_list.append(name)

        table_sizes = dataset.get_table_sizes(args.recount)
        
        num_compounds = table_sizes.num_compounds
        num_compounds_width = max(num_compounds_width, len(str(num_compounds)))
        num_compounds_list.append(num_compounds)

        num_rules = table_sizes.num_rules
        num_rules_width = max(num_rules_width, len(str(num_rules)))
        num_rules_list.append(num_rules)

        num_pairs = table_sizes.num_pairs
        num_pairs_width = max(num_pairs_width, len(str(num_pairs)))
        num_pairs_list.append(num_pairs)

        num_envs = table_sizes.num_rule_environments
        num_envs_width = max(num_envs_width, len(str(num_envs)))
        num_envs_list.append(num_envs)
        
        num_stats = table_sizes.num_rule_environment_stats
        num_stats_width = max(num_stats_width, len(str(num_stats)))
        num_stats_list.append(num_stats)
    
        title = dataset.title
        title_width = max(title_width, len(title))
        titles.append(title)

        prop_names = dataset.get_property_names()
        if prop_names:
            s = " ".join(prop_names)
        else:
            s = "<none>"
        property_names.append(s)

        all_fragment_options.append(dataset.fragment_options_str)
        all_index_options.append(dataset.index_options_str)

    fmt = "%-{}s %-{}s %-{}s %-{}s %-{}s %-{}s  %-{}s Properties".format(
        name_width, num_compounds_width, num_rules_width, num_pairs_width,
        num_envs_width, num_stats_width, title_width)
    fancy_title = " Title ".center(title_width, "-")
    fancy_title = "|" + fancy_title[1:-1] + "|"
    print(fmt % ("Name".center(name_width), "#cmpds".center(num_compounds_width), "#rules".center(num_rules_width),
                 "#pairs".center(num_pairs_width), "#envs".center(num_envs_width), "#stats".center(num_stats_width),
                 fancy_title))
    
    fmt = "%{}s %{}d %{}d %{}d %{}d %{}d  %-{}s %s".format(
        name_width, num_compounds_width, num_rules_width, num_pairs_width, num_envs_width, num_stats_width, title_width)
    
    prefix = " "*num_compounds_width
    for (name, num_compounds, num_rules, num_pairs, num_envs, num_stats, title,
         names_and_counts, fragment_options, index_options) in zip(
            name_list, num_compounds_list, num_rules_list, num_pairs_list, num_envs_list, num_stats_list, titles,
            property_names, all_fragment_options, all_index_options):
        print(fmt % (name, num_compounds, num_rules, num_pairs, num_envs, num_stats, title, names_and_counts))
        if args.all:
            creation_date = dataset.creation_date
            creation_date_str = creation_date.isoformat(" ")
            print(prefix + "Created:", creation_date_str)
            
            s = " "  # Always have a trailing space
            for property_name, count in dataset.get_property_names_and_counts():
                s += "%s/%s " % (count, property_name)
            if s == " ":
                s = "(no properties)"
            else:
                s = s[:-1]  # strip the trailing space
            print(prefix + "  #compounds/property:", s)
            
            print(prefix + "  #smiles for rules: %d  for constants: %d"
                  % (dataset.get_num_rule_smiles(), dataset.get_num_constant_smiles()))
            
            options = json.loads(fragment_options)
            print(prefix + "  Fragment options:")
            for k, v in sorted(options.items()):
                print(prefix + "    %s: %s" % (k, v))
    
            options = json.loads(index_options)
            print(prefix + "  Index options:")
            for k, v in sorted(options.items()):
                print(prefix + "    %s: %s" % (k, v))
    

# mmpdb create_index <filename>
def create_index_command(parser, args):
    mmpdb = dbutils.open_database_from_args_or_exit(args)
    with mmpdb.atomic():
        schema.create_index(mmpdb)
        

# mmpdb drop_index <filename>
def drop_index_command(parser, args):
    mmpdb = dbutils.open_database_from_args_or_exit(args)
    with mmpdb.atomic():
        schema.drop_index(mmpdb)
    mmpdb.execute("VACUUM")

##### loadprops/reaggregate

def reaggregate_properties(dataset, property_name_ids, compound_values_for_property_name_id,
                           cursor, reporter):
    # Mapping from rule environment id to rule environment statistics id
    ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print("Starting reaggregate_properties at " + ts)
    print("Starting reaggregate_properties at " + ts, file=sys.stderr)

    reporter.update("Computing aggregate statistics")
    t1 = time.time()
    num_pairs = dataset.get_num_pairs(cursor=cursor)
    
    all_pairs = dataset.iter_pairs(cursor=cursor)
    all_pairs = reporter.progress(all_pairs, "Computing aggregate statistics", num_pairs)

    stats_info = []
    seen_rule_environment_ids = set()
    for rule_environment_id, rule_environment_pairs in itertools.groupby(
            all_pairs, (lambda pair: pair.rule_environment_id)):

        seen_rule_environment_ids.add(rule_environment_id)
        rule_environment_pairs = list(rule_environment_pairs)  # now a list, not iterator
        
        for property_name_id in property_name_ids:
            deltas = []
            compound_values = compound_values_for_property_name_id[property_name_id]
            for pair in rule_environment_pairs:
                value1 = compound_values.get(pair.compound1_id, None)
                if value1 is None:
                    continue
                value2 = compound_values.get(pair.compound2_id, None)
                if value2 is None:
                    continue
                deltas.append(value2-value1)
            if deltas:
                stats = index_algorithm.compute_aggregate_values(deltas)
                stats_info.append( (rule_environment_id, property_name_id, stats) )

    t2 = time.time()
    diff = t2 - t1
    ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print("Finished aggregating statistics at " + ts)
    print("Finished aggregating statistics at " + ts, file=sys.stderr)
    print("Aggregate statistics computing took (in seconds): " + str(diff))
    # Need to figure out if the statistics exist or need to be created
    reporter.report("Generated %d rule statistics (%d rule environments, %d properties)"
                    % (len(stats_info), len(seen_rule_environment_ids), len(property_name_ids)))
    reporter.update("Getting information about which rule statistics exist...")

    ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print("Reconnecting in reaggregate_properties at: " + ts, file=sys.stderr)
    # Reconnect, because connection timeouts have occurred if the statistics aggregation takes longer than 3 hours
    dataset.mmpa_db.db.close()
    time_last_closed = time.time()
    dataset.mmpa_db.db.connect()
    cursor = dataset.mmpa_db.get_cursor()

    existing_stats_ids = dataset.get_rule_environment_statistics_mapping(
        property_name_ids, cursor=cursor)

    stats_info_progress = reporter.progress(
        stats_info, "Updating statistics table", len(stats_info))
    seen_stats_ids = set()
    num_updated = num_added = 0
    # Write in batches when using a DB server
    if type(dataset.mmpa_db.db) in [peewee.OracleDatabase, peewee.CustomPostgresqlDatabase]:
        ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print("Started batch statistics writing at " + ts)
        print("Started batch statistics writing at " + ts, file=sys.stderr)
        batch_size = 100000
        for (rule_environment_id, property_name_id, stats) in stats_info_progress:
            key = (rule_environment_id, property_name_id)
            stats_id = existing_stats_ids.get(key, None)
            # Need to reformat very large numbers, and very small floating point numbers, otherwise can get errors with Oracle / Postgres
            stats = dataset.validate_statistics(stats)
            if stats_id is not None:
                dataset.update_rule_environment_statistics(stats_id, stats)
                seen_stats_ids.add(stats_id)
                num_updated += 1
                if num_updated % batch_size == 0:
                    dataset.batch_update_re_stats()
                    # Commit after each batch, trying to avoid enormous inserts that can take up all temp space before final commit
                    dataset.mmpa_db.db.commit()
                    if time.time() - time_last_closed > 3600:
                        dataset.mmpa_db.db.close()
                        time_last_closed = time.time()
                        dataset.mmpa_db.db.connect()
            else:
                dataset.add_rule_environment_statistics(rule_environment_id, property_name_id, stats)
                num_added += 1
                if num_added % batch_size == 0:
                    dataset.batch_insert_re_stats()
                    # Commit after each batch, trying to avoid enormous inserts that can take up all temp space before final commit
                    dataset.mmpa_db.db.commit()
                    if time.time() - time_last_closed > 3600:
                        dataset.mmpa_db.db.close()
                        time_last_closed = time.time()
                        dataset.mmpa_db.db.connect()
        if len(dataset.re_stats_updates) > 0:
            dataset.batch_update_re_stats()
        if len(dataset.re_stats_inserts) > 0:
            dataset.batch_insert_re_stats()
        # Commit after each batch, trying to avoid enormous inserts that can take up all temp space before final commit
        dataset.mmpa_db.db.commit()
        if time.time() - time_last_closed > 3600:
            dataset.mmpa_db.db.close()
            time_last_closed = time.time()
            dataset.mmpa_db.db.connect()
    else:
        for (rule_environment_id, property_name_id, stats) in stats_info_progress:
            key = (rule_environment_id, property_name_id)
            stats_id = existing_stats_ids.get(key, None)
            if stats_id is not None:
                dataset.update_rule_environment_statistics(stats_id, stats)
                seen_stats_ids.add(stats_id)
                num_updated += 1
            else:
                dataset.add_rule_environment_statistics(rule_environment_id, property_name_id, stats)
                num_added += 1

    t3 = time.time()
    diff = t3 - t2
    ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print("Finished batch statistics writing at " + ts)
    print("Finished batch statistics writing at " + ts, file=sys.stderr)
    print("Updating statistics took (in seconds): " + str(diff))
    
    to_delete = set(existing_stats_ids.values()) - seen_stats_ids
    num_deleted = len(to_delete)
    if to_delete:
        delete_report = reporter.progress(
            iter(to_delete), "Deleting rule statistics", num_deleted)
        while 1:
            ids = list(itertools.islice(delete_report, 0, 1000))
            if not ids:
                break
            dataset.delete_rule_environment_statistics(ids)

    t4 = time.time()
    diff = t4 - t3
    print("Deleting old statistics took (in seconds): " + str(diff))

    reporter.report("Number of rule statistics added: %d updated: %d deleted: %d"
                    % (num_added, num_updated, num_deleted))
        

def reaggregate_properties_v2(dataset, property_name_ids, compound_values_for_property_name_id,
                           cursor, reporter, new_cps):
    
    # This is only designed to work with the current way mmpdb is set up:
    # When adding properties, submit all properties for everything in the database.
    # Do not submit subsets of the total set of intended properties

    # Reworking reaggregate_properties function to be more efficient

    # reaggregate_properties (original function) scans through every pair for every property
    # This is probably wasteful given that in a large multiproperty database,
    # Many pairs will not have all (or even more than a few?) of the properties recorded

    reporter.update("Computing aggregate statistics")

    stats_info = []
    seen_rule_environment_ids = set()

    # For each property, gather compound_ids in list of tuples with 0: [(cmpd1_id, 0), (cmpd2_id, 0) ...)]
    # Oracle has strange tastes, requiring tuples for lists > 1000 elements


    for property_name_id, new_cp_cluster in itertools.groupby(new_cps, lambda x: x[1]):
        print("Starting to find pairs for property_name_id " + str(property_name_id))
        print("time = " + str(time.time()))
        compound_id_0_tuples = []
        for new_cp in new_cp_cluster:
            compound_id_0_tuples.append((new_cp[0], 0))
        tupleized_string = '(' + str(compound_id_0_tuples)[1:-1] + ')'
        # Now we can query DB with these compound_ids, to get pairs where we will be able to calculate
        # a difference in the property. Aside from pairs selected below, there should be no other pairs
        # in the database where we can calculate a difference in the property, because all property data
        # for everything in the database was submitted to the loadprops_command.

        # pairs = dataset.iter_pairs_with_property(tupleized_string)
        print("Fetching pairs from database at " + str(time.time()))
        pairs = dataset.get_pairs_with_property(tupleized_string)
        num_pairs = len(pairs)
        
        print("Finished fetching pairs at " + str(time.time()))
        pairs = reporter.progress(pairs, "Computing statistics for property_name_id " + str(property_name_id))
        for rule_environment_id, rule_environment_pairs in itertools.groupby(
            pairs, (lambda pair: pair.rule_environment_id)):

            seen_rule_environment_ids.add(rule_environment_id)
            rule_environment_pairs = list(rule_environment_pairs)  # now a list, not iterator
            
            deltas = []
            # Current algorithm requires complete list of properties for all compounds in DB, to be submitted in loadprops command
            # Not efficient for when we want to add properties for just a few new compounds, when many compounds are already in DB
            # Could probably update algorithm to combine a DB query for existing properties, with new/updated properties submitted to loadprops
            compound_values = compound_values_for_property_name_id[property_name_id]
            for pair in rule_environment_pairs:
                value1 = compound_values.get(pair.compound1_id, None)
                value2 = compound_values.get(pair.compound2_id, None)
                try:
                    deltas.append(value2-value1)
                except TypeError:
                    print("Something is wrong with the code: property values expected but not found")
                    raise TypeError
            if deltas:
                stats = index_algorithm.compute_aggregate_values(deltas)
                stats_info.append( (rule_environment_id, property_name_id, stats) )

    # Need to figure out if the statistics exist or need to be created
    reporter.report("Generated %d rule statistics (%d rule environments, %d properties)"
                    % (len(stats_info), len(seen_rule_environment_ids), len(property_name_ids)))
    reporter.update("Getting information about which rule statistics exist...")
    existing_stats_ids = dataset.get_rule_environment_statistics_mapping(
        property_name_ids, cursor=cursor)

    stats_info_progress = reporter.progress(
        stats_info, "Updating statistics table", len(stats_info))
    seen_stats_ids = set()
    num_updated = num_added = 0
    for (rule_environment_id, property_name_id, stats) in stats_info_progress:
        key = (rule_environment_id, property_name_id)
        stats_id = existing_stats_ids.get(key, None)
        if stats_id is not None:
            dataset.update_rule_environment_statistics(stats_id, stats)
            seen_stats_ids.add(stats_id)
            num_updated += 1
        else:
            dataset.add_rule_environment_statistics(rule_environment_id, property_name_id, stats)
            num_added += 1

    # If writing to an Oracle DB, we need to executemany with lists generated in the above for loop.
    from .peewee import OracleDatabase
    if isinstance(dataset.mmpa_db.db, OracleDatabase):
        dataset.batch_write_re_stats_to_oracle()

    # The code is currently set up so that properties should be loaded all at once, not incrementally.
    # To tell the DB to remove compound property values currently in the DB,
    # omit these values to be deleted from the property file passed into loadprops command.

    # The below code will delete stats that were in DB before, but which were not recalculated this time,
    # because the compound property values for these stats were not passed into the loadprops command
    to_delete = set(existing_stats_ids.values()) - seen_stats_ids
    num_deleted = len(to_delete)
    if to_delete:
        delete_report = reporter.progress(
            iter(to_delete), "Deleting rule statistics", num_deleted)
        while 1:
            ids = list(itertools.islice(delete_report, 0, 1000))
            if not ids:
                break
            dataset.delete_rule_environment_statistics(ids)

    reporter.report("Number of rule statistics added: %d updated: %d deleted: %d"
                    % (num_added, num_updated, num_deleted))

# mmpdb loadprops <filename> -p <properties_filename>

# When adding properties, submit all properties for everything in the database, in the <properties_filename> file
# Do not submit subsets of the total set of intended properties - this will cause all other properties not in the file to be deleted
# To delete properties, leave them out of the file

def loadprops_command(parser, args):
    ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print("Started loadprops_command at " + ts)
    print("Started loadprops_command at " + ts, file=sys.stderr)
    print("args", file=sys.stderr)
    print(args, file=sys.stderr)
    from . import properties_io
    reporter = command_support.get_reporter(args.quiet)

    db = dbutils.open_database_from_args_or_exit(args)
    c = db.get_cursor()
    dataset = db.get_dataset()

    # Testing reconnection
    ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print("Reconnecting in loadprops_command at: " + ts, file=sys.stderr)
    db.db.close()
    time_last_closed = time.time()
    db.db.connect()
    c = db.get_cursor()

    reporter.report("Using dataset: %s" % (dataset.title,))
        
    if args.properties is None:
        reporter.report("Reading properties from stdin")
        properties_file = sys.stdin
        close = None
        source = "<stdin>"
    else:
        reporter.report("Reading properties from %r" % (args.properties,))
        try:
            properties_file = open(args.properties, "U")
        except IOError as err:
            parser.error("Cannot open properties file: %s" % (err,))
        close = properties_file.close
        source = args.properties

    try:
        try:
            with properties_file:
                properties = properties_io.load_properties(properties_file, reporter)
        except ValueError as err:
            parser.error("Problem reading --properties file %r: %s"
                         % (args.properties, err))
    finally:
        if close is not None:
            close()

    reporter.report("Read %d properties for %d compounds from %r"
                    % (len(properties.property_columns), len(properties.id_column),
                       source))

    # mmpdb loadprops -p <properties_filename> --metadata <metadata_filename>
    if args.metadata is None:
        reporter.report("""
No metadata file was passed with the option --metadata <file>.
In this case, all property values displayed in the UI will be identical to those stored in the database.
All property changes will be represented as deltas by default in the UI.
""")
        metadata = None
    else:
        try:
            metadata_file = open(args.metadata, "U")
        except IOError as err:
            parser.error("Cannot open metadata file: %s" % (err,))
            metadata=None
        close = metadata_file.close
        source = args.metadata

        try:
            try:
                with metadata_file:
                    metadata = properties_io.load_metadata(metadata_file, reporter)
            except ValueError as err:
                parser.error("Problem reading --metadata file %r: %s"
                            % (args.metadata, err))
                metadata=None
        finally:
            if close is not None:
                close()

        reporter.report("Read metadata for %d properties from %r" % (len(metadata), source))

    public_id_to_id = dataset.get_public_id_to_id_table(c)
    
    compound_ids = [public_id_to_id.get(public_id, None) for public_id in properties.id_column]
    num_missing = compound_ids.count(None)
    if num_missing:
        reporter.report("%d compounds from %r are not in the dataset at %r"
                        % (num_missing, source, args.databases[0]))
        ## missing = [public_id for public_id in properties.id_column if public_id not in public_id_to_id]
        ## del missing[6:]
        ## if len(missing) > 5:
        ##     reporter.warning("First five missing records: %r" % (missing[:5],))
        ## else:
        ##     reporter.warning("Missing records: %r" % (missing,))
            

    UPDATE_COMPOUND_PROPERTY_SQL = db.SQL(
        "UPDATE compound_property "
        "      SET value = ?"
        "       WHERE compound_id = ?"
        "         AND property_name_id = ?")
    INSERT_COMPOUND_PROPERTY_SQL = db.SQL(
        "INSERT INTO compound_property (compound_id, property_name_id, value) "
        " VALUES (?, ?, ?)")
    
    with db.atomic():
        # Remember which compound properties exist, so I can tell if a
        # value should replace an existing value or create a new value.
        c.execute(db.SQL(
            "SELECT compound_id, property_name_id from compound_property"))
        seen_properties = dict((key, False) for key in c)
        ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print("Finished downloading compound_property rows at " + ts)
        print("Finished downloading compound_property rows at " + ts, file=sys.stderr)

        compound_values_for_property_name_id = {}
        property_name_ids = []
        
        # Enabling batch statement execution, for Oracle and Postgres

        #from .peewee import OracleDatabase, CustomPostgresqlDatabase
        if type(db.db) in [peewee.OracleDatabase, peewee.CustomPostgresqlDatabase]:
        #if isinstance(db.db, OracleDatabase):
            cp_updates = []
            cp_inserts = []

            #def update_prop(val, cmpd_id, p_name_id, update_list):
                #update_list.append((val, cmpd_id, p_name_id))
            #def insert_prop(val, cmpd_id, p_name_id, insert_list):
                #insert_list.append((cmpd_id, p_name_id, val))

            #def batch_write(cursor, statement, batch_list):
                #start_batch_write_time = time.time()
                #print("Starting batch property writing at " + str(start_batch_write_time))
                #batch_execute(cursor, statement, batch_list)
                #finish_batch_write_time = time.time()
                #print(str(len(batch_list)) + " writes took " + str(finish_batch_write_time - start_batch_write_time) + " seconds to write")

            if type(db.db) == peewee.OracleDatabase:
                insert_statement = INSERT_COMPOUND_PROPERTY_SQL
                update_statement = UPDATE_COMPOUND_PROPERTY_SQL
                def batch_execute(cursor, statement, batch_list):
                    start_batch_write_time = time.time()
                    print("Starting batch property writing at " + str(start_batch_write_time))

                    cursor.executemany(statement, batch_list)

                    finish_batch_write_time = time.time()
                    print(str(len(batch_list)) + " writes took " + str(finish_batch_write_time - start_batch_write_time) + " seconds to write")
                batch_insert = batch_execute
                batch_update = batch_execute

            elif type(db.db) == peewee.CustomPostgresqlDatabase:
                from psycopg2.extras import execute_values
                update_statement = """
UPDATE compound_property 
SET value = update_batch.value
FROM (VALUES %s) AS update_batch (value, compound_id, property_name_id)
WHERE compound_property.compound_id = update_batch.compound_id
AND compound_property.property_name_id = update_batch.property_name_id
"""
                def batch_update(cursor, statement, batch_list):
                    #start_batch_write_time = time.time()
                    #print("Starting batch property writing at " + str(start_batch_write_time))

                    execute_values(cursor, statement, batch_list)

                    #finish_batch_write_time = time.time()
                    #print(str(len(batch_list)) + " writes took " + str(finish_batch_write_time - start_batch_write_time) + " seconds to write")

                insert_statement = b"INSERT INTO compound_property (compound_id, property_name_id, value) VALUES "
                def batch_insert(cursor, statement, batch_list):
                    #start_batch_write_time = time.time()
                    #print("Starting batch property writing at " + str(start_batch_write_time))

                    insert_values = b','.join(cursor.mogrify("(%s,%s,%s)", x) for x in batch_list)
                    cursor.execute(statement + insert_values)

                    #finish_batch_write_time = time.time()
                    #print(str(len(batch_list)) + " writes took " + str(finish_batch_write_time - start_batch_write_time) + " seconds to write")

            for property_name, property_values in properties.iter_properties():
                property_name_id = dataset.get_or_add_property_name(property_name)
                property_name_ids.append(property_name_id)
                #print(str(property_name_id))
                #reporter.report("Loading property %r (id %d)" % (property_name.name, property_name.id))
                ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                print("Starting compound_property writes for " + property_name + " at " + ts)
                print("Starting compound_property writes for " + property_name + " at " + ts, file=sys.stderr)

                num_created = num_updated = 0
                compound_values_for_property_name_id[property_name_id] = compound_values = {}
            
                for compound_id, value in zip(compound_ids, property_values):
                    # Prevent the batch execution lists from getting too large and taking too much memory
                    if len(cp_inserts) > 0 and len(cp_inserts) % 100000 == 0:
                        batch_insert(c, insert_statement, cp_inserts)
                        cp_inserts = []
                        # Commit after each batch, trying to avoid enormous inserts that can take up all temp space before final commit
                        db.db.commit()
                        if time.time() - time_last_closed > 3600:
                            db.db.close()
                            time_last_closed = time.time()
                            db.db.connect()
                            c = db.get_cursor()
                    if len(cp_updates) > 0 and len(cp_updates) % 100000 == 0:
                        batch_update(c, update_statement, cp_updates)
                        cp_updates = []
                        # Commit after each batch, trying to avoid enormous inserts that can take up all temp space before final commit
                        db.db.commit()
                        if time.time() - time_last_closed > 3600:
                            db.db.close()
                            time_last_closed = time.time()
                            db.db.connect()
                            c = db.get_cursor()

                    if compound_id is None:
                        # Property specified but not in database
                        continue
                    if value is None:
                        # Property value is "*", meaning it isn't available
                        continue
                    key = (compound_id, property_name_id)
                    if key in seen_properties:
                        seen_properties[key] = True
                        num_updated += 1
                        #update_prop(value, compound_id, property_name_id)
                        cp_updates.append((value, compound_id, property_name_id))
                    else:
                        num_created += 1
                        #insert_prop(value, compound_id, property_name_id)
                        cp_inserts.append((compound_id, property_name_id, value))
                    compound_values[compound_id] = value
                reporter.report(
                    "Imported %d %r records (%d new, %d updated)."
                    % (num_updated + num_created, property_name, num_created, num_updated))

            if len(cp_inserts) > 0:
                batch_insert(c, insert_statement, cp_inserts)
                cp_inserts = []
            if len(cp_updates) > 0:
                batch_update(c, update_statement, cp_updates)
                cp_updates = []
            # Commit after each batch, trying to avoid enormous inserts that can take up all temp space before final commit
            db.db.commit()

        else:
            def update_prop(val, cmpd_id, p_name_id, cursor=c, statement=UPDATE_COMPOUND_PROPERTY_SQL):
                cursor.execute(statement, (val, cmpd_id, p_name_id))
            def insert_prop(val, cmpd_id, p_name_id, cursor=c, statement=INSERT_COMPOUND_PROPERTY_SQL):
                cursor.execute(statement, (cmpd_id, p_name_id, val))

            for property_name, property_values in properties.iter_properties():
                property_name_id = dataset.get_or_add_property_name(property_name)
                property_name_ids.append(property_name_id)
                #print(str(property_name_id))
                #reporter.report("Loading property %r (id %d)" % (property_name.name, property_name.id))

                num_created = num_updated = 0
                compound_values_for_property_name_id[property_name_id] = compound_values = {}
            
                for compound_id, value in zip(compound_ids, property_values):
                    if compound_id is None:
                        # Property specified but not in database
                        continue
                    if value is None:
                        # Property value is "*", meaning it isn't available
                        continue
                    key = (compound_id, property_name_id)
                    if key in seen_properties:
                        seen_properties[key] = True
                        num_updated += 1
                        update_prop(value, compound_id, property_name_id)
                    else:
                        num_created += 1
                        insert_prop(value, compound_id, property_name_id)
                    compound_values[compound_id] = value
                reporter.report(
                    "Imported %d %r records (%d new, %d updated)."
                    % (num_updated + num_created, property_name, num_created, num_updated))
        
        ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print("Finished compound_property writes at " + ts)
        print("Finished compound_property writes at " + ts, file=sys.stderr)
        # The code is currently set up so that properties should be loaded all at once, not incrementally.
        # To tell the DB to remove compound property values currently in the DB,
        # omit these values to be deleted from the property file passed into loadprops command.

        # Remove existing compound properties where the property name was in the
        # properties file but the where the file did not specify a value.
        properties_to_delete = [key for key, was_updated in seen_properties.items()
                                if not was_updated and key[1] in property_name_ids]
        if properties_to_delete:
            dataset.delete_compound_properties(properties_to_delete)
        
        reaggregate_properties(dataset, property_name_ids, compound_values_for_property_name_id,
                                cursor=c, reporter=reporter)

        #db.db.close()
        #db.db.connect()
        # We may have reset the db connection during reaggregate_properties, therefore get a fresh cursor here
        c = db.get_cursor()

        # Check if any of the properties are completely gone
        if properties_to_delete:
            for property_name_id in property_name_ids:
                n = dataset.get_num_compound_properties(property_name_id, cursor=c)
                if n == 0:
                    dataset.delete_property_name_id(property_name_id, cursor=c)

        # Update the environment statistics
        reporter.update("Updating environment statistics count ...")
        if type(db.db) == peewee.CustomPostgresqlDatabase:
            c.execute("SELECT count(*) from rule_environment_statistics")
            num = schema._get_one(c)
            c.execute("UPDATE dataset SET num_rule_environment_stats=%s", (num,))
        else:
            num = schema._get_one(c.execute("SELECT count(*) from rule_environment_statistics"))
            c.execute(db.SQL("UPDATE dataset SET num_rule_environment_stats=?"), (num,))

        seen_props = set()
        if metadata is not None:
            for prop in metadata:
                seen_props.add(prop)
                insert_values = [metadata[prop][key] for key in ('base', 'unit', 'display_name', 'display_base', 'display_unit', 'change_displayed')]
                c.execute("""UPDATE property_name SET (base, unit, display_name, display_base, display_unit, change_displayed) = 
(%s, %s, %s, %s, %s, %s) WHERE property_name.name = %s""", insert_values + [prop])

        for prop in set(properties.property_names) - seen_props:
            c.execute("""UPDATE property_name SET (base, display_name, display_base, change_displayed) = 
('raw', %s, 'raw', 'delta') WHERE property_name.name = %s""", (prop, prop))
        
        reporter.update("Commiting changed ...")
        # It seems like the actual commit happens when peewee.py exits the context manager for the DB at this point
        # If an exception occurred earlier, then it seems like it would have instead rolled back the writes
        ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print("Exiting db.atomic and committing writes at " + ts)
        print("Exiting db.atomic and committing writes at " + ts, file=sys.stderr)

    ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print("Exited db.atomic at " + ts)
    print("Exited db.atomic at " + ts, file=sys.stderr)
                
    reporter.report("Loaded all properties and re-computed all rule statistics.")

            
def reaggregate_command(parser, args):
    reporter = command_support.get_reporter(args.quiet)

    db = dbutils.open_database_from_args_or_exit(args)
    c = db.get_cursor()
    dataset = db.get_dataset()
    
    with db.atomic():
        # Which properties to reaggregate?
        property_names = command_support.get_property_names_or_error(parser, args, dataset)
        reporter.report("Reaggregating %d properties" % (len(property_names),))

        # Convert them into ids
        property_name_ids = [dataset.get_property_name_id(name) for name in property_names]

        # Get their values
        compound_values_for_property_name_id = dict(
            (property_name_id, dataset.get_property_values(property_name_id))
             for property_name_id in property_name_ids)

        # Reaggregate
        reaggregate_properties(dataset, property_name_ids, compound_values_for_property_name_id,
                               cursor=c, reporter=reporter)
                

# mmpdb smicat [--input-smiles] [-o filename] filename
def smicat_command(parser, args):
    reporter = command_support.get_reporter(args.quiet)

    db = dbutils.open_database_from_args_or_exit(args)
    c = db.get_cursor()
    dataset = db.get_dataset()
    input_smiles = args.input_smiles

    output_filename = args.output
    with fileio.open_output(output_filename, output_filename) as outfile:
        for compound_row in dataset.iter_compounds(cursor=c):
            if input_smiles:
                smiles = compound_row.input_smiles
            else:
                smiles = compound_row.clean_smiles
            outfile.write("%s %s\n" % (smiles, compound_row.public_id))

# mmpdb propcat [--property <name>]* [--all] [-o filename] filename
def propcat_command(parser, args):
    reporter = command_support.get_reporter(args.quiet)

    db = dbutils.open_database_from_args_or_exit(args)
    c = db.get_cursor()
    dataset = db.get_dataset()

    show_all = args.all

    property_names = command_support.get_property_names_or_error(parser, args, dataset)

    property_values_list = []
    for property_name in property_names:
        property_name_id = dataset.get_property_name_id(property_name)
        property_values_list.append(dataset.get_property_values(property_name_id))

    output_filename = args.output
    with fileio.open_output(output_filename, output_filename) as outfile:
        print("ID", *property_names, sep="\t", file=outfile)
        for compound_row in dataset.iter_compounds(cursor=c):
            columns = [compound_row.public_id]
            is_empty = True
            for property_values in property_values_list:
                value = property_values.get(compound_row.id, None)
                if value is None:
                    columns.append("*")
                else:
                    columns.append(value)
                    is_empty = False
            if show_all or not is_empty:
                print(*columns, sep="\t", file=outfile)
        


# mmpdb rulecat? or would paircat be more useful? both?
