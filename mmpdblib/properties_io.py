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


class Properties(object):
    def __init__(self, id_header, id_column, property_names, property_columns,
                 property_table):
        self.id_header = id_header
        self.id_column = id_column
        self.property_names = property_names
        self.property_columns = property_columns
        self.property_table = property_table
        
    def get_ids(self):
        return self.id_column
    
    def get_property_values(self, id):
        return self.property_table[id]
        
    def iter_properties(self):
        for name, column in zip(self.property_names, self.property_columns):
            yield name, column


# If a line contains a tab then it's tab-delimited.  Otherwise it's
# whitespace delimited.  (Originally it was whitespace delimited. Then
# I decided to support identifiers which contained a space in
# them. This seemed like a hacky-but-good-enough solution.)
def _split(line):
    if "\t" in line:
        return line.rstrip("\n").split("\t")
    return line.split()


def load_properties(properties_file, reporter):
    try:
        header_line = next(properties_file)
    except StopIteration:
        raise ValueError("First line of the properties file must contain the header")
    header_names = _split(header_line)
    if not header_names:
        raise ValueError("The properties file must contain at least one column name, for the id")
    if header_names[0] not in ("id", "ID", "Name", "name"):
        reporter.warning("the identifier column in the properties file (column 1) has "
                         "a header of %r; should be 'id', 'ID', 'Name', or 'name'" % (header_names[0],))

    seen = set()
    for header_name in header_names:
        if header_name in seen:
            raise ValueError(
                "Duplicate header %r found. A property name may not be listed more than once."
                % (header_name,))
        seen.add(header_name)
        
    n = len(header_names)

    id_column = []
    property_table = {}
    property_rows = []
    for lineno, line in enumerate(properties_file, 2):
        fields = _split(line)
        if len(fields) != n:
            raise ValueError("Line %d has %d fields but the header has %d"
                             % (lineno, len(fields), n))
        float_fields = []
        try:
            for field in fields[1:]:
                if field == "*":
                    float_fields.append(None)
                else:
                    float_fields.append(float(field))
        except ValueError:
            raise ValueError("Line %d value %r cannot be converted to a float"
                             % (lineno, field))

        id = fields[0]
        id_column.append(id)
        property_table[id] = float_fields
        property_rows.append(float_fields)

    if property_rows:
        property_columns = list(zip(*property_rows))
    else:
        property_columns = [[] for _ in header_names]
    return Properties(header_names[0], id_column, header_names[1:], property_columns,
                      property_table)

def load_metadata(metadata_file, reporter):

    header_names = ['property_name', 'base', 'unit', 'display_name', 'display_base', 'display_unit', 'change_displayed']
    try:
        header_line = next(metadata_file)
    except StopIteration:
        raise ValueError("First line of the metadata file must contain the header")
    input_header_names = _split(header_line)
    if not input_header_names or input_header_names != header_names:
        raise ValueError("""The header values must be, in this order:
property_name	base	unit	display_name	display_base	display_unit	change_displayed""")

    metadata = {}
    for line in metadata_file:
        fields = _split(line)
        if len(fields) < 7 or len(fields) > 7:
            reporter.warning(f"The line for {fields[0]} needs 7 fields, but has {len(fields)}, skipping this line...")
            continue
        elif fields == header_names:
            continue
        fields = [str(x) for x in fields]
        metadata[fields[0]] = dict(zip(header_names[1:], fields[1:]))

    for prop in metadata:
        if metadata[prop]['base'] not in ('raw', 'log', 'negative_log'):
            metadata[prop]['base'] = 'raw'

        elif metadata[prop]['base'] in ('log', 'negative_log'):
            if metadata[prop]['display_base'] in ('log', 'negative_log'):

                reporter.warning(f"""The only base conversions supported are log or negative_log to raw.
However, conversion from {metadata[prop]['base']} to {metadata[prop]['display_base']} was attempted.
{prop} values displayed in the matcher UI will not undergo any changes relative to values in the properties file.
""")
                metadata[prop]['base'], metadata[prop]['display_base'] = 'raw', 'raw'
                metadata[prop]['unit'], metadata[prop]['display_unit'] = None, None

        if metadata[prop]['unit'] not in ('M', 'uM') or metadata[prop]['display_unit'] not in ('M', 'uM'):

            reporter.warning(f"""The only unit conversions supported are from M to uM, or uM to M.
However, conversion from {metadata[prop]['unit']} to {metadata[prop]['display_unit']} (or vice-versa) was attempted.
{prop} will not undergo unit conversion relative to values in the properties file.
""")
            metadata[prop]['unit'], metadata[prop]['display_unit'] = None, None

        if metadata[prop]['display_name'] in ('', '*'):
            reporter.warning(f"Invalid display_name entered for {prop}, setting display_name = property_name")
            metadata[prop]['display_name'] = prop

        if metadata[prop]['display_base'] not in ('raw'):
            metadata[prop]['display_base'] = 'raw'

        if metadata[prop]['change_displayed'] not in ('fold-change', 'delta'):
            reporter.warning(f"Only fold-change and delta are supported as values of change_displayed, setting default value of delta for {prop}")
            metadata[prop]['change_displayed'] = 'delta'

    return metadata