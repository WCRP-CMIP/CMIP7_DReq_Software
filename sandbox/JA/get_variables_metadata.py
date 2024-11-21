#!/usr/bin/env python
'''
Extract metadata of CMOR variables and write to json
'''

import sys
import json
add_paths = ['../MS/dreq_api/', '../GR/']
for path in add_paths:
    if path not in sys.path:
        sys.path.append(path)
import dreq_content as dc
import dreq_query as dq
import dreq_classes

from collections import OrderedDict

# from importlib import reload
# reload(dq)
# reload(dc)

###############################################################################
# Load data request content

use_dreq_version = 'v1.0beta'

# Download specified version of data request content (if not locally cached)
dc.retrieve(use_dreq_version)
# Load content into python dict
content = dc.load(use_dreq_version)

###############################################################################
# Retrive info about variables

base = dq.create_dreq_tables_for_variables(content)

Vars = base['Variables']

# Choose which table to use for freqency
# freq_table_name = 'Frequency'  # not available in v1.0beta release export, need to use CMIP7 or CMIP6 one instead
# freq_table_name = 'CMIP7 Frequency'
# freq_table_name = 'CMIP6 Frequency (legacy)'

try_freq_table_name = []
try_freq_table_name.append('Frequency')
try_freq_table_name.append('CMIP7 Frequency')
try_freq_table_name.append('CMIP6 Frequency (legacy)')

for freq_table_name in try_freq_table_name:
    freq_attr_name = dreq_classes.format_attribute_name(freq_table_name)
    # assert freq_attr_name in Vars.attr2field, 'attribute not found: ' + freq_attr_name
    if freq_attr_name not in Vars.attr2field:
        continue
    if 'frequency' not in Vars.attr2field:
        # code below assumes a variable's frequency is given by its "frequency" 
        Vars.rename_attr(freq_attr_name, 'frequency')
    if freq_table_name in base:
        Frequency = base[freq_table_name]
    break

SpatialShape = base['Spatial Shape']
Dimensions = base['Coordinates and Dimensions']
TemporalShape = base['Temporal Shape']
CellMethods = base['Cell Methods']
PhysicalParameter = base['Physical Parameters']

CFStandardName = None
if 'CF Standard Names' in base:
    CFStandardName = base['CF Standard Names']

# Use compound name to look up record id of each variable in the Vars table
var_name_map = {record.compound_name : record_id for record_id, record in Vars.records.items()}
assert len(var_name_map) == len(Vars.records), 'compound names do not uniquely map to variable record ids'

# Dicts to store the results
all_var_info = {}

for var in Vars.records.values():

    var_info = {}

    if isinstance(var.frequency[0], str):
        # retain this option for non-consolidated raw export?
        assert isinstance(var.frequency, list)
        frequency = var.frequency[0]
    else:
        link = var.frequency[0]
        freq = Frequency.get_record(link)
        frequency = freq.name

    link = var.temporal_shape[0]
    temporal_shape = TemporalShape.get_record(link)

    if hasattr(var, 'cell_methods'):
        assert len(var.cell_methods) == 1
        link = var.cell_methods[0]
        cell_methods = CellMethods.get_record(link).cell_methods
    else:
        cell_methods = ''

    # get the 'Spatial Shape' record, which contains info about dimensions
    assert len(var.spatial_shape) == 1
    link = var.spatial_shape[0]
    spatial_shape = SpatialShape.get_record(link)

    var_dims = []
    dims = None
    if hasattr(spatial_shape, 'dimensions'):
        for link in spatial_shape.dimensions:
            dims = Dimensions.get_record(link)
            var_dims.append(dims.name)

    # Get CF standard name, if it exists
    # record_id = var.cf_standard_name_from_physical_parameter[0]  # not a real link! 
    # phys_param = PhysicalParameter.get_record(record_id)
    link = var.physical_parameter[0]
    phys_param = PhysicalParameter.get_record(link)
    if hasattr(phys_param, 'cf_standard_name'):
        if isinstance(phys_param.cf_standard_name, str):
            # retain this option for non-consolidated raw export?
            var_info.update({
                'CF standard name' : phys_param.cf_standard_name,
            })
        else:
            link = phys_param.cf_standard_name[0]
            cfsn = CFStandardName.get_record(link)
            var_info.update({
                'CF standard name' : cfsn.name,
            })
    else:
        var_info.update({
            'CF standard name (proposed)' : phys_param.proposed_cf_standard_name,
        })

    var_info.update({
        'units' : phys_param.units,
        'cell_methods' : cell_methods,
        'dimensions' : ' '.join(var_dims),
        'frequency' : frequency,
        'spatial_shape' : spatial_shape.name,
        'temporal_shape' : temporal_shape.name,


        # 'hor_label_dd' : spatial_shape.hor_label_dd,
        # 'vertical_label_dd' : spatial_shape.vertical_label_dd,
        # 'temporal_brand' : temporal_shape.brand,
    })

    # var_name = var.compound_name
    var_name = dq.get_unique_var_name(var)
    assert var_name not in all_var_info, 'non-unique variable name: ' + var_name
    all_var_info[var_name] = var_info


# Sort the all-variables dict
d = OrderedDict()
for var_name in sorted(all_var_info.keys(), key=str.lower):
    d[var_name] = all_var_info[var_name]
all_var_info = d
del d

d = OrderedDict({
    'Header' : OrderedDict({
        "dreq version": use_dreq_version,
        "Comment" : "Provisional selected information on all variables requested by at least one Opportunity. To be superseded by a systematic collection of all metadata needed by CMOR and other relevant information on each variable (likely using a single file for each variable). As an interim partial solution this file provides a very basic subset of information about the requested variables. Each variable is idenfied by a compound name comprised of a CMIP6-era table name and a short variable name."
    }),
    dq.UNIQUE_VAR_NAME : all_var_info,
})

filepath = 'all_var_info.json'
with open(filepath, 'w') as f:
    json.dump(d, f, indent=4)
    print(f'wrote {filepath} for {len(all_var_info)} variables')

