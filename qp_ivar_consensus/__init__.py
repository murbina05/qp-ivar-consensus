# -----------------------------------------------------------------------------
# Copyright (c) 2020, Qiita development team.
#
# Distributed under the terms of the BSD 3-clause License License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from qiita_client import QiitaPlugin, QiitaCommand
from .qp_ivar_consensus import get_dbs_list, ivar_consensus
from .utils import plugin_details
from os.path import splitext


THREADS = 15


# Initialize the plugin
plugin = QiitaPlugin(**plugin_details)

# Define the command
dbs = get_dbs_list()
dbs_without_extension = [splitext(db)[0] for db in dbs]
dbs_defaults = ', '.join([f'"{x}"' for x in dbs_without_extension])
req_params = {'input': ('artifact', ['per_sample_FASTQ'])}
opt_params = {
    'reference': [
        f'choice:["None", {dbs_defaults}]', dbs_without_extension],
    'threads': ['integer', f'{THREADS}']}

outputs = {'Filtered files': 'bam'}
default_params = {
    'auto-detect adapters only filtering [not recommended]': {
        'reference': "None", 'threads': THREADS}}
for db in dbs_without_extension:
    name = f'auto-detect adapters and {db} + phix filtering'
    default_params[name] = {'reference': db, 'threads': THREADS}

IVAR_consensus_cmd = QiitaCommand(
    'consensus reads in aligned BAM ', "",
    ivar_consensus, req_params, opt_params, outputs, default_params)
# i have no idea what to put here :) ^^^^
plugin.register_command(IVAR_consensus_cmd)
