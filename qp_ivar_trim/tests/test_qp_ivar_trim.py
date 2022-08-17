# -----------------------------------------------------------------------------
# Copyright (c) 2020--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
from unittest import main
from qiita_client.testing import PluginTestCase
from os import remove, environ, dirname
from os.path import exists, isdir, join
from shutil import rmtree, copyfile
from tempfile import mkdtemp
from json import dumps
from itertools import zip_longest
from functools import partial

from qp_ivar_trim import plugin
from qp_ivar_trim.utils import plugin_details
from qp_ivar_trim.qp_ivar_trim import (
    get_dbs_list, _generate_commands,
    ivar_trim_to_array, QC_REFERENCE, IVAR_TRIM_CMD)


class IvarTrimTests(PluginTestCase):
    def setUp(self):
        plugin("https://localhost:21174", 'register', 'ignored')

        out_dir = mkdtemp()
        self.maxDiff = None
        self.out_dir = out_dir
        self.dbs = get_dbs_list()
        self.db_path = QC_REFERENCE  # need to make envrionment variable
        self.params = {'reference': 'artifacts', 'threads': 2}
        self._clean_up_files = []
        self._clean_up_files.append(out_dir)

    def tearDown(self):
        for fp in self._clean_up_files:
            if exists(fp):
                if isdir(fp):
                    rmtree(fp)
                else:
                    remove(fp)

#    def test_get_dbs_list(self):
#        dbs = get_dbs_list()
#        self.assertCountEqual(dbs, ['artifacts.mmi', 'empty.mmi'])
        # might need to change this

    def test_generate_commands(self):
        params = {'database': 'artifacts', 'nprocs': 2,
                  'out_dir': '/foo/bar/output'}
        # need to change these to bam
        BAM_file = ['CALM_SEP_001970_03_S265_L001.sorted.bam',
                    'CALM_SEP_001970_03_S265_L002.sorted.bam']
        obs = _generate_commands(BAM_file, params['database'],
                                 params['nprocs'], params['out_dir'])
        cmd = IVAR_TRIM_CMD.format(**params)
        ecmds = [cmd % (bam, bam)
                 for bam in zip_longest(BAM_file)]
        eof = [(f'{params["out_dir"]}/{f}', 'bam')
               for bam in sorted(BAM_file)]
        # for f in sorted(rev_seqs):
        #    eof.append((f'{params["out_dir"]}/{f}',
        #  'raw_reverse_seqs'))
        self.assertCountEqual(obs[0], ecmds)
        self.assertCountEqual(obs[1], eof)

        # cmd = COMBINED_CMD_SINGLE.format(**params)
        # ecmds = [cmd % (f, f) for
        # f in fwd_seqs]
        # eof = [(f'{params["out_dir"]}/{f}',
        #  'raw_forward_seqs')
        #       for f in sorted(fwd_seqs)]
        # self.assertCountEqual(obs[0], ecmds)
        # self.assertCountEqual(obs[1], eof)

if __name__ == '__main__':
            main()

