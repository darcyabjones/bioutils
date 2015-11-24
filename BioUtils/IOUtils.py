
__program__ = "IOUtils"
__version__ = "0.1.0"
__author__ = "Darcy Jones"
__date__ = "17 May 2014"
__author_email__ = "darcy.ab.jones@gmail.com"
__license__ = """
##############################################################################

    Copyright (C) 2014  Darcy Jones

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.

##############################################################################
"""

############################ Import all modules ##############################

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import os
import subprocess

from Bio.Application import _Option
from Bio.Application import _Switch
from Bio.Application import _ArgumentList
from Bio.Application import AbstractCommandline

class MugsyCommandline(AbstractCommandline):
    """ Command line wrapper for the multiple alignment program MUGSY.

    http://mugsy.sourceforge.net/
    """

    def __init__(self, *args, cmd='mugsy', **kwargs):
        """ . """
        self.parameters = [
            _Option(
                ['--prefix', 'prefix'],
                'prefix for output files',
                filename=True,
                equate=False,
                ),
            _Option(
                ['--directory', 'directory'],
                ('directory used to store output and temporary'
                 'files. Must be a absolute path'),
                filename=True,
                equate=False,
                ),
            _Option(
                ['--minlength', 'minlength'],
                ('minimum span of an aligned region in a colinear\n'
                 'block (bp). This is used by the segmentation step\n'
                 'synchain-mugsy. Default is 30bp.'),
                checker_function=lambda x: isinstance(x, int),
                equate=False,
                ),
            _Option(
                ['--duplications', 'duplications'],
                '1 - Detect and report duplications. 0 - Skip. Default is 0.',
                checker_function=lambda x: x in (0, 1),
                equate=False,
                ),
            _Option(
                ['--nucmeropts', 'nucmeropts'],
                ('options passed through to the Nucmer\n'
                 'package. Eg. -nucmeropts "-l 15" sets the minimum MUM length\n'
                 'in NUCmer to 15. See the Nucmer documentation at\n'
                 'http://mummer.sf.net for more information. Default is -l 15.'),
                checker_function=lambda x: "'" in x or '"' in x,
                equate=False,
                ),
            _Switch(
                ['-allownestedlcbs', 'allownestedlcbs'],
                ('Default=false. Places each multi-genome anchor\n'
                 'in exactly one LCB; the longest spanning LCB'),
                ),
            _Switch(
                ['-plot', 'plot'],
                ('output genome dot plots in GNUplot format. Overlays LCBS\n'
                 'onto pairwise plots from mummerplot. Display of draft\n'
                 'genomes in these plots is not supported.'),
                ),
            _Switch(
                ['-fullsearch', 'fullsearch'],
                ('Run a complete all pairs Nucmer search with each \n'
                 'sequence as a reference and query (n^2-1 total searches).\n'
                 'Default is one direction only (n^2-1/2 searches).'),
                ),
            _Option(
                ['-refine', 'refine'],
                ('run an second iteration of Mugsy on each LCB to refine the \n'
                 'alignment using either Mugsy (--refine mugsy), FSA (--refine\n'
                 'fsa), Pecan (--refine pecan), MLAGAN (--refine mlagan).\n'
                 'Requires necessary tools are in your path:\n'
                 'fsa: fsa\n'
                 'pecan: muscle,exonerate, in the path. classpath set '
                 'for bp.pecan.Pecan.\n'
                 'mlagan: mlagan.sh'),
                checker_function=lambda x: x in ('mugsy', 'fsa',
                                                 'pecan', 'mlagan'),
                equate=False,
                ),
            _ArgumentList(
                ['input'],
                ('Input is one or more (multi)FASTA files, one per genome.\n'
                 'Each file should contain all the sequences for a single\n'
                 'organism/species. The filename is used as the genome name.\n'
                 '\n'
                 'Limitations on FASTA input: input FASTA headers must\n'
                 "not contain ':' or '-' ambiguity characters are converted \n"
                 'to N in output'),
                filename=True,
                is_required=True,
                )
            ]
        AbstractCommandline.__init__(self, cmd, **kwargs)
