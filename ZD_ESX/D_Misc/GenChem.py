#!/usr/bin/env python
"""
GenChem.py - a Python script to read in equations and output production and
loss terms for Fortran programs.
"""

import logging
import csv
import re
import collections
import itertools
from textwrap import dedent
import sys


def indent(data, amount=1, string='  '):
    return ''.join(string * amount + line for line in data.splitlines(True))


class IndentingLogger(logging.LoggerAdapter):
    """Stateful logger adapter that indents messages.

    Provides :meth:`indent` and :meth:`outdent` to increase and decrease the
    indent level.  All messages have the indentation prepended to them when
    they pass through the adapter.

    >>> log = IndentingLogger(logging.getLogger('foo'))
    >>> log.debug('hello world')
    >>> log.indent()
    >>> log.debug('I am indented')
    >>> log.outdent()
    >>> log.debug('and now I am not')
    """
    def __init__(self, log):
        super(IndentingLogger, self).__init__(log, None)
        self._indent_level = 0

    def indent(self):
        self._indent_level += 1

    def outdent(self):
        self._indent_level = max(0, self._indent_level - 1)

    def process(self, msg, kwargs):
        return (indent(msg, self._indent_level), kwargs)


class IndentingStreamWriter(object):
    """Stateful stream wrapper that indents written strings.

    Provides :meth:`indent` and :meth:`outdent` to increase and decrease the
    indent level.  All messages have the indentation prepended to them when
    they pass through :meth:`write` unless ``indent=False`` is present.
    """
    def __init__(self, stream):
        self._stream = stream
        self._indent_level = 0

    def indent(self):
        self._indent_level += 1

    def outdent(self):
        self._indent_level = max(0, self._indent_level - 1)

    def write(self, data, skip_indent=False):
        if not skip_indent:
            data = indent(data, self._indent_level)
        self._stream.write(data)


# A horrible horrible hack to work around http://bugs.python.org/issue21172
old_LogRecord_init = logging.LogRecord.__init__
# This argument list is only correct for Python 2.7
def new_LogRecord_init(self, name, level, pathname, lineno,
                       msg, args, exc_info, func=None):
    # Initialise LogRecord as normal
    old_LogRecord_init(self, name, level, pathname, lineno,
                       msg, args, exc_info, func)
    # Perform same check as the original constructor, but replace
    # isinstance(args[0], dict) check with hasattr(args[0], '__getitem__'), to
    # match the expectations of the % operator
    if args and len(args) == 1 and hasattr(args[0], '__getitem__') and args[0]:
        # Don't re-do the special case if it succeeded the first time
        if self.args is not args[0]:
            self.args = args[0]
logging.LogRecord.__init__ = new_LogRecord_init


def is_numeric(n):
    """Check if *n* is numeric.

    This replaces ``is_integer`` and ``is_float`` from ``GenChem.pl``.  Valid
    floats are a superset of valid integers, so just uses ``float()`` to check.
    """
    try:
        x = float(n)
        return True
    except (TypeError, ValueError):
        return False


class ShorthandMap(object):
    """Read shorthands from *stream*.

    Given a file with a format like::

        * Ignore this line
        XT           temp(iq)
        H2O          H2O(iq)        ! A comment
        FH2O         (1.+1.4e-21*h2o*exp(2200./XT))

    creates a mapping from shorthand to expansions, which can be applied to
    other equations with the :meth:`expand` method.  Shorthand that appears in
    other expansions is expanded during reading, so in the above example::

        self.mapping['FH2O'] = '(1.+1.4e-21*H2O(IQ)*EXP(2200./TEMP(IQ)))'
    """
    log = IndentingLogger(logging.getLogger('shorthand'))

    def __init__(self, stream):
        self.mapping = collections.OrderedDict()

        self.log.info('Processing shorthands...')
        self.log.indent()
        for line in stream:
            line = line.strip().upper()
            # Skip empty lines and comments
            if line and not line.startswith('*'):
                # Split into pattern, expansion, comments
                parts = line.split(None, 3)

                pattern = parts[0]
                # Expand any shorthand that appears in the expansion
                expanded = self.expand(parts[1])

                self.mapping[pattern] = expanded
                self.log.debug('%-12s  =>  %s', pattern, expanded)

        self.log.outdent()
        self.log.info('%s shorthands processed.', len(self.mapping))

    def expand(self, eqn):
        """Expand shorthand in *eqn*."""
        if len(self.mapping) == 0:
            return eqn

        pattern = r'\b' + '|'.join(self.mapping) + r'\b'

        def replace(matchobj):
            return self.mapping[matchobj.group(0)]

        return re.sub(pattern, replace, eqn)


class Species(object):
    # Species type constants
    SHORT_LIVED = 0
    ADVECTED = 1
    SEMIVOL = 2
    SLOW = 3

    # Atomic weights of atoms we care about
    ATOMS = {
        'C': 12,
        'H': 1,
        'N': 14,
        'O': 16,
        'S': 32,
    }

    def __init__(self, **kwargs):
        self.name = None
        self.type = None
        self.formula = None
        self.molwt = None
        self.groups = []
        self.extinc = None
        self.cstar = None
        self.DeltaH = None
        self.counts = collections.Counter()
        self.NMHC = False

        for k, v in kwargs.iteritems():
            if not hasattr(self, k):
                raise KeyError(k)
            setattr(self, k, v)

    def __getitem__(self, key):
        """Redirect dict-like access to attribute access.

        Particularly useful for the old-style string formatting still used by
        the logging module, e.g. ``"%(name)s" % spec``.
        """
        return getattr(self, key)

    def process_formula(self):
        """Calculate molwt, NMHC flag and atom counts from formula."""
        formula = self.formula or ''

        # Count atoms in the formula
        self.counts = collections.Counter()
        for atom, n in re.findall('([A-Z][a-z]?)(\d*)', formula):
            n = int(n) if n else 1
            # Only store counts for atoms we care about
            if atom in self.ATOMS:
                self.counts.update({atom: n})

        # Sum the molecular weight for the formula
        self.molwt = float(sum(self.ATOMS[atom] * n for atom, n in self.counts.iteritems()))

        # Test for non-methane hydrocarbon (NMHC) - contains only C and H atoms
        # but excluding CH4
        self.NMHC = set(self.counts) == {'C', 'H'} and self.counts['C'] >= 2


class SpeciesReader(object):
    FIELDS = ('Spec', 'type', 'formula', 'in_rmm', 'dry', 'wet', 'extinc',
              'cstar', 'DeltaH', None, 'groups', None, 'comment')
    log = IndentingLogger(logging.getLogger('species'))

    def __init__(self):
        self._species = collections.OrderedDict()
        self._groups = collections.defaultdict(list)

    def read(self, stream):
        """Read species from *stream*."""
        slow = False

        reader = csv.DictReader(stream, self.FIELDS)
        self.log.info('Processing species...')
        self.log.indent()
        new_species = []
        for row in reader:
            # After encountering #SLOW, mark species as Species.SLOW
            if row['Spec'] == '#SLOW':
                slow = True
                continue

            # Skip empty/comment rows
            if not row['Spec'] or row['Spec'] == 'Spec' or row['Spec'][0] in {'*', '#'}:
                continue

            # Strip out ignored field(s)
            del row[None]
            # Replace "blank" fields with None
            for k, v in row.iteritems():
                if v == 'xx' or v == '':
                    row[k] = None
            # Store 'slow' value
            row['slow'] = slow

            spec = Species(name=row['Spec'],
                           type=Species.SLOW if slow else int(row['type']),  # TODO: stop using #SLOW
                           formula=row['formula'],
                           extinc=None if row['extinc'] == '0' else row['extinc'],
                           cstar=row['cstar'],
                           DeltaH=row['DeltaH'])

            if spec.name in self._species:
                raise ValueError("SPEC %s already defined!" % spec.name)

            self.log.debug('SPEC %(name)s', spec)
            self.log.indent()

            #process_groups
            if row['groups'] is not None:
                groups = row['groups'].upper().split(';')
                wet_groups = [] if row['wet'] is None else ['WDEP_' + g for g in groups]
                dry_groups = [] if row['dry'] is None else ['DDEP_' + g for g in groups]
                spec.groups = groups + wet_groups + dry_groups
                self.log.debug('In groups: ' + ', '.join(spec.groups))

            # Get molecular weight, NMHC flag and atom counts from formula
            spec.process_formula()
            self.log.debug('process_formula: %s  =>  MOLWT=%s, NHMC=%s, COUNTS=%r',
                           spec.formula, spec.molwt, spec.NMHC, dict(spec.counts))

            # If molecular weight has been specified, use that instead
            if row['in_rmm'] is not None:
                spec.molwt = float(row['in_rmm'])
                self.log.debug('INPUT MOLWT: %(molwt)s', spec)

            self._species[spec.name] = spec
            new_species.append(spec.name)
            self.log.outdent()

        self.log.outdent()
        self.log.info('%s species processed.', len(new_species))

        self.log.info('Processing groups...')
        self.log.indent()

        # Collect new group memberships from species
        new_groups = collections.defaultdict(list)
        for s in new_species:
            for g in self._species[s].groups:
                new_groups[g].append(s)

        # Merge (and log) group changes
        for g in sorted(new_groups):
            self._groups[g].extend(new_groups[g])
            self.log.debug('%-11s  =>  %s', g, ', '.join(new_groups[g]))

        self.log.outdent()
        self.log.info('%s groups processed.', len(new_groups))

    def species_list(self):
        """Get a list of Species objects ordered by type."""
        return sorted(self._species.values(), key=lambda spec: spec.type)


class CodeGenerator(object):
    MODULE = dedent("""\
    module {module}
    """)



class SpeciesWriter(CodeGenerator):
    log = IndentingLogger(logging.getLogger('write_species'))

    TYPES = collections.OrderedDict([
        (None, {
            'txt': 'tot',
            'desc': 'All reacting species',
            'nspec': 'NSPEC_TOT',
            'ixlab': '',
        }),
        (Species.ADVECTED, {
            'txt': 'adv',
            'desc': 'Advected species',
            'nspec': 'NSPEC_ADV',
            'ixlab': 'IXADV_',
        }),
        (Species.SHORT_LIVED, {
            'txt': 'shl',
            'desc': 'Short-lived (non-advected) species',
            'nspec': 'NSPEC_SHL',
            'ixlab': 'IXSHL_',
        }),
    ])

    INDICES_HEADER = dedent("""
    !+ Defines indices and NSPEC for {txt} : {desc}

    integer, public, parameter :: {nspec}={count}
    """)

    AEROSOL_BLOCK = dedent("""
    integer, public, parameter :: &
      NAEROSOL={count:<5},      & ! Number of aerosol species
      FIRST_SEMIVOL={first:<5}, & ! First aerosol species
      LAST_SEMIVOL={last:<5}     ! Last aerosol species
    """)

    def write(self, stream, all_species):
        self.log.info('Writing species')
        self.log.indent()

        # Wrap stream in indenting writer
        stream = IndentingStreamWriter(stream)

        # Group species into (type, species) pairs
        groups = [(t, list(g)) for t, g in itertools.groupby(all_species, lambda s: s.type)]

        # Find first and last aerosol species
        offset = 0
        for type, species in groups:
            if type == Species.SEMIVOL:
                count = len(species)
                first = offset + 1
                last = offset + count
                break
            else:
                offset += len(species)
        else:
            # If we didn't find aerosol type (didn't hit "break") then
            # set default values
            count = 0
            first = -999
            last = -999

        # Write information for all species
        info = self.TYPES[None]
        self.log.info('PROCESS %s NSPEC %s', info['txt'], len(all_species))
        stream.write(self.INDICES_HEADER.format(count=len(all_species), **info))
        stream.write(self.AEROSOL_BLOCK.format(count=count, first=first, last=last))
        stream.write(self._gen_indices(all_species, info['ixlab']))

        # Write information for individual types
        for type, species in groups:
            # Skip types we don't have information for
            if type not in self.TYPES:
                continue
            info = self.TYPES[type]
            self.log.info('PROCESS %s NSPEC %s', info['txt'], len(species))
            species = list(species)
            stream.write(self.INDICES_HEADER.format(count=len(species), **info))
            stream.write(self._gen_indices(species, info['ixlab']))

        self.log.outdent()

    def _gen_indices(self, species, prefix=''):
        INDEX = '{prefix}{spec.name:<12}={i:4d}'
        INDEX_GROUP = dedent("""
        integer, public, parameter :: &
          {}
        """)
        indices = [INDEX.format(prefix=prefix, i=i, spec=spec) for i, spec in enumerate(species, 1)]
        chunks = (indices[i:i+10] for i in xrange(0, len(indices), 10))
        groups = (INDEX_GROUP.format(self.indent('   &\n, '.join(chunk))) for chunk in chunks)
        return ''.join(groups)


class ReactionsReader(object):
    def __init__(self, reactions_file):
        pass


if __name__ == '__main__':
    # Send logging output to stderr
    stream_handler = logging.StreamHandler()
    #stream_handler.setLevel(logging.INFO)
    stream_handler.setFormatter(logging.Formatter('[%(levelname).1s] %(message)s'))
    # Send logging output to Log.GenOut
    file_handler = logging.FileHandler('Log.GenOut', 'w')
    file_handler.setFormatter(logging.Formatter('%(message)s'))
    # Attach log handlers
    rootlogger = logging.getLogger('')
    rootlogger.setLevel(logging.DEBUG)
    rootlogger.addHandler(stream_handler)
    rootlogger.addHandler(file_handler)

    shorthand = ShorthandMap(open('GenIn.shorthand', 'r'))
    species_reader = SpeciesReader()
    species_reader.read(open('GenIn.species', 'r'))

    species_writer = SpeciesWriter()
    species_writer.write(open('CM_ChemSpecs.f90', 'w'), species_reader.species_list())
