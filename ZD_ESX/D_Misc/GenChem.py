#!/usr/bin/env python
"""
GenChem.py - a Python script to read in equations and output production and
loss terms for Fortran programs.
"""

import logging
import csv
import re
import collections


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
    INDENT_STRING = '  '
    _indent_level = 0

    def __init__(self, log):
        super(IndentingLogger, self).__init__(log, None)

    def indent(self):
        self._indent_level += 1

    def outdent(self):
        self._indent_level = max(0, self._indent_level - 1)

    def process(self, msg, kwargs):
        return ((self.INDENT_STRING * self._indent_level + msg), kwargs)


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

        ATOMS = {'C': 12, 'H': 1, 'N': 14, 'O': 14, 'S': 32}
        self.counts = collections.Counter()

        # Count atoms in the formula
        for atom, n in re.findall('([A-Z][a-z]?)(\d*)', formula):
            n = int(n) if n else 1
            # Only store counts for atoms we care about
            if atom in ATOMS:
                self.counts.update({atom: n})

        # Sum the molecular weight for the formula
        self.molwt = float(sum(ATOMS[atom] * n for atom, n in self.counts.iteritems()))

        # Test for non-methane hydrocarbon (NMHC) - contains only C and H atoms
        # but excluding CH4
        self.NMHC = set(self.counts) == {'C', 'H'} and self.counts['C'] >= 2


class SpeciesReader(object):
    """Read species from *stream*.
    """
    FIELDS = ('Spec', 'type', 'formula', 'in_rmm', 'dry', 'wet', 'extinc',
              'cstar', 'DeltaH', None, 'groups', None, 'comment')
    log = IndentingLogger(logging.getLogger('species'))

    def __init__(self, stream):
        self.species = collections.OrderedDict()

        slow = False

        reader = csv.DictReader(stream, self.FIELDS)
        self.log.info('Processing species...')
        self.log.indent()
        for row in reader:
            # After encountering #SLOW, use self.slow_species
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
                           type=3 if slow else int(row['type']),  # TODO: stop using #SLOW
                           formula=row['formula'],
                           extinc=row['extinc'],
                           cstar=row['cstar'],
                           DeltaH=row['DeltaH'])

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

            self.species[spec.name] = spec
            self.log.outdent()

        self.log.outdent()
        self.log.info('%s species processed.', len(self.species))

        self.log.info('Processing groups...')
        self.log.indent()

        # Build groups of species from the groups declared for each species
        self.groups = collections.defaultdict(list)
        for s in self.species.itervalues():
            for g in s.groups:
                self.groups[g].append(s.name)

        for g in sorted(self.groups):
            self.log.debug('%-11s  =>  %s', g, ', '.join(self.groups[g]))

        self.log.outdent()
        self.log.info('%s groups created.', len(self.groups))


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
    species = SpeciesReader(open('GenIn.species', 'r'))
