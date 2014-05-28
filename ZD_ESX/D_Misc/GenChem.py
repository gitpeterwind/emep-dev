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


def split(s, sep=None, maxsplit=-1):
    """Call ``s.split(...)`` and strip whitespace.

    Returns an empty list if *s* is empty.  This differs from ``str.split(...)``
    behaviour which would return ``['']``."""
    if s:
        return [x.strip() for x in s.split(sep, maxsplit)]
    else:
        return []


def ichunk(iterable, size):
    i = iter(iterable)
    while True:
        chunk = list(itertools.islice(i, size))
        if len(chunk) == 0:
            break
        else:
            yield chunk


class DefaultListOrderedDict(collections.OrderedDict):
    """An ordered dictionary that also acts like ``defaultdict(list)``.

    We have a lot of key -> list mappings in this script, it's nice to keep
    insertion order, so this is a special case that does what we want.
    """
    def __getitem__(self, key):
        try:
            return collections.OrderedDict.__getitem__(self, key)
        except KeyError:
            self[key] = new_list = list()
            return new_list


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

    def is_advected(self):
        """Is this an advected species?"""
        return self.type > 0


class SpeciesReader(object):
    FIELDS = ('Spec', 'type', 'formula', 'in_rmm', 'dry', 'wet', 'extinc',
              'cstar', 'DeltaH', None, 'groups', None, 'comment')
    log = IndentingLogger(logging.getLogger('species'))

    def __init__(self):
        self._species = collections.OrderedDict()
        self._groups = DefaultListOrderedDict()

    def read(self, stream):
        """Read species from *stream*."""
        slow = False

        reader = csv.DictReader(stream, self.FIELDS)
        self.log.info('Processing species...')
        self.log.indent()
        new_species = []
        for row in reader:
            # After encountering #SLOW, mark species as Species.SLOW
            if row['Spec'].startswith('#SLOW'):
                self.log.debug('FOUND #SLOW, REST OF SPECIES ARE ADV=3')
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
                           cstar=float(row['cstar']),
                           DeltaH=float(row['DeltaH']))

            if spec.name in self._species:
                raise ValueError("SPEC %s already defined!" % spec.name)

            self.log.debug('SPEC %(name)s', spec)
            self.log.indent()

            #process_groups
            if row['groups'] is not None:
                for g in row['groups'].upper().split(';'):
                    spec.groups.append(g)
                    if row['wet'] is not None:
                        spec.groups.append('WDEP_' + g)
                    if row['dry'] is not None:
                        spec.groups.append('DDEP_' + g)
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
        new_groups = DefaultListOrderedDict()
        for s in new_species:
            for g in self._species[s].groups:
                new_groups[g].append(s)

        # Merge (and log) group changes
        for g in new_groups:
            self._groups[g].extend(new_groups[g])
            self.log.debug('%-11s  =>  %s', g, ', '.join(new_groups[g]))

        self.log.outdent()
        self.log.info('%s groups processed.', len(new_groups))

    def species_list(self):
        """Get a list of Species objects ordered by type."""
        return sorted(self._species.values(), key=lambda spec: spec.type)

    def group_list(self):
        return self._groups.items()


class ReactionsReader(object):
    log = IndentingLogger(logging.getLogger('reactions'))

    def __init__(self, shorthand):
        self.shorthand = shorthand
        self.emisfiles = set()
        self.rate_coefficients = []
        self.photol_specs = set()
        self.emis_specs = set()
        self.bionat_specs = set()

    def read(self, stream):
        """Read reactions from *stream*."""
        self.log.info('Processing reactions...')
        self.log.indent()
        for linenum, line in enumerate(stream, 1):
            # Strip whitespace, convert to uppercase
            line = line.strip().upper()
            # Skip empty lines and comments
            if not line or line.startswith('*'):
                continue
            # Strip comment from end of line
            line = line.partition(';')[0].strip()

            self.log.info('Line %3d: %s', linenum, line)
            self.log.indent()

            if line.startswith('EMISFILES:'):
                self._read_emisfiles(split(line, ':', 1)[1].lower())
            else:
                self._read_reaction(line)

            self.log.outdent()

        self.log.outdent()

    def _read_emisfiles(self, emisfiles):
        """Handle the "foo,bar,baz" from an "emisfiles:foo,bar,baz" line."""
        # Split comma-separated files and add them to the set of emisfiles
        emisfiles = split(emisfiles, ',')
        self.emisfiles.update(emisfiles)
        self.log.info('EMISFILES added %s, now: %s', emisfiles, self.emisfiles)

    def _read_reaction(self, reaction):
        """Handle a "<rate> <LHS> = <RHS>" definition of a reaction."""
        assert '=' in reaction
        # Split up "<rate> <reactants> = <products>"
        rate, terms = reaction.split(None, 1)
        lhs, rhs = split(terms, '=', 1)
        # Parse terms
        lhs = [self._read_reaction_term(_) for _ in split(lhs, '+')]
        rhs = [self._read_reaction_term(_) for _ in split(rhs, '+')]
        # Make sure we don't have [] terms on the RHS
        assert all(type != 'catalyst' for type, _, _ in rhs), '[] terms not allowed in RHS'

        self.log.debug('rate: %s,  LHS: %s,  RHS: %s', rate, lhs, rhs)

        rate = self._read_reaction_rate(rate)

        self.log.debug('processed rate: %s', rate)

    def _read_reaction_term(self, term):
        """reactant or product -> (type, factor, species)"""
        # Split "<factor> <species>" pair
        factor, _, species = term.rpartition(' ')
        # Default factor if not supplied
        factor = factor or '1.0'

        # Figure out type of term
        if species.startswith('[') and species.endswith(']'):
            type = 'catalyst'
            species = species[1:-1]
        elif species.startswith('{') and species.endswith('}'):
            type = 'tracer'
            species = species[1:-1]
        else:
            type = None

        return (type, factor, species)

    def _read_reaction_rate(self, rate):
        """Process a rate, defining new coefficients if necessary, and return
        a usable Fortran expression to refer to the rate."""
        expanded_rate = self.shorthand.expand(rate)
        if expanded_rate != rate:
            self.log.debug('expanded rate: %s', expanded_rate)
        rate = expanded_rate

        # TODO: why do we care?
        if 'TROE' in rate:
            self.log.debug('is TROE')

        # Clean up the rate a bit
        # Lowercase EXP() (TODO: do we need to?)
        #rate = re.sub(r'\bEXP\b', 'exp', rate)
        # Normalise e-notation
        rate = re.sub(r'([\d\.])[EdD]([\+\-]?\d)', r'\1e\2', rate)
        # Replace 1. -> 1.0
        rate = re.sub(r'\.(?=\D)', r'.0', rate)
        self.log.debug('cleaned rate: %s', rate)

        # Resolve rates to Fortran expressions.  In genchem.pl, appears to
        # define rct (the rate coefficient expression), rcttext (a reference to
        # all values for a coefficient in rct) and k/rate_label (the expression
        # to return from define_rates and use in building the overall rate).
        if is_numeric(rate):
            # Use numeric rate as-is
            return rate
        elif 'RCEMIS' in rate or 'RCBIO' in rate:
            # Turn rcemis:foo into rcemis(foo,k)
            return re.sub(r'(RCEMIS|RCBIO):(\w+)', self._get_emis_rate, rate)
        elif 'DJ' in rate:
            # Turn DJ(foo) into rcphot(foo,k)
            return re.sub(r'DJ\((\w+)\)', self._get_Jrate, rate)
        elif 'AQRCK' in rate:
            # Use rate as-is
            return rate
        elif rate.startswith('_FUNC_'):
            # Strip _func_ and then get/create rate coefficient as below
            # TODO: remove _func_ instances from input files
            rate = rate[6:]
            rct = self._get_rate_coefficient(rate)
            return rct['rate_label']
        else:
            # Get/create rate coefficient and return expression referring to it
            rct = self._get_rate_coefficient(rate)
            return rct['rate_label']

    def _get_emis_rate(self, match):
        """Get expression for emission rate.

        *match* is a :class:`re.MatchObject` which should have the type (RCEMIS
        or RCBIO) in ``match.group(1)`` and emission species in
        ``match.group(2)``.  Use with :func:`re.sub`.
        """
        type, spec = match.group(1, 2)
        self.log.debug('PROCESS EMIS: %s : %s', type, spec)
        self.log.indent()

        if spec in self.emis_specs:
            self.log.debug('RCEMIS FOUND: %s, rate = 0', spec)
            rate = '0'
        else:
            self.log.debug('RCEMIS NEW: %s', spec)
            self.emis_specs.add(spec)
            rate = 'rcemis({},k)'.format(spec)

        if type == 'RCBIO':
            if spec in self.bionat_specs:
                self.log.debug('BIONAT FOUND: %s', spec)
            else:
                self.log.debug('BIONAT NEW: %s', spec)
                self.bionat_specs.add(spec)

        self.log.outdent()
        return rate

    def _get_Jrate(self, match):
        """Get expression for photolysis rate.

        *match* is a :class:`re.MatchObject` which should have the photolysis
        species in ``match.group(1)``.  Use with :func:`re.sub`.
        """
        id = match.group(1)
        if id in self.photol_specs:
            self.log.debug('PHOTOL FOUND: %s', id)
        else:
            self.log.debug('PHOTOL NEW: %s', id)
            self.photol_specs.add(id)
        return 'rcphot({},k)'.format(id)

    def _get_rate_coefficient(self, rate):
        """Get (or create, if necessary) a rate coefficient from *rate*."""
        # Bail out early if this rate already has a rate coefficient defined
        for rct in self.rate_coefficients:
            if rct['rct'] == rate:
                self.log.debug('Found rate coefficient for %r: %r', rate, rct)
                return rct

        # Otherwise, define a new rate coefficient
        # (Most of these names are based on those used in GenChem.pl)
        # TODO: handle this whole DIMFLAG/ESXFLAG crazyness from GenChem.pl
        index = len(self.rate_coefficients) + 1
        rct = {
            'index': index,
            'rct': rate,
            'rcttext': 'rct(%d,:)' % index,
            'rate_label': 'rct(%d,k)' % index,
        }
        self.log.info('NEW RCT: %s', rct)
        self.rate_coefficients.append(rct)
        return rct

class CodeGenerator(object):
    def write_module_header(self, stream, module, use=None):
        """Write module header to *stream*, leaving *stream* indented."""
        stream.write(dedent("""\
        ! Generated by GenChem.py - DO NOT EDIT
        """))
        stream.write('module {}\n\n'.format(module))
        stream.indent()
        if use is not None:
            for line in use:
                stream.write(line + '\n')
            stream.write('\n')
        stream.write('implicit none\n')
        # Duplicate behaviour of GenChem.pl (TODO: can we *alway* add this?)
        if use is not None:
            stream.write('private\n')

    def write_contains(self, stream):
        stream.outdent()
        stream.write('\ncontains\n')
        stream.indent()

    def write_module_footer(self, stream, module):
        stream.outdent()
        stream.write('\nend module {}\n'.format(module))


class SpeciesWriter(CodeGenerator):
    log = IndentingLogger(logging.getLogger('write_species'))

    # Groups of species indices to write out
    INDEX_GROUPS = [
        {
            'filter': None,
            'tag': 'TOT',
            'desc': 'All reacting species',
            'ixlab': '',
        },
        {
            'filter': lambda s: s.is_advected(),
            'tag': 'ADV',
            'desc': 'Advected species',
            'ixlab': 'IXADV_',
        },
        {
            'filter': lambda s: s.type == Species.SHORT_LIVED,
            'tag': 'SHL',
            'desc': 'Short-lived (non-advected) species',
            'ixlab': 'IXSHL_',
        },
        {
            'filter': lambda s: s.type == Species.SEMIVOL,
            'tag': 'SEMIVOL',
            'desc': 'Semi-volatile organic aerosols',
            'ixlab': 'IXSOA_',
        },
    ]

    INDICES_HEADER = dedent("""
    !+ Defines indices and NSPEC for {tag} : {desc}

    integer, public, parameter :: NSPEC_{tag}={count}
    """)
    INDICES_HEADER_WITH_EXTENT = INDICES_HEADER.rstrip() + dedent("""
    integer, public, parameter :: FIRST_{tag}={first}, &
                                   LAST_{tag}={last}
    """)

    DECLS = dedent("""
    !/--   Characteristics of species:
    !/--   Number, name, molwt, carbon num, nmhc (1) or not(0)

    public :: define_chemicals    ! Sets names, molwts, carbon num, advec, nmhc

    type, public :: Chemical
         character(len=20) :: name
         real              :: molwt
         integer           :: nmhc      ! nmhc (1) or not(0)
         integer           :: carbons   ! Carbon-number
         real              :: nitrogens ! Nitrogen-number
         integer           :: sulphurs  ! Sulphur-number
         real              :: CiStar     ! VBS param
         real              :: DeltaH    ! VBS param
    endtype Chemical
    type(Chemical), public, dimension(NSPEC_TOT), target :: species

    ! Pointers to parts of species (e.g. short-lived, advected)
    """)

    DEFINE_HEADER = 'subroutine define_chemicals()\n'

    DEFINE_FOOTER = 'end subroutine define_chemicals\n'

    def write(self, stream, all_species):
        self.log.info('Writing species')
        self.log.indent()

        # Wrap stream in indenting writer
        stream = IndentingStreamWriter(stream)

        self.write_module_header(stream, 'ChemSpecs')

        # Write each defined group of indices
        for info in self.INDEX_GROUPS:
            if info['filter'] is None:
                species = all_species
                stream.write(self.INDICES_HEADER.format(count=len(species), **info))
            else:
                species, first, last = self._find_extent(all_species, info['filter'])
                stream.write(self.INDICES_HEADER_WITH_EXTENT.format(
                    count=len(species), first=first, last=last, **info))

            self.log.info('PROCESS %s NSPEC %s', info['tag'], len(species))
            self._write_indices(stream, species, info['ixlab'])

        # Compatibility with existing code
        # TODO: remove this
        stream.write(dedent("""\
        ! Compatibility with existing code
        integer, public, parameter :: NAEROSOL=NSPEC_SEMIVOL
        """))

        stream.write(self.DECLS)

        for info in self.INDEX_GROUPS:
            if info['filter'] is not None:
                stream.write('type(Chemical), public, dimension(:), pointer :: species_{tag}=>null()\n'
                             .format(tag=info['tag'].lower()))

        self.write_contains(stream)

        stream.write(self.DEFINE_HEADER)
        stream.indent()

        stream.write(dedent("""\
        !+
        ! Pointers to parts of species (e.g. short-lived, advected), only assigned if
        ! non-empty.
        !
        """))
        for info in self.INDEX_GROUPS:
            if info['filter'] is not None:
                stream.write(dedent("""\
                if (NSPEC_{tag_uc} > 0) then
                  species_{tag_lc} => species(FIRST_{tag_uc}:LAST_{tag_uc})
                end if
                """).format(tag_uc=info['tag'], tag_lc=info['tag'].lower()))

        stream.write(dedent("""
        !+
        ! Assigns names, mol wts, carbon numbers, advec,  nmhc to user-defined Chemical
        ! array, using indices from total list of species (advected + short-lived).
        !                                                        MW  NM   C    N   S       C*      dH
        """))
        for spec in all_species:
            self._write_spec(stream, spec)

        stream.outdent()
        stream.write(self.DEFINE_FOOTER)

        self.write_module_footer(stream, 'ChemSpecs')

        self.log.outdent()

    def _write_indices(self, stream, species, prefix):
        INDEX_HEADER = '\ninteger, public, parameter :: &\n  '
        INDEX = '{prefix}{spec.name:<12}={i:4d}'
        indices = (INDEX.format(prefix=prefix, i=i, spec=spec) for i, spec in enumerate(species, 1))
        for chunk in ichunk(indices, 10):
            stream.write(INDEX_HEADER)
            stream.write('  &\n  , '.join(chunk) + '\n')

    def _write_spec(self, stream, spec):
        SPEC_FORMAT = ('species({s.name:<12}) = Chemical("{s.name:<12}",{s.molwt:9.4f},'
                       '{nmhc:3d},{s.counts[C]:3d},{s.counts[N]:4d},{s.counts[S]:3d},'
                       '{s.cstar:8.4f},{s.DeltaH:7.1f} )\n')
        stream.write(SPEC_FORMAT.format(s=spec, nmhc=(1 if spec.NMHC else 0)))

    def _find_extent(self, all_species, predicate):
        """Find ``(species, first_ix, last_ix)`` for species matching *predicate*.

        Assumes *all_species* is ordered in such a way that everything matching
        *predicate* is in one contiguous group.  Returns ``([], -999, -999)``
        if nothing is matched.
        """
        offset = 0
        for check, species in itertools.groupby(all_species, predicate):
            if check:
                species = list(species)
                return (species, offset + 1, offset + len(species))
            else:
                offset += len(list(species))
        else:
            return ([], -999, -999)


class GroupsWriter(CodeGenerator):
    log = IndentingLogger(logging.getLogger('write_groups'))

    DECLS = dedent("""
    ! Assignment of groups from GenIn.species
    public :: Init_ChemGroups

    type(typ_sp), dimension({ngroups}), public, save :: chemgroups
    """)

    DECLARE_GROUP = dedent("""
    integer, public, target, save, dimension ({size}) :: &
      {group}_GROUP = (/ {specs} /)
    """)

    DEFINE_GROUP = dedent("""
    chemgroups({i})%name="{group}"
    chemgroups({i})%ptr=>{group}_GROUP
    """)

    def write(self, stream, groups):
        self.log.info('Writing groups')
        self.log.indent()

        # Wrap stream in indenting writer
        stream = IndentingStreamWriter(stream)

        self.write_module_header(stream, 'ChemGroups',
                                 ['use ChemSpecs        ! => species indices',
                                  'use OwnDataTypes     ! => typ_sp'])

        stream.write(self.DECLS.format(ngroups=len(groups)))

        for g, specs in groups:
            self.log.info('PROCESS %-12s N=%2d  =>  %s', g, len(specs), ', '.join(specs))
            stream.write(self.DECLARE_GROUP.format(size=len(specs), group=g,
                                                   specs=','.join(specs)))

        # TODO: RO2 pool species
        stream.write(dedent("""
        ! ------- RO2 Pool     species ------------------
        integer, public, parameter :: SIZE_RO2_POOL = 1
        integer, public, parameter, dimension(1) :: &
          RO2_POOL = (/ -99 /)
        """))

        self.write_contains(stream)

        stream.write('\nsubroutine Init_ChemGroups()\n')
        stream.indent()
        for i, (g, specs) in enumerate(groups, 1):
            stream.write(self.DEFINE_GROUP.format(i=i, group=g))
        stream.outdent()
        stream.write('\nend subroutine Init_ChemGroups\n')

        self.write_module_footer(stream, 'ChemGroups')

        self.log.outdent()


class PrettyStreamHandler(logging.StreamHandler):
    """A :class:`logging.StreamHandler` that wraps log messages with
    severity-dependent ANSI colours."""
    #: Mapping from logging levels to ANSI colours.
    COLOURS = {
        logging.DEBUG: '\033[36m',      # Cyan
        logging.WARNING: '\033[33m',    # Yellow foreground
        logging.ERROR: '\033[31m',      # Red foreground
        logging.CRITICAL: '\033[31;7m'  # Red foreground, inverted
    }
    #: ANSI code for resetting the terminal to default colour.
    COLOUR_END = '\033[0m'

    def format(self, record):
        """Call :meth:`logging.StreamHandler.format`, and apply a colour to the
        message if output stream is a TTY."""
        msg = super(PrettyStreamHandler, self).format(record)
        if self.stream.isatty():
            colour = self.COLOURS.get(record.levelno, '')
            return colour + msg + self.COLOUR_END
        else:
            return msg


if __name__ == '__main__':
    # Send logging output to stderr
    stream_handler = PrettyStreamHandler()
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

    species_reader = SpeciesReader()
    with open('GenIn.species', 'r') as f:
        species_reader.read(f)

    shorthand = ShorthandMap(open('GenIn.shorthand', 'r'))
    reactions_reader = ReactionsReader(shorthand)
    with open('GenIn.reactions', 'r') as f:
        reactions_reader.read(f)

    species_writer = SpeciesWriter()
    species_writer.write(open('CM_ChemSpecs.f90', 'w'), species_reader.species_list())

    groups_writer = GroupsWriter()
    groups_writer.write(open('CM_ChemGroups.f90', 'w'), species_reader.group_list())
