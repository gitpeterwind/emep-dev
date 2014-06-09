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

from ordered_set import OrderedSet


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


def element_remainder_pairs(elements):
    """'ABC' -> [('A', 'BC'), ('B', 'AC'), ('C', 'AB')]"""
    for i, x in enumerate(elements):
        yield (x, elements[:i] + elements[i+1:])


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


LOG = IndentingLogger(logging.getLogger(''))


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
    def __init__(self, stream):
        self.mapping = collections.OrderedDict()

        LOG.info('Processing shorthands...')
        LOG.indent()
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
                LOG.debug('%-12s  =>  %s', pattern, expanded)

        LOG.outdent()
        LOG.info('%s shorthands processed.', len(self.mapping))

    def expand(self, eqn):
        """Expand shorthand in *eqn*."""
        if len(self.mapping) == 0:
            return eqn

        pattern = r'\b(' + '|'.join(self.mapping) + r')\b'

        def replace(matchobj):
            return self.mapping[matchobj.group(0)]

        return re.sub(pattern, replace, eqn)


#: Atomic weights of atoms we care about
ATOMS = {
    'C': 12,
    'H': 1,
    'N': 14,
    'O': 16,
    'S': 32,
}


def count_atoms(formula):
    """formula -> (counts, molwt)

    Process a formula to find out how many of each atom in *ATOMS* it contains
    and the molecular weight from only those atoms.
    """
    counts = collections.Counter()
    # Count atoms in the formula
    for atom, n in re.findall('([A-Z][a-z]?)(\d*)', formula):
        n = int(n) if n else 1
        if atom in ATOMS:
            counts.update({atom: n})
    # Sum the molecular weight
    molwt = float(sum(ATOMS[atom] * n for atom, n in counts.iteritems()))
    return (counts, molwt)


class ChemicalScheme(object):
    def __init__(self):
        self.species = collections.OrderedDict()
        self.groups = DefaultListOrderedDict()
        self.emis_files = OrderedSet()
        self.reactions = []

    def add_species(self, spec):
        """Add *spec* to species, updating groups as necessary."""
        if spec.name in self.species:
            raise ValueError("SPEC %s already defined!" % spec.name)

        self.species[spec.name] = spec

        for group in spec.groups:
            self.groups[group].append(spec.name)

    def add_emis_files(self, files):
        """Add *files* to emis_files."""
        self.emis_files |= files
        LOG.info('EMISFILES added %s, now: %s', files, self.emis_files)

    def add_reaction(self, reaction):
        """Add *reaction* to reactions."""
        self.reactions.append(reaction)

    def get_species_list(self):
        """Get a list of species sorted by type."""
        return sorted(self.species.values(), key=lambda spec: spec.type)

    def get_group_list(self):
        """Get a list of groups."""
        return self.groups.items()


class Species(object):
    # Species type constants
    SHORT_LIVED = 0
    ADVECTED = 1
    SEMIVOL = 2
    SLOW = 3

    def __init__(self, name, type, formula, extinc, cstar, DeltaH, molwt=None, groups=None):
        # Arguments simply copied to attributes
        self.name = name
        self.type = type
        self.formula = formula
        self.extinc = extinc
        self.cstar = float(cstar)
        self.DeltaH = float(DeltaH)

        # Calculate atom counts and molecular weight from formula
        self.counts, self.molwt = count_atoms(self.formula or '')
        # Is non-methane hydrocarbon?
        self.NMHC = set(self.counts) == {'C', 'H'} and self.counts['C'] >= 2

        LOG.debug('processed formula: %s  =>  MOLWT=%s, NMHC=%s, COUNTS=%r',
                  self.formula, self.molwt, self.NMHC, dict(self.counts))

        # Override molecular weight?
        if molwt is not None:
            self.molwt = float(molwt)
            LOG.debug('override MOLWT: %s', self.molwt)

        # Groups, default to empty list
        self.groups = groups or []

    def __getitem__(self, key):
        """Redirect dict-like access to attribute access.

        Particularly useful for the old-style string formatting still used by
        the logging module, e.g. ``"%(name)s" % spec``.
        """
        return getattr(self, key)

    def is_advected(self):
        """Is this an advected species?"""
        return self.type > 0


class SpeciesReader(object):
    FIELDS = ('Spec', 'type', 'formula', 'in_rmm', 'dry', 'wet', 'extinc',
              'cstar', 'DeltaH', None, 'groups', None, 'comment')
    def __init__(self, scheme):
        self.scheme = scheme

    def read(self, stream):
        """Read species from *stream*."""
        LOG.info('Processing species...')
        LOG.indent()

        reader = csv.DictReader(stream, self.FIELDS)
        for row in reader:
            # The #SLOW line is the old way of doing adv=3, and is inflexible.
            # We don't allow it here.
            if row['Spec'].startswith('#SLOW'):
                raise ValueError('#SLOW not allowed, use adv=3')

            # Skip empty/comment rows
            if not row['Spec'] or row['Spec'] == 'Spec' or row['Spec'][0] in {'*', '#'}:
                continue

            # Strip out ignored field(s)
            del row[None]
            # Replace "blank" fields with None
            for k, v in row.iteritems():
                if v == 'xx' or v == '':
                    row[k] = None

            LOG.debug('SPEC %s', row['Spec'])
            LOG.indent()

            # Process groups
            groups = []
            if row['groups'] is not None:
                for g in row['groups'].upper().split(';'):
                    groups.append(g)
                    if row['wet'] is not None:
                        groups.append('WDEP_' + g)
                    if row['dry'] is not None:
                        groups.append('DDEP_' + g)
                LOG.debug('In groups: ' + ', '.join(groups))

            spec = Species(name=row['Spec'],
                           type=int(row['type']),
                           formula=row['formula'],
                           extinc=None if row['extinc'] == '0' else row['extinc'],
                           cstar=row['cstar'],
                           DeltaH=row['DeltaH'],
                           molwt=row['in_rmm'],
                           groups=groups)

            self.scheme.add_species(spec)

            LOG.outdent()

        LOG.outdent()


Term = collections.namedtuple('Term', ['species', 'factor', 'type'])


class Reaction(object):
    def __init__(self, rate, LHS, RHS):
        self.rate = rate
        # Make sure we don't have factors on the LHS
        assert all(_.factor == '1.0' for _ in LHS), 'factors not allowed in LHS'
        self.LHS = LHS
        # Make sure we don't have [] terms on the RHS
        assert all(_.type != 'catalyst' for _ in RHS), '[] terms not allowed in RHS'
        self.RHS = RHS

    def __repr__(self):
        return 'Reaction(rate=%r, LHS=%r, RHS=%r)' % (self.rate, self.LHS, self.RHS)

    def __str__(self):
        return self.__repr__()

    @property
    def reactants(self):
        """Reaction reactants (from LHS)."""
        return [_ for _ in self.LHS if _.type is None]

    @property
    def catalysts(self):
        """Reaction catalysts (from LHS)."""
        return [_ for _ in self.LHS if _.type == 'catalyst']

    @property
    def products(self):
        """Reaction products (from RHS)."""
        return [_ for _ in self.RHS if _.type is None]

    @property
    def tracers(self):
        """Reaction tracers (from RHS)."""
        return [_ for _ in self.RHS if _.type is not None]

    def get_full_rate(self):
        """Full rate, including rates contributed by LHS catalysts."""
        rate = self.rate[:]
        for term in self.catalysts:
            rate.extend(['*', ('amount', term.species)])
        return rate

    def get_loss_rates(self):
        """Iterate over ``(species, rate)`` pairs for this reaction's losses."""
        base_rate = self.get_full_rate()
        for term, others in element_remainder_pairs(self.reactants):
            # Multiply rate by amounts of all other reactants
            term_rate = []
            for t in others:
                term_rate.extend(['*', ('amount', t.species)])
            yield (term.species, base_rate + term_rate)

    def get_prod_rates(self):
        """Iterate over ``(species, rate)`` pairs for this reaction's production."""
        base_rate = self.get_full_rate()
        # Multiply rate by amounts of all reactants
        for term in self.reactants:
            base_rate.extend(['*', ('amount', term.species)])
        for term in self.products:
            # Include product factor if present
            if term.factor == '1.0':
                term_rate = []
            else:
                term_rate = [term.factor + '*']
            yield (term.species, term_rate + base_rate)


class ReactionsReader(object):
    def __init__(self, scheme, shorthand):
        self.scheme = scheme
        self.shorthand = shorthand

    def read(self, stream):
        """Read reactions from *stream*."""
        LOG.info('Processing reactions...')
        LOG.indent()
        for linenum, line in enumerate(stream, 1):
            # Strip whitespace, convert to uppercase
            line = line.strip().upper()
            # Skip empty lines and comments
            if not line or line.startswith('*'):
                continue
            # Strip comment from end of line
            line = line.partition(';')[0].strip()

            LOG.info('Line %3d: %s', linenum, line)
            LOG.indent()

            if line.startswith('EMISFILES:'):
                files = split(line, ':', 1)[1].lower()
                self.scheme.add_emis_files(split(files, ','))
            else:
                self.scheme.add_reaction(self._read_reaction(line))

            LOG.outdent()

        LOG.outdent()

    def _read_reaction(self, reaction):
        """Turn a "<rate> <LHS> = <RHS>" definition into a Reaction object."""
        assert '=' in reaction
        # Split up "<rate> <reactants> = <products>"
        rate, terms = reaction.split(None, 1)
        lhs, rhs = split(terms, '=', 1)
        # Parse LHS and RHS terms
        lhs = [self._read_reaction_term(_) for _ in split(lhs, '+')]
        rhs = [self._read_reaction_term(_) for _ in split(rhs, '+')]

        # Process rate and create reaction
        rate = self._read_reaction_rate(rate)
        reaction = Reaction(rate, lhs, rhs)
        LOG.debug('%s', reaction)

        # TODO: check that specs exist, where appropriate
        # TODO: check reaction balance
        # TODO: what is check_multipliers()?
        # TODO: RHS tracers?

        LOG.debug('full rate including catalysts: %r', reaction.get_full_rate())

        LOG.debug('loss rates:')
        LOG.indent()
        for spec, rate in reaction.get_loss_rates():
            LOG.debug('%s: %s', spec, rate)
        LOG.outdent()

        LOG.debug('prod rates:')
        LOG.indent()
        for spec, rate in reaction.get_prod_rates():
            LOG.debug('%s: %s', spec, rate)
        LOG.outdent()

        return reaction

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

        return Term(species, factor, type)

    def _read_reaction_rate(self, rate):
        """Process *rate* into a list of rate parts.

        The rate first has shorthand expansions and some basic cleanup applied
        to it.  It's then split by "top-level" operators (operators that don't
        appear inside parentheses) and each "part" is checked for special forms
        that need special handling instead of being included as-is in the rate.
        The special forms are replaced with ``(kind, arg)`` tuples, and
        everything else is left as literal strings.

        For example, ``DJ(IDNO2)*0.222`` becomes ``[('photol', 'IDNO2'), '*0.222']``.
        The responsibility for turning this into a Fortran expression describing
        the rate (e.g. ``rcphot(IDNO2,k)*0.222``) belongs to the code generator.

        Non-trivial rates that do not cause any special handling become
        ``[('coeff', rate)]``, which gets used in the code generator to re-use
        rate calculations.
        """
        # Expand shorthands in the rate
        expanded_rate = self.shorthand.expand(rate)
        if expanded_rate != rate:
            LOG.debug('expanded rate: %s', expanded_rate)
        rate = expanded_rate

        # TODO: why do we care?
        if 'TROE' in rate:
            LOG.debug('is TROE')

        # Clean up the rate a bit
        # Lowercase EXP() (TODO: do we need to?)
        #rate = re.sub(r'\bEXP\b', 'exp', rate)
        # Normalise e-notation
        rate = re.sub(r'([\d\.])[EdD]([\+\-]?\d)', r'\1e\2', rate)
        # Replace 1. -> 1.0
        rate = re.sub(r'\.(?=\D)', r'.0', rate)
        LOG.debug('cleaned rate: %s', rate)

        # Short-circuit some cases where we just want the literal rate
        if is_numeric(rate) or 'AQRCK' in rate:
            return [rate]
        elif rate.startswith('_FUNC_'):
            return [rate[6:]]

        parts = []
        def handle_part(part, op=None):
            # Recognise and replace special forms with (kind, arg) tuples
            if part.startswith('RCEMIS:'):
                # rcemis:foo -> ('emis', 'foo')
                part = ('emis', part[7:])
            elif part.startswith('RCBIO:'):
                # rcbio:foo -> ('emis', 'foo')
                # TODO: remove this, after fixing input files
                LOG.warning('rcbio:spec deprecated, replace with rcemis:spec')
                part = ('emis', part[6:])
            elif part.startswith('DJ('):
                # DJ(foo) -> ('photol', 'foo')
                part = ('photol', part[3:-1])

            # If this is the first part, a special case, or follows a special
            # case, save as a new part, otherwise append it to the previous one.
            # (Keeps it down to the most concise form that still represents the
            # rate properly.)
            if isinstance(part, tuple) or len(parts) == 0 or isinstance(parts[-1], tuple):
                parts.append(part)
            else:
                parts[-1] += part

            # If the part has an associated operator (rather than terminated by
            # the end of the string), add that to the most recent part if it
            # was a string, or start a new part if it was a special case.
            if op is not None:
                if isinstance(parts[-1], tuple):
                    parts.append(op)
                else:
                    parts[-1] += op

        # Process the rate, splitting on operators, but never breaking inside
        # parentheses.  Use the handle_part function above to handle each part.
        OPS = set('+-/*')
        part_start = 0
        paren_depth = 0
        for i, c in enumerate(rate):
            if c == '(':
                paren_depth += 1
            elif c == ')':
                paren_depth -= 1
            elif paren_depth == 0 and c in OPS:
                # Found an operator outside of paretheses, end the current part
                # and process it
                handle_part(rate[part_start:i], rate[i])
                # Start a new part
                part_start = i + 1
            else:
                # Proceed to next character, accumulating the current part
                pass
        assert paren_depth == 0, 'mismatched parentheses'
        # Handle final part
        handle_part(rate[part_start:])

        # If we just have just a single literal expression by this point, it
        # should probably be a rate coefficient (see the short-circuit section
        # at the top of this function for cases where this won't happen).  The
        # code generator can build a database of rate coefficients and reuse
        # calculations.
        if parts == [rate]:
            parts = [('coeff', rate)]

        return parts


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

    def __init__(self, scheme):
        self.scheme = scheme

    def write(self, stream):
        LOG.info('Writing species')
        LOG.indent()

        all_species = self.scheme.get_species_list()

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

            LOG.info('PROCESS %s NSPEC %s', info['tag'], len(species))
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

        LOG.outdent()

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

    def __init__(self, scheme):
        self.scheme = scheme

    def write(self, stream):
        LOG.info('Writing groups')
        LOG.indent()

        groups = self.scheme.get_group_list()

        # Wrap stream in indenting writer
        stream = IndentingStreamWriter(stream)

        self.write_module_header(stream, 'ChemGroups',
                                 ['use ChemSpecs        ! => species indices',
                                  'use OwnDataTypes     ! => typ_sp'])

        stream.write(self.DECLS.format(ngroups=len(groups)))

        for g, specs in groups:
            LOG.info('PROCESS %-12s N=%2d  =>  %s', g, len(specs), ', '.join(specs))
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

        LOG.outdent()


class ReactionsWriter(CodeGenerator):
    EMIS_REF = 'rcemis({spec},k)'
    PHOTOL_REF = 'rcphot({spec},k)'
    COEFF_REF = 'rct({i},k)'
    COEFF_ALL = 'rct({i},:)'
    AMOUNT_REF = 'xnew({spec})'

    CHEMEQN_FULL = 'xnew({spec}) = (xold({spec}) + dt2 * P) / (1.0 + dt2 * L)\n'
    CHEMEQN_NOPROD = 'xnew({spec}) = xold({spec}) / (1.0 + dt2 * L)\n'
    CHEMEQN_NOLOSS = 'xnew({spec}) = xold({spec}) + dt2 * P\n'
    CHEMEQN_NONE = '! Nothing to do for {spec}! xnew({spec}) = max(0.0, xold({spec}))\n'
    # Map (has_prod, has_loss) to an equation string
    CHEMEQN_LOOKUP = {
        (True, True): CHEMEQN_FULL,
        (False, True): CHEMEQN_NOPROD,
        (True, False): CHEMEQN_NOLOSS,
        (False, False): CHEMEQN_NONE,
    }

    CHEMRATES_USE = [
        'use Zmet_ml        ! => tinv, h2o, m, Fgas !ESX',
        'use ZchemData      ! => rct',
        'use ChemSpecs      ! => NSPEC_TOT, PINALD, .... for FgasJ08',
        'use DefPhotolysis  ! => IDNO2 etc.',
    ]

    CHEMRATES_DECLS = dedent("""
    public :: setChemRates
    public :: setPhotolUsed

    integer, parameter, public :: NCHEMRATES = {coeff_count}  ! No. coefficients

    ! Photolysis rates
    integer, parameter, public :: NPHOTOLRATES = {photol_count}  ! No. DJ vals used
    integer, save, public, dimension(NPHOTOLRATES) :: photol_used
    """)

    def __init__(self, scheme):
        self.scheme = scheme
        self.coefficients = OrderedSet()

        self.prod = DefaultListOrderedDict()
        self.loss = DefaultListOrderedDict()
        self.emis_specs = OrderedSet()
        self.photol_specs = OrderedSet()
        self.coefficients = OrderedSet()

        LOG.info('Extracting prod/loss from reactions...')
        LOG.indent()

        def process_rate_part(part):
            if isinstance(part, tuple):
                kind, arg = part
                if kind == 'emis':
                    if arg in self.emis_specs:
                        LOG.warning('RCEMIS duplicate: %s (using rate = 0)', arg)
                        return '0'
                    else:
                        self.emis_specs.add(arg)
                        LOG.debug('RCEMIS new: %s', arg)
                        return self.EMIS_REF.format(spec=arg)
                elif kind == 'photol':
                    if arg in self.photol_specs:
                        LOG.debug('PHOTOL found: %s', arg)
                    else:
                        LOG.debug('PHOTOL new: %s', arg)
                        self.photol_specs.add(arg)
                    return self.PHOTOL_REF.format(spec=arg)
                elif kind == 'coeff':
                    if arg not in self.coefficients:
                        LOG.debug('RCT NEW %2d: %s', len(self.coefficients), arg)
                    index = self.coefficients.add(arg) + 1
                    return self.COEFF_REF.format(i=index)
                elif kind == 'amount':
                    return self.AMOUNT_REF.format(spec=arg)
                else:
                    raise ValueError('unhandled rate part', part)
            else:
                return part

        for reaction in self.scheme.reactions:
            LOG.info('PROCESS %r', reaction)
            LOG.indent()
            for spec, rate in reaction.get_prod_rates():
                LOG.info('PROD  %-12s: %s', spec, rate)
                LOG.indent()
                processed_rate = ' '.join(process_rate_part(r) for r in rate)
                self.prod[spec].append(processed_rate)
                LOG.info('RATE: %s', processed_rate)
                LOG.outdent()
            for spec, rate in reaction.get_loss_rates():
                LOG.info('LOSS  %-12s: %s', spec, rate)
                LOG.indent()
                processed_rate = ' '.join(process_rate_part(r) for r in rate)
                self.loss[spec].append(processed_rate)
                LOG.info('RATE: %s', processed_rate)
                LOG.outdent()
            LOG.outdent()

        LOG.outdent()

    def write_prod_loss(self, normal, slow):
        LOG.info('Writing prod/loss files...')
        LOG.indent()

        # Wrap streams in indenting writers
        normal = IndentingStreamWriter(normal)
        slow = IndentingStreamWriter(slow)

        # TODO: handle "RO2POOL"
        for spec in self.scheme.get_species_list():
            prod = self.prod[spec.name]
            loss = self.loss[spec.name]
            LOG.info('SPEC %-12s: nprod=%2d, nloss=%2d', spec.name, len(prod), len(loss))
            # Choose appropriate output stream
            stream = slow if spec.type == Species.SLOW else normal

            stream.write('!-> {spec}\n\n'.format(spec=spec.name))
            stream.indent()

            if prod:
                stream.write('P = ' + '  &\n  + '.join(prod) + '\n')
            else:
                stream.write('! P = 0.0\n')
            stream.write('\n')

            if loss:
                stream.write('L = ' + '  &\n  + '.join(loss) + '\n')
            else:
                stream.write('! L = 0.0\n')
            stream.write('\n')

            eqn = self.CHEMEQN_LOOKUP[(bool(prod), bool(loss))].format(spec=spec.name)
            stream.write(eqn)

            stream.outdent()
            stream.write('\n\n')

        LOG.outdent()

    def write_rates(self, stream):
        LOG.info("Writing rates...")
        LOG.indent()

        stream = IndentingStreamWriter(stream)

        self.write_module_header(stream, 'ChemRates', self.CHEMRATES_USE)
        stream.write(self.CHEMRATES_DECLS.format(coeff_count=len(self.coefficients),
                                                 photol_count=len(self.photol_specs)))
        self.write_contains(stream)

        stream.write('\nsubroutine setPhotolUsed()\n')
        stream.indent()
        stream.write('photol_used = (/ &\n    ' +
                     '  &\n  , '.join(self.photol_specs) +
                     '  &\n/)')
        stream.outdent()
        stream.write('\nend subroutine setPhotolUsed\n')

        stream.write(dedent("""
        subroutine setChemRates(debug_level)
          integer, intent(in) :: debug_level

        """))
        stream.indent()

        for i, rate in enumerate(self.coefficients, 1):
            # Replace XT with temp
            # TODO: shouldn't this be in shorthands or ... the input?
            rate = re.sub(r'\b(XT|xt)\b', 'temp', rate)
            # TODO: wrap equations so the lines aren't too long for the compiler
            stream.write((self.COEFF_ALL + ' = {rate}\n').format(i=i, rate=rate))

        stream.outdent()
        stream.write('\nend subroutine setChemRates\n')

        self.write_module_footer(stream, 'ChemRates')

        LOG.outdent()


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

    # Hook unhandled exception handling into logging
    def excepthook(type, value, tb):
        import traceback
        for part in traceback.format_exception(type, value, tb):
            for line in part.splitlines():
                LOG.critical(line)
    sys.excepthook = excepthook

    scheme = ChemicalScheme()

    species_reader = SpeciesReader(scheme)
    with open('GenIn.species', 'r') as f:
        species_reader.read(f)

    shorthand = ShorthandMap(open('GenIn.shorthand', 'r'))
    reactions_reader = ReactionsReader(scheme, shorthand)
    with open('GenIn.reactions', 'r') as f:
        reactions_reader.read(f)

    species_writer = SpeciesWriter(scheme)
    species_writer.write(open('CM_ChemSpecs.f90', 'w'))

    groups_writer = GroupsWriter(scheme)
    groups_writer.write(open('CM_ChemGroups.f90', 'w'))

    reactions_writer = ReactionsWriter(scheme)
    with open('CM_Reactions1.inc', 'w') as f1, open('CM_Reactions2.inc', 'w') as f2:
        reactions_writer.write_prod_loss(f1, f2)
    with open('CM_ChemRates.f90', 'w') as f:
        reactions_writer.write_rates(f)
