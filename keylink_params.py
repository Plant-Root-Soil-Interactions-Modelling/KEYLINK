"""
Schema-driven, validated reading of the KEYLINK parameter files (KL_*.txt).

The file formats are unchanged: whitespace separated numbers, one entry per line (or one
functional group per column), followed by a '# comment'. What is new is that every entry is
identified by NAME instead of by its line number: the schema below lists, for each file, the
expected entries in order with a keyword that must occur in the line's comment, the unit, and the
valid range. Loading therefore fails with an explicit message when a line is missing, reordered,
mislabelled or out of range, instead of silently shifting all following values (the error class
behind several bugs of earlier versions).

    >>> from keylink_params import load_parameter_file
    >>> soil = load_parameter_file('KL_initSoil')
    >>> soil['pH'], soil['PVstruct']

Group-wise files (KL_FaunalParams, KL_densityParams) have a header line with the nine functional
groups, which is checked against keylink_names.GROUP_KEYS (aliases allowed), and return one
array of length 9 per parameter.
"""
import os
import re

import numpy as np

from keylink_names import NGROUPS, NPOOLS, GROUP_KEYS

INF = float('inf')


class ParameterFileError(ValueError):
    """A KL_*.txt file is missing, has the wrong number of lines, a mislabelled line or a value
    outside its valid range."""


class Entry:
    """One parameter of a file: `key` in the returned dict, `keyword` (regex, case-insensitive)
    that must occur in the line's comment, unit, valid (min, max) range, `n` values per line
    (1, or NGROUPS for group-wise files, or 12 for monthly values)."""

    def __init__(self, key, keyword, unit='', lo=-INF, hi=INF, n=1, description=''):
        self.key, self.keyword, self.unit, self.lo, self.hi, self.n, self.description = \
            key, keyword, unit, lo, hi, n, description


# names accepted in the header of group-wise files, per functional group (case-insensitive prefix)
GROUP_ALIASES = {
    'bact': ('bact',), 'fungi': ('fung',), 'myc': ('em', 'myc'), 'bvores': ('bacvor', 'bvores', 'bacterivor'),
    'fvores': ('funvor', 'fvores', 'fungivor'), 'sap': ('detrvor', 'sap', 'detritivor'),
    'eng': ('engin', 'eng'), 'hvores': ('herbvor', 'hvores', 'herbivor'), 'pred': ('pred',)}

SCHEMA = {
    'KL_FaunalParams': dict(groupwise=True, entries=[
        Entry('gmax', r'gmax', 'g C g C-1 day-1', 0, 20, NGROUPS, 'maximal growth rate'),
        Entry('ks', r'\bks\b', 'g C m-3', 1e-9, INF, NGROUPS, 'half saturation constant'),
        Entry('death', r'death', 'g C g C-1 day-1', 0, 1, NGROUPS, 'turnover rate d'),
        Entry('resp', r'resp', 'g C g C-1 day-1', 0, 1, NGROUPS, 'respiration rate at Topt'),
        Entry('faec', r'f(a)?ec', 'fraction', 0, 1, NGROUPS, 'fraction of consumed C to faeces'),
        Entry('cn', r'\bcn\b', '-', 0.5, 1000, NGROUPS, 'C:N of the biomass'),
        Entry('rec', r'\brec\b', 'fraction', 0, 1, NGROUPS, 'recalcitrance of the biomass (unused)'),
        Entry('mcn', r'(pm|m)cn', '-', 0, 5, NGROUPS, 'sensitivity of gmax to the C:N of the food'),
        Entry('mrec', r'(pm|m)rec', '-', 0, 1, NGROUPS, 'sensitivity of gmax to litter recalcitrance'),
        Entry('tmin', r'tmin', 'degC', -50, 60, NGROUPS, 'minimum temperature for growth'),
        Entry('topt', r'topt', 'degC', -50, 60, NGROUPS, 'optimum temperature'),
        Entry('tmax', r'tmax', 'degC', -50, 80, NGROUPS, 'maximum temperature'),
        Entry('q10', r'q10', '-', 1, 10, NGROUPS, 'Q10 of respiration')]),
    'KL_densityParams': dict(groupwise=True, entries=[
        Entry('beta', r'beta', '-', 1, 5, NGROUPS, 'exponent of the density dependent turnover'),
        Entry('bref', r'bref', 'g C m-3', -INF, INF, NGROUPS, 'reference biomass (<= 0: initial biomass)')]),
    'KL_initSoil': dict(groupwise=False, entries=[
        Entry('d', r'depth', 'm', 1e-3, 100, 1, 'depth of the soil layer'),
        Entry('BD', r'bulk', 'kg m-3', 0, 3000, 1, 'bulk density (not used by the equations)'),
        Entry('alpha', r'alpha', 'kPa-1', 0, INF, 1, 'van Genuchten alpha'),
        Entry('n', r'^\s*n\b', '-', 1, 10, 1, 'van Genuchten n'),
        Entry('m', r'^\s*m\b', '-', 0, 1, 1, 'van Genuchten m'),
        Entry('Ksat', r'ksat', 'mm day-1', 0, INF, 1, 'saturated conductivity = maximal matrix infiltration rate'),
        Entry('pH', r'\bph\b', '-', 0, 14, 1, 'soil pH'),
        Entry('litterCN0', r'litter\s*cn', '-', 1, 1000, 1, 'initial C:N of the litter pool'),
        Entry('SOMCN0', r'som\s*cn', '-', 1, 1000, 1, 'initial C:N of SOM'),
        Entry('drainmax', r'drainmax', 'mm day-1', 0, INF, 1, 'infiltration capacity of the layer below'),
        Entry('PV_inac', r'inac', 'l m-3', 1e-9, 1000, 1, 'textural volume of inaccessible pores'),
        Entry('PV_bact', r'bact', 'l m-3', 1e-9, 1000, 1, 'textural volume of bacterial pores'),
        Entry('PV_micro', r'micro', 'l m-3', 1e-9, 1000, 1, 'textural volume of micropores'),
        Entry('PV_meso', r'meso', 'l m-3', 1e-9, 1000, 1, 'textural volume of mesopores'),
        Entry('PV_macro', r'macro', 'l m-3', 1e-9, 1000, 1, 'textural volume of macropores')]),
    'KL_engineerParams': dict(groupwise=False, entries=[
        Entry('ratioPVBeng', r'ratio', 'l g C-1', 0, INF, 1, 'burrow volume per engineer biomass'),
        Entry('fPVB', r'fpvb', 'fraction', 0, 1, 1, 'fraction of burrow volume that is additional porosity'),
        Entry('tPVB', r'turnover', 'day-1', 0, 1, 1, 'daily turnover of burrows'),
        Entry('PVBmax', r'pvbmax', 'l m-3', 0, 1000, 1, 'maximal burrow volume'),
        Entry('frag', r'fragment', 'fraction', 0, 1, 1, 'fraction of litter eaten by engineers that becomes SOM'),
        Entry('pfaec_sap', r'm?faec.*(sap|non.?eng)', '-', 0, 5, 1, 'sensitivity of the detritivore faeces fraction to C:N'),
        Entry('pfaec_eng', r'm?faec', '-', 0, 5, 1, 'sensitivity of the engineer faeces fraction to C:N'),
        Entry('bioturbRate', r'bioturb', '% g C-1 day-1', 0, 100, 1, '% of SOM mixed downwards per g C engineer per day'),
        Entry('moveRate', r'litter mov', '% g C-1 day-1', 0, 100, 1, '% of litter moved downwards per g C engineer per day')]),
    'KL_runparams': dict(groupwise=False, entries=[
        Entry('tStop', r'days', 'day', 1, INF, 1, 'number of simulated days'),
        Entry('initWater', r'water', '%', 0, 100, 1, 'initial water content, % of pore volume'),
        Entry('Nmin0', r'n mineral|mineral n', 'g N m-3', 0, INF, 1, 'initial mineral N'),
        Entry('rrg', r'root growth', 'g C m-3 day-1', 0, INF, 1, 'daily root growth'),
        Entry('rootTO', r'root turnover', 'day-1', 0, 1, 1, 'root turnover'),
        Entry('inputLit', r'litter input', 'g C m-3 day-1', 0, INF, 1, 'daily litter input'),
        Entry('CNlit', r'cn of litter', '-', 1, 1000, 1, 'C:N of the litter input'),
        Entry('recLit', r'recalcitrant', '%', 0, 100, 1, 'recalcitrance of the litter input'),
        Entry('CtoMyc', r'c to myc', 'g C m-3 day-1', 0, INF, 1, 'plant C input to mycorrhiza'),
        Entry('NmyctoPlant', r'n from myc', 'fraction', 0, 1, 1, 'fraction of plant N uptake delivered by mycorrhiza'),
        Entry('ee', r'evapotranspiration', 'fraction', 0, 1, 1, 'effective / potential evapotranspiration')]),
    'KL_initC_Pools': dict(groupwise=False, entries=[
        Entry('bact', r'bact', 'g C m-3', 0, INF), Entry('fungi', r'fung', 'g C m-3', 0, INF),
        Entry('myc', r'myc', 'g C m-3', 0, INF), Entry('bvores', r'bacterivor', 'g C m-3', 0, INF),
        Entry('fvores', r'fungivor', 'g C m-3', 0, INF), Entry('sap', r'saprotroph.*non|detritivor', 'g C m-3', 0, INF),
        Entry('eng', r'engineer', 'g C m-3', 0, INF), Entry('hvores', r'herbivor', 'g C m-3', 0, INF),
        Entry('pred', r'predator', 'g C m-3', 0, INF), Entry('litter', r'litter', 'g C m-3', 0, INF),
        Entry('SOM', r'\bsom\b', 'g C m-3', 1e-9, INF), Entry('roots', r'root', 'g C m-3', 0, INF),
        Entry('CO2', r'respiration|co2', 'g C m-3', 0, INF)]),
    'KL_climateParams': dict(groupwise=False, header=True, entries=[
        Entry('TEMP', r'temperature', 'degC', -50, 60, 12, 'monthly mean temperature'),
        Entry('SUNH', r'sun', 'h', 0, 744, 12, 'monthly sunshine hours')]),
}

# entries of the 'PVstruct' array (kept together for convenience)
PV_KEYS = ['PV_inac', 'PV_bact', 'PV_micro', 'PV_meso', 'PV_macro']


def _read_lines(path):
    """(values, comment, line number) for every non-empty line, values as list of floats."""
    out = []
    with open(path, encoding='latin1') as f:
        for k, raw in enumerate(f, start=1):
            line = raw.strip()
            if not line:
                continue
            text, _, comment = line.partition('#')
            values = []
            for tok in text.split():
                try:
                    values.append(float(tok))
                except ValueError:
                    raise ParameterFileError('%s line %d: cannot read number %r' % (path, k, tok))
            out.append((values, comment.strip(), k))
    return out


def _check_group_header(path, comment, lineno):
    tokens = comment.lower().split()
    if len(tokens) != NGROUPS:
        raise ParameterFileError('%s line %d: header must name the %d functional groups (%s), found %r'
                                 % (path, lineno, NGROUPS, ' '.join(GROUP_KEYS), comment))
    for key, tok in zip(GROUP_KEYS, tokens):
        if not any(tok.startswith(a) for a in GROUP_ALIASES[key]):
            raise ParameterFileError('%s line %d: column %r found where %r was expected; the group order must be %s'
                                     % (path, lineno, tok, key, ' '.join(GROUP_KEYS)))


def load_parameter_file(name, directory='.'):
    """Read KL_<name>.txt according to the schema and return {key: value or array}."""
    if name not in SCHEMA:
        raise ParameterFileError('no schema for parameter file %r' % name)
    spec = SCHEMA[name]
    path = os.path.join(directory, name + '.txt')
    if not os.path.exists(path):
        raise ParameterFileError('parameter file %s not found' % path)
    lines = _read_lines(path)
    # header line (group names / month names): a comment without values
    if spec.get('groupwise') or spec.get('header'):
        if not lines or lines[0][0]:
            raise ParameterFileError('%s: first line must be the header (a # comment)' % path)
        header = lines.pop(0)
        if spec.get('groupwise'):
            _check_group_header(path, header[1], header[2])
    entries = spec['entries']
    if len(lines) != len(entries):
        raise ParameterFileError('%s: expected %d parameter lines (%s), found %d'
                                 % (path, len(entries), ', '.join(e.key for e in entries), len(lines)))
    result = {}
    for entry, (values, comment, lineno) in zip(entries, lines):
        if len(values) != entry.n:
            raise ParameterFileError('%s line %d (%s): expected %d value(s), found %d'
                                     % (path, lineno, entry.key, entry.n, len(values)))
        if not re.search(entry.keyword, comment, flags=re.IGNORECASE):
            raise ParameterFileError('%s line %d should be %r (%s) but is labelled %r; the order of the lines must be: %s'
                                     % (path, lineno, entry.key, entry.description or entry.unit, comment,
                                        ', '.join(e.key for e in entries)))
        arr = np.array(values, dtype=float)
        if np.any(arr < entry.lo) or np.any(arr > entry.hi):
            raise ParameterFileError('%s line %d (%s): value(s) %s outside the valid range [%g, %g] %s'
                                     % (path, lineno, entry.key, values, entry.lo, entry.hi, entry.unit))
        result[entry.key] = arr if entry.n > 1 else float(arr[0])
    return result


def load_all(directory='.', density_dependent_turnover=True):
    """All parameter files as one dict of dicts, keyed by file name (without KL_ and .txt)."""
    names = ['KL_FaunalParams', 'KL_initSoil', 'KL_engineerParams', 'KL_runparams', 'KL_initC_Pools', 'KL_climateParams']
    if density_dependent_turnover:
        names.append('KL_densityParams')
    return {n[3:]: load_parameter_file(n, directory) for n in names}


def describe(name):
    """Human readable description of the expected content of a parameter file."""
    spec = SCHEMA[name]
    lines = ['%s.txt (%s):' % (name, 'one column per functional group: ' + ' '.join(GROUP_KEYS)
                                if spec.get('groupwise') else 'one value per line')]
    for i, e in enumerate(spec['entries'], start=1):
        rng = '' if (e.lo == -INF and e.hi == INF) else ' [%g, %g]' % (e.lo, e.hi)
        lines.append('  %2d. %-12s %s%s  %s' % (i, e.key, e.unit, rng, e.description))
    return '\n'.join(lines)
