"""
Structure of the code base: the model is split into modules by concern and the dependency
direction between them is fixed. hydrology and plant depend on nothing else of the model,
structure only on the names, foodweb on names + plant, core on all of them; keylink_functions is a
pure re-export facade. Nothing below core may import core, matplotlib or read files.

Run with:  python run_tests.py   (all tests) from the model directory.
"""
import ast
import os
import re
import unittest

from tests import MODEL_DIR   # noqa: F401  (sets cwd / sys.path)
import keylink_functions as mf
import keylink_hydrology as hydro
import keylink_structure as struct
import keylink_foodweb as foodweb
import keylink_plant as plant

ALLOWED_IMPORTS = {
    'keylink_names': set(),
    'keylink_params': {'keylink_names'},
    'keylink_plant': set(),
    'keylink_hydrology': set(),
    'keylink_structure': {'keylink_names'},
    'keylink_foodweb': {'keylink_names', 'keylink_plant'},
    'keylink_core': {'keylink_names', 'keylink_params', 'keylink_hydrology', 'keylink_structure',
                     'keylink_foodweb', 'keylink_plant'},
    'keylink_functions': {'keylink_hydrology', 'keylink_structure', 'keylink_foodweb', 'keylink_plant'},
}
FORBIDDEN_BELOW_CORE = {'matplotlib', 'scipy'}   # equations must stay independent of plotting/integration


def imported_modules(filename):
    tree = ast.parse(open(filename, encoding='utf-8').read())
    names = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            names.update(a.name.split('.')[0] for a in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            names.add(node.module.split('.')[0])
    return names


class TestModuleLayering(unittest.TestCase):

    def test_dependency_direction(self):
        for module, allowed in ALLOWED_IMPORTS.items():
            imports = {n for n in imported_modules(module + '.py') if n.startswith('keylink_')}
            self.assertTrue(imports <= allowed, '%s imports %s, allowed: %s' % (module, sorted(imports - allowed), sorted(allowed)))

    def test_equation_modules_do_not_depend_on_plotting_or_solver(self):
        for module in ('keylink_hydrology', 'keylink_structure', 'keylink_foodweb', 'keylink_plant', 'keylink_names', 'keylink_params'):
            imports = imported_modules(module + '.py')
            self.assertFalse(imports & FORBIDDEN_BELOW_CORE, '%s imports %s' % (module, imports & FORBIDDEN_BELOW_CORE))

    def test_facade_reexports_the_same_objects(self):
        for name in mf.__all__:
            obj = getattr(mf, name)
            home = {'keylink_hydrology': hydro, 'keylink_structure': struct, 'keylink_foodweb': foodweb,
                    'keylink_plant': plant}[obj.__module__]
            self.assertIs(obj, getattr(home, name))
        self.assertEqual(len(mf.__all__), 22)

    def test_archive_is_not_imported_by_the_model(self):
        """archive/ holds superseded code and unused functions; no model module may depend on it."""
        for module in list(ALLOWED_IMPORTS):
            self.assertNotIn('keylink_unused', imported_modules(module + '.py'), module)
        self.assertFalse(any(n.startswith('archive') for m in ALLOWED_IMPORTS for n in imported_modules(m + '.py')))

    def test_each_function_lives_in_the_module_of_its_concern(self):
        self.assertEqual(hydro.calcPW.__module__, 'keylink_hydrology')
        self.assertEqual(hydro.PET.__module__, 'keylink_hydrology')
        self.assertEqual(struct.calcAvail.__module__, 'keylink_structure')
        self.assertEqual(struct.calcPVD.__module__, 'keylink_structure')
        self.assertEqual(struct.calcBioturb.__module__, 'keylink_structure')
        self.assertEqual(foodweb.fCompSpecies.__module__, 'keylink_foodweb')
        self.assertEqual(foodweb.calcDeath.__module__, 'keylink_foodweb')
        self.assertEqual(plant.mycNtoPlant.__module__, 'keylink_plant')


if __name__ == '__main__':
    unittest.main()


class TestDocumentation(unittest.TestCase):
    """Every public function/class of the model must have a docstring, and the equation map in
    docs/equations.md must only refer to functions that exist."""

    MODULES = ['keylink_names', 'keylink_params', 'keylink_plant', 'keylink_hydrology',
               'keylink_structure', 'keylink_foodweb', 'keylink_core']

    def test_all_functions_and_classes_have_docstrings(self):
        import ast
        for module in self.MODULES:
            tree = ast.parse(open(module + '.py', encoding='utf-8').read())
            self.assertTrue(ast.get_docstring(tree), '%s has no module docstring' % module)
            for node in ast.walk(tree):
                if isinstance(node, (ast.FunctionDef, ast.ClassDef)):
                    self.assertTrue(ast.get_docstring(node), '%s.%s has no docstring' % (module, node.name))

    def test_units_are_documented(self):
        """The docstrings of the flux/rate functions must state their units."""
        import keylink_structure as st
        import keylink_foodweb as fw
        import keylink_hydrology as hy
        for fn, unit in [(fw.calcgrowth, 'g C m-3 day-1'), (fw.calcDeath, 'g C m-3 day-1'),
                         (fw.calcresp, 'g C g C-1 day-1'), (st.calcBioturb, 'g C m-3 day-1'),
                         (st.calcPVD, 'l m-3'), (hy.PET, 'mm day-1')]:
            self.assertIn(unit, fn.__doc__, fn.__name__)

    def test_equation_map_refers_to_existing_functions(self):
        import importlib
        text = open(os.path.join('docs', 'equations.md'), encoding='utf-8').read()
        refs = set(re.findall(r'`(keylink_\w+)\.(\w+)`', text))
        self.assertGreater(len(refs), 15)
        for module, attr in sorted(refs):
            if module.startswith('keylink_') and module not in ('keylink_unused',):
                mod = importlib.import_module(module)
                self.assertTrue(hasattr(mod, attr), '%s.%s in docs/equations.md does not exist' % (module, attr))

    def test_equation_map_covers_all_equations(self):
        """Every equation 1-37 of Flores et al. (2021) must appear in the map, either with its own
        row or inside a range such as '4-8'."""
        text = open(os.path.join('docs', 'equations.md'), encoding='utf-8').read()
        covered = set()
        for spec in re.findall(r'^\|\s*([0-9\u2013\-]+)\s*\|', text, flags=re.M):
            parts = re.split(r'[\u2013\-]', spec)
            if len(parts) == 2 and all(p.isdigit() for p in parts):
                covered.update(range(int(parts[0]), int(parts[1])+1))
            elif spec.isdigit():
                covered.add(int(spec))
        self.assertEqual(sorted(set(range(1, 38)) - covered), [], 'equations missing in docs/equations.md')
