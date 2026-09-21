"""
Structure of the code base: the model is split into modules by concern and the dependency
direction between them is fixed. hydrology and plant depend on nothing else of the model,
structure only on the names, foodweb on names + plant, core on all of them; keylink_functions is a
pure re-export facade. Nothing below core may import core, matplotlib or read files.

Run with:  python run_tests.py   (all tests) from the model directory.
"""
import ast
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
        self.assertEqual(len(mf.__all__), 28)

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
